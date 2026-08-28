"""
Sri Lanka Peaks Atlas - Step 20: detect every peak from the DEM.

OSM lists 213 named peaks for the whole island - nowhere near the real count,
and biased toward whatever somebody bothered to tag. This derives peaks from
the terrain itself and ranks them by TOPOGRAPHIC PROMINENCE, which is what
separates an independent mountain from a bump on the shoulder of a bigger one.

ALGORITHM (level-set flooding)

  Imagine the sea rising to the top of the island, then falling. At any water
  level each dry patch is one connected component. As the water drops the
  components grow and eventually merge; the level at which two of them merge is
  the key saddle between their summits, so the lower summit has

      prominence = its elevation - that merge level

  and is thereafter absorbed by the higher one. Sweeping the level from the
  island high point down to MIN_ELEV therefore prices every summit in one pass.

  A summit enters the sweep simply when the falling level crosses its
  elevation, which is why the local maxima are found ONCE up front. The obvious
  alternative - asking scipy for the maximum inside every component at every
  level - rescans all 76 M cells 243 times and takes hours; this takes minutes.

  Prominence is exact to +/- STEP. The island high point never merges with
  anything, so its prominence is its full height above the sea.

Peaks are matched against the OSM natural=peak layer so the known ones keep
their real names.
"""
import csv
import json
import os
import time
import numpy as np
from osgeo import gdal, osr
from scipy import ndimage

gdal.UseExceptions()

SRIPADA = r"C:\Users\thush\OneDrive\Desktop\SriPada"
ROOT = r"C:\Users\thush\OneDrive\Desktop\SriLankaPeaks"
DEM = os.path.join(SRIPADA, "DEM", "SriLanka_DEM_30m.tif")
LM = os.path.join(SRIPADA, "DEM", "SriLanka_landmask_30m.tif")
OSM_PEAKS = os.path.join(SRIPADA, "OSM", "SriLanka_peaks.geojson")
OUT = os.path.join(ROOT, "Results")

MIN_ELEV = 80.0         # flood floor (was 100; lower = bigger working array)
STEP = 10.0             # prominence resolution
MIN_PROM = 50.0         # a summit worth listing (was 100 m)
MIN_PEAK_ELEV = 200.0   # ... and tall enough to call a hill (was 300 m)
NAME_TOL_M = 900.0


def main():
    t0 = time.time()
    ds, lm_ds = gdal.Open(DEM), gdal.Open(LM)
    gt = ds.GetGeoTransform()
    nx, ny = ds.RasterXSize, ds.RasterYSize

    print("loading DEM ...", flush=True)
    z = np.zeros((ny, nx), np.int16)
    BL = 2048
    for r0 in range(0, ny, BL):
        h = min(BL, ny - r0)
        b = ds.GetRasterBand(1).ReadAsArray(0, r0, nx, h)
        l = lm_ds.GetRasterBand(1).ReadAsArray(0, r0, nx, h)
        z[r0:r0 + h] = np.clip(np.where(l == 1, b, -100.0) * 10.0,
                               -32000, 32000).astype(np.int16)
    del b, l

    hi = z >= int(MIN_ELEV * 10)
    rows = np.where(hi.any(axis=1))[0]
    cols = np.where(hi.any(axis=0))[0]
    R0, C0 = int(rows.min()), int(cols.min())
    zc = np.ascontiguousarray(z[R0:int(rows.max()) + 1, C0:int(cols.max()) + 1])
    del z, hi
    zmax = int(zc.max())
    print(f"  highland bbox {zc.shape[1]}x{zc.shape[0]} "
          f"({zc.size/1e6:.0f} M cells), high point {zmax/10:.1f} m "
          f"({time.time()-t0:.0f}s)", flush=True)

    # ---- local maxima, once
    t1 = time.time()
    mx = ndimage.maximum_filter(zc, size=3, mode="nearest")
    cand = (zc == mx) & (zc >= int(MIN_PEAK_ELEV * 10))
    del mx
    lab_c, nc = ndimage.label(cand, structure=np.ones((3, 3), bool))
    pos = ndimage.maximum_position(zc, lab_c, np.arange(1, nc + 1))
    pos = np.array(pos, dtype=np.int64).reshape(-1, 2)
    pel = zc[pos[:, 0], pos[:, 1]].astype(np.int32)
    del cand, lab_c
    order = np.argsort(-pel)
    pos, pel = pos[order], pel[order]
    print(f"  {len(pel):,} local maxima >= {MIN_PEAK_ELEV:g} m "
          f"({time.time()-t1:.0f}s)", flush=True)

    # ---- flood
    struct = np.ones((3, 3), bool)
    mask = np.empty(zc.shape, bool)
    labo = np.empty(zc.shape, np.int32)

    prom = np.full(len(pel), -1, np.int32)      # decimetres
    active = np.zeros(len(pel), bool)
    nxt = 0

    levels = np.arange(zmax, int(MIN_ELEV * 10) - 1, -int(STEP * 10))
    print(f"  flooding {len(levels)} levels ...", flush=True)
    t2 = time.time()

    for li, L in enumerate(levels):
        while nxt < len(pel) and pel[nxt] >= L:
            active[nxt] = True
            nxt += 1
        if not active.any():
            continue

        np.greater_equal(zc, L, out=mask)
        n = ndimage.label(mask, structure=struct, output=labo)
        if n == 0:
            continue

        idx = np.flatnonzero(active)
        comp = labo[pos[idx, 0], pos[idx, 1]]

        o = np.argsort(comp, kind="stable")
        idx_s, comp_s = idx[o], comp[o]
        starts = np.flatnonzero(np.r_[True, np.diff(comp_s) != 0])
        ends = np.r_[starts[1:], len(comp_s)]
        for s, e in zip(starts, ends):
            if e - s < 2:
                continue
            grp = idx_s[s:e]
            keep = grp[np.argmax(pel[grp])]
            for g in grp:
                if g != keep:
                    prom[g] = pel[g] - int(L)
                    active[g] = False

        if li % 20 == 0:
            print(f"    {L/10:7.1f} m  comps {n:6d}  active "
                  f"{int(active.sum()):5d}  {time.time()-t2:5.0f}s", flush=True)

    prom[active] = pel[active]          # never absorbed -> rises from the sea
    print(f"  flooding done ({time.time()-t2:.0f}s)", flush=True)

    # ---- assemble
    sld = osr.SpatialReference(); sld.ImportFromEPSG(5235)
    sld.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    wgs = osr.SpatialReference(); wgs.ImportFromEPSG(4326)
    wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    to_wgs = osr.CoordinateTransformation(sld, wgs)

    keep_i = np.flatnonzero(prom >= int(MIN_PROM * 10))
    recs = []
    for i in keep_i:
        x = gt[0] + (C0 + pos[i, 1] + 0.5) * gt[1]
        y = gt[3] + (R0 + pos[i, 0] + 0.5) * gt[5]
        recs.append({"x": x, "y": y, "elev_m": round(pel[i] / 10.0, 1),
                     "prominence_m": round(prom[i] / 10.0, 1)})
    pts = to_wgs.TransformPoints([(r["x"], r["y"]) for r in recs])
    for r, p in zip(recs, pts):
        r["lon"], r["lat"] = round(p[0], 6), round(p[1], 6)
    recs.sort(key=lambda r: -r["prominence_m"])
    print(f"  {len(recs):,} peaks with prominence >= {MIN_PROM:g} m", flush=True)

    named = 0
    if os.path.exists(OSM_PEAKS):
        of = json.load(open(OSM_PEAKS, encoding="utf-8"))["features"]
        ol = [(o["geometry"]["coordinates"][0], o["geometry"]["coordinates"][1])
              for o in of]
        op = osr.CoordinateTransformation(wgs, sld).TransformPoints(ol)
        ox = np.array([p[0] for p in op])
        oy = np.array([p[1] for p in op])
        onm = [o["properties"].get("name", "") for o in of]
        used = set()
        for r in recs:
            d = np.hypot(ox - r["x"], oy - r["y"])
            i = int(np.argmin(d))
            if d[i] <= NAME_TOL_M and onm[i] and i not in used:
                r["name"] = onm[i]
                used.add(i)
                named += 1
            else:
                r["name"] = ""
    print(f"  {named} matched an OSM name", flush=True)
    for r in recs:
        r.pop("x", None)
        r.pop("y", None)

    os.makedirs(OUT, exist_ok=True)
    cols_out = ["name", "lat", "lon", "elev_m", "prominence_m"]
    fp = os.path.join(OUT, "SriLanka_peaks_detected.csv")
    with open(fp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols_out)
        w.writeheader()
        for r in recs:
            w.writerow({k: r.get(k, "") for k in cols_out})
    json.dump({"type": "FeatureCollection", "features": [
        {"type": "Feature",
         "geometry": {"type": "Point", "coordinates": [r["lon"], r["lat"]]},
         "properties": {k: r.get(k, "") for k in cols_out}} for r in recs]},
        open(os.path.join(OUT, "SriLanka_peaks_detected.geojson"), "w",
             encoding="utf-8"), ensure_ascii=False)

    print("\ntop 25 by prominence:")
    print(f"  {'#':>3} {'name':26s}{'elev':>8}{'prom':>8}   lat, lon")
    for i, r in enumerate(recs[:25], 1):
        print(f"  {i:3d} {(r['name'] or '(unnamed)')[:25]:26s}"
              f"{r['elev_m']:8.0f}{r['prominence_m']:8.0f}   "
              f"{r['lat']:.4f}, {r['lon']:.4f}")
    print(f"\nwrote {fp}  ({time.time()-t0:.0f}s total)")


if __name__ == "__main__":
    main()
