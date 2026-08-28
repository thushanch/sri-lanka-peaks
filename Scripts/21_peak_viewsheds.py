"""
Sri Lanka Peaks Atlas - Step 21: a real 30 m viewshed for every listed peak.

Same physics as the Sri Pada engine - curvature and refraction via
c(d) = (1-k)d^2/(2R), k = 0.13 - run once per peak.

Two things keep this tractable for ~100 peaks instead of ~10 hours:

  * Each peak only reaches as far as its own height allows. A 2192 m summit
    carries 184 km; a 600 m one carries 101 km. Rays stop at that peak's own
    geometric limit, and only the raster window inside it is touched.
  * The DEM is loaded once as int16 decimetres and reused across every peak.

Writes, per peak: a Web-Mercator PNG of the ground-level viewshed (for the
app) plus exact 30 m land-area statistics (for the record - the PNG is
downsampled for display only, the numbers are not).
"""
import base64
import csv
import json
import os
import time
import numpy as np
from osgeo import gdal, osr
from scipy.ndimage import map_coordinates

gdal.UseExceptions()

SRIPADA = r"C:\Users\thush\OneDrive\Desktop\SriPada"
ROOT = r"C:\Users\thush\OneDrive\Desktop\SriLankaPeaks"
DEM = os.path.join(SRIPADA, "DEM", "SriLanka_DEM_30m.tif")
LM = os.path.join(SRIPADA, "DEM", "SriLanka_landmask_30m.tif")
PEAKS = os.path.join(ROOT, "Results", "SriLanka_peaks_detected.geojson")
VSD = os.path.join(ROOT, "Viewsheds")

K, RE, EYE = 0.13, 6371008.8, 1.7
TILE = 1024
PNG_MAXDIM = 1500
TOP_N = 45          # ~85 min; the whole famous set, and resumable if raised
VIS_RGBA = (255, 209, 102, 190)


def geometric_range(z_peak):
    """How far a 1.7 m observer could see this summit over a smooth earth."""
    return 3.86 * (np.sqrt(max(z_peak, 1.0)) + np.sqrt(EYE)) * 1000.0 * 1.05


def png_from_mask(mask, gt, wkt, out_path):
    """Binary viewshed -> Web Mercator RGBA PNG + latlon bounds.

    The MEM band is Float32, not Byte: averaging a 0/1 Byte band rounds each
    output cell straight back to 0 or 1 and the coverage fraction - which is
    what gives the overlay a soft, non-aliased edge - is lost.
    """
    h, w = mask.shape
    mem = gdal.GetDriverByName("MEM").Create("", w, h, 1, gdal.GDT_Float32)
    mem.SetGeoTransform(gt)
    mem.SetProjection(wkt)
    mem.GetRasterBand(1).WriteArray(mask.astype(np.float32))
    warp = gdal.Warp("", mem, format="MEM", dstSRS="EPSG:3857",
                     resampleAlg="average", outputType=gdal.GDT_Float32)
    sc = max(warp.RasterXSize, warp.RasterYSize) / PNG_MAXDIM
    if sc > 1:
        warp = gdal.Warp("", warp, format="MEM",
                         width=max(int(warp.RasterXSize / sc), 1),
                         height=max(int(warp.RasterYSize / sc), 1),
                         resampleAlg="average", outputType=gdal.GDT_Float32)
    a = np.clip(warp.GetRasterBand(1).ReadAsArray(), 0.0, 1.0)
    rgba = np.zeros((a.shape[0], a.shape[1], 4), np.uint8)
    rgba[:, :, 0] = VIS_RGBA[0]
    rgba[:, :, 1] = VIS_RGBA[1]
    rgba[:, :, 2] = VIS_RGBA[2]
    rgba[:, :, 3] = (a * VIS_RGBA[3]).astype(np.uint8)

    png = gdal.GetDriverByName("MEM").Create("", a.shape[1], a.shape[0], 4,
                                             gdal.GDT_Byte)
    for i in range(4):
        png.GetRasterBand(i + 1).WriteArray(rgba[:, :, i])
    gdal.GetDriverByName("PNG").CreateCopy(out_path, png, options=["ZLEVEL=9"])

    g = warp.GetGeoTransform()
    m = osr.SpatialReference(); m.ImportFromEPSG(3857)
    m.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ll = osr.SpatialReference(); ll.ImportFromEPSG(4326)
    ll.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(m, ll)
    x1 = g[0] + g[1] * warp.RasterXSize
    y1 = g[3] + g[5] * warp.RasterYSize
    lo1, la1, _ = tr.TransformPoint(g[0], y1)
    lo2, la2, _ = tr.TransformPoint(x1, g[3])
    return [[la1, lo1], [la2, lo2]]


def main():
    t0 = time.time()
    ds = gdal.Open(DEM)
    lm_ds = gdal.Open(LM)
    gt, wkt = ds.GetGeoTransform(), ds.GetProjection()
    nx, ny = ds.RasterXSize, ds.RasterYSize
    px = gt[1]

    # Only the DEM is held resident. The land mask is read per tile instead:
    # keeping it in RAM costs 150 MB that the polar array needs more.
    print("loading DEM ...", flush=True)
    z16 = np.empty((ny, nx), np.int16)
    BL = 2048
    for r0 in range(0, ny, BL):
        h = min(BL, ny - r0)
        blk = ds.GetRasterBand(1).ReadAsArray(0, r0, nx, h)
        z16[r0:r0 + h] = np.clip(blk * 10.0, -32000, 32000).astype(np.int16)
    del blk
    print(f"  {z16.nbytes/1e6:.0f} MB ({time.time()-t0:.0f}s)", flush=True)

    with open(PEAKS, encoding="utf-8") as f:
        feats = json.load(f)["features"]
    feats = feats[:TOP_N]
    print(f"{len(feats)} peaks to process\n", flush=True)

    wgs = osr.SpatialReference(); wgs.ImportFromEPSG(4326)
    wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    sld = osr.SpatialReference(); sld.ImportFromEPSG(5235)
    sld.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    to_sld = osr.CoordinateTransformation(wgs, sld)

    os.makedirs(VSD, exist_ok=True)
    index = []
    for pi, ft in enumerate(feats):
        pr = ft["properties"]
        lon, lat = ft["geometry"]["coordinates"]
        sx, sy, _ = to_sld.TransformPoint(lon, lat)
        pid = f"p{pi:03d}"
        name = pr.get("name") or f"Peak {pr['elev_m']:.0f} m"
        png_path = os.path.join(VSD, f"{pid}.png")
        meta_path = os.path.join(VSD, f"{pid}.json")

        if os.path.exists(png_path) and os.path.exists(meta_path):
            index.append(json.load(open(meta_path, encoding="utf-8")))
            continue

        tp = time.time()
        col = int((sx - gt[0]) / px)
        row = int((gt[3] - sy) / px)
        zs = float(z16[row, col]) / 10.0
        rmax = geometric_range(zs)
        nr = max(int(rmax / px), 10)
        # capped at 20000: at the longest range this is ~65 m of azimuthal
        # spacing, and the polar array stays near 550 MB rather than 830 MB,
        # which this machine cannot hold alongside the DEM
        naz = int(np.clip(2 * np.pi * rmax / px, 3600, 20000))

        # --- polar sweep from this summit
        r = ((np.arange(nr, dtype=np.float32) + 1.0) * px)
        c_r = (1.0 - K) * r * r / (2.0 * RE)
        Aprev = np.empty((naz + 1, nr), np.float32)
        az = np.arange(naz, dtype=np.float64) * (2.0 * np.pi / naz)
        CH = 512
        for a0 in range(0, naz, CH):
            a1 = min(a0 + CH, naz)
            sa = np.sin(az[a0:a1])[:, None]
            ca = np.cos(az[a0:a1])[:, None]
            cc = (sx + r[None, :] * sa - gt[0]) / px - 0.5
            rr = (gt[3] - (sy + r[None, :] * ca)) / px - 0.5
            zr = map_coordinates(z16, [rr, cc], order=1, mode="constant",
                                 cval=0.0, output=np.float32)
            zr *= 0.1
            al = (zr - c_r[None, :] - zs) / r[None, :]
            np.maximum.accumulate(al, axis=1, out=al)
            Aprev[a0:a1, 0] = -10.0
            Aprev[a0:a1, 1:] = al[:, :-1]
            del cc, rr, zr, al
        Aprev[naz] = Aprev[0]

        # --- solve only the window this peak can actually reach
        c_lo = max(int((sx - rmax - gt[0]) / px), 0)
        c_hi = min(int((sx + rmax - gt[0]) / px) + 1, nx)
        r_lo = max(int((gt[3] - (sy + rmax)) / px), 0)
        r_hi = min(int((gt[3] - (sy - rmax)) / px) + 1, ny)
        W, H = c_hi - c_lo, r_hi - r_lo
        vis = np.zeros((H, W), np.uint8)
        naz_f = naz / (2.0 * np.pi)
        nvis_land = 0

        for ry in range(r_lo, r_hi, TILE):
            th = min(TILE, r_hi - ry)
            yy = gt[3] - (np.arange(ry, ry + th, dtype=np.float64)[:, None] + 0.5) * px
            for rx in range(c_lo, c_hi, TILE):
                tw = min(TILE, c_hi - rx)
                xx = gt[0] + (np.arange(rx, rx + tw, dtype=np.float64)[None, :] + 0.5) * px
                dx = xx - sx
                dy = yy - sy
                d = np.hypot(dx, dy)
                if d.min() > rmax:
                    continue
                a = np.mod(np.arctan2(dx, dy), 2 * np.pi)
                ai = a * naz_f
                ri = np.clip(d / px - 1.0, 0.0, nr - 1.0)
                lo, hi = ai.min(), ai.max()
                if hi - lo > naz * 0.5:
                    sub, off = Aprev, 0.0
                else:
                    i0 = max(int(np.floor(lo)) - 1, 0)
                    i1 = min(int(np.ceil(hi)) + 2, naz + 1)
                    sub, off = Aprev[i0:i1], float(i0)
                astar = map_coordinates(sub, [(ai - off).ravel(), ri.ravel()],
                                        order=1, mode="nearest").reshape(d.shape)
                cd = (1.0 - K) * d * d / (2.0 * RE)
                zg = z16[ry:ry + th, rx:rx + tw].astype(np.float32) * 0.1
                hreq = zs + astar * d + cd - zg
                v = ((hreq <= EYE) & (d <= rmax)).astype(np.uint8)
                vis[ry - r_lo:ry - r_lo + th, rx - c_lo:rx - c_lo + tw] = v
                lmt = lm_ds.GetRasterBand(1).ReadAsArray(rx, ry, tw, th)
                nvis_land += int((v & lmt).sum())
        del Aprev

        sub_gt = (gt[0] + c_lo * px, px, 0, gt[3] - r_lo * px, 0, -px)
        bounds = png_from_mask(vis, sub_gt, wkt, png_path)
        area = nvis_land * (px * px) / 1e6
        rec = {"id": pid, "name": name, "lat": lat, "lon": lon,
               "elev_m": pr["elev_m"], "prominence_m": pr["prominence_m"],
               "dem_elev_m": round(zs, 1),
               "range_km": round(rmax / 1000, 1),
               "visible_land_km2": round(area, 1),
               "visible_pct_of_island": round(100 * area / 65896.0, 2),
               "png": f"Viewsheds/{pid}.png", "bounds": bounds}
        json.dump(rec, open(meta_path, "w", encoding="utf-8"))
        index.append(rec)
        # write the index after every peak, not just at the end: this job runs
        # for hours and the app must be buildable from whatever is finished
        json.dump(index, open(os.path.join(ROOT, "Results",
                                           "peak_viewshed_index.json"), "w",
                              encoding="utf-8"))
        print(f"  [{pi+1:3d}/{len(feats)}] {name[:26]:28s} "
              f"{zs:6.0f} m  range {rmax/1000:5.0f} km  "
              f"sees {area:7,.0f} km2  ({time.time()-tp:4.0f}s)", flush=True)

    idx_path = os.path.join(ROOT, "Results", "peak_viewshed_index.json")
    json.dump(index, open(idx_path, "w", encoding="utf-8"))

    with open(os.path.join(ROOT, "Results", "peak_visibility.csv"), "w",
              newline="", encoding="utf-8-sig") as f:
        cols = ["id", "name", "lat", "lon", "elev_m", "prominence_m",
                "range_km", "visible_land_km2", "visible_pct_of_island"]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(index)

    print(f"\nwrote {idx_path}   ({time.time()-t0:.0f}s total)")
    top = sorted(index, key=lambda r: -r["visible_land_km2"])[:12]
    print("\nmost widely visible peaks:")
    for r in top:
        print(f"  {r['name'][:28]:30s}{r['dem_elev_m']:7.0f} m  "
              f"{r['visible_land_km2']:8,.0f} km2  "
              f"{r['visible_pct_of_island']:5.1f}% of the island")


if __name__ == "__main__":
    main()
