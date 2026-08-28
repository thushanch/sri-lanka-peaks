"""
Sri Lanka Peaks Atlas - Step 22: build the app.

Writes peaks_atlas.html next to a Viewsheds/ folder. Viewshed PNGs stay as
separate files rather than base64 blobs: there are dozens of them, and an
<img> loads fine from file:// (only fetch/XHR is blocked there), so the page
still works by double-clicking it.
"""
import base64
import json
import os
import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()

SRIPADA = r"C:\Users\thush\OneDrive\Desktop\SriPada"
ROOT = r"C:\Users\thush\OneDrive\Desktop\SriLankaPeaks"
DEMDIR = os.path.join(SRIPADA, "DEM")
RES = os.path.join(ROOT, "Results")
GRID_DEG = 0.002
MAXDIM = 2600


def png_bytes(rgba):
    h, w, nb = rgba.shape
    mem = gdal.GetDriverByName("MEM").Create("", w, h, nb, gdal.GDT_Byte)
    for i in range(nb):
        mem.GetRasterBand(i + 1).WriteArray(rgba[:, :, i])
    tmp = f"/vsimem/a_{id(rgba)}.png"
    gdal.GetDriverByName("PNG").CreateCopy(tmp, mem, options=["ZLEVEL=9"])
    f = gdal.VSIFOpenL(tmp, "rb")
    gdal.VSIFSeekL(f, 0, 2)
    n = gdal.VSIFTellL(f)
    gdal.VSIFSeekL(f, 0, 0)
    data = gdal.VSIFReadL(1, n, f)
    gdal.VSIFCloseL(f)
    gdal.Unlink(tmp)
    return data


def latlon_bounds(w):
    g = w.GetGeoTransform()
    m = osr.SpatialReference(); m.ImportFromEPSG(3857)
    m.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ll = osr.SpatialReference(); ll.ImportFromEPSG(4326)
    ll.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(m, ll)
    x1 = g[0] + g[1] * w.RasterXSize
    y1 = g[3] + g[5] * w.RasterYSize
    lo1, la1, _ = tr.TransformPoint(g[0], y1)
    lo2, la2, _ = tr.TransformPoint(x1, g[3])
    return [[la1, lo1], [la2, lo2]]


def packed_grid(path, scale=10.0):
    w = gdal.Warp("", gdal.Open(path), format="MEM", dstSRS="EPSG:4326",
                  xRes=GRID_DEG, yRes=GRID_DEG, resampleAlg="bilinear")
    a = w.GetRasterBand(1).ReadAsArray().astype(np.float64)
    a = np.where(np.isfinite(a), a, -3276.0)
    v = np.clip(np.round(a * scale), -32768, 32767).astype(np.int32) + 32768
    rgb = np.zeros((a.shape[0], a.shape[1], 3), np.uint8)
    rgb[:, :, 0] = (v >> 8) & 0xFF
    rgb[:, :, 1] = v & 0xFF
    gt = w.GetGeoTransform()
    data = png_bytes(rgb)
    print(f"  terrain grid {a.shape[1]}x{a.shape[0]}  {len(data)/1e6:5.2f} MB")
    return {"b64": base64.b64encode(data).decode(), "x0": gt[0], "y0": gt[3],
            "dx": gt[1], "dy": gt[5], "scale": scale}


def relief_overlay():
    hs = gdal.Open(os.path.join(DEMDIR, "SriLanka_hillshade_30m_masked.tif"))
    w = gdal.Warp("", hs, format="MEM", dstSRS="EPSG:3857",
                  resampleAlg="bilinear")
    sc = max(w.RasterXSize, w.RasterYSize) / MAXDIM
    if sc > 1:
        w = gdal.Warp("", w, format="MEM", width=int(w.RasterXSize / sc),
                      height=int(w.RasterYSize / sc), resampleAlg="bilinear")
    a = w.GetRasterBand(1).ReadAsArray()
    wgt = w.GetGeoTransform()
    lm = gdal.Warp("", gdal.Open(os.path.join(DEMDIR, "SriLanka_landmask_30m.tif")),
                   format="MEM", dstSRS="EPSG:3857",
                   width=w.RasterXSize, height=w.RasterYSize,
                   outputBounds=[wgt[0], wgt[3] + wgt[5] * w.RasterYSize,
                                 wgt[0] + wgt[1] * w.RasterXSize, wgt[3]],
                   resampleAlg="near")
    lma = lm.GetRasterBand(1).ReadAsArray()
    rgba = np.zeros((a.shape[0], a.shape[1], 4), np.uint8)
    t = (120 + (a.astype(np.int32) - 140) * 0.5).clip(45, 235).astype(np.uint8)
    rgba[:, :, 0] = (t * 0.90).astype(np.uint8)
    rgba[:, :, 1] = (t * 0.97).astype(np.uint8)
    rgba[:, :, 2] = t
    rgba[:, :, 3] = np.where(lma == 1, 255, 0).astype(np.uint8)
    data = png_bytes(rgba)
    print(f"  relief       {a.shape[1]}x{a.shape[0]}  {len(data)/1e6:5.2f} MB")
    return {"b64": base64.b64encode(data).decode(), "bounds": latlon_bounds(w)}


if __name__ == "__main__":
    peaks = json.load(open(os.path.join(RES, "SriLanka_peaks_detected.geojson"),
                           encoding="utf-8"))["features"]
    ALL = []
    for ft in peaks:
        pr = ft["properties"]
        lon, lat = ft["geometry"]["coordinates"]
        ALL.append({"name": pr.get("name", ""), "lat": lat, "lon": lon,
                    "elev_m": pr["elev_m"], "prominence_m": pr["prominence_m"]})
    print(f"{len(ALL)} peaks detected")

    vs_index = []
    vsp = os.path.join(RES, "peak_viewshed_index.json")
    if os.path.exists(vsp):
        vs_index = json.load(open(vsp, encoding="utf-8"))

    # Attach viewsheds BY POSITION, not by list index. The index-based pairing
    # this replaced only held while the peak list never changed; re-running
    # detection reorders it and would hand a mountain someone else's viewshed.
    def attach(plist):
        hits = 0
        for v in vs_index:
            best, bd = None, 1e9
            for p in plist:
                d = (p["lat"] - v["lat"]) ** 2 + (p["lon"] - v["lon"]) ** 2
                if d < bd:
                    bd, best = d, p
            if best is not None and bd < (0.004) ** 2:      # ~400 m
                best["vsid"] = v["id"]
                hits += 1
        return hits

    # Display filter: named summits, plus unnamed ones tall enough to matter.
    # Peaks carrying a precomputed viewshed are kept regardless - dropping them
    # would orphan the most expensive data in the project.
    MIN_UNNAMED_ELEV = 900.0
    attach(ALL)
    P = [p for p in ALL
         if p["name"] or p["elev_m"] >= MIN_UNNAMED_ELEV or p.get("vsid")]
    named = sum(1 for p in P if p["name"])
    withvs = sum(1 for p in P if p.get("vsid"))
    print(f"  shown: {len(P)}  ({named} named, "
          f"{len(P)-named} unnamed >= {MIN_UNNAMED_ELEV:g} m or with a viewshed)")
    print(f"  viewsheds matched to a peak: {withvs} of {len(vs_index)}")

    print("assets ...")
    # keep only viewshed records whose peak survived, so the app never offers
    # a footprint for a mountain that is not in its list
    kept = {p.get("vsid") for p in P if p.get("vsid")}
    vs_index = [v for v in vs_index if v["id"] in kept]

    brand = None
    bp = os.path.join(ROOT, "brand.json")
    if os.path.exists(bp):
        brand = json.load(open(bp, encoding="utf-8"))
        print("  brand assets embedded (wordmark + mark)")
    else:
        print("  ! brand.json missing - run 24_branding.py")

    payload = {"peaks": P, "viewsheds": vs_index, "brand": brand,
               "relief": relief_overlay(),
               "grids": {"z": packed_grid(
                   os.path.join(DEMDIR, "SriLanka_DEM_30m.tif"))}}

    sd = os.path.join(ROOT, "Scripts")
    tpl = open(os.path.join(sd, "peaks_app_template.html"),
               encoding="utf-8").read()

    # Inline Leaflet rather than pulling it from a CDN. Without this the app
    # is only "offline" in the sense that its basemap is - the map library
    # itself would still need the network, and the whole page would be blank.
    # (The app uses divIcon/circleMarker only, so Leaflet's default marker
    # PNGs are never requested.)
    for token, fn in (("/*__LEAFLET_CSS__*/", "leaflet.css"),
                      ("/*__LEAFLET_JS__*/", "leaflet.js")):
        p = os.path.join(sd, fn)
        if not os.path.exists(p):
            raise SystemExit(f"missing {p} - fetch it from unpkg first")
        tpl = tpl.replace(token, open(p, encoding="utf-8").read())

    out = os.path.join(ROOT, "peaks_atlas.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(tpl.replace("/*__DATA__*/", json.dumps(payload)))
    print(f"\nwrote {out}  ({os.path.getsize(out)/1e6:.1f} MB)")
