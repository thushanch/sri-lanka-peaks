"""
Prepare Thushan Chamika brand assets for embedding in both web apps.

Source of truth is the Desktop Branding folder. Per the brand guide this uses
REVERSE (white) lockup, because both apps are dark-surfaced, and keeps the
branding in the app chrome rather than over the map or the panorama.

Writes brand.json (base64 PNGs) into both project folders:
  wordmark  white horizontal lockup, trimmed to its ink and scaled down
  mark      the peak-over-water glyph with its off-white card knocked out,
            for the favicon
"""
import base64
import json
import os
import numpy as np
from osgeo import gdal

gdal.UseExceptions()

BRAND = r"C:\Users\thush\OneDrive\Desktop\Branding"
TARGETS = [r"C:\Users\thush\OneDrive\Desktop\SriLankaPeaks",
           r"C:\Users\thush\OneDrive\Desktop\SriPada"]


def read_rgba(path):
    ds = gdal.Open(path)
    a = ds.ReadAsArray()
    if a.ndim == 2:
        a = np.stack([a, a, a])
    a = np.transpose(a, (1, 2, 0))
    if a.shape[2] == 3:
        a = np.dstack([a, np.full(a.shape[:2], 255, np.uint8)])
    return a.astype(np.uint8)


def to_png_b64(rgba):
    h, w, _ = rgba.shape
    mem = gdal.GetDriverByName("MEM").Create("", w, h, 4, gdal.GDT_Byte)
    for i in range(4):
        mem.GetRasterBand(i + 1).WriteArray(rgba[:, :, i])
    tmp = f"/vsimem/b_{id(rgba)}.png"
    gdal.GetDriverByName("PNG").CreateCopy(tmp, mem, options=["ZLEVEL=9"])
    f = gdal.VSIFOpenL(tmp, "rb")
    gdal.VSIFSeekL(f, 0, 2)
    n = gdal.VSIFTellL(f)
    gdal.VSIFSeekL(f, 0, 0)
    data = gdal.VSIFReadL(1, n, f)
    gdal.VSIFCloseL(f)
    gdal.Unlink(tmp)
    return base64.b64encode(data).decode(), n


def trim_alpha(rgba, pad=2):
    al = rgba[:, :, 3]
    ys, xs = np.where(al > 8)
    if not len(ys):
        return rgba
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad + 1, rgba.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad + 1, rgba.shape[1])
    return rgba[y0:y1, x0:x1]


def resize(rgba, out_h):
    """Plain pixel resize. Translate, not Warp: these are images with no
    geotransform, and Warp refuses to work without one."""
    h, w, _ = rgba.shape
    out_w = max(1, int(round(w * out_h / h)))
    mem = gdal.GetDriverByName("MEM").Create("", w, h, 4, gdal.GDT_Byte)
    for i in range(4):
        mem.GetRasterBand(i + 1).WriteArray(rgba[:, :, i])
    out = gdal.Translate("", mem, format="MEM", width=out_w, height=out_h,
                         resampleAlg="lanczos")
    a = np.transpose(out.ReadAsArray(), (1, 2, 0))
    return np.clip(a, 0, 255).astype(np.uint8)


if __name__ == "__main__":
    # ---- wordmark: white lockup, trimmed, small
    wm = read_rgba(os.path.join(BRAND, "wm_horizontal_white.png"))
    wm = resize(trim_alpha(wm), 72)
    wm_b64, wm_n = to_png_b64(wm)
    print(f"wordmark  {wm.shape[1]}x{wm.shape[0]}  {wm_n/1024:6.1f} KB")

    # ---- mark: knock out the off-white card so it sits on any surface
    mk = read_rgba(os.path.join(BRAND, "markB.png"))
    rgb = mk[:, :, :3].astype(np.int16)
    # brand paper is #F6F7F4 / #F4F5F1; treat anything near it as background
    near_paper = (np.abs(rgb - np.array([246, 247, 244])).max(axis=2) < 26)
    mk[:, :, 3] = np.where(near_paper, 0, 255).astype(np.uint8)
    mk = resize(trim_alpha(mk), 128)
    mk_b64, mk_n = to_png_b64(mk)
    kept = 100 * np.count_nonzero(mk[:, :, 3] > 8) / mk[:, :, 3].size
    print(f"mark      {mk.shape[1]}x{mk.shape[0]}  {mk_n/1024:6.1f} KB  "
          f"{kept:.0f}% opaque after knockout")

    payload = {"wordmark": wm_b64, "mark": mk_b64,
               "name": "Thushan Chamika",
               "palette": {"deep": "#14416B", "water": "#1E78B0",
                           "green": "#2E6B4F", "offwhite": "#F6F7F4",
                           "ink": "#1C2A33"}}
    for t in TARGETS:
        p = os.path.join(t, "brand.json")
        json.dump(payload, open(p, "w", encoding="utf-8"))
        print("wrote", p)
