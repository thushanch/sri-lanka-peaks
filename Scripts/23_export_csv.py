"""
Sri Lanka Peaks Atlas - Step 23: export peak_visibility.csv from the index.

21_peak_viewsheds.py only writes this CSV when the whole run finishes, but the
index is written after every peak. This regenerates the CSV from whatever is
done, so the table is never missing just because the long job is still going.
Safe to re-run at any time.
"""
import csv
import json
import os

ROOT = r"C:\Users\thush\OneDrive\Desktop\SriLankaPeaks"
RES = os.path.join(ROOT, "Results")

idx = json.load(open(os.path.join(RES, "peak_viewshed_index.json"),
                     encoding="utf-8"))
idx = sorted(idx, key=lambda v: -v["visible_land_km2"])

cols = ["id", "name", "lat", "lon", "elev_m", "dem_elev_m", "prominence_m",
        "range_km", "visible_land_km2", "visible_pct_of_island"]
out = os.path.join(RES, "peak_visibility.csv")
with open(out, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    w.writerows(idx)

print(f"wrote {out}  ({len(idx)} peaks with a full 30 m viewshed)")
print(f"\n{'peak':32s}{'elev':>7}{'prom':>7}{'sees km2':>11}{'% island':>10}")
for v in idx:
    print(f"  {v['name'][:29]:30s}{v['dem_elev_m']:7.0f}{v['prominence_m']:7.0f}"
          f"{v['visible_land_km2']:11,.0f}{v['visible_pct_of_island']:9.1f}%")
