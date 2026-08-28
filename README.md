# Sri Lanka Peaks Atlas

A separate app from the Sri Pada study, built on the same 30 m terrain and the
same physics. It answers visibility in **both directions**:

* **From a point** — click anywhere and see the skyline you would actually
  look at, with every catalogued peak that breaks it numbered and named.
* **From a peak** — select a summit and see every place in Sri Lanka it can be
  seen from.

Open `peaks_atlas.html` in a browser. It needs no network at all — Leaflet is
inlined and the default basemap is baked from the study's own hillshade; the
online basemaps are optional extras. It works on a phone.

Built for general use, not just for GIS people: a first-run explainer, a
locate-me button, search by name or number, shareable links that restore the
exact spot, bearing and eye height, a drag-to-turn panorama with a sunrise-light
mode, and a mobile layout where the view fills the screen and the peak list
becomes a swipeable strip. Every peak carries one number, used identically on
the map, in the list and on the horizon, so nothing is buried under overlapping
name labels.

## What is in it

| | |
|---|---|
| Peaks detected from the DEM | **4,406** with prominence ≥ 50 m |
| Matched to an OSM name | 146 |
| Shown in the app | **809** — every named summit, plus unnamed ones ≥ 900 m |
| Peaks with a full 30 m viewshed | see `Results/peak_visibility.csv` |
| Island high point | Pidurutalagala, 2,521 m |

Every peak carries a **number** — the same number on the map, in the peak list,
and on the panorama — so the map stays readable without hundreds of overlapping
name labels.

## How the peaks were found

Elevation alone is a poor way to rank mountains: a 1,500 m bump on the shoulder
of a 2,000 m massif is not a mountain. The atlas ranks by **topographic
prominence** instead, computed by level-set flooding (`20_peaks_detect.py`):

> Imagine the sea rising to the island high point and then falling. At any
> water level each dry patch is one connected component. As the water drops,
> components merge; the level at which two merge is the key saddle between
> their summits, so the lower summit has `prominence = its elevation − that
> merge level`, and is thereafter absorbed by the higher one.

Sweeping from 2,521 m down to 80 m in 10 m steps prices every summit in one
pass. Local maxima are found **once** up front, so a summit enters the sweep
simply when the level crosses its elevation — asking scipy for the maximum
inside every component at every level instead rescans 76 M cells 243 times and
turns 8 minutes into 4½ hours.

Prominence is exact to ±10 m. The island high point never merges with anything,
so its prominence equals its full height above the sea — which is the check
that the algorithm is wired up correctly.

## Visibility

Identical to the Sri Pada engine. Rays are cast from the summit, and along each
one the running maximum of the vertical angle

    alpha(r) = ( z(r) − c(r) − Zs ) / r ,   c(d) = (1−k) d² / (2R)

gives the grazing angle of the highest obstruction. **Earth curvature and
atmospheric refraction are both modelled** (k = 0.13, R = 6,371,008.8 m),
applied to every terrain sample and to the target cell. This is not cosmetic:
it moves the sea-level horizon for a 2,192 m peak from 172 km to 184 km.

Each peak only reaches as far as its own height allows — `3.86·(√h + √1.7)` km —
so rays stop there and only the raster window inside that radius is touched.

## The panorama

Terrain is drawn by walking each column outward and painting a strip **only
where the running-maximum angle rises**. That gives correct occlusion for free,
and lets each strip be shaded by the real DEM gradient (sun at 315°, 45°) and
outlined at its crest. Nearer ridges paint over farther ones, and colour blends
toward the sky with distance.

**True eye view.** One pixel is the same number of degrees across as it is up,
so the picture matches what the eye actually receives. An earlier version fitted
the vertical axis to whatever was on screen, which came out ~6.6× taller than
life — mountains felt far too close — and clipped the view at about −3°, so from
a summit you could not see your own slopes fall away. At 1:1 the view spans
~35° vertically instead of 6.7°, and reaches 17° below the horizon. **Height
scale** breaks the 1:1 on purpose when distant hills need to be legible, and
**Tilt** looks up or down.

Controls: eye height (0–200 m, so the tall-building effect can be felt rather
than read), lens width, bearing, tilt, height scale, a Fit button, and a 360°
sweep. Drag the view to turn and tilt.

Two things keep it responsive: sampling steps along a bearing use a local
planar step instead of spherical trigonometry (~30 m drift over 190 km, far
inside the 220 m grid), and sub-pixel strips are coalesced before painting.
Together those took a sea-level frame from 15.3 s to under 1 s.

## Honest limits

- **The DEM is a DSM.** Canopy is included: it blocks realistically, but it
  also lifts forested observers to treetop height.
- **Panorama ridges are drawn from a ~220 m resample** of the 30 m DEM, because
  the browser cannot hold 150 M cells. The *verdict* (visible / blocked /
  height needed) always comes from the full 30 m raster; only the silhouette is
  approximate. Peak markers use true 30 m summit elevations, so a marker can
  sit slightly above its resampled ridge.
- **Haze is not modelled.** Distance shading is a visual cue, not radiative
  transfer. Anything past ~100 km is geometric visibility only.
- **Names come from OSM** and only 146 of 4,406 peaks have one; the rest are
  genuine summits that nobody has tagged.

## Scripts

| | |
|---|---|
| `20_peaks_detect.py` | prominence by level-set flooding → 4,406 peaks |
| `21_peak_viewsheds.py` | one 30 m viewshed per peak → PNG + area stats (resumable) |
| `22_build_app.py` | bundles everything into `peaks_atlas.html` |

Run under the QGIS python:
`"C:\Program Files\QGIS 3.40.8\bin\python-qgis-ltr.bat"`
