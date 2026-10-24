# Alpha Shape

A Fusion 360 add-in that wraps the selected bodies — typically the margin cylinders
of screws and pins — in a **lofted solid**: sections perpendicular to a direction,
each one the outline of all the objects with the gaps between them closed, joined
by a loft.

**Solid › Create › Alpha Shape** (also Mesh › Create).

| Input | |
|---|---|
| Bodies | the objects to wrap (solid, surface or mesh) |
| Direction | optional: a planar face or plane (its normal), an axis, a straight edge or a cylinder. Empty = the mean axis of the cylinders among the bodies |
| Closing radius | in each section, gaps narrower than about twice this are filled with round arcs of this radius (0 = just the union) |
| Planes | number of sections |
| Range | **Where every object is cut** (default) or **Whole extent of the objects** |
| Bottom / Top offset | shorten the range at either end |
| Join the objects to the loft | off (default): the result is just the loft of the sections. On: the selected solids (kept) are joined to it, so what sticks out of the range is part of one body |
| Bridge planes that fall apart | see below (on) |
| Preview | the section contours and their matching lines, live |

## How

1. N planes perpendicular to the direction span the range.
2. Each plane cuts every body exactly (`TemporaryBRepManager.planeIntersection`; mesh
   bodies are sliced in the engine).
3. In each plane the cuts are united and **closed** with a disk of radius r (dilate
   r, erode r; shapely), and the outer contour is kept. Closing is the alpha hull
   of the cuts: it fills the gaps a disk of radius r cannot pass, with round arcs
   instead of the straight chords and corners of an alpha shape, so a spline through
   it does not ripple. Arcs are subdivided to a 0.002 mm chord error at any r.
4. Every contour is resampled to the same number of points by arc length, and each
   one is rotated to the cyclic shift that best matches the section below, so equal
   indices correspond and the loft does not twist.
5. The contours become closed fitted splines and a solid loft joins them; optionally
   the selected solids are joined in.

## Why "where every object is cut"

Objects that end at different heights leave planes that cut only some of them,
too far apart for the closing radius. On the Dlab_155 margins the three AK pins
(Ø 4 × 53 mm) run 13–20 mm past the screw margins (Ø 7.2 × 24–26 mm): over the
whole extent those planes needed 36 mm of bridging, the loft creased, and with 32
planes it failed outright. Over the 22.6 mm every margin crosses, the bridges
were 3.5–6 mm and 12 and 24 planes gave the same volume within 0.5 %; joining the
cylinders back in puts the protruding pins into the one body.

**Bridge planes that fall apart:** a plane whose cuts stay separate islands at the
closing radius gets its own radius — the least that joins it, plus 15 % so the
bridge is not a pinch — and the report lists them. Off, Create refuses and reports
the radius that would join every plane.

End planes that only graze a tilted cap (under 5 % of the median section area) are
dropped; the loft cannot use slivers.

## Measured

| | |
|---|---|
| two Ø 2.4 circles 6 mm apart | apart at r = 1; joined with a waist at r = 4; at r = 100, the hull less the waist's sag (0.36 mm²) |
| three test pins, one tilted 1.5° | solid of 3 faces after aligning the matching points (before: a crease where the start point jumped) |
| 12 Dlab_155 margins, shared range, 16 planes | one joined solid, 11 232 mm³; sections in 0.65 s |

## Fusion traps found here

- `planeIntersection` with a plane that misses the body **raises**
  (`InternalValidationError : result`), and the error stays pending and resurfaces
  on the next API call (`TemporaryBRepManager.get()` included). Only planes inside
  each body's own extent are asked for.
- `displayMeshes.bestMesh` raises `count > 0` on a body not drawn yet (hidden, or
  its document just activated); a coarse `createMeshCalculator` is used instead.
- `MeshBodies` has no `itemByName`; `BaseFeature.meshBodies` is a plain list.

## Install

```
powershell -ExecutionPolicy Bypass -File install.ps1
```

then in Fusion: Utilities › Add-Ins › Scripts and Add-Ins › Add-Ins › Alpha Shape ›
Run (tick *Run on Startup*). The engine needs a Python with numpy and shapely
(`pip install shapely`); the add-in looks for `ALPHASHAPE_PYTHON`, then
`python_path.txt` beside it, then the usual install folders, then `py -3.10`.

## Tests

```
py -3.10 tests/test_engine.py
```
