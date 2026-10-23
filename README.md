# Alpha Shape

A Fusion 360 add-in that wraps the selected bodies in their **alpha shape**: the
solid that a ball of radius *alpha* cannot carve out of the space between them. A
small alpha hugs each object; a large one tends to the convex hull of them all. The
result is a mesh body.

**Mesh › Create › Alpha Shape** (also Solid › Create).

| Input | |
|---|---|
| Bodies | any number of solid, surface or mesh bodies |
| Alpha (ball radius) | gaps narrower than this are filled |
| Resolution | sampling step on the surfaces |
| Preview | the shape is drawn as you change the values |
| Options | repair non-manifold spots (on), drop pieces under 1 % of the volume (on), hide the inputs afterwards |

## How

1. The bodies are tessellated (a tenth of the resolution, never finer than 0.05 mm)
   and their surfaces sampled every *resolution* mm. Sharp edges, such as a
   cylinder's rims, are walked at half the resolution and never thinned, so a rim
   survives any resolution.
2. A 3D Delaunay of the samples (Qhull, joggled: a cylinder's rim points are
   cospherical, and unjoggled Qhull leaves flat tetrahedra that punch holes) and the
   circumradius of every tetrahedron. This is cached on disk per selection and
   resolution, so **changing alpha only re-filters**: about a second on the Dlab_155
   margin cylinders at 0.5 mm.
3. The tetrahedra whose circumscribed ball is no larger than alpha make the solid;
   its boundary, oriented outward, is the shape.
4. Where two parts touch only along an edge or at a point, the tetrahedra there are
   removed until the boundary is a closed 2-manifold. Removing beat adding: adding
   the tetrahedra round a bad edge cascaded, put up to 23 % on the volume and never
   converged; removing converged in 1-5 rounds for -0.2 % at most.

The report gives the sample and triangle counts, volume, area, pieces and whether
the result is watertight.

## Measured

| | |
|---|---|
| one Ø 2.4 × 16 mm pin, alpha 3, resolution 0.4 | closed, volume within 3 % of the 24-gon |
| same, alpha 1000 | the convex hull, within 1 % |
| two such pins 6 mm apart, alpha 1.5 / 8 / 1000 | two pieces / one bridged solid / the hull within 2 % |
| the 12 margin cylinders of Dlab_155, resolution 0.5 | 32 663 points; alpha 3, 6, 15 mm all watertight, one piece; ~1 s per new alpha |
| same at resolution 0.25 | 135 258 points, 7.5 s for the new Delaunay |

Below an object's own radius the shape falls apart into shreds: there is no ball
small enough to fit inside it. Keep alpha above the thinnest object's radius and well
above the resolution.

## Install

```
powershell -ExecutionPolicy Bypass -File install.ps1
```

then in Fusion: Utilities › Add-Ins › Scripts and Add-Ins › Add-Ins › Alpha Shape ›
Run (tick *Run on Startup*). `as_core.py` reloads on every click and `as_engine.py`
runs fresh per call; only `AlphaShape.py` needs a Stop/Run.

The engine needs a Python with numpy and scipy. The add-in looks for
`ALPHASHAPE_PYTHON`, then `python_path.txt` beside it, then the usual install folders,
then `py -3.10`.

## Tests

```
py -3.10 tests/test_engine.py
```
