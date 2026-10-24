"""Alpha Shape engine - section contours for a lofted envelope. numpy + shapely.

    python as_engine.py job.json out.json

For every section plane: the cuts of all the objects (closed loops in the plane's
own u, v coordinates, mm) are united, CLOSED with a disk of radius r - dilated by
r, then eroded by r - and the outer contour of the result is kept. Closing fills
every gap narrower than 2r, like an alpha shape with alpha = r, but leaves round
arcs of radius r where an alpha shape leaves straight chords and corners, so a
spline through it does not ripple and the loft does not wrinkle.

All contours are then resampled to the SAME number of points, by arc length,
counter-clockwise, starting where a ray from the contour's centroid along +u
leaves it: equal indices correspond from section to section, and the loft has
nothing to twist about.

Job:
    sections   [{'h': mm, 'loops': [[[u, v], ...], ...]}, ...]   (BRep cuts)
    meshes     optional [{'coords', 'idx', 'origin', 'e1', 'e2', 'd'}] - sliced here
    r_mm       closing radius
    spacing_mm target spacing of the resampled points
"""

import json
import math
import sys
import time

import numpy as np
from shapely.geometry import Polygon, MultiPolygon, LineString
from shapely.ops import unary_union, polygonize
from shapely.geometry.polygon import orient


def slice_mesh(m, heights):
    """Closed loops (u, v) where each plane h cuts a triangle mesh."""
    v = np.asarray(m['coords'], float).reshape(-1, 3)
    f = np.asarray(m['idx'], np.int64).reshape(-1, 3)
    o, e1, e2, d = (np.asarray(m[k], float) for k in ('origin', 'e1', 'e2', 'd'))
    hv = (v - o) @ d
    uv = np.c_[(v - o) @ e1, (v - o) @ e2]
    out = []
    for h in heights:
        s = hv - h
        segs = []
        for tri in f:
            pts = []
            for p, q in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                if (s[p] < 0) != (s[q] < 0):
                    t = s[p] / (s[p] - s[q])
                    pts.append(uv[p] + t * (uv[q] - uv[p]))
            if len(pts) == 2:
                segs.append(LineString(pts))
        out.append([list(p.exterior.coords) for p in polygonize(segs)] if segs else [])
    return out


def close_section(loops, r):
    polys = []
    for lp in loops:
        if len(lp) >= 3:
            p = Polygon(lp).buffer(0)
            if not p.is_empty:
                polys.append(p)
    if not polys:
        return None, 0
    u = unary_union(polys)
    # arcs subdivided so the chord error stays at 0.002 mm whatever r is: a fixed
    # 24 segments per quarter is 0.02 mm at r = 36 but 0.5 mm at r = 1000, and
    # there the closing fell apart
    q = 24 if r <= 0 else int(min(2000, max(24, math.ceil(
        (math.pi / 2) / (2 * math.acos(max(-1.0, 1.0 - 0.002 / r)))))))
    c = u.buffer(r, quad_segs=q).buffer(-r, quad_segs=q) if r > 0 else u
    if c.is_empty:
        c = u
    parts = list(c.geoms) if isinstance(c, MultiPolygon) else [c]
    parts.sort(key=lambda g: -g.area)
    return orient(Polygon(parts[0].exterior), 1.0), len(parts)


def resample(poly, n):
    ring = np.asarray(poly.exterior.coords)[:-1]
    cx, cy = poly.centroid.x, poly.centroid.y
    # start: where the ray from the centroid along +u leaves the contour
    seg = np.vstack([ring, ring[:1]])
    best = None
    for i in range(len(ring)):
        (x0, y0), (x1, y1) = seg[i], seg[i + 1]
        if (y0 - cy) * (y1 - cy) <= 0 and y0 != y1:
            t = (cy - y0) / (y1 - y0)
            x = x0 + t * (x1 - x0)
            if x > cx and (best is None or x > best[0]):
                best = (x, i, t)
    if best is None:
        k = int(np.argmax(ring[:, 0]))
        ring = np.roll(ring, -k, axis=0)
    else:
        _x, i, t = best
        start = seg[i] + t * (seg[i + 1] - seg[i])
        ring = np.vstack([start, ring[i + 1:], ring[:i + 1]])
    closed = np.vstack([ring, ring[:1]])
    L = np.r_[0, np.cumsum(np.linalg.norm(np.diff(closed, axis=0), axis=1))]
    s = np.linspace(0, L[-1], n, endpoint=False)
    return np.c_[np.interp(s, L, closed[:, 0]), np.interp(s, L, closed[:, 1])], L[-1]


def main():
    job = json.load(open(sys.argv[1], encoding='utf-8'))
    t0 = time.time()
    r = float(job['r_mm'])
    secs = job['sections']
    heights = [s['h'] for s in secs]
    loops = [list(s['loops']) for s in secs]
    for m in job.get('meshes') or []:
        for k, lp in enumerate(slice_mesh(m, heights)):
            loops[k].extend(lp)
    closed, warnings, multi = [], [], []
    for k, (h, lp) in enumerate(zip(heights, loops)):
        poly, islands = close_section(lp, r)
        if poly is None:
            warnings.append('plane %d (%.2f mm) cuts nothing - skipped' % (k + 1, h))
            continue
        if islands > 1:
            multi.append(k)
        closed.append((h, poly, lp))
    r_needed, bridged = None, []
    if multi:
        # Closing joins two objects of radius a, centres D apart, only when
        # D/2 <= sqrt(a^2 + 2 a r): for thin pins that takes far more than half
        # the gap. Objects that end at different heights leave planes that cut
        # only some of them, far apart - on the Dlab_155 margins a global radius
        # of 31-41 mm would have been needed. Lofting through whichever island is
        # largest changes object from plane to plane and the loft jumps.
        def r_single(lps, start):
            def single(rr):
                return close_section(lps, rr)[1] <= 1
            lo, hi = start, max(2.0 * start, 1.0)
            while not single(hi) and hi < 1000.0:
                lo, hi = hi, 2.0 * hi
            for _ in range(30):
                mid = 0.5 * (lo + hi)
                lo, hi = (lo, mid) if single(mid) else (mid, hi)
            return hi
        if job.get('bridge', True):
            # each plane that falls apart gets its own radius - the least that
            # joins it, plus 15 % so the bridge is not a pinch - and only it
            for k in multi:
                h, _p, lp = closed[k]
                rk = 1.15 * r_single(lp, r)
                poly, _n = close_section(lp, rk)
                closed[k] = (h, poly, lp)
                bridged.append((h, rk))
            warnings.append('%d plane(s) fell apart at r = %.2f mm and were bridged with '
                            'their own radius: %s'
                            % (len(bridged), r, ', '.join('%.1f mm at h %.1f' % (rk, h)
                                                          for h, rk in bridged)))
        else:
            r_needed = max(r_single(lp, r) for _h, _p, lp in (closed[k] for k in multi))
            warnings.append('%d plane(s) split into separate islands at r = %.2f mm; one '
                            'contour in every plane needs r >= %.2f mm'
                            % (len(multi), r, r_needed))
    # end planes that only graze a tilted cap give slivers the loft cannot use
    if closed:
        med = sorted(p.area for _h, p, _l in closed)[len(closed) // 2]
        keep = [c for c in closed if c[1].area >= 0.05 * med]
        if len(keep) < len(closed):
            warnings.append('%d end plane(s) dropped: they only graze the objects '
                            '(under 5 %% of the median section area)'
                            % (len(closed) - len(keep)))
        closed = keep
    closed = [(h, p) for h, p, _l in closed]
    if len(closed) < 2:
        raise ValueError('fewer than two planes cut the objects')
    perim = max(p.exterior.length for _h, p in closed)
    n = int(max(64, min(int(job.get('max_points', 160)),
                        math.ceil(perim / float(job.get('spacing_mm', 1.0))))))
    out, prev = [], None
    for h, p in closed:
        pts, L = resample(p, n)
        if prev is not None:
            # each section starts where it best matches the one below: the cyclic
            # shift with the least total distance between matching points. The
            # ray-from-centroid start alone jumped between sections wherever the
            # shape changed, and the loft creased along the jump.
            cost = [((np.roll(pts, -k, axis=0) - prev) ** 2).sum() for k in range(n)]
            pts = np.roll(pts, -int(np.argmin(cost)), axis=0)
        prev = pts
        out.append({'h': h, 'uv': pts.tolist(), 'area_mm2': p.area, 'perimeter_mm': L})
    json.dump({'ok': True, 'sections': out, 'points': n, 'warnings': warnings,
               'r_needed_mm': r_needed,
               'bridged': [[h, rk] for h, rk in bridged],
               'secs': round(time.time() - t0, 2)}, open(sys.argv[2], 'w', encoding='utf-8'))


if __name__ == '__main__':
    try:
        main()
    except Exception:
        import traceback
        json.dump({'ok': False, 'error': traceback.format_exc()},
                  open(sys.argv[2], 'w', encoding='utf-8'))
