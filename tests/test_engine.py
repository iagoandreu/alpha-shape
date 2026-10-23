"""Offline tests for the alpha-shape engine; every answer is known in advance.

    py -3.10 tests/test_engine.py
"""
import math
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import as_engine as E                                       # noqa: E402

FAILED = []


def check(name, cond, detail=''):
    print('  %s %s%s' % ('ok  ' if cond else 'FAIL', name, '' if cond else '   %s' % (detail,)))
    if not cond:
        FAILED.append(name)


def cylinder(n, r, L, axis, origin):
    ax = np.array(axis, float) / np.linalg.norm(axis)
    h = np.array([1, 0, 0.]) if abs(ax[0]) < 0.9 else np.array([0, 1, 0.])
    e1 = np.cross(ax, h); e1 /= np.linalg.norm(e1); e2 = np.cross(ax, e1)
    V, T = [], []
    for z in (0, L):
        for k in range(n):
            a = 2 * math.pi * k / n
            V.append(origin + z * ax + r * (math.cos(a) * e1 + math.sin(a) * e2))
    for k in range(n):
        j = (k + 1) % n
        T += [(k, j, n + j), (k, n + j, n + k)]
    for k in range(1, n - 1):
        T += [(0, k + 1, k), (n, n + k, n + k + 1)]
    return np.array(V), np.array(T)


cache_dir = tempfile.mkdtemp(prefix='alpha_test_')
log = lambda m: None

print('one cylinder')
v, t = cylinder(24, 1.2, 16.0, (0, 0, 1), np.zeros(3))
c = E.prepare({'coords': v.ravel().tolist(), 'idx': t.ravel().tolist(),
               'resolution_mm': 0.4, 'cache_dir': cache_dir}, log)
poly = 0.5 * 24 * 1.2 ** 2 * math.sin(2 * math.pi / 24) * 16.0
r = E.shape(c, 3.0, log)
check('alpha above the radius gives one closed solid', r['watertight'] and r['pieces'] == 1, r)
check('its volume is the cylinder\'s within 3 %', abs(r['volume_mm3'] / poly - 1) < 0.03,
      (r['volume_mm3'], poly))
r = E.shape(c, 1000.0, log)
check('a huge alpha is the convex hull (the cylinder itself)', abs(r['volume_mm3'] / poly - 1) < 0.01,
      (r['volume_mm3'], poly))
r = E.shape(c, 0.3, log)
check('alpha below the radius loses the solid', r['volume_mm3'] < 0.2 * poly, r['volume_mm3'])
c2 = E.prepare({'coords': v.ravel().tolist(), 'idx': t.ravel().tolist(),
                'resolution_mm': 0.4, 'cache_dir': cache_dir}, log)
check('the same input and resolution reuse the cache', c2 == c)

print('two cylinders 6 mm apart')
v1, t1 = cylinder(24, 1.2, 16.0, (0, 0, 1), np.zeros(3))
v2, t2 = cylinder(24, 1.2, 16.0, (0, 0, 1), np.array([6.0, 0, 0]))
v = np.vstack([v1, v2]); t = np.vstack([t1, t2 + len(v1)])
c = E.prepare({'coords': v.ravel().tolist(), 'idx': t.ravel().tolist(),
               'resolution_mm': 0.4, 'cache_dir': cache_dir}, log)
r = E.shape(c, 1.5, log)
check('alpha smaller than the gap keeps them apart', r['pieces'] == 2, r['pieces'])
r = E.shape(c, 8.0, log)
check('alpha larger than the gap bridges them into one', r['pieces'] == 1 and r['watertight'], r)
# hull of the pair: one pin's section (two half 24-gons) plus the 6 x 2.4 web
hull = poly + 6.0 * 2.4 * 16.0
check('bridged, it holds the web: more than the two pins, at most the hull',
      2 * poly * 1.2 < r['volume_mm3'] <= hull * 1.005, (r['volume_mm3'], 2 * poly, hull))
r = E.shape(c, 1000.0, log)
check('a huge alpha is the hull of the pair (within 2 %)',
      abs(r['volume_mm3'] / hull - 1) < 0.02, (r['volume_mm3'], hull))

print('repair')
P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1.]])
T = np.array([[0, 1, 2, 3], [0, 4, 5, 6]])
check('two tetrahedra touching at a point are a pinch',
      list(E._pinched_vertices(E.Topology(P, T).faces[E.Topology(P, T).boundary(np.array([True, True]))])) == [0])

print()
if FAILED:
    print('%d FAILED: %s' % (len(FAILED), ', '.join(FAILED)))
    raise SystemExit(1)
print('all engine tests passed')
