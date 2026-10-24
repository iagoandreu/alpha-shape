"""Offline tests for the section engine; every answer is known in advance.

    py -3.10 tests/test_engine.py
"""
import json
import math
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.path.join(os.path.dirname(HERE), 'as_engine.py')
sys.path.insert(0, os.path.dirname(HERE))
import as_engine as E                                       # noqa: E402

FAILED = []


def check(name, cond, detail=''):
    print('  %s %s%s' % ('ok  ' if cond else 'FAIL', name, '' if cond else '   %s' % (detail,)))
    if not cond:
        FAILED.append(name)


def circle(cx, cy, r, n=96):
    return [[cx + r * math.cos(2 * math.pi * k / n), cy + r * math.sin(2 * math.pi * k / n)]
            for k in range(n)]


def run(job):
    d = tempfile.mkdtemp()
    jp, op = os.path.join(d, 'j.json'), os.path.join(d, 'o.json')
    json.dump(job, open(jp, 'w'))
    subprocess.run([sys.executable, ENGINE, jp, op], check=True)
    return json.load(open(op))


print('closing one section')
poly, n = E.close_section([circle(0, 0, 1.2)], 3.0)
check('one circle closes to itself (a 96-gon: within 0.5 %)',
      n == 1 and abs(poly.area / (math.pi * 1.44) - 1) < 0.005, poly.area)
poly, n = E.close_section([circle(0, 0, 1.2), circle(6, 0, 1.2)], 1.0)
check('two circles 6 mm apart stay apart at r = 1', n == 2, n)
poly, n = E.close_section([circle(0, 0, 1.2), circle(6, 0, 1.2)], 4.0)
stadium = math.pi * 1.44 + 6 * 2.4
check('at r = 4 they join, with a waist: more than two circles, less than the hull',
      n == 1 and 2 * math.pi * 1.44 < poly.area < stadium, (n, poly.area, stadium))
areas = [E.close_section([circle(0, 0, 1.2), circle(6, 0, 1.2)], rr)[0].area
         for rr in (5.0, 20.0, 100.0)]
check('a larger r fills more, towards the hull', areas[0] < areas[1] < areas[2] < stadium,
      (areas, stadium))
# at r = 100 the waist is an arc of radius 100 over the 6 mm span: its sag costs
# about 2 x (2/3) x 6 x 6^2 / (8 x 100) = 0.36 mm2 of the hull
check('and at r = 100 it is the hull less the waist sag (within 0.15 mm2)',
      abs((stadium - areas[2]) - 0.36) < 0.15, stadium - areas[2])

print('resampling and correspondence')
pts, L = E.resample(E.close_section([circle(0, 0, 5)], 0)[0], 64)
check('64 points, evenly spaced on the contour',
      len(pts) == 64 and abs(L - 2 * math.pi * 5) < 0.05, (len(pts), L))
check('the first point is on +u from the centroid', pts[0][0] > 4.9 and abs(pts[0][1]) < 0.1,
      pts[0])
res = run({'r_mm': 0.0, 'spacing_mm': 0.75,
           'sections': [{'h': h, 'loops': [circle(0.3 * h, 0, 5)]} for h in range(6)]})
first = [s['uv'][0] for s in res['sections']]
check('matching points follow a drifting section without twisting',
      all(abs(p[1]) < 0.3 for p in first), first)

print('planes that fall apart')
secs = [{'h': 0, 'loops': [circle(0, 0, 1.2), circle(4, 0, 1.2)]},
        {'h': 5, 'loops': [circle(0, 0, 1.2), circle(20, 0, 1.2)]}]
res = run({'r_mm': 2.0, 'bridge': False, 'sections': secs})
check('without bridging, the radius that joins every plane is reported',
      res['r_needed_mm'] and res['r_needed_mm'] > 2.0, res.get('r_needed_mm'))
res = run({'r_mm': 2.0, 'bridge': True, 'sections': secs})
check('with bridging, only the plane that fell apart gets its own radius',
      len(res['bridged']) == 1 and res['bridged'][0][0] == 5, res['bridged'])

print()
if FAILED:
    print('%d FAILED: %s' % (len(FAILED), ', '.join(FAILED)))
    raise SystemExit(1)
print('all engine tests passed')
