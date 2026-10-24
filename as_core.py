"""Alpha Shape - the work. Re-imported on every command event and dev job.

Wraps the selected bodies in a LOFTED SOLID, not a mesh:
  1. a direction - picked (planar face, plane, axis, edge, cylinder) or, by
     default, the mean axis of the cylinders among the bodies;
  2. N planes across the bodies' extent along it (minus optional end offsets);
  3. each plane cuts every body exactly (TemporaryBRepManager.planeIntersection;
     mesh bodies are sliced in the engine);
  4. in each plane the cuts are united and closed with a disk of radius r
     (dilate r, erode r): gaps narrower than 2r fill with round arcs of radius
     r, and the outer contour is kept - as_engine.py, numpy + shapely;
  5. the contours, resampled with matching points, become closed splines and a
     solid loft joins them.

The preview draws the contours. All geometry crossing into Fusion is in cm; the
engine works in mm in each plane's own (u, v).
"""

import adsk.core
import adsk.fusion

import hashlib
import json
import math
import os
import subprocess
import sys
import tempfile
import time
import traceback
import types

ADDIN_DIR = os.path.dirname(os.path.realpath(__file__))
ENGINE = os.path.join(ADDIN_DIR, 'as_engine.py')
WORK_DIR = os.path.join(tempfile.gettempdir(), 'alpha_shape_work')
SETTINGS_PATH = os.path.join(ADDIN_DIR, 'last_settings.json')
LOG_PATH = os.path.join(ADDIN_DIR, 'AlphaShape.log')
RES_DIR = os.path.join(ADDIN_DIR, 'resources', 'AlphaShape')
CM = 10.0
MM = 0.1
NO_WINDOW = 0x08000000 if os.name == 'nt' else 0

CMD_ID = 'AlphaShape_Create'
CMD_NAME = 'Alpha Shape'
PANELS = ['SolidCreatePanel', 'ParaMeshCreatePanel']

RANGE_SHARED = 'Where every object is cut'
RANGE_ALL = 'Whole extent of the objects'
RANGES = [RANGE_SHARED, RANGE_ALL]

DEFAULTS = {'r_mm': 3.0, 'planes': 12, 'bottom_mm': 0.0, 'top_mm': 0.0,
            'range': RANGE_SHARED, 'join': False,
            'bridge': True, 'preview': True, 'hide_inputs': False}
CONTOUR_COLOUR = (40, 140, 230)
END_EPS_MM = 0.02          # a plane exactly on an end cap cuts a face, not a solid


def log(msg):
    try:
        with open(LOG_PATH, 'a', encoding='utf-8') as fh:
            fh.write('%s  %s\n' % (time.strftime('%H:%M:%S'), msg))
    except Exception:
        pass


def short_error(tb):
    lines = [l for l in tb.strip().splitlines() if l.strip()]
    return lines[-1] if lines else tb


def state():
    st = sys.modules.get('_as_state')
    if st is None:
        st = types.ModuleType('_as_state')
        st.python = None
        st.cuts = {}
        st.last = None
        st.report = ''
        sys.modules['_as_state'] = st
    return st


def load_settings():
    o = dict(DEFAULTS)
    try:
        with open(SETTINGS_PATH, encoding='utf-8') as fh:
            o.update({k: v for k, v in json.load(fh).items() if k in DEFAULTS})
    except Exception:
        pass
    return o


def save_settings(o):
    try:
        with open(SETTINGS_PATH, 'w', encoding='utf-8') as fh:
            json.dump({k: o[k] for k in DEFAULTS if k in o}, fh, indent=1)
    except Exception:
        pass


def is_parametric(design):
    return design.designType == adsk.fusion.DesignTypes.ParametricDesignType


# --------------------------------------------------------------------------- #
# small vector helpers (Fusion's Python has no numpy)
# --------------------------------------------------------------------------- #

def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _unit(a):
    n = math.sqrt(_dot(a, a))
    return (a[0] / n, a[1] / n, a[2] / n) if n > 0 else a


def _tup(p):
    return (p.x, p.y, p.z)


# --------------------------------------------------------------------------- #
# the engine
# --------------------------------------------------------------------------- #

def python_cmd():
    """An interpreter with numpy and shapely: ALPHASHAPE_PYTHON, python_path.txt
    beside the add-in, the usual install folders, then py -3.10."""
    st = state()
    if st.python:
        return st.python
    cands = []
    if os.environ.get('ALPHASHAPE_PYTHON'):
        cands.append([os.environ['ALPHASHAPE_PYTHON']])
    try:
        with open(os.path.join(ADDIN_DIR, 'python_path.txt'), encoding='utf-8') as fh:
            line = fh.read().strip()
            if line:
                cands.append([line])
    except Exception:
        pass
    local = os.environ.get('LOCALAPPDATA', '')
    for v in ('310', '311', '312', '313'):
        cands.append([os.path.join(local, 'Programs', 'Python', 'Python%s' % v, 'python.exe')])
    cands.append(['py', '-3.10'])
    for c in cands:
        try:
            r = subprocess.run(c + ['-c', 'import numpy, shapely'], capture_output=True,
                               timeout=60, creationflags=NO_WINDOW)
            if r.returncode == 0:
                st.python = c
                return c
        except Exception:
            pass
    raise RuntimeError('No Python with numpy and shapely was found (pip install shapely). '
                       'Set ALPHASHAPE_PYTHON or put the interpreter path in '
                       'python_path.txt beside the add-in.')


def run_engine(job, timeout=600):
    os.makedirs(WORK_DIR, exist_ok=True)
    jp = os.path.join(WORK_DIR, 'job.json')
    op = os.path.join(WORK_DIR, 'out.json')
    with open(jp, 'w', encoding='utf-8') as fh:
        json.dump(job, fh)
    if os.path.exists(op):
        os.remove(op)
    r = subprocess.run(python_cmd() + [ENGINE, jp, op], capture_output=True, text=True,
                       timeout=timeout, creationflags=NO_WINDOW)
    if not os.path.exists(op):
        raise RuntimeError('the section engine produced nothing:\n%s' % r.stderr[-800:])
    with open(op, encoding='utf-8') as fh:
        res = json.load(fh)
    if not res.get('ok'):
        err = res.get('error', '')
        raise RuntimeError(short_error(err) if 'ValueError' in err else
                           'the section engine failed:\n%s' % err[-800:])
    return res


# --------------------------------------------------------------------------- #
# direction and range
# --------------------------------------------------------------------------- #

def selected(sel):
    return [sel.selection(i).entity for i in range(sel.selectionCount)]


def direction_of(entity):
    """A unit direction from a picked plane, face, axis or edge; None if unusable."""
    for cast, get in (
            (adsk.fusion.ConstructionPlane, lambda e: e.geometry.normal),
            (adsk.fusion.ConstructionAxis, lambda e: e.geometry.direction),
            (adsk.fusion.SketchLine, lambda e: _sub(_tup(e.worldGeometry.endPoint),
                                                    _tup(e.worldGeometry.startPoint)))):
        e = cast.cast(entity)
        if e:
            g = get(e)
            return _unit(g if isinstance(g, tuple) else _tup(g))
    f = adsk.fusion.BRepFace.cast(entity)
    if f:
        g = f.geometry
        pl = adsk.core.Plane.cast(g)
        if pl:
            return _unit(_tup(pl.normal))
        cy = adsk.core.Cylinder.cast(g)
        if cy:
            return _unit(_tup(cy.axis))
        cn = adsk.core.Cone.cast(g)
        if cn:
            return _unit(_tup(cn.axis))
    e = adsk.fusion.BRepEdge.cast(entity)
    if e:
        ln = adsk.core.Line3D.cast(e.geometry)
        if ln:
            return _unit(_sub(_tup(ln.endPoint), _tup(ln.startPoint)))
    return None


def body_points(b):
    """A cloud of points on a body, cm: the nodes of a coarse (0.1 mm) tessellation."""
    mb = adsk.fusion.MeshBody.cast(b)
    if mb:
        c = mb.displayMesh.nodeCoordinatesAsDouble
    else:
        # not displayMeshes.bestMesh: a body that has not been drawn yet (hidden,
        # or its document just activated) has none and it raises "count > 0"
        calc = adsk.fusion.BRepBody.cast(b).meshManager.createMeshCalculator()
        calc.surfaceTolerance = 0.01
        c = calc.calculate().nodeCoordinatesAsDouble
    return [(c[i], c[i + 1], c[i + 2]) for i in range(0, len(c), 3)]


def auto_direction(bodies, pts):
    """Mean axis of the cylindrical faces, weighted by area and sign-aligned; with
    none, the principal axis of the points."""
    acc, ref = [0.0, 0.0, 0.0], None
    for b in bodies:
        bb = adsk.fusion.BRepBody.cast(b)
        if not bb:
            continue
        for f in bb.faces:
            cy = adsk.core.Cylinder.cast(f.geometry)
            if not cy:
                continue
            a = _unit(_tup(cy.axis))
            if ref is None:
                ref = a
            s = 1.0 if _dot(a, ref) >= 0 else -1.0
            for k in range(3):
                acc[k] += s * f.area * a[k]
    if ref is not None and _dot(acc, acc) > 0:
        return _unit(tuple(acc)), 'mean axis of the cylinders'
    # principal axis by power iteration on the covariance
    n = float(len(pts))
    c = tuple(sum(p[k] for p in pts) / n for k in range(3))
    C = [[sum((p[i] - c[i]) * (p[j] - c[j]) for p in pts) / n for j in range(3)] for i in range(3)]
    v = (1.0, 0.7, 0.3)
    for _ in range(100):
        v = _unit(tuple(sum(C[i][j] * v[j] for j in range(3)) for i in range(3)))
    return v, 'principal axis of the bodies'


def frame(d):
    h = (1.0, 0.0, 0.0) if abs(d[0]) < 0.9 else (0.0, 1.0, 0.0)
    e1 = _unit(_cross(d, h))
    return e1, _cross(d, e1)


# --------------------------------------------------------------------------- #
# the cuts
# --------------------------------------------------------------------------- #

def _chain(polylines, tol=1e-4):
    """Join open polylines end to end into closed loops (cm)."""
    loops, open_ = [], [list(p) for p in polylines if len(p) >= 2]
    d = lambda a, b: math.sqrt(_dot(_sub(a, b), _sub(a, b)))
    while open_:
        cur = open_.pop(0)
        grown = True
        while grown and d(cur[0], cur[-1]) > tol:
            grown = False
            for i, p in enumerate(open_):
                for rev in (False, True):
                    q = p[::-1] if rev else p
                    if d(cur[-1], q[0]) <= tol:
                        cur.extend(q[1:])
                        open_.pop(i)
                        grown = True
                        break
                if grown:
                    break
        if d(cur[0], cur[-1]) <= tol and len(cur) >= 4:
            loops.append(cur[:-1])
    return loops


def plane_cuts(body, origin, d, e1, e2, heights_mm, tol_cm=0.001):
    """Closed cut loops of a BRep body at each height, as (u, v) mm.

    Only planes inside the body's own extent are asked for: a plane that misses
    RAISES (InternalValidationError : result), and the error then stays pending
    and surfaces on the next API call - TemporaryBRepManager.get() included."""
    tbm = adsk.fusion.TemporaryBRepManager.get()
    hs = [_dot(_sub(p, origin), d) * CM for p in body_points(body)]
    lo, hi = min(hs) + END_EPS_MM, max(hs) - END_EPS_MM
    out = []
    for h in heights_mm:
        if not lo <= h <= hi:
            out.append([])
            continue
        o = tuple(origin[k] + h * MM * d[k] for k in range(3))
        wire = tbm.planeIntersection(body, adsk.core.Plane.create(
            adsk.core.Point3D.create(*o), adsk.core.Vector3D.create(*d)))
        polys = []
        if wire:
            for e in wire.edges:
                ev = e.evaluator
                ok, p0, p1 = ev.getParameterExtents()
                ok, pts = ev.getStrokes(p0, p1, tol_cm)
                polys.append([_tup(p) for p in pts])
        loops = []
        for lp in _chain(polys):
            loops.append([[_dot(_sub(p, origin), e1) * CM, _dot(_sub(p, origin), e2) * CM]
                          for p in lp])
        out.append(loops)
    return out


def setup(bodies, dir_entity, o):
    """Direction, frame, range and the planes' heights (mm from the origin)."""
    pts = [p for b in bodies for p in body_points(b)]
    d, why = None, ''
    if dir_entity is not None:
        d = direction_of(dir_entity)
        why = 'picked'
        if d is None:
            raise ValueError('that selection gives no direction; pick a planar face, '
                             'a plane, an axis, a straight edge or a cylinder')
    if d is None:
        d, why = auto_direction(bodies, pts)
    e1, e2 = frame(d)
    n = float(len(pts))
    origin = tuple(sum(p[k] for p in pts) / n for k in range(3))
    # Where every object is cut, by default. Over the whole extent, objects that
    # end at different heights leave planes cutting only some of them, far
    # apart: on the Dlab_155 margins the AK pins run 13-20 mm past the screw
    # margins, those planes needed 36 mm bridges and the loft creased, or failed
    # outright with 32 planes. Over the shared range the bridges were 3.5-6 mm and
    # 12 and 24 planes gave the same volume within 0.5 %.
    ext = []
    for b in bodies:
        hb = [_dot(_sub(p, origin), d) * CM for p in body_points(b)]
        ext.append((min(hb), max(hb)))
    if o.get('range', RANGE_SHARED) == RANGE_SHARED:
        lo0, hi0 = max(e[0] for e in ext), min(e[1] for e in ext)
        if hi0 <= lo0:
            raise ValueError('no stretch of the direction cuts every object; choose '
                             '"%s"' % RANGE_ALL)
    else:
        lo0, hi0 = min(e[0] for e in ext), max(e[1] for e in ext)
    lo = lo0 + float(o['bottom_mm']) + END_EPS_MM
    hi = hi0 - float(o['top_mm']) - END_EPS_MM
    if hi <= lo:
        raise ValueError('the end offsets leave no length between the planes')
    N = max(2, int(o['planes']))
    heights = [lo + (hi - lo) * k / (N - 1) for k in range(N)]
    return {'d': d, 'why': why, 'e1': e1, 'e2': e2, 'origin': origin,
            'heights': heights, 'length_mm': hi - lo,
            'outside_mm': max(0.0, lo0 - min(e[0] for e in ext)) +
            max(0.0, max(e[1] for e in ext) - hi0)}


def sections(bodies, dir_entity, o):
    st = state()
    toks = sorted(getattr(b, 'entityToken', b.name) for b in bodies)
    dtok = getattr(dir_entity, 'entityToken', '') if dir_entity is not None else ''
    fr = setup(bodies, dir_entity, o)
    key = hashlib.sha1(('|'.join(toks) + dtok + repr(fr['heights']) + repr(fr['d'])
                        ).encode()).hexdigest()
    if key not in st.cuts:
        loops = [[] for _ in fr['heights']]
        meshes = []
        for b in bodies:
            bb = adsk.fusion.BRepBody.cast(b)
            if bb:
                for k, lp in enumerate(plane_cuts(bb, fr['origin'], fr['d'], fr['e1'],
                                                  fr['e2'], fr['heights'])):
                    loops[k].extend(lp)
            else:
                dm = adsk.fusion.MeshBody.cast(b).displayMesh
                meshes.append({'coords': [v * CM for v in dm.nodeCoordinatesAsDouble],
                               'idx': list(dm.nodeIndices),
                               'origin': [v * CM for v in fr['origin']],
                               'e1': list(fr['e1']), 'e2': list(fr['e2']), 'd': list(fr['d'])})
        st.cuts = {key: (loops, meshes)}
    loops, meshes = st.cuts[key]
    return fr, loops, meshes


def compute(bodies, dir_entity, o):
    st = state()
    t0 = time.time()
    fr, loops, meshes = sections(bodies, dir_entity, o)
    # the engine's own date too: an edited engine must not be answered from cache
    sig = (id(loops), float(o['r_mm']), bool(o.get('bridge', True)),
           os.path.getmtime(ENGINE))
    if st.last and st.last[0] == sig:
        return fr, st.last[1]
    res = run_engine({'r_mm': float(o['r_mm']), 'spacing_mm': 0.75, 'max_points': 160,
                      'bridge': bool(o.get('bridge', True)),
                      'sections': [{'h': h, 'loops': lp} for h, lp in zip(fr['heights'], loops)],
                      'meshes': meshes})
    res['wall_s'] = round(time.time() - t0, 2)
    st.last = (sig, res)
    return fr, res


def to_world(fr, h, uv):
    o, d, e1, e2 = fr['origin'], fr['d'], fr['e1'], fr['e2']
    return [tuple(o[k] + (h * d[k] + u * e1[k] + v * e2[k]) * MM for k in range(3))
            for u, v in uv]


# --------------------------------------------------------------------------- #
# preview graphics
# --------------------------------------------------------------------------- #

GROUP_ID = 'AlphaShapePreview'


def clear_graphics(app):
    try:
        des = adsk.fusion.Design.cast(app.activeProduct)
        if not des:
            return
        cg = des.rootComponent.customGraphicsGroups
        for i in range(cg.count - 1, -1, -1):
            if cg.item(i).id == GROUP_ID:
                cg.item(i).deleteMe()
    except Exception:
        pass


def draw(app, fr, res):
    clear_graphics(app)
    secs = res.get('sections') or []
    if not secs:
        return
    des = adsk.fusion.Design.cast(app.activeProduct)
    g = des.rootComponent.customGraphicsGroups.add()
    g.id = GROUP_ID
    coords, idx = [], []
    for s in secs:
        base = len(coords) // 3
        pts = to_world(fr, s['h'], s['uv'])
        for p in pts:
            coords.extend(p)
        n = len(pts)
        for i in range(n):
            idx.extend([base + i, base + (i + 1) % n])
    # the matching points, joined along the stack: the loft's own lines
    n = len(secs[0]['uv'])
    for j in range(0, n, max(1, n // 12)):
        for k in range(len(secs) - 1):
            idx.extend([k * n + j, (k + 1) * n + j])
    lines = g.addLines(adsk.fusion.CustomGraphicsCoordinates.create(coords), idx, False)
    c = CONTOUR_COLOUR
    lines.color = adsk.fusion.CustomGraphicsSolidColorEffect.create(
        adsk.core.Color.create(c[0], c[1], c[2], 255))
    lines.weight = 2.0
    lines.depthPriority = 1
    app.activeViewport.refresh()


# --------------------------------------------------------------------------- #
# the solid
# --------------------------------------------------------------------------- #

def build_loft(app, fr, res, o):
    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.rootComponent
    sk = comp.sketches.add(comp.xYConstructionPlane)
    sk.name = 'Alpha shape sections'
    sk.isComputeDeferred = True
    splines = []
    for s in res['sections']:
        coll = adsk.core.ObjectCollection.create()
        for p in to_world(fr, s['h'], s['uv']):
            coll.add(sk.modelToSketchSpace(adsk.core.Point3D.create(*p)))
        spl = sk.sketchCurves.sketchFittedSplines.add(coll)
        spl.isClosed = True
        splines.append(spl)
    sk.isComputeDeferred = False
    sk.isVisible = False
    lin = comp.features.loftFeatures.createInput(
        adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
    lin.isSolid = True
    for spl in splines:
        lin.loftSections.add(comp.features.createPath(spl, False))
    try:
        loft = comp.features.loftFeatures.add(lin)
    except Exception:
        sk.deleteMe()
        raise RuntimeError('the loft refused these sections: %s'
                           % short_error(traceback.format_exc()))
    loft.name = 'Alpha shape loft'
    body = loft.bodies.item(0)
    body.name = 'Alpha shape r%.2f n%d' % (float(o['r_mm']), len(splines))
    return body


def join_objects(app, body, bodies):
    """Join the loft with the selected solids, keeping them: the parts outside
    the lofted range become part of one body."""
    comp = body.parentComponent
    tools = adsk.core.ObjectCollection.create()
    for b in bodies:
        bb = adsk.fusion.BRepBody.cast(b)
        if bb and bb.isSolid:
            tools.add(bb)
    if not tools.count:
        return 0
    ci = comp.features.combineFeatures.createInput(body, tools)
    ci.operation = adsk.fusion.FeatureOperations.JoinFeatureOperation
    ci.isKeepToolBodies = True
    f = comp.features.combineFeatures.add(ci)
    f.name = 'Alpha shape join'
    return tools.count


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #

def describe(fr, res, o):
    secs = res.get('sections') or []
    areas = [s['area_mm2'] for s in secs]
    rows = ['direction: %s (%.3f, %.3f, %.3f)' % ((fr['why'],) + tuple(fr['d'])),
            '%d sections over %.2f mm, %d matching points each'
            % (len(secs), fr['length_mm'], res.get('points', 0))]
    if fr.get('outside_mm', 0) > 0.05:
        rows.append('%.1f mm of the objects lies outside the lofted range%s'
                    % (fr['outside_mm'], ' (joined back in)' if o.get('join', False) else ''))
    if areas:
        rows.append('section area %.1f - %.1f mm2' % (min(areas), max(areas)))
    for w in res.get('warnings', []):
        rows.append('<b>%s</b>' % w)
    rows.append('%.1f s' % res.get('wall_s', 0))
    return '<br>'.join(rows)


# --------------------------------------------------------------------------- #
# the dialog
# --------------------------------------------------------------------------- #

def on_created(app, cmd):
    o = load_settings()
    inputs = cmd.commandInputs
    cmd.setDialogInitialSize(360, 420)
    cmd.okButtonText = 'Create'
    cmd.isExecutedWhenPreEmpted = False
    sel = inputs.addSelectionInput('bodies', 'Bodies', 'The objects to wrap')
    sel.addSelectionFilter('SolidBodies')
    sel.addSelectionFilter('SurfaceBodies')
    sel.addSelectionFilter('MeshBodies')
    sel.setSelectionLimits(1, 0)
    ds = inputs.addSelectionInput('direction', 'Direction',
                                  'Optional: a planar face or plane (its normal), an axis, '
                                  'a straight edge or a cylinder')
    for f in ('PlanarFaces', 'CylindricalFaces', 'ConstructionPlanes', 'ConstructionLines',
              'LinearEdges', 'SketchLines'):
        ds.addSelectionFilter(f)
    ds.setSelectionLimits(0, 1)
    ds.tooltip = ('The sections are perpendicular to this. Leave it empty to use the mean '
                  'axis of the cylinders among the bodies.')
    # capped: above ~300 mm GEOS loses precision on millimetre objects
    v = inputs.addFloatSpinnerCommandInput('r', 'Closing radius', 'mm', 0.0, 20.0, 0.5,
                                           o['r_mm'] * MM)
    v.tooltip = ('In every section, gaps between the objects narrower than twice this are '
                 'filled, with round arcs of this radius. 0 just unites the cuts.')
    inputs.addIntegerSpinnerCommandInput('planes', 'Planes', 2, 400, 1, int(o['planes']))
    dd = inputs.addDropDownCommandInput('range', 'Range',
                                        adsk.core.DropDownStyles.TextListDropDownStyle)
    for label in RANGES:
        dd.listItems.add(label, label == o.get('range', RANGE_SHARED))
    dd.tooltip = ('Where the planes go. "%s": only the stretch every object crosses, '
                  'which keeps the sections alike; the parts that stick out stay outside '
                  'the loft unless "Join" is on. "%s": end to end.'
                  % (RANGE_SHARED, RANGE_ALL))
    inputs.addFloatSpinnerCommandInput('bottom', 'Bottom offset', 'mm', 0.0, 1000.0, 0.5,
                                       o['bottom_mm'] * MM)
    inputs.addFloatSpinnerCommandInput('top', 'Top offset', 'mm', 0.0, 1000.0, 0.5,
                                       o['top_mm'] * MM)
    v = inputs.addBoolValueInput('join', 'Join the objects to the loft', True, '',
                                 bool(o.get('join', False)))
    v.tooltip = ('Join the loft with the selected solids (which are kept), so the parts '
                 'outside the lofted range are part of one body.')
    v = inputs.addBoolValueInput('bridge', 'Bridge planes that fall apart', True, '',
                                 bool(o.get('bridge', True)))
    v.tooltip = ('Where objects end at different heights, a plane can cut only some of '
                 'them, too far apart for the closing radius to join. On: that plane alone '
                 'gets the radius that joins it (the report lists them). Off: Create '
                 'refuses and says which radius would join every plane.')
    inputs.addBoolValueInput('preview', 'Preview', True, '', bool(o['preview']))
    inputs.addBoolValueInput('hide_inputs', 'Hide the selected bodies afterwards', True, '',
                             bool(o['hide_inputs']))
    inputs.addTextBoxCommandInput('info', '', state().report or
                                  'Select the bodies to wrap.', 7, True)


def read_options(inputs):
    return {'r_mm': inputs.itemById('r').value * CM,
            'planes': inputs.itemById('planes').value,
            'bottom_mm': inputs.itemById('bottom').value * CM,
            'top_mm': inputs.itemById('top').value * CM,
            'bridge': inputs.itemById('bridge').value,
            'range': inputs.itemById('range').selectedItem.name,
            'join': inputs.itemById('join').value,
            'preview': inputs.itemById('preview').value,
            'hide_inputs': inputs.itemById('hide_inputs').value}


def _dir_entity(inputs):
    ds = inputs.itemById('direction')
    return ds.selection(0).entity if ds.selectionCount else None


def on_input_changed(app, inputs, changed_id):
    if changed_id == 'preview' and not inputs.itemById('preview').value:
        clear_graphics(app)


def on_validate(inputs):
    return (inputs.itemById('bodies').selectionCount > 0
            and inputs.itemById('planes').value >= 2)


def on_preview(app, inputs, args):
    if not inputs.itemById('preview').value:
        clear_graphics(app)
        return
    o = read_options(inputs)
    try:
        fr, res = compute(selected(inputs.itemById('bodies')), _dir_entity(inputs), o)
        text = describe(fr, res, o)
    except Exception:
        log(traceback.format_exc())
        fr, res, text = None, {}, 'Preview failed: %s' % short_error(traceback.format_exc())
    # any write to an input wipes custom graphics: text first, contours last
    info = inputs.itemById('info')
    if info.formattedText != text:
        info.formattedText = text
    if fr:
        draw(app, fr, res)
    else:
        clear_graphics(app)


def on_execute(app, ui, inputs):
    o = read_options(inputs)
    save_settings(o)
    bodies = selected(inputs.itemById('bodies'))
    clear_graphics(app)
    try:
        fr, res = compute(bodies, _dir_entity(inputs), o)
        if res.get('r_needed_mm'):
            raise ValueError('With a closing radius of %.2f mm some planes fall apart into '
                             'separate islands, so there is no single envelope to loft. '
                             'Raise the closing radius to at least %.2f mm.'
                             % (o['r_mm'], res['r_needed_mm']))
        body = build_loft(app, fr, res, o)
        if o.get('join', False):
            join_objects(app, body, bodies)
    except (ValueError, RuntimeError) as exc:
        ui.messageBox(str(exc), CMD_NAME)
        return
    if o['hide_inputs']:
        for b in bodies:
            try:
                b.isLightBulbOn = False
            except Exception:
                pass
    state().report = ('Last: "%s", %.1f mm3<br>' % (body.name, body.volume * 1000.0)
                      + describe(fr, res, o))


def on_destroy(app):
    clear_graphics(app)


# --------------------------------------------------------------------------- #
# toolbar
# --------------------------------------------------------------------------- #

TOOLTIP = ('Wrap the selected bodies in a lofted solid.\n\n'
           'Planes perpendicular to a direction (by default the mean axis of the '
           'cylinders) cut every body; in each plane the cuts are joined, gaps '
           'narrower than twice the closing radius are filled with round arcs, and '
           'the outer contours become the sections of a loft.')


def install(app, ui, created_handler_factory):
    old = ui.commandDefinitions.itemById(CMD_ID)
    if old:
        old.deleteMe()
    folder = RES_DIR if os.path.isdir(RES_DIR) and os.listdir(RES_DIR) else ''
    cd = ui.commandDefinitions.addButtonDefinition(CMD_ID, CMD_NAME, TOOLTIP, folder)
    cd.commandCreated.add(created_handler_factory())
    placed = []
    for pid in PANELS:
        panel = ui.allToolbarPanels.itemById(pid)
        if panel and not panel.controls.itemById(CMD_ID):
            panel.controls.addCommand(cd)
            placed.append(pid)
    return placed


def uninstall(ui):
    for pid in PANELS:
        panel = ui.allToolbarPanels.itemById(pid)
        if panel:
            c = panel.controls.itemById(CMD_ID)
            if c:
                c.deleteMe()
    d = ui.commandDefinitions.itemById(CMD_ID)
    if d:
        d.deleteMe()


# --------------------------------------------------------------------------- #
# dev bridge
# --------------------------------------------------------------------------- #

def _design(app):
    return adsk.fusion.Design.cast(app.activeProduct)


def _bodies_named(app, names):
    comp = _design(app).rootComponent
    out = []
    for n in names:
        b = comp.bRepBodies.itemByName(n)
        if b is None:                     # MeshBodies has no itemByName
            b = next((m for m in comp.meshBodies if m.name == n), None)
        if b is None:
            raise ValueError('no body named %r' % n)
        out.append(b)
    return out


def dev_ping(app, job):
    return {'doc': app.activeDocument.name if app.activeDocument else None,
            'engine': python_cmd()}


def dev_run(app, job):
    """Sections (and with create=True, the loft) for named bodies."""
    o = dict(DEFAULTS)
    o.update(job.get('options') or {})
    bodies = _bodies_named(app, job['bodies'])
    fr, res = compute(bodies, None, o)
    out = {'direction': fr['d'], 'why': fr['why'], 'length_mm': fr['length_mm'],
           'sections': len(res['sections']), 'points': res['points'],
           'areas': [round(s['area_mm2'], 2) for s in res['sections']],
           'warnings': res['warnings'], 'wall_s': res['wall_s']}
    if job.get('draw'):
        draw(app, fr, res)
    out['r_needed_mm'] = res.get('r_needed_mm')
    out['bridged'] = res.get('bridged')
    if job.get('create') and not res.get('r_needed_mm'):
        body = build_loft(app, fr, res, o)
        if o.get('join', False):
            out['joined'] = join_objects(app, body, bodies)
        out.update({'body': body.name, 'solid': body.isSolid, 'faces': body.faces.count,
                    'volume_mm3': round(body.volume * 1000.0, 2)})
    return out


def dev_clear(app, job):
    clear_graphics(app)
    return {'cleared': True}


def dev_snapshot(app, job):
    """Save the current view as it is - the camera is not touched."""
    vp = app.activeViewport
    if job.get('fit'):                    # only in a throwaway document
        cam = vp.camera
        cam.isSmoothTransition = False
        cam.viewOrientation = adsk.core.ViewOrientations.IsoTopRightViewOrientation
        cam.isFitView = True
        vp.camera = cam
        adsk.doEvents()
        vp.refresh()
    ok = vp.saveAsImageFile(job['path'], int(job.get('width', 1400)),
                            int(job.get('height', 1000)))
    return {'saved': bool(ok), 'path': job['path']}


def dev_make_test(app, job):
    """A new document with three pins, so a build never touches the user's."""
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    des = _design(app)
    if not is_parametric(des):
        des.designType = adsk.fusion.DesignTypes.ParametricDesignType
    comp = des.rootComponent
    tbm = adsk.fusion.TemporaryBRepManager.get()
    base = comp.features.baseFeatures.add()
    base.startEdit()
    for x, y, tilt in job.get('pins', [(0, 0, 0), (0.8, 0.2, 0.15), (0.3, 0.9, -0.1)]):
        cyl = tbm.createCylinderOrCone(adsk.core.Point3D.create(x, y, 0), 0.12,
                                       adsk.core.Point3D.create(x + tilt, y, 1.6), 0.12)
        comp.bRepBodies.add(cyl, base)
    base.finishEdit()
    names = []
    for i, b in enumerate(base.bodies):
        b.name = 'Pin %d' % (i + 1)
        names.append(b.name)
    return {'doc': doc.name, 'bodies': names}


def dev_close(app, job):
    doc = app.activeDocument
    name = doc.name
    doc.close(False)
    return {'closed': name, 'now': app.activeDocument.name if app.activeDocument else None}


def dev_activate(app, job):
    for i in range(app.documents.count):
        doc = app.documents.item(i)
        if doc.name.startswith(job['name']):
            doc.activate()
            return {'active': app.activeDocument.name}
    raise ValueError('no open document named %r' % job['name'])


def dev_cylinders(app, job):
    """Axis, radius, length and centre of each named body's cylindrical face."""
    out = []
    for b in _bodies_named(app, job['bodies']):
        bb = adsk.fusion.BRepBody.cast(b)
        row = {'name': b.name}
        for f in bb.faces:
            cy = adsk.core.Cylinder.cast(f.geometry)
            if cy:
                pts = body_points(bb)
                a = _unit(_tup(cy.axis))
                o = _tup(cy.origin)
                hs = [_dot(_sub(p, o), a) for p in pts]
                c = tuple(sum(p[k] for p in pts) / len(pts) * CM for k in range(3))
                row.update({'axis': [round(v, 4) for v in a], 'diameter_mm': round(2 * cy.radius * CM, 3),
                            'length_mm': round((max(hs) - min(hs)) * CM, 2),
                            'centre_mm': [round(v, 1) for v in c]})
        out.append(row)
    return out


def dev_clone(app, job):
    """Copy named bodies of the active document into a NEW document (temp BRep
    copies), so a build can be tried on real geometry without touching the user's."""
    tbm = adsk.fusion.TemporaryBRepManager.get()
    copies = [(b.name, tbm.copy(b)) for b in _bodies_named(app, job['bodies'])]
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    des = _design(app)
    if not is_parametric(des):
        des.designType = adsk.fusion.DesignTypes.ParametricDesignType
    comp = des.rootComponent
    base = comp.features.baseFeatures.add()
    base.startEdit()
    for _n, c in copies:
        comp.bRepBodies.add(c, base)
    base.finishEdit()
    names = []
    for (n, _c), b in zip(copies, base.bodies):
        b.name = n
        names.append(b.name)
    return {'doc': doc.name, 'bodies': names}


def dev_extents(app, job):
    """Each body's extent along the stack direction (mm from the stack origin)."""
    bodies = _bodies_named(app, job['bodies'])
    fr = setup(bodies, None, dict(DEFAULTS))
    out = []
    for b in bodies:
        hs = [_dot(_sub(p, fr['origin']), fr['d']) * CM for p in body_points(b)]
        out.append({'name': b.name, 'from': round(min(hs), 2), 'to': round(max(hs), 2)})
    return {'d': fr['d'], 'bodies': out}


DEV_COMMANDS = {
    'extents': dev_extents,
    'clone': dev_clone,
    'cylinders': dev_cylinders,
    'activate': dev_activate,
    'ping': dev_ping,
    'run': dev_run,
    'clear': dev_clear,
    'snapshot': dev_snapshot,
    'make_test': dev_make_test,
    'close': dev_close,
}
