"""Alpha Shape - the work. Re-imported on every command event and dev job.

Wraps the selected bodies (solid, surface or mesh) in their alpha shape: the
surfaces are sampled every `resolution` mm, a 3D Delaunay is built once, and the
tetrahedra whose circumscribed ball is no larger than `alpha` make the solid; its
boundary becomes a mesh body. A small alpha hugs each object, a large one tends
to the convex hull of the lot. The maths runs in as_engine.py, in a Python with
numpy and scipy, because Fusion's own has neither.

The preview follows the dialog: the tessellation is cached per selection and the
Delaunay per resolution (on disk), so moving alpha only re-filters.

All geometry crossing into Fusion is in centimetres; the engine works in mm.
"""

import adsk.core
import adsk.fusion

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
import traceback
import types

ADDIN_DIR = os.path.dirname(os.path.realpath(__file__))
ENGINE = os.path.join(ADDIN_DIR, 'as_engine.py')
CACHE_DIR = os.path.join(tempfile.gettempdir(), 'alpha_shape_cache')
SETTINGS_PATH = os.path.join(ADDIN_DIR, 'last_settings.json')
LOG_PATH = os.path.join(ADDIN_DIR, 'AlphaShape.log')
RES_DIR = os.path.join(ADDIN_DIR, 'resources', 'AlphaShape')
CM = 10.0
MM = 0.1
NO_WINDOW = 0x08000000 if os.name == 'nt' else 0

CMD_ID = 'AlphaShape_Create'
CMD_NAME = 'Alpha Shape'
PANELS = ['ParaMeshCreatePanel', 'SolidCreatePanel']

DEFAULTS = {'alpha_mm': 3.0, 'resolution_mm': 0.5, 'preview': True,
            'regular': True, 'drop_small': True, 'hide_inputs': False}
PREVIEW_COLOUR = (64, 160, 220)


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
        st.tess = {}          # selection key -> (coords mm, idx)
        st.caches = {}        # (selection key, resolution) -> npz path
        st.last = None        # (selection key, resolution, alpha, options) -> result
        st.report = ''
        sys.modules['_as_state'] = st
    return st


def load_settings():
    o = dict(DEFAULTS)
    try:
        with open(SETTINGS_PATH, encoding='utf-8') as fh:
            o.update(json.load(fh))
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
# the engine
# --------------------------------------------------------------------------- #

def python_cmd():
    """An interpreter with numpy and scipy: ALPHASHAPE_PYTHON, python_path.txt
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
            r = subprocess.run(c + ['-c', 'import numpy, scipy'], capture_output=True,
                               timeout=60, creationflags=NO_WINDOW)
            if r.returncode == 0:
                st.python = c
                return c
        except Exception:
            pass
    raise RuntimeError('No Python with numpy and scipy was found. Set ALPHASHAPE_PYTHON '
                       'or put the interpreter path in python_path.txt beside the add-in.')


def run_engine(job, timeout=900):
    os.makedirs(CACHE_DIR, exist_ok=True)
    jp = os.path.join(CACHE_DIR, 'job.json')
    op = os.path.join(CACHE_DIR, 'out.json')
    with open(jp, 'w', encoding='utf-8') as fh:
        json.dump(job, fh)
    if os.path.exists(op):
        os.remove(op)
    r = subprocess.run(python_cmd() + [ENGINE, jp, op], capture_output=True, text=True,
                       timeout=timeout, creationflags=NO_WINDOW)
    if not os.path.exists(op):
        raise RuntimeError('the alpha-shape engine produced nothing:\n%s' % r.stderr[-800:])
    with open(op, encoding='utf-8') as fh:
        res = json.load(fh)
    if not res.get('ok'):
        raise RuntimeError('the alpha-shape engine failed:\n%s' % res.get('error', '')[-800:])
    return res


# --------------------------------------------------------------------------- #
# the input bodies
# --------------------------------------------------------------------------- #

def selected_bodies(sel):
    return [sel.selection(i).entity for i in range(sel.selectionCount)]


def body_key(bodies, res_mm):
    toks = []
    for b in bodies:
        try:
            toks.append(b.entityToken)
        except Exception:
            toks.append(b.name)
    return hashlib.sha1(('|'.join(sorted(toks)) + '@%.4f' % res_mm).encode()).hexdigest()[:16]


def body_triangles(b, res_mm):
    """World triangles of a body, in mm."""
    mb = adsk.fusion.MeshBody.cast(b)
    if mb:
        pm = mb.mesh
        if pm and pm.triangleCount and not pm.quadCount and not pm.polygonCount:
            c, i = pm.nodeCoordinatesAsDouble, pm.triangleNodeIndices
        else:
            dm = mb.displayMesh
            c, i = dm.nodeCoordinatesAsDouble, dm.nodeIndices
        return [v * CM for v in c], list(i)
    bb = adsk.fusion.BRepBody.cast(b)
    calc = bb.meshManager.createMeshCalculator()
    # a tenth of the sampling step, but never finer than 0.05 mm: 0.01 mm on a
    # 30 mm ball once took Fusion down
    calc.surfaceTolerance = max(0.005, min(0.1, 0.1 * res_mm * MM))
    m = calc.calculate()
    return [v * CM for v in m.nodeCoordinatesAsDouble], list(m.nodeIndices)


def tessellate(bodies, res_mm):
    st = state()
    key = body_key(bodies, res_mm)
    if key not in st.tess:
        coords, idx = [], []
        for b in bodies:
            c, i = body_triangles(b, res_mm)
            base = len(coords) // 3
            coords.extend(c)
            idx.extend(base + k for k in i)
        st.tess = {key: (coords, idx)}          # one selection at a time
    return key, st.tess[key]


def compute(bodies, o):
    """Alpha shape of the bodies; cached at every level that can be."""
    st = state()
    res_mm, alpha_mm = float(o['resolution_mm']), float(o['alpha_mm'])
    key, (coords, idx) = tessellate(bodies, res_mm)
    sig = (key, res_mm, alpha_mm, bool(o['regular']), bool(o['drop_small']))
    if st.last and st.last[0] == sig:
        return st.last[1]
    job = {'alpha_mm': alpha_mm, 'regular': bool(o['regular']),
           'min_piece': 0.01 if o['drop_small'] else 0.0}
    cache = st.caches.get((key, res_mm))
    if cache and os.path.exists(cache):
        job['cache'] = cache
    else:
        job.update({'coords': coords, 'idx': idx, 'resolution_mm': res_mm,
                    'cache_dir': CACHE_DIR})
    t0 = time.time()
    res = run_engine(job)
    res['wall_s'] = round(time.time() - t0, 2)
    st.caches[(key, res_mm)] = res['cache']
    st.last = (sig, res)
    return res


# --------------------------------------------------------------------------- #
# preview graphics
# --------------------------------------------------------------------------- #

GROUP_ATTR = 'AlphaShapePreview'


def clear_graphics(app):
    try:
        des = adsk.fusion.Design.cast(app.activeProduct)
        if not des:
            return
        cg = des.rootComponent.customGraphicsGroups
        for i in range(cg.count - 1, -1, -1):
            g = cg.item(i)
            if g.id == GROUP_ATTR:
                g.deleteMe()
    except Exception:
        pass


def draw(app, res):
    clear_graphics(app)
    if not res.get('tris'):
        return
    des = adsk.fusion.Design.cast(app.activeProduct)
    g = des.rootComponent.customGraphicsGroups.add()
    g.id = GROUP_ATTR
    coords = adsk.fusion.CustomGraphicsCoordinates.create([v * MM for v in res['coords']])
    # per-vertex normals: without them the shell renders flat and the long
    # triangles an alpha shape spans between objects are hard to read
    nrm = res.get('normals') or []
    mesh = g.addMesh(coords, res['idx'], nrm, res['idx'] if nrm else [])
    c = PREVIEW_COLOUR
    mesh.color = adsk.fusion.CustomGraphicsBasicMaterialColorEffect.create(
        adsk.core.Color.create(c[0], c[1], c[2], 255), None, None, None, 5.0, 0.55)
    try:
        mesh.isSelectable = False
    except Exception:
        pass
    app.activeViewport.refresh()


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #

def describe(res, o):
    if not res.get('tris'):
        return ('<b>Nothing at alpha %.3f mm.</b> The smallest ball through the samples '
                'has radius %.3f mm: raise alpha above it, or refine the resolution.'
                % (o['alpha_mm'], res.get('r_min_mm', 0)))
    rows = [
        '%d sample points, %d tetrahedra kept, %d triangles'
        % (res['points'], res['tets'], res['tris']),
        'volume %.1f mm3, area %.1f mm2, %d piece(s)'
        % (res['volume_mm3'], res['area_mm2'], res['pieces']),
        ('<b>watertight</b>' if res['watertight'] else
         '<b>NOT watertight</b>: %d non-manifold spot(s), %d open edge(s)'
         % (res['nonmanifold_edges'], res['open_edges'])),
        '%.1f s' % res.get('wall_s', res.get('secs', 0)),
    ]
    for l in res.get('log', []):
        if l.startswith('regularised'):
            rows.append(l)
    return '<br>'.join(rows)


# --------------------------------------------------------------------------- #
# the dialog
# --------------------------------------------------------------------------- #

def on_created(app, cmd):
    o = load_settings()
    inputs = cmd.commandInputs
    cmd.setDialogInitialSize(360, 380)
    cmd.okButtonText = 'Create'
    cmd.isExecutedWhenPreEmpted = False
    sel = inputs.addSelectionInput('bodies', 'Bodies', 'The objects to wrap')
    sel.addSelectionFilter('SolidBodies')
    sel.addSelectionFilter('SurfaceBodies')
    sel.addSelectionFilter('MeshBodies')
    sel.setSelectionLimits(1, 0)
    a = inputs.addFloatSpinnerCommandInput('alpha', 'Alpha (ball radius)', 'mm',
                                           0.001, 1000.0, 0.5, o['alpha_mm'] * MM)
    a.tooltip = ('The radius of the probe ball. Space narrower than this between the '
                 'objects is filled; a small alpha hugs each object, a large one tends '
                 'to the convex hull of them all. Below an object\'s own radius it '
                 'falls apart into shreds.')
    r = inputs.addFloatSpinnerCommandInput('resolution', 'Resolution', 'mm',
                                           0.002, 2.0, 0.1, o['resolution_mm'] * MM)
    r.tooltip = ('Sampling step on the surfaces. Finer follows the objects more '
                 'closely and takes longer (a new resolution rebuilds the Delaunay; '
                 'a new alpha does not). Keep alpha well above it.')
    inputs.addBoolValueInput('preview', 'Preview', True, '', bool(o['preview']))
    adv = inputs.addGroupCommandInput('advanced', 'Options')
    adv.isExpanded = False
    ch = adv.children
    v = ch.addBoolValueInput('regular', 'Repair non-manifold spots', True, '',
                             bool(o['regular']))
    v.tooltip = ('Where two parts of the shape touch only along an edge or at a point, '
                 'remove the tetrahedra there, so the result is a closed manifold.')
    v = ch.addBoolValueInput('drop_small', 'Drop pieces under 1 % of the volume', True,
                             '', bool(o['drop_small']))
    ch.addBoolValueInput('hide_inputs', 'Hide the selected bodies afterwards', True, '',
                         bool(o['hide_inputs']))
    inputs.addTextBoxCommandInput('info', '', state().report or
                                  'Select the bodies to wrap.', 6, True)


def read_options(inputs):
    return {'alpha_mm': inputs.itemById('alpha').value * CM,
            'resolution_mm': inputs.itemById('resolution').value * CM,
            'preview': inputs.itemById('preview').value,
            'regular': inputs.itemById('regular').value,
            'drop_small': inputs.itemById('drop_small').value,
            'hide_inputs': inputs.itemById('hide_inputs').value}


def on_input_changed(app, inputs, changed_id):
    if changed_id == 'preview' and not inputs.itemById('preview').value:
        clear_graphics(app)


def on_validate(inputs):
    return (inputs.itemById('bodies').selectionCount > 0
            and inputs.itemById('alpha').value > 0
            and inputs.itemById('resolution').value > 0)


def on_preview(app, inputs, args):
    if not inputs.itemById('preview').value:
        clear_graphics(app)
        return
    o = read_options(inputs)
    bodies = selected_bodies(inputs.itemById('bodies'))
    try:
        res = compute(bodies, o)
        text = describe(res, o)
    except Exception:
        log(traceback.format_exc())
        res, text = {}, 'Preview failed: %s' % short_error(traceback.format_exc())
    # any write to an input wipes custom graphics: text first, mesh last
    info = inputs.itemById('info')
    if info.formattedText != text:
        info.formattedText = text
    draw(app, res)


def create_body(app, bodies, o, res):
    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.rootComponent
    base = comp.features.baseFeatures.add() if is_parametric(des) else None
    if base:
        base.startEdit()
    try:
        mb = comp.meshBodies.addByTriangleMeshData([v * MM for v in res['coords']],
                                                   res['idx'], [], [])
    finally:
        if base:
            base.finishEdit()
    if base:
        base.name = 'Alpha shape'
        mbs = base.meshBodies               # a plain list, not a collection
        mb = mbs[0] if mbs else mb
    mb.name = 'Alpha shape a%.2f r%.2f mm' % (o['alpha_mm'], o['resolution_mm'])
    if o.get('hide_inputs'):
        for b in bodies:
            try:
                b.isLightBulbOn = False
            except Exception:
                pass
    return mb


def on_execute(app, ui, inputs):
    o = read_options(inputs)
    save_settings(o)
    bodies = selected_bodies(inputs.itemById('bodies'))
    clear_graphics(app)
    res = compute(bodies, o)
    if not res.get('tris'):
        ui.messageBox(describe(res, o).replace('<br>', '\n').replace('<b>', '')
                      .replace('</b>', ''), CMD_NAME)
        return
    mb = create_body(app, bodies, o, res)
    state().report = ('Last: "%s"<br>' % mb.name) + describe(res, o)
    log(json.dumps({k: v for k, v in res.items() if k not in ('coords', 'idx', 'normals')}))


def on_destroy(app):
    clear_graphics(app)


# --------------------------------------------------------------------------- #
# toolbar
# --------------------------------------------------------------------------- #

TOOLTIP = ('Wrap the selected bodies in their alpha shape, as a mesh body.\n\n'
           'The surfaces are sampled every "Resolution" mm and every gap a ball of '
           'radius "Alpha" cannot pass through is filled. Small alpha hugs each '
           'object; large alpha tends to the convex hull. The preview follows the '
           'values as you change them.')


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
    out = []
    comp = _design(app).rootComponent
    for n in names:
        b = comp.bRepBodies.itemByName(n) or comp.meshBodies.itemByName(n)
        if b is None:
            raise ValueError('no body named %r' % n)
        out.append(b)
    return out


def dev_ping(app, job):
    return {'doc': app.activeDocument.name if app.activeDocument else None,
            'engine': python_cmd()}


def dev_run(app, job):
    """Compute (and with create=True, build) the alpha shape of named bodies."""
    o = dict(DEFAULTS)
    o.update(job.get('options') or {})
    bodies = _bodies_named(app, job['bodies'])
    res = compute(bodies, o)
    out = {k: v for k, v in res.items() if k not in ('coords', 'idx', 'normals')}
    if job.get('draw'):
        draw(app, res)
        out['drawn'] = True
    if job.get('create') and res.get('tris'):
        mb = create_body(app, bodies, o, res)
        out['body'] = mb.name
        out['body_tris'] = mb.mesh.triangleCount
    return out


def dev_clear(app, job):
    clear_graphics(app)
    return {'cleared': True}


def dev_snapshot(app, job):
    """Save the current view as it is - the camera is not touched."""
    vp = app.activeViewport
    ok = vp.saveAsImageFile(job['path'], int(job.get('width', 1400)), int(job.get('height', 1000)))
    return {'saved': bool(ok), 'path': job['path']}


def dev_make_test(app, job):
    """A new document with three cylinders, so a build never touches the user's."""
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    des = _design(app)
    if not is_parametric(des):
        des.designType = adsk.fusion.DesignTypes.ParametricDesignType
    comp = des.rootComponent
    tbm = adsk.fusion.TemporaryBRepManager.get()
    base = comp.features.baseFeatures.add()
    base.startEdit()
    for k, (x, y) in enumerate(((0, 0), (0.8, 0.2), (0.3, 0.9))):
        cyl = tbm.createCylinderOrCone(adsk.core.Point3D.create(x, y, 0), 0.12,
                                       adsk.core.Point3D.create(x, y, 1.6), 0.12)
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


DEV_COMMANDS = {
    'make_test': dev_make_test,
    'close': dev_close,
    'snapshot': dev_snapshot,
    'ping': dev_ping,
    'run': dev_run,
    'clear': dev_clear,
}
