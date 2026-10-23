"""Alpha Shape engine - runs in a Python with numpy and scipy.

    python as_engine.py job.json out.json

Two steps, so the preview can follow the alpha slider quickly:
  prepare   sample the selected surfaces at `resolution_mm`, 3D Delaunay, and
            cache the points, tetrahedra and circumradii (npz) under a key made of
            the input and the resolution; skipped when the cache exists
  shape     keep the tetrahedra whose circumradius is <= alpha and return the
            boundary triangles of that set, oriented outward

The alpha shape here is the classic one (Edelsbrunner & Muecke 1994) on a
sampled surface: alpha is the radius of the probe ball, so a small alpha hugs
every object and a large one tends to the convex hull. Everything is in mm.

Job:
    triangles   flat [x,y,z,...] per vertex and flat indices, every selected body
    resolution_mm, alpha_mm, cache_dir
"""

import hashlib
import json
import os
import sys
import time

import numpy as np
from scipy.spatial import Delaunay


def sample_surface(v, f, res, seed=20261023):
    """Points every ~res mm over the triangles, plus their sharp edges densely.

    Surface: area-weighted random samples at 4x the target density, thinned to
    one per res-sized voxel. Sharp edges (dihedral over 30 degrees, e.g. a
    cylinder's rim) are walked at res/2 and never thinned, so a rim survives
    any resolution."""
    rng = np.random.default_rng(seed)
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    cr = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(cr, axis=1)
    total = area.sum()
    n = int(min(4e6, max(1000, 4.0 * total / (res * res))))
    pick = rng.choice(len(f), size=n, p=area / total)
    r1 = np.sqrt(rng.random(n))
    r2 = rng.random(n)
    pts = ((1 - r1)[:, None] * a[pick] + (r1 * (1 - r2))[:, None] * b[pick]
           + (r1 * r2)[:, None] * c[pick])
    key = np.floor(pts / res).astype(np.int64)
    _u, first = np.unique(key, axis=0, return_index=True)
    surf = pts[first]

    # sharp edges
    nrm = cr / np.maximum(np.linalg.norm(cr, axis=1, keepdims=True), 1e-30)
    e = np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    fid = np.tile(np.arange(len(f)), 3)
    s = np.sort(e, axis=1)
    order = np.lexsort((s[:, 1], s[:, 0]))
    s, fid = s[order], fid[order]
    same = np.all(s[1:] == s[:-1], axis=1)
    i0 = np.where(same)[0]
    dots = (nrm[fid[i0]] * nrm[fid[i0 + 1]]).sum(1)
    sharp = s[i0][dots < np.cos(np.radians(30.0))]
    edge_pts = [v[np.unique(sharp)] if len(sharp) else np.zeros((0, 3))]
    for p, q in sharp:
        L = np.linalg.norm(v[q] - v[p])
        k = int(np.ceil(L / (0.5 * res)))
        if k > 1:
            t = np.linspace(0, 1, k + 1)[1:-1, None]
            edge_pts.append(v[p] + t * (v[q] - v[p]))
    edges = np.vstack(edge_pts)
    allp = np.vstack([surf, edges])
    # merge exact and near duplicates (Delaunay dislikes them)
    key = np.round(allp / (0.01 * res)).astype(np.int64)
    _u, first = np.unique(key, axis=0, return_index=True)
    return allp[np.sort(first)], len(surf), len(edges)


def circumradii(P, T):
    a = P[T[:, 0]]
    u, w_, z = P[T[:, 1]] - a, P[T[:, 2]] - a, P[T[:, 3]] - a
    vw, wu, uv = np.cross(w_, z), np.cross(z, u), np.cross(u, w_)
    den = 2.0 * (u * vw).sum(1)
    off = ((u * u).sum(1)[:, None] * vw + (w_ * w_).sum(1)[:, None] * wu
           + (z * z).sum(1)[:, None] * uv)
    with np.errstate(divide='ignore', invalid='ignore'):
        R = np.linalg.norm(off, axis=1) / np.abs(den)
    R[~np.isfinite(R)] = np.inf
    vol = np.abs(den) / 12.0
    return R, vol


def prepare(job, log):
    v = np.asarray(job['coords'], float).reshape(-1, 3)
    f = np.asarray(job['idx'], np.int64).reshape(-1, 3)
    res = float(job['resolution_mm'])
    h = hashlib.sha1(v.tobytes() + f.tobytes() + repr(res).encode()).hexdigest()[:16]
    path = os.path.join(job['cache_dir'], 'alpha_%s.npz' % h)
    if os.path.exists(path):
        log('cache hit %s' % os.path.basename(path))
        return path
    t0 = time.time()
    P, n_surf, n_edge = sample_surface(v, f, res)
    log('sampled %d points (%d surface every %.3f mm, %d on sharp edges) in %.2f s'
        % (len(P), n_surf, res, n_edge, time.time() - t0))
    t0 = time.time()
    # joggled: a cylinder's rim points are cospherical, and without QJ Qhull
    # leaves flat tetrahedra whose huge circumradii punch holes at any alpha
    T = Delaunay(P, qhull_options='QJ Qbb Qc').simplices.astype(np.int64)
    R, vol = circumradii(P, T)
    log('Delaunay: %d tetrahedra in %.2f s' % (len(T), time.time() - t0))
    os.makedirs(job['cache_dir'], exist_ok=True)
    np.savez(path, P=P, T=T, R=R, vol=vol)
    return path


class Topology:
    """Face identities of the Delaunay tetrahedra, computed ONCE per alpha.

    Recomputing them for every kept set (np.unique over all faces, rows at a
    time) was 80 % of a 5 s regularisation; with integer keys and the interior
    face pairs known in advance, a boundary is a bincount and the pieces a
    connected-components call."""

    def __init__(self, P, T):
        self.P, self.T = P, T
        n, N = len(P), len(T)
        self.faces = np.vstack([T[:, [1, 2, 3]], T[:, [0, 3, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 1]]])
        self.opp = np.concatenate([T[:, 0], T[:, 1], T[:, 2], T[:, 3]])
        self.tet = np.tile(np.arange(N), 4)
        s = np.sort(self.faces, axis=1).astype(np.int64)
        key = (s[:, 0] * n + s[:, 1]) * n + s[:, 2]
        _u, self.fid = np.unique(key, return_inverse=True)
        self.fid = self.fid.ravel()
        self.nf = int(self.fid.max()) + 1
        order = np.argsort(self.fid, kind='stable')
        sf = self.fid[order]
        same = np.where(sf[1:] == sf[:-1])[0]
        self.pair = np.c_[self.tet[order][same], self.tet[order][same + 1]]
        pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        te = np.sort(np.stack([T[:, [a, b]] for a, b in pairs], axis=1), axis=2).astype(np.int64)
        self.tkey = te[..., 0] * n + te[..., 1]

    def boundary(self, keep):
        """Indices (into self.faces) of the faces with exactly one kept tetrahedron."""
        m = keep[self.tet]
        cnt = np.bincount(self.fid[m], minlength=self.nf)
        return np.where(m & (cnt[self.fid] == 1))[0]

    def pieces(self, keep):
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        idx = np.where(keep)[0]
        loc = -np.ones(len(keep), np.int64)
        loc[idx] = np.arange(len(idx))
        pr = self.pair[keep[self.pair[:, 0]] & keep[self.pair[:, 1]]]
        g = coo_matrix((np.ones(len(pr)), (loc[pr[:, 0]], loc[pr[:, 1]])),
                       shape=(len(idx), len(idx)))
        n, lab = connected_components(g, directed=False)
        return idx, n, lab


def _bad_edges(F, npt):
    e = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1).astype(np.int64)
    ek, c = np.unique(e[:, 0] * npt + e[:, 1], return_counts=True)
    return ek[c > 2]


def _pinched_vertices(F):
    """Vertices whose boundary faces form more than one fan (two regions touching
    at a point). Vectorised: one node per (face, corner); the two corners of a
    vertex on either side of an edge shared by exactly two faces are joined; a
    vertex whose nodes fall in more than one component is pinched."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    nf = len(F)
    if not nf:
        return np.zeros(0, np.int64)
    # each directed edge (corner k -> k+1) of each face, with its corner nodes
    ea = np.concatenate([F[:, 0], F[:, 1], F[:, 2]])
    eb = np.concatenate([F[:, 1], F[:, 2], F[:, 0]])
    na = np.concatenate([np.arange(nf) * 3 + 0, np.arange(nf) * 3 + 1, np.arange(nf) * 3 + 2])
    nb = np.concatenate([np.arange(nf) * 3 + 1, np.arange(nf) * 3 + 2, np.arange(nf) * 3 + 0])
    key = np.minimum(ea, eb) * (int(F.max()) + 1) + np.maximum(ea, eb)
    order = np.argsort(key, kind='stable')
    k = key[order]
    _u, start, cnt = np.unique(k, return_index=True, return_counts=True)
    two = start[cnt == 2]
    i1, i2 = order[two], order[two + 1]
    # match corners by vertex: the same vertex sits at a or b of either half-edge
    rows, cols = [], []
    for x1, y1 in ((ea, na), (eb, nb)):
        for x2, y2 in ((ea, na), (eb, nb)):
            m = x1[i1] == x2[i2]
            rows.append(y1[i1][m])
            cols.append(y2[i2][m])
    r, c = np.concatenate(rows), np.concatenate(cols)
    g = coo_matrix((np.ones(len(r)), (r, c)), shape=(3 * nf, 3 * nf))
    _n, lab = connected_components(g, directed=False)
    vert = F.ravel()                      # node 3f+k is vertex F[f, k]
    pairs = np.unique(np.c_[vert, lab], axis=0)
    v, fans = np.unique(pairs[:, 0], return_counts=True)
    return v[fans > 1]


def regularise(topo, vol, keep, min_piece, log, max_iter=30):
    """Make the kept set bound a 2-manifold. Where two kept regions touch only
    along an edge (four boundary faces on it) or at a vertex (two fans), the kept
    tetrahedra on that edge or vertex are removed, and that is repeated until
    none is left. Removing beat adding: adding the tetrahedra round a bad edge
    cascaded, put up to 23 % on the volume and never converged; removing
    converged in 1-5 rounds for -0.2 % at most on 10 pins. Then pieces under
    `min_piece` of the kept volume are dropped."""
    T, npt = topo.T, len(topo.P)
    v0 = vol[keep].sum()
    removed, rounds = 0, 0
    for rounds in range(1, max_iter + 1):
        F = topo.faces[topo.boundary(keep)]
        if not len(F):
            break
        be = _bad_edges(F, npt)
        if len(be):
            hit = np.isin(topo.tkey, be).any(axis=1)
        else:
            pv = _pinched_vertices(F)
            if not len(pv):
                break
            hit = np.isin(T, pv).any(axis=1)
        hit &= keep
        if not hit.any():
            break
        keep = keep & ~hit
        removed += int(hit.sum())
    idx, n, lab = topo.pieces(keep)
    dropped = 0
    if n:
        pv = np.bincount(lab, weights=vol[idx], minlength=n)
        small = pv < min_piece * pv.sum()
        dropped = int(small.sum())
        if dropped:
            keep = keep.copy()
            keep[idx[small[lab]]] = False
    if removed or dropped:
        log('regularised: %d tetrahedra removed at non-manifold edges/vertices in '
            '%d round(s) (volume %+.2f %%), %d small piece(s) dropped'
            % (removed, rounds, 100.0 * (vol[keep].sum() / v0 - 1.0), dropped))
    return keep


def shape(cache, alpha, log, min_piece=0.01, regular=True):
    d = np.load(cache)
    P, T, R, vol = d['P'], d['T'], d['R'], d['vol']
    keep = R <= alpha
    if not keep.any():
        return {'coords': [], 'idx': [], 'tris': 0, 'tets': 0, 'volume_mm3': 0.0,
                'area_mm2': 0.0, 'watertight': False, 'nonmanifold_edges': 0,
                'open_edges': 0, 'pieces': 0, 'points': int(len(P)),
                'r_min_mm': float(R.min())}
    topo = Topology(P, T)
    if regular:
        keep = regularise(topo, vol, keep, min_piece, log)
    b = topo.boundary(keep)
    F = topo.faces[b].copy()
    o = topo.opp[b]
    # outward: the normal points away from the tetrahedron's fourth vertex
    a, bb, c = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]
    flip = (np.cross(bb - a, c - a) * (P[o] - a)).sum(1) > 0
    F[flip] = F[flip][:, [0, 2, 1]]
    used, Fi = np.unique(F, return_inverse=True)
    Fi = Fi.reshape(-1, 3)
    V = P[used]
    e = np.sort(np.vstack([Fi[:, [0, 1]], Fi[:, [1, 2]], Fi[:, [2, 0]]]), axis=1).astype(np.int64)
    _ue, ecnt = np.unique(e[:, 0] * len(V) + e[:, 1], return_counts=True)
    nonman = int((ecnt > 2).sum()) + int(len(_pinched_vertices(Fi)))
    open_e = int((ecnt == 1).sum())
    _idx, pieces, _lab = topo.pieces(keep)
    fn = np.cross(V[Fi[:, 1]] - V[Fi[:, 0]], V[Fi[:, 2]] - V[Fi[:, 0]])   # area-weighted
    vn = np.zeros_like(V)
    for k in range(3):
        np.add.at(vn, Fi[:, k], fn)
    vn /= np.maximum(np.linalg.norm(vn, axis=1, keepdims=True), 1e-30)
    return {'coords': V.ravel().tolist(), 'idx': Fi.ravel().tolist(),
            'normals': np.round(vn, 5).ravel().tolist(),
            'tris': int(len(Fi)), 'tets': int(keep.sum()), 'points': int(len(P)),
            'volume_mm3': float(vol[keep].sum()),
            'area_mm2': float(0.5 * np.linalg.norm(np.cross(V[Fi[:, 1]] - V[Fi[:, 0]],
                                                             V[Fi[:, 2]] - V[Fi[:, 0]]),
                                                    axis=1).sum()),
            'watertight': bool(open_e == 0 and nonman == 0),
            'nonmanifold_edges': nonman, 'open_edges': open_e, 'pieces': int(pieces),
            'r_min_mm': float(R.min())}


def main():
    job = json.load(open(sys.argv[1], encoding='utf-8'))
    lines = []

    def log(m):
        lines.append(m)

    t0 = time.time()
    cache = job.get('cache') or prepare(job, log)
    out = {'ok': True, 'cache': cache}
    if job.get('alpha_mm') is not None:
        t1 = time.time()
        out.update(shape(cache, float(job['alpha_mm']), log,
                         float(job.get('min_piece', 0.01)), bool(job.get('regular', True))))
        log('alpha %.3f mm: %d boundary triangles from %d tetrahedra in %.2f s'
            % (job['alpha_mm'], out['tris'], out['tets'], time.time() - t1))
    out['log'] = lines
    out['secs'] = round(time.time() - t0, 2)
    json.dump(out, open(sys.argv[2], 'w', encoding='utf-8'))


if __name__ == '__main__':
    try:
        main()
    except Exception:
        import traceback
        json.dump({'ok': False, 'error': traceback.format_exc()},
                  open(sys.argv[2], 'w', encoding='utf-8'))
