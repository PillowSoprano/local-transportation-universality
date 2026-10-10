"""Build slim arrays for positive 1-in-3 formulas (Theorem 1.3).

Each clause t_i+t_j+t_k=1 is encoded by the matching gadgets of Lemma 3.2
(Prusa-Werner), and each hyperedge is placed on a cycle of six cells
(Lemma 3.1). The output has format 3 x 12q x 9q, line sums in {0,1,2},
and the fractional witness obtained from t_i = 1/3.
"""
import argparse
import json
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

FOUR = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
SATISFIABLE = [[0, 1, 2], [2, 3, 4]]


def hypergraph(clauses):
    q = len(clauses)
    occ = [(s, j) for s in range(q) for j in (1, 2, 3)]
    var = {(s, j): clauses[s][j - 1] for s, j in occ}
    by_var = {}
    for o in occ:
        by_var.setdefault(var[o], []).append(o)
    pi = {}
    for lst in by_var.values():
        for k, o in enumerate(lst):
            pi[o] = lst[(k + 1) % len(lst)]

    def vid(cls, s, j):
        return 4 * q * cls + 4 * s + j

    edges, kinds = [], []
    for s, j in occ:
        ps, pj = pi[(s, j)]
        edges.append((vid(0, s, 0), vid(1, s, j), vid(2, ps, pj))); kinds.append(('P', var[(s, j)]))
        edges.append((vid(0, s, j), vid(1, s, 0), vid(2, s, 0))); kinds.append(('Q', var[(s, j)]))
        edges.append((vid(0, s, j), vid(1, s, j), vid(2, s, j))); kinds.append(('R', var[(s, j)]))
    return 12 * q, edges, kinds


def array(nv, edges, weights):
    """Six cells per hyperedge; weights are the matching weights t_e."""
    cells = {}
    for k, ((a, b, c), w) in enumerate(zip(edges, weights)):
        for i, v, val in [(0, a, 1 - w), (0, c, w), (1, b, w), (1, c, 1 - w), (2, a, w), (2, b, 1 - w)]:
            cells[(i, v, k)] = val
    return cells


def line_sums(nv, edges):
    deg = [0] * nv
    for e in edges:
        for v in e:
            deg[v] += 1
    q4 = nv // 3
    iv = {}
    for v in range(nv):
        d = deg[v]
        cls = v // q4
        iv[v] = [[d - 1, 0, 1], [0, 1, d - 1], [1, d - 1, 0]][cls]
    ve = {v: [0] * len(edges) for v in range(nv)}
    for k, e in enumerate(edges):
        for v in e:
            ve[v][k] = 1
    return iv, ve, deg


def build(clauses, name):
    nv, edges, kinds = hypergraph(clauses)
    iv, ve, deg = line_sums(nv, edges)
    # Matching weights at t_i = 1/3, as numerators over 3: p_o = q_o = 1/3, r_o = 2/3.
    weight_third = [1 if kind in ('P', 'Q') else 2 for kind, _ in kinds]
    cells = sorted(array(nv, edges, [0] * len(edges)).keys())
    wit = {}
    for k, ((a, b, c), w) in enumerate(zip(edges, weight_third)):
        for i, v, val in [(0, a, 3 - w), (0, c, w), (1, b, w), (1, c, 3 - w), (2, a, w), (2, b, 3 - w)]:
            wit[(i, v, k)] = val
    nvars = 1 + max(max(c) for c in clauses)
    sat = []
    for bits in product([0, 1], repeat=nvars):
        if all(sum(bits[x] for x in cl) == 1 for cl in clauses):
            ws = [bits[x] if kind in ('P', 'Q') else 1 - bits[x] for kind, x in kinds]
            arr = array(nv, edges, ws)
            sat.append({'assignment': list(bits), 'ones': sorted([list(c) for c, val in arr.items() if val == 1])})
    record = {
        'name': name,
        'clauses': clauses,
        'grid': [3, nv, len(edges)],
        'line_sums_IV': [iv[v] for v in range(nv)],
        'line_sums_VE': [ve[v] for v in range(nv)],
        'line_sums_IE': [[1] * len(edges) for _ in range(3)],
        'allowed_cells': [list(c) for c in cells],
        'witness_denominator': 3,
        'witness_numerators': [wit[c] for c in cells],
        'integer_points_from_satisfying_assignments': sat,
    }
    out = ROOT / f'data/slim_hardness_{name}.json'
    out.write_text(json.dumps(record) + '\n')
    print(f'{name}: q={len(clauses)} format 3x{nv}x{len(edges)}, '
          f'{len(cells)} allowed cells, {len(sat)} satisfying assignments -> {out.name}')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--clauses', type=str, help='JSON list of clauses, e.g. [[0,1,2],[2,3,4]]')
    p.add_argument('--name', type=str, default='custom')
    a = p.parse_args()
    if a.clauses:
        build(json.loads(a.clauses), a.name)
    else:
        build(FOUR, 'four')
        build(SATISFIABLE, 'satisfiable')


if __name__ == '__main__':
    main()
