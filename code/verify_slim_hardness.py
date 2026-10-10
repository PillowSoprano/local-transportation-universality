"""Independently check the slim certificates for Theorem 1.3.

Reads only the stored line sums, allowed cells, witness, and integer points;
it does not import the construction code.
"""
import json
from collections import defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def rank_Q(rows, ncols):
    rows = [dict(r) for r in rows if r]
    r = 0
    pivots = {}
    for row in rows:
        row = {k: Fraction(v) for k, v in row.items() if v}
        for col, prow in pivots.items():
            if col in row:
                f = row[col]
                for k, v in prow.items():
                    row[k] = row.get(k, 0) - f * v
                    if row[k] == 0:
                        del row[k]
        if row:
            col = min(row)
            f = row[col]
            row = {k: v / f for k, v in row.items()}
            for c2, prow in pivots.items():
                if col in prow:
                    g = prow[col]
                    for k, v in row.items():
                        prow[k] = prow.get(k, 0) - g * v
                        if prow[k] == 0:
                            del prow[k]
            pivots[col] = row
            r += 1
    return r


def check(name):
    d = json.loads((ROOT / f'data/slim_hardness_{name}.json').read_text())
    I, V, E = d['grid']
    iv, ve, ie = d['line_sums_IV'], d['line_sums_VE'], d['line_sums_IE']
    assert I == 3 and V == 12 * len(d['clauses']) and E == 9 * len(d['clauses'])
    vals_iv = {x for row in iv for x in row}
    vals_rest = {x for row in ve for x in row} | {x for row in ie for x in row}
    assert vals_iv <= {0, 1, 2} and vals_rest <= {0, 1}
    # allowed cells are exactly the triples whose three line sums are positive
    allowed = {(i, v, e) for i in range(I) for v in range(V) for e in range(E)
               if iv[v][i] > 0 and ie[i][e] > 0 and ve[v][e] > 0}
    stored = {tuple(c) for c in d['allowed_cells']}
    assert allowed == stored
    # line bounds on the allowed cells
    lines = {'IV': defaultdict(list), 'IE': defaultdict(list), 'VE': defaultdict(list)}
    for idx, (i, v, e) in enumerate(d['allowed_cells']):
        lines['IV'][(i, v)].append(idx); lines['IE'][(i, e)].append(idx); lines['VE'][(v, e)].append(idx)
    bounds = tuple(max(len(x) for x in lines[k].values()) for k in ('IV', 'IE', 'VE'))
    assert bounds[0] <= 3 and bounds[1] <= 2 and bounds[2] <= 2
    # the witness: values 1/3 and 2/3, positive on every allowed cell, exact line sums
    den = d['witness_denominator']
    x = [Fraction(n, den) for n in d['witness_numerators']]
    assert set(x) <= {Fraction(1, 3), Fraction(2, 3)}
    target = {'IV': lambda k: iv[k[1]][k[0]], 'IE': lambda k: ie[k[0]][k[1]], 'VE': lambda k: ve[k[0]][k[1]]}
    for fam in lines:
        for key, idxs in lines[fam].items():
            assert sum(x[j] for j in idxs) == target[fam](key)
    # lines with positive sum but no allowed cell would make the polytope empty
    for v in range(V):
        for i in range(I):
            assert iv[v][i] == 0 or (i, v) in lines['IV']
    # dimension of the polytope from the rank of its equations
    rows = [{j: 1 for j in idxs} for fam in lines for idxs in lines[fam].values()]
    rk = rank_Q(rows, len(x))
    dim = len(x) - rk
    # integer points coming from satisfying assignments
    n = 1 + max(max(c) for c in d['clauses'])
    sat = [bits for bits in product([0, 1], repeat=n)
           if all(sum(bits[t] for t in cl) == 1 for cl in d['clauses'])]
    pts = d['integer_points_from_satisfying_assignments']
    assert sorted(tuple(p['assignment']) for p in pts) == sorted(sat)
    index = {tuple(c): j for j, c in enumerate(d['allowed_cells'])}
    for p in pts:
        y = [0] * len(x)
        for c in p['ones']:
            y[index[tuple(c)]] = 1
        for fam in lines:
            for key, idxs in lines[fam].items():
                assert sum(y[j] for j in idxs) == target[fam](key)
    if dim == 0:
        assert not pts, 'a singleton with a fractional point has no integer point'
    summary = {'name': name, 'format': d['grid'], 'allowed_cells': len(x), 'line_bounds': bounds,
               'line_sum_values': {'IV': sorted(vals_iv), 'IE_VE': sorted(vals_rest)},
               'witness_values': ['1/3', '2/3'], 'dimension': dim,
               'satisfying_assignments': len(sat), 'integer_points_checked': len(pts)}
    print(f"{name}: format {d['grid'][0]}x{d['grid'][1]}x{d['grid'][2]}, line bounds {bounds}, "
          f"witness in {{1/3,2/3}}, dimension {dim}, {len(sat)} satisfying assignments, "
          f"{len(pts)} integer points verified")
    return summary


def main():
    out = [check('four'), check('satisfiable')]
    four = out[0]
    assert four['format'] == [3, 48, 36] and four['dimension'] == 0 and four['integer_points_checked'] == 0
    (ROOT / 'results/slim_hardness_replay.json').write_text(json.dumps(out, indent=2) + '\n')


if __name__ == '__main__':
    main()
