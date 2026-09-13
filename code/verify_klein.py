"""Reconstruct the Klein certificate, exact ranks, homology, and LP vertex."""
import json
from pathlib import Path
from topology import incidence, rank, smith_diagonal, surface_checks

ROOT = Path(__file__).resolve().parents[1]


def main():
    data = json.loads((ROOT / 'data/klein_4x4x3.json').read_text())
    v, e, t, ef, rows = incidence(data['cells'])
    checks = surface_checks(v, t, ef)
    rq = rank(rows)
    field_ranks = {str(p): rank(rows, p) for p in [2, 3, 5, 7, 101]}
    snf = smith_diagonal([[row.get(j, 0) for j in range(len(t))] for row in rows])
    assert (len(v), len(e), len(t), rq) == (11, 33, 22, 22)
    assert checks['edge_degrees'] == [2]
    assert all(checks['vertex_links_single_cycles'].values())
    assert checks['dual_components'] == 1 and not checks['dual_bipartite']
    assert field_ranks == {'2': 21, '3': 22, '5': 22, '7': 22, '101': 22}
    assert snf == [1] * 21 + [2]
    # Signed boundary is obtained by multiplying AC rows by -1.
    # Torsion(coker d2) = torsion(H1), since im d1 is free.
    b1 = len(e) - len(v) + 1 - rq
    all_margins = []
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        for a in range(data['grid'][i]):
            for b in range(data['grid'][j]):
                count = sum(s[i] == a and s[j] == b for s in data['cells'])
                assert count % 2 == 0
                all_margins.append({'axes': [i, j], 'indices': [a, b], 'value': count // 2})
    output = {'counts': {'vertices': len(v), 'edges': len(e), 'triangles': len(t)},
              'rank_Q': rq, 'rank_Fp': field_ranks, 'smith_diagonal': snf,
              'H1': {'free_rank': b1, 'torsion_invariant_factors': [2]},
              'H2_Z_rank': len(t) - rq, 'surface_checks': checks,
              'half_integral_vertex': {'grid': data['grid'], 'support': data['cells'],
                                       'value_on_support': '1/2', 'pair_margins': all_margins,
                                       'extremality': '22 positive-support columns independent over Q'}}
    (ROOT / 'results/klein_verified.json').write_text(json.dumps(output, indent=2) + '\n')
    print('Klein: (V,E,F)=(11,33,22); ranks Q/F2=22/21; H1=Z + Z/2; LP vertex verified.')


if __name__ == '__main__':
    main()
