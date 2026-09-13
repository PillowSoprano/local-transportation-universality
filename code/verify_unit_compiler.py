"""Verify affine compilation on a simplex, square, and pentagon.

Also barycentrically encode each entire polytope using margins 0,1,2 and
line bounds (3,2,2); verify that positive pair marginals alone enforce flags.
"""
from fractions import Fraction
from itertools import product,permutations,combinations
from collections import defaultdict
from pathlib import Path
import json
from construct_unit_torsion import compile_system,parity
from topology import incidence,rank
ROOT=Path(__file__).resolve().parents[1]


def subdivide(record):
    original=[tuple(t) for t in record['cells']]
    vars0=record['coordinate_variables'];signs0=record['coordinate_signs']
    vertices,edges,faces,ef,rows=incidence(original)
    vi={v:i for i,v in enumerate(vertices)};ei={e:i for i,e in enumerate(edges)}
    cells=[];variables=[];signs=[]
    for fi,face in enumerate(faces):
        for p in permutations(range(3)):
            cell=(vi[face[p[0]]],ei[tuple(sorted((face[p[0]],face[p[1]])))],fi)
            cells.append(cell);variables.append(vars0[fi]);signs.append(signs0[fi]*parity(p))
    values=list(map(Fraction,record['variable_values']))
    weights=[values[v] if z==1 else 1-values[v] for v,z in zip(variables,signs)]
    vv,ee,_,_,rr=incidence(cells)
    margins=[sum(weights[j] for j in row) for row in rr]
    assert all(v in (1,2) for v in margins)
    bypair={pair:max(len(row) for edge,row in zip(ee,rr) if (edge[0][0],edge[1][0])==pair) for pair in [(0,1),(0,2),(1,2)]}
    assert bypair[(0,1)]<=3 and bypair[(0,2)]==bypair[(1,2)]==2
    rq=rank(rr);assert len(cells)-rq==record['dimension']
    # Enumerate every triple consistent with its three nonzero pair margins.
    AB=defaultdict(set);AC=defaultdict(set);BC=set()
    for a,b,c in cells:AB[a].add(b);AC[a].add(c);BC.add((b,c))
    allowed={(a,b,c) for a in AB for b in AB[a] for c in AC[a] if (b,c) in BC}
    assert allowed==set(cells),'Structural zeros were not implied by margins'
    return {'cells':cells,'coordinate_variables':variables,'coordinate_signs':signs,
            'dimension':record['dimension'],'rank_Q':rq,'grid':[len(vertices),len(edges),len(faces)],
            'line_degrees':list(bypair.values()),'margin_values':sorted(set(map(int,margins))),
            'structural_zeros_implied_by_pair_margins':True,
            'margins':{str(edge):int(v) for edge,v in zip(ee,margins)}}


def verify_vertices(r,assignments):
    out=[]
    for values in assignments:
        weights=[values[v] if sign==1 else 1-values[v] for v,sign in zip(r['coordinate_variables'],r['coordinate_signs'])]
        assert all(0<=x<=1 for x in weights)
        allrows=incidence(r['cells'])[-1]
        margins=[sum(weights[j] for j in row) for row in allrows]
        assert all(v==1 for v in margins)
        support=[c for c,x in zip(r['cells'],weights) if x]
        assert rank(incidence(support)[-1])==len(support)
        out.append({'root_values':list(map(str,values)),'positive_entries':len(support),'vertex':True})
    return out


def run():
    F=Fraction
    cases=[]
    eq=[[(0,1),(1,1),(2,1)]]
    cases.append(('simplex',[F(1,3)]*3,eq,[[F(i==j) for i in range(3)] for j in range(3)]))
    eq=[[(0,1),(0,-1)],[(1,1),(1,-1)]]
    cases.append(('square',[F(1,2)]*2,eq,[list(map(F,p)) for p in product([0,1],repeat=2)]))
    # roots u,v,h,a,b,c,q,s; inequalities 0<=u,v<=1 and u+v<=3/2.
    eq=[[(2,1),(2,1)],[(3,1),(3,1),(0,-1)],[(4,1),(4,1),(1,-1)],
        [(3,1),(4,1),(5,-1)],[(6,1),(6,1),(2,-1)],[(5,1),(7,1),(6,1)]]
    def lift(u,v):
        u,v=F(u),F(v);return [u,v,F(1,2),u/2,v/2,(u+v)/2,F(1,4),F(3,4)-(u+v)/2]
    cases.append(('pentagon',lift(F(1,2),F(1,2)),eq,[lift(u,v) for u,v in [(0,0),(1,0),(1,F(1,2)),(F(1,2),1),(0,1)]]))
    summaries=[]
    for name,values,eq,verts in cases:
        r=compile_system(values,eq,smith_limit=0,expected_dimension=2)
        vertexchecks=verify_vertices(r,verts)
        s=subdivide(r)
        (ROOT/f'data/unit_compiler_{name}.json').write_text(json.dumps(r,indent=2)+'\n')
        (ROOT/f'data/integer_322_{name}.json').write_text(json.dumps(s,indent=2)+'\n')
        summary={'case':name,'unit_grid':r['grid'],'unit_support':r['counts'][2],'unit_rank':r['rank_Q'],
                 'vertices':vertexchecks,'subdivided_grid':s['grid'],'subdivided_support':len(s['cells']),
                 'subdivided_rank':s['rank_Q'],'subdivided_degree':s['line_degrees'],
                 'implicit_structural_zeros':True}
        summaries.append(summary);print(summary,flush=True)
    (ROOT/'results/unit_compiler_polytope_checks.json').write_text(json.dumps(summaries,indent=2)+'\n')
if __name__=='__main__':run()
