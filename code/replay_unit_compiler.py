"""Replay compiler certificates from arrays alone, without importing construction.

Two-entry elimination here computes the integer presentation from actual pair
rows, including odd components, and checks its full reduced Smith form.
"""
from collections import defaultdict
from pathlib import Path
from fractions import Fraction
from math import gcd,lcm
import json
from topology import rank,smith_diagonal
ROOT=Path(__file__).resolve().parents[1]

def actual_rows(cells):
    rr=defaultdict(dict)
    for j,(a,b,c) in enumerate(cells):
        for e in [(0,a,b),(1,a,c),(2,b,c)]:rr[e][j]=1
    return list(rr.values())

def reduce_rows(rows,s):
    adj=[[] for _ in range(s)]
    for row in rows:
        if len(row)==2:
            a,b=row;adj[a].append(b);adj[b].append(a)
    roots=[None]*s;sign=[0]*s;count=0
    for v in range(s):
        if roots[v] is not None:continue
        roots[v]=count;sign[v]=1;stack=[v]
        while stack:
            a=stack.pop()
            for b in adj[a]:
                if roots[b] is None:roots[b]=count;sign[b]=-sign[a];stack.append(b)
        count+=1
    reduced=set()
    for row in rows:
        r=[0]*count
        for j in row:r[roots[j]]+=sign[j]
        if any(r):
            if next(z for z in r if z)<0:r=[-z for z in r]
            reduced.add(tuple(r))
    return count,sorted(reduced)

reports=[]
for path in sorted((ROOT/'data').glob('unit_torsion_*.json')):
    base=json.loads(path.read_text());cells=base['cells'];nums=base['numerators'];m=base['denominator'];s=len(cells)
    assert len({tuple(c) for c in cells})==s and all(0<z<m for z in nums)
    rows=actual_rows(cells)
    assert all(sum(nums[j] for j in row)==m for row in rows)
    assert max(map(len,rows))<=3
    rq=rank(rows);assert rq==s
    assert lcm(*(m//gcd(z,m) for z in nums))==m
    nr,reduced=reduce_rows(rows,s);diag=smith_diagonal(reduced)
    assert sum(z!=0 for z in diag)==nr and [z for z in diag if z>1]==[m]
    L=m.bit_length()-1;r=L+m.bit_count()-1-(m%2)
    assert s==189*r-(71 if m%2==0 else 0)
    report={'file':path.name,'m':m,'triangles':s,'rank_Q':rq,'actual_two_line_components':nr,
            'actual_reduced_rows':len(reduced),'actual_torsion_smith':[z for z in diag if z>1],
            'all_margins_one':True,'result':'PASS'}
    reports.append(report);print(report,flush=True)
(ROOT/'results/unit_compiler_independent_replay.json').write_text(json.dumps(reports,indent=2)+'\n')

completions=[]
for path in sorted((ROOT/'data').glob('full_unit_*.json')):
    b=json.loads(path.read_text());cells=b['original_cells'];nums=b['original_numerators'];m=b['source_denominator'];N=b['order']
    sums=defaultdict(int);degs=defaultdict(int)
    for (a,c,d),z in zip(cells,nums):
        for e in [(0,a,c),(1,a,d),(2,c,d)]:sums[e]+=z;degs[e]+=1
    assert all(z==m for z in sums.values())
    added=0
    for a,row in enumerate(b['integer_completion_by_AB']):
        assert len(row)==N
        for c,d in enumerate(row):
            if d<0:continue
            assert 0<=d<N;added+=1
            for e in [(0,a,c),(1,a,d),(2,c,d)]:
                assert e not in sums;sums[e]=m;degs[e]=1
    assert len(sums)==3*N*N and max(degs.values())<=3
    rq=rank(actual_rows(cells));dimension=len(cells)-rq
    assert dimension==b.get('dimension',0) and rq+added==b['rank_Q']
    assert lcm(*(m//gcd(z,m) for z in nums))==m
    result={'file':path.name,'order':N,'m':m,'dimension':dimension,'added_private_columns':added,
            'all_margins_one':True,'denominator_preserved':True,'result':'PASS'}
    completions.append(result);print(result,flush=True)
(ROOT/'results/full_unit_independent_replay.json').write_text(json.dumps(completions,indent=2)+'\n')
