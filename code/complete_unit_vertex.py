"""Complete a partial unit-margin polytope to a full tristochastic face.

Explicit order-2n construction using two bipartite edge colorings.
Existing fractional entries are unchanged; every new entry is 0 or 1.
The affine dimension is preserved, including dimension-zero vertices.
The method is a Latin-framework embedding; no asymptotic decomposition bound.
"""
from collections import defaultdict,deque
from pathlib import Path
from fractions import Fraction
from math import gcd,lcm
import argparse,json
from topology import incidence,rank
ROOT=Path(__file__).resolve().parents[1]


def matching(counts,size):
    # Hopcroft--Karp on the support of a regular bipartite multigraph.
    adj=[list(row) for row in counts];left=[-1]*size;right=[-1]*size
    while True:
        dist=[-1]*size;q=deque()
        for a in range(size):
            if left[a]<0:dist[a]=0;q.append(a)
        found=False
        while q:
            a=q.popleft()
            for b in adj[a]:
                aa=right[b]
                if aa<0:found=True
                elif dist[aa]<0:dist[aa]=dist[a]+1;q.append(aa)
        if not found:break
        def dfs(a):
            for b in adj[a]:
                aa=right[b]
                if aa<0 or (dist[aa]==dist[a]+1 and dfs(aa)):
                    left[a]=b;right[b]=a;return True
            dist[a]=-1;return False
        for a in range(size):
            if left[a]<0:dfs(a)
    assert all(b>=0 for b in left)
    return left


def edge_colors(rows,degree,size):
    counts=[{b:1 for b in row} for row in rows]+[{} for _ in range(size-len(rows))]
    coldeg=[0]*size
    for row in counts:
        for b,z in row.items():coldeg[b]+=z
    deficits=[degree-z for z in coldeg]
    cursor=0
    for row in counts:
        need=degree-sum(row.values())
        while need:
            while deficits[cursor]==0:cursor+=1
            z=min(need,deficits[cursor]);row[cursor]=row.get(cursor,0)+z;need-=z;deficits[cursor]-=z
    assert all(z==0 for z in deficits)
    colors=[]
    for c in range(degree):
        mat=matching(counts,size);colors.append(mat)
        for a,b in enumerate(mat):
            counts[a][b]-=1
            if not counts[a][b]:del counts[a][b]
    assert not any(counts)
    return colors


def complete(base,out_path):
    cells=[tuple(t) for t in base['cells']];m=base['denominator'];nums=base['numerators']
    n=max(max(t) for t in cells)+1;N=2*n
    AB=defaultdict(set);AC=defaultdict(set);BC=defaultdict(set);marg=defaultdict(int)
    for (a,b,c),z in zip(cells,nums):
        AB[a].add(b);AC[a].add(c);BC[b].add(c)
        for e in [(0,a,b),(1,a,c),(2,b,c)]:marg[e]+=z
    assert all(z==m for z in marg.values())
    assert all(len(AB[a])==len(AC[a]) for a in range(n))
    assert all(len(BC[b])==sum(b in AB[a] for a in range(n)) for b in range(n))
    assert all(sum(c in AC[a] for a in range(n))==sum(c in BC[b] for b in range(n)) for c in range(n))
    table=[[-1]*N for _ in range(N)]
    # Old row/old column cells outside the active graph receive new symbols.
    for a in range(n):
        for b in range(n):
            if b not in AB[a]:table[a][b]=n+(a+b)%n
    # Fill all old-row/new-column cells by edge-coloring available symbols.
    available=[]
    for a in range(n):
        used={z for z in table[a] if z>=0}
        available.append(set(range(n))-AC[a] | (set(range(n,N))-used))
        assert len(available[-1])==n
    coloring=edge_colors(available,n,N)
    for j,match in enumerate(coloring):
        for a in range(n):table[a][n+j]=match[a]
    # Each column has n available symbols after excluding reserved old BC edges.
    available=[]
    for b in range(N):
        used={table[a][b] for a in range(n) if table[a][b]>=0}
        reserved=BC[b] if b<n else set()
        available.append(set(range(N))-used-reserved)
        assert len(available[-1])==n
    coloring=edge_colors(available,n,N)
    for i,match in enumerate(coloring):
        for b in range(N):table[n+i][b]=match[b]
    # Independent final verification of all pair sums and degree bounds.
    new_count=0;seen=set(marg);original=set(marg)
    for a in range(N):
        for b in range(N):
            c=table[a][b]
            if c<0:
                assert a<n and b<n and b in AB[a];continue
            new_count+=1
            for e in [(0,a,b),(1,a,c),(2,b,c)]:
                assert e not in seen,'A completion entry overlaps another entry or an original pair'
                seen.add(e)
    assert len(seen)==3*N*N
    assert new_count==N*N-len({(a,b) for a in AB for b in AB[a]})
    dimension=base.get('dimension',0)
    rq=rank(incidence(cells)[-1]);assert rq==len(cells)-dimension
    assert lcm(*(m//gcd(z,m) for z in nums))==m
    # Each added coordinate has three private rows, so rank increases by one.
    result={'source_denominator':m,'order':N,'original_cells':cells,'original_numerators':nums,
            'integer_completion_by_AB':table,'integer_entries_added':new_count,'positive_entries':len(cells)+new_count,
            'rank_Q':rq+new_count,'dimension':dimension,'rank_proof':'original support independently ranked over Q; added columns have disjoint private pair rows',
            'max_line_degree':max(len(row) for row in incidence(cells)[-1]),
            'all_3N2_margins_one':True,'denominator_preserved':True}
    out_path.write_text(json.dumps(result,separators=(',',':'))+'\n')
    print({k:v for k,v in result.items() if k not in ('original_cells','original_numerators','integer_completion_by_AB')},flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('--output');a=p.parse_args()
    base=json.loads(Path(a.input).read_text());out=Path(a.output) if a.output else ROOT/f'data/full_unit_torsion_{base["denominator"]}.json'
    complete(base,out)
