"""Compile literal sum-one equations into unit-margin tripartite supports.

Variable pieces are barycentric subdivisions of triangulated spheres with
separated triangular holes. A port exposes either t or 1-t. Gluing two or
three ports implements a sum-one equation in all three pair directions.
Exact verification includes absence of unintended edge identifications.
"""
from collections import defaultdict,Counter
from itertools import permutations,combinations
from fractions import Fraction
from pathlib import Path
import argparse,json
from topology import incidence,rank,smith_diagonal
ROOT=Path(__file__).resolve().parents[1]

class DSU:
    def __init__(self):self.p={}
    def find(self,x):
        if x not in self.p:self.p[x]=x
        if self.p[x]!=x:self.p[x]=self.find(self.p[x])
        return self.p[x]
    def union(self,a,b):self.p[self.find(a)]=self.find(b)


def parity(p):return (-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))


def sphere(k):
    # k separated port faces, with distinguished vertices in even-numbered rings.
    # Two rings per port ensure distinguished vertices are pairwise nonadjacent.
    rings=2*k
    faces=[]
    for j in range(rings-1):
        for a in range(3):
            b=(a+1)%3
            faces.extend([(3*j+a,3*j+b,3*(j+1)+b),
                          (3*j+a,3*(j+1)+b,3*(j+1)+a)])
    faces.extend([(0,2,1),(3*(rings-1),3*(rings-1)+1,3*(rings-1)+2)])
    # Orient coherently using the dual graph, independent of initial ordering.
    ef=defaultdict(list)
    for i,f in enumerate(faces):
        for a,b in [(f[0],f[1]),(f[1],f[2]),(f[2],f[0])]:ef[tuple(sorted((a,b)))].append((i,1 if a<b else -1))
    assert all(len(z)==2 for z in ef.values())
    adj=defaultdict(list)
    for z in ef.values():
        (i,si),(j,sj)=z;adj[i].append((j,-si*sj));adj[j].append((i,-si*sj))
    orient={0:1};stack=[0]
    while stack:
        i=stack.pop()
        for j,factor in adj[i]:
            new=orient[i]*factor
            if j in orient:assert orient[j]==new
            else:orient[j]=new;stack.append(j)
    assert len(orient)==len(faces)
    chosen=[(6*i,6*i+1,3*(2*i+1)+1) for i in range(k)]
    fset={tuple(sorted(f)):i for i,f in enumerate(faces)}
    return faces,orient,[(fset[tuple(sorted(f))],f[0]) for f in chosen]


def variable_piece(var,signs):
    # signs +1 expose t, -1 expose 1-t along all three edges of the port.
    faces,ori,chosen=sphere(len(signs))
    triangles=[];colors=[];face_flags={}
    for i,f in enumerate(faces):
        for p in permutations(range(3)):
            vv=f[p[0]];edge=tuple(sorted((f[p[0]],f[p[1]])))
            tri=((0,var,'v',vv),(1,var,'e',*edge),(2,var,'f',i))
            face_flags[(i,vv,edge)]=len(triangles)
            triangles.append(tri);colors.append(ori[i]*parity(p))
    edges=defaultdict(list)
    for i,tri in enumerate(triangles):
        for e in combinations(tri,2):edges[e].append(i)
    assert all(len(z)==2 and colors[z[0]]==-colors[z[1]] for z in edges.values())
    ports=[];deleted=set()
    for (fi,v),sign in zip(chosen,signs):
        options=[j for (i,vv,e),j in face_flags.items() if i==fi and vv==v and colors[j]==-sign]
        assert len(options)==1
        j=options[0];deleted.add(j);ports.append(triangles[j])
    assert len({v for p in ports for v in p})==3*len(ports)
    kept=[j for j in range(len(triangles)) if j not in deleted]
    # The dual graph of the punctured sphere is connected.
    adj=defaultdict(set)
    for z in edges.values():
        zz=[j for j in z if j not in deleted]
        if len(zz)==2:a,b=zz;adj[a].add(b);adj[b].add(a)
    seen={kept[0]};stack=[kept[0]]
    while stack:
        for j in adj[stack.pop()]-seen:seen.add(j);stack.append(j)
    assert seen==set(kept)
    return [triangles[j] for j in kept],[colors[j] for j in kept],ports


def binary_system(m):
    assert m>=2
    values=[Fraction(1,m)];coefficients=[1];equations=[]
    current=0
    bits=bin(m)[3:]
    for i,bit in enumerate(bits):
        last=i==len(bits)-1
        if last:
            equations.append([(current,1),(current,1)]+([(0,1)] if bit=='1' else []));break
        doubled=len(values);coefficients.append(2*coefficients[current]);values.append(2*values[current])
        equations.append([(current,1),(current,1),(doubled,-1)])
        if bit=='1':
            added=len(values);coefficients.append(coefficients[doubled]+1);values.append(values[doubled]+values[0])
            equations.append([(doubled,1),(0,1),(added,-1)]);current=added
        else:current=doubled
    assert len(values)==len(equations) and all(0<x<1 for x in values)
    assert all(sum(values[v] if sign==1 else 1-values[v] for v,sign in eq)==1 for eq in equations)
    return values,equations,coefficients


def compile_system(values,equations,smith_limit=300,expected_dimension=0):
    occurrences=defaultdict(list)
    for ei,eq in enumerate(equations):
        assert len(eq) in (2,3)
        for oi,(v,sign) in enumerate(eq):occurrences[v].append((ei,oi,sign))
    triangles=[];weights=[];coordinate_variables=[];coordinate_signs=[];ports={};dsu=DSU();pieces=[]
    for v,value in enumerate(values):
        occ=occurrences[v]
        tt,cc,pp=variable_piece(v,[sign for _,_,sign in occ])
        triangles.extend(tt);weights.extend(value if sign==1 else 1-value for sign in cc)
        coordinate_variables.extend([v]*len(tt));coordinate_signs.extend(cc)
        for (ei,oi,sign),port in zip(occ,pp):ports[ei,oi]=port
        pieces.append({'variable':v,'ports':len(pp),'triangles':len(tt)})
    # Identify boundary vertices by their three colors.
    for ei,eq in enumerate(equations):
        for oi in range(1,len(eq)):
            for a,b in zip(ports[ei,0],ports[ei,oi]):dsu.union(a,b)
    converted=[tuple(dsu.find(v) for v in tri) for tri in triangles]
    assert len(set(converted))==len(converted),'Duplicate triangle after gluing'
    assert all(len(set(tri))==3 for tri in converted)
    # All identifications of pre-gluing edges must be intentional port identifications.
    old_edges=defaultdict(list)
    for j,tri in enumerate(triangles):
        for e in combinations(tri,2):old_edges[e].append(j)
    merged=defaultdict(list)
    for e,fs in old_edges.items():merged[tuple(dsu.find(v) for v in e)].append((e,fs))
    expected={}
    for ei,eq in enumerate(equations):
        for ab in combinations(range(3),2):
            key=tuple(dsu.find(ports[ei,0][a]) for a in ab)
            assert key not in expected,'Different junctions collided'
            expected[key]=len(eq)
    for edge,old in merged.items():
        if edge in expected:
            assert len(old)==expected[edge] and all(len(fs)==1 for _,fs in old)
        else:assert len(old)==1 and len(old[0][1])==2,'Unintended edge merge'
    labels=[{v:i for i,v in enumerate(sorted({tri[a] for tri in converted}))} for a in range(3)]
    cells=[tuple(labels[a][tri[a]] for a in range(3)) for tri in converted]
    verts,edges,_,_,rows=incidence(cells)
    assert all(sum(weights[j] for j in row)==1 for row in rows)
    assert all(2<=len(row)<=3 for row in rows)
    rq=rank(rows);assert rq==len(cells)-expected_dimension
    from math import lcm
    denominator=lcm(*(x.denominator for x in weights))
    rr=[];rhs=[]
    for eq in equations:
        row=[0]*len(values)
        for v,sign in eq:row[v]+=sign
        rr.append(row);rhs.append(1-sum(sign==-1 for v,sign in eq))
    reduced_snf=smith_diagonal(rr)
    actual_snf=None
    if len(cells)<=smith_limit:
        diag=smith_diagonal([[row.get(j,0) for j in range(len(cells))] for row in rows])
        actual_snf=[z for z in diag if z>1]
        assert actual_snf==[z for z in reduced_snf if z>1]
    return {'grid':list(map(len,labels)),'counts':[len(verts),len(edges),len(cells)],
            'denominator':denominator,'cells':cells,'numerators':[int(x*denominator) for x in weights],
            'rank_Q':rq,'dimension':expected_dimension,'coordinate_variables':coordinate_variables,'coordinate_signs':coordinate_signs,'rank_at_denominator':rank(rows,denominator) if denominator in [2,3,5,7,11,13,17,19,31,97,257,65537] else None,
            'max_line_degrees':[max(len(row) for e,row in zip(edges,rows) if (e[0][0],e[1][0])==ab) for ab in [(0,1),(0,2),(1,2)]],
            'variable_values':list(map(str,values)),'literal_equations':equations,'reduced_matrix':rr,'reduced_rhs':rhs,
            'reduced_smith':reduced_snf,'full_torsion_smith':actual_snf,'pieces':pieces,
            'checked':'exact unit margins, injective triangle map, only intended edge merges, exact predicted Q rank'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--orders',type=int,nargs='+',default=[2,3,4,5,7,8,11,13]);p.add_argument('--smith-limit',type=int,default=300);a=p.parse_args()
    reports=[]
    for m in a.orders:
        values,eq,coeff=binary_system(m);r=compile_system(values,eq,a.smith_limit);assert r['denominator']==m
        assert [z for z in r['reduced_smith'] if z>1]==[m]
        (ROOT/f'data/unit_torsion_{m}.json').write_text(json.dumps(r,indent=2)+'\n')
        report={k:v for k,v in r.items() if k not in ['cells','numerators','coordinate_variables','coordinate_signs']};reports.append(report)
        print('m',m,'grid',r['grid'],'s',r['counts'][2],'rank',r['rank_Q'],'degrees',r['max_line_degrees'],'SNF',r['full_torsion_smith'],flush=True)
    (ROOT/'results/unit_torsion_construction.json').write_text(json.dumps(reports,indent=2)+'\n')
if __name__=='__main__':main()
