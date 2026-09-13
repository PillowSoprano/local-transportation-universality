"""Exact, dependency-free tools for tripartite two-complex certificates."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations


def incidence(cells):
    triangles = [tuple((i, v) for i, v in enumerate(t)) for t in cells]
    vertices = sorted({v for t in triangles for v in t})
    edge_faces = defaultdict(list)
    for j, t in enumerate(triangles):
        for edge in combinations(t, 2):
            edge_faces[edge].append(j)
    edges = sorted(edge_faces)
    rows = [{j: 1 for j in edge_faces[e]} for e in edges]
    return vertices, edges, triangles, edge_faces, rows


def rank(rows, prime=None):
    """Sparse Gaussian elimination over Q or a specified prime field.

    Python integers avoid overflow in modular multiplication. A modular
    rank is never labelled a rational rank.
    """
    basis = {}
    for source in rows:
        row = {j: (v % prime if prime else Fraction(v))
               for j, v in source.items() if v}
        row = {j: v for j, v in row.items() if v}
        while row:
            pivot = min(row)
            if pivot not in basis:
                d = row[pivot]
                inv = pow(d, -1, prime) if prime else 1 / d
                basis[pivot] = {j: (v * inv % prime if prime else v * inv)
                                for j, v in row.items()}
                break
            scale = row[pivot]
            for j, value in basis[pivot].items():
                v = row.get(j, 0) - scale * value
                if prime:
                    v %= prime
                if v:
                    row[j] = v
                else:
                    row.pop(j, None)
    return len(basis)


def components(adjacency):
    remaining = set(adjacency)
    answer = []
    while remaining:
        seen = set()
        stack = [min(remaining)]
        while stack:
            v = stack.pop()
            if v not in seen:
                seen.add(v)
                stack.extend(adjacency[v] - seen)
        remaining -= seen
        answer.append(sorted(seen))
    return answer


def surface_checks(vertices, triangles, edge_faces):
    links = {}
    for v in vertices:
        adjacency = defaultdict(set)
        for t in triangles:
            if v in t:
                a, b = [u for u in t if u != v]
                adjacency[a].add(b)
                adjacency[b].add(a)
        links[str(v)] = (len(components(adjacency)) == 1
                         and all(len(n) == 2 for n in adjacency.values()))
    dual = {i: set() for i in range(len(triangles))}
    for faces in edge_faces.values():
        for a, b in combinations(faces, 2):
            dual[a].add(b)
            dual[b].add(a)
    colors = {}
    bipartite = True
    for root in dual:
        if root in colors:
            continue
        colors[root] = 0
        stack = [root]
        while stack:
            a = stack.pop()
            for b in dual[a]:
                if b not in colors:
                    colors[b] = 1 - colors[a]
                    stack.append(b)
                elif colors[b] == colors[a]:
                    bipartite = False
    return {'edge_degrees': sorted({len(x) for x in edge_faces.values()}),
            'vertex_links_single_cycles': links,
            'dual_components': len(components(dual)),
            'dual_bipartite': bipartite}


def smith_diagonal(matrix):
    """Integer Smith invariants using unimodular Euclidean row/column steps.

    Intended for small certificates. Does not return transformation matrices.
    """
    a = [list(map(int, row)) for row in matrix]
    if not a:
        return []
    m, n = len(a), len(a[0])
    diagonal = []
    for k in range(min(m, n)):
        candidates = [(abs(a[i][j]), i, j) for i in range(k, m)
                      for j in range(k, n) if a[i][j]]
        if not candidates:
            break
        _, i, j = min(candidates)
        a[k], a[i] = a[i], a[k]
        for row in a:
            row[k], row[j] = row[j], row[k]
        while True:
            changed = False
            for i in range(k + 1, m):
                if a[i][k]:
                    q = a[i][k] // a[k][k]
                    a[i] = [x - q * y for x, y in zip(a[i], a[k])]
                    if a[i][k]:
                        a[k], a[i] = a[i], a[k]
                    changed = True
                    break
            if changed:
                continue
            for j in range(k + 1, n):
                if a[k][j]:
                    q = a[k][j] // a[k][k]
                    for row in a:
                        row[j] -= q * row[k]
                    if a[k][j]:
                        for row in a:
                            row[k], row[j] = row[j], row[k]
                    changed = True
                    break
            if changed:
                continue
            bad = next(((i, j) for i in range(k + 1, m)
                        for j in range(k + 1, n) if a[i][j] % a[k][k]), None)
            if bad is None:
                break
            i, _ = bad
            a[k] = [x + y for x, y in zip(a[k], a[i])]
        diagonal.append(abs(a[k][k]))
    assert all(b % a == 0 for a, b in zip(diagonal, diagonal[1:]))
    return diagonal


def sign_permutation(values):
    return (-1) ** sum(values[i] > values[j]
                       for i in range(len(values)) for j in range(i + 1, len(values)))
