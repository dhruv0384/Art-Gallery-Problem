# GroupID-22 (22114030_22114057_22114073)
# 24th Sept, 2025
# EdgeGuardSolver.py: mobile guard solver based on the recursive proof in "Art Galleries and Mobile Guards: Revisiting O'Rourke's Proof" by A. Biniaz.

from collections import defaultdict
import copy

DEBUG = False


def _coords_of_point(p):
    """Return a plain tuple coordinate for a DCEL.Point or a coord-tuple."""
    try:
        return tuple(p.coords)
    except Exception:
        return tuple(p)

def polygon_boundary_coords(d):
    """Return polygon boundary coords list (interior face outer boundary)."""
    faces = d.getFaces()
    if len(faces) >= 2:
        return [tuple(x) for x in faces[1].getOuterBoundaryCoords()]
    return [tuple(v.coords) for v in d.getVertices()]

def unordered_edge(a, b):
    """Order-invariant edge key (tuple of tuples)."""
    a_tup, b_tup = _coords_of_point(a), _coords_of_point(b)
    return (a_tup, b_tup) if a_tup <= b_tup else (b_tup, a_tup)

def small_side_vertices(poly_order, u, v):
    """Return the smaller side (inclusive endpoints) when splitting poly_order by u and v."""
    try:
        iu = poly_order.index(u)
        iv = poly_order.index(v)
    except ValueError:
        return poly_order[:]
    if iu > iv: iu, iv = iv, iu
    side1 = poly_order[iu:iv+1]
    side2 = poly_order[iv:] + poly_order[:iu+1]
    return side1 if len(side1) <= len(side2) else side2

def cut_off_count(poly_order, u, v):
    """Return number of boundary edges on the smaller side of (u,v)."""
    side = small_side_vertices(poly_order, u, v)
    return max(0, len(side) - 1)

def build_triangulation_edges(all_triangles):
    """Builds a set of all unique edges from a list of triangles."""
    tri_edges = set()
    for tri in all_triangles:
        coords = tuple(_coords_of_point(p) for p in tri)
        tri_edges.add(unordered_edge(coords[0], coords[1]))
        tri_edges.add(unordered_edge(coords[1], coords[2]))
        tri_edges.add(unordered_edge(coords[2], coords[0]))
    return tri_edges

def polygon_edge_set(poly_order):
    """Return set of unordered polygon boundary edges for current poly_order"""
    n = len(poly_order)
    s = set()
    for i in range(n):
        s.add(unordered_edge(poly_order[i], poly_order[(i+1) % n]))
    return s

def _solve_base_case(poly_order, V, all_tri_edges):
    """Handles the base cases for n=4 and n=5 as per the paper."""
    guards = set()
    n = len(poly_order)
    poly_edges = polygon_edge_set(poly_order)

    if n == 4:
        # For |V| >= 2, we need 2 guards. Two opposite edges cover all vertices.
        if len(V) >= 2:
            e1 = unordered_edge(poly_order[0], poly_order[1])
            e2 = unordered_edge(poly_order[2], poly_order[3])
            guards.add(e1)
            guards.add(e2)
        # For |V| <= 1, one guard suffices.
        else:
            if V:
                v_must_cover = list(V)[0]
                # Find an edge incident to the vertex in V
                for edge in poly_edges:
                    if v_must_cover in edge:
                        guards.add(edge)
                        break
            else: # V is empty, add any edge
                guards.add(next(iter(poly_edges)))
        return guards

    if n == 5:
        # Find a vertex with two diagonals incident to it.
        # This vertex and its diagonals form the triangulation of a 5-gon.
        diagonals = all_tri_edges - poly_edges
        v_counts = defaultdict(int)
        for u, v in diagonals:
            v_counts[u] += 1
            v_counts[v] += 1
        
        center_v = max(v_counts, key=v_counts.get)
        
        # Guard placement logic based on |V|
        if len(V) <= 1: # 1 guard
            v_to_cover = list(V)[0] if V else center_v
            # Add a diagonal/edge incident to v_to_cover
            for edge in all_tri_edges:
                if v_to_cover in edge:
                    guards.add(edge)
                    break
        elif len(V) <= 3: # 2 guards
            # Add two edges that cover all vertices in V.
            # A simple strategy is to pick edges incident to V vertices.
            covered = set()
            for v_must_cover in V:
                if v_must_cover in covered: continue
                for edge in poly_edges:
                    if v_must_cover in edge:
                        guards.add(edge)
                        covered.update(edge)
                        break
            # Ensure we have at least 2 guards if needed
            while len(guards) < 2 and len(poly_edges) > len(guards):
                guards.add(list(poly_edges - guards)[0])

        else: # |V| is 4 or 5, needs 3 guards
            # 3 edges can cover all 5 vertices.
            guards.add(unordered_edge(poly_order[0], poly_order[1]))
            guards.add(unordered_edge(poly_order[2], poly_order[3]))
            guards.add(unordered_edge(poly_order[4], poly_order[0]))

        return guards
    
    return None # Should not be reached

def _solve_guards(poly_order, all_tri_edges, V):
    """The core recursive function that implements the inductive proof."""
    n = len(poly_order)
    
    # Base Cases: n <= 5
    if n <= 5:
        return _solve_base_case(poly_order, V, all_tri_edges)

    # Inductive Step: n >= 6
    poly_edges = polygon_edge_set(poly_order)
    diagonals = all_tri_edges - poly_edges

    # Find a diagonal 'd' that cuts off 3 or 4 edges (Lemma 2)
    d = None
    cut_k = 0
    
    # Prefer a 3-edge cut
    for diag in diagonals:
        k = cut_off_count(poly_order, diag[0], diag[1])
        if k == 3:
            d = diag
            cut_k = 3
            break
    # Otherwise find a 4-edge cut
    if not d:
        for diag in diagonals:
            k = cut_off_count(poly_order, diag[0], diag[1])
            if k == 4:
                d = diag
                cut_k = 4
                break
    
    if not d:
        # Fallback, this shouldn't happen in a valid triangulation per Lemma 2
        # but as a safeguard, return a simple solution.
        return {next(iter(poly_edges))}

    # Get the vertices on the small side to identify x, y, w, z etc.
    small = small_side_vertices(poly_order, d[0], d[1])
    
    # --- Case 1: Diagonal d cuts off 3 edges ---
    if cut_k == 3:
        # Here small side has 4 vertices: w, x, y, z. d = (w,z)
        w, x, y, z = small[0], small[1], small[2], small[3]
        if d != unordered_edge(w,z): w,z = z,w # Ensure d = (w,z)
        
        P_prime = [p for p in poly_order if p not in {x, y}]
        
        # Case 1(a): x and y are NOT in V
        if x not in V and y not in V:
            V_prime = V.copy()
            V_prime.add(z) # Ensure z is guarded in the subproblem
            return _solve_guards(P_prime, all_tri_edges, V_prime)
        
        # Case 1(b): x or y (or both) are in V
        else:
            V_prime = V - {x, y}
            sub_guards = _solve_guards(P_prime, all_tri_edges, V_prime)
            # Add edge 'zy' as per paper to guard triangle xyz
            sub_guards.add(unordered_edge(z, y))
            return sub_guards

    # --- Case 2: Diagonal d cuts off 4 edges ---
    elif cut_k == 4:
        # Here small side has 5 vertices: v, w, x, y, z. d = (v,z)
        v, w, x, y, z = small[0], small[1], small[2], small[3], small[4]
        if d != unordered_edge(v,z): v,z = z,v
        
        P_prime = [p for p in poly_order if p not in {w, y}]

        # Case 2(a): w and y are NOT in V
        if w not in V and y not in V:
            V_prime = V.copy()
            V_prime.add(z)
            return _solve_guards(P_prime, all_tri_edges, V_prime)
        
        # Case 2(b): Exactly one of w, y is in V
        elif (w in V and y not in V) or (y in V and w not in V):
            V_prime = V - {w, y}
            sub_guards = _solve_guards(P_prime, all_tri_edges, V_prime)
            guard_to_add = unordered_edge(v, w) if w in V else unordered_edge(z, y)
            sub_guards.add(guard_to_add)
            return sub_guards

        # Case 2(c): Both w and y are in V
        else:
            V_prime = V - {w, y}
            V_prime.add(x) # Ensure x is guarded in subproblem
            sub_guards = _solve_guards(P_prime, all_tri_edges, V_prime)
            # Paper's logic: remove a guard on xz and add wx, yz.
            # A simpler, correct implementation is to just add wx and yz.
            sub_guards.add(unordered_edge(v, w))
            sub_guards.add(unordered_edge(z, y))
            return sub_guards
            
    return set() # Should not be reached

def mobile_edge_guards(d, listOfTriangles, V=None, debug=False):

    global DEBUG
    DEBUG = debug

    poly_order = polygon_boundary_coords(d)
    all_tri_edges = build_triangulation_edges(listOfTriangles)
    
    # Convert V from points/tuples to a canonical coordinate tuple set
    V_coords = set()
    if V:
        for v_item in V:
            V_coords.add(_coords_of_point(v_item))

    if DEBUG:
        print(f"Starting mobile guard solver for polygon with {len(poly_order)} vertices.")
        if V_coords:
            print(f"Required to cover {len(V_coords)} vertices in V.")

    guards = _solve_guards(poly_order, all_tri_edges, V_coords)
    
    if DEBUG:
        print(f"\nFinished. Found {len(guards)} guards.")
        print(guards)
        
    return guards, poly_order
    
def plot_mobile_guards(d, guards, nowDraw, newline):
    toDraw = []
    # polygon boundary
    for e in d.getEdges():
        p1 = list(e.origin.coords); p2 = list(e.getTwin().origin.coords)
        X, Y = newline(p1, p2)
        toDraw.append([X, Y, 'k-'])
    # guards (dashed + midpoint)
    for g in guards:
        p1, p2 = g
        X, Y = newline(list(p1), list(p2))
        toDraw.append([X, Y, 'b--'])
        mx = (p1[0] + p2[0]) / 2.0
        my = (p1[1] + p2[1]) / 2.0
        toDraw.append([[mx], [my], 'ro'])
    nowDraw(toDraw)