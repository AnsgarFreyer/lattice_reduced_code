from sage.geometry.polyhedron.backend_cdd import Polyhedron_RDF_cdd
from sage.misc.prandom import randrange

def diff_polar(P):
    #
    # returns (P-P)*, that is, the polar of the difference body of the polytope P
    #
    return (P + (-P)).polar()

def interior_integral_points(P):
    #
    # returns the integral points in the interior of the polytope P
    #
    L = []
    for v in P.integral_points():
        if P.interior_contains(v):
            L.append(v)
    return L

def nonzero_integral_points(P):
    #
    # returns the integral points in the interior of the polytope P
    #
    L = []
    d = P.dim()
    o = vector([0 for i in range(d)])
    for v in P.integral_points():
        if v != o:
            L.append(v)
    return L

def interior_integral_points_count(P):
    #
    # returns the number of integral points in the interior of the polytope P
    #
    return len(interior_integral_points(P))

def poly_norm(P,x):
    #
    # determines the norm of the vector x with respect to the distance function induced by the polytope P
    # assumes that the origin is interior to P
    #
    d = P.dim()
    o = vector([0 for i in range(d)])
    if not(P.interior_contains(o)):
        return "The origin is not contained in the interior of the polytope."
    facts = []
    for ieq in P.inequalities_list():
        q = (vector(x)*vector(ieq[1:d+1])) / (-ieq[0])
        facts.append(q)
    return max(facts)

def successive_minimum(P):
    #
    # determines the first successive minimum of the polytope P
    # assumes that the origin is interior to P
    #
    d = P.dim()
    o = vector([0 for i in range(d)])
    if not(P.interior_contains(o)):
        return "The origin is not contained in the interior of the polytope."

    s = 1
    while (s*P).integral_points_count() == 1:
        s += 1
    depths = []
    for x in (s*P).integral_points():
        if vector(x) != o:
            depths.append(poly_norm(P,x))
    return min(depths)

def lattice_width(P):
    #
    # determines the lattice-width of the polytope P
    #
    return successive_minimum(diff_polar(P))

def lattice_directions(P):
    # determines the directions where lattice width of the polytope P are achieved
    return nonzero_integral_points(successive_minimum(diff_polar(P))*diff_polar(P))



def is_reduced(P):
    #
    # returns True if the polytope P is lattice reduced; False otherwise
    #
    maximisers = []
    V = []
    W = []
    D = []
    for v in P.vertices():
        V.append(vector(v))
        W.append(tuple(v))
    for f in lattice_directions(P):
        D.append(vector(f))
    for l in D:
        single_max= false
        max = -oo
        for w in V:
            if l*w == max:
                single_max = false
            elif l*w > max:
                max = l*w
                extreme = tuple(w)
                single_max = true
        if single_max == true:
            if extreme not in maximisers:
                maximisers.append(extreme)
    if set(maximisers) == set(W):
        return True
    else:
        return False

def lattice_diameter(P):
    #
    # determines the lattice-diameter of the full dimensional polytope P
    #
    return successive_minimum(P+(-P))


def diameter_directions(P):
    # determines the directions where lattice diameter of the polytope P are achieved
    return nonzero_integral_points(lattice_diameter(P)*(P+(-P)))

def is_complete_simplex(S):
    #
    # returns True if the simplex S is complete; False otherwise
    # does NOT work for polytopes other than simplices.
    #
    V = [v.vector() for v in S.vertices()]
    TCs = [Polyhedron(rays = [v-w for v in V], base_ring=AA) for w in V]
    directions = diameter_directions(S)
    # print("Testing vertices...")
    for K in TCs:
        vertex_is_complete = false
        for v in directions:
            if K.relative_interior_contains(v):
                vertex_is_complete = true
        if not vertex_is_complete:
            return false
    return true
