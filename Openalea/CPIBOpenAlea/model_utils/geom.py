
"""Geometrical utility functions"""

import numpy as np
from numpy.linalg import eigh
from utils import pairs, sgn
from math import sqrt



def ear_clip(poly):
    """ 
    Ear clipping triangulation of a simple polygon
    Track of reflex angles and those corners known not
    to be ears to reduce computational complexity.

    Args:
        poly: list of the vertices of a polygon in counterclockwise order

    Returns:
        list. list of tuples of the three indices for each triangle

    """
    tris=[]
    N = len(poly)
    # set of indices of vertices (in poly) which are reflex
    reflex = set()
    # List of vertex indices; elements will be removed from v
    # when that vertex is at the tip of a clipped ear.
    v = list(range(N))
    # Check that polygon is given in counterclockwise order.
    # If not, reverse the order of the indices in v
    if area(poly)<0:
        v.reverse()
    for i in range(N):
        # Test all vertices to see if reflex
        x = poly[v[(i-1)%N]]
        y = poly[v[i]]
        z = poly[v[(i+1)%N]]
        if ccw(x, y, z) < 0:
            reflex.add(v[i])
    if len(reflex) == 0:
        # Convex polygon - special case
        return [ (0, i, i+1) for i in range(1,N-1) ]
    # Now loop over polygon looking for ears
    not_ears = set()
    while len(v) > 3 :
#        print v
        M = len(v)
        # Loop over all vertices in v. Use an index, as we need to
        # be able to access the previous and next (remaining) vertices.
        # Would be easier if v was a doubly-linked list
        for i in range(M):
            if v[i] not in reflex and v[i] not in not_ears:
                is_ear = True
                idx = v[i]
                prev = v[(i-1)%M]
                next = v[(i+1)%M]
                # a, b, c are the positions of the vertices of
                # the triangle with tip v[i]
                a=poly[prev]
                b=poly[idx]
                c=poly[next]
                # check to see if any reflex vertices lie in a, b, c
                for j in reflex:
                    if j!=prev and j!=next:
                        if pt_in_tri((a,b,c), poly[j]):
                            # mark as not ear, and add to set of 
                            # vertices not to test until neighbour changed
                            is_ear = False
                            not_ears.add(idx)
                            break
                # No convex vertices found inside a, b, c -> ear
                if is_ear:
                    # add triangle to output list
                    tris.append((prev, idx, next))
                    # update convexity of prev, next 
                    if prev in reflex:        
                        x = poly[v[(i-2)%M]]
                        if not ccw(x,a,c)<0:
                            reflex.remove(prev)
                    if next in reflex:
                        z = poly[v[(i+2)%M]]
                        if not ccw(a,c,z)<0:
                            reflex.remove(next)
                    # remove prev and next from not_ears, as it's
                    # now possible for them to be ear tips.
                    not_ears.discard(prev)
                    not_ears.discard(next)
                    # remove v from the working polygon
                    del v[i]
                    break
    tris.append(tuple(v))
    return tris

def triangulate_polygon(pts):
    """ 
    Ear clipping triangulation of a simple polygon

    Args:
        poly: list of the vertices of a polygon in counterclockwise order

    Returns:
        list. list of tuples of the three indices for each triangle
    """

    v=range(len(pts))
    if area(pts)<0:
        v.reverse()
    tris=[]
    is_convex = True
    N = len(v)
    for i in range(N):
        a=pts[(i-1)%N]
        b=pts[i]
        c=pts[(i+1)%N]
        if ccw(a,b,c)<0:
            is_convex = False
    if is_convex:
        return [ (0, i, i+1) for i in range(1,N-1) ]

    while(len(v)>2):
        N=len(v)
        for i in range(N-1, -1, -1):
            a=pts[v[(i-1)%N]]
            b=pts[v[i]]
            c=pts[v[(i+1)%N]]
            if ccw(a,b,c)>0:
                is_ear=True
                for j in range(2, N-1):
                    p=pts[v[(i+j)%N]]
                    if pt_in_tri((a, b, c), p):
                        is_ear=False
                        break
                if is_ear:
                    tris.append((v[(i-1)%N], v[i%N], v[(i+1)%N]))
                    del v[i]
                    break
    return tris
    

def pt_in_polygon(pts, x):
    """    
    Test whether a point lies within a polygon. Adapted from
    http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

    int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy)
    {
    int i, j, c = 0;
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
         c = !c;
    }
    return c;
    }

    Args:
        pts: list of the vertices of a polygon in counterclockwise order
        x: Point

    Returns:
        bool. Whether x is within polygon (unsure about corner cases).

    """
    c=False
    for pj, pi in pairs(pts):
        if (((pj[1]>x[1]) != (pi[1]>x[1])) and
            (x[0]< ((pj[0]-pi[0])*(x[1]-pi[1])/(pj[1]-pi[1]) + pi[0]))):
            c=not c
    return c

def pt_in_tri(tri, x):
    """
    Test whether a point lies within a triangle
    Args:
        tri: list of vertex positions (tuple/Vector2/numpy.ndarray) 
             in counterclockwise order
        x: Point

    Returns:
        bool. Whether x is within triangle (unsure about corner cases).
    """

    s=[ccw(tri[0],tri[1],x), ccw(tri[1], tri[2], x), ccw(tri[2], tri[0], x)]
 #   print 'TRI', s, tri, x
    return not ((1 in s) and (-1 in s))

def dist2(p0, p1):
    """
    Squared istance between points with positions p1 and p2
    
    Args:
        p1, p2: points (tuple/Vector2/list/numpy.ndarray)
    
    Returns:
        float. Square of distance from p1 to p2
    """
    dx=p1[0]-p0[0]
    dy=p1[1]-p0[1]
    return dx*dx+dy*dy

def dist(p0, p1):
    """
    Distance between points with positions p1 and p2
    
    Args:
        p1, p2: points (tuple/Vector2/list/numpy.ndarray)
    
    Returns:
        float. Distance from p1 to p2
    """
    dx=p1[0]-p0[0]
    dy=p1[1]-p0[1]
    return sqrt(dx*dx+dy*dy)

def dist_segment_pt(v, w, p):
    """ 
    Calculates the distance of the point p from the line segment
    with ends v and w.
    
    Args:
        v, w - positions of segment endpoint (tuple / Vector2)
        p - position of point

    Returns:
        float. distance from p to nearest point on v-w
    """

    l2=dist2(v,w)
    t = ((p[0] - v[0])*(w[0] - v[0])+(p[1] - v[1])*(w[1] - v[1]))/l2
    if t < 0.0: 
        return dist(p, v)
    elif t > 1.0:
        return dist(p, w)
    projection=(v[0]+t*(w[0]-v[0]), v[1]+t*(w[1]-v[1]))
    return dist(p, projection);



def ccw(a, b, c) :
    """ Test whether three points are in a counter
        clockwise order in the plane

        Args:
            a, b, c: Positions of the three points (tuple / numpy.ndarray)

        Returns:    
            int. Order of points in the plane:: 
              1  --  counter-clockwise 
              0  --  colinear 
              -1  --  clockwise

        It is probably unwise to use this as a test 
        for colinearity, owing to the inexactness of floating
        point arithmetic.
    """ 
    # test for some of the degenerate cases (two point coninciding)
    # ignore collinear but not coincident points
    if np.all(a==b) or np.all(a==c) or np.all(b==c) :
        return 0
    # rule from 
    # http://compgeom.cs.uiuc.edu/~jeffe/teaching/
    #        algorithms/notes/xn-convexhull.pdf
    else :
        s=(b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0])
        if abs(s)<1e-10:
            return 0
        return cmp(s, 0.0)

def area(poly):
    """ 
    Calculates the area of a polygon
    
    Args:
        poly: list of the vertices of a polygon in counterclockwise order

    Returns:
        float. Area of the polygon
    """
    return sum((0.5*(pp[0][0]*pp[1][1]-pp[1][0]*pp[0][1]) 
                for pp in pairs(poly)))

     
def centroid(poly):
    """ 
    Calculates the centroid of a polygon, defined by
    a list of points

    Args:
        poly: list of the vertices of a polygon in counterclockwise order

    Returns:
        numpy.ndarray. Centroid of the polygon
    
    >>> centroid([(0,0),(1,0),(1,1),(0,1)])
    array([ 0.5,  0.5])
    """
    A = area(poly)
    cx = sum(((pp[0][0]+pp[1][0])*(pp[0][0]*pp[1][1]-pp[1][0]*pp[0][1]) 
            for pp in pairs(poly)))/(6.0*A)
    cy = sum(((pp[0][1]+pp[1][1])*(pp[0][0]*pp[1][1]-pp[1][0]*pp[0][1]) 
            for pp in pairs(poly)))/(6.0*A)
    return (cx, cy)

def det(a, b):
    """
    Calculates dot(a, rot(b)) for two vectors a, b, or alternatively
    the determinant of the 2x2 matrix with columns [a, b]

    Args:
        a, b: The two vectors (tuples / numpy.ndarray)

    Returns:
        float. a[0]*b[1]-a[1]*b[0]

    """
    return a[0]*b[1]-a[1]*b[0]

def moments(poly):
    """ 
    Calculates the first moments of area of a polygon

    Args:
        poly:  List of positions of polygon vertices, in ccw order

    Returns:
        float, float, float. The first moments A, B and C.
    """
    A = 0
    B = 0
    C = 0
    centre = centroid(poly)
    for s, e in pairs(poly):
        vs = (s[0]-centre[0], s[1]-centre[1])
        ve = (e[0]-centre[0], e[1]-centre[1])
        a = 0.5*det(vs, ve)
        C += (1.0/6.0)*a*(vs[0]*vs[0]+vs[0]*ve[0]+ve[0]*ve[0])
        A += (1.0/6.0)*a*(vs[1]*vs[1]+vs[1]*ve[1]+ve[1]*ve[1])
        B -= (1.0/12.0)*a*(2*vs[0]*vs[1]+vs[0]*ve[1]+vs[1]*ve[0]+2*ve[0]*ve[1])
    return A, B, C

def long_axis(poly):
    """
    Calculates the long axis of a polygon, described
    by a list of points

    Args:
        poly:  List of positions of polygon vertices, in ccw order

    Returns:
        numpy.ndarray. direction of the long axis
    
    """
    A, B, C = moments(poly)
    v = eigh([[A, B], [B, C]])[1][0]
    return v
