# -*- python -*-
#
#       mesh: function used to deal with mesh geometry
#
#       Copyright 2006 INRIA - CIRAD - INRA  
#
#       File author(s): Jerome Chopard <revesansparole@gmail.com>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
# 
#       OpenAlea WebSite : http://sa_oa.gforge.inria.fr
#

__doc__ = """
This module defines algos to apply on a mesh
"""

__license__ = "Cecill-C"
__revision__ = " $Id: mesh_algo.py 14917 2013-09-27 12:28:55Z pradal $ "

try :
    from numpy import sqrt, subtract, cross, dot
    from numpy.linalg import norm
    
    from topomesh_algo import flip_edge
    
    __all__ = ["center_mesh"
             , "triangle_quality"
             , "is_flip_better"
             , "triangulate"]
except ImportError :
    __all__ = []


def center_mesh (mesh) :
    """modify all position in the mesh to center it
    around its barycenter
    """
    pos = tuple(mesh.position(did) for did in mesh.darts(0) )
    if len(pos) == 0 :
        return
    
    if None in pos :
        raise UserWarning("One position is not defined")
    
    bary = reduce(lambda x, y: x + y, pos) / len(pos)
    
    def center (vec) :
        return vec - bary
    
    mesh.apply_geom_transfo(center)


def triangle_quality (pt1, pt2, pt3) :
    """Construct a quality index for a triangle
    
    The value of the index varie between 1 and infinity,
    1 being an equilateral triangle.
    
    Q = hmax P / (4 sqrt(3) S)
    where :
     - hmax, length of longest edge
     - P, perimeter
     - S, surface
    
    :Parameters:
     - `pt1` (Vector) - corner
     - `pt2` (Vector) - corner
     - `pt3` (Vector) - corner
    
    :Returns: (float) - triangle quality
    """
    e1 = norm(subtract(pt3,pt2) )
    e2 = norm(subtract(pt3,pt1) )
    e3 = norm(subtract(pt2,pt1) )
    
    P = (e1 + e2 + e3) / 2.
    hmax = max(e1,e2,e3)
    S2 = P * (P - e1) * (P - e2) * (P - e3)
    if S2 <= 0 :
        return 1e6
    else :
        return hmax * P / 2. / sqrt(3.) / sqrt(S2)


def is_flip_better (mesh, eid) :
    """Test wether flipping the edge gain something
    in terms of triangles quality.
    
    .. warning:: mesh must be planar with triangle faces only
    
    :Parameters:
     - `mesh` (Mesh)
     - `eid` (eid) - id of edge to test, must be of degree 1
    
    :Returns Type: bool
    """
    #test wether edge is between two triangles
    if mesh.nb_regions(eid) != 2 :
        return False
    
    #find points
    lpids = set()
    for fid in mesh.regions(eid) :
        lpids.update(mesh.borders(fid,2) )

    seg = set(mesh.borders(eid) )
    pt1, pt2 = (mesh.position(pid) for pid in seg)
    pta, ptb = (mesh.position(pid) for pid in (lpids - seg) )
    
    #test wether flipped edge is inside the quadrangle
    if dot(cross(subtract(pt2, pt1), subtract(pta, pt1) ),
           cross(subtract(pt2, pt1), subtract(ptb, pt1) ) ) > 0 :
        return False
    
    if dot(cross(subtract(ptb, pta), subtract(pt1, pta) ),
           cross(subtract(ptb, pta), subtract(pt2, pta) ) ) > 0 :
        return False
    
    #test the quality of triangles
    cur_shape_qual = triangle_quality(pt1, pt2, pta) \
                   + triangle_quality(pt1, pt2, ptb)
    flp_shape_qual = triangle_quality(pta, ptb, pt1) \
                   + triangle_quality(pta, ptb, pt2)
    
    return flp_shape_qual < cur_shape_qual


def triangulate (mesh, fid) :
    """Triangulate the given face of a mesh
    
    :Parameters:
     - `mesh` (Mesh)
     - `fid` (did) - id of face to triangulate. Must be of
                     degree 2
    
    :Warning: Modify mesh in place.
    
    :Note: This method does not introduce new points in the
    mesh.
    """
    #store info
    edges = tuple(mesh.borders(fid) )
    if len(edges) == 3 :#face is already a triangle
        return
    
    if len(edges) < 3 :
        raise UserWarning("face is not geometricaly defined")
    
    cells = tuple(mesh.regions(fid) )
    
    #remove face
    mesh.remove_dart(fid)
    
    #walk along edges
    edge_pool = list(edges)
    front_eid = edge_pool.pop()
    ref_pid, front_pid = mesh.borders(front_eid)
    ref_pos = mesh.position(ref_pid)
    
    new_edges = []
    while len(edge_pool) > 2 :#remaining two edges will
                              #naturally form a triangle
        #find next edge
        next_eid, = (eid for eid in edge_pool \
                     if front_pid in mesh.borders(eid) )
        next_pid, = (pid for pid in mesh.borders(next_eid) \
                     if pid != front_pid)
        
        #triangulate, create an edge between ref_pid and
        #next_pid
        new_eid = mesh.add_dart(1)
        mesh.link(new_eid, ref_pid)
        mesh.link(new_eid, next_pid)
        l = norm(mesh.position(next_pid) - ref_pos)
        new_edges.append( (l, new_eid) )
        
        #create new facet
        fid = mesh.add_dart(2)
        for cid in cells :
            mesh.link(cid, fid)
        
        for eid in (front_eid, next_eid, new_eid) :
            mesh.link(fid, eid)
        
        #increment step
        edge_pool.remove(next_eid)
        front_eid = new_eid
        front_pid = next_pid
    
    #create last facet
    fid = mesh.add_dart(2)
    for cid in cells :
        mesh.link(cid, fid)
    
    for eid in edge_pool + [front_eid] :
        mesh.link(fid, eid)
    
    #improve triangulation quality
    new_edges.sort()
    
    for l, eid in new_edges :
        if is_flip_better(mesh, eid) :
            flip_edge(mesh, eid)





















