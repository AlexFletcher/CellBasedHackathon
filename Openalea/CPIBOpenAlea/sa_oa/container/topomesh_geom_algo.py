# -*- python -*-
# -*- coding: utf-8 -*-
#
#       Topomesh : container package
#
#       Copyright or  or Copr. 2006 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       VPlants WebSite : https://gforge.inria.fr/projects/vplants/
#

__doc__="""
This module provide a simple pure python implementation
for a some topomesh algorithms that require the position of points in space
"""

__license__= "Cecill-C"
__revision__=" $Id: topomesh_geom_algo.py 14917 2013-09-27 12:28:55Z pradal $ "

from math import sqrt
from topomesh import Topomesh
from topomesh_algo import topo_triangulate_polygon,ordered_pids,flip_edge

try :
    from numpy import subtract,cross,dot
    from numpy.linalg import norm
    __all__ = ["triangle_quality","flip_necessary",
               "triangulate_polygon","triangulate_face"]
except ImportError :
    __all__ = []


###########################################################
#
#       mesh edition
#
###########################################################
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
    
    :Returns Type: float
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

def flip_necessary (mesh, eid, pos) :
    """Test wether flipping the edge gain something
    in terms of triangles quality.
    
    .. warning:: mesh must be planar with triangle faces only
    
    :Parameters:
     - `mesh` (Topomesh)
     - `eid` (eid) - id of edge to test
     - `pos` (dict of (pid|Vector) ) - position of points
    
    :Returns Type: bool
    """
    #test wether edge is between two triangles
    if mesh.nb_regions(1,eid) != 2 :
        return False
    
    #find points
    lpids = set()
    for fid in mesh.regions(1,eid) :
        lpids.update(mesh.borders(2,fid,2) )
    
    pt1,pt2 = (pos[pid] for pid in mesh.borders(1,eid) )
    pta,ptb = (pos[pid] for pid in (lpids - set(mesh.borders(1,eid) ) ) )
    
    #test wether flipped edge is inside the quadrangle
    if dot(cross(subtract(pt2,pt1),subtract(pta,pt1) ),
           cross(subtract(pt2,pt1),subtract(ptb,pt1) ) ) > 0 :
        return False
    
    if dot(cross(subtract(ptb,pta),subtract(pt1,pta) ),
           cross(subtract(ptb,pta),subtract(pt2,pta) ) ) > 0 :
        return False
    
    #test the quality of triangles
    cur_shape_qual = triangle_quality(pt1,pt2,pta) \
                   + triangle_quality(pt1,pt2,ptb)
    flp_shape_qual = triangle_quality(pta,ptb,pt1) \
                   + triangle_quality(pta,ptb,pt2)
    
    return flp_shape_qual < cur_shape_qual

###############################################
#
#        triangulation
#
###############################################
def triangulate_polygon (pids, pos) :
    """Sort pids in tuples of 3 points
    
    .. seealso:: Compared to `topo_triangulate_polygon`, this function will use
                 geometrical informations to produce nice triangles
    :Parameters:
     - `pids` (list of pid) - ordered list of points that form a closed polygon
     - `pos` (dict of (pid|Vector) ) - position of points in the mesh.
    
    :Returns Type: list of (pid,pid,pid)
    """
    triangles = topo_triangulate_polygon(pids)
    if len(triangles) == 1 :
        return triangles
    
    #change triangles to tends towards equilateral one
    #local triangulated mesh
    mesh = Topomesh(2)
    for pid in pids :
        mesh.add_wisp(0,pid)
    edge = {}
    for pid1,pid2,pid3 in triangles :
        tid = mesh.add_wisp(2)
        for pida,pidb in [(pid1,pid2),(pid2,pid3),(pid3,pid1)] :
            key = (min(pida,pidb),max(pida,pidb) )
            try :
                eid = edge[key]
            except KeyError :
                eid = mesh.add_wisp(1)
                mesh.link(1,eid,pida)
                mesh.link(1,eid,pidb)
                edge[key] = eid
            mesh.link(2,tid,eid)
    
    #flip edges
    test = True
    while test :
        test = False
        for eid in tuple(mesh.wisps(1) ) :
            if flip_necessary(mesh,eid,pos) :
                flip_edge(mesh,eid)
                test = True
    
    return [tuple(mesh.borders(2,fid,2) ) for fid in mesh.wisps(2)]

def triangulate_face (mesh, fid, pos) :
    """Triangulate a face of a mesh
    
    Return a list of triangles for this face.
    
    .. seealso:: A pure topological implementation of this function is defined
                 in `topo_triangulate_face`
    
    :Parameters:
     - `mesh` (:class:`sa_oa.container.Topomesh`)
     - `fid` (fid) - id of the face to triangulate
     - `pos` (dict of (pid|Vector) ) - position of points in the mesh.
    
    :Returns Type: list of (pid,pid,pid)
    """
    pids = ordered_pids(mesh,fid)
    return triangulate_polygon(pids,pos)






