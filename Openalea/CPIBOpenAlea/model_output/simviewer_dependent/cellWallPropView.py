

"""
Class to view properties associated with a (cell, wall) pair.
[ e.g. PINs ]
Based on  ScalarPropView (JC)
"""

from sa_vp.plantgl.algo import GLRenderer
from sa_vp.plantgl.scenegraph import Material,Scene,Shape
from sa_oa.pglviewer import SceneView
from sa_vp.plantgl.math import Vector2,Vector3, dot, norm
from sa_vp.plantgl.scenegraph import Scene, Shape, Material,\
              Translated, Scaled, TriangleSet
from sa_oa.tissueshape import centroid
from sa_oa.container import ordered_pids
from model_utils.geom import area

def mat_map_func (prop, cmap, default_mat) :
    def mat_map (wid) :
        try :
            dat = prop[wid]
            if dat is None :
                return None
            col = cmap(dat)
            if col is None :
                return None
            else :
                mat =  Material(col.i3tuple())
            return mat
        except KeyError :
            return default_mat
	
    return mat_map

def draw(mesh, graph, wall, pos, mat_map, shrink, width, offset):

    
    pos3d = dict( (pid,Vector3(*(float(v) for v in vec) ) ) \
                    for pid,vec in pos.iteritems() )

    offset_vec = Vector3(0, 0, offset)
    
    k = Vector3(0, 0, 1)
    
    sc = Scene()
    for cid in mesh.wisps(2):
        for wid in mesh.borders(2, cid):
            w_idx = list(mesh.borders(1, wid))
            w_pts = [pos3d[pid] for pid in mesh.borders(1, wid)]
            mat = mat_map((cid, wid))
            if mat is not None:
                bary = centroid(mesh, pos, 2, cid)

                o_idx = ordered_pids(mesh, cid)
                if area([pos[pid] for pid in o_idx]) < 0:
                    o_idx.reverse()
                i0=o_idx.index(w_idx[0])
                i1=o_idx.index(w_idx[1])
                if (i1-i0)%len(o_idx)!=1:
                    w_pts.reverse()
                    i0, i1 = i1, i0
                
                i2=(i1+1)%len(o_idx)
                im1=(i0-1)%len(o_idx)
            
                d2 = pos3d[o_idx[i2]] - pos3d[o_idx[i1]]
                d2 = d2 / norm(d2)
                d1 = pos3d[o_idx[i1]] - pos3d[o_idx[i0]]
                d1 = d1 / norm(d1)
                d0 = pos3d[o_idx[i0]] - pos3d[o_idx[im1]]
                d0 = d0 / norm(d0)
                
                n1 = (d1 + d2) ^ k
                n1 = n1 / norm(n1)
                n0 = (d0 + d1) ^ k
                n0 = n0 / norm(n0)

                normal = (d1 ^ k) 

                pts = [w_pts[0], 
                       w_pts[1], 
                       w_pts[1]-width*n1/dot(normal, n1),
                       w_pts[0]-width*n0/dot(normal, n0)]

#            N = (pts[1] - pts[0]) ^ (pts[2] - pts[1])
#            if N.z < 0:
#                pts.reverse()
                geom = TriangleSet(pts, [(0, 1, 2), (0, 2, 3)])  


                scale = 1. - shrink / 100.
        
                if shrink > 0: 
                    bary = Vector3(*(float(v) for v in bary) )
                    scale = 1. - shrink / 100.
                    geom = Translated(bary,
                                      Scaled((scale, scale, scale),
                                             Translated(-bary, geom)))
                if offset !=0:
                    geom = Translated(offset_vec, geom)
            
                shp = Shape(geom, mat)
               # shp.id = (int(cid), int(wid))
                sc.add(shp)
    
    return sc



class CellWallPropView2D(SceneView):

    def __init__(self, mesh, graph, wall, pos, 
                 prop, shrink, offset, width, cmap):
        SceneView.__init__(self)
        self.idmode = GLRenderer.ShapeId
        self.set_alpha_threshold(0.1)
        self.set_name("edge_prop")
        self._mesh = mesh
        self._graph = graph
        self._wall = wall
        self._pos = pos
        self._prop = prop
        self._shrink = shrink
        self._offset = offset
        self._width = width
        self._cmap = cmap
        self._default_mat = Material()
        self._default_mat.transparency = 0.5

    def redraw(self, send_signal=True):
        sc = draw(self._mesh,
                  self._graph,
                  self._wall,
                  self._pos, 
                  mat_map_func(self._prop, self._cmap, self._default_mat),
                  self._shrink,
                  self._width,
                  self._offset)
        self.clear(False)
        self.merge(sc, send_signal)
