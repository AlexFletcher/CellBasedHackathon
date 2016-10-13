
# Based on ScalarPropView (JC)


from sa_vp.plantgl.algo import GLRenderer
from sa_vp.plantgl.scenegraph import Material, Scene, Shape
from sa_oa.pglviewer import SceneView
from sa_vp.plantgl.math import Vector2,Vector3, dot, norm
from sa_vp.plantgl.scenegraph import Scene, Shape, Material,\
              Translated, Scaled, TriangleSet
from sa_oa.tissueshape import centroid
from sa_oa.container import ordered_pids
from model_utils.geom import area
from PyQt4.QtCore import Qt, SIGNAL
from time import clock
from sa_oa.tissueview.mesh_pgl import face_geom3D


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
#    start= clock()
    
    pos3d = dict( (pid,Vector3(*(float(v) for v in vec) ) ) \
                    for pid,vec in pos.iteritems() )

    offset_vec = Vector3(0, 0, offset)
    
    k = Vector3(0, 0, 1)
    
    sc = Scene()
    if mesh.degree()==2:
	    cell_pts = {}
	    cell_inner_pts = {}
	    for cid in mesh.wisps(mesh.degree()):
		o_idx = ordered_pids(mesh, cid)
		if area([pos[pid] for pid in o_idx]) < 0:
		        o_idx.reverse()
		cell_pts[cid] = o_idx
		N = len(o_idx)
		i_pts = []
		for i in range(N):
		    pi = pos3d[o_idx[i]]
		    ip1 = (i+1)%N
		    im1 = (i-1)%N
		    dp = pos3d[o_idx[ip1]] - pi
		    dp = dp / norm(dp)
		    dm = pi - pos3d[o_idx[im1]]
		    dm = dm / norm(dm)
		    ni = (dp + dm) ^ k
		    ni = ni / norm(ni)
		    normal = (dp ^ k)
		    i_pts.append(pi - width*ni/dot(normal, ni))
		cell_inner_pts[cid] = i_pts

	    for eid in graph.edges():
		wid = wall[eid]
		cid = graph.source(eid)
		w_idx = list(mesh.borders(mesh.degree()-1, wid))
		w_pts = [pos3d[pid] for pid in mesh.borders(mesh.degree()-1, wid)]
		mat = mat_map(eid)
		if mat is not None:

		    o_idx = cell_pts[cid]

		    i0=o_idx.index(w_idx[0])
		    i1=o_idx.index(w_idx[1])
		    if (i1-i0)%len(o_idx)!=1:
		        w_pts.reverse()
		        i0, i1 = i1, i0

		    pi0 = cell_inner_pts[cid][i0]
		    pi1 = cell_inner_pts[cid][i1]
		    

		    pts = [w_pts[0], 
		           w_pts[1], 
		           pi1,
		           pi0 ]

		    geom = TriangleSet(pts, [(0, 1, 2), (0, 2, 3)])  


		    scale = 1. - shrink / 100.
		
		    if shrink > 0: 
		        bary = centroid(mesh, pos, mesh.degree(), cid)
		        bary = Vector3(*(float(v) for v in bary) )
		        scale = 1. - shrink / 100.
		        geom = Translated(bary,
		                          Scaled((scale, scale, scale),
		                                 Translated(-bary, geom)))
		    if offset !=0:
		        geom = Translated(offset_vec, geom)
		    
		    shp = Shape(geom, mat)
		    shp.id = int(eid)
		    sc.add(shp)
    elif mesh.degree()==3:
	    cell_pts = {}
	    cell_inner_pts = {}
	    for cid in mesh.wisps(mesh.degree()-1):
		o_idx = ordered_pids(mesh, cid)
		if area([pos[pid] for pid in o_idx]) < 0:
		        o_idx.reverse()
		cell_pts[cid] = o_idx
		N = len(o_idx)
		i_pts = []
		for i in range(N):
		    pi = pos3d[o_idx[i]]
		    ip1 = (i+1)%N
		    im1 = (i-1)%N
		    dp = pos3d[o_idx[ip1]] - pi
		    dp = dp / norm(dp)
		    dm = pi - pos3d[o_idx[im1]]
		    dm = dm / norm(dm)
		    ni = (dp + dm) ^ k
		    ni = ni / norm(ni)
		    normal = (dp ^ k)
		    i_pts.append(pi - width*ni/dot(normal, ni))
		cell_inner_pts[cid] = o_idx

	    for eid in graph.edges():
		wid = wall[eid]
		
		cid = graph.source(eid)
		w_idx = list(mesh.borders(mesh.degree()-1, wid))
		#w_pts = [pos3d[pid] for pid in mesh.borders(mesh.degree()-1, wid)]
		
		mat = mat_map(eid)
		if mat is not None:
		    """
		    o_idx = cell_pts[cid]

		    i0=o_idx.index(w_idx[0])
		    i1=o_idx.index(w_idx[1])
		    if (i1-i0)%len(o_idx)!=1:
		        w_pts.reverse()
		        i0, i1 = i1, i0

		    pi0 = cell_inner_pts[cid][i0]
		    pi1 = cell_inner_pts[cid][i1]
		    

		    pts = [w_pts[0], 
		           w_pts[1], 
		           pi1,
		           pi0 ]
		    """
		    up = centroid(mesh,pos,2,wid)
		    geom = face_geom3D(mesh,pos,wid,up,'topo')

		    """
		    scale = 1. - shrink / 100.
		
		    if shrink > 0: 
		        bary = centroid(mesh, pos, mesh.degree(), cid)
		        bary = Vector3(*(float(v) for v in bary) )
		        scale = 1. - shrink / 100.
		        geom = Translated(bary,
		                          Scaled((scale, scale, scale),
		                                 Translated(-bary, geom)))
		    if offset !=0:
		        geom = Translated(offset_vec, geom)
		    """
		    shp = Shape(geom, mat)
		    shp.id = int(eid)
		    sc.add(shp)
#    print 'time: ', clock()-start
    return sc



class EdgePropView2D(SceneView):

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
	self._cache_geometry = None
	self._deg=mesh.degree()
    def redraw(self, send_signal=True):

		if self._cache_geometry is None :
			sc = draw(self._mesh,
				  self._graph,
				  self._wall,
				  self._pos, 
				  mat_map_func(self._prop, self._cmap, self._default_mat),
				  self._shrink,
				  self._width,
				  self._offset)
		else :
			sc = Scene()
			mat_map = mat_map_func(self._prop,
			                        self._cmap,
			                        self._default_mat)
			for wid,geom in self._cache_geometry.iteritems() :
				mat = mat_map(wid)
				if mat is not None :
					shp = Shape(geom,mat)
					shp.id = wid
					sc.add(shp)
		
		self.clear(False)
		self.merge(sc,send_signal)



    def cache_geometry (self, cache = True) :
	"""Precompute the geometry.

	Use this function to accelerate
	the redraw of a scene if the
	underlying geometry do not change.

	:Parameters:
	 - `cache` (bool) - if True, compute
	   the actual geometry of a scene and
	   store it. Else, reset any precomputed
	   geometry

	:Returns Type: False
	"""
	if cache :
		sc = draw(self._mesh,
				  self._graph,
				  self._wall,
				  self._pos, 
				  mat_map_func(self._prop, self._cmap, self._default_mat),
				  self._shrink,
				  self._width,
				  self._offset)
	
		self._cache_geometry = dict( (shp.id,shp.geometry) \
			                      for shp in sc)

	else :
		self._cache_geometry = None

        
    def prop_type(self):
        """
        Retrieve the type of value stored in this quantity.
		
        :Return:
        - a string defining the type
        - None if no type is defined
		
        :Returns Type:
        - str
        - None
        """
        try :
            return self._prop.type()
        except AttributeError :
            return 'float'

    def value (self, elmid) :
        """
        Returns the value of the property.
        """
        return self._prop.get(elmid, None)

    def set_value (self, elmid, val, send_signal = True) :
        """
        Set the value of the property.
        """
        if val is None :
            self._prop.pop(elmid, None)
        else :
            self._prop[elmid] = val
        if send_signal :
            self.emit(SIGNAL("set_value"),elmid,val)

    def values (self) :
        """
        Iterator on all (wid,val).
        """
        return self._prop.iteritems()
	
    def set_values (self, items) :
        """
        Change all values.
        """
        self._prop.clear()
        self._prop.update(items)
        self.emit(SIGNAL("set_values"))
	
    def clear_values (self) :
        """
        Clear property
        """
        self._prop.clear()
        self.emit(SIGNAL("clear_values"))

    def fill (self, elmid, val) :
        """
        Fill a zone around the given element
        with the given value.
        """
        self.set_value(elmid, val)
        self.emit(SIGNAL("fill"), elmid, val)
