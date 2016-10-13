from sa_oa.tissueview import ScalarPropView
from model_output.edgePropView import EdgePropView2D
from sa_oa.tissueview.prop_display import draw_scalar_prop
from sa_vp.plantgl.scenegraph import Material

default_mat = Material()
default_mat.transparency = 0.5

from model_utils.db_utilities import get_mesh,get_graph

#can only use with fixed geometry

class CustomScalarPropView_edge(EdgePropView2D):

	def __init__(self, db, prop_name, cmap) :
		mesh=get_mesh(db)
		graph=get_graph(db)
		wall=db.get_property('wall')
		pos=db.get_property('position')
		deg=mesh.degree()
		self._prop_name=prop_name
		prop=db.get_property(prop_name)
		shrink=0
		offset=0
		width=5
		#ScalarPropView.__init__(self, mesh, pos, deg, prop, shrink, cmap)
		EdgePropView2D.__init__(self,mesh, graph, wall, pos, prop, shrink, offset, width, cmap)
		self.sc = draw_scalar_prop(self._mesh,
			                      self._pos,
			                      2,
			                      self._prop,
			                      self._shrink,
			                      self._cmap,
			                      default_mat,
			                      'topo')
	def cache_geometry_subset (self,type_map,active_types) :


		self._cache_geometry={}
		for shp in self.sc:
			if type_map[shp.id] in active_types:
				self._cache_geometry[shp.id]=shp.geometry


	def change_property (self,prop,propname):
		self._prop=prop
		self._prop_name=propname


	def renormalise(self):
		new_max=max(list(self._prop.itervalues()))
		def_max=0.01
		from sa_vp.plantgl.ext.color import JetMap
		if new_max>0:
			self._cmap=JetMap(0,new_max, outside_values=True)
		else:
			self._cmap=JetMap(0,def_max, outside_values=True)
