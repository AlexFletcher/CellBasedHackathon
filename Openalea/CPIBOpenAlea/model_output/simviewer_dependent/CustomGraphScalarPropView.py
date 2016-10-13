from model_output.graph_prop_view import GraphScalarPropView
from sa_oa.tissueview import ScalarPropView
from model_output.graph_prop_display import draw_graph_prop
from sa_vp.plantgl.scenegraph import Material

default_mat = Material()
default_mat.transparency = 0.5

from model_utils.db_utilities import get_graph

#can only use with fixed geometry

class CustomGraphScalarPropView(GraphScalarPropView):

	def __init__(self, db, prop_name, cmap) :
		graph=get_graph(db)
		pos=db.get_property('position')
		self._prop_name=prop_name
		prop=db.get_property(prop_name)
		shrink=0
		deg=3
		GraphScalarPropView.__init__(self, graph, pos, deg, prop, shrink, cmap)
		self.sc = draw_graph_prop (self._graph, self._pos, self._deg, self._prop, self._shrink, self._cmap,
                      None, 'topo')
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
