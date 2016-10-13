# -*- python -*-
#
#       tissueview: function used to display tissue properties
#
#       Copyright 2006 INRIA - CIRAD - INRA  
#
#       File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
# 
#       OpenAlea WebSite : http://sa_oa.gforge.inria.fr
#

__doc__="""
This module defines functions to display a property on a topomesh
"""

__license__= "Cecill-C"
__revision__=" $Id: $ "

from PyQt4.QtCore import Qt,SIGNAL
from sa_vp.plantgl.algo import GLRenderer
from sa_vp.plantgl.scenegraph import Material,Scene,Shape
from sa_oa.pglviewer import SceneView
from model_output.graph_prop_display import (draw_graph_prop,
                          _mat_map_func)

default_mat = Material()
default_mat.transparency = 0.5

class GraphScalarPropView (SceneView) :
	"""View on a scalar prop defined on a mesh.
	"""
	def __init__ (self, graph, pos, deg, prop, shrink, cmap) :
		SceneView.__init__(self)
		self.idmode = GLRenderer.ShapeId
		self.set_alpha_threshold(0.1)
		
		self.set_name("scalar_prop")
		
		self._graph = graph
		self._pos = pos
		self._deg = deg
		self._prop = prop
		self._shrink = shrink
		self._cmap = cmap
		
		self._triangulation_method = 'topo'
		self._cache_geometry = None
	
	def redraw (self, send_signal = True) :
		if self._cache_geometry is None :
			sc = draw_graph_prop (self._graph, self._pos, self._deg, self._prop, self._shrink, self._cmap,
                      None, 'topo')
		else :
			sc = Scene()
			mat_map = _mat_map_func(self._prop,
			                        self._cmap,
			                        default_mat)
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
			sc = draw_scalar_prop(self._graph,
			                      self._pos,
			                      self._deg,
			                      self._prop,
			                      self._shrink,
			                      self._cmap,
			                      default_mat,
			                      self._triangulation_method)
			
			self._cache_geometry = dict( (shp.id,shp.geometry) \
			                              for shp in sc)

		else :
			self._cache_geometry = None


	
	##############################################
	#
	#	interaction
	#
	##############################################
	def prop_type (self) :
		"""Retrieve the type of value
		stored in this quantity.
		
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
			return None
	
	def value (self, elmid) :
		"""Returns the value of the property.
		"""
		return self._prop.get(elmid,None)
	
	def set_value (self, elmid, val, send_signal = True) :
		"""Set the value of the property.
		"""
		if val is None :
			self._prop.pop(elmid,None)
		else :
			self._prop[elmid] = val
		if send_signal :
			self.emit(SIGNAL("set_value"),elmid,val)
	
	def values (self) :
		"""Iterator on all (wid,val).
		"""
		return self._prop.iteritems()
	
	def set_values (self, items) :
		"""Change all values.
		"""
		self._prop.clear()
		self._prop.update(items)
		self.emit(SIGNAL("set_values") )
	
	def clear_values (self) :
		"""Clear property
		"""
		self._prop.clear()
		self.emit(SIGNAL("clear_values") )
	
	def fill (self, elmid, val) :
		"""Fill a zone around the given element
		with the given value.
		"""
		mesh = self._mesh
		deg = self._deg
		prop = self._prop
		old_val = prop.get(elmid,None)
		#temporary
		assert deg > 0
		#
		front = set([elmid])
		while len(front) > 0 :
			wid = front.pop()
			self.set_value(wid,val,False)
			for nid in mesh.border_neighbors(deg,wid) :
				if (self.value(nid) == old_val) \
				   and (self.value(nid) != val) :
					front.add(nid)
		self.emit(SIGNAL("fill"),elmid,val)


