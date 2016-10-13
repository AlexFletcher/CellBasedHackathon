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

from numpy.linalg import eig
from sa_vp.plantgl.math import Vector2,Vector3
from sa_vp.plantgl.scenegraph import Scene,Shape,Material,\
              Polyline
from sa_vp.plantgl.ext.color import red,green,blue,black
from sa_oa.container import Graph,Topomesh
from sa_oa.tissueshape import centroid,gcentroid
from sa_oa.tissueview.graph_display import draw_graph


axes_mat = [Material(col.i3tuple() ) for col in (red,green,blue,black)]

def _mat_map_func (prop, cmap, default_mat) :
	def mat_map (wid) :
		try :
			#print 'prop',prop
			dat = prop[wid]
			if dat is None :
				return None
			
			col = cmap(dat)
			if col is None :
				return None
			else :
				mat =  Material(col.i3tuple() )
				mat.transparency = col.transparency / 255.
#				mat.diffuse = 0.5
#				mat.specular = (0,0,0)#tuple(int(v / 3.) for v in col.i3tuple() )
#				mat.emission = (0,0,0)#tuple(int(v / 3.) for v in col.i3tuple() )
				return mat
		except KeyError :
			#print 'KeyError',wid
			return default_mat
	
	return mat_map

def draw_graph_prop (graph, pos, deg, prop, shrink, cmap,
                      default_mat = None, method = 'topo') :
	"""Draw elements of the mesh for which the
	property is defined using the specified color
	map
	"""
	mat_map = _mat_map_func(prop,cmap,default_mat)
	
	return draw_graph (graph, pos, "edge",mat_map, mat_map)



