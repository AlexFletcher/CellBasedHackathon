# -*- python -*-
#
#       tissueshape: function used to deal with tissue geometry
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

__doc__ = """
This module defines a set of functions to access graph geometrical properties
"""

__license__ = "Cecill-C"
__revision__ = " $Id: $ "

from math import pi,atan2
from numpy import zeros,add,subtract,multiply,divide,cross,dot,outer
from numpy.linalg import eig,norm
from sa_oa.container import triangulate_face

__all__ = ["gcentroid","gedge_length"]

####################################################
#
#		geometrical direct properties
#
####################################################
def gcentroid (graph, pos, elm_type, eid) :
	"""Compute centroid of a element of the graph
	
	:Parameters:
	 - `graph` (:class:`sa_oa.container.Graph`)
	 - `pos` (dict of pid|array) - geometrical position of vertices in space
	 - `elm_type` (str) - type of element, either 'vertex' or 'edge'
	 - `eid` (eid) - id of the element to consider
	
	:Returns Type: array
	"""
	if elm_type == "vertex" :
		return pos[eid]
	
	pts = pos[graph.source(eid)]
	ptt = pos[graph.target(eid)]
	return divide(add(pts,ptt),2.)

def gedge_length (graph, pos, eid) :
	"""Compute length of a segment edge
	
	:Parameters:
	 - `graph` (:class:`sa_oa.container.Graph`)
	 - `pos` (dict of pid|array) - geometrical position of vertices in space
	 - `eid` (eid) - id of the edge to consider
	
	:Returns Type: float
	"""
	pid1 = graph.source(eid)
	pid2 = graph.target(eid)
	return norm(subtract(pos[pid2],pos[pid1]) )

