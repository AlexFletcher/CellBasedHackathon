# -*- python -*-
#
#       shapes: function used to create tissues
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
This module defines a function to extrude a mesh
"""

__license__= "Cecill-C"
__revision__=" $Id: $ "

from math import radians,degrees
from sa_oa.plantgl.math import Vector2,norm
from sa_oa.physics.mechanics import LinearSpring2D,CircularSpring2D,ForwardEuler2D

K = 1.
L = 1.
dt = 0.3

def roundify_2D (mesh, pos, boundary_func) :
	"""
	change edges orientation to obtain more spherical cells
	empirical rule is that this algo is working correctly
	if each edge is around 1 in length
	mesh: a topomesh
	pos: a dict of {pid:Vector2)
	fixed_points : a list of points that must not move
	return : None, modify pos in place
	"""
	springs = []
	for cid in mesh.wisps(2) :
		edges = set(mesh.borders(2,cid))
		points = list(mesh.borders(2,cid,2))
		bary = sum((pos[pid] for pid in points),Vector2())/len(points)
		tot_length = 0
		for eid in edges :
			pid1,pid2 = mesh.borders(1,eid)
			tot_length += norm(pos[pid1] - pos[pid2])
		ref_length = tot_length/len(edges)
		for eid in edges :
			pid1,pid2 = mesh.borders(1,eid)
			spring = LinearSpring2D(pid1,pid2,K,ref_length)
			springs.append(spring)
		ref_angle = radians(360./len(edges))
		for pid in points :
			local_edges = set(mesh.regions(0,pid)) & edges
			if len(local_edges) == 2 :
				eid1,eid2 = local_edges
				pid1,pid2 = ( set(mesh.borders(1,eid1))|set(mesh.borders(1,eid2)) ) - set([pid])
				if (pos[pid]-bary)^(pos[pid1]-bary)>0 :
					spring = CircularSpring2D(pid,pid1,pid2,L,ref_angle)
				else :
					spring = CircularSpring2D(pid,pid2,pid1,L,ref_angle)
				springs.append(spring)
	#resolution
	weight = dict( (pid,1.) for pid in mesh.wisps(0) )
	algo = ForwardEuler2D(weight,springs,boundary_func)
	algo.deform(pos,dt,300)
	return springs
