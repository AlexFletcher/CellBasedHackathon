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
This module defines a set of functions to access mesh geometrical properties
"""

__license__ = "Cecill-C"
__revision__ = " $Id: $ "

from math import pi,atan2
from numpy import zeros,add,subtract,multiply,divide,cross,dot,outer
from numpy.linalg import eig,norm
from sa_oa.container import triangulate_face

__all__ = ["centroid","barycenter",
           "edge_length",
           "face_surface_2D","face_surface_3D",
           "cell_volume",
           "edge_normal","face_normal",
           "face_curvature",
           "face_main_axes_2D","face_main_axes_3D",
           "cell_main_axes",
           "find_nearby_wisp",
           "edge_loop_around"]

####################################################
#
#		geometrical direct properties
#
####################################################
def centroid (mesh, pos, degree, wid) :
	"""Compute centroid of a wisp
	
	.. seealso: this function considers only set of points with the same weight.
	            For an exact barycenter see :func:`barycenter`
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `degree` (int) - degree of the considered wisp
	 - `wid` (wid) - id of the wisps to consider
	
	:Returns Type: array
	"""
	if degree == 0 :
		return pos[wid]
	
	pts = [pos[pid] for pid in mesh.borders(degree,wid,degree)]
	return divide(reduce(add,pts),len(pts) )

def barycenter (mesh, pos, degree, wid) :
	"""Compute the barycenter of a wisp
	
	.. seealso: this function compute the exact barycenter
	            for an approximation faster to compute
	            see :func:`centroid`
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `degree` (int) - degree of the considered wisp
	 - `wid` (wid) - id of the wisps to consider
	
	:Returns Type: array
	"""
	print "warning TODO"
	return centroid(mesh,pos,degree,wid)

def edge_length (mesh, pos, eid) :
	"""Compute length of a segment edge
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `eid` (eid) - id of the edge to consider
	
	:Returns Type: float
	"""
	pid1,pid2 = mesh.borders(1,eid)
	return norm(subtract(pos[pid2],pos[pid1]) )

def face_surface_2D (mesh, pos, fid, return_barycenter = False) :
	"""Compute surface of a polygonal convex face
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `fid` (fid) - id of the face to consider
	 - `return_barycenter` (bool) - tells wether the function will return
	                                the barycenter of the face too
	
	:Returns Type: float or (float,array)
	"""
	bary = centroid(mesh,pos,2,fid)
	
	#compute triangle for each edge
	surface = 0.
	for eid in mesh.borders(2,fid) :
		pid1,pid2 = mesh.borders(1,eid)
		surface += abs(cross(subtract(pos[pid1],bary),
		                     subtract(pos[pid2],bary) ) )
	
	#return
	if return_barycenter :
		return surface / 2.,bary
	else :
		return surface / 2.

def face_surface_3D (mesh, pos, fid, return_barycenter = False) :
	"""Compute surface of a polygonal convex face
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `fid` (fid) - id of the face to consider
	 - `return_barycenter` (bool) - tells wether the function will return 
	                                the barycenter of the face too
	
	:Returns Type: float or (float,array)
	"""
	bary = centroid(mesh,pos,2,fid)
	
	#compute triangle for each edge
	surface = 0.
	for eid in mesh.borders(2,fid) :
		pid1,pid2 = mesh.borders(1,eid)
		normal = cross(subtract(pos[pid1],bary),subtract(pos[pid2],bary) )
		surface += norm(normal)
	
	#return
	if return_barycenter :
		return surface / 2.,bary
	else :
		return surface / 2.

def cell_volume (mesh, pos, cid, return_barycenter = False) :
	"""Compute the volume of a given polyhedral cell
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `cid` (cid) - id of the cell to consider
	 - `return_barycenter` (bool) - tells wether the function will return 
	                                the barycenter of the cell too
	
	:Returns Type: float or (float,array)
	"""
	cell_center = centroid(mesh,pos,3,cid)
	
	#compute the pyramid for each face
	V = 0.
	for fid in mesh.borders(3,cid) :
		S,face_center = face_surface_3D(mesh,pos,fid,True)
		V += norm(subtract(face_center,cell_center) ) * S / 3.
	
	#return
	if return_barycenter :
		return V,cell_center
	else :
		return V

def edge_normal (mesh, pos, eid, ref_point, return_barycenter = False) :
	"""Compute the outer normal of an edge in a 2D space
	
	The normal of an edge is perpendicular to
	the segment that represent this edge
	and its length is equal to the surface of the polygon
	
	The returned normal is leading away from the given ref_point
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `eid` (eid) - id of the edge to consider
	 - `ref_point` (array) - point located below
	     the edge (typically the center of the cell)
	 - `return_barycenter` (bool) - tells wether the function will return 
	                                the barycenter of the cell too
	
	:Returns Type: array or (array,array)
	"""
	pt1,pt2 = (pos[pid] for pid in mesh.borders(1,eid) )
	
	#barycenter of the edge
	bary = divide(add(pt1,pt2),2.)
	
	#outside direction
	outside_dir = subtract(bary,ref_point)
	
	#normal
	seg_dir = subtract(pt2,pt1)
	normal = (-seg_dir[1],seg_dir[0])
	
	if dot(normal,(outside_dir[0],outside_dir[1]) ) < 0 :
		normal = multiply(normal,-1.)
	
	#return
	if return_barycenter :
		return normal,bary
	else :
		return normal

def triangle_normal (pt1, pt2, pt3, ref_point) :
	"""Compute the outer normal of a triangle
	
	The length of the normal is equal to the surface of the triangle
	
	:Returns Type: array
	"""
	bary = divide(add(add(pt1,pt2),pt3),3.)
	
	normal = divide(cross(subtract(pt2,pt1),subtract(pt3,pt1) ),2.)
	
	if dot(normal,subtract(bary,ref_point) ) < 0 :
		normal = multiply(normal,-1.)
	
	return normal,bary

def face_normal (mesh, pos, fid, ref_point, return_barycenter = False) :
	"""Compute the outer normal of a face in a 3D space
	
	The normal of a polygonal face is perpendicular to
	the subdivision of the polygon in triangles
	and its length is equal to the surface of the polygon
	
	The returned normal is leading away from the given ref_point
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `fid` (fid) - id of the face to consider
	 - `ref_point` (array) - point located below
	     the face (typically the center of the cell)
	 - `return_barycenter` (bool) - tells wether the function will return 
	                                the barycenter of the cell too
	
	:Returns Type: array or (array,array)
	"""
	pids = tuple(mesh.borders(2,fid,2) )
	if len(pids) == 3 :
		#triangle face
		normal,bary = triangle_normal(*tuple(pos[pid] for pid in pids) )
	else :
		#polygonal face
		
		#barycenter of the face
		bary = centroid(mesh,pos,2,fid)
		
		#normal for each triangle
		normal_list = []
		for eid in mesh.borders(2,fid) :
			pt1,pt2 = (pos[pid] for pid in mesh.borders(1,eid) )
			normal_list.append(triangle_normal(bary,pt1,pt2,ref_point)[0])
		#normal
		if len(normal_list) == 0 :
			raise UserWarning("face not geometrically defined")
		normal = divide(reduce(add,normal_list),2.)
	
	#return
	if return_barycenter :
		return normal,bary
	else :
		return normal

def face_curvature (mesh, pos, fid) :
	"""Compute a measure of the curvature of the face
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `fid` (fid) - id of the face to consider
	
	:Returns Type: float
	"""
	ref_point = (0,0,0)
	N_list = []
	for pid1,pid2,pid3 in triangulate_face(mesh,fid) :
		N,bary = triangle_normal(pos[pid1],pos[pid2],pos[pid3],ref_point)
		N_list.append(divide(N,norm(N) ) )
	
	Nref = N_list[0]
	
	return sum(abs(dot(N,Nref) ) for N in N_list) / len(N_list)


####################################################
#
#		main axes
#
####################################################
def face_main_axes_2D (mesh, pos, fid) :
	"""Compute main axis of a face in 2D
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `fid` (fid) - id of the face to consider
	
	:Returns: barycenter and axes
	
	:Returns Type: array,array,array
	"""
	#compute local weighted point according to
	#barycenter and connected edge length
	L = dict( (eid,edge_length(mesh,pos,eid) ) for eid in mesh.borders(2,fid) )
	tot = sum(L.itervalues() )
	L = dict( (eid,val / tot) for eid,val in L.iteritems() )
	
	bary = centroid(mesh,pos,2,fid)
	
	pts = [multiply(subtract(pos[pid],bary),
	                sum(L.get(eid,0.) for eid in mesh.regions(0,pid) ) ) \
	       for pid in mesh.borders(2,fid,2)]
	
	#compute inertia
	I = divide(reduce(add,[outer(pt,pt) for pt in pts]),float(len(pts) ) )
	
	#compute eigen values
	w,v = eig(I)
	
	vecs = [(w[i],multiply(v[:,i],w[i]) ) for i in (0,1)]
	vecs.sort(reverse=True)
	
	#return
	return bary,vecs[0][1],vecs[1][1]

def face_main_axes_3D (mesh, pos, fid) :
	"""Compute main axis of a face in 3D
	
	.. warning:: assume the face is almost flat.
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `fid` (fid) - id of the face to consider
	
	:Returns: barycenter and axes
	
	:Returns Type: array,array,array
	"""
	#compute local weighted point according to
	#barycenter and connected edge length
	L = dict( (eid,edge_length(mesh,pos,eid) ) for eid in mesh.borders(2,fid) )
	tot = sum(L.itervalues() )
	L = dict( (eid,val / tot) for eid,val in L.iteritems() )
	bary = centroid(mesh,pos,2,fid)
	
	pts = [multiply(subtract(pos[pid],bary),
	                sum(L.get(eid,0.) for eid in mesh.regions(0,pid) ) ) \
	       for pid in mesh.borders(2,fid,2)]
	
	#compute inertia
	I = divide(reduce(add,[outer(pt,pt) for pt in pts]),float(len(pts) ) )
	
	#compute eigen values
	w,v = eig(I)
	
	vecs = [(w[i],multiply(v[:,i],w[i]) ) for i in (0,1,2)]
	vecs.sort(reverse=True)
	
	#return
	return bary,vecs[0][1],vecs[1][1]#,vecs[2][1]

def cell_main_axes (mesh, pos, cid) :
	"""Compute main axes of a cell
	
	.. warning:: assume the geometry of the cell is a polyhedron
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `cid` (cid) - id of the cell to consider
	
	:Returns: barycenter and 3 main axes
	
	:Returns Type: array,array,array,array
	"""
	#compute local weighted point according to
	#barycenter and connected edge length
	L = dict( (eid,edge_length(mesh,pos,eid) ) \
	          for eid in mesh.borders(3,cid,2) )
	tot = sum(L.itervalues() )
	L = dict( (eid,val / tot) for eid,val in L.iteritems() )
	bary = centroid(mesh,pos,3,cid)
	
	pts = [multiply(subtract(pos[pid],bary),
	                sum(L.get(eid,0.) for eid in mesh.regions(0,pid) ) ) \
	       for pid in mesh.borders(3,cid,3)]
	
	#compute inertia
	I = divide(reduce(add,[outer(pt,pt) for pt in pts]),float(len(pts) ) )
	
	#compute eigen values
	w,v = eig(I)
	
	vecs = [(w[i],multiply(v[:,i],w[i]) ) for i in (0,1,2)]
	vecs.sort(reverse=True)
	
	#return
	return bary,vecs[0][1],vecs[1][1],vecs[2][1]

####################################################
#
#		geometrical advanced properties
#
####################################################
def find_nearby_wisp (mesh, pos, degree, point) :
	"""Find the wisp at the given degree whose centroid
	is most proximal to the given point.
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `degree` (int) - degree of wisps to consider
	 - `point` (array) - reference point
	
	:Returns: distance and id of the most proximal wisp
	
	:Returns Type: float,wid
	"""
	dist = [(norm(subtract(centroid(mesh,pos,degree,wid),point) ),wid) \
	         for wid in mesh.wisps(degree)]
	dist.sort()
	return dist[0]

def edge_loop_around (mesh, pos, point) :
	"""Find a loop of edges around a point
	
	.. warning:: the dimension of space must be 2
	
	.. warning:: the loop must exists!
	
	:Parameters:
	 - `mesh` (:class:`Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `point` (array) - point inside the future loop
	
	:Returns: the list of edge id that form a closed loop around the point
	
	:Returns Type: list of eid
	"""
	#find most proximal edge
	edge_prox = []
	for eid in mesh.wisps(1) :
		pt1,pt2 = (pos[pid] for pid in mesh.borders(1,eid) )
		seg_dir = subtract(pt2,pt1)
		l = norm(seg_dir)
		seg_dir = divide(seg_dir,l)
		
		proj = dot(subtract(point,pt1),seg_dir)
		if 0 <= proj <= l :
			d = norm(subtract(subtract(point,pt1),
			                  multiply(seg_dir,proj) ) )
			edge_prox.append( (d,eid) )
	
	edge_prox.sort()
	starting_eid = edge_prox[0][1]
	
	#rotation clockwise
	pid1,pid2 = mesh.borders(1,starting_eid)
	prod = cross(subtract(pos[pid1],point),subtract(pos[pid2],point) )
	if prod < 0 :
		cycle = [pid1]
		front = (starting_eid,pid2,pid1)
	else :
		cycle = [pid2]
		front = (starting_eid,pid1,pid2)
	
	#follow edge up to close the loop
	edge_cycle = [starting_eid]
	while front[1] != cycle[0] :
		last_eid,last_pid,prev_pid = front
		cycle.append(last_pid)
		#sort all possible paths to find the one
		#that makes the smallest angle
		angles = []
		ref_dir = subtract(pos[prev_pid],pos[last_pid])
		for eid in mesh.regions(0,last_pid) :
			if eid != last_eid :
				pid, = set(mesh.borders(1,eid) ) - set([last_pid])
				cur_dir = subtract(pos[pid],pos[last_pid])
				alpha = atan2(cross(ref_dir,cur_dir),
				              dot(ref_dir,cur_dir) )
				if alpha < 0 :
					alpha += 2 * pi
				angles.append( (alpha,eid) )
		
		angles.sort()
		last_eid = angles[0][1]
		prev_pid = last_pid
		last_pid, = set(mesh.borders(1,last_eid) ) - set([last_pid])
		front = (last_eid,last_pid,prev_pid)
		edge_cycle.append(last_eid)
	
	#return
	return edge_cycle


