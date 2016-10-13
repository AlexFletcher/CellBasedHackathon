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
This module defines a set of functions to modify complex meshes
"""

__license__ = "Cecill-C"
__revision__ = " $Id: $ "

from numpy import add,subtract,multiply,divide,dot
from numpy.linalg import norm
from sa_oa.container import topo_divide_edge,\
                               topo_divide_face,\
                               topo_divide_cell
from mesh_geometry import centroid,edge_length

__all__ = ["subdivide_triangle",
           "divide_segment",
           "divide_edge",
           "divide_face",
           "divide_cell",
           "merge_faces_2D",
           "extrude_2D"]

###############################################
#
#		mesh
#
###############################################
def subdivide_triangle (mesh, pos, fid, method = "centroid") :
	"""Subdivide a triangle
	
	:Parameters:
	 - `mesh` (Topomesh)
	 - `pos` (dict of (pid,array) ) - position of points
	 - `fid` (fid) - id of face to divide, must be a triangle
	 - `method` (str) - method used to divide the face. Must be one of:
	
	     - 'centroid': divide triangle in 3 around the centroid
	     - 'median': divide triangle in 2 around the longest mediatrice
	     - 'inscribed': divide triangle in 4 subtriangles
	
	:Returns: for each degree of the mesh, return a mapping between old elements
	 and newly created ones where:
	
	   - None is the key for newly created elements without parent (separation)
	   - unchanged elements are not recorded
	   - removed elements are associated with None
	
	:Returns Type: list of dict of {wid:list of wid}
	"""
	lineage = [{} for i in xrange(mesh.degree() + 1)]
	
	if method == "centroid" :
		cent_pid = mesh.add_wisp(0)
		lineage[0][None] = [cent_pid]
		pos[cent_pid] = centroid(mesh,pos,2,fid)
		
		ray = {}
		lineage[1][None] = []
		for pid in mesh.borders(2,fid,2) :
			eid = mesh.add_wisp(1)
			mesh.link(1,eid,pid)
			mesh.link(1,eid,cent_pid)
			ray[pid] = eid
			lineage[1][None].append(eid)
		
		lineage[2][fid] = []
		for eid in mesh.borders(2,fid) :
			tid = mesh.add_wisp(2)
			lineage[2][fid].append(tid)
			
			if mesh.degree() > 2 :
				for cid in mesh.regions(2,fid) :
					mesh.link(3,cid,tid)
			mesh.link(2,tid,eid)
			for pid in mesh.borders(1,eid) :
				mesh.link(2,tid,ray[pid])
		
		mesh.remove_wisp(2,fid)
	
	elif method == "median" :
		#find longest edge
		d = [(edge_length(mesh,pos,eid),eid) for eid in mesh.borders(2,fid)]
		d.sort()
		eid = d[-1][1]
		
		#subdivide it
		pid1,pid2 = mesh.borders(1,eid)
		eid1,eid2,cent_pid = topo_divide_edge(mesh,eid)
		lineage[0][None] = [cent_pid]
		lineage[1][eid] = (eid1,eid2)
		if pid1 not in mesh.borders(1,eid1) :
			eid1,eid2 = eid2,eid1
		pos[cent_pid] = (pos[pid1] + pos[pid2]) / 2.
		
		#subdivide faces
		med_pids = set([pid1,pid2,cent_pid])
		lineage[1][None] = []
		for fid in tuple(mesh.regions(1,eid1) ) :
			#create median edge
			pid, = set(mesh.borders(2,fid,2) ) - med_pids
			meid = mesh.add_wisp(1)
			lineage[1][None].append(meid)
			mesh.link(1,meid,pid)
			mesh.link(1,meid,cent_pid)
			
			#create first subtriangle
			tid1 = mesh.add_wisp(2)
			if mesh.degree() > 2 :
				for cid in mesh.regions(2,fid) :
					mesh.link(3,cid,tid1)
			
			seid, = set(mesh.regions(0,pid) ) & set(mesh.regions(0,pid1) )
			mesh.link(2,tid1,seid)
			mesh.link(2,tid1,eid1)
			mesh.link(2,tid1,meid)
			
			#create second subtriangle
			tid2 = mesh.add_wisp(2)
			if mesh.degree() > 2 :
				for cid in mesh.regions(2,fid) :
					mesh.link(3,cid,tid2)
			
			seid, = set(mesh.regions(0,pid) ) & set(mesh.regions(0,pid2) )
			mesh.link(2,tid2,seid)
			mesh.link(2,tid2,eid2)
			mesh.link(2,tid2,meid)
			
			#lineage
			lineage[2][fid] = (tid1,tid2)
			mesh.remove_wisp(2,fid)
	
	elif method == "inscribed" :
		pass
	else :
		raise UserWarning("method '%s' unrecognized" % method)
	
	return lineage


###############################################
#
#		division
#
###############################################
def divide_segment (extr1, extr2, plane_normal, point) :
	"""Compute the intersection between a segment (extr1,extr2)
	and a line or plan specified by a point and a normal.
	
	Returns None if the intersection do not exist.
	
	:Parameters:
	 - `extr1` (array) - first extremity of the segment
	 - `extr2` (array) - second extremity of the segment
	 - `plane_normal` (array) - vector normal to the plane
	 - `point` (array) - a point belonging to the plane
	
	:Returns: Position of intersection
	
	:Returns Type: array
	"""
	seg = subtract(extr2,extr1)
	#no test for superposed points
	#assert norm(seg) > 0
	
	dnd = dot(plane_normal,seg)
	if abs(dnd) / norm(seg) / norm(plane_normal) < 1e-6 :
		#plane parallel to the segment
		return None#TODO
	
	u = divide(dot(plane_normal,subtract(point,extr1) ),dnd)
	if 0. <= u <= 1. :#intersection inside the segment
		return add(extr1,multiply(seg,u) )

def divide_edge (mesh, pos, eid, point, axis, ratio_min = 1e-4) :
	"""Divide an edge
	
	Modify mesh in place.
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `fid` (fid) : id of the face to divide
	 - `point` (array) - a point inside the division plane
	 - `axis` (array) - a direction perpendicular to the division plane
	 - `ratio_min` (float) - a number between 0. and 0.5 that tell 
	    the minimum ratio between a daugther segment and the length of the edge
	
	:Returns: for each degree of the mesh, return a mapping between old elements
	 and newly created ones where:
	
	  - None is the key for newly created elements without parent (separation)
	  - unchanged elements are not recorded
	  - removed elements are associated with None
	
	:Returns Type: list of dict of {wid:list of wid}
	"""
	lineage = [{} for i in xrange(mesh.degree() + 1)]
	
	pt1,pt2 = (pos[pid] for pid in mesh.borders(1,eid) )
	
	eid1,eid2,pid = topo_divide_edge(mesh,eid)
	
	pos[pid] = divide_segment(pt1,pt2,axis,point)
	
	#ensure no small wall are created
	d1 = norm(subtract(pt1,pos[pid]) )
	d2 = norm(subtract(pt2,pos[pid]) )
	if divide(d1,d1 + d2) < ratio_min :
		pos[pid] = add(pt1,multiply(subtract(pt2,pt1),ratio_min) )
	elif divide(d2,d1 + d2) < ratio_min :
		pos[pid] = add(pt2,multiply(subtract(pt1,pt2),ratio_min) )
	
	lineage[0][None] = (pid,)
	lineage[1][eid] = (eid1,eid2)
	
	return lineage

def divide_face (mesh, pos, fid, point, axis, ratio_min = 1e-4) :
	"""Divide a face into two faces
	
	Modify mesh in place.
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `fid` (fid) : id of the face to divide
	 - `point` (array) - a point inside the division plane
	 - `axis` (array) - a direction perpendicular to the division plane
	 - `ratio_min` (float) - a number between 0. and 0.5 that tell the 
	    minimum ratio between a daugther segment and the length of the edge
	    that have been divided
	
	:Returns: for each degree of the mesh, return a mapping between old elements
	 and newly created ones where:
	
	  - None is the key for newly created elements without parent (separation)
	  - unchanged elements are not recorded
	  - removed elements are associated with None
	
	:Returns Type: list of dict of {wid:list of wid}
	"""
	#create lineage
	lineage = [{} for i in xrange(mesh.degree() + 1)]
	
	#divide edges
	lineage[0][None] = []
	for eid in tuple(mesh.borders(2,fid) ) :
		if len(set(dot(pos[pid] - point,axis) > 0 \
		            for pid in mesh.borders(1,eid) ) ) > 1 :
			edge_lineage = divide_edge(mesh,pos,eid,point,axis,ratio_min)
			lineage[0][None].extend(edge_lineage[0][None])
			lineage[1][eid] = edge_lineage[1][eid]
	
	#divide face
	pid1,pid2 = lineage[0][None]
	
	fid1,fid2,eid = topo_divide_face(mesh,fid,pid1,pid2)
	lineage[1][None] = (eid,)
	lineage[2][fid] = (fid1,fid2)
	
	#return
	return lineage

def divide_cell (mesh, pos, cid, point, axis) :
	"""Divide a cell into two cells
	
	Modify mesh in place.
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `cid` (cid) : id of the cell to divide
	 - `point` (array) - a point inside the division plane
	 - `axis` (array) - a direction perpendicular to the division plane
	
	:Returns: for each degree of the mesh, return a mapping between old elements
	 and newly created ones where:
	
	  - None is the key for newly created elements without parent (separation)
	  - unchanged elements are not recorded
	  - removed elements are associated with None
	
	:Returns Type: list of dict of {wid:list of wid}
	"""
	#create lineage
	lineage = [{} for i in xrange(mesh.degree() + 1)]
	
	#divide edges
	lineage[0][None] = []
	for eid in tuple(mesh.borders(3,cid,2) ) :
		if len(set(dot(pos[pid] - point,axis) > 0 \
		            for pid in mesh.borders(1,eid) ) ) > 1 :
			edge_lineage = divide_edge(mesh,pos,eid,point,axis)
			lineage[0][None].extend(edge_lineage[0][None])
			lineage[1][eid] = edge_lineage[1][eid]
	
	#divide faces
	lineage[1][None] = []
	new_pids = set(lineage[0][None])
	for fid in tuple(mesh.borders(3,cid) ) :
		face_new_pids = set(mesh.borders(2,fid,2) ) & new_pids
		if len(face_new_pids) == 2 :
			pid1,pid2 = face_new_pids
			fid1,fid2,eid = topo_divide_face(mesh,fid,pid1,pid2)
			lineage[1][None].append(eid)
			lineage[2][fid] = (fid1,fid2)
	
	#divide cell
	cid1,cid2,fid = topo_divide_cell(mesh,cid,lineage[1][None])
	
	lineage[2][None] = (fid,)
	lineage[3][cid] = (cid1,cid2)
	
	#return
	return lineage

###############################################
#
#		edition
#
###############################################
def merge_faces_2D (mesh, pos, fid1, fid2) :
	"""Merge two cells into a single one
	
	Relink borders to keep a closed geometry if possible
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical position of points in space
	 - `fid1` (fid) - id of the first face to merge
	 - `fid2` (fid) - id of the second face to merge
	
	:Returns: id of the newly created face
	
	:Returns Type: fid
	"""
	raise NotImplementedError
	
	frontier_edges = set(mesh.borders(2,cid1)) & set(mesh.borders(2,cid2))
	frontier_points = set()
	for eid in frontier_edges :
		frontier_points.update(mesh.borders(1,eid))
	#remove edges
	for eid in frontier_edges :
		mesh.remove_wisp(1,eid)
	#remove points
	extremal_points = []
	for pid in frontier_points :
		if mesh.nb_regions(0,pid) == 0 :
			mesh.remove_wisp(0,pid)
			del pos[pid]
		else :
			extremal_points.append(pid)
	#simplification of external edges if possible
	for pid in extremal_points :
		if mesh.nb_regions(0,pid) == 2 :
			eid1,eid2 = mesh.regions(0,pid)
			pid1,pid2 = ( set(mesh.borders(1,eid1))|set(mesh.borders(1,eid2)) ) - set([pid])
			mesh.remove_wisp(0,pid)
			del pos[pid]
			neighbor_regions = set()
			for eid in (eid1,eid2) :
				neighbor_regions |= set(mesh.regions(1,eid))
				mesh.remove_wisp(1,eid)
			eid = mesh.add_wisp(1)
			mesh.link(1,eid,pid1)
			mesh.link(1,eid,pid2)
			for cid in neighbor_regions :
				mesh.link(2,cid,eid)
	#simplification of cells
	border_edges = set()
	for cid in (cid1,cid2) :
		border_edges |= set(mesh.borders(2,cid))
		mesh.remove_wisp(2,cid)
	cid = mesh.add_wisp(2)
	for eid in border_edges :
		mesh.link(2,cid,eid)
	#return id of new cell
	return cid

def extrude_2D (mesh, pos, main_direction, point_list, volumic_extrude=True) :
	"""Extrude the edges of a mesh
	
	Extrude perpendicularly to each edge.
	Modify mesh and pos in place.
	
	:Parameters:
	 - `mesh` (:class:`sa_oa.container.Topomesh`)
	 - `pos` (dict of (pid|array) ) - geometrical
	         position of points in space
	 - `main_direction` (array) - general direction of extrusion
	          the norm of this vector will be the extrusion length
	 - `point_list` (list of pid) - a list of points to extrude
	 - `volumic_extrude` (bool) - if true create faces from extruded edges
	
	:Returns: None, modify mesh and pos in place
	"""
	raise NotImplementedError
	
	#compute list of edges linked to points
	edge_list = set()
	for pid in point_list :
		for eid in mesh.regions(0,pid) :
			other_pid, = (bid for bid in mesh.borders(1,eid) if bid != pid)
			if other_pid in point_list :
				edge_list.add(eid)
	#compute normals to each edge
	normals = {}
	for eid in edge_list :
		pid1,pid2 = mesh.borders(1,eid)
		ori = pos[pid2] - pos[pid1]
		n = array([ori[1],-ori[0] ])
		n.normalize()
		if mesh.nb_regions(1,eid) == 1 :
			cid, = mesh.regions(1,eid)
			bary = cell_barycenter_2D(mesh,pos,cid)
			extrude_direction = (pos[pid1]+pos[pid2])/2. - bary
		else :
			extrude_direction = main_direction
		if (n*extrude_direction)<0 :
			n *= -1
		normals[eid] = n
	#compute twin point of each point
	l = norm(main_direction)
	twin = {}
	for pid in point_list :
		n_list = [normals[eid] for eid in mesh.regions(0,pid) if eid in edge_list]
		if len(n_list) == 0 :
			n = main_direction
		else :
			n = reduce(add,n_list)/len(n_list)
		tid = mesh.add_wisp(0)
		eid = mesh.add_wisp(1)
		mesh.link(1,eid,pid)
		mesh.link(1,eid,tid)
		pos[tid] = pos[pid] + n*l
		twin[pid] = tid
	#link twin points
	for eid in edge_list :
		tid = mesh.add_wisp(1)
		if volumic_extrude :
			cid = mesh.add_wisp(2)
			mesh.link(2,cid,eid)
			mesh.link(2,cid,tid)
		for pid in mesh.borders(1,eid) :
			mesh.link(1,tid,twin[pid])
			if volumic_extrude :
				rid, = set(mesh.regions(0,twin[pid])) & set(mesh.regions(0,pid))
				mesh.link(2,cid,rid)

