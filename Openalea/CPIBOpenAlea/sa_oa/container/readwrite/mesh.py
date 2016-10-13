# -*- python -*-
# -*- coding: utf-8 -*-
#
#       Topomesh : container package
#
#       Copyright or  or Copr. 2006 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       VPlants WebSite : https://gforge.inria.fr/projects/vplants/
#
"""
This module provide a function to read .mesh files

see http://www.ann.jussieu.fr/~frey/publications/RT-0253.pdf
for a rought description of this format.
"""

__license__= "Cecill-C"
__revision__=" $Id: mesh.py 14917 2013-09-27 12:28:55Z pradal $ "

from ..topomesh import Topomesh

def _create_face (mesh, pids, edge) :
	"""Internal function used to create a face from an ordered list of points
	
	:Parameters:
	 - mesh (Topomesh) - the mesh in which the face will be creates
	 - `pids` (list of int) - an ordered list of pids
	 - `edge` (dict of ( (pid1,pid2)|eid ) - a map that associate a sorted pair
	          of pids with the eid of the edge that links these two points
	
	:Returns: fid of the created face
	
	:Returns Type: fid
	"""
	fid = mesh.add_wisp(2)
	
	#link points to face and create edges
	nb_pids = len(pids)
	for j in xrange(len(pids) ) :
		pid1 = pids[j]
		pid2 = pids[(j + 1) % nb_pids]
		key = (min(pid1,pid2),max(pid1,pid2) )
		try :
			eid = edge[key]
		except KeyError :
			eid = mesh.add_wisp(1)
			mesh.link(1,eid,pid1)
			mesh.link(1,eid,pid2)
			edge[key] = eid
		
		mesh.link(2,fid,eid)
	
	#return
	return fid

def read (filename) :
	"""Read a .mesh file
	
	:Parameters:
	 - `filename` - (str) - name of the file to read
	
	:Returns:
	 - the constructed mesh
	 - the position of points
	 - the reference of each element as defined in the file.
	   if the element is not explicitely defined in the file, its reference
	   is set to None
	
	:Returns type:
	 - Topomesh
	 - dict of (pid,(float,float,float) )
	 - list of (dict of (wid,int) )
	"""
	#open file
	f = open(filename,'r')
	line = ""
	
	#create mesh
	mesh = Topomesh(3)
	point = []
	edge = {}
	face = {}
	
	pos = {}
	ref = [{} for i in xrange(mesh.degree() + 1)]
	
	#read vertices
	while line != "Vertices" :
		line = f.readline().strip()
	
	nb_pts = int(f.readline().strip() )
	for i in xrange(nb_pts) :
		gr = f.readline().strip().split()
		pid = mesh.add_wisp(0)
		point.append(pid)
		pos[pid] = tuple(float(val) for val in gr[:3])
		ref[0][pid] = int(gr[3])
	
	#read other higher degree elements
	while line != "End" :
		while line == "" :
			line = f.readline().strip()
		
		#read faces
		if line in ("Triangles","Quadrilaterals") :
			nb_faces = int(f.readline().strip() )
			for i in xrange(nb_faces) :
				gr = f.readline().strip().split()
				pids = [point[int(val) - 1] for val in gr[:-1] ]
				fid = _create_face(mesh,pids,edge)
				ref[2][fid] = int(gr[-1])
				pids.sort()
				face[tuple(pids)] = fid
		
		#read cells
		if line == "Hexaedra" :
			raise UserWarning("don't know how to handle Hexaedra")
		
		if line == "Tetrahedra" :
			nb_cells = int(f.readline().strip() )
			for i in xrange(nb_cells) :
				gr = f.readline().strip().split()
				cid = mesh.add_wisp(3)
				pids = [point[int(val) - 1] for val in gr[:-1] ]
				ref[3][cid] = int(gr[-1])
				
				#link face to cells
				for inds in [(0,1,2),(1,2,3),(2,3,0),(3,0,1)] :
					lpids = [pids[j] for j in inds]
					lpids.sort()
					key = tuple(lpids)
					try :
						fid = face[key]
					except KeyError :
						fid = _create_face(mesh,lpids,edge)
						ref[2][fid] = 0
					
					mesh.link(3,cid,fid)
		
		line = f.readline().strip()
	
	#return
	return mesh,pos,ref


