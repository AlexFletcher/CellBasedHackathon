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
This module defines functions to create a tissue from a segmented 3D image
"""

__license__ = "Cecill-C"
__revision__ = " $Id: $ "

from sys import stdout
from numpy import array,add,outer
from numpy.linalg import eig
from sa_oa.celltissue import Tissue,TissueDB,Config,ConfigItem

def extract_graph_tissue (img, debug_info = True) :
	"""Construct a tissue from a segmented image
	
	.. warning:: by convention index 0 correspond to the external part of the
	             image and index 1 correspond to the external part of the
	             tissue. Hence, the constructed tissue will have a cell 0
	             and 1 without volume
	
	.. warning:: the geometrical center of voxel i,j,k is
	             (i + 0.5,j + 0.5,k + 0.5) * im.resolution
	
	.. warning:: in the image i correspond to y and j to x
	
	:Parameters:
	 - `img` (LxMxN array of int) - a 3D array where each voxel contains the
	                                the id of the cell it belongs to.
	 - `debug_info` (bool) - tells the algo to display evolution
	
	:Returns: a tissue that contains :
	 - a graph from cell neighborhood
	 - the volume of each cell
	 - the contact surface of each wall (edge) between two cells
	 - the barycenter of each cell
	 - the principal directions of each cell
	 - the umber of voxels in each cell
	 - the number of voxel faces in each wall
	
	:Returns Type: TissueDB
	"""
	imax,jmax,kmax = img.shape
	try :
		res = array(img.resolution)
	except AttributeError :
		res = array([1.,1.,1.])
	
	Sres = [res[1] * res[2],res[2] * res[0],res[0] * res[1] ]
	Vres = res[0] * res[1] * res[2]
	
	#extract infos
	vox = [[] for i in xrange(img.max() + 1)]
	facets = {}
	facets_nb = {}
	
	if debug_info :
		print "infos",img.shape
	
	for k in range(kmax) :
		if debug_info :
			print k,
			stdout.flush()
		
		for i in range(imax) :
			for j in range(jmax) :
				cid = img[i,j,k]
				if cid > 1 :
					#cell voxels
					vox[cid].append( (i + 0.5,j + 0.5,k + 0.5) * res)
					
					#surface facets
					if i == 0 :
						nid = 0
						key = (0,cid)
						facets[key] = facets.get(key,0) + Sres[1]
						facets_nb[key] = facets_nb.get(key,0) + 1
					elif i == (imax - 1) :
						nid = 0
						key = (0,cid)
						facets[key] = facets.get(key,0) + Sres[1]
						facets_nb[key] = facets_nb.get(key,0) + 1
					else :
						nid = img[i + 1,j,k]
						if nid != cid :
							key = (min(cid,nid),max(cid,nid) )
							facets[key] = facets.get(key,0) + Sres[1]
							facets_nb[key] = facets_nb.get(key,0) + 1
					
					if j == 0 :
						nid = 0
						key = (0,cid)
						facets[key] = facets.get(key,0) + Sres[0]
						facets_nb[key] = facets_nb.get(key,0) + 1
					elif j == (jmax - 1) :
						nid = 0
						key = (0,cid)
						facets[key] = facets.get(key,0) + Sres[0]
						facets_nb[key] = facets_nb.get(key,0) + 1
					else :
						nid = img[i,j + 1,k]
						if nid != cid :
							key = (min(cid,nid),max(cid,nid) )
							facets[key] = facets.get(key,0) + Sres[0]
							facets_nb[key] = facets_nb.get(key,0) + 1
					
					if k == 0 :
						nid = 0
						key = (0,cid)
						facets[key] = facets.get(key,0) + Sres[2]
						facets_nb[key] = facets_nb.get(key,0) + 1
					elif k == (kmax - 1) :
						nid = 0
						key = (0,cid)
						facets[key] = facets.get(key,0) + Sres[2]
						facets_nb[key] = facets_nb.get(key,0) + 1
					else :
						nid = img[i,j,k + 1]
						if nid != cid :
							key = (min(cid,nid),max(cid,nid) )
							facets[key] = facets.get(key,0) + Sres[2]
							facets_nb[key] = facets_nb.get(key,0) + 1
					
	if debug_info :
		print "#end of war"
	
	#create tissue
	t = Tissue()
	CELL = t.add_type("CELL")
	WALL = t.add_type("WALL")
	
	#properties
	V = {} #cell volume
	S = {} #wall surface
	nb = {} #number of elements that geometricaly define an element
	bary = {} #cell barycenter
	axes = {} #cell main axes
	
	#create graph
	graph_id = t.add_relation("graph",(CELL,WALL) )
	graph = t.relation(graph_id)
	
	for cid in range(len(vox) ) :
		graph.add_vertex(cid)
		V[cid] = len(vox[cid]) * Vres
		nb[cid] = len(vox[cid])
	
	for key in facets :
		eid = graph.add_edge(*key)
		S[eid] = facets[key]
		nb[eid] = facets_nb[key]
	
	#geometrical descriptors
	if debug_info :
		print "geometrical descriptors"
	
	for cid in graph.vertices() :
		pts = vox[cid]
		
		if len(pts) == 0 :
			bary[cid] = None
			axes[cid] = None
		else :
			b = reduce(add,pts) / len(pts)
			
			#compute inertia
			I = reduce(add,[outer(pt - b,pt - b) for pt in pts]) / float(len(pts) )
			
			#compute eigen values
			w,v = eig(I)
			
			vecs = [(w[i],tuple(float(a) for a in v[:,i] * w[i]) ) \
			         for i in (0,1,2)]
			vecs.sort(reverse=True)
			
			#return
			bary[cid] = b
			axes[cid] = [v for w,v in vecs]
	
	#write result
	db = TissueDB()
	db.set_tissue(t)
	
	cfg = Config("topology")
	cfg.add_item(ConfigItem("CELL",CELL) )
	cfg.add_item(ConfigItem("WALL",WALL) )
	cfg.add_item(ConfigItem("graph_id",graph_id) )
	
	db.set_config("config",cfg)
	
	db.set_property("bary",bary)
	db.set_description("bary","barycenter of cells")
	
	db.set_property("V",V)
	db.set_description("V","cell volume")
	
	db.set_property("S",S)
	db.set_description("S","wall surface")
	
	db.set_property("nb",nb)
	db.set_description("nb",
	                   "number of geom elements used to define this element")
	
	db.set_property("axes",axes)
	db.set_description("axes","main axes of the cells")
	
	#return
	return db





