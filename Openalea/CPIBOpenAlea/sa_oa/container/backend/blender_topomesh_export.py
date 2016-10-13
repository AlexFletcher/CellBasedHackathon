#!BPY
""" 
Name: 'topomesh'
Blender: 244
Group: 'Export'
Tooltip: 'Export into topomesh file format (.msh)'
"""

# -*- coding: utf-8 -*-
# -*- python -*-
#
#       OpenAlea.Container
#
#       Copyright 2008-2009 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <jerome.chopard.at.sophia.inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://sa_oa.gforge.inria.fr
#
###############################################################################

'''
This module export a mesh from Blender
'''

__docformat__ = "restructuredtext"
__license__ = "Cecill-C"
__revision__ = " $Id: blender_topomesh_export.py 14917 2013-09-27 12:28:55Z pradal $ "

from random import randint
import Blender
from Blender import Window,Mesh,Object,Scene
from sa_oa.container import Topomesh,Quantity,write_topomesh

def save_topomesh (filename, bldmesh):
	"""Write the given blender mesh
	into a file
	"""
	
	mesh = Topomesh(2)
	X = Quantity({},"pix","float","x coordinate of vertices")
	Y = Quantity({},"pix","float","y coordinate of vertices")
	Z = Quantity({},"pix","float","z coordinate of vertices")
	
	#vertices
	ptrans = {}
	for vtx in bldmesh.verts :
		pid = mesh.add_wisp(0)
		X[pid],Y[pid],Z[pid] = vtx.co
		ptrans[vtx.index] = pid
	
	#edges
	etrans = {}
	for edge in bldmesh.edges :
		pids = [ptrans[edge.v1.index],
		        ptrans[edge.v2.index] ]
		pids.sort()
		try :
			eid = etrans[tuple(pids)]
		except KeyError :
			eid = mesh.add_wisp(1)
			etrans[tuple(pids)] = eid
		
		for pid in pids :
			mesh.link(1,eid,pid)
	
	#faces
	for face in bldmesh.faces :
		fid = mesh.add_wisp(2)
		pids = [ptrans[vtx.index] for vtx in face.v]
		nb = len(pids)
		for i in xrange(nb) :
			bids = [pids[i],pids[(i + 1) % nb]]
			bids.sort()
			mesh.link(2,fid,etrans[tuple(bids)])
	
	#write
	props = [[('X',X),('Y',Y),('Z',Z)],
	         [],
	         [] ]
	write_topomesh(filename,mesh,"from blender",props)

def create_topomesh (filename) :
	Window.WaitCursor(True)
	
	#current object
	sc = Scene.GetCurrent()
	ref_obj = sc.getActiveObject()
	
	save_topomesh(filename,ref_obj.data)
	
	Window.RedrawAll()


if __name__ == '__main__' :
	Window.FileSelector(create_topomesh, 'Export topomesh', '*.msh')
