#!BPY
""" 
Name: 'topomesh'
Blender: 244
Group: 'Import'
Tooltip: 'Import from topomesh file format (.msh)'
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
This module import a mesh in Blender
'''

__docformat__ = "restructuredtext"
__license__ = "Cecill-C"
__revision__ = " $Id: blender_topomesh_import.py 7977 2010-02-15 15:27:44Z chopard $ "

from random import randint
import Blender
from Blender import Window,Mesh,Object,Scene
from sa_oa.container import read_topomesh

def load_topomesh (filename, obj):
	"""read a mesh and create it in Blender.
	
	associate the mesh to the given object
	"""
	#read mesh
	mesh,descr,props = read_topomesh(filename)
	
	#create blender mesh
	me = Mesh.New('myMesh')
	
	#create points
	X, = (prop[1] for prop in props[0] \
	      if prop[0] == 'X')
	Y, = (prop[1] for prop in props[0] \
	      if prop[0] == 'Y')
	Z = [prop[1] for prop in props[0] \
	      if prop[0] == 'Z']
	if len(Z) == 1 :
		Z, = Z
	else :
		Z = dict( (pid,0) for pid in X)
	
	trans = {}
	coords = []
	for pid in X :
		trans[pid] = len(coords)
		coords.append( (X[pid],Y[pid],Z[pid]) )
	
	me.verts.extend(coords)
	for i in xrange(len(coords) ) :
		me.verts[i].sel = 0
	
	#create edges
	edges = [tuple(trans[pid] for pid in mesh.borders(1,eid) ) \
	         for eid in mesh.wisps(1) ]
	me.edges.extend(edges)
	
	#link to obj
	obj.link(me)
	
	#fill faces
	me.vertexColors = 1
	for wid in mesh.wisps(2) :
		R = randint(0,255)
		G = randint(0,255)
		B = randint(0,255)
		for pid in mesh.borders(2,wid,2) :
			me.verts[trans[pid] ].sel = 1
		ind_first_face = len(me.faces)
		me.fill()
		ind_last_face = len(me.faces)
		for i in xrange(ind_first_face,ind_last_face) :
			for j in xrange(len(me.faces[i].verts) ) :
				me.faces[i].col[j].r = R
				me.faces[i].col[j].g = G
				me.faces[i].col[j].b = B
		
		for pid in mesh.borders(2,wid,2) :
			me.verts[trans[pid] ].sel = 0
	
	#return
	return me

def create_topomesh (filename) :
	Window.WaitCursor(True)
	
	#link to object
	obj = Object.New('Mesh','myObj')
	
	#link to scene
	sc = Scene.GetCurrent()
	sc.link(obj)
	
	me = load_topomesh(filename,obj)
	
	Window.RedrawAll()


if __name__ == '__main__' :
	Window.FileSelector(create_topomesh, 'Import topomesh', '*.msh')
