#!BPY
""" 
Name: 'pointset'
Blender: 244
Group: 'Import'
Tooltip: 'Import from pickled point set (.msh)'
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
This module import a point set in Blender
'''

__docformat__ = "restructuredtext"
__license__ = "Cecill-C"
__revision__ = " $Id: blender_pointset.py 14917 2013-09-27 12:28:55Z pradal $ "

from random import randint
import Blender
from Blender import Window,Object,Scene
from pickle import load

def load_point_set (filename):
	"""Read a point set and duplicate
	obj for each location.
	
	Insert every created object in sc
	"""
	Window.WaitCursor(True)
	
	sc = Scene.GetCurrent()
	ref_obj = sc.getActiveObject()
	
	#unselect everything except ref
	for obj in Object.GetSelected() :
		obj.sel = 0
	
	ref_obj.sel = 1
	
	#read points
	pts = load(open(filename,'rb') )
	
	for pid,vec in pts.iteritems() :
		Object.Duplicate()
		obj = sc.getActiveObject()
		obj.setLocation(*vec)
	
	#return
	Window.RedrawAll()

if __name__ == '__main__' :
	Window.FileSelector(load_point_set, 'Import point set', '*.pkl')
