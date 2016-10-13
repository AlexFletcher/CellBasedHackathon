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

__doc__="""
This module defines a set of functions to draw the visual representation
associated to a tissue
"""

__license__= "Cecill-C"
__revision__=" $Id: $ "

from sa_oa.svgdraw import SVGScene,SVGLayer,\
                             SVGText,SVGSphere,SVGBox,\
                             SVGConnector,SVGPath,\
                             expand_path

def planar_mesh () :
	"""Visual representation of (CELL,WALL,POINT)
	
	return an SVGScene
	"""
	sc = SVGScene(650,400)
	
	lay = SVGLayer("elements",sc.width(),sc.height(),"lay1")
	sc.append(lay)
	
	#cell1
	cell1 = SVGSphere(200,200,50,50,"cell1")
	cell1.set_fill(None)
	cell1.set_stroke( (255,0,255) )
	cell1.set_stroke_width(3)
	lay.append(cell1)
	cell_txt = SVGText(200 - 32,200 + 8,"CELL",24,"cell_txt")
	cell_txt.set_fill( (255,0,255) )
	lay.append(cell_txt)
	#cell2
	cell2 = SVGSphere(400,200,50,50,"cell2")
	cell2.set_fill(None)
	cell2.set_stroke( (255,0,255) )
	cell2.set_stroke_width(3)
	lay.append(cell2)
	
	#points
	point_txt = SVGText(300 - 35,70,"POINT",24,"point_txt")
	point_txt.set_fill( (0,255,0) )
	lay.append(point_txt)
	pts = []
	for i in xrange(3) :
		for j in xrange(2) :
			pt = SVGSphere(i * 200 + 100,j * 200 + 100,10,10,"point%d" % (i * 2 + j) )
			pt.set_fill( (0,255,0) )
			lay.append(pt)
			pts.append(pt)
	
	#walls
	wall_txt = SVGText(530,200 + 10,"WALL",24,"wall_txt")
	wall_txt.set_fill( (0,0,255) )
	lay.append(wall_txt)
	walls = []
	for i in xrange(3) :
		wall = SVGBox(i * 200 + 100 + (i - 2) * 8,150,16,100,"wall%d" % len(walls) )
		wall.set_fill( (0,0,255) )
		lay.append(wall)
		walls.append(wall)
	
	for i in xrange(2) :
		for j in xrange(2) :
			wall = SVGBox(i * 200 + 150,j * 200 + 100 + (j - 1) * 16,100,16,"wall%d" % len(walls) )
			wall.set_fill( (0,0,255) )
			lay.append(wall)
			walls.append(wall)
	
	#connectors
	lay = SVGLayer("mesh_id",sc.width(),sc.height(),"lay2")
	sc.append(lay)
	
	#cell1 walls
	for i,wall_ind in enumerate( (0,1,3,4) ) :
		con = SVGConnector(cell1.id(),walls[wall_ind].id(),"c1w_con%d" % i)
		con.set_stroke( (0,0,0) )
		con.set_stroke_width(2)
		lay.append(con)
	
	#cell2 walls
	for i,wall_ind in enumerate( (1,2,5,6) ) :
		con = SVGConnector(cell2.id(),walls[wall_ind].id(),"c2w_con%d" % i)
		con.set_stroke( (0,0,0) )
		con.set_stroke_width(2)
		lay.append(con)
	
	#walls pts
	for i,(pt_ind,wall_ind) in enumerate([(0,0),(0,3),
	                            (1,0),(1,4),
	                            (2,1),(2,3),(2,5),
	                            (3,1),(3,4),(3,6),
	                            (4,2),(4,5),
	                            (5,2),(5,6)]) :
		con = SVGConnector(walls[wall_ind].id(),pts[pt_ind].id(),"wp_con%d" % i)
		con.set_stroke( (0,0,0) )
		con.set_stroke_width(2)
		lay.append(con)
	
	#return
	expand_path(sc)
	return sc

def add_graph_layer (sc) :
	"""Add a layer representing a bidirectional
	graph between cells.
	"""
	lay = sc.get_layer("elements")
	
	edge1_txt = SVGText(320,260 + 20,"EDGE",24,"edge1_txt")
	edge1_txt.set_fill( (255,115,0) )
	lay.append(edge1_txt)
	
	edge2_txt = SVGText(320,140,"EDGE",24,"edge2_txt")
	edge2_txt.set_fill( (255,115,0) )
	lay.append(edge2_txt)
	
	lay = SVGLayer("graph_id",sc.width(),sc.height(),"lay3")
	sc.append(lay)
	
	edge1 = SVGPath("edge1")
	edge1.set_fill(None)
	edge1.set_stroke( (255,115,0) )
	edge1.set_stroke_width(3)
	edge1.move_to(230,230)
	edge1.curve_to( (280,270),(320,270),(370,230) )
	lay.append(edge1)
	
	edge2 = SVGPath("edge2")
	edge2.set_fill(None)
	edge2.set_stroke( (255,115,0) )
	edge2.set_stroke_width(3)
	edge2.move_to(230,170)
	edge2.curve_to( (280,130),(320,130),(370,170) )
	lay.append(edge2)
	



