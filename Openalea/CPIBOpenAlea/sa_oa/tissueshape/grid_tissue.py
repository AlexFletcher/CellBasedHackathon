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
This module defines functions to create predefined tissues
with regular geometry
"""

__license__ = "Cecill-C"
__revision__ = " $Id: $ "

from math import sqrt,sin,cos,pi
from sa_oa.container import Quantity,Grid
from sa_oa.celltissue import Tissue,Config,ConfigItem,TissueDB
from sa_oa.svgdraw import to_xml
from tissue_visual_descr import planar_mesh

def _single_cell1D () :
	#create tissue
	t = Tissue()
	CELL = t.add_type("cell")
	WALL = t.add_type("wall")
	mesh_id = t.add_relation("mesh",(WALL,CELL) )
	mesh = t.relation(mesh_id)
	
	#create cells
	cid = mesh.add_wisp(1)
	
	#create walls
	walls = [mesh.add_wisp(0) for i in (0,1)]
	
	#link walls and cells
	for wid in walls :
		mesh.link(1,cid,wid)
	
	#create property position of walls
	pos = dict( (wid,i) for i,wid in enumerate(walls) )
	
	#create config
	cfg = Config("topology")
	cfg.add_item(ConfigItem("CELL",CELL) )
	cfg.add_item(ConfigItem("WALL",WALL) )
	cfg.add_item(ConfigItem("mesh_id",mesh_id) )
	
	#return
	return t,cfg,pos

def _single_cell2D (nb_sides) :
	#create tissue
	t = Tissue()
	CELL = t.add_type("cell")
	WALL = t.add_type("wall")
	POINT = t.add_type("point")
	mesh_id = t.add_relation("mesh",(POINT,WALL,CELL) )
	mesh = t.relation(mesh_id)
	
	#create cells
	cid = mesh.add_wisp(2)
	
	#create walls
	walls = [mesh.add_wisp(1) for i in xrange(nb_sides)]
	
	#create points
	point = [mesh.add_wisp(0) for i in xrange(nb_sides)]
	
	#link walls and cells
	for wid in walls :
		mesh.link(2,cid,wid)
	
	#link walls and points
	for i,wid in enumerate(walls) :
		mesh.link(1,wid,point[i])
		mesh.link(1,wid,point[(i + 1) % nb_sides])
	
	#create property position of points
	pos = dict( (pid,(cos(2 * pi * i / nb_sides),
	                  sin(2 * pi * i / nb_sides) ) ) \
	             for i,pid in enumerate(point) )
	
	#create config
	cfg = Config("topology")
	cfg.add_item(ConfigItem("CELL",CELL) )
	cfg.add_item(ConfigItem("WALL",WALL) )
	cfg.add_item(ConfigItem("POINT",POINT) )
	cfg.add_item(ConfigItem("mesh_id",mesh_id) )
	
	#return
	return t,cfg,pos

def single_cell (nb_sides = 2) :
	"""Create a single polygonal cell
	
	:Parameters:
	 - `nb_sides` (int) - number of sides of the cell. If less than 3, the
	   returned tissue will be 1D
	
	:Return: a tissuedb containing :
	            - a mesh
	            - a config
	            - a pos property
	
	:Returns Type: TissueDB
	"""
	if nb_sides < 3 :
		t,cfg,pos = _single_cell1D()
	else :
		t,cfg,pos = _single_cell2D(nb_sides)
	
	#create tissueDB
	db = TissueDB()
	db.set_tissue(t)
	db.set_property("position",Quantity(pos,"None") )
	db.set_description("position","geometrical position of points")
	db.set_config("config",cfg)
	
	return db

def _regular_grid1D (shape) :
	if type(shape) == int :
		nb_cells = shape
	else :
		nb_cells, = shape
	
	#create tissue
	t = Tissue()
	CELL = t.add_type("cell")
	WALL = t.add_type("wall")
	mesh_id = t.add_relation("mesh",(WALL,CELL) )
	mesh = t.relation(mesh_id)
	
	#create cells
	cells = [mesh.add_wisp(1) for i in xrange(nb_cells)]
	
	#create walls
	wall = [mesh.add_wisp(0) for i in xrange(nb_cells + 1)]
	
	#link walls and cells
	for i,cid in enumerate(cells) :
		mesh.link(1,cid,wall[i])
		mesh.link(1,cid,wall[i + 1])
	
	#create property position of walls
	pos = dict( (wid,i) for i,wid in enumerate(wall) )
	
	#create config
	cfg = Config("topology")
	cfg.add_item(ConfigItem("CELL",CELL) )
	cfg.add_item(ConfigItem("WALL",WALL) )
	cfg.add_item(ConfigItem("mesh_id",mesh_id) )
	
	#return
	return t,cfg,pos

def _regular_grid2D (shape) :
	cell_grid = Grid(shape)
	point_grid = Grid( (i + 1 for i in shape) )
	imax,jmax = shape
	
	#create tissue
	t = Tissue()
	CELL = t.add_type("cell")
	WALL = t.add_type("wall")
	POINT = t.add_type("point")
	mesh_id = t.add_relation("mesh",(POINT,WALL,CELL) )
	mesh = t.relation(mesh_id)
	
	#create cells
	cell = [mesh.add_wisp(2) for ind in cell_grid]
	
	#create points
	point = [mesh.add_wisp(0) for ind in point_grid]
	
	#create vertical walls
	for i in xrange(imax + 1) :
		for j in xrange(jmax) :
			wid = mesh.add_wisp(1)
			mesh.link(1,wid,point[point_grid.index( (i,j) )])
			mesh.link(1,wid,point[point_grid.index( (i,j + 1) )])
			if i > 0 :
				mesh.link(2,cell[cell_grid.index( (i - 1,j) )],wid)
			if i < imax :
				mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
	
	#create horizontal walls
	for i in xrange(imax) :
		for j in xrange(jmax + 1) :
			wid = mesh.add_wisp(1)
			mesh.link(1,wid,point[point_grid.index( (i,j) )])
			mesh.link(1,wid,point[point_grid.index( (i + 1,j) )])
			if j > 0 :
				mesh.link(2,cell[cell_grid.index( (i,j - 1) )],wid)
			if j < jmax :
				mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
	
	#create property position of points
	pos = dict( (pid,point_grid.coordinates(i) ) for i,pid in enumerate(point) )
	
	#create config
	cfg = Config("topology")
	cfg.add_item(ConfigItem("CELL",CELL) )
	cfg.add_item(ConfigItem("WALL",WALL) )
	cfg.add_item(ConfigItem("POINT",POINT) )
	cfg.add_item(ConfigItem("mesh_id",mesh_id) )
	
	#return
	return t,cfg,pos

def _regular_grid3D (shape) :
	cell_grid = Grid(shape)
	point_grid = Grid( (i + 1 for i in shape) )
	imax,jmax,kmax = shape
	
	#create tissue
	t = Tissue()
	CELL = t.add_type("cell")
	WALL = t.add_type("wall")
	EDGE = t.add_type("edge")
	POINT = t.add_type("point")
	mesh_id = t.add_relation("mesh",(POINT,EDGE,WALL,CELL) )
	mesh = t.relation(mesh_id)
	
	#create cells
	cell = [mesh.add_wisp(3) for ind in cell_grid]
	
	#create points
	point = [mesh.add_wisp(0) for ind in point_grid]
	
	#create edges
	point_to_edge = {}
	#between i and i + 1
	for i in xrange(imax) :
		for j in xrange(jmax + 1) :
			for k in xrange(kmax + 1) :
				eid = mesh.add_wisp(1)
				pid1 = point[point_grid.index( (i,j,k) )]
				pid2 = point[point_grid.index( (i + 1,j,k) )]
				mesh.link(1,eid,pid1)
				mesh.link(1,eid,pid2)
				key = (min(pid1,pid2),max(pid1,pid2) )
				point_to_edge[key] = eid
	#between j and j + 1
	for i in xrange(imax + 1) :
		for j in xrange(jmax) :
			for k in xrange(kmax + 1) :
				eid = mesh.add_wisp(1)
				pid1 = point[point_grid.index( (i,j,k) )]
				pid2 = point[point_grid.index( (i,j + 1,k) )]
				mesh.link(1,eid,pid1)
				mesh.link(1,eid,pid2)
				key = (min(pid1,pid2),max(pid1,pid2) )
				point_to_edge[key] = eid
	#between k and k + 1
	for i in xrange(imax + 1) :
		for j in xrange(jmax + 1) :
			for k in xrange(kmax) :
				eid = mesh.add_wisp(1)
				pid1 = point[point_grid.index( (i,j,k) )]
				pid2 = point[point_grid.index( (i,j,k + 1) )]
				mesh.link(1,eid,pid1)
				mesh.link(1,eid,pid2)
				key = (min(pid1,pid2),max(pid1,pid2) )
				point_to_edge[key] = eid
	
	#create walls
	#in the (i,j) plane, k constant
	for k in xrange(kmax + 1) :
		for i in xrange(imax) :
			for j in xrange(jmax) :
				wid = mesh.add_wisp(2)
				pids = (point[point_grid.index( (i,j,k) )],
				        point[point_grid.index( (i + 1,j,k) )],
				        point[point_grid.index( (i + 1,j + 1,k) )],
				        point[point_grid.index( (i,j + 1,k) )] )
				for ind in xrange(4) :
					pid1 = pids[ind]
					pid2 = pids[(ind + 1) % 4]
					eid = point_to_edge[(min(pid1,pid2),max(pid1,pid2) )]
					mesh.link(2,wid,eid)
				if k > 0 :
					mesh.link(3,cell[cell_grid.index( (i,j,k - 1) )],wid)
				if k < kmax :
					mesh.link(3,cell[cell_grid.index( (i,j,k) )],wid)
	
	#in the (j,k) plane, i constant
	for i in xrange(imax + 1) :
		for j in xrange(jmax) :
			for k in xrange(kmax) :
				wid = mesh.add_wisp(2)
				pids = (point[point_grid.index( (i,j,k) )],
				        point[point_grid.index( (i,j + 1,k) )],
				        point[point_grid.index( (i,j + 1,k + 1) )],
				        point[point_grid.index( (i,j,k + 1) )] )
				for ind in xrange(4) :
					pid1 = pids[ind]
					pid2 = pids[(ind + 1) % 4]
					eid = point_to_edge[(min(pid1,pid2),max(pid1,pid2) )]
					mesh.link(2,wid,eid)
				if i > 0 :
					mesh.link(3,cell[cell_grid.index( (i - 1,j,k) )],wid)
				if i < imax :
					mesh.link(3,cell[cell_grid.index( (i,j,k) )],wid)
	
	#in the (k,i) plane, j constant
	for j in xrange(jmax + 1) :
		for k in xrange(kmax) :
			for i in xrange(imax) :
				wid = mesh.add_wisp(2)
				pids = (point[point_grid.index( (i,j,k) )],
				        point[point_grid.index( (i,j,k + 1) )],
				        point[point_grid.index( (i + 1,j,k + 1) )],
				        point[point_grid.index( (i + 1,j,k) )] )
				for ind in xrange(4) :
					pid1 = pids[ind]
					pid2 = pids[(ind + 1) % 4]
					eid = point_to_edge[(min(pid1,pid2),max(pid1,pid2) )]
					mesh.link(2,wid,eid)
				if j > 0 :
					mesh.link(3,cell[cell_grid.index( (i,j - 1,k) )],wid)
				if j < jmax :
					mesh.link(3,cell[cell_grid.index( (i,j,k) )],wid)
	
	#create property position of points
	pos = dict( (pid,point_grid.coordinates(i) ) for i,pid in enumerate(point) )
	
	#create config
	cfg = Config("topology")
	cfg.add_item(ConfigItem("CELL",CELL) )
	cfg.add_item(ConfigItem("WALL",WALL) )
	cfg.add_item(ConfigItem("EDGE",EDGE) )
	cfg.add_item(ConfigItem("POINT",POINT) )
	cfg.add_item(ConfigItem("mesh_id",mesh_id) )
	
	#return
	return t,cfg,pos

def regular_grid (shape) :
	"""Create a tissue with a regular grid shape
	
	:Parameters:
	 - `shape` (int or tuple of int) - nb of cells in each dimension
	
	:Return: a tissuedb containing :
	            - a mesh
	            - a config
	            - a pos property
	
	:Returns Type: TissueDB
	"""
	if type(shape) == int or len(shape) == 1 :
		t,cfg,pos = _regular_grid1D(shape)
	elif len(shape) == 2 :
		t,cfg,pos = _regular_grid2D(shape)
	elif len(shape) == 3 :
		t,cfg,pos = _regular_grid3D(shape)
	else :
		raise NotImplementedError("nD grid still in the box, %s" % str(shape) )
	
	#create tissueDB
	db = TissueDB()
	db.set_tissue(t)
	db.set_property("position",Quantity(pos,"None") )
	db.set_description("position","geometrical position of points")
	db.set_config("config",cfg)
	
	#add visual descrition
	data = to_xml(planar_mesh() )
	db.set_external_data("visual_descr.svg",data)
	return db

def _hexagonal_grid2D (shape, shape_geom = "hexa") :
	cell_grid = Grid(shape)
	imax,jmax = shape
	assert imax > 1
	assert jmax > 1
	point_grid = Grid( (2 * imax + 1, 3 + (jmax - 1) ) )
	
	#create tissue
	t = Tissue()
	CELL = t.add_type("cell")
	WALL = t.add_type("wall")
	POINT = t.add_type("point")
	mesh_id = t.add_relation("mesh",(POINT,WALL,CELL) )
	mesh = t.relation(mesh_id)
	
	#create cells
	cell = [mesh.add_wisp(2) for ind in cell_grid]
	
	#create points
	point = [mesh.add_wisp(0) for ind in point_grid]
	
	#create horizontal walls
	for i in xrange(imax) :
		for j in xrange(jmax + 1) :
			wid = mesh.add_wisp(1)
			if j % 2 == 0 :
				mesh.link(1,wid,point[point_grid.index( (2 * i,j) )])
				mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j) )])
			else :
				mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j) )])
				mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j) )])
			
			if j > 1 :
				mesh.link(2,cell[cell_grid.index( (i,j - 2) )],wid)
			if j < jmax :
				mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
	
	j = jmax
	if j % 2 != 0 :
		for i in xrange(imax) :
			wid = mesh.add_wisp(1)
			mesh.link(1,wid,point[point_grid.index( (2 * i,j + 1) )])
			mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
			
			mesh.link(2,cell[cell_grid.index( (i,j - 1) )],wid)
	else :
		for i in xrange(imax) :
			wid = mesh.add_wisp(1)
			mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
			mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 1) )])		
			
			mesh.link(2,cell[cell_grid.index( (i,j - 1) )],wid)
	
	#create vertical walls
	for i in xrange(imax) :
		for j in xrange(jmax) :
			#bottom left wall
			wid = mesh.add_wisp(1)
			if j % 2 == 0 :
				mesh.link(1,wid,point[point_grid.index( (2 * i,j) )])
				mesh.link(1,wid,point[point_grid.index( (2 * i,j + 1) )])
			else :
				mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j) )])
				mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
			
			mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
			if j > 0 :
				if j % 2 == 0 :
					if i > 0 :
						mesh.link(2,cell[cell_grid.index( (i - 1,j - 1) )],wid)
				else :
					mesh.link(2,cell[cell_grid.index( (i,j - 1) )],wid)
			#bottom right wall
			wid = mesh.add_wisp(1)
			if j % 2 == 0 :
				mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j) )])
				mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
			else :
				mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j) )])
				mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 1) )])
			
			mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
			if j > 0 :
				if j % 2 == 0 :
					mesh.link(2,cell[cell_grid.index( (i,j - 1) )],wid)
				else :
					if i < (imax - 1) :
						mesh.link(2,cell[cell_grid.index( (i + 1,j - 1) )],wid)
	
	#top left, first column
	i = 0
	for j in (2 * k for k in xrange( (jmax + 1) / 2) ) :
		wid = mesh.add_wisp(1)
		mesh.link(1,wid,point[point_grid.index( (i,j + 1) )])
		mesh.link(1,wid,point[point_grid.index( (i,j + 2) )])
		mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
		
	#top right, last column
	i = imax - 1
	for j in (2 * k + 1 for k in xrange(jmax / 2) ) :
		wid = mesh.add_wisp(1)
		mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 1) )])
		mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 2) )])
		mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
	
	#top left and right for top line
	j = jmax - 1
	if j % 2 == 0 :
		#top left
		for i in xrange(1,imax) :
			wid = mesh.add_wisp(1)
			mesh.link(1,wid,point[point_grid.index( (2 * i,j + 1) )])
			mesh.link(1,wid,point[point_grid.index( (2 * i,j + 2) )])
			mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
		#top right
		for i in xrange(imax) :
			wid = mesh.add_wisp(1)
			mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
			mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 2) )])
			mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
	else :
		#top left
		for i in xrange(imax) :
			wid = mesh.add_wisp(1)
			mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
			mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 2) )])
			mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
		#top right
		for i in xrange(imax - 1) :
			wid = mesh.add_wisp(1)
			mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 1) )])
			mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 2) )])
			mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
		
	#create property position of points
	if shape_geom == 'box' :
		pos = dict( (pid,point_grid.coordinates(i) ) for i,pid in enumerate(point) )
	elif shape_geom == 'hexa' :
		pos = {}
		for ind,pid in enumerate(point) :
			i,j = point_grid.coordinates(ind)
			if j % 2 == 0 :
				x = 0.5 + (i / 2) * 3 + (i % 2)
				y = j * sqrt(3) / 2.
			else :
				x = - 1 + ( (i + 1) / 2) * 3 + ( (i + 1) % 2)
				y = j * sqrt(3) / 2.
			pos[pid] = (x,y)
	else :
		raise NotImplementedError("shape_geom: %s not recognized. only 'box' or 'hexa'" % str(shape_geom) )
	
	#create config
	cfg = Config("topology")
	cfg.add_item(ConfigItem("CELL",CELL) )
	cfg.add_item(ConfigItem("WALL",WALL) )
	cfg.add_item(ConfigItem("POINT",POINT) )
	cfg.add_item(ConfigItem("mesh_id",mesh_id) )
	
	#return
	return t,cfg,pos

def hexagonal_grid (shape, shape_geom = 'hexa') :
	"""Create a tissue with a regular hexagonal grid topology
	
	The shape of each cell may either be a box or an hexagon
	
	:Parameters:
	 - `shape` (int or tuple of int) - nb of cells in each dimension
	 - `shape_geom` (str) : either 'hexa' or 'box'
	
	:Return: a tissuedb containing :
	            - a mesh
	            - a config
	            - a pos property
	
	:Returns Type: TissueDB
	"""
	if type(shape) == int or len(shape) == 1 :
		t,cfg,pos = _regular_grid1D(shape)
	elif len(shape) == 2 :
		imax,jmax = shape
		if imax == 1 and jmax == 1 :
			t,cfg,pos = _single_cell2D(6)
		else :
			t,cfg,pos = _hexagonal_grid2D(shape,shape_geom)
	
	#create tissueDB
	db = TissueDB()
	db.set_tissue(t)
	db.set_property("position",Quantity(pos,"None") )
	db.set_description("position","geometrical position of points")
	db.set_config("config",cfg)
	
	return db






