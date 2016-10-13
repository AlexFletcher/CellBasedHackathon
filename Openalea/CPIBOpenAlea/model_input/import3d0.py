from import_xml import import_xml
from model_utils.db_utilities import get_mesh,get_graph
from sa_oa.container import ordered_pids
from model_utils.db_geom import *
import time
from model_output.vtk_output import *
from numpy import array, around
from import_utils import *
from sa_oa.celltissue import Tissue,TissueDB
from sa_oa.tissueshape import face_surface_3D,cell_volume,divide_cell,centroid
from numpy.random import normal
from model_utils.celltissue_util import def_property

from sa_oa.celltissue import ConfigFormat,ConfigItem
from sa_oa.celltissue import topen
from sa_oa.tissueshape import totup

def import3d(filename):
	print 'importing', filename
	db = import_xml(filename)
	#db=TissueDB()
	#db.read("twocells.zip")	
	print 'TODO: make sure centred and scaled'
	ncells=1
	print 'extruding geometry:',ncells,'cells'
	mesh = get_mesh(db)



	# define z values: dict by cell in template with associated z values
	zeds=[0.0,50.0]#,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0]
	#zeds=zeds[:ncells+1]

	zvals={}
	for cid in mesh.wisps(2):
		zvals[cid]=[]
		for zed in zeds:
			zvals[cid].append(round(normal(zed,2.5),4))
			#zvals[cid].append(zed)
	#zvals={0:[0.0,50.0],1:[25.0,75.0]}


	set_3d_db(db,ncells,zvals)

	newdb=TissueDB()
	newdb.read("testext.zip")

	wall_dict = get_wall_strings(newdb)
	cell_dict = get_cell_strings(newdb)
	edge_dict = get_edge_strings(newdb)




	newdb,PIN = def_property(newdb,'PIN7',0.0,'EDGE',"config","")


	graph=get_graph(newdb)

	cell_type=newdb.get_property("cell_type")

	
	for eid in graph.edges():
		sid=graph.source(eid)
		tid=graph.target(eid)
		if cell_type[sid]==22 and cell_type[tid]!=16:
			PIN[eid]=1.0
		if cell_type[sid]==23 and cell_type[tid]!=16:
			PIN[eid]=1.0
		if cell_type[sid]==16 and cell_type[tid]!=3:
			PIN[eid]=1.0		
		if cell_type[sid]==22 and cell_type[tid]==21:
			PIN[eid]=0.0
		if cell_type[sid]==20:
			PIN[eid]=1.0
		if cell_type[sid]==19 and cell_type[tid]==16:
			PIN[eid]=0.0
		if cell_type[sid]==21 and cell_type[tid]!=3:
			PIN[eid]=1.0



	species_desc=newdb.get_property("species_desc")
	species_desc['PIN7']="edge"
	

	newdb.set_property('vtk_strings',{'cell':cell_dict})
	newdb.set_description('vtk_strings','description of cells, walls, edges for vtk format')

	vtk_strings=newdb.get_property('vtk_strings')
	vtk_strings['edge']=edge_dict
	vtk_strings['wall']=wall_dict

		
	db_to_vtu(newdb,"../","ccccc.vtu",["cell","edge","wall"])



def set_3d_db(db,ncells,zvals):

	min_zed=min([zval[0] for zval in zvals.itervalues()])
	max_zed=max([zval[-1] for zval in zvals.itervalues()])

	mesh0 = get_mesh(db)
	cfg = db.get_config('config')
	cell_types = cfg.cell_types


	tissue = Tissue()
	POINT = tissue.add_type("point")
	LINE = tissue.add_type("line")
	WALL = tissue.add_type("wall")
	CELL = tissue.add_type("cell")
	EDGE = tissue.add_type("edge")

	mesh_id = tissue.add_relation("mesh",(POINT,LINE,WALL,CELL) )
	graph_id = tissue.add_relation("graph",(CELL,EDGE) )
	position = db.get_property('position')
	mesh3d = tissue.relation(mesh_id)
	
	cell_type=db.get_property('cell_type')
	cell_type_3d={}
	db_cids={}



	print 'converting to 3d'
	#add cells

	for cid in mesh0.wisps(2):
		ncid = mesh3d.add_wisp(3)
		db_cids[cid]=ncid
		cell_type_3d[ncid]=cell_type[cid]


		

	#add horiz walls and link
	#top walls
	db_wids={}
	for cid in mesh0.wisps(2):
		nwid = mesh3d.add_wisp(2)
		db_wids[-cid-1]=nwid
		mesh3d.link(3,db_cids[cid],nwid)

	#bottom walls
	for cid in mesh0.wisps(2):
		nwid = mesh3d.add_wisp(2)
		db_wids[cid]=nwid
		mesh3d.link(3,db_cids[cid],nwid)

	#vert walls
	for wid in mesh0.wisps(1):
		nwid = mesh3d.add_wisp(2)
		db_wids[wid]=nwid
		for cid in mesh0.regions(1,wid):
			mesh3d.link(3,db_cids[cid],nwid)
	db_lids={}
	#top lines
	for wid in mesh0.wisps(1):
		nlid = mesh3d.add_wisp(1)
		db_lids[-wid-1]=nlid	
		for cid in mesh0.regions(1,wid):
			mesh3d.link(2,db_wids[-cid-1],nlid)
		mesh3d.link(2,db_wids[wid],nlid)
	#bottom lines
	for wid in mesh0.wisps(1):
		nlid = mesh3d.add_wisp(1)
		db_lids[wid]=nlid	
		for cid in mesh0.regions(1,wid):
			mesh3d.link(2,db_wids[cid],nlid)
		mesh3d.link(2,db_wids[wid],nlid)
	#vert lines
	for pid in mesh0.wisps(0):
		nlid = mesh3d.add_wisp(1)
		db_lids[pid]=nlid
		for wid in mesh0.regions(0,pid):
			mesh3d.link(2,db_wids[wid],nlid)
	db_pids={}
	#top points
	for pid in mesh0.wisps(0):
		npid = mesh3d.add_wisp(0)
		db_pids[-pid-1]=npid
		for wid in mesh0.regions(0,pid):
			mesh3d.link(1,db_lids[-wid-1],npid)
		mesh3d.link(1,db_lids[pid],npid)						

	#bottom points			
	for pid in mesh0.wisps(0):
		npid = mesh3d.add_wisp(0)
		db_pids[pid]=npid
		for wid in mesh0.regions(0,pid):
			mesh3d.link(1,db_lids[wid],npid)
		mesh3d.link(1,db_lids[pid],npid)

	#add positions
	newpos3d={}
	for pid in mesh0.wisps(0):
		opos=position[pid]
		newpos3d[db_pids[-pid-1]]=[opos[0],opos[1],min_zed]
		newpos3d[db_pids[pid]]=[opos[0],opos[1],max_zed]
		
	
	print 'making z divisions'
	#divide horizontally by zvals
	zcidrefs={}
	for ocid,cid in db_cids.iteritems():
		div_cid=cid
		zcidrefs[ocid]=[]
		for zed in zvals[ocid]:
			if zed!=min_zed and zed!=max_zed:
				if div_cid in zcidrefs[ocid]:
					zcidrefs[ocid].remove(div_cid)
				lineage = divide_cell (mesh3d, newpos3d, div_cid, array([0,0,zed]), array([0,0,1]))
				ncids = lineage[3][div_cid]
				cell_type_3d[ncids[0]]=cell_type_3d[div_cid]
				cell_type_3d[ncids[1]]=cell_type_3d[div_cid]
				cell_type_3d.pop(div_cid, None)
				centres=[centroid(mesh3d,newpos3d,3,ncids[0])[2],centroid(mesh3d,newpos3d,3,ncids[1])[2]]
				maxz_index=centres.index(max(centres))
				div_cid=ncids[maxz_index]#one of ncids
				zcidrefs[ocid].append(ncids[centres.index(min(centres))])
				zcidrefs[ocid].append(ncids[maxz_index])

						
	print 'removing part cells'
	#remove part cells
	for ocid,cid in db_cids.iteritems():
		topcid=zcidrefs[ocid][0]
		botcid=zcidrefs[ocid][-1]
		if zvals[ocid][0]!=min_zed:
			zcidrefs[ocid].remove(topcid)
			wids=list(mesh3d.borders(3,topcid))
			mesh3d.remove_wisp(3,topcid)			
			for wid in wids:
				if mesh3d.nb_regions(2,wid)==0:
					lids=list(mesh3d.borders(2,wid))
					mesh3d.remove_wisp(2,wid)
					for lid in lids:
						if mesh3d.nb_regions(1,lid)==0:
							pids=list(mesh3d.borders(1,lid))
							mesh3d.remove_wisp(1,lid)
							for pid in pids:
								if mesh3d.nb_regions(0,pid)==0:
									mesh3d.remove_wisp(0,pid)
									del newpos3d[pid]
		
		if zvals[ocid][-1]!=max_zed:
			zcidrefs[ocid].remove(botcid)
			wids=list(mesh3d.borders(3,botcid))
			mesh3d.remove_wisp(3,botcid)			
			for wid in wids:
				if mesh3d.nb_regions(2,wid)==0:
					lids=list(mesh3d.borders(2,wid))
					mesh3d.remove_wisp(2,wid)
					for lid in lids:
						if mesh3d.nb_regions(1,lid)==0:
							pids=list(mesh3d.borders(1,lid))
							mesh3d.remove_wisp(1,lid)
							for pid in pids:
								if mesh3d.nb_regions(0,pid)==0:
									mesh3d.remove_wisp(0,pid)
									del newpos3d[pid]
	
	newpos3d={pid:list(around(array(pos),decimals=3)) for pid,pos in newpos3d.iteritems()}

	#reference cells by original cell id and position in z stack	
	orig_cid={}
	vindex={}
	for ocid,ncids in zcidrefs.iteritems():
		for i,ncid in enumerate(ncids):
			orig_cid[ncid]=ocid
			vindex[ncid]=i

	print 'making graph'
	graph3d = tissue.relation(graph_id)

	wall = {}

	for wid in mesh3d.wisps(2) :
		if mesh3d.nb_regions(2,wid) == 2 :
			cid1,cid2 = mesh3d.regions(2,wid)
			eid1 = graph3d.add_edge(cid1,cid2)
			wall[eid1] = wid
			eid2 = graph3d.add_edge(cid2,cid1)
			wall[eid2] = wid




	#compute geometrical various properties
	if 'scale_px_to_micron' not in db.properties():
		db.set_property('scale_px_to_micron',(1.0,'microns'))
		db.set_description('scale_px_to_micron','scale: pixels to microns')
	scale_px_to_micron = db.get_property('scale_px_to_micron')
	V = dict( (cid,scale_px_to_micron[0]*scale_px_to_micron[0]*scale_px_to_micron[0]*cell_volume(mesh3d,newpos3d,cid) ) for cid in mesh3d.wisps(3) )
	S = dict( (wid,scale_px_to_micron[0]*scale_px_to_micron[0]*face_surface_3D(mesh3d,newpos3d,wid) ) for wid in mesh3d.wisps(2) )


	#configuration
	if 'species_desc' not in db.properties():
		species_desc={'V':'cell','S':'wall','cell_type':'cell','vindex':'cell','orig_cid':'cell'}
		db.set_property('species_desc', species_desc)
    		db.set_description('species_desc', 'type of property, i.e. cell,edge or wall')



	cfg = ConfigFormat(vars() )
	cfg.add_section("elements types")
	cfg.add("POINT")
	cfg.add("LINE")
	cfg.add("WALL")
	cfg.add("CELL")
	cfg.add("EDGE")
	cfg.add("mesh_id")
	cfg.add("graph_id")
	cfg.add("cell_types")


	cfg = cfg.config()

	##########################################
	#
	print "writing tissue file"
	#
	##########################################

	f = topen('testext.zip','w')
	f.write(tissue)
	f.write_config(cfg,"config")
	f.write(newpos3d,"position","position of points")
	f.write(species_desc,"species_desc","association of properties with mesh degree")
	f.write(vindex,"vindex","position in z direction")
	f.write(orig_cid,"orig_cid","cell id from orig 2d template")
	f.write(cell_type_3d,"cell_type","type of each cell")
	f.write(wall,"wall","wall corresponding to a given edge")
	f.write(V,"V","Volume of each cell in m2")
	f.write(S,"S","Surface of each wall in m")
	f.write(wall,"wall","wall")
	f.write(scale_px_to_micron,"scale_px_to_micron","conversion factor [0] from pixels to given unit [1]")
	f.close()

	
