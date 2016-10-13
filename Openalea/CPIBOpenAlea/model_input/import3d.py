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

from trim_walls import *
import os,sys

def import3d(filename):
	dirname = filename[:-4]
	os.system("mkdir %s" % dirname)
	os.system("cp %s %s/%s" % (filename,dirname,filename))
	print 'importing', filename
	db = import_xml(filename)
	print 'TODO: make sure centred and scaled'
	dbt = trim_walls(db)
	set_vtk_strings(db)
	set_vtk_strings(dbt)
	db.write("%s/%s_2d.zip" % (dirname,dirname))
	dbt.write("%s/%s_2dt.zip" % (dirname,dirname))

	zvals=get_zvals(db)

	db3d = set_3d_db(db,zvals)
	db3dt = get_trim_3d(db,db3d,zvals)
	set_vtk_strings(db3d)
	set_vtk_strings(db3dt)
	db3d.write("%s/%s_3d.zip" % (dirname,dirname))
	db3dt.write("%s/%s_3dt.zip" % (dirname,dirname))

	zeds=[0.0,10.0,20.0,30.0]#,100.0,125.0,150.0,175.0,200.0,225.0,250.0]

	dbu3d = set_3d_db(db,zeds)
	dbu3dt = get_trim_3d(db,dbu3d,zeds)
	set_vtk_strings(dbu3d)
	set_vtk_strings(dbu3dt)
	dbu3d.write("%s/%s_u_3d.zip" % (dirname,dirname))
	dbu3dt.write("%s/%s_u_3dt.zip" % (dirname,dirname))
		
	#db_to_vtu(db,dirname,"2d.vtu",["cell","wall","edge"])
	#db_to_vtu(dbt,dirname,"2dt.vtu",["cell","wall","edge"])
	#db_to_vtu(db3d,dirname,"3d.vtu",["cell","wall","edge"])
	#db_to_vtu(db3dt,dirname,"3dt.vtu",["cell","wall","edge"])
	#db_to_vtu(dbu3d,dirname,"u_3d.vtu",["cell","wall","edge"])
	#db_to_vtu(dbu3dt,dirname,"u_3dt.vtu",["cell","wall","edge"])

def set_3d_db(db,zvals):
	if type(zvals)==dict:
		return set_3d_irregular(db,zvals)
		
	elif type(zvals)==list:
		return set_3d_uniform(db,zvals)		

def set_3d_uniform(db,zvals):
	mesh0 = get_mesh(db)
	cfg = db.get_config('config')
	cell_types = cfg.cell_types
	nzeds=len(zvals)

	tissue = Tissue()
	POINT = tissue.add_type("point")
	LINE = tissue.add_type("line")
	WALL = tissue.add_type("wall")
	CELL = tissue.add_type("cell")
	EDGE = tissue.add_type("edge")

	mesh_id = tissue.add_relation("mesh",(POINT,LINE,WALL,CELL) )
	graph_id = tissue.add_relation("graph",(CELL,EDGE) )
	position = db.get_property('position')

	offset_vectors=db.get_property('offset_vectors')

	mesh3d = tissue.relation(mesh_id)
	
	cell_type=db.get_property('cell_type')
	cell_type_3d={}
	db_cids={}
	db_wids={}
	db_lids={}
	db_pids={}


	print 'converting to 3d'

	for i in range(nzeds-1):
		#cells
		for cid in mesh0.wisps(2):
			ncid = mesh3d.add_wisp(3)
			db_cids[(cid,i)]=ncid
			cell_type_3d[ncid]=cell_type[cid]
	for i in range(nzeds):
		#horiz walls
		for cid in mesh0.wisps(2):
			nwid = mesh3d.add_wisp(2)
			db_wids[(cid,i)]=nwid
			if i!=0:
				mesh3d.link(3,db_cids[(cid,i-1)],nwid)
			if i!=nzeds-1:
				mesh3d.link(3,db_cids[(cid,i)],nwid)
	for i in range(nzeds-1):
		#vert walls
		for wid in mesh0.wisps(1):
			nwid = mesh3d.add_wisp(2)
			db_wids[(wid,i)]=nwid
			for cid in mesh0.regions(1,wid):
				mesh3d.link(3,db_cids[(cid,i)],nwid)
	for i in range(nzeds):
		#horiz lines
		for wid in mesh0.wisps(1):
			nlid = mesh3d.add_wisp(1)
			db_lids[(wid,i)]=nlid	
			for cid in mesh0.regions(1,wid):
				mesh3d.link(2,db_wids[(cid,i)],nlid)
			if i!=0:
				mesh3d.link(2,db_wids[(wid,i-1)],nlid)
			if i!=nzeds-1:
				mesh3d.link(2,db_wids[(wid,i)],nlid)
	for i in range(nzeds-1):
		#vert lines
		for pid in mesh0.wisps(0):
			nlid = mesh3d.add_wisp(1)
			db_lids[(pid,i)]=nlid
			for wid in mesh0.regions(0,pid):
				mesh3d.link(2,db_wids[(wid,i)],nlid)
	for i in range(nzeds):
		#points
		for pid in mesh0.wisps(0):
			npid = mesh3d.add_wisp(0)
			db_pids[(pid,i)]=npid
			for wid in mesh0.regions(0,pid):
				mesh3d.link(1,db_lids[(wid,i)],npid)
			if i!=0:
				mesh3d.link(1,db_lids[(pid,i-1)],npid)
			if i!=nzeds-1:
				mesh3d.link(1,db_lids[(pid,i)],npid)	
	newpos3d={}	
	for i in range(nzeds):			
		#add positions		
		for pid in mesh0.wisps(0):
			opos=position[pid]
			newpos3d[db_pids[(pid,i)]]=[opos[0],opos[1],zvals[i]]

	newpos3d={pid:[round(pos[0],4),round(pos[1],4),round(pos[2],4)] for pid,pos in newpos3d.iteritems()}
	print 'referencing point ids'
	#print newpos3d

	position={pid:[round(pos[0],4),round(pos[1],4)] for pid,pos in position.iteritems()}
	
	orig_pid={}
	for pid,pos in newpos3d.iteritems():
		x=pos[0]
		y=pos[1]
		for opid,opos in position.iteritems():
			ox=opos[0]
			oy=opos[1]
			if x==ox and y==oy:
				orig_pid[pid]=opid
	print len(newpos3d),len(orig_pid)

	#reference cells by original cell id and position in z stack	
	orig_cid={}
	vindex={}
	for (ocid,i),ncid in db_cids.iteritems():
		orig_cid[ncid]=ocid
		vindex[ncid]=i

	print 'making graph'
	#make graph
	graph3d = tissue.relation(graph_id)

	wall = {}

	for wid in mesh3d.wisps(2) :
		if mesh3d.nb_regions(2,wid) == 2 :
			cid1,cid2 = mesh3d.regions(2,wid)
			eid1 = graph3d.add_edge(cid1,cid2)
			wall[eid1] = wid
			eid2 = graph3d.add_edge(cid2,cid1)
			wall[eid2] = wid




	#compute volumes and surface areas
	if 'scale_px_to_micron' not in db.properties():
		db.set_property('scale_px_to_micron',(1.0,'microns'))
		db.set_description('scale_px_to_micron','scale: pixels to microns')
	scale_px_to_micron = db.get_property('scale_px_to_micron')
	
	#calculate V and S from zvals and orig mesh to allow 'true' sizes to transferred from 3d mesh with full detail to one with minimal walls and edges
	#V = {cid:V0[orig_cid[cid]]*(zvals[vindex[cid]+1]-zvals[vindex[cid]]) for cid in mesh3d.wisps(3)}
 
	V = dict( (cid,scale_px_to_micron[0]*scale_px_to_micron[0]*scale_px_to_micron[0]*cell_volume(mesh3d,newpos3d,cid) ) for cid in mesh3d.wisps(3) )
	S = dict( (wid,scale_px_to_micron[0]*scale_px_to_micron[0]*face_surface_3D(mesh3d,newpos3d,wid) ) for wid in mesh3d.wisps(2) )


	#configuration
	species_desc={'V':'cell','S':'wall','cell_type':'cell','vindex':'cell','orig_cid':'cell'}

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

	dbprops={"position":(newpos3d,"position of points"),"species_desc":(species_desc,"association of properties with mesh degree"),"vindex":(vindex,"position in z direction"),\
		"orig_cid":(orig_cid,"cell id from orig 2d template"),"cell_type":(cell_type_3d,"type of each cell"),"wall":(wall,"wall corresponding to a given edge"),\
		"V":(V,"Volume of each cell in m2"),"S":(S,"Surface of each wall in m"),"wall":(wall,"dictionary giving wall from edge id"),\
		"scale_px_to_micron":(scale_px_to_micron,"conversion factor [0] from pixels to given unit [1]"),\
		"orig_pid":(orig_pid,"point id from orig 2d template"),"offset_vectors":(offset_vectors, "horizontal offset vectors to shrink cells etc.")}

	db3d=TissueDB()
	db3d.set_tissue(tissue)
	db3d.set_config('config',cfg)
	for pname,(prop,desc) in dbprops.iteritems():
		db3d.set_property(pname,prop)
		db3d.set_description(pname,desc)

	return db3d


def set_3d_irregular(db,zvals):

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

	offset_vectors=db.get_property('offset_vectors')

	mesh3d = tissue.relation(mesh_id)
	
	cell_type=db.get_property('cell_type')
	cell_type_3d={}
	db_cids={}
	db_wids={}
	db_lids={}
	db_pids={}

	horiz_walls=[]

	print 'converting to 3d'
	#add cells
	for cid in mesh0.wisps(2):
		ncid = mesh3d.add_wisp(3)
		db_cids[cid]=ncid
		cell_type_3d[ncid]=cell_type[cid]
	#top walls
	for cid in mesh0.wisps(2):
		nwid = mesh3d.add_wisp(2)
		db_wids[-cid-1]=nwid
		mesh3d.link(3,db_cids[cid],nwid)
		horiz_walls.append(nwid)
	#bottom walls
	for cid in mesh0.wisps(2):
		nwid = mesh3d.add_wisp(2)
		db_wids[cid]=nwid
		mesh3d.link(3,db_cids[cid],nwid)
		horiz_walls.append(nwid)
	#vert walls
	for wid in mesh0.wisps(1):
		nwid = mesh3d.add_wisp(2)
		db_wids[wid]=nwid
		for cid in mesh0.regions(1,wid):
			mesh3d.link(3,db_cids[cid],nwid)
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
				horiz_walls.append(lineage[2][None][0])

						
	print 'removing part cells'
	#remove part cells
	for ocid,cid in db_cids.iteritems():
		if len(zcidrefs[ocid])!=0:
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
						if wid in horiz_walls:
							horiz_walls.remove(wid)
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
						if wid in horiz_walls:
							horiz_walls.remove(wid)
						for lid in lids:
							if mesh3d.nb_regions(1,lid)==0:
								pids=list(mesh3d.borders(1,lid))
								mesh3d.remove_wisp(1,lid)
								for pid in pids:
									if mesh3d.nb_regions(0,pid)==0:
										mesh3d.remove_wisp(0,pid)
										del newpos3d[pid]
	
	newpos3d={pid:[round(pos[0],4),round(pos[1],4),round(pos[2],4)] for pid,pos in newpos3d.iteritems()}
	print 'referencing point ids'
	#print newpos3d

	position={pid:[round(pos[0],4),round(pos[1],4)] for pid,pos in position.iteritems()}
	
	orig_pid={}
	for pid,pos in newpos3d.iteritems():
		x=pos[0]
		y=pos[1]
		for opid,opos in position.iteritems():
			ox=opos[0]
			oy=opos[1]
			if x==ox and y==oy:
				orig_pid[pid]=opid
	print len(newpos3d),len(orig_pid)
			
	#reference cells by original cell id and position in z stack	
	orig_cid={}
	vindex={}
	for ocid,ncids in zcidrefs.iteritems():
		for i,ncid in enumerate(ncids):
			orig_cid[ncid]=ocid
			vindex[ncid]=i
		if len(ncids)==0:
			orig_cid[db_cids[ocid]]=ocid			
			vindex[db_cids[ocid]]=0
	print 'making graph'
	#make graph
	graph3d = tissue.relation(graph_id)

	wall = {}

	for wid in mesh3d.wisps(2) :
		if mesh3d.nb_regions(2,wid) == 2 :
			cid1,cid2 = mesh3d.regions(2,wid)
			eid1 = graph3d.add_edge(cid1,cid2)
			wall[eid1] = wid
			eid2 = graph3d.add_edge(cid2,cid1)
			wall[eid2] = wid




	#compute volumes and surface areas
	if 'scale_px_to_micron' not in db.properties():
		db.set_property('scale_px_to_micron',(1.0,'microns'))
		db.set_description('scale_px_to_micron','scale: pixels to microns')
	scale_px_to_micron = db.get_property('scale_px_to_micron')
	V = dict( (cid,scale_px_to_micron[0]*scale_px_to_micron[0]*scale_px_to_micron[0]*cell_volume(mesh3d,newpos3d,cid) ) for cid in mesh3d.wisps(3) )
	S = dict( (wid,scale_px_to_micron[0]*scale_px_to_micron[0]*face_surface_3D(mesh3d,newpos3d,wid) ) for wid in mesh3d.wisps(2) )
	print 'TODO: retain S,V in trimmed db'

	#configuration
	species_desc={'V':'cell','S':'wall','cell_type':'cell','vindex':'cell','orig_cid':'cell'}

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

	dbprops={"position":(newpos3d,"position of points"),"species_desc":(species_desc,"association of properties with mesh degree"),"vindex":(vindex,"position in z direction"),\
		"orig_cid":(orig_cid,"cell id from orig 2d template"),"cell_type":(cell_type_3d,"type of each cell"),"wall":(wall,"wall corresponding to a given edge"),\
		"V":(V,"Volume of each cell in m2"),"S":(S,"Surface of each wall in m"),"wall":(wall,"dictionary giving wall from edge id"),\
		"scale_px_to_micron":(scale_px_to_micron,"conversion factor [0] from pixels to given unit [1]"),"horiz_walls":(horiz_walls,"list of all horizontal walls"),\
		"orig_pid":(orig_pid,"point id from orig 2d template"),"offset_vectors":(offset_vectors, "horizontal offset vectors to shrink cells etc.")}

	db3d=TissueDB()
	db3d.set_tissue(tissue)
	db3d.set_config('config',cfg)
	for pname,(prop,desc) in dbprops.iteritems():
		db3d.set_property(pname,prop)
		db3d.set_description(pname,desc)

	return db3d

def get_zvals(db):
	mesh = get_mesh(db)
	# define z values: dict by cell in template with associated z values
	#zeds=[0.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0]
	zeds=[0.0,10.0,20.0,30.0]#,40.0,50.0,60.0,70.0,80.0,90.0,100.0]
	#zeds=zeds[:ncells+1]

	zvals={}
	for cid in mesh.wisps(2):
		zvals[cid]=[]
		for i,zed in enumerate(zeds):
			if i==0 or i==len(zeds)-1:
				zvals[cid].append(zed)
			else:
				zvals[cid].append(round(normal(zed,1.0),4))
				#zvals[cid].append(zed)
	#zvals={0:[0.0,50.0,100.0],1:[0.0,50.0,100.0]}
	#zvals=[0.0,50.0,100.0]
	return zvals



def get_trim_3d(db,db3d,zvals):
	db_trim=trim_walls(db)
	cell_trans=db_trim.get_property('cell_trans')

	if type(zvals)==dict:
		zval_trans={}
		for ncid,ocid in cell_trans.iteritems():
			zval_trans[ncid]=zvals[ocid]
	else:
		zval_trans=zvals

	db3dt = set_3d_db(db_trim,zval_trans)

	ocell_trans=cell_trans
	cell_trans={}
	t_orig_cid=db3dt.get_property('orig_cid')
	orig_cid=db3d.get_property('orig_cid')
	t_vindex=db3dt.get_property('vindex')
	vindex=db3d.get_property('vindex')

	
	for tcid,tocid in t_orig_cid.iteritems():
		tvin=t_vindex[tcid]
		xcid=ocell_trans[tocid]
		for cid,ocid in orig_cid.iteritems():
			xvin=vindex[cid]
			if ocid==xcid and xvin==tvin:
				cell_trans[cid]=tcid

	mesh=get_mesh(db3d)
	mesh_t=get_mesh(db3dt)

	
	graph=get_graph(db3d)
	graph_t=get_graph(db3dt)
	wall=db3d.get_property('wall')
	wall_t=db3dt.get_property('wall')
	edge_trans={}
	wall_trans={}
	for cid in graph.vertices():
		for eid in graph.edges(cid):
			tid=graph.target(eid)
			for eidt in graph_t.edges(cell_trans[cid]):
				tidt=graph_t.target(eidt)
				if cell_trans[tid]==tidt:
					edge_trans[eid]=eidt
					wall_trans[wall[eid]]=wall_t[eidt]

	for wid in mesh.wisps(2):
		rcids = list(mesh.regions(2,wid))
		if len(rcids)==1:
			cidt=cell_trans[rcids[0]]
			for twid in mesh_t.borders(3,cidt):
				if len(list(mesh_t.regions(2,twid)))==1:
					wall_trans[wid]=twid
			
	db3d.set_property('cell_trans',cell_trans)
	db3d.set_description('cell_trans','reference to cell in minimal graph version of db')	
	db3d.set_property('wall_trans',wall_trans)
	db3d.set_description('wall_trans','reference to wall in minimal graph version of db')	
	db3d.set_property('edge_trans',edge_trans)
	db3d.set_description('edge_trans','reference to edge in minimal graph version of db')
	#db3d.set_property('trim_db',db3dt)
	#db3d.set_description('trim_db','minimal graph version of db')

	return db3dt
