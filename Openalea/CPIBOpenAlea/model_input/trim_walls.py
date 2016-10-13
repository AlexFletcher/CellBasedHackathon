from sa_oa.container import Topomesh,graph
from sa_oa.container import read_topomesh,write_topomesh,Quantity
from sa_oa.tissueshape import centroid
from sa_oa.celltissue import TissueDB
from sa_oa.tissueshape import face_surface_3D, edge_length,cell_volume,face_surface_2D
from model_utils.db_utilities import get_mesh,get_graph
from model_utils.db_geom import updateVS
from import_utils import get_2d_offset_vectors



def trim_walls(db):


	# check for scale of tissue - assume =1.0 if not given - check this if you have problems
	if 'scale_px_to_micron' not in db.properties():
		db.set_property('scale_px_to_micron',(1.0,'microns'))
		db.set_description('scale_px_to_micron','scale: pixels to microns')
	scale_px_to_micron  = db.get_property('scale_px_to_micron')


	mesh0=get_mesh(db)
	graph0=get_graph(db)
	position = db.get_property('position')
	updateVS(db)
	V0=db.get_property('V')
	S0=db.get_property('S')
	Vtemp={}
	Stemp={}
	mesh=Topomesh(2)

	pos_new = {}
	trans = {}
	rev={}
	cells={}

	raw_cell_type = db.get_property('cell_type')
	init_cell_type={}
	cfg = db.get_config('config')
	cell_types = cfg.cell_types


	for pid in mesh0.wisps(0):
		if len(list(mesh0.regions(0,pid)))>2:
			npid = mesh.add_wisp(0)
			pos_new[npid] = position[pid]
			trans[pid] = npid


	wall_cells={}
	for wid in mesh0.wisps(1):
		wall_cells[wid]=sorted(list(mesh0.regions(1,wid)))
	cwdict={}
	for cids in wall_cells.itervalues():
		if cids not in cwdict.itervalues():
			cwdict[len(list(cwdict.iterkeys()))+1]=cids
	cell_walls={}
	for aid,cids in cwdict.iteritems():
		cell_walls[aid]=[]
		for wid,ncids in wall_cells.iteritems():
			if ncids==cids:
				cell_walls[aid].append(wid)


	
	end_points={}
	for aid,wids in cell_walls.iteritems():
		end_points[aid]=[]
		for wid in wids:
			for pid in mesh0.borders(1,wid):
				if len(list(mesh0.regions(0,pid)))>2:
					end_points[aid].append(pid)

	wtrans={}
	for aid,pts in end_points.iteritems():
		nwid = mesh.add_wisp(1)
		mesh.link(1,nwid,trans[pts[0]])
		mesh.link(1,nwid,trans[pts[1]])
		wtrans[aid]=nwid
		Stemp[nwid]=0
		for wid in cell_walls[aid]:
			Stemp[nwid]+=S0[wid]

	ctrans={}
	for cid in mesh0.wisps(2):
		ncid = mesh.add_wisp(2)
		init_cell_type[ncid] = raw_cell_type[cid]
		Vtemp[ncid]=V0[cid]
		ctrans[cid]=ncid

	for aid,cids in cwdict.iteritems():
		for cid in cids:
			mesh.link(2,ctrans[cid],wtrans[aid])


	raw_mesh=mesh
	raw_pos=pos_new


	##########################################
	#
	print "create tissue"
	#
	##########################################
	from sa_oa.celltissue import Tissue

	tissue = Tissue()
	POINT = tissue.add_type("point")
	LINE = tissue.add_type("line")
	WALL = tissue.add_type("wall")
	CELL = tissue.add_type("cell")
	EDGE = tissue.add_type("edge")

	mesh_id = tissue.add_relation("mesh",(POINT,WALL,CELL) )
	graph_id = tissue.add_relation("graph",(CELL,EDGE) )

	##########################################
	#
	print "fill mesh"
	#
	##########################################
	mesh = tissue.relation(mesh_id)

	trans0 = [{} for deg in xrange(3)]
	cell_type={}
	vindex={}
	V={}
	S={}
	cell_trans={}
	wall_trans={}
	#add wisps
	for deg in xrange(3) :
		for wid in raw_mesh.wisps(deg) :
			trans0[deg][wid] = mesh.add_wisp(deg)
			if deg==2:
				cell_type[trans0[deg][wid]]=init_cell_type[wid]
				vindex[trans0[deg][wid]]=0
				V[trans0[deg][wid]]=Vtemp[wid]
				for ocid,ncid in ctrans.iteritems():
					if wid == ncid:
						xcid=ocid
				cell_trans[trans0[deg][wid]]=xcid
			if deg==1:
				S[trans0[deg][wid]]=Stemp[wid]
				for owid,nwid in wtrans.iteritems():
					if wid == nwid:
						xwid=owid		
				wall_trans[trans0[deg][wid]]=cell_walls[xwid]
		
	#link wisps
	for deg in xrange(1,3) :
		for wid in raw_mesh.wisps(deg) :
			for bid in raw_mesh.borders(deg,wid) :
				mesh.link(deg,trans0[deg][wid],trans0[deg - 1][bid])

	pos = Quantity(dict( (trans0[0][pid],vec) for pid,vec in raw_pos.iteritems() ),
		       [],
		           "Vector2",
		           "position of points" )

	##########################################
	#
	print "fill graph"
	#
	##########################################
	graph = tissue.relation(graph_id)

	wall = {}

	for wid in mesh.wisps(1) :
		if mesh.nb_regions(1,wid) == 2 :
			cid1,cid2 = mesh.regions(1,wid)
			eid1 = graph.add_edge(cid1,cid2)
			wall[eid1] = wid
			eid2 = graph.add_edge(cid2,cid1)
			wall[eid2] = wid
	#TODO fix this so it only runs if necessary
	eidtrans={}
	rev_trans={}
	for a,b in cell_trans.iteritems():
		rev_trans[b]=a
	for eid0 in graph0.edges():
		sid0=graph0.source(eid0)
		tid0=graph0.target(eid0)
		sidtrans=rev_trans[sid0]
		tidtrans=rev_trans[tid0]
		for eid in graph.edges():
			sid=graph.source(eid)
			tid=graph.target(eid)
			if sid==sidtrans and tid==tidtrans:
				eidtrans[eid0]=eid
	print len(eidtrans)
	print len(list(graph0.edges()))		
	AUX1=db.get_property('AUX1')
	PIN=db.get_property('PIN')
	LAX=db.get_property('LAX')
	AUXnew={}
	PINnew={}
	LAXnew={}
	for eid0,eid in eidtrans.iteritems():
		AUXnew[eid]=AUX1[eid0]
		PINnew[eid]=PIN[eid0]
		LAXnew[eid]=LAX[eid0]
	##########################################
	#
	print "geometry properties"
	#
	##########################################
	from sa_oa.tissueshape import face_surface_2D,edge_length


	pos = dict( (pid,vec) for pid,vec in pos.iteritems() )

	offset_vectors=get_2d_offset_vectors(mesh,pos)
	##########################################
	#
	print "create config"
	#
	##########################################
	#if 'species_desc' not in db.properties():
	#	species_desc={'V':'cell','S':'wall','cell_type':'cell'}
	#	db.set_property('species_desc', species_desc)
    	#	db.set_description('species_desc', 'type of property, i.e. cell,edge or wall')
	from sa_oa.celltissue import ConfigFormat,ConfigItem

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

	species_desc={'V':'cell','S':'wall','cell_type':'cell'}
	##########################################
	#
	print "return new tissue file"
	#
	##########################################

	dbprops={"position":(pos,"position of points"),\
		"cell_type":(cell_type,"type of each cell"),"wall":(wall,"wall corresponding to a given edge"),"species_desc":(species_desc,"association of properties with mesh degree"),\
		"V":(V,"Volume of each cell in m2"),"S":(S,"Surface of each wall in m"),"wall":(wall,"dictionary giving wall from edge id"),\
		"scale_px_to_micron":(scale_px_to_micron,"conversion factor [0] from pixels to given unit [1]"),\
		"cell_trans":(cell_trans,"translation of new cell ids to original cell ids"),"wall_trans":(wall_trans,"translation of new wall ids to list of original wall ids"),\
		"offset_vectors":(offset_vectors, "horizontal offset vectors to shrink cells etc.")}

	db_trim=TissueDB()
	db_trim.set_tissue(tissue)
	db_trim.set_config('config',cfg)
	for pname,(prop,desc) in dbprops.iteritems():
		db_trim.set_property(pname,prop)
		db_trim.set_description(pname,desc)
	#TODO sort this out
	db_trim.set_property('AUX1',AUXnew)
	db_trim.set_description('AUX1','AUX localisation')
	db_trim.set_property('PIN',PINnew)
	db_trim.set_description('PIN','PIN localisation')
	db_trim.set_property('LAX',LAXnew)
	db_trim.set_description('LAX','LAX localisation')
	return db_trim
