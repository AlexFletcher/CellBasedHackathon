from model_utils.db_utilities import get_mesh,get_graph,get_wall_decomp
from sa_oa.container import ordered_pids
from model_utils.db_geom import get_offset_vec2
from numpy import array, around

def set_vtk_strings(db):
	mesh=get_mesh(db)
	position=db.get_property('position')

	if mesh.degree()==3:
		wall_dict = get_wall_strings(db)
		cell_dict = get_cell_strings(db)
		edge_dict = get_edge_strings(db)
	if mesh.degree()==2:
		wall_dict,cell_dict,edge_dict = get_2d_strings(db)

	db.set_property('vtk_strings',{'cell':cell_dict,'edge':edge_dict,'wall':wall_dict})
	db.set_description('vtk_strings','description of cells, walls, edges for vtk format')

def get_2d_strings(db):
	offset_vectors=db.get_property('offset_vectors')
	width=0.5 # offset for edges

	position=db.get_property('position')
	mesh=get_mesh(db)
	graph=get_graph(db)

	#cells
	offsets=[]
	connections=[]
	pts=[]
	current_offset=0
	cid_trans={}
	pid_trans={}

	for pid in mesh.wisps(0):
		pts.append([position[pid][0],position[pid][1],0.0])
		pid_trans[pid]=len(pid_trans)

	for cid in mesh.wisps(2):
		pids=ordered_pids(mesh,cid)
		for pid in pids:
			connections.append(pid_trans[pid])
		current_offset+=len(pids)
		offsets.append(current_offset)
		cid_trans[len(cid_trans)]=cid

	position_string=' '.join([' '.join([str(coord) for coord in point]) for point in pts])		
	connect_string=' '.join([str(val) for val in connections])
	offset_string=' '.join([str(val) for val in offsets])
	total_points=len(pid_trans)
	total_cells=len(cid_trans)
	type_string = ' '.join(["7" for num in xrange(total_cells)])

	cell_dict={'position_string':position_string,'connect_string':connect_string,\
			'offset_string':offset_string,'type_string':type_string,\
			'n_cells':str(total_cells),'n_points':str(total_points),'cid_trans':cid_trans}

	#walls
	offsets=[]
	connections=[]
	current_offset=0
	wid_trans={}

	for wid in mesh.wisps(1):
		for pid in mesh.borders(1,wid):
			connections.append(pid_trans[pid])
		current_offset+=2
		offsets.append(current_offset)	
		wid_trans[len(wid_trans)]=wid

	connect_string=' '.join([str(val) for val in connections])
	offset_string=' '.join([str(val) for val in offsets])
	total_points=len(pid_trans)
	total_cells=len(wid_trans)
	type_string = ' '.join(["3" for num in xrange(total_cells)])

	wall_dict={'position_string':position_string,'connect_string':connect_string,\
			'offset_string':offset_string,'type_string':type_string,\
			'n_cells':str(total_cells),'n_points':str(total_points),'cid_trans':wid_trans}

	#edges

	wall=db.get_property('wall')
	offsets=[]
	connections=[]
	pts=[]
	current_offset=0
	eid_trans={}
	pid_trans={}
	for eid in graph.edges():
		cid=graph.source(eid)
		wid=wall[eid]
		pids=list(mesh.borders(1,wid))
		for pid in pids:
			if (cid,pid) not in pid_trans:
				vec=offset_vectors[cid][pid]
				pt=array([position[pid][0],position[pid][1],0.0])+width*vec
				pts.append(pt)
				pid_trans[(cid,pid)]=len(pid_trans)
			connections.append(pid_trans[(cid,pid)])
		current_offset+=2
		offsets.append(current_offset)
		eid_trans[len(eid_trans)]=eid	
		
	connect_string=' '.join([str(val) for val in connections])
	offset_string=' '.join([str(val) for val in offsets])		
	position_string=' '.join([' '.join([str(coord) for coord in point]) for point in pts])
	total_points=len(pid_trans)
	total_cells=len(eid_trans)
	type_string = ' '.join(["3" for num in xrange(total_cells)])
	
	edge_dict={'position_string':position_string,'connect_string':connect_string,\
			'offset_string':offset_string,'type_string':type_string,\
			'n_cells':str(total_cells),'n_points':str(total_points),'cid_trans':eid_trans}

	return wall_dict,cell_dict,edge_dict

def get_wall_strings(db):
	print 'setting wall strings'

	position=db.get_property('position')
	mesh=get_mesh(db)
	pid_trans={}
	pts=[]
	for pid in mesh.wisps(0):
		pts.append(list(position[pid]))
		pid_trans[pid]=len(pid_trans)
	total_points=len(pts)
	position_string=' '.join([' '.join([str(coord) for coord in point]) for point in pts])

	wid_trans={}
	connections=[]
	offsets=[]
	current_offset=0
	for wid in mesh.wisps(2):
		pids = ordered_pids(mesh,wid)
		for pid in pids:
			connections.append(pid_trans[pid])
		current_offset+=len(pids)
		offsets.append(current_offset)
		wid_trans[len(wid_trans)]=wid
	offset_string=' '.join([str(val) for val in offsets])
	connect_string=' '.join([str(val) for val in connections])
	total_cells=len(wid_trans)
	type_string = ' '.join(["7" for num in xrange(total_cells)])		
	wall_dict={'position_string':position_string,'connect_string':connect_string,\
			'offset_string':offset_string,'type_string':type_string,\
			'n_cells':str(total_cells),'n_points':str(total_points),'cid_trans':wid_trans}

	return wall_dict

def get_cell_strings(db):
	print 'setting cell strings'

	position=db.get_property('position')
	mesh=get_mesh(db)
	
	position_string=""
	connect_string=""
	offset_string=""
	
	current_offset=0
	cid_trans={}
	point_count=0
	for cid in mesh.wisps(3):
		pid_trans={}	
		pts=[]
		offsets=[]
		connections=[]
		pid_pos,wids=get_cell_points(db,mesh,cid,position)#{pid:new_pos for all points of cell}+wids needed

		for pid in pid_pos:
			pts.append(list(pid_pos[pid]))
			pid_trans[pid]=point_count
			point_count+=1
			
		for wid in wids:
			xpids=ordered_pids(mesh,wid)
			for xpid in xpids:
				connections.append(pid_trans[xpid])
			current_offset+=len(xpids)
			offsets.append(current_offset)
			cid_trans[len(cid_trans)]=cid

		temp_str=' '.join([' '.join([str(coord) for coord in point]) for point in pts])
		position_string =' '.join([position_string,temp_str])		
		temp_str=' '.join([str(val) for val in connections])
		connect_string =' '.join([connect_string,temp_str])
		temp_str=' '.join([str(val) for val in offsets])
		offset_string =' '.join([offset_string,temp_str])
	total_points=point_count
	total_cells=len(cid_trans)

	type_string = ' '.join(["7" for num in xrange(total_cells)])
		
	cell_dict={'position_string':position_string,'connect_string':connect_string,\
			'offset_string':offset_string,'type_string':type_string,\
			'n_cells':str(total_cells),'n_points':str(total_points),'cid_trans':cid_trans}

	return cell_dict

def get_edge_strings(db):
	print 'setting edge strings'

	position=db.get_property('position')
	mesh=get_mesh(db)
	graph=get_graph(db)
	
	position_string=""
	connect_string=""
	offset_string=""
	
	current_offset=0
	cid_trans={}

	wall_decomp=get_wall_decomp(db)
	point_count=0
	for cid in mesh.wisps(3):
		pid_trans={}
		pts=[]
		offsets=[]
		connections=[]
		pid_pos,wids=get_edge_points(db,mesh,cid,position)#{pid:new_pos for all points of cell}+wids needed

		for pid in pid_pos: 
			pts.append(list(pid_pos[pid]))
			pid_trans[pid]=point_count
			point_count+=1

			
		for wid in wids:
			xpids=ordered_pids(mesh,wid)
			for xpid in xpids:
				connections.append(pid_trans[xpid])
			current_offset+=len(xpids)
			offsets.append(current_offset)
			for eid in wall_decomp[wid]:
				if graph.source(eid)==cid:
					xeid=eid

			cid_trans[len(cid_trans)]=xeid

		temp_str=' '.join([' '.join([str(coord) for coord in point]) for point in pts])
		position_string =' '.join([position_string,temp_str])		
		temp_str=' '.join([str(val) for val in connections])
		connect_string =' '.join([connect_string,temp_str])
		temp_str=' '.join([str(val) for val in offsets])
		offset_string =' '.join([offset_string,temp_str])
	total_points=point_count
	total_cells=len(cid_trans)

	type_string = ' '.join(["7" for num in xrange(total_cells)])
		
	edge_dict={'position_string':position_string,'connect_string':connect_string,\
			'offset_string':offset_string,'type_string':type_string,\
			'n_cells':str(total_cells),'n_points':str(total_points),'cid_trans':cid_trans}

	return edge_dict


def get_cell_points(db,mesh,cid,position,width=0.2):
	offset_vectors=db.get_property("offset_vectors")	
	orig_pid=db.get_property('orig_pid')
	orig_cid=db.get_property('orig_cid')
	pid_pos={}
	wids=[]
	all_pids=[]
	for wid in mesh.borders(3,cid):
		wids.append(wid)
		for lid in mesh.borders(2,wid):
			for pid in mesh.borders(1,lid):
				all_pids.append(pid)
	all_pids=list(set(all_pids))

	pid_pos	={pid:list(array(position[pid])+width*offset_vectors[orig_cid[cid]][orig_pid[pid]]) for pid in all_pids}
	pid_pos={pid:[round(pos[0],4),round(pos[1],4),round(pos[2],4)] for pid,pos in pid_pos.iteritems()}
	zmin=min([pos[2] for pos in pid_pos.itervalues()])
	zmax=max([pos[2] for pos in pid_pos.itervalues()])
	for pid,pos in pid_pos.iteritems():
		if pos[2]==zmin:
			pid_pos[pid]=[pos[0],pos[1],pos[2]+width]
		if pos[2]==zmax:
			pid_pos[pid]=[pos[0],pos[1],pos[2]-width]
	return pid_pos,wids



def get_edge_points(db,mesh,cid,position,width=0.2):
	offset_vectors=db.get_property("offset_vectors")
	orig_pid=db.get_property('orig_pid')
	orig_cid=db.get_property('orig_cid')
	pid_pos={}
	wids=[]
	all_pids=[]
	for wid in mesh.borders(3,cid):
		if mesh.nb_regions(2,wid)>1:
			wids.append(wid)
			for lid in mesh.borders(2,wid):
				for pid in mesh.borders(1,lid):
					all_pids.append(pid)
	all_pids=list(set(all_pids))
	
	pid_pos	={pid:list(array(position[pid])+width*offset_vectors[orig_cid[cid]][orig_pid[pid]]) for pid in all_pids}
	pid_pos={pid:[round(pos[0],4),round(pos[1],4),round(pos[2],4)] for pid,pos in pid_pos.iteritems()}

	zmin=min([pos[2] for pos in pid_pos.itervalues()])
	zmax=max([pos[2] for pos in pid_pos.itervalues()])

	for pid,pos in pid_pos.iteritems():
		if pos[2]==zmin:
			pid_pos[pid]=[pos[0],pos[1],pos[2]+width]
		if pos[2]==zmax:
			pid_pos[pid]=[pos[0],pos[1],pos[2]-width]

	return pid_pos,wids

def get_2d_offset_vectors(mesh,position):

	offset_vectors={}
	for cid in mesh.wisps(2):
		offset_vectors[cid]={}
		pids = ordered_pids(mesh,cid)
		cell_pts = [position[i] for i in pids]
		xvals=[i[0] for i in cell_pts]
		yvals=[i[1] for i in cell_pts]
		xlen=len(cell_pts)
		xy_centre=[sum(xvals)/xlen,sum(yvals)/xlen,0.0]
		N=len(pids)
		for i,pid in enumerate(pids):
			
			pvals = [list(position[pid]),list(position[pids[(i+1)%N]]),list(position[pids[i-1]])]
			for val in pvals:
				if len(val)==2:
					val.append(0.0)
			pvals.append(xy_centre) 
			offset_vectors[cid][pid]=get_offset_vec2(pvals)

	return offset_vectors

