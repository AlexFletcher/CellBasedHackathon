import lxml.etree as ET
from model_utils.db_utilities import get_mesh, get_graph, get_wall_decomp, get_cell_pids
from model_utils.db_geom import *
from sa_oa.container import ordered_pids
import glob
from numpy import array,dot,cross
from scipy.linalg import norm
from model_utils.geom import area

def set_vtk_strings(db):
	import time
	start = time.time()
	print 'setting vtk strings'
	master_count=0
	cell_wall_width=db.get_property('cell_wall_width')
	position=db.get_property('position')
	wall_decomp=get_wall_decomp(db)
	mesh=get_mesh(db)
	graph=get_graph(db)

	cell_ps=[]
	cell_ts=[]
	cell_cs=[]
	cell_os=[]

	wall_ps=[]
	wall_ts=[]
	wall_cs=[]
	wall_os=[]

	edge_ps=[]
	edge_ts=[]
	edge_cs=[]
	edge_os=[]

	wid_trans={}
	wid_pos=0
	wall_pp=0
	wall_co=0
	
	cell_pp=0
	cell_co=0		
	cid_trans={}
	cid_pos=0

	eid_trans={}
	eid_pos=0
	edge_pp=0
	edge_co=0

	cell_pids=get_cell_pids(mesh)
	moved_cid_pids={}
	moved_ecid_pids={}
	wall_pid_refs={}
	for wid in mesh.wisps(2):
		pids = ordered_pids(mesh, wid)

		for xpid in pids:
			if xpid not in wall_pid_refs.iterkeys():
				pt = position[xpid]
				wall_ps.append("%s %s %s " % (str(pt[0]),str(pt[1]),str(pt[2])))
				wall_pid_refs[xpid]=wall_pp
				wall_cs.append("%s " % str(wall_pp))
				wall_pp +=1
			else:
				wall_cs.append("%s " % str(wall_pid_refs[xpid]))
			wall_co +=1
				
		wall_os.append("%s " % str(wall_co))
		wid_trans[wid_pos]=wid
		wid_pos+=1

		if len(list(mesh.regions(2,wid)))==1:
			outer_wall=True
		else:
			outer_wall=False

		for cid in mesh.regions(2,wid):


			for i,pid in enumerate(pids):
				p_pids=[]
				for lid in mesh.regions(0,pid):
					for xpid in mesh.borders(1,lid):
						p_pids.append(xpid)
				p_pids=list(set(p_pids))
				for xpid in p_pids:
					if xpid not in cell_pids[cid]:
						p_pids.remove(xpid)
				
				for xpid in p_pids:
					if position[xpid][2]!=position[pid][2]:
						nv=xpid
				p_pids.remove(nv)
				p_pids.append(nv)
				p_pids.remove(pid)
				p_pids.insert(0,pid)
				


				o_pt = array(position[pid])
				if (cid,pid) not in moved_cid_pids.iterkeys():

					vec = get_offset_vec([position[xpid] for xpid in p_pids])
					new_pt = o_pt + cell_wall_width*vec
					cell_ps.append("%s %s %s " % (str(new_pt[0]),str(new_pt[1]),str(new_pt[2])))
					master_count+=1
					if master_count%500==0:
						print 'count',master_count
					cell_cs.append("%s " % str(cell_pp))
					#moved_cid_pids[(cid,pid)]=cell_pp
					cell_pp +=1

					if not outer_wall:
						new_pt= o_pt + 0.5*cell_wall_width*vec
						edge_ps.append("%s %s %s " % (str(new_pt[0]),str(new_pt[1]),str(new_pt[2])))
						#moved_ecid_pids[(cid,pid)]=edge_pp
						edge_cs.append("%s " % str(edge_pp))					
						edge_co +=1
						edge_pp +=1

					else:
						outer_point=True
						for xlid in mesh.regions(0,pid):
							for xwid in mesh.regions(1,xlid):
								if len(list(mesh.regions(2,xwid)))>1:
									outer_point=False
						if not outer_point:
							new_pt= o_pt + 0.5*cell_wall_width*vec
							edge_ps.append("%s %s %s " % (str(new_pt[0]),str(new_pt[1]),str(new_pt[2])))
							#moved_ecid_pids[(cid,pid)]=edge_pp
							edge_pp +=1


				else:
					cell_cs.append("%s " % str(moved_cid_pids[(cid,pid)]))
					if not outer_wall:				
						edge_cs.append("%s " % str(moved_ecid_pids[(cid,pid)]))
						edge_co +=1
				cell_co +=1
				


			cell_os.append("%s " % str(cell_co))
			cid_trans[cid_pos]=cid
			cid_pos+=1
			if not outer_wall:
				edge_os.append("%s " % str(edge_co))
				eid1,eid2=wall_decomp[wid]
				if graph.source(eid1)==cid:
					eid=eid1
				else:
					eid=eid2

				eid_trans[eid_pos]=eid
				eid_pos+=1					

	cell_ts = ''.join(["7 " for num in xrange(cid_pos)])
	wall_ts = ''.join(["7 " for num in xrange(wid_pos)])
	edge_ts = ''.join(["7 " for num in xrange(eid_pos)])
	scell_ps=''.join([val for val in cell_ps])
	scell_ts=''.join([val for val in cell_ts])
	scell_cs=''.join([val for val in cell_cs])
	scell_os=''.join([val for val in cell_os])

	swall_ps=''.join([val for val in wall_ps])
	swall_ts=''.join([val for val in wall_ts])
	swall_cs=''.join([val for val in wall_cs])
	swall_os=''.join([val for val in wall_os])

	sedge_ps=''.join([val for val in edge_ps])
	sedge_ts=''.join([val for val in edge_ts])
	sedge_cs=''.join([val for val in edge_cs])
	sedge_os=''.join([val for val in edge_os])	
	
	edge_dict={'position_string':sedge_ps,'connect_string':sedge_cs,\
			'offset_string':sedge_os,'type_string':sedge_ts,'n_cells':str(eid_pos),'n_points':str(edge_pp),'cid_trans':eid_trans}
	wall_dict={'position_string':swall_ps,'connect_string':swall_cs,\
			'offset_string':swall_os,'type_string':swall_ts,'n_cells':str(wid_pos),'n_points':str(wall_pp),'cid_trans':wid_trans}			
	cell_dict={'position_string':scell_ps,'connect_string':scell_cs,\
			'offset_string':scell_os,'type_string':scell_ts,'n_cells':str(cid_pos),'n_points':str(cell_pp),'cid_trans':cid_trans}

	db.set_property('vtk_strings',{'cell':cell_dict,'wall':wall_dict,'edge':edge_dict})
	db.set_description('vtk_strings','description of cells, walls, edges for vtk format')
	print time.time()-start
	print 'done'


def vtu_series(directory,types):
	"""
	converts series of TissueDB .zip files (named 'iter_*****.zip') in given directory into individual .vtu files then time series .pvd 		format suitable for viewing in paraview. If overwrite = False the function will not write the .vtu files but paraview will expect a 		set of .vtu files with the correct naming convention
	"""

	dblist = glob.glob("%s/iter_*.zip" % directory)
	iters=[]
	for name in dblist:
		iters.append(int(name[-9:-4]))
	for vtype in types:
		TCfile = ET.Element("VTKFile")
		TCfile.set("type","Collection")
		TCfile.set("version","0.1") 
		Collection = ET.SubElement(TCfile, "Collection")
		for iteration in iters:
			DataSet = ET.SubElement(Collection, "DataSet")
			DataSet.set("timestep",str(iteration))
			DataSet.set("group","")
			DataSet.set("part","0")
			DataSet.set("file","%s_iter_%05d.vtu" % (vtype,iteration))
		final = ET.tostring(TCfile, pretty_print=True)
		text_file = open("%s/%s_time_series.pvd" % (directory,vtype), "w")
		text_file.write('<?xml version="1.0"?>%s' % final)
		text_file.close()


def db_to_vtu(db,directory,savename,types):
	""" converts a TissueDB into a .vtu xml file. requires a output filename"""
	Division=db.get_property('Division')
	Growth=db.get_property('Growth')
	try:
		if Division or Growth:
			set_vtk_strings(db)
	except KeyError:
		set_vtk_strings(db)

	for vtype in types:
		vtk_strings=db.get_property('vtk_strings')
		vtk_string=vtk_strings[vtype]
		position_string=vtk_string['position_string']
		connect_string=vtk_string['connect_string']
		offset_string=vtk_string['offset_string']
		type_string=vtk_string['type_string']
		n_cells=vtk_string['n_cells']
		n_points=vtk_string['n_points']
		cid_trans=vtk_string['cid_trans']
	
	
		VTKFile = ET.Element("VTKFile")
		VTKFile.set("type","UnstructuredGrid")
		VTKFile.set("type","UnstructuredGrid")
		VTKFile.set("version","0.1") 
		VTKFile.set("byte_order","LittleEndian")

		UG = ET.SubElement(VTKFile, "UnstructuredGrid")

		piece_args = {"NumberOfPoints":n_points,"NumberOfCells":n_cells}
		Piece = add_vtu_subelement(UG,"Piece",piece_args,"")

		Points = ET.SubElement(Piece, "Points")

		position_args = {"type":"Float32","Name":"position","NumberOfComponents":"3","format":"ascii"}
		PosDA = add_vtu_subelement(Points,"DataArray",position_args,position_string)

		Cells = ET.SubElement(Piece, "Cells")

		connect_args = {"type":"Int32","Name":"connectivity","NumberOfComponents":"1","format":"ascii"}
		ConDA = add_vtu_subelement(Cells,"DataArray",connect_args,connect_string)

		off_args = {"type":"Int32","Name":"offsets","NumberOfComponents":"1","format":"ascii"}
		OffDA = add_vtu_subelement(Cells,"DataArray",off_args,offset_string)

		type_args = {"type":"UInt8","Name":"types","NumberOfComponents":"1","format":"ascii"}
		TypeDA = add_vtu_subelement(Cells,"DataArray",type_args,type_string)


		CellData = ET.SubElement(Piece, "CellData")
	
		species_desc=db.get_property('species_desc')
		for propname,ntype in species_desc.iteritems():
			if ntype==vtype.upper() or ntype==vtype:
				prop_string=""
				prop=db.get_property(propname)
				for ncid,ocid in cid_trans.iteritems():
					prop_string +="%s " % str(prop[ocid])
				prop_args={"type":"Float32","Name":propname,"NumberOfComponents":"1","format":"ascii"}
				propDA = add_vtu_subelement(CellData, "DataArray",prop_args,prop_string)			
	
		final = ET.tostring(VTKFile, pretty_print=True)

		text_file = open("%s/%s_%s" % (directory,vtype,savename), "w")
		text_file.write('<?xml version="1.0"?>%s' % final)
		text_file.close()


def add_vtu_subelement(parent,name,args,text_string):
	SE = ET.SubElement(parent,name)
	for a1,a2 in args.iteritems():
		SE.set(a1,a2)
	SE.text=text_string
	return SE

