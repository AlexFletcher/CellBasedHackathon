import lxml.etree as ET
#from model_input.import_tissue import import_tissue
from model_utils.db_utilities import get_mesh, get_graph, get_wall_decomp
from model_utils.db_geom import interior_offset, face_interior_offset, horiz_face, init_ordered_pts, face_interior_offset_3d_vert
from sa_oa.container import ordered_pids
import glob
from numpy import array,dot,cross
from scipy.linalg import norm

def set_vtk_strings(db):
	import time
	start_time=time.time()
	print 'setting vtk strings'
	cell_wall_width=db.get_property('cell_wall_width')
	position=db.get_property('position')
	#wall_decomp=get_wall_decomp(db)
	mesh=get_mesh(db)
	#graph=get_graph(db)

	cell_ps=""
	cell_ts=""
	cell_cs=""
	cell_os=""

	wall_ps=""
	wall_ts=""
	wall_cs=""
	wall_os=""

	edge_ps=""
	edge_ts=""
	edge_cs=""
	edge_os=""

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
	import time
	start=time.time()
	centroids={}
	for cid in mesh.wisps(3):
		cpids=[]
		for fid in mesh.borders(3,cid):
			for lid in mesh.borders(2,fid):
				for pid in mesh.borders(1,lid):
					if pid not in cpids:
						cpids.append(pid)
		

		xvals=[position[i][0] for i in cpids]
		yvals=[position[i][1] for i in cpids]
		zvals=[position[i][2] for i in cpids]
		xlen=len(cpids)
		centroids[cid]=[sum(xvals)/xlen,sum(yvals)/xlen,sum(zvals)/xlen]
	print 'centroid time',time.time()-start



	for wid in mesh.wisps(2):
		pids = ordered_pids(mesh, wid)
		pts = [position[pid] for pid in pids]
		for pt in pts:
			wall_ps += "%s %s %s " % (str(pt[0]),str(pt[1]),str(pt[2]))
			wall_cs += "%s " % str(wall_pp)
			wall_co +=1
			wall_pp +=1
		wall_os+="%s " % str(wall_co)
		wid_trans[wid_pos]=wid
		wid_pos+=1



		xvals=[i[0] for i in pts]
		yvals=[i[1] for i in pts]
		zvals=[i[2] for i in pts]
		xlen=len(pts)
		face_centre=array([sum(xvals)/xlen,sum(yvals)/xlen,sum(zvals)/xlen])
		
		k,offset_vecs,xpts=face_interior_offset_3d_vert(pts)
		
		
		fc_k=k-face_centre

		for cid in mesh.regions(2,wid):
			fc_cc=array(centroids[cid])-face_centre
			fc_cc=fc_cc/norm(fc_cc)
	
			k_direction=cmp(dot(k,fc_cc),0)
			if not horiz_face(wid,db):	
				
				for i,pt in enumerate(xpts):
					svec=array(xpts[(i+1)%len(xpts)])-array(pt)
					svec=svec/norm(svec)
					x=cross(svec,array([0,0,1]))
					new_pt=(cell_wall_width*k_direction)*k/(dot(x,k)) + pt - cell_wall_width*offset_vecs[i]				

					cell_ps += "%s %s %s " % (str(new_pt[0]),str(new_pt[1]),str(new_pt[2]))
					cell_cs += "%s " % str(cell_pp)
					cell_co +=1
					cell_pp +=1
					new_pt=0.5*k*cell_wall_width*k_direction + pt - 0.5*cell_wall_width*offset_vecs[i]
					edge_ps += "%s %s %s " % (str(new_pt[0]),str(new_pt[1]),str(new_pt[2]))
					edge_cs += "%s " % str(edge_pp)
					edge_co +=1
					edge_pp +=1

						
						
					

				cell_os+="%s " % str(cell_co)
				cid_trans[cid_pos]=cid
				cid_pos+=1
				edge_os+="%s " % str(edge_co)
				eid_trans[eid_pos]=1880
				eid_pos+=1					





	for i in range(cid_pos):
		cell_ts +="7 "
	for i in range(wid_pos):
		wall_ts +="7 "
	for i in range(eid_pos):
		edge_ts +="7 "

	db.set_property('vtk_edge_strings',{'position_string':edge_ps,'connect_string':edge_cs,\
			'offset_string':edge_os,'type_string':edge_ts,'n_cells':str(eid_pos),'n_points':str(edge_pp),'cid_trans':eid_trans})
	db.set_description('vtk_edge_strings','description of edges for vtk format')			
	db.set_property('vtk_strings',{'position_string':cell_ps,'connect_string':cell_cs,\
			'offset_string':cell_os,'type_string':cell_ts,'n_cells':str(cid_pos),'n_points':str(cell_pp),'cid_trans':cid_trans})
	db.set_description('vtk_strings','description of cells for vtk format')
	db.set_property('vtk_wall_strings',{'wall_ps':wall_ps,'wall_cs':wall_cs,\
			'wall_os':wall_os,'wall_ts':wall_ts,'wall_nc':str(wid_pos),'wall_np':str(wall_pp),'wid_trans':wid_trans})
	db.set_description('vtk_wall_strings','description of walls for vtk format')
	print 'done'
	print 'time:',time.time()-start_time

def vtu_series(directory):
	"""
	converts series of TissueDB .zip files (named 'iter_*****.zip') in given directory into individual .vtu files then time series .pvd 		format suitable for viewing in paraview. If overwrite = False the function will not write the .vtu files but paraview will expect a 		set of .vtu files with the correct naming convention
	"""

	dblist = glob.glob("%s/iter_*.zip" % directory)
	iters=[]
	for name in dblist:
		iters.append(int(name[-9:-4]))

	TCfile = ET.Element("VTKFile")
	TCfile.set("type","Collection")
	TCfile.set("version","0.1") 
	Collection = ET.SubElement(TCfile, "Collection")
	for iteration in iters:
		DataSet = ET.SubElement(Collection, "DataSet")
		DataSet.set("timestep",str(iteration))
		DataSet.set("group","")
		DataSet.set("part","0")
		DataSet.set("file","iter_%05d.vtu" % iteration)
	final = ET.tostring(TCfile, pretty_print=True)
	text_file = open("%s/time_series.pvd" % directory, "w")
	text_file.write('<?xml version="1.0"?>%s' % final)
	text_file.close()


def zip_to_vtu(db,savename):
	""" converts a TissueDB into a .vtu xml file. requires a output filename"""
	Division=db.get_property('Division')
	Growth=db.get_property('Growth')
	try:
		if Division or Growth:
			set_vtk_strings(db)
	except KeyError:
		set_vtk_strings(db)


	vtk_strings=db.get_property('vtk_strings')
	position_string=vtk_strings['position_string']
	connect_string=vtk_strings['connect_string']
	offset_string=vtk_strings['offset_string']
	type_string=vtk_strings['type_string']
	n_cells=vtk_strings['n_cells']
	n_points=vtk_strings['n_points']
	cid_trans=vtk_strings['cid_trans']
	
	
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
		if ntype=='CELL' or ntype=='cell':
			prop_string=""
			prop=db.get_property(propname)
			for ncid,ocid in cid_trans.iteritems():
				prop_string +="%s " % str(prop[ocid])
			prop_args={"type":"Float32","Name":propname,"NumberOfComponents":"1","format":"ascii"}
			propDA = add_vtu_subelement(CellData, "DataArray",prop_args,prop_string)			
	
	final = ET.tostring(VTKFile, pretty_print=True)

	text_file = open(savename, "w")
	text_file.write('<?xml version="1.0"?>%s' % final)
	text_file.close()


def add_vtu_subelement(parent,name,args,text_string):
	SE = ET.SubElement(parent,name)
	for a1,a2 in args.iteritems():
		SE.set(a1,a2)
	SE.text=text_string
	return SE


def wall_to_vtu(db,savename):
	""" converts a TissueDB into a .vtu xml file. requires a output filename"""
	Division=db.get_property('Division')
	Growth=db.get_property('Growth')
	try:
		if Division or Growth:
			set_vtk_strings(db)
	except KeyError:
		set_vtk_strings(db)


	vtk_strings=db.get_property('vtk_wall_strings')
	position_string=vtk_strings['wall_ps']
	connect_string=vtk_strings['wall_cs']
	offset_string=vtk_strings['wall_os']
	type_string=vtk_strings['wall_ts']
	n_cells=vtk_strings['wall_nc']
	n_points=vtk_strings['wall_np']
	cid_trans=vtk_strings['wid_trans']
	
	
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
		if ntype=='WALL' or ntype=='wall':
			prop_string=""
			prop=db.get_property(propname)
			for ncid,ocid in cid_trans.iteritems():
				prop_string +="%s " % str(prop[ocid])
			prop_args={"type":"Float32","Name":propname,"NumberOfComponents":"1","format":"ascii"}
			propDA = add_vtu_subelement(CellData, "DataArray",prop_args,prop_string)			
	
	final = ET.tostring(VTKFile, pretty_print=True)

	text_file = open(savename, "w")
	text_file.write('<?xml version="1.0"?>%s' % final)
	text_file.close()

def edge_to_vtu(db,savename):
	""" converts a TissueDB into a .vtu xml file. requires a output filename"""
	Division=db.get_property('Division')
	Growth=db.get_property('Growth')
	try:
		if Division or Growth:
			set_vtk_strings(db)
	except KeyError:
		set_vtk_strings(db)


	vtk_strings=db.get_property('vtk_edge_strings')
	position_string=vtk_strings['position_string']
	connect_string=vtk_strings['connect_string']
	offset_string=vtk_strings['offset_string']
	type_string=vtk_strings['type_string']
	n_cells=vtk_strings['n_cells']
	n_points=vtk_strings['n_points']
	cid_trans=vtk_strings['cid_trans']
	
	
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
		if ntype=='EDGE' or ntype=='edge':
			prop_string=""
			prop=db.get_property(propname)
			for ncid,ocid in cid_trans.iteritems():
				prop_string +="%s " % str(prop[ocid])
			prop_args={"type":"Float32","Name":propname,"NumberOfComponents":"1","format":"ascii"}
			propDA = add_vtu_subelement(CellData, "DataArray",prop_args,prop_string)
				
	
	final = ET.tostring(VTKFile, pretty_print=True)

	text_file = open(savename, "w")
	text_file.write('<?xml version="1.0"?>%s' % final)
	text_file.close()
	
