import lxml.etree as ET
#from model_input.import_tissue import import_tissue
from model_utils.db_utilities import get_mesh, get_graph
from model_utils.db_geom import interior_offset, face_interior_offset, horiz_face, init_ordered_pts
from sa_oa.container import ordered_pids
import glob

def set_vtk_strings(db):
	print 'setting vtk strings'
	cell_wall_width=db.get_property('cell_wall_width')
	position=db.get_property('position')
	mesh=get_mesh(db)

	position_string=""
	type_string=""
	connect_string=""
	offset_string=""
	n_cells=0
	n_points=0

	wall_ps=""
	wall_ts=""
	wall_cs=""
	wall_os=""
	wall_nc=0
	wall_np=0

	if mesh.degree()==2:
		n_cells=len(list(mesh.wisps(mesh.degree())))

		
		for i in range(n_cells):
			type_string +="7 "

		pid_trans={}
		
		interior_points = interior_offset(db,cell_wall_width)
		str_pos=0
		for cid,pts in interior_points.iteritems():
			pid_trans[cid]=[]
			for pt in pts:
				position_string += "%s %s %s " % (str(pt[0]),str(pt[1]),str(pt[2]))
				pid_trans[cid].append(str_pos)
				str_pos+=1
		n_points=str_pos


		current_offset=0
		cid_trans={}
		cid_pos=0
		for cid,pts in interior_points.iteritems():
			current_offset +=len(pts)
			pidpos=0
			for pt in pts:
				connect_string += "%s " % str(pid_trans[cid][pidpos])
				pidpos+=1
			offset_string+="%s " % str(current_offset)
			cid_trans[cid_pos]=cid
			cid_pos+=1
	elif mesh.degree()==3:

		cell_faces_pids = {}
		zvals={}
		interior_points = {}

		for cid in mesh.wisps(3):
			faces = list(mesh.borders(3,cid))
			zvals[cid]=[]
			cell_faces_pids[cid]={}
			for fid in faces:
				cell_faces_pids[cid][fid]=ordered_pids(mesh, fid)
				interior_points[fid]=[]
				if horiz_face(fid,db):					
					interior_points[fid] = face_interior_offset(fid,db,cell_wall_width)
					zvals[cid].append(interior_points[fid][0][2])
			zvals[cid]=sorted(zvals[cid])

		
		pt_pos=0
		current_offset=0
		wall_pp=0
		wall_co=0
		cid_trans={}
		cid_pos=0
		wid_trans={}
		wid_pos=0

		for cid,fids in cell_faces_pids.iteritems():


			for fid,pids in fids.iteritems():
				if len(interior_points[fid])>0:
					
					
					#add top
					if interior_points[fid][0][2]==zvals[cid][0]:			
						for pt in interior_points[fid]:
							ps="%s %s %s " % (str(pt[0]),str(pt[1]),\
									str(pt[2]+cell_wall_width))#move in z direction
							#cell
							position_string += ps
							connect_string += "%s " % str(pt_pos)
							current_offset +=1
							pt_pos+=1
							#wall
							#wall_ps += ps
							#wall_cs += "%s " % str(wall_pp)
							#wall_co +=1
							#wall_pp +=1
						#wall_os+="%s " % str(wall_co)
						#wid_trans[wid_pos]=fid
						#wid_pos+=1
							
						
						opos=init_ordered_pts([position[pid] for pid in cell_faces_pids[cid][fid]])
						for pt in opos:
							wall_ps += "%s %s %s " % (str(pt[0]),str(pt[1]),str(pt[2]))
							wall_cs += "%s " % str(wall_pp)
							wall_pp +=1
							wall_co +=1
						
						wall_os+="%s " % str(wall_co)
						wid_trans[wid_pos]=fid
						wid_pos+=1						
					
						#for i,pt in enumerate(opos):
						#	N=len(opos)
						#	wall_cs += "%s " % str(wall_pp-2*N + i%N)
						#	wall_cs += "%s " % str(wall_pp-2*N + (i+1)%N)
						#	wall_cs += "%s " % str(wall_pp-N + (i+1)%N)
						#	wall_cs += "%s " % str(wall_pp-N + i%N)
						#	wall_co+=4
						#	wall_os+="%s " % str(wall_co)
						#	wid_trans[wid_pos]=fid
						#	wid_pos+=1
					#add bottom		
					elif interior_points[fid][0][2]==zvals[cid][1]:
						for pt in interior_points[fid]:
							ps="%s %s %s " % (str(pt[0]),str(pt[1]),\
									str(pt[2]-cell_wall_width))#move in z direction
							#cell
							position_string += ps
							connect_string += "%s " % str(pt_pos)
							current_offset +=1
							pt_pos+=1
							#wall
							#wall_ps += ps
							#wall_cs += "%s " % str(wall_pp)
							#wall_co +=1
							#wall_pp +=1
						#wall_os+="%s " % str(wall_co)
						#wid_trans[wid_pos]=fid
						#wid_pos+=1
							
						
						opos=init_ordered_pts([position[pid] for pid in cell_faces_pids[cid][fid]])
						for pt in opos:
							wall_ps += "%s %s %s " % (str(pt[0]),str(pt[1]),str(pt[2]))
							wall_cs += "%s " % str(wall_pp)
							wall_pp +=1
							wall_co +=1
						
						wall_os+="%s " % str(wall_co)
						wid_trans[wid_pos]=fid
						wid_pos+=1						
					
						#for i,pt in enumerate(opos):
						#	N=len(opos)
						#	wall_cs += "%s " % str(wall_pp-2*N + i%N)
						#	wall_cs += "%s " % str(wall_pp-2*N + (i+1)%N)
						#	wall_cs += "%s " % str(wall_pp-N + (i+1)%N)
						#	wall_cs += "%s " % str(wall_pp-N + i%N)
						#	wall_co+=4
						#	wall_os+="%s " % str(wall_co)
						#	wid_trans[wid_pos]=fid
						#	wid_pos+=1


					offset_string+="%s " % str(current_offset)
					cid_trans[cid_pos]=cid
					cid_pos+=1

					wall_count=len(interior_points[fid])

				
			#side wall rectangles
			for pt in range(wall_count):
				connect_string += "%s " % str(pt_pos-2*wall_count + pt%wall_count)
				connect_string += "%s " % str(pt_pos-2*wall_count + (pt+1)%wall_count)
				connect_string += "%s " % str(pt_pos-wall_count + (pt+1)%wall_count)
				connect_string += "%s " % str(pt_pos-wall_count + pt%wall_count)
				current_offset += 4
				offset_string+="%s " % str(current_offset)
				cid_trans[cid_pos]=cid
				cid_pos+=1
				wall_cs += "%s " % str(wall_pp-2*wall_count + pt%wall_count)
				wall_cs += "%s " % str(wall_pp-2*wall_count + (pt+1)%wall_count)
				wall_cs += "%s " % str(wall_pp-wall_count + (pt+1)%wall_count)
				wall_cs += "%s " % str(wall_pp-wall_count + pt%wall_count)
				wall_co += 4
				wall_os+="%s " % str(wall_co)
				wid_trans[wid_pos]=fid
				wid_pos+=1
					

		n_cells=cid_pos
		n_points=pt_pos		
		wall_nc=wid_pos
		wall_np=wall_pp

		for i in range(n_cells):
			type_string +="7 "
		for i in range(wall_nc):
			wall_ts +="7 "





			
	db.set_property('vtk_strings',{'position_string':position_string,'connect_string':connect_string,\
			'offset_string':offset_string,'type_string':type_string,'n_cells':str(n_cells),'n_points':str(n_points),'cid_trans':cid_trans})
	db.set_description('vtk_strings','description of cells for vtk format')
	db.set_property('vtk_wall_strings',{'wall_ps':wall_ps,'wall_cs':wall_cs,\
			'wall_os':wall_os,'wall_ts':wall_ts,'wall_nc':str(wall_nc),'wall_np':str(wall_np),'wid_trans':wid_trans})
	db.set_description('vtk_wall_strings','description of walls for vtk format')
	print 'done'


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
	
