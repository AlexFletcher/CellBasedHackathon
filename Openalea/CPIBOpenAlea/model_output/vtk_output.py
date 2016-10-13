import lxml.etree as ET
from model_utils.db_utilities import get_mesh, get_graph, get_wall_decomp, get_cell_pids
from model_utils.db_geom import *
from sa_oa.container import ordered_pids
import glob
from numpy import array,dot,cross
from scipy.linalg import norm
from model_utils.geom import area


def vtu_series(directory,types):
	"""
	converts series of TissueDB .zip files (named 'iter_*****.zip') in given directory into individual .vtu files then time series .pvd 		format suitable for viewing in paraview. If overwrite = False the function will not write the .vtu files but paraview will expect a 		set of .vtu files with the correct naming convention
	"""


	for vtype in types:
		dblist = glob.glob("%s/%s_iter_*.vtu" % (directory,vtype))
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
			DataSet.set("file","%s_iter_%05d.vtu" % (vtype,iteration))
		final = ET.tostring(TCfile, pretty_print=True)
		text_file = open("%s/%s_time_series.pvd" % (directory,vtype), "w")
		text_file.write('<?xml version="1.0"?>%s' % final)
		text_file.close()


def db_to_vtu(db,directory,savename,types):
	""" converts a TissueDB into a .vtu xml file. requires a output filename"""
	"""
	try:
		Division=db.get_property('Division')
		Growth=db.get_property('Growth')
		if Division or Growth:
			set_vtk_strings(db)
	except KeyError:
		set_vtk_strings(db)
	"""
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

