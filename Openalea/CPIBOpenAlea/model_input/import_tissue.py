from model_utils.db_utilities import get_mesh,get_graph
from model_utils.celltissue_util import def_property
from import_xml import import_xml
from sa_oa.celltissue import TissueDB
from sa_oa.tissueshape import centroid
from import_utils import *
import sys
from trim_walls import *

def import_tissue(fname):

	#print 'importing tissue',fname

	if fname.split('.')[-1]=='zip':
		db=TissueDB()	
		db.read(fname)

	elif fname.split('.')[-1]=='xml':
		db = import_xml(fname)
		set_vtk_strings(db)
		db.write("%s.zip" % fname.split('.')[0])
		dbt = trim_walls(db)
		set_vtk_strings(dbt)
		dbt.write("%s_trim.zip" % fname.split('.')[0])
	else:
		print "invalid filename or extension '%s'" % fname.split('.')[-1]
		sys.exit()
		

	if 'scale_px_to_micron' not in db.properties():
		db.set_property('scale_px_to_micron',(1.0,'microns'))
		db.set_description('scale_px_to_micron','scale: pixels to microns')

	
	scale = db.get_property('scale_px_to_micron')
	print 'scale: 1 pixel =',scale[0],scale[1]

	position=db.get_property('position')
	xspan=scale[0]*(max([pos[0] for pos in position.itervalues()])-min([pos[0] for pos in position.itervalues()]))
	yspan=scale[0]*(max([pos[1] for pos in position.itervalues()])-min([pos[1] for pos in position.itervalues()]))
	print 'x span =',xspan,'um'
	print 'y span =',yspan,'um'

	"""
	# cell type codes - may need to edit this
	cell_strings = {0: 'undefined',
		   2: 'epiderm',
                   4: 'cortex',
                   3: 'endoderm',
                   5: 'vasculature',
                   6: 'LRC1',
                   7: 'LRC2',
                   8: 'LRC3',
		   9: 'LRC4',
                   17: 'QC',
                   10: 'initials',
                   11: 'S1',
                   12: 'S2',
                   13: 'S3',
                   14: 'S4',
                   15: 'S5',
		   16: 'pericycle',
		   19: 'protoxylem',
		   20: 'metaxylem',
		   21: 'xylem pole pericycle',
		   22: 'procambium',
		   23: 'phloem',
		   24: 'meta sieve element',
		   25: 'peri_adj_pc'}

	db.set_property('cell_type_strings',cell_strings)
	db.set_description('cell_type_strings','cell types by code number')
	"""


	if 'vindex' not in db.properties():
		def_property(db,'vindex',0,'CELL',"config","")

	if 'border' not in db.properties():
		print 'border reset'
		def_property(db,'border',0,'CELL',"config","")

	if 'time' not in db.properties():
    		db.set_property('time', 0)
   		db.set_description('time', 'current simulation time')

	if 'iteration' not in db.properties():
    		db.set_property('iteration', 0)
   		db.set_description('iteration', 'current iteration')

	if 'parameters' in db.properties():
		parameters = db.get_property('parameters')
		parameters.clear()

	if 'species_desc' not in db.properties():
		species_desc={'V':'cell','S':'wall','cell_type':'cell'}
		db.set_property('species_desc', species_desc)
    		db.set_description('species_desc', 'type of property, i.e. cell,edge or wall')

	if 'divided_props' not in db.properties():
		divided_props={'cell':['cell_type','orig_cid'],'wall':[],'edge':[],'point':[]}
		db.set_property('divided_props', divided_props)
    		db.set_description('divided_props', 'properties to reassign on cell division')

	if 'diluted_props' not in db.properties():
		diluted_props={}
		db.set_property('diluted_props', diluted_props)
    		db.set_description('diluted_props', 'properties to dilute as cell grows')

	if 'vtk_strings' not in db.properties():
		set_vtk_strings(db)

	if 'cell_centres' not in db.properties():
		db,centres = def_property(db,'cell_centres',0,'CELL',"config","")
		mesh = get_mesh(db)
		for cid in centres.iterkeys():
			centres[cid] = centroid (mesh, position, mesh.degree(), cid)	

	if 'edge_vindex' not in db.properties():
		db,edge_vindex = def_property(db,'edge_vindex',0,'EDGE',"config","")
		vindex=db.get_property('vindex')
		graph = get_graph(db)
		for eid in graph.edges():
			sid=graph.source(eid)


             
	return db


