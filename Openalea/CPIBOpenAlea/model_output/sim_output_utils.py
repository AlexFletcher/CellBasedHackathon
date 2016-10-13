import glob
from model_input.import_tissue import import_tissue
import os,sys
from model_output.pylab_plot import plot_pylab_figure
from sa_oa.scheduler import Scheduler,Task
from model_input.import_tissue import import_tissue
from model_utils.db_utilities import get_parameters, set_parameters, get_graph, get_mesh, get_wall_decomp
from model_utils.celltissue_util import get_property
from model_utils.celltissue_util import def_property
from sa_oa.celltissue import TissueDB
def csv_from_zips(directory):
	import csv

	start_db=TissueDB()
	start_db.read('%s/iter_%05d.zip' % (directory,0))
	mesh=get_mesh(start_db)
	species_desc=start_db.get_property('species_desc')
	#if mesh.degree()==3:
	#	vindex=start_db.get_property('vindex')
	times=[]
	tcprops={}
	for prop_name in species_desc.iterkeys():
		tcprops[prop_name]={}		
		prop=start_db.get_property(prop_name)
		for cid,value in prop.iteritems():
			tcprops[prop_name][cid]={}

	dblist = glob.glob("%s/iter_*.zip" % directory)

	iters=[]
	for name in dblist:
		iters.append(int(name[-9:-4]))

	for i in range (0,max(iters)+1):
		new_db=TissueDB()
		new_db.read('%s/iter_%05d.zip' % (directory,i))
		times.append(new_db.get_property('time'))

		for prop_name in species_desc.iterkeys():	
			prop=new_db.get_property(prop_name)
			for cid,value in prop.iteritems():
				tcprops[prop_name][cid][new_db.get_property('time')]=value


	filename='%s/sim_results.csv' % directory
	spreadsheet = csv.writer(open(filename, 'wb'), delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
	cell_type=start_db.get_property('cell_type')
	
	headrow=['prop_name']
	headrow.append('cell_type')
	#if mesh.degree()==3:
	#	headrow.append('vindex')
	headrow.append('id')
	for tp in times:
	    headrow.append(tp)
	spreadsheet.writerow(headrow)


	graph=get_graph(start_db)
	wd = get_wall_decomp(start_db)

	for prop_name,cdict in tcprops.iteritems():
	    for cid,tcourse in cdict.iteritems():
		newrow = [prop_name]
		if species_desc[prop_name]=='cell' or species_desc[prop_name]=='CELL':
			newrow.append(cell_type[cid])
			#if mesh.degree()==3:
			#	newrow.append(vindex[cid])
		elif species_desc[prop_name]=='edge' or species_desc[prop_name]=='EDGE':
			newrow.append('%s-%05d-%05d' % (species_desc[prop_name],graph.source(cid),graph.target(cid)))
			#if mesh.degree()==3:
			#	newrow.append('%02d-%02d' % (vindex[graph.source(cid)],vindex[graph.target(cid)]))
		elif species_desc[prop_name]=='wall' or species_desc[prop_name]=='WALL':
			if cid in wd.iterkeys():
				cella=graph.source(wd[cid][0])
				cellb=graph.source(wd[cid][1])				
				newrow.append('%s-%05d-%05d' % (species_desc[prop_name],cella,cellb))
				#if mesh.degree()==3:
				#	newrow.append('%02d-%02d' % (vindex[cella],vindex[cellb]))
			else:
				
				newrow.append('%s-%s' % (species_desc[prop_name],'outer'))
				#if mesh.degree()==3:
				#	cella=list(mesh.regions(2,cid))[0]
				#	newrow.append('%02d-%s' % (vindex[cella],'outer'))
		newrow.append(cid)
		for tp in times:
			newrow.append(tcourse[tp])
		spreadsheet.writerow(newrow)


def draw_from_zips(directory,view_names,cmap_lims=None):
	
	if cmap_lims is None:		
		cmap_lims=[]
		for name in view_names:
			cmap_lims.append(get_range(directory,name))

	

	db=import_tissue('%s/iter_00000.zip' % directory)
	plot_pylab_figure(db,0,directory,view_names,cmap_lims)

	dblist = glob.glob("%s/iter_*.zip" % directory)

	iters=[]
	for name in dblist:
		iters.append(int(name[-9:-4]))

	for i in range (1,max(iters)+1):
		new_db=import_tissue('%s/iter_%05d.zip' % (directory,i))
		#for prop in db.properties():
		#	nprop=new_db.get_property(prop)
		#	eprop=self.db.get_property(prop)
		#	if type(nprop)==dict:
		#		for cid,val in nprop.iteritems():
		#			eprop[cid]=val	
		plot_pylab_figure(new_db,i,directory,view_names,cmap_lims)

def draw_fluxes(directory,maxlength):

	db=import_tissue('%s/iter_00000.zip' % directory)
	from model_output.flux_plot import flux_plot
	flux_plot(db,0,directory,maxlength)

	dblist = glob.glob("%s/iter_*.zip" % directory)

	iters=[]
	for name in dblist:
		iters.append(int(name[-9:-4]))

	for i in range (1,max(iters)+1):
		new_db=import_tissue('%s/iter_%05d.zip' % (directory,i))
		flux_plot(new_db,i,directory,maxlength)


def make_ffmpeg(directory):
	command='ffmpeg -f image2 -r 1.0 -i %s/cmap_%%05d.png -s 1280x960 %s/%s.mp4' % (directory,directory,directory)
	os.system(command)



def get_range(directory,prop_name):
	from math import log10, floor
	maxval=(0.0,0.0)		

	dblist = glob.glob("%s/iter_*.zip" % directory)

	iters=[]
	for name in dblist:
		iters.append(int(name[-9:-4]))

	for i in range (0,max(iters)+1):
		new_db=import_tissue('%s/iter_%05d.zip' % (directory,i))

		prop=new_db.get_property(prop_name)
		for cid,value in prop.iteritems():
			if value>maxval[1]:
				maxval=(0,value)
	if maxval[1]==0.0:
		maxval=(0.0,1.0)
	else:
		maxval=(0.0,round(maxval[1], -int(floor(log10(maxval[1])))))

	return maxval

def get_range_db(db,prop_name):
	from math import log10, floor
	maxval=(0.0,0.0)		


	prop=db.get_property(prop_name)
	for cid,value in prop.iteritems():
		if value>maxval[1]:
			maxval=(0,value)
	if maxval[1]==0.0:
		maxval=(0.0,1.0)
	else:
		maxval=(0.0,round(maxval[1], -int(floor(log10(maxval[1])))))

	return maxval

def writeParameters(db,directory):
	# Write parameters of db to filename
	plist = glob.glob('%s/parameters*.dat' % directory)
	fname='%s/parameters%02d.dat' % (directory,len(plist))
	f = open(fname,"w")
	params = db.get_property('parameters')
	for param_name, prop in params.iteritems():
		f.write(param_name + '\n')
		if isinstance(prop, dict):
	    		for key, value in prop.iteritems():
				f.write('    ' + str(key) + ' : ' + str(value) + '\n')
		else:
	    		f.write('    '+ str(prop) + '\n')
	f.write('\n')


def convert_from_trimmed(directory,original):
	new_dir='xpts_%s' % directory
	os.system('mkdir %s' % new_dir)	
	db=import_tissue('%s/iter_00000.zip' % directory)
	cell_trans=db.get_property('cell_trans')
	wall_trans=db.get_property('wall_trans')
	orig_db=import_tissue('%s.zip' % original)

	species_desc=db.get_property('species_desc')
	orig_db.set_property('species_desc',species_desc)
	orig_db.set_description('species_desc','species_desc')
	for propname,proptype in species_desc.iteritems():
		transprop={}
		prop=db.get_property(propname)
		if proptype=='cell' or proptype=='CELL':
			for ncid,ocid in cell_trans.iteritems():
				transprop[ocid]=prop[ncid]
		if proptype=='wall' or proptype=='WALL':
			for nwid,owids in wall_trans.iteritems():
				for owid in owids:
					transprop[owid]=prop[nwid]
		orig_db.set_property(propname,transprop)
		orig_db.set_description(propname, propname)
		

	orig_db.write('%s/iter_00000.zip' % new_dir)

	dblist = glob.glob("%s/iter_*.zip" % directory)

	iters=[]
	for name in dblist:
		iters.append(int(name[-9:-4]))

	for i in range (1,max(iters)+1):
		new_db=import_tissue('%s/iter_%05d.zip' % (directory,i))
		for propname,proptype in species_desc.iteritems():
			transprop={}
			prop=new_db.get_property(propname)
			if proptype=='cell' or proptype=='CELL':
				for ncid,ocid in cell_trans.iteritems():
					transprop[ocid]=prop[ncid]
			if proptype=='wall' or proptype=='WALL':
				for nwid,owids in wall_trans.iteritems():
					for owid in owids:
						transprop[owid]=prop[nwid]
			orig_db.set_property(propname,transprop)
			orig_db.set_description(propname, propname)
			

		orig_db.write('%s/iter_%05d.zip' % (new_dir,i))

