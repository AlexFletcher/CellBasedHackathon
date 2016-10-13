import model_output.sim_output_utils as out
from model_utils.db_utilities import get_parameters, set_parameters
from model_utils.celltissue_util import def_property, get_property
from model_structure.genenetwork import CombinedModel
from model_structure.genenetwork import GeneNetwork
from model_input.import_tissue import import_tissue

from sa_oa.scheduler import Scheduler,Task
from sa_oa.celltissue import TissueDB

from model_output.vtk_output import *

import os,sys
from sa_oa.container import Topomesh,graph

import time


#from model_structure.parameters import writeParameters

class Sim_dg_unif(object):
	def __init__(self,tissuename,model_list):
		
		self.sch = Scheduler()

		#default values
		#self.savedb_from_viewer=True
		self.initial_values={}
		self.timestep=1
		self.fixed_dict={}
		self.param_values={}
		
		self.tissuename=tissuename.replace('/','_').split('.')[0]
		self.db=import_tissue(tissuename)



		set_parameters(self.db, 'initial_values', self.initial_values)
		set_parameters(self.db, 'timestep', self.timestep)
		set_parameters(self.db, 'fixed', self.fixed_dict)
		if type(model_list)==list:
			mlist=[model() for model in list(model_list)]
		else:
			mlist=[model_list()]
		self.combo_model=CombinedModel(mlist)
		self.vtk_output_types=['cell','wall','edge']


		############
		
		names=self.combo_model.get_model_names()
		cnames=''
		for name in names:
			cnames='%s_%s' % (cnames,name)
		self.set_sim_name(cnames)

		self.gene_network=GeneNetwork(self.db, self.combo_model)
		self.sch.register(Task(self.gene_network.step, 1, 8, "genenetwork") )

		from gr_unif import Gr_unif
		self.gr=Gr_unif(self.db)
		self.sch.register(Task(self.gr.step, 1, 7, "growth") )

		from dv_unif import Dv_unif
		self.dv=Dv_unif(self.db)
		self.sch.register(Task(self.dv.step, 1, 6, "division") )

		def update_iteration():
			t=self.db.get_property('time')
			t+=get_parameters(self.db, 'timestep')
			self.db.set_property('time', t)

			iteration=self.db.get_property('iteration')
			iteration+=1
			self.db.set_property('iteration', iteration)		


		self.sch.register(Task(update_iteration, 1, 2, "add iteration"))


	#######run simulations#############
	
	def run_simulation(self,timesteps):

		os.system('mkdir %s' % self.dir_name)		
		out.writeParameters(self.db,self.dir_name)
		#self.db.write('%s/iter_%05d.zip' % (self.dir_name,0))
		#db_to_vtu(self.db,self.dir_name,'iter_%05d.vtu' % 0,self.vtk_output_types)
		if type(timesteps)==int:
			ts_list=range(timesteps+1)
		else:
			ts_list=range(timesteps[0],timesteps[1]+1)
		steptime=self.db.get_property('steptime')
		import time
		for i in ts_list:
			print 'iteration',i,'of',ts_list[-1]
			self.db.write('%s/iter_%05d.zip' % (self.dir_name,i))
			stime=time.time()
			db_to_vtu(self.db,self.dir_name,'iter_%05d.vtu' % i,self.vtk_output_types)
			next ( self.sch . run ())
			steptime.append(time.time()-stime)
		print 'writing vtu time series'
		vtu_series(self.dir_name,self.vtk_output_types)
		print 'done'
		print 'NOT writing .csv file'	
		#self.csv_from_zips()
		print 'done'


	#########change settings ###############

	def set_sim_name(self,name):
		self.simname=name
		self.dir_name= '%s%s' % (self.tissuename,self.simname)

	def get_db(self):
		return self.db

	def set_initial_values(self,ivs):
		for prop,vals in ivs.iteritems():
			self.initial_values[prop]=vals
		set_parameters(self.db, 'initial_values', self.initial_values)

	def set_timestep(self,tstep):
		self.timestep=tstep
		set_parameters(self.db, 'timestep', self.timestep)

	def set_fixed_dict(self,fixed_dict):
		self.fixed_dict=fixed_dict
		set_parameters(self.db, 'fixed', self.fixed_dict)

	def set_param_values(self,param_values):
		self.param_values=param_values

		names = self.combo_model.get_model_names()
		for name in names:
			p = get_parameters(self.db, name)
			for par,value in self.param_values.iteritems():
				if par in p.iterkeys():
					p[par]=value
					print 'par',par,'val',value,'overwritten'
			set_parameters(self.db, name, p)

	#def set_save_db_from_viewer(self,savedb):
	#	self.savedb_from_viewer=savedb

	def set_scale(self,px_to_micron):
		self.db.set_property('scale_px_to_micron',(px_to_micron,'microns'))
		self.db.set_description('scale_px_to_micron','scale: pixels to microns')


	def set_cell_wall_width(self,cell_wall_width):
		self.db.set_property('cell_wall_width',cell_wall_width)
		set_vtk_strings(self.db)

	def set_vtk_output_types(self,types):
		self.vtk_output_types=types

	############output methods#############


	def draw_from_zips(self,view_names,cmap_lims=None):
		out.draw_from_zips(self.dir_name,view_names,cmap_lims)

	def make_ffmpeg(self):
		out.make_ffmpeg(self.dir_name)

	def get_range(self,prop_name):
		out.get_range(self.dir_name,prop_name)

	def csv_from_zips(self):
		out.csv_from_zips(self.dir_name)

	def draw_fluxes(self,maxlength):
		out.draw_fluxes(self.dir_name,maxlength)


	def translate_db(self,db):
		species_desc=db.get_property('species_desc')
		cell_trans=self.master_db.get_property('cell_trans')
		wall_trans=self.master_db.get_property('wall_trans')
		edge_trans=self.master_db.get_property('edge_trans')
		for propname,ntype in species_desc.iteritems():
			mprop={}
			prop=db.get_property(propname)
			if ntype=='cell':
				for mcid,cid in cell_trans.iteritems():
					mprop[mcid]=prop[cid]
			if ntype=='wall':
				for mwid,wid in wall_trans.iteritems():
					mprop[mwid]=prop[wid]
			if ntype=='edge':
				for meid,eid in edge_trans.iteritems():
					mprop[meid]=prop[eid]
			self.master_db.set_property(propname,mprop)
			self.master_db.set_description(propname,propname)
			m_species_desc=	self.master_db.get_property('species_desc')
			m_species_desc[propname]=ntype		
			

