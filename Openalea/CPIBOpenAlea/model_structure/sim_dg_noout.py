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
from model_input.import_utils import set_vtk_strings

#from model_structure.parameters import writeParameters

class Sim_dg_noout(object):
	def __init__(self,*args,**kwargs):

		self.step_start_time=time.time()
		self.steptimes=[]
		
		self.sch = Scheduler()
		self.division=False
		self.growth=False

		self.initial_values={}
		self.timestep=1
		self.fixed_dict={}
		self.param_values={}
		self.vtk_output_types=['cell','wall','edge']
		
		self.tissuename=args[0].replace('/','_').split('.')[0]
		self.db=import_tissue(args[0])


		set_parameters(self.db, 'initial_values', self.initial_values)
		set_parameters(self.db, 'timestep', self.timestep)
		set_parameters(self.db, 'fixed', self.fixed_dict)


		m_init_funcs={'genenetwork':self.init_genenetwork,'division':self.init_division,'growth':self.init_growth}
		for mtype,model in kwargs.iteritems():
			m_init_funcs[mtype](model)
			

		
		names=self.combo_model.get_model_names()
		cnames=''
		for name in names:
			cnames='%s_%s' % (cnames,name)
		self.modelname=cnames
		self.simname=int(time.time())
		self.dir_name= '%s/%s/t%10d' % (self.tissuename,self.modelname,self.simname)
		#os.system('mkdir %s' % (self.tissuename))
		#os.system('mkdir %s/%s' % (self.tissuename,self.modelname))
		#os.system('mkdir %s' % self.dir_name)
		#os.system('cp %s %s/%s' % (sys.argv[0],self.dir_name,sys.argv[0]))
		self.db.set_property('iteration', 0)
		#self.output()
		#self.sch.register(Task(self.output, 1, 1, "output"))


		def update_iteration():
			t=self.db.get_property('time')
			t+=get_parameters(self.db, 'timestep')
			self.db.set_property('time', t)

			iteration=self.db.get_property('iteration')
			iteration+=1
			self.db.set_property('iteration', iteration)
			self.steptimes.append(time.time()-self.step_start_time)
			self.step_start_time=time.time()
		self.sch.register(Task(update_iteration, 1, 10, "add iteration"))



		def update_vtk():
			set_vtk_strings(self.db)
		if self.division or self.growth:
			self.sch.register(Task(update_vtk, 1, 3, "set vtk strings"))


	#######run simulations#############
	
	def run_simulation(self,first_step,last_step):
			
		self.db.set_property('iteration', first_step)
		
		for i in range(first_step,last_step):
			#print 'iteration',i+1,'of',last_step
			next ( self.sch . run ())


		#final output

		#out.writeParameters(self.db,self.dir_name)
		#print 'writing vtu time series'
		#vtu_series(self.dir_name,self.vtk_output_types)
		#print 'done'
		#print 'NOT writing .csv file'	
		#self.csv_from_zips()

		#print 'simulation time:',sum(self.steptimes)




	#########change settings ###############

	def set_sim_name(self,name):
		old_name = self.dir_name
		self.dir_name= '%s/%s/%s' % (self.tissuename,self.modelname,name)
		os.system('mv %s %s' % (old_name,self.dir_name))

	def get_db(self):
		return self.db

	def set_initial_values(self,ivs):
		for pname,val in ivs.iteritems():
			self.initial_values[pname]=val
			prop=self.db.get_property(pname)
			for cid in prop.iterkeys():
				prop[cid]=val
		set_parameters(self.db, 'initial_values', self.initial_values)
		self.output()

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
		if self.growth:
			p = get_parameters(self.db, self.gr.name)
			for par,value in self.param_values.iteritems():
				if par in p.iterkeys():
					p[par]=value
					print 'par',par,'val',value,'overwritten'
			set_parameters(self.db, self.gr.name, p)
		if self.division:
			p = get_parameters(self.db, self.dv.name)
			for par,value in self.param_values.iteritems():
				if par in p.iterkeys():
					p[par]=value
					print 'par',par,'val',value,'overwritten'
			set_parameters(self.db, self.dv.name, p)			

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

	def init_genenetwork(self,model_list):
		if type(model_list)==list:
			mlist=[model() for model in model_list]
		else:
			mlist=[model_list()]
		self.combo_model=CombinedModel(mlist)

		self.gene_network=GeneNetwork(self.db, self.combo_model)
		self.sch.register(Task(self.gene_network.step, 1, 6, "genenetwork") )

	def init_growth(self,model):
		self.gr=model(self.db)
		self.sch.register(Task(self.gr.step, 1, 8, "growth") )
		self.growth=True

	def init_division(self,model):
		self.dv=model(self.db)
		self.sch.register(Task(self.dv.step, 1, 7, "division") )		
		self.division=True

	def output(self):
		i=self.db.get_property('iteration')
		self.db.write('%s/iter_%05d.zip' % (self.dir_name,i))
		db_to_vtu(self.db,self.dir_name,'iter_%05d.vtu' % i,self.vtk_output_types)

