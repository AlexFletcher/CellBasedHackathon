from numpy import zeros
from math import exp
from model_utils.db_utilities import get_mesh, get_graph, get_parameters, set_parameters, get_wall_decomp, get_tissue_maps
from model_utils.celltissue_util import def_property
from numpy import size, sqrt
from model_structure.gn_utils import *

class delta_notch(object):
    """

    Implements AbstractModel
    """
    name = 'delta_notch'
    mobile_species = {}
    
    diluted_props = []
    divided_props = ['delta', 'notch']  

    def get_species(self):

        return [('delta', 'cell'), ('notch', 'cell')]

    def get_steady_species(self):

        return []

    def set_default_parameters(self, db):

	p={}

	p['a']=0.01
	p['k']=2
	p['b']=100
	p['h']=2
	p['v']=1.0

        set_parameters(db, self.name, p)

    def deriv(self, y, t, db, tissue_maps, ssl, wall_decomposition, fixed_idx):

   	self.set_ss(db, tissue_maps, wall_decomposition)
	
	p,var,c,steady=deriv_init(db,ssl,tissue_maps,y,self.get_species(),self.get_steady_species(),self.name)

	ret = deriv_diff(db,y,self.mobile_species,p,tissue_maps,ssl)


	# model equations

	ret[ssl['delta']] += 1/(1+p['b']*pow(var['notch'],p['h'])) - var['delta']

        ret[ssl['notch']] += - p['v']*var['notch']
	

	
	cell_map=tissue_maps['cell']
	graph=get_graph(db)
	ret_notch=ret[ssl['notch']]

	V=db.get_property('V')
	S=db.get_property('S')
	wall=db.get_property('wall')
	cell_type=db.get_property('cell_type')
	mesh=get_mesh(db)
	delta=db.get_property('delta')
	V=db.get_property('V')
	for cid in mesh.wisps(2):
		cidx = cell_map[cid]
		d_avs = [delta[nid] for nid in mesh.border_neighbors(2,cid)]
		#Vtots = [V[nid] for nid in mesh.border_neighbors(2,cid)]
		d_av = sum(d_avs)/len(d_avs)#/sum(Vtots)
		ret_notch[cidx]+=p['v']*(pow(d_av,p['k'])/(p['a']+pow(d_av,p['k'])))
		#if cell_type[cid]==3:
		#	ret_notch[cidx]+=p['nu']
			


        for i in fixed_idx:
            ret[i] = 0.0
        return ret



    def set_ss(self, db, tissue_maps, fixed_steady):

	pass
	

