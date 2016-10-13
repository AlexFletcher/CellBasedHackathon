from numpy import zeros
from model_utils.db_utilities import get_parameters,get_graph,get_mesh,set_parameters,get_wall_decomp
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

def deriv_init(db,ssl,tissue_maps,y,species,steady_species,name):
	p = get_parameters(db, name)
	set_parameters(db, name, p)
	var={}
	for sname,stype in species:
		var[sname]=y[ssl[sname]]


	cell_map = tissue_maps['cell']
	wall_map = tissue_maps['wall']
	edge_map = tissue_maps['edge']
	c={}
	#convert dictionaries of cell type dependent parameters to slices
	for pname,ptypes in p.iteritems():
		if type(ptypes)==tuple:
			adict = db.get_property(pname)
			tdict = dict((cell_map[cid],val) for cid,val in adict.iteritems())
			alist=[]
			for i in range(0,len(cell_map)):
				alist.append(tdict[i])
			c[pname]=alist[slice(0,len(cell_map))]

	steady={}


	for sname,stype in steady_species:
		nprop=db.get_property(sname)
		if stype=='cell':
			steady[sname]=zeros(len(cell_map))
			for cid,val in nprop.iteritems():
				steady[sname][cell_map[cid]]=nprop[cid]			
		if stype=='wall':
			steady[sname]=zeros(len(wall_map))
			for wid in nprop.iterkeys():
				steady[sname][wall_map[wid]]=nprop[wid]		
		if stype=='edge':
			steady[sname]=zeros(len(edge_map))
			for eid in nprop.iterkeys():
				steady[sname][edge_map[eid]]=nprop[eid]
	return p,var,c,steady

def deriv_diff(db,y,mobile_species,p,tissue_maps,ssl):
	ret = zeros((len(y),))
	graph = get_graph(db)      
	S = db.get_property('S')
	V = db.get_property('V')
	cell_map = tissue_maps['cell']

	wall=db.get_property('wall')

	for eid in graph.edges():
		sid=graph.source(eid)
		tid=graph.target(eid)
		sidx = cell_map[sid]
		tidx = cell_map[tid]

		Vs = float(V[sid])
		Vt = float(V[tid])
		Swall = float(S[wall[eid]])
		for sname,pval in mobile_species.iteritems():
			ret_sp=ret[ssl[sname]]
			s_sp = y[ssl[sname]][sidx]
			ret_sp[sidx] += -Swall/Vs * p[pval]*s_sp			 
			ret_sp[tidx] += Swall/Vt * p[pval]*s_sp


	return ret
	

def deriv_pars(db,tissue_maps,name):
	p = get_parameters(db, name)
	set_parameters(db, name, p)

	cell_map = tissue_maps['cell']
	wall_map = tissue_maps['wall']
	edge_map = tissue_maps['edge']
	c={}

	#convert dictionaries of cell type dependent parameters to slices
	for pname,ptypes in p.iteritems():
		if type(ptypes)==tuple:
			adict = db.get_property(pname)
			tdict = dict((cell_map[cid],val) for cid,val in adict.iteritems())
			alist=[]
			for i in range(0,len(cell_map)):
				alist.append(tdict[i])
			c[pname]=alist[slice(0,len(cell_map))]
	return p,c


def calc_ss_walls_cells(db,cellpropname,wallpropname,tissue_maps,InFluxes,OutFluxes,InFluxes_noedge,OutFluxes_noedge,production,degradation,shoot_boundary,tip_boundary,lamb, fixed_steady,wall_diff):

        S = db.get_property('S')
        V = db.get_property('V')
	graph=get_graph(db)
        mesh = get_mesh(db)
	centres = db.get_property('cell_centres')	
	vindex=db.get_property('vindex')
	position=db.get_property('position')

 	cellprop=db.get_property(cellpropname)
	wallprop=db.get_property(wallpropname)

        cell_map = tissue_maps['cell']
        wall_map = tissue_maps['wall']	
	try:
		horiz_walls=db.get_property('horiz_walls')
	except KeyError:
		horiz_walls=[]
	wall_decomposition=get_wall_decomp(db)
	
	c_off = len(list(cell_map.iterkeys()))
	w_off = len(list(wall_map.iterkeys()))


        N = c_off+w_off
        J = lil_matrix((N,N))
        r = zeros((N,))

	walldeg=mesh.degree()-1
	#set fluxes	
        for wid in mesh.wisps(walldeg):
		widx = c_off + wall_map[wid]
		if wid in wall_decomposition:
			eid1 = wall_decomposition[wid][0]
			eid2 = wall_decomposition[wid][1]		
			sid = graph.source(eid1)
			tid = graph.target(eid1)
			sidx = cell_map[sid]
			tidx = cell_map[tid]

		
			J[sidx, sidx] += -OutFluxes[eid1]*S[wid]/V[sid]
			J[sidx, widx] += InFluxes[eid1]*S[wid]/V[sid]
			J[widx, sidx] += OutFluxes[eid1]/lamb
			J[widx, widx] += -InFluxes[eid1]/lamb

			J[tidx, tidx] += -OutFluxes[eid2]*S[wid]/V[tid]
			J[tidx, widx] += InFluxes[eid2]*S[wid]/V[tid]
			J[widx, tidx] += OutFluxes[eid2]/lamb
			J[widx, widx] += -InFluxes[eid2]/lamb

		else:
			sid = list(mesh.regions(walldeg,wid))[0]
			sidx = cell_map[sid]
			J[sidx, sidx] += -OutFluxes_noedge[wid]*S[wid]/V[sid]
			J[sidx, widx] += InFluxes_noedge[wid]*S[wid]/V[sid]
			J[widx, sidx] += OutFluxes_noedge[wid]/lamb
			J[widx, widx] += -InFluxes_noedge[wid]/lamb
	
	

	#production / degradation
        for cid in mesh.wisps(mesh.degree()):
            cidx = cell_map[cid]      
            J[cidx, cidx] -= degradation
            r[cidx] += production[cidx]


	for wid in horiz_walls:
		cids=list(mesh.regions(walldeg,wid))
		if len(cids)==1:
			cid=cids[0]
			zval = position[list(mesh.borders(walldeg-1,list(mesh.borders(walldeg,wid))[0]))[0]][2] # this is just the z position of any point of the wall
			cidx = cell_map[cid]
			widx = c_off + wall_map[wid]


			if zval>centres[cid][2]:
				flux=shoot_boundary[cidx]
			else:
				flux=tip_boundary[cidx]
				
			if flux>0.0:
				r[widx] +=  flux/lamb
			if flux<0.0:
				J[widx, widx] += flux/lamb
	


	# wall-wall diffusion

	if wall_diff>0.0:
		if walldeg==1:
			print 'wall diffusion in 2D mesh TODO'
		for wid in mesh.wisps(walldeg):
			widx = c_off + wall_map[wid]
			for  nwid in mesh.border_neighbors(walldeg,wid):
				nwidx = c_off + wall_map[nwid]
				J[widx, widx] -= wall_diff/S[wid]
				J[nwidx, widx] += wall_diff/S[nwid]	


	# fixed properties	
	for prop_name,fcids in fixed_steady.iteritems():		
		if prop_name==cellpropname:
			for cid in fcids:				
				i=cell_map[cid]
				J[i, :] *=0
				J[i, i] = 1
				r[i] = -cellprop[cid]

	for prop_name,fcids in fixed_steady.iteritems():		
		if prop_name==wallpropname:
			for cid in fcids:				
				i=c_off + wall_map[cid]
				J[i, :] *=0
				J[i, i] = 1
				r[i] = -wallprop[cid]

	#solve
        y = -spsolve(J.tocsr(), r)

	#assign values
	for cid in cellprop.iterkeys():
		cellprop[cid]=y[cell_map[cid]]		

	for wid in wallprop.iterkeys():
		wallprop[wid]=y[c_off+wall_map[wid]]	

