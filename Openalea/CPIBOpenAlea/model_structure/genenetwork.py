from sa_oa.tissueshape import edge_length, face_surface_2D

from model_utils.celltissue_util import def_property
from model_utils.db_utilities import refresh_property, get_wall_decomp, get_tissue_maps, \
    get_parameters,get_graph,get_mesh
from model_utils.db_geom import updateVS

from odesparse import odeints
from scipy.sparse import lil_matrix
from numpy import zeros, concatenate
from scipy.integrate import odeint
class AbstractModel(object):
    """
    Abstract object used to document Model interface -
    implement this interface rather than subclass
    Here for documentation purposes
    """
    
    name = 'AbstractModel'
    # Name of this model - used for storing parameters in TissueDB


    def get_species(self):
        """
        Species used in model
        
        :returns: list of (species_name, type) pairs.
                  e.g. [('auxin', 'cell'), ('PIN', 'edge')]
        """
        pass

    def set_default_parameters(self, db):
        """
        Set default parameters in TissueDB
        :param db: Tissue database
        :type db: TissueDB
        """
        pass

        
    def deriv(y, t, db, tissue_maps, species_slices, wall_decomposition, 
              fixed_idx):
        """
        Right hand side for the system of equations.
        :param y: current state vector
        :type y: numpy.ndarray
        :param t: current simulation time 
        :type t: float
        :param db:  tissue database
        :type db: TissueDB
        :param tissue_maps: dictionary mapping property types ('cell', 'wall',
                            'edge') to a dictionary which takes id numbers to
                            positions in slices of the y for properties of 
                            that type.
                            For example, tissue_maps['cell'] gives a dictionary
                            of the form { 101: 0, 123: 1, 121: 2, ...}
                            so the cell with cid 101 maps to the first element
                            of the slice.
                            y[species_slice['auxin']][tissue_maps['cell'][28]]
                            is the auxin concentration in cell with cid 28
        :type tissue_maps: dict
        :param species_slices: dictionary mapping property names to slices of
                               y. e.g. auxin = y[species_slice['auxin']]
        :type species_slices: dict
        :param wall_decomposition: dictionary which takes a wall id to the
                                   ids of the two corresponding edges.
              
        :param fixed_idx: indices of the elements of y[] which are to be kept
                          constant
        :type fixed_idx: list
        :returns: - numpy.array which is the right hand side of
                        the system of ODEs              
        """
        pass

class CombinedModel(object):
    """
    Additive combination of GRN models
    """

    name = 'CombinedModel'

    def __init__(self, models):
        """
        :param models: List of models to combine (each implementing
                       the AbstractModel interface)
        """
        self.models=models


    def set_default_parameters(self, db):
        for model in self.models:
            model.set_default_parameters(db)

    def get_transporters(self):

        return set.union(*[set(m.transporters) for m in self.models])

    def get_diluted_props(self):

        return set.union(*[set(m.diluted_props) for m in self.models])

    def get_mobile_species(self):

        return set.union(*[set(m.mobile_species) for m in self.models])

    def get_species(self):
        """
        :returns: Set of names of species in the combined model
        """
        return set.union(*[set(m.get_species()) for m in self.models])

    def get_steady_species(self):

        return set.union(*[set(m.get_steady_species()) for m in self.models])

    def get_model_names(self):
	names = []
	for model in self.models:
		names.append(model.name)
	return names
    
    def deriv(self, y, t, *args):
        """
        Derivative function - sum of the derivatives for 
        each component model
        :param y: Current state vector
        :type y: numpy.array (or list)
        :param t: Current time
        :type t: float
        :param *args: additional arguments for derivative function
        """
        return sum((m.deriv(y, t, *args) for m in self.models), \
                       zeros(y.shape))


    def set_ss(self, db, tissue_maps, wall_decomposition):
	for model in self.models:
		model.set_ss(db, tissue_maps, wall_decomposition)

		 

    def get_Jpat(self,species_slices,offset,db):
	from scipy.sparse import lil_matrix
	Jpat=lil_matrix((offset,offset))
	for model in self.models:
		tJpat=model.get_Jpat(species_slices,offset,db)
		Jpat=Jpat+tJpat
	return Jpat

class GeneNetwork(object):
    """
    Class to handle integration of GRN and cell-cell transport
    """
    def __init__(self, db, model):
        """ 
        Initialise the networks 
        :param db: Tissue database
        :type db: TissueDB
        :param model: implementation of AbstractModel


        """
        

        self.db = db
        # Construct a model for each of the model classes
        self.model=model
        self.model_species=model.get_species()
	self.steady_species=model.get_steady_species()
	species_desc=self.db.get_property('species_desc')
	divided_props=self.db.get_property('divided_props')
	diluted_props=self.db.get_property('diluted_props')


        # Initialise all of the species needed which are not already
        # contained in the TissueDB
        for species_name, species_type  in self.model_species:
		species_desc[species_name]=species_type
		if species_name not in db.properties():
		    self.db,prop = def_property(
		            db, 
		            species_name, 
		            0.0,
		            species_type.upper(), 
		            "config",
		            "")
		if species_name not in divided_props[species_type]:
			divided_props[species_type].append(species_name)

        for species_name, species_type  in self.steady_species:
		species_desc[species_name]=species_type
		if species_name not in db.properties():
		    self.db,prop = def_property(
		            db, 
		            species_name, 
		            0.0,
		            species_type.upper(), 
		            "config",
		            "")
		if species_name not in divided_props[species_type]:
			divided_props[species_type].append(species_name)
	for prop in model.get_diluted_props():
		diluted_props[prop]=species_desc[prop]


        model.set_default_parameters(db)

        set_fixed_properties(db)

	updateVS(db)
	
    def step(self):
        """
        Evolve the gene network for one timestep
        """



	#set fixed properties, may use updated dictionaries from previous step
        set_fixed_properties(self.db)
        
        
        t=self.db.tissue()



        # Recalculate the wall to edge map, and the tissue maps
        wall_decomposition = get_wall_decomp(self.db)
        tissue_maps = get_tissue_maps(self.db)
     
        # For each model species, calculate the slice (indices) of
        # the y[] array which will store property values during 
        # integration
        offset=0 # current offset into the y[] array
        species_slices={} # map from species name to the slice of y[]
        for species_name, species_type in self.model_species:
            # 'size' is how many values are needed to store a property
            # of that type 
            size = len(tissue_maps[species_type])
            species_slices[species_name] = slice(offset, offset+size)
            offset+=size
	tJpat=self.get_Jpat(species_slices,offset,self.db,self.model)



	# Extract mesh
        cfg = self.db.get_config('config')
        mesh = t.relation(cfg.mesh_id)
        # Type of each cell
        cell_type = self.db.get_property('cell_type')
        
        # Call routine to find fixed property values
        fixed_idx,fixed_steady=self.get_fixed_idx(tissue_maps, species_slices)


    	# assign steady state values before integration step
	self.model.set_ss(self.db, tissue_maps, fixed_steady)
	#

        #initialise with zero
        y0=zeros((offset,))

        ### SETUP INITIAL VALUES
        # Loop over all species referred to in these models
        for species_name, species_type in self.model_species:
            # extract property dict from TissueDB
            prop=self.db.get_property(species_name)
	    prop = self.db.get_property(species_name)
            # obtain appropriate section of y0[]
            y0_slice=y0[species_slices[species_name]]
            for tid, idx in tissue_maps[species_type].iteritems() :
                y0_slice[idx]=prop[tid]

 
        ###INTEGRATE OVER TIMESTEP
        dt=get_parameters(self.db, 'timestep')
        #integrated = odeint(self.model.deriv, y0, [0., dt], \
        #                        (self.db, tissue_maps, species_slices, 
        #                         wall_decomposition, fixed_idx))
	

	integrated = odeints(self.model.deriv, y0, [0., dt], (self.db, tissue_maps, species_slices,wall_decomposition, fixed_idx)\
				,rtol=1e-8,lrw=1000000,JPat=tJpat)
	#integrated = odeint(self.model.deriv, y0, [0., dt],  (self.db, tissue_maps, species_slices,wall_decomposition, fixed_idx))



        ###SET SPECIES TO NEW VALUE (AFTER TIMESTEP)
        yend=integrated[-1,:]
        
        for species_name, species_type in self.model_species: 
            prop=self.db.get_property(species_name)
            yend_slice=yend[species_slices[species_name]]
            for tid, idx in tissue_maps[species_type].iteritems() :
                prop[tid]=yend_slice[idx]





    def get_fixed_idx(self, tissue_maps, species_slices):
        """
        Obtain a list of those elements in the ODE state vector y
        which are constant during the simulation
        (From the property 'fixed' in the tissue database.)
        This currently only works for cell-based properties
        
        The dictionary fixed maps a property name to a list
        of (expression, value) rules; for each cell, if the expression
        is satisfied, then the property is fixed
        This is a Python expression, and may depend on the cell
        id number and the values of all the properties of the tissue.

        e.g. fixed = { 'auxin':[ ('cell_type[cid]=='pericycle', 2) ] }
        would fix the auxin concentration in the pericycle.

        :param tissue_maps:
        :param species_slices:
        :returns: list of the offsets into the state vector y[] which
                  are to be held constant.

        """
        fixed = get_parameters(self.db, 'fixed')
        
        cell_map = tissue_maps['cell']
	edge_map = tissue_maps['edge']
	wall_map = tissue_maps['wall']
        fixed_idx = [] # list of fixed indices
	fixed_steady={} # dist of steady props and cids
        # Extract all the properties from the tissue data,
        # so these can be used in the eval statement
        all_properties=dict((name, self.db.get_property(name)) \
                                for name in self.db.properties())
        # Loop over all properties with rules

        for prop_name, rules in fixed.iteritems():	    
            # Ignore properties not in model
            if prop_name in species_slices:		
                s=species_slices[prop_name]
                prop=self.db.get_property(prop_name)
                for (expr, value) in rules:
                    # Parse expression and convert it to python
                    # bytecode
                    code=compile(expr, '', 'eval')
		    if dict(self.model_species)[prop_name]=='cell':
		            # Loop over all cells
			    for cid, i in cell_map.iteritems():
				        # Test the expression for this cid
				        if eval(code, all_properties, { 'cid':cid }):
					    #check cid 
					    if cid in prop.iterkeys():
						    idx=s.start+i
						    fixed_idx.append(idx)
		    if dict(self.model_species)[prop_name]=='wall':
		            # Loop over all walls
			    for cid, i in wall_map.iteritems():
				        # Test the expression for this cid
				        if eval(code, all_properties, { 'cid':cid }):
					    #check cid
					    if cid in prop.iterkeys():
						    idx=s.start+i
						    fixed_idx.append(idx)
		    if dict(self.model_species)[prop_name]=='edge':
		            # Loop over all walls
			    for cid, i in edge_map.iteritems():
				        # Test the expression for this cid
				        if eval(code, all_properties, { 'cid':cid }):
					    #check cid
					    if cid in prop.iterkeys():
						    idx=s.start+i
						    fixed_idx.append(idx)
	    
            elif prop_name in dict(self.steady_species):
                prop=self.db.get_property(prop_name)
		fixed_steady[prop_name]=[]
                for (expr, value) in rules:
                    # Parse expression and convert it to python
                    # bytecode
                    code=compile(expr, '', 'eval')
		    if dict(self.steady_species)[prop_name]=='cell':
		            # Loop over all cells
			    for cid, i in cell_map.iteritems():
				        # Test the expression for this cid
				        if eval(code, all_properties, { 'cid':cid }):
					    #check cid 
					    if cid in prop.iterkeys():
						    fixed_steady[prop_name].append(cid)
		    if dict(self.steady_species)[prop_name]=='wall':
		            # Loop over all walls
			    for cid, i in wall_map.iteritems():
				        # Test the expression for this cid
				        if eval(code, all_properties, { 'cid':cid }):
					    #check cid
					    if cid in prop.iterkeys():
						    fixed_steady[prop_name].append(cid)
		    if dict(self.steady_species)[prop_name]=='edge':
		            # Loop over all walls
			    for cid, i in edge_map.iteritems():
				        # Test the expression for this cid
				        if eval(code, all_properties, { 'cid':cid }):
					    #check cid
					    if cid in prop.iterkeys():
						    fixed_steady[prop_name].append(cid)
	

        return fixed_idx,fixed_steady

    def get_Jpat(self,species_slices,offset,db,model):

	JPat=lil_matrix((offset,offset))
	tissue_maps = get_tissue_maps(db)
	cell_map = tissue_maps['cell']
	wall_map = tissue_maps['wall']
	graph=get_graph(db)
	mesh=get_mesh(db)
	wall_deg=mesh.degree()-1
	if 'auxin_wall' in species_slices:

		asl = species_slices['auxin']
		wsl = species_slices['auxin_wall']

		for wid in mesh.wisps(wall_deg):
			cells=list(mesh.regions(wall_deg,wid))
			for cid in cells:
				JPat[asl.start+cell_map[cid],wsl.start+wall_map[wid]]=1.0
				JPat[wsl.start+wall_map[wid],asl.start+cell_map[cid]]=1.0			

		for cid in cell_map.iterkeys():
			JPat[asl.start+cell_map[cid],asl.start+cell_map[cid]]=1.0

		for pid in mesh.wisps(wall_deg-1):
			walls=list(mesh.regions(wall_deg-1,pid))
			idx=[wsl.start+wall_map[wid] for wid in walls]
			for i in idx:
				JPat[idx,i]=1.0


 
	
	return JPat


    def set_fixed_template(self,value):
	if value==False:
		self.fixed_template=False
	else:
		self.fixed_template=False


def set_fixed_properties(db):
    """
    The dictionary fixed maps a property name to a list
    of (expression, value) rules; for each cell, if the expression
    is satisfied, then the property is set to the value.
    This is a Python expression, and may depend on the cell
    id number and the values of all the properties of the tissue.
    If multiple expressions are satisfied, then all of these
    are applied, so the last in the list is used.

    e.g. fixed = { 'auxin':[ ('cell_type[cid]=='pericycle', 2) ] }
    would set the auxin concentration in the pericycle to 2.

    :param db:  tissue database
    :type db: TissueDB
    """
    # Extract mesh from TissueDB
    t = db.tissue()
    cfg = db.get_config('config')
    mesh = t.relation(cfg.mesh_id)
    graph = t.relation(cfg.graph_id)
    # Extract all the properties from the tissue data,
    # so these can be used in the eval statement
    all_properties=dict((name, db.get_property(name)) \
                            for name in db.properties())
    # Loop over all properties with rules 
    for prop_name, rules in get_parameters(db, 'fixed').iteritems():

        # Check if property currently in tissue database
        if prop_name in db.properties():
            prop = db.get_property(prop_name)
            for (expr, value) in rules:
                # compile expression to python bytecode
                code=compile(expr, '', 'eval')

                # loop over all cells
		try:
			for cid in mesh.wisps(mesh.degree()):
			    # test whether expr is satisfied for this cell
			    if eval(code, all_properties, {'cid':cid}):
				#check cid is already a key of property
				if cid in prop.iterkeys():
			        	prop[cid] = value
		#if rule does not apply to cells pass
		except KeyError:
			pass
		try:
			for cid in mesh.wisps(mesh.degree()-1):
			    # test whether expr is satisfied for this wall
			    if eval(code, all_properties, {'cid':cid}):
				#check cid is already a key of property
				if cid in prop.iterkeys():
			        	prop[cid] = value
		#if rule does not apply to walls pass
		except KeyError:
			pass
		try:
			for cid in graph.edges():
			    # test whether expr is satisfied for this edge
			    if eval(code, all_properties, {'cid':cid}):
				#check cid is already a key of property
				if cid in prop.iterkeys():
			        	prop[cid] = value
		#if rule does not apply to walls pass			
		except KeyError:
			pass


