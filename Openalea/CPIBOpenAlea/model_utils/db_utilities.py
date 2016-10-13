from model_utils.celltissue_util import def_property

def set_parameters(db, name, value):
    """
    Set parameter set with name 'name' from TissueDB
    :param db: The tissue database
    :type db: TissueDB
    :param name: Name of parameter set
    :param value: Parameter values (dictionary of name-value pairs)
    :type value: dict
    :type name: str
    """
    
    if 'parameters' in db.properties():
        parameters = db.get_property('parameters')
    else:
        parameters = {}
        db.set_property('parameters', parameters)
        db.set_description('parameters', "Model parameters")
    parameters[name]=value

    
    if type(value)==dict:
	    for pname,ptypes in value.iteritems():
	    	if type(ptypes)==tuple:
			db,par_dict=def_property(db,pname,ptypes[1]['else'],'CELL',"config","")
			prop=db.get_property(ptypes[0])
			for ctype,val in ptypes[1].iteritems():
				for cid in par_dict.iterkeys():
					if prop[cid]==ctype:
						par_dict[cid]=val

def get_parameters(db, name):
    """
    Extract parameter set with name 'name' from TissueDB
    :param db: The tissue database
    :type db: TissueDB
    :param name: Name of parameter
    :type name: str
    :returns: Parameter set dictionary
    """
    return db.get_property('parameters').get(name)

def get_mesh(db):
    """
    Extract graph from TissueDB
    :param db: The tissue database
    :type db: TissueDB
    :returns: Mesh
    """
    t = db.tissue()
    cfg = db.get_config('config')
    return t.relation(cfg.mesh_id)

def get_graph(db):
    """
    Extract graph from TissueDB
    :param db: The tissue database
    :type db: TissueDB
    :returns: graph
    """
    t = db.tissue()
    cfg = db.get_config('config')
    return t.relation(cfg.graph_id)


def refresh_property(db, prop_name):
    """ 
    Look in the db - if there is a property called prop_name, clear it.
    Otherwise, create an empty property of that name
    We don't want to call db.set_property with a new dictionary, as
    other references to that property (particularly in the viewer)
    will still point to the old dictionary
    :param db: The tissue database
    :type db: TissueDB
    :param prop_name: Name of the property
    :type prop_name: str
    :returns: Empty property
    """
    try:
        prop = db.get_property(prop_name)
        prop.clear() # Empty the dictionary
    except KeyError:
        prop = {}
        db.set_property(prop_name, prop) # Create an empty property of that name
    return prop

def get_wall_decomp(db):
    """ 
    Calculate the 'wall decomposition', which maps each wall in the mesh
    seperating a pair of cells to its corresponding pair of edges in the
    graph
    :param db: Tissue Database
    :type db: TissueDB
    :returns dict: Wall decomposition. Map from wall ids to the
                   corresponding pair of edge ids.
    """
    wall_decomposition = {}
    wall=db.get_property('wall')
    for eid in wall:
        if wall[eid] not in wall_decomposition:
            wall_decomposition[wall[eid]] = [eid]
        else:
            wall_decomposition[wall[eid]].append(eid)
    return wall_decomposition


def get_tissue_maps(db):
    """
    Calculate the 'cell_map' and 'edge_map':
    maps between the id numbers in the graph and their
    position in the list of the cells and edges
    :param db: Tissue Database (modified in place)
    :type db: TissueDB
    :returns: {'cell': cell_map, 'edge': edge_map}
              cell_map is a dictionary mapping cell ids
              to offsets. Similarly for edge_map.
    """
    t = db.tissue()
    cfg = db.get_config('config')
    graph = t.relation(cfg.graph_id) # directed cell to cell graph
    mesh = t.relation(cfg.mesh_id)

    cell_map={}  
    for ind, cid in enumerate(mesh.wisps(mesh.degree())) :
	cell_map[cid] = ind
    edge_map={}
    for ind, eid in enumerate(graph.edges()) :
	edge_map[eid] = ind
    wall_map={}
    for ind, wid in enumerate(mesh.wisps(mesh.degree()-1)) :
	wall_map[wid] = ind

    return {'cell': cell_map, 'edge': edge_map, 'wall': wall_map }

def add_cell_property(db, name, init_value):
    """
    Add a property to the tissue database, with every cell
    having the value init_value
    :param db: Tissue Database
    :type db: TissueDB
    :param name: Name of parameter
    :type name: str
    :param init_value: Default initial value for the parameter
    :type init_value: float, int (any?)
    """
    mesh=get_mesh(db)
    p = dict( (cid, init_value) for cid in mesh.wisps(2) )
    db.set_property(name, p)
    db.set_description(name, '')
    return p


def iter_prop_type(db, prop_type):
    #needs updating for use with 3D
    if prop_type == 'cell':
        mesh=get_mesh(db)
        return mesh.wisps(2)
    if prop_type == 'wall':
        mesh = get_mesh(db)
        return mesh.wisps(1)
    if prop_type == '(cell, wall)':
        mesh = get_mesh(db)
        return cell_walls(mesh)
    if prop_type == 'point':
        mesh = get_mesh(db)
        return mesh.wisps(0)
    if prop_type == 'edge':
        graph = get_graph(db)
        return graph.edges()
    if prop_type == '(cell, tri)':
        mesh = get_mesh(db)
        cell_triangles = db.get_property('cell_triangles')
        return ((cid, t) for cid, l in cell_triangles.iteritems() for t in l)
        


def add_property(db, prop_name, prop_type, default_value):
    """
    Add a property to the tissue database, with every cell
    having the value init_value
    :param db: Tissue Database
    :type db: TissueDB
    :param prop_name: Name of parameter
    :type prop_name: str
    :param prop_type: Type of mesh element property is defined upon
    :type prop_type: str
    :param init_value: Default initial value for the parameter
    :type init_value: float, int (any?)
    """

    p = {}
    for wid in iter_prop_type(db, prop_type):
        p[wid] = default_value
    db.set_property(prop_name, p)
    db.set_description(prop_name, '')
    return p

def cell_walls(mesh):
    #needs updating for use with 3D
    """
    Extract all (cell, wall) pairs in a Topomesh 
    :param mesh: The mesh
    :type mesh: Topomesh
    :returns: generator of (cell, wall) id pairs
    """
    for wid in mesh.wisps(1):
        for cid in mesh.regions(1, wid):
            yield (cid, wid)


def translate_cell_types_and_borders(db):
    """
    converts imported/given cell types and border types to standard values
    to add/change standard values change the dictionarys below
    ** also changes config for saving db again **
    """
    
    cfg = db.get_config('config')
    print 'imported cell types',cfg.cell_types
    standard_cell_types={'epiderm':0,'cortex':1,'endoderm':2,'pericycle':3,'vasculature':4,'columella':7, 'cap':8, 'QC':9,'primordium':5}
    print 'standard_cell_types',standard_cell_types
    ## convert cell_types to 'standard' values
    cell_type=db.get_property('cell_type')
    for cid,ctype in cell_type.iteritems():
	if ctype is not None:
    		cell_type[cid]=standard_cell_types[cfg.cell_types[ctype]]
    ## change config to match
    oldcfg=cfg.cell_types
    cfg.cell_types={}
    for ctype in oldcfg.itervalues():
	cfg.cell_types[standard_cell_types[ctype]]=ctype
    
    #convert border types
    print 'imported border types',cfg.border
    standard_border_types={'apical':5,'basal':6,'primed':7}
    print 'standard border types',standard_border_types

    border = db.get_property('border')
    for cid,btype in border.iteritems():
	if btype is not None:
    		border[cid]=standard_border_types[cfg.border[btype]]
    #amend config
    oldcfg=cfg.border
    cfg.border={}
    for btype in oldcfg.itervalues():
	cfg.border[standard_border_types[btype]]=btype   

def get_outer_walls(db):
    """
    returns dictionary of outer walls of tissue with wall ids as keys and adjoining cell as value
    """
    mesh=get_mesh(db)
    wall=db.get_property('wall')
    outer_walls={}
    for wid in mesh.wisps(mesh.degree()-1):
	if wid not in wall.itervalues():
		out_cells=[]
		for cid in mesh.regions(mesh.degree()-1,wid):
			out_cells.append(cid)
		outer_walls[wid]=out_cells[0]
    return outer_walls

def get_cell_pids(mesh):
	""" returns dictionary of cell id and all associated points for given 3D mesh"""
	cell_pids={}
	for cid in mesh.wisps(3):
		cell_pids[cid]=[]
		for fid in mesh.borders(3,cid):
			for lid in mesh.borders(2,fid):
				for pid in mesh.borders(1,lid):
					if pid not in cell_pids[cid]:
						cell_pids[cid].append(pid)
	return cell_pids

