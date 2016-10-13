
from scipy.linalg import norm
from numpy import array,cross,dot
from geom import area, centroid
from sa_oa.tissueshape import face_surface_3D, edge_length,cell_volume,face_surface_2D
from db_utilities import get_mesh, get_graph, refresh_property
from sa_oa.container import ordered_pids

def wall_angle(db, wid):
    """
    Determine whether a wall is nearer to vertical than it is to
    being horizontal (i.e. making an angle with the horizontal of
    between pi/4 and 3pi/4, or between -3pi/4 and -pi/4)
    :param db: Tissue database
    :type db: TissueDB
    :param wid: Wall Id
    :type wid: int
    """
    pos=db.get_property('position')
    mesh=get_mesh(db)
    pts=[pos[pid] for pid in mesh.borders(1,wid)]
    return abs(pts[1][1]-pts[0][1])>abs(pts[1][0]-pts[0][0])

def ordered_wids_pids(mesh, pos, cid):
    wid_set = set(mesh.borders(2,cid))
    wid = wid_set.pop()
    o_wids = [wid]
    end = list(mesh.borders(1, wid))[-1]
    o_pids= [end]
    while wid_set:
        for wid in wid_set:
            if end in mesh.borders(1, wid):
                wid_set.discard(wid)
                o_wids.append(wid)
                end = (set(mesh.borders(1, wid)) - set([end])).pop()
                o_pids.append(end)
                break
    if area([pos[pid] for pid in o_pids])<0:
        o_wids.reverse()
        o_pids=o_pids[-2::-1]+[o_pids[-1]]
    return o_wids, o_pids



def updateVS(db):
    """ 
    Calculate cell volumes and surface area between cells
    (updating these in the TissueDB properties V and S). 
    :param db: Tissue Database (modified in place)
    :type db: TissueDB
    """
#    graph = get_graph(db)
    mesh = get_mesh(db)
    pos = db.get_property('position')
    wall = db.get_property('wall')
    scale = db.get_property('scale_px_to_micron')

    if mesh.degree()==3:
	    # Area of each cell
	    V = refresh_property(db, 'V')
	    scale3=scale[0]*scale[0]*scale[0]
	    for cid in mesh.wisps(3):
		V[cid]=scale3*cell_volume(mesh, pos, cid)
    
	    # Lengths of each wall
	    scale2=scale[0]*scale[0]
	    S = refresh_property(db, 'S')
	    for wid in mesh.wisps(2):
		S[wid] = scale2*face_surface_3D(mesh, pos, wid)

    if mesh.degree()==2:
	    # Area of each cell
	    V = refresh_property(db, 'V')
	    for cid in mesh.wisps(2):
		V[cid]=scale[0]*scale[0]*face_surface_2D(mesh, pos, cid)

	    
	    # Lengths of each wall
	    S = refresh_property(db, 'S')
	    for wid in mesh.wisps(1):
		S[wid] = scale[0]*edge_length(mesh, pos, wid)
	    

def getSbyCell(db):
    """
    gives sum of lengths of walls surrounding each cell
    """	
    mesh = get_mesh(db)

    # Lengths of each wall
    from model_utils.celltissue_util import get_property
    S = get_property(db, 'S')

    Scell={}
    for cid in mesh.wisps(2):
	Slist=[]
	for wids in mesh.borders(2,cid):
		Slist.append(S[1][wids])
		Scell[cid]=sum(Slist)
    db.set_property('Scell',Scell)



def get_offset_vec2(pvals):
	#points: pt, pt+1, pt-1, centroid

	dp = array(pvals[1]) - array(pvals[0]) 
    	dp = dp / norm(dp) # unit pi->p1
	dm = array(pvals[0]) - array(pvals[2])
	dm = dm / norm(dm) # unit pl->pi

	to_centre=array(pvals[3])-array(pvals[0])
	k = cross(to_centre,dm)
	k = k/norm(k)

	ni = cross((dp + dm),k) # normal to k, pt 1 to -1

	ni = ni / norm(ni)
	normal = cross(dp,k) # normal to k, pt 0 to 1
	dot_p = dot(normal, ni)

	vec = ni/dot_p
	# assumes z direction is vertical
	#vec = vec - array([0.0,0.0,1.0])
	
	return vec
