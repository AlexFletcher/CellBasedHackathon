
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
	    for cid in mesh.wisps(3):
		V[cid]=scale[0]*scale[0]*scale[0]*cell_volume(mesh, pos, cid)
    
	    # Lengths of each wall
	    S = refresh_property(db, 'S')
	    for wid in mesh.wisps(2):
		S[wid] = scale[0]*scale[0]*face_surface_3D(mesh, pos, wid)

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


def interior_offset(db,width):


	pos=db.get_property('position')
    	
	#pos3d = dict( (pid,array(list(float(v) for v in vec))) for pid,vec in pos.iteritems() )
	    
	k = array([0, 0, 1])
	    
	mesh=get_mesh(db)
	cell_pts = {}
	cell_inner_pts = {}
	for cid in mesh.wisps(2):
		o_idx = ordered_pids(mesh, cid)
		if area([pos[pid] for pid in o_idx]) < 0:
			o_idx.reverse()
		cell_pts[cid] = o_idx
		N = len(o_idx)
		i_pts = []
		for i in range(N):
			pi = array(pos[o_idx[i]])
		    	ip1 = (i+1)%N
		    	im1 = (i-1)%N
		    	dp = array(pos[o_idx[ip1]]) - pi # pi->p1
		    	dp = dp / norm(dp) # unit pi->p1
			dm = pi - array(pos[o_idx[im1]]) #pl->pi
			dm = dm / norm(dm) # unit pl->pi
			ni = cross((dp + dm),k) # normal to 
			ni = ni / norm(ni)
			normal = cross(dp,k)
			i_pts.append(list(pi - width*ni/dot(normal, ni)))
		cell_inner_pts[cid] = i_pts

	return cell_inner_pts


def face_interior_offset(fid,db,width):


	pos=db.get_property('position')
    	
	#pos3d = dict( (pid,array(list(float(v) for v in vec))) for pid,vec in pos.iteritems() )
	    
	k = array([0, 0, 1])
	    
	mesh=get_mesh(db)


	o_idx = ordered_pids(mesh, fid)
	if area([pos[pid] for pid in o_idx]) < 0:
		o_idx.reverse()
	#p1 = array(pos[o_idx[2]]) - array(pos[o_idx[1]]) 
	#p2 = array(pos[o_idx[0]]) - array(pos[o_idx[1]])  
	#k = cross(p1,p2)
	#k = k/norm(k)
	#print 'k',k
	N = len(o_idx)
	i_pts = []
	for i in range(N):
		pi = array(pos[o_idx[i]])
	    	ip1 = (i+1)%N
	    	im1 = (i-1)%N
	    	dp = array(pos[o_idx[ip1]]) - pi # pi->p1
	    	dp = dp / norm(dp) # unit pi->p1
		dm = pi - array(pos[o_idx[im1]]) #pl->pi
		dm = dm / norm(dm) # unit pl->pi
		ni = cross((dp + dm),k) # normal to 
		ni = ni / norm(ni)
		normal = cross(dp,k)
		i_pts.append(list(pi - width*ni/dot(normal, ni)))
		
	
	
	return init_ordered_pts(i_pts)


	
def horiz_face(fid,db):


  

	pos=db.get_property('position')
	    
	mesh=get_mesh(db)


	idx = ordered_pids(mesh, fid)
	"""
	pt0=array(pos[idx[0]])
	pt1=array(pos[idx[1]])
	ptl=array(pos[idx[-1]])
	nm = cross(pt1-pt0,ptl-pt0)
	if abs(nm[0])==0.0 and abs(nm[1])==0.0:
		return True
	else:
		return False
	"""
	
	zvals = [pos[i][2] for i in idx]

	if zvals.count(min(zvals))==len(zvals):
		return True
	else:
		return False

def init_ordered_pts(i_pts):
	#keep order, but make min x first point. if >1 equal min x point, use min x with min y (z vals should be identical anyway)
	xvals = [i[0] for i in i_pts]
	
	if xvals.count(min(xvals))==1:
		st_pt = xvals.index(min(xvals))

	else:
		yvals = [i[1] for i in i_pts]
		ymind=[]
		yval_xmin=[]
		for i,val in  enumerate(xvals):
			if val==min(xvals):
				ymind.append(i)
				yval_xmin.append(yvals[i])
		if yval_xmin.count(min(yval_xmin))>1:
			print 'coincident points'
		else:
			st_pt = yvals.index(min(yval_xmin))

	ret_pts = i_pts[st_pt:]
	ret_pts.extend(i_pts[:st_pt])		
	
	
	return ret_pts


def face_interior_offset_3d_vert(o_idx):
	"""assumes all faces are horizontal or vertical""" 

	if area(o_idx) < 0:
		o_idx.reverse()
	p1 = array(o_idx[2]) - array(o_idx[1]) 
	p2 = array(o_idx[0]) - array(o_idx[1])  
	k = cross(p1,p2)
	k = k/norm(k)

	N = len(o_idx)
	i_pts = []
	for i in range(N):
		pi = array(o_idx[i])
	    	ip1 = (i+1)%N
	    	im1 = (i-1)%N
	    	dp = array(o_idx[ip1]) - pi # pi->p1
	    	dp = dp / norm(dp) # unit pi->p1
		dm = pi - array(o_idx[im1]) #pl->pi
		dm = dm / norm(dm) # unit pl->pi
		ni = cross((dp + dm),k) # normal to k, pt 1 to -1
		ni = ni / norm(ni)
		normal = cross(dp,k) # normal to k, pt 0 to 1
		i_pts.append(ni/dot(normal, ni))
				
	return k,i_pts,o_idx

def orig_face_interior_offset_3d_vert(o_idx,widths,cell_centre):
	"""assumes all faces are horizontal or vertical""" 

	if area(o_idx) < 0:
		o_idx.reverse()
	p1 = array(o_idx[2]) - array(o_idx[1]) 
	p2 = array(o_idx[0]) - array(o_idx[1])  
	k = cross(p1,p2)
	k = k/norm(k)

	xvals=[i[0] for i in o_idx]
	yvals=[i[1] for i in o_idx]
	zvals=[i[2] for i in o_idx]
	xlen=len(o_idx)
	f_centre=array(sum(xvals)/xlen,sum(yvals)/xlen,sum(zvals)/xlen)
	fc_k=k-face_centre
	fc_cc=array(cell_centre)-face_centre
	k_direction=cmp(norm(fc_k,fc_cc),0)

	N = len(o_idx)
	i_pts = [[] for width in widths]
	for i in range(N):
		pi = array(o_idx[i])
	    	ip1 = (i+1)%N
	    	im1 = (i-1)%N
	    	dp = array(o_idx[ip1]) - pi # pi->p1
	    	dp = dp / norm(dp) # unit pi->p1
		dm = pi - array(o_idx[im1]) #pl->pi
		dm = dm / norm(dm) # unit pl->pi
		ni = cross((dp + dm),k) # normal to 
		ni = ni / norm(ni)
		normal = cross(dp,k)
		for j,width in enumerate(widths):
			i_pts[j].append(list(k*width*k_direction + pi - width*ni/dot(normal, ni)))
			
	
	return i_pts

def get_components(position,pid,n_pids,norm_pid):
	
	dp = array(position[n_pids[1]]) - array(position[pid]) 
    	dp = dp / norm(dp) # unit pi->p1
	dm = array(position[pid]) - array(position[n_pids[0]])
	dm = dm / norm(dm) # unit pl->pi
	k = cross(dp,dm)
	k = k/norm(k)

	ni = cross((dp + dm),k) # normal to k, pt 1 to -1
	ni = ni / norm(ni)
	normal = cross(dp,k) # normal to k, pt 0 to 1
	planar_comp = ni/dot(normal, ni)
	dn = array(position[norm_pid]) - array(position[pid])
	if dp[0]==0.0 and dp[1]==0.0:
		uvec=dm
	else:
		uvec=dp
	xnormal = cross(dn,uvec)
	normal_comp = k#/dot(xnormal,k)


	normal_comp=planar_comp			
	return planar_comp, normal_comp


def get_offset_vec(pos):

	dp = array(pos[1]) - array(pos[0]) 
    	dp = dp / norm(dp) # unit pi->p1
	dm = array(pos[0]) - array(pos[2])
	dm = dm / norm(dm) # unit pl->pi
	k = cross(dp,dm)
	k = k/norm(k)

	ni = cross((dp + dm),k) # normal to k, pt 1 to -1
	ni = ni / norm(ni)
	normal = cross(dp,k) # normal to k, pt 0 to 1
	vec = ni/dot(normal, ni)

	# assumes z direction is vertical
	vec = vec - array([0.0,0.0,1.0])
			
	return vec


