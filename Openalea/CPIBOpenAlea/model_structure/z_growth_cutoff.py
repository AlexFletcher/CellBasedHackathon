import time
from math import sqrt
from model_utils.db_geom import updateVS
from model_utils.db_utilities import get_mesh,get_graph,set_parameters,get_parameters
from sa_oa.tissueshape import centroid
class z_growth_cutoff(object):

	name = 'z_growth_cutoff'
	
	def __init__(self,db):

		# set basic information
		self.db=db
		self.mesh=get_mesh(self.db)
		self.position=self.db.get_property('position')
		
		# default parameters
		self.p={}
		self.p['maxL']=100.0
		self.p['rate_dz']=0.1
		self.p['rate_ez']=0.5
		self.p['dz']=450.0
		self.p['end_gz']=200.0

		set_parameters(self.db, self.name, self.p)

		intime=time.time()

	# return minimum and maximum z position for a given cell
	def get_zrange(self,cid):
		cellpids=[]
		for wid in self.mesh.borders(3,cid):
			for lid in self.mesh.borders(2,wid):
				for pid in self.mesh.borders(1,lid):
					cellpids.append(pid)
			
		return [min([self.position[pid][2] for pid in cellpids]),max([self.position[pid][2] for pid in cellpids])]


	# function to remove a cell from tissue db and tidy up walls, edges, properties etc.
	def remove_cell(self,cid):
		divided_props=self.db.get_property('divided_props')
		vindex=self.db.get_property('vindex')
		V=self.db.get_property('V')	
		wall_props=['S']
		point_props=['position','orig_pid']
		edge_props=[]
		mesh=self.mesh
		pos=self.position
		wids=list(mesh.borders(3,cid))
		mesh.remove_wisp(3,cid)
		for pname in divided_props['cell']:
			prop=self.db.get_property(pname)
			del prop[cid]
		del V[cid]
		del vindex[cid]

		horiz_walls=self.db.get_property('horiz_walls')
		for wid in wids:
			if mesh.nb_regions(2,wid)==0:
				lids=list(mesh.borders(2,wid))
				mesh.remove_wisp(2,wid)
				if wid in horiz_walls:
					horiz_walls.remove(wid)
				for pname in divided_props['wall']:
					prop=self.db.get_property(pname)
					del prop[wid]			
				for lid in lids:
					if mesh.nb_regions(1,lid)==0:
						pids=list(mesh.borders(1,lid))
						mesh.remove_wisp(1,lid)
						for pid in pids:
							if mesh.nb_regions(0,pid)==0:
								mesh.remove_wisp(0,pid)
								for pname in point_props:
									prop=self.db.get_property(pname)
									del prop[pid]

		graph=get_graph(self.db)
		etd=[]
		wall=self.db.get_property('wall')
		for eid,wid in wall.iteritems():
			try:
				test=graph.source(eid)
			except:
				etd.append(eid)
		for eid in etd:
			del wall[eid]
		for pname in divided_props['edge']:
			prop=self.db.get_property(pname)
			for eid in etd:
				del prop[eid]
		
		
	def step(self):

		print 'growth'

		# get basic info from db. every cell effectively has a 2D coordinate (orig_cell,vindex)
		# i.e. (cell file, position in file)
		vindex=self.db.get_property('vindex')
		orig_cid=self.db.get_property('orig_cid')
		p = get_parameters(self.db, self.name)
		dt=get_parameters(self.db, 'timestep')


		# DEFINE MAPS ETC. FOR NAVIGATION OF TISSUE

		# dict of original cell id:{vindex : current cell id} i.e. the files of cells
		pos_map={ocid:{} for ocid in orig_cid.itervalues()}
		# max and min z value for every cell
		all_czeds={}
		for cid in self.mesh.wisps(3):
			pos_map[orig_cid[cid]][vindex[cid]]=cid
			all_czeds[cid]=self.get_zrange(cid)
		
		# current vindex position in each cell file (initially zero in every file)		
		cvin={ocid:0 for ocid in pos_map.iterkeys()}
		# current z range for each cell file
		czeds={ocid:self.get_zrange(pos_map[ocid][cvin[ocid]]) for ocid in pos_map.iterkeys()}
		# 
		maxvins={ocid:max(list(pos_map[ocid].iterkeys())) for ocid in pos_map.iterkeys()}

		# all unique z positions ordered
		allzeds=[]		
		for pid,pos in self.position.iteritems():
			allzeds.append(pos[2])			
		allzeds=sorted(list(set(allzeds)))

		# current z position
		cpos=allzeds[0]
		all_segs=[]
		cum_growth=0.0


		# MOVE UP ROOT, CALCULATE GROWTH OF EACH SEGMENT, AND NOTE CUMULATIVE GROWTH FOR EACH Z POSITION
		while cpos<allzeds[-1]:

			#list of 'optimal' growth rates
			Gopt=[]

			# find width of current segment
			ind = allzeds.index(cpos)
			gap = allzeds[ind+1]-allzeds[ind]
			"""
			# find Gopt for every cell in segment
			for cid,zv in czeds.iteritems():

				#find cell length and middle z position
				cell_length=zv[1]-zv[0]
				mid_cell=zv[0]+cell_length/2.0

				# get growth parameter depending on middle position of cell (TODO: could also base on cell type here)
				if mid_cell>p['dz']:
					rate=p['rate_ez']*dt
				else:
					rate=p['rate_dz']*dt

				# rate then depends on growth parameter (scaled as proportion of gap to total cell length)
				# and difference between total cell length and max length
				Gopt.append(rate*(gap/(cell_length))*(p['maxL']-(cell_length)))

			# actual growth of segment is the minimum calculated growth of all the files
			# assumption: each cell is constrained by all its neighbours
			# could just as easily take the maximum, or an average, or etc.etc.....
			if len(Gopt)>0:
				seg_growth=min(Gopt)
			else:
				seg_growth=0.0
			"""
			if allzeds[ind+1]<p['dz']:
				seg_growth=p['rate_dz']*dt*gap
			else:
				seg_growth=p['rate_ez']*dt*gap

			
			# add to previous cumulative growth, append to list of cumulative growths fro each segment
			cum_growth+=seg_growth
			all_segs.append(cum_growth)

			# move to next z position and update all relevant dicts etc.
			cpos+=gap
			end_cids=[]
			for cid,zv in czeds.iteritems():
				if zv[1]==cpos:
					if cvin[cid]<maxvins[cid]:
						cvin[cid]+=1
						czeds[cid]=all_czeds[pos_map[cid][cvin[cid]]]
					else:
						end_cids.append(cid)
			for cid in end_cids:
				del czeds[cid]


		# UPDATE POSITION OF EVERY POINT BASED ON CALCULATED GROWTH
		gr_dict=dict(zip(allzeds[1:],all_segs))
	
		for pid,pos in self.position.iteritems():
			if pos[2]>0.0:
				pos[2]+=round(gr_dict[pos[2]],4)



		# TEST IF MIDDLE OF CELL IS OUTSIDE GROWTH ZONE AND IF SO REMOVE IT
		allwisps3=list(self.mesh.wisps(3))
		count=0
		for xcid in allwisps3:
			zr=all_czeds[xcid]
			if zr[0]+(zr[1]-zr[0])*0.5>p['end_gz']:
				self.remove_cell(xcid)
				count+=1


		# update volumes and surface areas of cells and walls, and cell centres - as with all this, it assumes growth in z direction only
		# plus update props that should dilute with growth e.g. auxin, or anything non-nuclear localised
		

		horiz_walls=self.db.get_property('horiz_walls')
		S=self.db.get_property('S')
		V=self.db.get_property('V')
		centres=self.db.get_property('cell_centres')

		diluted_props=self.db.get_property('diluted_props')

		for wid in self.mesh.wisps(2):
			if wid not in horiz_walls:
				cellpids=[]
				for lid in self.mesh.borders(2,wid):
					for pid in self.mesh.borders(1,lid):
						cellpids.append(pid)
				sx=max([self.position[pid][0] for pid in cellpids])-min([self.position[pid][0] for pid in cellpids])
				sy=max([self.position[pid][1] for pid in cellpids])-min([self.position[pid][1] for pid in cellpids])
				sz=max([self.position[pid][2] for pid in cellpids])-min([self.position[pid][2] for pid in cellpids])
				Snew=sz*sqrt(sx*sx+sy*sy)
				for pname,ptype in diluted_props.iteritems():
					if ptype=='wall' or ptype=='WALL':
						prop=self.db.get_property(pname)
						prop[wid]=(prop[wid]*S[wid])/Snew
				S[wid]=Snew

		for cid in self.mesh.wisps(3):
			new_czed=self.get_zrange(cid)
			Vnew=(V[cid]*(new_czed[1]-new_czed[0]))/(all_czeds[cid][1]-all_czeds[cid][0])
			for pname,ptype in diluted_props.iteritems():
				if ptype=='cell' or ptype=='CELL':
					prop=self.db.get_property(pname)
					prop[cid]=(prop[cid]*V[cid])/Vnew
			V[cid]=Vnew
			centres[cid] = centroid (self.mesh, self.position, self.mesh.degree(), cid)
			

