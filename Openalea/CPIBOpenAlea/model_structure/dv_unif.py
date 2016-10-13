from numpy import array, around
from model_input.import_utils import set_vtk_strings
from model_utils.db_utilities import get_mesh,get_graph
from model_utils.db_geom import updateVS
from sa_oa.tissueshape import divide_cell
from numpy.random import normal
from model_input.import_tissue import import_tissue
from sa_oa.tissueshape import face_surface_3D, edge_length,cell_volume,face_surface_2D
from sa_oa.celltissue import ConfigFormat,ConfigItem



class Dv_unif(object):
	def __init__(self,db):
		self.db=db
		self.mesh=get_mesh(self.db)
		self.position=self.db.get_property('position')
		db0=import_tissue('full_root_2dt.zip')
		self.mesh0=get_mesh(db0)
		self.pos0=db0.get_property('position')
		

	def get_zrange(self,cid):
		cellpids=[]
		for wid in self.mesh.borders(3,cid):
			for lid in self.mesh.borders(2,wid):
				for pid in self.mesh.borders(1,lid):
					cellpids.append(pid)
		cellpids=list(set(cellpids))
			
		return [min([self.position[pid][2] for pid in cellpids]),max([self.position[pid][2] for pid in cellpids])],cellpids
		
		
	def step(self):
		import time
		print 'dv'
		vtime=time.time()
		div_zone=450.0
		div_thresh=15.0
		cells_to_divide={}
		divided_props=self.db.get_property('divided_props')
		V=self.db.get_property('V')
		S=self.db.get_property('S')
		scale = self.db.get_property('scale_px_to_micron')
		vindex=self.db.get_property('vindex')
		orig_cid=self.db.get_property('orig_cid')
		orig_pid=self.db.get_property('orig_pid')				
		graph=get_graph(self.db)
		wall=self.db.get_property('wall')

		# get list of cells to divide
		vins_to_divide=[]
		vzeds=[]

		maxv=max([vin for vin in vindex.itervalues()])
		zeds=sorted(list(set([pos[2] for pos in self.position.itervalues()])))

		for vin in range(1,maxv):
			cell_length=zeds[vin+1]-zeds[vin]
			mid_cell=zeds[vin]+cell_length/2.0
			if mid_cell<=div_zone and cell_length>=div_thresh:
				vins_to_divide.append(vin)
				vzeds.append([zeds[vin],zeds[vin+1]])
		vins_to_divide.reverse()
		vzeds.reverse()
		

		pos_map={ocid:{} for ocid in orig_cid.itervalues()}
		for cid in self.mesh.wisps(3):
			pos_map[orig_cid[cid]][vindex[cid]]=cid
		pos_map2={vin:{} for vin in range(maxv+1)}

		for cid in self.mesh.wisps(3):
			pos_map2[vindex[cid]][orig_cid[cid]]=cid		

		# divide cells on list
		for i,vin in enumerate(vins_to_divide):
			zrange=vzeds[i]
			newz=zrange[0]+0.5*(zrange[1]-zrange[0])


			the_cells=[cid for cid in pos_map2[vin].itervalues()]
			up_cells=[cid for cid in pos_map2[vin+1].itervalues()]

			all_up_walls=[]
			all_up_lines=[]
			for cid in up_cells:
				for wid in self.mesh.borders(3,cid):
					if wid not in all_up_walls:
						all_up_walls.append(wid)
						for lid in self.mesh.borders(2,wid):
							if lid not in all_up_lines:
								all_up_lines.append(lid)


			all_the_walls=[]
			for cid in the_cells:
				for wid in self.mesh.borders(3,cid):
					if wid not in all_the_walls:
						all_the_walls.append(wid)

			up_walls=[wid for wid in all_up_walls if wid in all_the_walls]

			unlinked={}

			for cid in the_cells:
				for wid in self.mesh.borders(3,cid):
					if wid in up_walls:
						self.mesh.unlink(3,cid,wid)
						unlinked[wid]=[cid,3]
						for lid in self.mesh.borders(2,wid):
							if lid not in unlinked:	
								xwid=[xw for xw in self.mesh.regions(1,lid) if xw not in all_up_walls][0]
								self.mesh.unlink(2,xwid,lid)
								unlinked[lid]=[xwid,2]
							for pid in self.mesh.borders(1,lid):
								if pid not in unlinked:
									xlid=[xl for xl in self.mesh.regions(0,pid)\
											 if xl not in all_up_lines][0]
									self.mesh.unlink(1,xlid,pid)
									unlinked[pid]=[xlid,1]
			trans={}

			for xid,info in unlinked.iteritems():
				nid = self.mesh.add_wisp(info[1])				
				trans[xid]=[nid]
				if info[1]==3:
					for pname in divided_props['cell']:
						prop = self.db.get_property(pname)
						prop[nid]=prop[info[0]]
					V[nid]=0.5*V[info[0]]
					V[info[0]]=0.5*V[info[0]]


				nid = self.mesh.add_wisp(info[1]-1)				
				trans[xid].append(nid)
				if info[1]-1==0:
					opos=self.position[xid]
					self.position[nid]=[opos[0],opos[1],newz]
					orig_pid[nid]=orig_pid[xid]


#trans
#cells and down walls based on up walls - link cells to up/down walls
#sidewalls and down lines based on up lines - link sidewalls to up/down lines
#sidelines and down points based on up points - link sidelines to up/down points

#link cells to sidewalls
#link downwalls to downlines

#link sidewalls to sidelines
#link downlines to downpoints


#link downcells to downwalls
#link downsidewalls to downlines
#link downsidelines to downpoints		

			for xid,info in unlinked.iteritems():
				self.mesh.link(info[1],trans[xid][0],xid)
				self.mesh.link(info[1],trans[xid][0],trans[xid][1])
				#info[0] is 3:downcell for upwall, 2:downsidewalls for upline, 1:downsidelines for uppoints
				#xid is 3:upwall, 2:upline, 1:uppoint
				#info[1] is degree
				#trans[xid] is 3:[cell,downwall],2:[sidewall,downline],1:[sideline,downpoint]
				self.mesh.link(info[1],info[0],trans[xid][1])
				# new wall props
				if info[1]==3:
					for pname in divided_props['wall']:
						prop =self.db.get_property(pname)
						prop[trans[xid][1]]=0.0
				if info[1]==2:
					for pname in divided_props['wall']:
						prop =self.db.get_property(pname)
						prop[trans[xid][0]]=prop[info[0]]				

				if info[1]>1:
					for xnid in self.mesh.borders(info[1]-1,xid):
						#3:uplines for upcell,2:uppoints for uplines
						self.mesh.link(info[1],trans[xid][0],trans[xnid][0])
						self.mesh.link(info[1]-1,trans[xid][1],trans[xnid][1])					
							
			# sort out graph
			rem_eids=[]
			for cid in the_cells:
				for eid in graph.out_edges(cid):
					if graph.target(eid) in up_cells:
						rem_eids.append(eid)
				for eid in graph.in_edges(cid):
					if graph.source(eid) in up_cells:
						rem_eids.append(eid)
			
			for eid in rem_eids:
				graph.remove_edge(eid)
				del wall[eid]
				for pname in divided_props['edge']:
					prop =self.db.get_property(pname)
					del prop[eid]
			done_wids=[]
			for xid,info in unlinked.iteritems():
				if info[1]==3:
					cid = trans[xid][0]
					for wid in self.mesh.borders(3,cid):
						S[wid] = scale[0]*scale[0]*face_surface_3D(self.mesh, self.position, wid)						
						if self.mesh.nb_regions(2,wid) == 2 and wid not in done_wids:
							done_wids.append(wid)							
							cid1,cid2 = self.mesh.regions(2,wid)
							eid1 = graph.add_edge(cid1,cid2)
							wall[eid1] = wid
							eid2 = graph.add_edge(cid2,cid1)
							wall[eid2] = wid
							for pname in divided_props['edge']:
								prop =self.db.get_property(pname)
								prop[eid1]=0.0
								prop[eid2]=0.0
					

			# add new vindex of cells
			for cid in self.mesh.wisps(3):
				try:
					if vindex[cid]>vin:
						vindex[cid]+=1
				except KeyError:
					vindex[cid]=vin+1



		print 'div time',time.time()-vtime
		vtime=time.time()	
		set_vtk_strings(self.db)
		print 'vtk time',time.time()-vtime


