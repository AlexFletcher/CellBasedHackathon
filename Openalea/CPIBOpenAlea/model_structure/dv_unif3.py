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
		divided_props=['orig_cid','cell_type','auxin']
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
		vindex=self.db.get_property('vindex')
		maxv=max([vin for vin in vindex.itervalues()])
		zeds=sorted(list(set([pos[2] for pos in self.position.itervalues()])))
		print 'zeds',zeds
		for vin in range(1,maxv):
			cell_length=zeds[vin+1]-zeds[vin]
			print 'cell',cell_length
			mid_cell=zeds[vin]+cell_length/2.0
			print 'mc',mid_cell
			if mid_cell<=div_zone and cell_length>=div_thresh:
				vins_to_divide.append(vin)
				vzeds.append([zeds[vin],zeds[vin+1]])
		vins_to_divide.reverse()
		vzeds.reverse()
		print vins_to_divide,'vins to divide'
		print vzeds,'vzeds'



		pos_map={ocid:{} for ocid in orig_cid.itervalues()}
		for cid in self.mesh.wisps(3):
			pos_map[orig_cid[cid]][vindex[cid]]=cid
		pos_map2={vin:{} for vin in range(maxv+1)}
		print 'cid',cid
		print pos_map2
		for cid in self.mesh.wisps(3):
			pos_map2[vindex[cid]][orig_cid[cid]]=cid		

		# divide cells on list
		for i,vin in enumerate(vins_to_divide):
			zrange=vzeds[i]
			#the_cells={ocid:cid for ocid,cid in pos_map2[vin].iteritems()}
			#up_cells={ocid:cid for ocid,cid in pos_map2[vin+1].iteritems()}
			#all_up_walls={ocid:list(self.mesh.borders(3,up_cells[ocid])) for ocid in up_cells.iterkeys()}
			#all_the_walls={ocid:list(self.mesh.borders(3,the_cells[ocid])) for ocid in the_cells.iterkeys()}
			#up_walls={ocid:[wid for wid in all_up_walls[ocid] if wid in all_the_walls[ocid]] for ocid in up_cells.iterkeys()}
			#up_lines={ocid:[lid for lid in self.mesh.borders(2,up_walls[ocid][0])] for ocid in up_cells.iterkeys()}


			the_cells=[cid for cid in pos_map2[vin].itervalues()]
			up_cells=[cid for cid in pos_map2[vin+1].itervalues()]

			all_up_walls=[]
			for cid in up_cells:
				for wid in self.mesh.borders(3,cid):
					if wid not in all_up_walls:
						all_up_walls.append(wid)
			all_the_walls=[]
			for cid in the_cells:
				for wid in self.mesh.borders(3,cid):
					if wid not in all_the_walls:
						all_the_walls.append(wid)


			up_walls=[wid for wid in all_up_walls if wid in all_the_walls]

			up_lids=[]
			for wid in up_walls:
				for lid in self.mesh.borders(2,wid):
					if lid not in up_lids:
						up_lids.append(lid)
			up_points=[]
			for lid in up_lids:
				for pid in self.mesh.borders(1,lid):
					if pid not in up_points:
						up_points.append(pid)
			unlinked={}

			for cid in the_cells:
				for wid in self.mesh.borders(3,cid):
					if wid in up_walls:
						self.mesh.unlink(3,cid,wid)
						unlinked[wid]=cid
						for lid in self.mesh.borders(2,wid):
							if lid not in unlinked:					
								self.mesh.unlink(2,wid,lid)
								unlinked[lid]=wid
							for pid in self.mesh.borders(1,lid):
								if pid not in unlinked:
									self.mesh.unlink(1,lid,pid)
									unlinked[pid]=lid

				
							
			
			"""
			down_cells={ocid:cid for ocid,cid in pos_map2[vin-1].iteritems()}



			all_down_walls={ocid:list(self.mesh.borders(3,down_cells[ocid])) for ocid in down_cells.iterkeys()}
			

			down_walls={ocid:[wid for wid in all_down_walls[ocid] if wid in all_the_walls[ocid]]\
					for ocid in down_cells.iterkeys()}

			
			for ocid,cid in the_cells.iteritems():
				self.mesh.unlink(3,cid,up_walls[ocid][0])

			the_side_walls={}
			for cid in the_cells.itervalues():
				the_side_walls[cid]=[]
				for wid in self.mesh.borders(3,cid):
					if wid not in [wids[0] for wids in up_walls.itervalues()] and wid not in [wids[0] for wids in down_walls.itervalues()]:
						the_side_walls[cid].append(wid)
			up_side_walls={}
			for cid in the_cells.itervalues():
				the_side_walls[cid]=[]
				for wid in self.mesh.borders(3,cid):
					if wid not in [wids[0] for wids in up_walls.itervalues()] and wid not in [wids[0] for wids in down_walls.itervalues()]:
						the_side_walls[cid].append(wid)
			print 'the_side_walls',the_side_walls
			"""
			# add new vindex of cells
			

		#haven't linked to lower cells

		"""	
		for cid,zv in cells_to_divide.iteritems():
			#pick division plane and divide cell
			cell_length=zv[1]-zv[0]
			mid_cell=zv[0]+cell_length/2
			zed = round(normal(mid_cell,1.0),4)

			lineage = divide_cell (self.mesh, self.position, cid, array([0,0,zed]), array([0,0,1]))

			#print 'lineage',lineage
			
			for xwid in lineage[2].iterkeys():
				if xwid is not None:
					del S[xwid]

			# get info from division lineage
			ocid=list(lineage[3].iterkeys())[0]

			ncell1=list(lineage[3].itervalues())[0][0]
			ncell2=list(lineage[3].itervalues())[0][1]
			zv1,cellpids1=self.get_zrange(ncell1)
			zv2,cellpids2=self.get_zrange(ncell2)


			#round position values and update orig_pid
			ptrans={}
			for pid in cellpids1: 
				pos=self.position[pid]
				self.position[pid]=[round(pos[0],4),round(pos[1],4),round(pos[2],4)]
				if pid not in lineage[0][None]:
					ptrans[orig_pid[pid]] = (self.position[pid][0],self.position[pid][1])

			for pid in cellpids2:
				pos=self.position[pid]
				self.position[pid]=[round(pos[0],4),round(pos[1],4),round(pos[2],4)]
		
			for pid in lineage[0][None]:
				px=self.position[pid][0]
				py=self.position[pid][1]
				for xpid,ref in ptrans.iteritems():
					if px==ref[0] and py==ref[1]:
						orig_pid[pid]=xpid

			#update vindex
			zv1,cellpids1=self.get_zrange(ncell1)
			zv2,cellpids2=self.get_zrange(ncell2)
			ovin=vindex[ocid]
			for xcid,vin in vindex.iteritems():
				if vin>ovin and orig_cid[xcid]==orig_cid[ocid]:
					vindex[xcid]+=1
			zs=sorted([zv1[0],zv1[1],zv2[0],zv2[1]])
			if zs[0] in zv1:
				vindex[ncell1]=vindex[ocid]
				vindex[ncell2]=vindex[ocid]+1
			else:
				vindex[ncell2]=vindex[ocid]
				vindex[ncell1]=vindex[ocid]+1

			#update divided props
			for pname in divided_props:
				prop =self.db.get_property(pname)
				prop[ncell1]=prop[ocid]
				prop[ncell2]=prop[ocid]
				del prop[ocid]
			del vindex[cid]
			del V[cid]

			#update graph, S
			for wid in self.mesh.borders(3,ncell1):
				S[wid] = scale[0]*scale[0]*face_surface_3D(self.mesh, self.position, wid)
				if self.mesh.nb_regions(2,wid) == 2 :
					cid1,cid2 = self.mesh.regions(2,wid)
					eid1 = graph.add_edge(cid1,cid2)
					wall[eid1] = wid
					eid2 = graph.add_edge(cid2,cid1)
					wall[eid2] = wid
			for wid in self.mesh.borders(3,ncell2):
				S[wid] = scale[0]*scale[0]*face_surface_3D(self.mesh, self.position, wid)
				if self.mesh.nb_regions(2,wid) == 2 :
					cid1,cid2 = self.mesh.regions(2,wid)
					if cid1!=ncell1 and cid2!=ncell1:
						eid1 = graph.add_edge(cid1,cid2)
						wall[eid1] = wid
						eid2 = graph.add_edge(cid2,cid1)
						wall[eid2] = wid

			# update V
			V[ncell1]=scale[0]*scale[0]*scale[0]*cell_volume(self.mesh, self.position, ncell1)
			V[ncell2]=scale[0]*scale[0]*scale[0]*cell_volume(self.mesh, self.position, ncell2)

		#tidy up
		wtd=[]
		for eid,wid in wall.iteritems():
			try:
				test=graph.source(eid)
			except:
				wtd.append(eid)
		for eid in wtd:
			del wall[eid]
		"""	
		
		set_vtk_strings(self.db)
		print 'div time',time.time()-vtime


