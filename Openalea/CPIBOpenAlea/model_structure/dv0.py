from numpy import array, around
from model_input.import_utils import set_vtk_strings
from model_utils.db_utilities import get_mesh,get_graph
from model_utils.db_geom import updateVS
from sa_oa.tissueshape import divide_cell
from numpy.random import normal
from sa_oa.tissueshape import face_surface_3D, edge_length,cell_volume,face_surface_2D
class Dv0(object):
	def __init__(self,db):
		self.db=db
		self.mesh=get_mesh(self.db)
		self.position=self.db.get_property('position')

	def get_zrange(self,cid):
		cellpids=[]
		for wid in self.mesh.borders(3,cid):
			for lid in self.mesh.borders(2,wid):
				for pid in self.mesh.borders(1,lid):
					cellpids.append(pid)
		cellpids=list(set(cellpids))
			
		return [min([self.position[pid][2] for pid in cellpids]),max([self.position[pid][2] for pid in cellpids])],cellpids
		
		
	def step(self):

		print 'dv'
		div_zone=200.0
		div_thresh=45.0
		cells_to_divide={}
		divided_props=['orig_cid','cell_type','auxin']
		V=self.db.get_property('V')
		S=self.db.get_property('S')
		scale = self.db.get_property('scale_px_to_micron')
		vindex=self.db.get_property('vindex')
		orig_cid=self.db.get_property('orig_cid')
		orig_pid=self.db.get_property('orig_pid')				
		old_cells=[]
		graph=get_graph(self.db)
		wall=self.db.get_property('wall')

		# get list of cells to divide
		for cid in self.mesh.wisps(3):
			zv,cellpids=self.get_zrange(cid)
			cell_length=zv[1]-zv[0]
			mid_cell=zv[0]+cell_length/2
			if mid_cell<=div_zone and cell_length>=div_thresh and vindex[cid]!=0:
				cells_to_divide[cid]=zv
		print len(cells_to_divide),'cells to divide'

		# divide cells on list
		for cid,zv in cells_to_divide.iteritems():
			#pick division plane and divide cell
			cell_length=zv[1]-zv[0]
			mid_cell=zv[0]+cell_length/2
			zed = round(normal(mid_cell,1.0),4)
			lineage = divide_cell (self.mesh, self.position, cid, array([0,0,zed]), array([0,0,1]))

			# get info from division lineage
			ocid=list(lineage[3].iterkeys())[0]
			old_cells.append(ocid)
			ncell1=list(lineage[3].itervalues())[0][0]
			ncell2=list(lineage[3].itervalues())[0][1]
			zv1,cellpids1=self.get_zrange(ncell1)
			zv2,cellpids2=self.get_zrange(ncell2)

			#update divided props
			for pname in divided_props:
				prop =self.db.get_property(pname)
				prop[ncell1]=prop[ocid]
				prop[ncell2]=prop[ocid]

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
			"""
			if zv1[1]==zv2[0]:
				vindex[ncell1]=vindex[ocid]
				vindex[ncell2]=vindex[ocid]+1
			elif zv2[1]==zv1[0]:
				vindex[ncell2]=vindex[ocid]
				vindex[ncell1]=vindex[ocid]+1
			else:
				print 'non-matching cells. severe difficulties ahead.'
				print 'zv1',zv1,'zv2',zv2
			"""
			zs=sorted([zv1[0],zv1[1],zv2[0],zv2[1]])
			if zs[0] in zv1:
				vindex[ncell1]=vindex[ocid]
				vindex[ncell2]=vindex[ocid]+1
			else:
				vindex[ncell2]=vindex[ocid]
				vindex[ncell1]=vindex[ocid]+1	




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
			if eid not in graph.edges():
				wtd.append(eid)
		for eid in wtd:
			del wall[eid]

		for cid in old_cells:
			for pname in divided_props:
				prop =self.db.get_property(pname)
				del prop[cid]
			del vindex[cid]
			del V[cid]
		std=[]
		for wid in S.iterkeys():
			if wid not in list(self.mesh.wisps(2)):
				std.append(wid)
		for wid in std:
			del S[wid]
		set_vtk_strings(self.db)
		#updateVS(self.db)
