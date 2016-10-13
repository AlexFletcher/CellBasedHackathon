from numpy import array, around
from model_input.import_utils import set_vtk_strings
from model_utils.db_utilities import get_mesh,get_graph,set_parameters
from model_utils.db_geom import updateVS
from sa_oa.tissueshape import divide_cell
from numpy.random import normal
from sa_oa.tissueshape import face_surface_3D, edge_length,cell_volume,face_surface_2D
from sa_oa.tissueshape import centroid
import time
class Dv(object):

	name = 'division'

	def __init__(self,db):
		self.db=db
		self.mesh=get_mesh(self.db)
		self.position=self.db.get_property('position')

		self.p={}
		self.p['div_zone']=250.0
		self.p['div_thresh']=('cell_type',{2:10.8,4:9.7,3:14.7,16:23.5,21:23.5,'else':24.4})

		set_parameters(self.db, self.name, self.p)	
		
	def step(self):
		print 'division'
		vtime=time.time()
		div_thresh=self.db.get_property('div_thresh')
		cells_to_divide={}
		# get list of cells to divide
		for cid in self.mesh.wisps(3):
			zv,cellpids=self.get_zrange(cid)
			cell_length=zv[1]-zv[0]
			mid_cell=zv[0]+cell_length/2
			if mid_cell<=self.p['div_zone'] and cell_length>=div_thresh[cid]:
				cells_to_divide[cid]=zv

		self.divide_cells(cells_to_divide)

		print 'div time',time.time()-vtime

	def divide_cells(self,cells_to_divide):
		divided_props=self.db.get_property('divided_props')
		V=self.db.get_property('V')
		S=self.db.get_property('S')
		centres=self.db.get_property('cell_centres')
		scale = self.db.get_property('scale_px_to_micron')
		vindex=self.db.get_property('vindex')
		orig_cid=self.db.get_property('orig_cid')
		orig_pid=self.db.get_property('orig_pid')				
		graph=get_graph(self.db)
		wall=self.db.get_property('wall')
		horiz_walls=self.db.get_property('horiz_walls')

		# divide cells on list
		for cid,zv in cells_to_divide.iteritems():
			#pick division plane and divide cell
			cell_length=zv[1]-zv[0]
			mid_cell=zv[0]+cell_length/2
			zed = round(normal(mid_cell,1.0),4)

			lineage = divide_cell (self.mesh, self.position, cid, array([0,0,zed]), array([0,0,1]))

			# update wall props
			new_wids=[]
			for xwid,ywid in lineage[2].iteritems():
				new_wids.append(ywid[0])
				try:
					new_wids.append(ywid[1])
				except:
					pass
			for xwid,ywid in lineage[2].iteritems():
				if xwid is not None:
					for pname in divided_props['wall']:
						prop = self.db.get_property(pname)
						assert len(ywid)==2
						prop[ywid[0]]=prop[xwid]
						prop[ywid[1]]=prop[xwid]
						if xwid not in new_wids:
							del prop[xwid]
					del S[xwid]
				else:
					for pname in divided_props['wall']:
						prop = self.db.get_property(pname)
						assert len(ywid)==1
						prop[ywid[0]]=0.0
					horiz_walls.append(ywid[0])
					

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
			for pname in divided_props['cell']:
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
					for pname in divided_props['edge']:
						prop =self.db.get_property(pname)
						prop[eid1]=0.0
						prop[eid2]=0.0
			for wid in self.mesh.borders(3,ncell2):
				S[wid] = scale[0]*scale[0]*face_surface_3D(self.mesh, self.position, wid)
				if self.mesh.nb_regions(2,wid) == 2 :
					cid1,cid2 = self.mesh.regions(2,wid)
					if cid1!=ncell1 and cid2!=ncell1:
						eid1 = graph.add_edge(cid1,cid2)
						wall[eid1] = wid
						eid2 = graph.add_edge(cid2,cid1)
						wall[eid2] = wid
						for pname in divided_props['edge']:
							prop =self.db.get_property(pname)
							prop[eid1]=0.0
							prop[eid2]=0.0

			# update V, cell_centres
			V[ncell1]=scale[0]*scale[0]*scale[0]*cell_volume(self.mesh, self.position, ncell1)
			V[ncell2]=scale[0]*scale[0]*scale[0]*cell_volume(self.mesh, self.position, ncell2)
			centres[ncell1]=centroid (self.mesh, self.position, self.mesh.degree(), ncell1)
			centres[ncell2]=centroid (self.mesh, self.position, self.mesh.degree(), ncell2)		

		#tidy up
		wtd=[]
		for eid,wid in wall.iteritems():
			try:
				test=graph.source(eid)
			except:
				wtd.append(eid)
		for eid in wtd:
			del wall[eid]
			for pname in divided_props['edge']:
				prop =self.db.get_property(pname)
				del prop[eid]
			
		
		set_parameters(self.db, self.name, self.p)
		


	def get_zrange(self,cid):
		cellpids=[]
		for wid in self.mesh.borders(3,cid):
			for lid in self.mesh.borders(2,wid):
				for pid in self.mesh.borders(1,lid):
					cellpids.append(pid)
		cellpids=list(set(cellpids))
			
		return [min([self.position[pid][2] for pid in cellpids]),max([self.position[pid][2] for pid in cellpids])],cellpids
