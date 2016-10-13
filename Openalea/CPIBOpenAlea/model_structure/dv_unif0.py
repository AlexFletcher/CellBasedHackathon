from numpy import array, around
from model_input.import_utils import set_vtk_strings
from model_utils.db_utilities import get_mesh,get_graph
from model_utils.db_geom import updateVS
from sa_oa.tissueshape import divide_cell
from numpy.random import normal
from sa_oa.celltissue import Tissue,TissueDB
from sa_oa.tissueshape import face_surface_3D, edge_length,cell_volume,face_surface_2D
from sa_oa.celltissue import ConfigFormat,ConfigItem
from model_output.vtk_output import *


class Dv_unif(object):
	def __init__(self,db):
		self.db=db
		self.mesh=get_mesh(self.db)
		self.position=self.db.get_property('position')

		tissue = Tissue()
		POINT = tissue.add_type("point")
		LINE = tissue.add_type("line")
		WALL = tissue.add_type("wall")
		CELL = tissue.add_type("cell")
		EDGE = tissue.add_type("edge")

		mesh_id = tissue.add_relation("mesh",(POINT,LINE,WALL,CELL) )
		graph_id = tissue.add_relation("graph",(CELL,EDGE) )
		
		orig_pid=self.db.get_property('orig_pid')

		self.mesh0 = tissue.relation(mesh_id)
		self.pos0 = {}
		self.opids = {}

		ptrans={}
		ltrans={}
		wtrans={}
		#points
		for pid,pos in self.position.iteritems():
			if pos[2]==0.0:
				npid = self.mesh0.add_wisp(0)
				self.pos0[npid]=pos
				ptrans[pid]=npid
				self.opids[npid]=orig_pid[pid]
		#walls
		for opid,npid in ptrans.iteritems():
			for lid in self.mesh.regions(0,opid):
				xzeds=[]
				for xpid in self.mesh.borders(1,lid):
					xzeds.append(self.position[xpid][2])
				if max(xzeds)==0.0:
					print 'k'					
					if lid not in ltrans:
						nwid = self.mesh0.add_wisp(1)
						self.mesh0.link(1,nwid,npid)
						ltrans[lid]=nwid
					else:
						self.mesh0.link(1,ltrans[lid],npid)
		#cells
		for olid,nwid in ltrans.iteritems():
			for wid in self.mesh.regions(1,olid):
				xzeds=[]
				for xlid in self.mesh.borders(2,wid):
					for xpid in self.mesh.borders(1,xlid):			
						xzeds.append(self.position[xpid][2])
				if max(xzeds)==0.0:
					print 'k'
					if wid not in wtrans:
						ncid = self.mesh0.add_wisp(2)
						self.mesh0.link(2,ncid,nwid)
						wtrans[wid]=ncid
					else:
						self.mesh0.link(2,wtrans[wid],nwid)			
						
		cfg = ConfigFormat(vars() )
		cfg.add_section("elements types")
		cfg.add("POINT")
		cfg.add("LINE")
		cfg.add("WALL")
		cfg.add("CELL")
		cfg.add("EDGE")
		cfg.add("mesh_id")
		cfg.add("graph_id")
		#cfg.add("cell_types")


		cfg = cfg.config()

		place={cid:0.0 for cid in self.mesh0.wisps(2)}
		ndb=TissueDB()
		ndb.set_tissue(tissue)
		ndb.set_config('config',cfg)
		dbprops={"position":(self.pos0,"position of points"),"species_desc":({},"association of properties with mesh degree"),"vindex":({},"position in z direction"),\
		"orig_cid":({},"cell id from orig 2d template"),"cell_type":({},"type of each cell"),"wall":({},"wall corresponding to a given edge"),\
		"V":({},"Volume of each cell in m2"),"S":({},"Surface of each wall in m"),\
		\
		"orig_pid":({},"point id from orig 2d template"),"place":(place, "horizontal offset vectors to shrink cells etc.")}
		for pname,(prop,desc) in dbprops.iteritems():
			ndb.set_property(pname,prop)
			ndb.set_description(pname,desc)
		from model_output.pylab_plot import *
		plot_pylab_figure(ndb,-1,'test.png',['place'],[0.0,1.0])

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
		for vin in range(1,maxv+2):
			cell_length=zeds[vin]-zeds[vin-1]
			print 'cell',cell_length
			mid_cell=zeds[vin-1]+cell_length/2.0
			print 'mc',mid_cell
			if mid_cell<=div_zone and cell_length>=div_thresh:
				vins_to_divide.append(vin)
				vzeds.append([zeds[vin-1],zeds[vin]])
		vins_to_divide.reverse()
		vzeds.reverse()
		print vins_to_divide,'vins to divide'
		print vzeds,'vzeds'



		pos_map={ocid:{} for ocid in orig_cid.itervalues()}
		for cid in self.mesh.wisps(3):
			pos_map[orig_cid[cid]][vindex[cid]]=cid

		# divide cells on list
		for i,vin in enumerate(vins_to_divide):
			zrange=vzeds[i]
			# shift points
			for pid,pos in self.position.iteritems():
				if pos[2]>zrange[0]:
					pos[2]-=round((zrange[1]-zrange[0])/2.0,4)

			# add new vindex of cells
			db_cids={}
			db_wids={}
			db_lids={}
			db_pids={}
			dwids=[]

			#cells
			for cid in self.mesh0.wisps(2):
				ncid = self.mesh.add_wisp(3)
				db_cids[cid]=ncid
			#horiz walls
			for cid in self.mesh0.wisps(2):
				nwid = self.mesh.add_wisp(2)
				db_wids[cid]=nwid
				self.mesh.link(3,db_cids[cid],nwid)
			#vert walls
			for wid in self.mesh0.wisps(1):
				nwid = self.mesh.add_wisp(2)
				db_wids[wid]=nwid
				for cid in self.mesh0.regions(1,wid):
					self.mesh.link(3,db_cids[cid],nwid)
			#horiz lines
			for wid in self.mesh0.wisps(1):
				nlid = self.mesh.add_wisp(1)
				db_lids[wid]=nlid	
				for cid in self.mesh0.regions(1,wid):
					self.mesh.link(2,db_wids[cid],nlid)
				self.mesh.link(2,db_wids[wid],nlid)

			#vert lines
			for pid in self.mesh0.wisps(0):
				nlid = self.mesh.add_wisp(1)
				db_lids[pid]=nlid
				for wid in self.mesh0.regions(0,pid):
					self.mesh.link(2,db_wids[wid],nlid)

			#points
			for pid in self.mesh0.wisps(0):
				npid = self.mesh.add_wisp(0)
				orig_pid[npid]=self.opids[pid]
				db_pids[pid]=npid
				for wid in self.mesh0.regions(0,pid):
					self.mesh.link(1,db_lids[wid],npid)
				
			#add positions		
			for pid in self.mesh0.wisps(0):
				opos=self.pos0[pid]
				self.position[db_pids[pid]]=[opos[0],opos[1],zeds[-1]]

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


