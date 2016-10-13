import time
from model_input.import_utils import set_vtk_strings
from model_utils.db_utilities import get_mesh,get_graph
from sa_oa.tissueshape import divide_cell
from numpy import array
from sa_oa.tissueshape import face_surface_3D, edge_length,cell_volume,face_surface_2D
class Gr3_cut(object):

	name = 'gr3_cut'
	
	def __init__(self,db):
		self.db=db

		self.mesh=get_mesh(self.db)

		self.position=self.db.get_property('position')

		intime=time.time()
	

	def get_zrange(self,cid):
		cellpids=[]
		for wid in self.mesh.borders(3,cid):
			for lid in self.mesh.borders(2,wid):
				for pid in self.mesh.borders(1,lid):
					cellpids.append(pid)
		cellpids=list(set(cellpids))
			
		return [min([self.position[pid][2] for pid in cellpids]),max([self.position[pid][2] for pid in cellpids])],cellpids			



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

	def divide_cell_full(self,cid,divz):
		divided_props=self.db.get_property('divided_props')
		V=self.db.get_property('V')
		S=self.db.get_property('S')
		scale = self.db.get_property('scale_px_to_micron')
		vindex=self.db.get_property('vindex')
		orig_cid=self.db.get_property('orig_cid')
		orig_pid=self.db.get_property('orig_pid')				
		graph=get_graph(self.db)
		wall=self.db.get_property('wall')
		horiz_walls=self.db.get_property('horiz_walls')

		#pick division plane and divide cell

		zed = round(divz,4)

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
			for pname in divided_props['edge']:
				prop =self.db.get_property(pname)
				del prop[eid]

		return lineage

		
		
	def step(self):
		pt=time.time()
		print 'gr'
		vindex=self.db.get_property('vindex')
		orig_cid=self.db.get_property('orig_cid')

		pos_map={ocid:{} for ocid in orig_cid.itervalues()}
		all_czeds={}
		for cid in self.mesh.wisps(3):
			pos_map[orig_cid[cid]][vindex[cid]]=cid
			all_czeds[cid],cps=self.get_zrange(cid)

		cvin={ocid:0 for ocid in pos_map.iterkeys()}
		czeds={}
		for ocid in pos_map.iterkeys():
			zv,cps=self.get_zrange(pos_map[ocid][cvin[ocid]])
			czeds[ocid]=zv
		maxvins={ocid:max(list(pos_map[ocid].iterkeys())) for ocid in pos_map.iterkeys()}
	
		allzeds=[]
		
		for pid,pos in self.position.iteritems():
			allzeds.append(pos[2])
			
		allzeds=sorted(list(set(allzeds)))
		gaps=[allzeds[i]-allzeds[i-1] for i in range(1,len(allzeds))]

		cpos=allzeds[0]
		i=0
		all_segs=[]
		cum_growth=0.0
		while cpos<allzeds[-1]:
			
			maxL=100.0
			rate_dz=0.0107
			rate_ez=0.1667
			dz=450.0
			end_gz=50.0

			Gopt=[]

			for cid,zv in czeds.iteritems():
				cell_length=zv[1]-zv[0]
				mid_cell=zv[0]+cell_length/2.0
				if mid_cell>dz:
					rate=rate_ez
				else:
					rate=rate_dz
				this_vin=cvin[cid]
				if 0<this_vin<maxvins[cid]:
					Gopt.append(rate*(gaps[i]/(cell_length))*(maxL-(cell_length)))
			if len(Gopt)>0:
				seg_growth=min(Gopt)
			else:
				seg_growth=0.0

			#zero growth in zero cells
			for cid,zv in czeds.iteritems():
				if cvin[cid]==maxvins[cid]:#cvin[cid]==0 or 
					seg_growth=0.0
			
			cum_growth+=seg_growth
			all_segs.append(cum_growth)

			cpos+=gaps[i]
			new_cid=None
			for cid,zv in czeds.iteritems():
				if zv[1]==cpos:
					if cvin[cid]<maxvins[cid]:
						cvin[cid]+=1
						czeds[cid]=all_czeds[pos_map[cid][cvin[cid]]]

	
			i+=1


		gr_dict=dict(zip(allzeds[1:],all_segs))
	
		for pid,pos in self.position.iteritems():
			if pos[2]>0.0:
				pos[2]+=round(gr_dict[pos[2]],4)
		print 'growth time',time.time()-pt
		allwisps3=list(self.mesh.wisps(3))
		count=0
		trim=0
		for xcid in allwisps3:
			zr=all_czeds[xcid]
			if zr[0]>=end_gz:
				self.remove_cell(xcid)
				count+=1
			elif zr[1]>end_gz:
				lineage = self.divide_cell_full(xcid,end_gz)
				ncell1=list(lineage[3].itervalues())[0][0]
				ncell2=list(lineage[3].itervalues())[0][1]
				if vindex[ncell1]==vindex[ncell2]:
					print 'equal vins'
				if vindex[ncell1]>vindex[ncell2]:
					self.remove_cell(ncell1)
				else:
					self.remove_cell(ncell2)
				trim+=1
		print 'count',count
		print 'trim',trim
		#set_vtk_strings(self.db)
