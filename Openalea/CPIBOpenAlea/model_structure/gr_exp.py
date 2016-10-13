import time
from model_input.import_utils import set_vtk_strings
from model_utils.db_utilities import get_mesh,get_graph
class Gr_exp(object):

	name = 'gr_exp'
	
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
			
		return [min([self.position[pid][2] for pid in cellpids]),max([self.position[pid][2] for pid in cellpids])]

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
		pt=time.time()
		print 'gr'
		vindex=self.db.get_property('vindex')
		orig_cid=self.db.get_property('orig_cid')

		pos_map={ocid:{} for ocid in orig_cid.itervalues()}
		all_czeds={}
		for cid in self.mesh.wisps(3):
			pos_map[orig_cid[cid]][vindex[cid]]=cid
			all_czeds[cid]=self.get_zrange(cid)

		cvin={ocid:0 for ocid in pos_map.iterkeys()}
		czeds={ocid:self.get_zrange(pos_map[ocid][cvin[ocid]]) for ocid in pos_map.iterkeys()}
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
			end_gz=150.0

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
		for xcid in allwisps3:
			zr=all_czeds[xcid]
			if zr[0]+(zr[1]-zr[0])*0.5>end_gz:
				self.remove_cell(xcid)
				count+=1
		print 'count',count
		#set_vtk_strings(self.db)
