import time
from model_input.import_utils import set_vtk_strings
from model_utils.db_utilities import get_mesh,get_graph

class Gr_unif(object):
	def __init__(self,db):
		self.db=db

		self.mesh=get_mesh(self.db)

		self.position=self.db.get_property('position')



	def remove_cell(self,cid):
		
		cell_props=['vindex','V','auxin','cell_type','orig_cid']
		wall_props=['S']
		point_props=['position','orig_pid']
		edge_props=[]
		mesh=self.mesh
		pos=self.position
		wids=list(mesh.borders(3,cid))
		mesh.remove_wisp(3,cid)
		for pname in cell_props:
			prop=self.db.get_property(pname)
			del prop[cid]			
		for wid in wids:
			if mesh.nb_regions(2,wid)==0:
				lids=list(mesh.borders(2,wid))
				mesh.remove_wisp(2,wid)
				for pname in wall_props:
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
		for pname in edge_props:
			prop=self.db.get_property(pname)
			for eid in etd:
				del prop[eid]
		
		
	def step(self):
		stime=time.time()
		maxL=100.0
		rate_dz=0.0107
		rate_ez=0.1667
		dz=450.0
		end_gz=165.0

		vindex=self.db.get_property('vindex')

		maxv=max([vin for vin in vindex.itervalues()])
		zeds=sorted(list(set([pos[2] for pos in self.position.itervalues()])))

		assert maxv+2==len(zeds)

		all_segs=[0.0]
		cum_growth=0.0
		for i in range(1,maxv):
			cell_length=zeds[i]-zeds[i-1]
			mid_cell=zeds[i-1]+cell_length/2.0
			if mid_cell>dz:
				rate=rate_ez
			else:
				rate=rate_dz
			cum_growth+=rate*(maxL-cell_length)
			all_segs.append(cum_growth)
		all_segs.append(cum_growth)	
		gr_dict=dict(zip(zeds[1:],all_segs))
	
		for pid,pos in self.position.iteritems():
			if pos[2]>0.0:
				pos[2]+=round(gr_dict[pos[2]],4)

		
		remove=[]
		vzeds={}

		for ozpos in sorted(gr_dict):
			vzeds[len(vzeds)]=ozpos+gr_dict[ozpos]
		print vzeds

		for xcid in self.mesh.wisps(3):		
			if vzeds[vindex[xcid]]>end_gz:
				remove.append(xcid)	
		for cid in remove:
			self.remove_cell(cid)

		print 'gtime',time.time()-stime

		#set_vtk_strings(self.db)
