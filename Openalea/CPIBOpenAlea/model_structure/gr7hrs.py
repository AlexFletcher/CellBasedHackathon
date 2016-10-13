
from model_input.import_utils import set_vtk_strings
from model_utils.db_utilities import get_mesh
class Gr(object):
	def __init__(self,db):
		self.db=db

		self.mesh=get_mesh(self.db)
		self.position=self.db.get_property('position')
		"""
		zero_pids=[]
		
		for pid,pos in self.position.iteritems():
			if pos[2]==0.0:
				zero_pids.append(pid)
		self.ccells={}
		for pid in zero_pids:
			for lid in self.mesh.regions(0,pid):
				for wid in self.mesh.regions(1,lid):
					for cid in self.mesh.regions(2,wid):
						if cid not in self.ccells:
							self.ccells[cid]=self.get_zrange(cid)
		"""
					
			


	def get_zrange(self,cid):
		cellpids=[]
		for wid in self.mesh.borders(3,cid):
			for lid in self.mesh.borders(2,wid):
				for pid in self.mesh.borders(1,lid):
					cellpids.append(pid)
			
		return [min([self.position[pid][2] for pid in cellpids]),max([self.position[pid][2] for pid in cellpids])]
		
		
	def step(self):
		
		print 'gr'
		vindex=self.db.get_property('vindex')
		orig_cid=self.db.get_property('orig_cid')

		pos_map={ocid:{} for ocid in orig_cid.itervalues()}
		for cid in self.mesh.wisps(3):
			pos_map[orig_cid[cid]][vindex[cid]]=cid
		
		cvin={ocid:0 for ocid in pos_map.iterkeys()}
		czeds={ocid:self.get_zrange(pos_map[ocid][cvin[ocid]]) for ocid in pos_map.iterkeys()}
		maxvins={ocid:max(list(pos_map[ocid].iterkeys())) for ocid in pos_map.iterkeys()}

		allzeds=sorted(list(set([pos[2] for pos in self.position.itervalues()])))
		gaps=[allzeds[i]-allzeds[i-1] for i in range(1,len(allzeds))]

				
		cpos=allzeds[0]
		i=0
		zshifts={pid:0.0 for pid in self.position.iterkeys()}

			
		while cpos<allzeds[-1]:
			maxL=200.0
			rate=0.1

			Gopt=[]

			for cid,zv in czeds.iteritems():
				cell_length=zv[1]-zv[0]
				mid_cell=zv[0]+cell_length/2
				this_vin=cvin[cid]
				if 0<this_vin<maxvins[cid]:
					Gopt.append(rate*(gaps[i]/(cell_length))*(maxL-(cell_length)))
			if len(Gopt)>0:
				seg_growth=round(min(Gopt),4)
			else:
				seg_growth=0.0

			for pid,pos in self.position.iteritems():
				if pos[2]>cpos:
					zshifts[pid]+=seg_growth
			cpos+=gaps[i]
			new_cid=None
			for cid,zv in czeds.iteritems():
				if zv[1]==cpos:
					if cvin[cid]<maxvins[cid]:
						cvin[cid]+=1
						czeds[cid]=self.get_zrange(pos_map[cid][cvin[cid]])	
			i+=1

		for pid,pos in self.position.iteritems():
			pos[2]+=zshifts[pid]
		
		#set_vtk_strings(self.db)
