from model_structure.sim_dg import Sim_dg
from delta_notch import delta_notch
from model_utils.db_utilities import get_mesh

import time
stime = time.time()

sim=Sim_dg('full_root_2dt.zip', genenetwork=delta_notch)
sim.set_timestep(100)
db=sim.get_db()

delta=db.get_property('delta')
cell_type=db.get_property('cell_type')
mesh=get_mesh(db)
for cid in mesh.wisps(2):
	if cell_type[cid] in [19,23,2]:
		delta[cid]=10
sim.output()
sim.run_simulation(0,50)



print 'TOTAL TIME:',time.time()-stime





