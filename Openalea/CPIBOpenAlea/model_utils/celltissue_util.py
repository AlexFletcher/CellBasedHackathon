
from sa_oa.celltissue import TissueDB
from numpy import array

class DummyFunc (object) :
        def __init__ (self, val) :
                self._val = val

        def __call__ (self) :
                return self._val

def tissue (db) :
        return db,db.tissue()

def get_property (db, name) :
        return db,db.get_property(name)

def get_config (db, name) :
        return db,db.get_config(name)

def get_topology (db, name, cfg) :
        return db,db.get_topology(name,cfg)

#########################################                                       
#                                                                               
#       edition                                                                 
#                                                                               
#########################################        
def def_property (db, name, val, elm_type, cfg, descr) :
        cfg = db.get_config(cfg)
        prop = {}

	if callable(val) :
                func = val
        else :
                func = DummyFunc(val)

	for elmid in db._tissue.elements(cfg[elm_type]) :
                prop[elmid] = func()

        db.set_property(name,prop)
        db.set_description(name,descr)
	#try:
	#	species_desc=db.get_property('species_desc')
	#	species_desc[name]=elm_type
	#except:
	#	species_desc={name:elm_type}
	#	db.set_property('species_desc',species_desc)
	#	db.set_description('species_desc', 'type of property, i.e. cell,edge or wall') 

        return db,prop

def scaled_property (db, name, ref_prop, scale) :
        prop = {}

        for key,val in ref_prop.iteritems() :
                prop[key] = val * scale

        db.set_property(name,prop)

        return db,prop

def merge_property (prop, ref_prop) :
        """Update the content of prop with ref_prop.                            
        """
        prop.update(ref_prop)

        return prop,
def read (filename) :
        #read                                                                   
        db = TissueDB()
        db.read(filename)
	#transform positions        
#        try :
#             	pos = db.get_property("position")
#                db.set_property("position",array(pos) )#
#	except KeyError :
#                pass
        #return                                                                 
	return db,

def write (db, filename) :
        try :
                #transform positions                                            
		pos = db.get_property("position")
                db.set_property("position",totup(pos) )
		print 'totup try'
		#write                                                          
                db.write(filename)
                #put back pos                                                   
                db.set_property("position",pos)
        except KeyError :
		#write                                                          
                db.write(filename)
        #return                                                                 
        return db,


def CellToCoordinates(mesh,cell_to_walls,border_to_coordinates):

    ptc=[]
    for ncid,cid in enumerate(mesh.wisps(2)):
        ptc.append([])
        patches=[]
  
        for i in range(len(cell_to_walls[cid])):
            patches.append(border_to_coordinates[cell_to_walls[cid][i]])

        ptc[ncid]=[]
        ptc[ncid].append(patches[0][0])
        ptc[ncid].append(patches[0][1])
        patches.remove(patches[0])
        i=0
        while i <len(patches):
            if patches[i].count(ptc[ncid][-1])==1:
                id1=patches[i].index(ptc[ncid][-1])
		id2=abs(-1+id1)
                ptc[ncid].append(patches[i][id2])
                patches.remove(patches[i])
                i=-1
            i=i+1

    cell_to_coordinates = dict((cid,ptc[ncid]) for ncid,cid in enumerate(mesh.wisps(2)))

    return cell_to_coordinates


