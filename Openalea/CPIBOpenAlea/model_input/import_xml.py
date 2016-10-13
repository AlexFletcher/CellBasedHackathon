# Program by John Fozard.
# Converts xml output of CellSeT to an OpenAlea zip data structure.
# Note, CellSeT uses an origin at the top left hand corner of the image, whereas OpenAlea has the origin at the bottom left hand corner.
# CellSeT gives geometries in terms of pixels, whereas OpenAlea works in microns - depending on the microscope settings, weprescribe 'scale' on line 100, to convert between the two.

import lxml.etree as ET

import numpy as np
import matplotlib.pylab as plt
from matplotlib.collections import PolyCollection, CircleCollection


from sa_oa.container import Topomesh, Quantity, ordered_pids
from sa_oa.celltissue import Tissue, topen, ConfigFormat, TissueDB, \
    ConfigItem
#from geom import perimeter, area
from sa_oa.tissueshape import face_surface_2D,edge_length

from import_utils import get_2d_offset_vectors

# Routine to convert xml attributes to strings or numbers,
# as appropriate
def to_val(s):
    s=str(s)
    if s=='None':
        return None
    else:
        try:
            s=int(s)
        except:
            try:
                s=float(s)
            except:
                pass
            pass
        return s


def import_xml(filename):
    wall_data={}
    cell_data=[]

    p_map={}

    t=Tissue()
    POINT = t.add_type('POINT')
    WALL = t.add_type('WALL')
    CELL = t.add_type('CELL')
    mesh_id = t.add_relation('mesh', (POINT,WALL,CELL) )
    mesh = t.relation(mesh_id)


    tree=ET.parse(filename)
    root=tree.getroot()
    cells=tree.find('cells')
    cell_list = []
    cell_group_list = []
    old_cell_id_list = []
    truncated_list = []
    for c in cells.iter('cell'):
        walls=c.find('walls')
        cell_walls=[]
        for w in walls.iter('wall'):
            cell_walls.append(int(w.get('id')))
        cell_data.append(cell_walls)
        cell_group_list.append(int(c.get('group')))
        old_cell_id_list.append(int(c.get('id')))
        truncated_list.append(c.get('truncated') == 'true')
        
  
    walls=tree.find('walls')
    wc = 0
    for w in walls.iter('wall'):
        pts=[]
        for point in w.find('points').iter('point'):
            pts.append((float(point.get('x')), float(point.get('y'))))
        wall_data[int(w.get('id'))]=pts
#        wall_data[int(w.get('id'))]=[pts[0], pts[-1]]
        wc+=1


    pt_map={}
    wall_map={}
    pos={}
    for i, w in wall_data.iteritems():
        wall_map[i]=[]
        for pp in zip(w[0:-1], w[1:]):
            wid=mesh.add_wisp(1)
            for p in pp:
                if p not in pt_map:
                    pid=mesh.add_wisp(0)
                    pos[pid]=p
                    pt_map[p]=pid
                mesh.link(1, wid, pt_map[p])
            wall_map[i].append(wid)

    y_min = min(v[1] for pid, v in pos.iteritems())
    y_max = max(v[1] for pid, v in pos.iteritems())
    x_min = min(v[0] for pid, v in pos.iteritems())
    x_max = max(v[0] for pid, v in pos.iteritems())
    print 'x',x_min,x_max
    print 'y',y_min,y_max
    scale = 0.0943
    for pid, v in pos.iteritems():
        pos[pid]=(v[0]*scale, (y_min+(y_max-v[1]))*scale)

    cid_list=[]
    for c in cell_data:
        cid=mesh.add_wisp(2)
        cid_list.append(cid)
        for w in c:
            for wid in wall_map.get(w,[]):
                mesh.link(2, cid, wid)

    EDGE = t.add_type('EDGE')
    graph_id = t.add_relation("graph",(CELL,EDGE) )
    graph = t.relation(graph_id)

    wall = {}

    for wid in mesh.wisps(1) :
	if mesh.nb_regions(1,wid) == 2 :
		cid1,cid2 = mesh.regions(1,wid)
		eid1 = graph.add_edge(cid1,cid2)
		wall[eid1] = wid
		eid2 = graph.add_edge(cid2,cid1)
		wall[eid2] = wid

    

    cell_type = dict((cid, cell_group_list[i]) 
                      for i, cid in enumerate(cid_list))

   # cell_types =dict((i,str(i)) for i in set(cell_group_list))
    cell_types = { 2: 'epiderm',
                   4: 'cortex',
                   3: 'endoderm',
                   5: 'vasculature',
                   6: 'LRC1',
                   7: 'LRC2',
                   8: 'LRC3',
		   9: 'LRC4',
                   17: 'QC',
                   10: 'initials',
                   11: 'S1',
                   12: 'S2',
                   13: 'S3',
                   14: 'S4',
                   15: 'S5',
		   16: 'pericycle',
		   19: 'protoxylem',
		   20: 'metaxylem',
		   21: 'xp_pericycle',
		   22: 'procambium',
		   23: 'phloem'}

                   

    old_cell_id = dict((cid, old_cell_id_list[i]) 
                      for i, cid in enumerate(cid_list))

    border_dict = dict((cid_list[i], 1) 
                  for i,t in enumerate(truncated_list) if t)


    file_cell_types = map(cell_types.get, ['epiderm', 'cortex', 'endoderm'])
    edge_type = {}
    for cid in mesh.wisps(2):
        oe = graph.out_edges(cid)
        nb_cells = set(graph.target(eid) for eid in oe)
        if cell_type[cid] in file_cell_types:
            nb_file = [nid for nid in nb_cells if cell_type[nid] == cell_type[cid]]
                
   
    V = {}
    for cid in mesh.wisps(2):
        V[cid] = face_surface_2D(mesh, pos, cid)
   
    S = {}
    for wid in mesh.wisps(1):
        S[wid] = edge_length(mesh, pos, wid)

    offset_vectors=get_2d_offset_vectors(mesh,pos)

    pos = Quantity(pos,'pix','tuple',"position of vertices")
    cfg = ConfigFormat(locals())
    cfg.add_section("tissue descr")
    cfg.add('POINT')
    cfg.add('WALL')
    cfg.add('CELL')
    cfg.add('EDGE')
    cfg.add('mesh_id')
    cfg.add('graph_id')
    cfg.add_section("prop types")
    cfg.add('cell_types')
    cfg.add_item(ConfigItem('border', {1:'truncated'}))

    species_desc={'V':'cell','S':'wall','cell_type':'cell'}

    #writing

    db=TissueDB()
    db.set_config('config', cfg.config())
    db.set_tissue(t)
    db.set_property('position', pos)
    db.set_description('position', "position")
    db.set_property('old_cell_id', old_cell_id)
    db.set_description('old_cell_id', "original cell id numbers in xml file")
    db.set_property('cell_type', cell_type)
    db.set_description('cell_type', "type of each cell")
    db.set_property('border', border_dict)
    db.set_description('border',"")
    db.set_property('wall', wall)
    db.set_description('wall', "wall corresponding to a given edge")
    db.set_property('V', V)
    db.set_description('V', "Volume of each cell")
    db.set_property('S', S)
    db.set_description('S', "Surface of each wall in m")
    db.set_property('offset_vectors', offset_vectors)
    db.set_description('offset_vectors', "horizontal offset vectors to shrink cells etc.")
    db.set_property("species_desc",species_desc)
    db.set_description("species_desc","association of properties with mesh degree")
    return db


