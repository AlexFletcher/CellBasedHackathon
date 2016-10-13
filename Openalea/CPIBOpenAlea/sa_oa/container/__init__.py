# CPL: Do not import anything in the __init__
# Otherelse it may breaks everything see commit  15337

from utils import IdDict
from data_prop import Quantity,DataProp
from graph import Graph
from property_graph import PropertyGraph
from temporal_property_graph import TemporalPropertyGraph
from tree import Tree, PropertyTree
from grid import Grid
from relation import Relation

################################
#
#       mesh
#
################################
from topomesh import *
from topomesh_txt import write_topomesh,read_topomesh
from topomesh_algo import *
from topomesh_geom_algo import *
