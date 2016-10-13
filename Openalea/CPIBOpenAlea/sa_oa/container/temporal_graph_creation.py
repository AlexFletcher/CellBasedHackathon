# -*- python -*-
#
#       OpenAlea.Container
#
#       Copyright 2012 INRIA - CIRAD - INRA
#
#       File author(s): Jonathan Legrand
#                       Vincent Mirabet
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://sa_oa.gforge.inria.fr
#
################################################################################
"""This module helps to create TemporalPropertyGraph from Spatial Images."""

__license__ = "Cecill-C"
__revision__ = " $Id: temporal_graph_creation.py 14917 2013-09-27 12:28:55Z pradal $ "

import numpy as np
import scipy.ndimage as nd

from sa_oa.container import PropertyGraph

from vplants.tissue_analysis.growth_analysis import dict_cells_slices


def cellNeighbours(im, labels=None, real=False, surf=False):
    """
    calculate cell to cell dictionnary, and optionaly returns a (cell,cell) -> surface dictionnary
    real let you choose between values in voxels or values in um
    :INPUT:
        .im: Spatial Image containing cells (segmented image)
        .labels: if given, return a dict containing only the cells labels provided; if None, will return a dict of all cells present in the spatial image (*Optional*)
    
    :OUTPUT:
    """
    cell_cell={}

    sl=dict_cells_slices(im,labels=labels)
        
    xmax,ymax,zmax=im.shape
    surface={}

    print 'Finding neighbours of cells...'
    for n,i in enumerate(sl.keys()):
        if n%20==0:
            print n,'/',len(sl.keys())
        if i!=1:
            x,y,z=sl[i]
            # On dilate la slice d'un voxel (si on touche pas le bord).
            xd=[x.start-1,x.start][x.start==0]
            xf=[x.stop+1,x.stop][x.stop==xmax-1]
            yd=[y.start-1,y.start][y.start==0]
            yf=[y.stop+1,y.stop][y.stop==ymax-1]
            zd=[z.start-1,z.start][z.start==0]
            zf=[z.stop+1,z.stop][z.stop==zmax-1]
            # On recupere l'info de la SpatialImage contenu dans la slice (creation d'une sous-SpImg)
            mlabel=im[xd:xf,yd:yf,zd:zf].copy()
            mlabel[mlabel!=i]=0
            mlabel[mlabel==i]=1
            m=nd.binary_dilation(mlabel)-mlabel
            res=m*im[xd:xf,yd:yf,zd:zf]
            l=list(np.unique(res))
            l.remove(0)
            cell_cell[i]=l
            if surf:
                for c in l:
                    if real:
                        surface[(i,c)]=len(np.where(res==c)[0])*im.resolution[0]*im.resolution[1]*im.resolution[2]
                        surface[(c,i)]=surface[(i,c)]
                    else:
                        surface[(i,c)]=len(np.where(res==c)[0])
                        surface[(c,i)]=surface[(i,c)]
    
    print 'done !!'
    if surf:
        return cell_cell, surface
    else:
        return cell_cell


def toPropertyGraph(self, cell_cell=None, labels=None):
    """
    creates a PropertyGraph
    """
    p=PropertyGraph()
    p.add_edge_property("source-target")
    # -- If the cell neighbours dictionnary is not provided, we create it.
    if not cell_cell:
        cell_cell=cellNeighbours(self, labels=labels)
    
    if labels!=None:
        l=labels
    else:
        l=cell_cell.keys()
    for c in l:
        if not p.has_vertex(c):
            p.add_vertex(c)
        for cv in cell_cell[c]:
            if (not p.has_vertex(cv)) and (cv in l):
                p.add_vertex(cv)
                eid=p.add_edge(c,cv)
                p.edge_property("source-target")[eid]=(c,cv)
    return p


def display_NXgraph(graph):
    import networkx as nx
    t=nx.Graph()
    
    if type(graph)!=type(t):
        nxGraph = graph.to_networkx()
    else:
        nxGraph=graph
    
    import matplotlib.pyplot as plt
    nx.draw(nxGraph)
    plt.show()


def add_prop2vrtx(graph, dic, prop_name=""):
    """
    Add a property to graph vertices.
    """
    if (prop_name=="") or ('old_label' not in graph.vertex_property_names()):
        # Be carefull : You DON'T have the right to exit a process based on an error. 
        #Just raise an error.
        #sys.exit(1)
        raise ValueError(prop_name)
    if prop_name not in graph.vertex_property_names():
        graph.add_vertex_property(str(prop_name))
    for k,v in graph.vertex_property('old_label').iteritems():
        if (v in dic.keys()) and (v!=1):
            graph.vertex_property(str(prop_name))[k]=dic[v]
    
    return graph


#~ def SegmentedImagetoPropertyGraph(cell_cell):
    #~ """
    #~ creates a PropertyGraph
    #~ """
    #~ p=PropertyGraph()
    #~ p.add_edge_property("source-target")
    #~ for c in cell_cell.keys():
        #~ if not p.has_vertex(c):
            #~ p.add_vertex(c)
        #~ for cv in cell_cell[c]:
            #~ if not p.has_vertex(cv):
                #~ p.add_vertex(cv)
            #~ eid=p.add_edge(c,cv)
            #~ p.edge_property("source-target")[eid]=(c,cv)
    #~ return p


#~ def createTemporalGraph(pgs, links):
    #~ """
    #~ takes propertygraphs and links and create
    #~ """
    #~ g = TemporalPropertyGraph()
    #~ for k,p in enumerate(pgs):
        #~ if not p.vertex_property("time_point"):
            #~ p.add_vertex_property("time_point")
        #~ for i in p.vertex():
            #~ p.vertex_property("time_point")[i]=k            
    #~ return g.extend(pgs,links)
    
