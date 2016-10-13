# -*- coding: utf-8 -*-
#
#       Graph : graph package
#
#       Copyright or Copr. 2006 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       VPlants WebSite : https://gforge.inria.fr/projects/vplants/
#
"""
This module provide a simple pure python implementation
for a graph interface
do not implement copy concept
"""

__license__= "Cecill-C"
__revision__=" $Id: graph.py 12282 2012-06-21 16:09:29Z boudon $ "

from interface.graph import InvalidEdge,InvalidVertex,IGraph,\
                                        IVertexListGraph,IEdgeListGraph,\
                                        IMutableVertexGraph,IMutableEdgeGraph,\
                                        IExtendGraph
from utils import IdDict

class Graph (IGraph,\
                        IVertexListGraph,IEdgeListGraph,\
                        IMutableVertexGraph,IMutableEdgeGraph,\
                        IExtendGraph):
    """
    directed graph with multiple links
    in this implementation :

        - vertices are tuple of edge_in,edge_out
        - edges are tuple of source,target
    """
    def __init__(self, graph=None, idgenerator = "set"):
        """constructor

        if graph is not none make a copy of the topological structure of graph
        (i.e. don't use the same id)

        :param graph: the graph to copy, default=None
        :type graph: Graph
        """
        self._vertices=IdDict(idgenerator = idgenerator)
        self._edges=IdDict(idgenerator = idgenerator)
        if graph is not None :
            dummy=self.extend(graph)

    # ##########################################################
    #
    # Graph concept
    #
    # ##########################################################
    def source(self, eid):
        try :
            return self._edges[eid][0]
        except KeyError :
            raise InvalidEdge(eid)
    source.__doc__=IGraph.source.__doc__

    def target(self, eid):
        try :
            return self._edges[eid][1]
        except KeyError :
            raise InvalidEdge(eid)
    target.__doc__=IGraph.target.__doc__

    def edge_vertices(self, eid):
        try :
            return self._edges[eid]
        except KeyError :
            raise InvalidEdge(eid)
    edge_vertices.__doc__=IGraph.edge_vertices.__doc__
    
    def edge(self, source, target) :
        link_in,link_out=self._vertices[source]
        for eid in link_in : 
            if self._edges[eid][0] == target: 
                return eid
        for eid in link_out :
            if self._edges[eid][1] == target: 
                return eid        
        return None
        
    edge.__doc__=IGraph.edge.__doc__

    def __contains__(self, vid):
        return self.has_vertex(vid)
    __contains__.__doc__=IGraph.__contains__.__doc__

    def has_vertex(self,vid):
        return self._vertices.has_key(vid)
    has_vertex.__doc__=IGraph.has_vertex.__doc__

    def has_edge(self,eid):
        return self._edges.has_key(eid)
    has_edge.__doc__=IGraph.has_edge.__doc__

    def is_valid(self):
        return True
    is_valid.__doc__=IGraph.is_valid.__doc__

    # ##########################################################
    #
    # Vertex List Graph Concept
    #
    # ##########################################################
    def vertices(self):
        return iter(self._vertices)
    vertices.__doc__=IVertexListGraph.vertices.__doc__

    def __iter__ (self) :
        return iter(self._vertices)
    __iter__.__doc__=IVertexListGraph.__iter__.__doc__

    def nb_vertices(self):
        return len(self._vertices)
    nb_vertices.__doc__=IVertexListGraph.nb_vertices.__doc__

    def __len__(self):
        return self.nb_vertices()
    __len__.__doc__=IVertexListGraph.__len__.__doc__

    def in_neighbors(self, vid):
        if vid not in self :
            raise InvalidVertex(vid)
        neighbors_list=[self.source(eid) for eid in self._vertices[vid][0] ]
        return iter(set(neighbors_list))
    in_neighbors.__doc__=IVertexListGraph.in_neighbors.__doc__

    def out_neighbors(self, vid):
        if vid not in self :
            raise InvalidVertex(vid)
        neighbors_list=[self.target(eid) for eid in self._vertices[vid][1] ]
        return iter(set(neighbors_list))
    out_neighbors.__doc__=IVertexListGraph.out_neighbors.__doc__

    def neighbors(self, vid):
        neighbors_list=list(self.in_neighbors(vid))
        neighbors_list.extend(self.out_neighbors(vid))
        return iter(set(neighbors_list))
    neighbors.__doc__=IVertexListGraph.neighbors.__doc__

    def nb_in_neighbors(self, vid):
        neighbors_set=list(self.in_neighbors(vid))
        return len(neighbors_set)
    nb_in_neighbors.__doc__=IVertexListGraph.nb_in_neighbors.__doc__

    def nb_out_neighbors(self, vid):
        neighbors_set=list(self.out_neighbors(vid))
        return len(neighbors_set)
    nb_out_neighbors.__doc__=IVertexListGraph.nb_out_neighbors.__doc__

    def nb_neighbors(self, vid):
        neighbors_set=list(self.neighbors(vid))
        return len(neighbors_set)
    nb_neighbors.__doc__=IVertexListGraph.nb_neighbors.__doc__

    # ##########################################################
    #
    # Edge List Graph Concept
    #
    # ##########################################################
    def _iter_edges (self, vid) :
        """
        internal function that perform 'edges' with vid not None
        """
        link_in,link_out=self._vertices[vid]
        for eid in link_in : yield eid
        for eid in link_out : yield eid

    def edges(self, vid=None):
        if vid is None :
            return iter(self._edges)
        if vid not in self :
            raise InvalidVertex(vid)
        return self._iter_edges(vid)
    edges.__doc__=IEdgeListGraph.edges.__doc__

    def nb_edges(self, vid=None):
        if vid is None :
            return len(self._edges)
        if vid not in self :
            raise InvalidVertex(vid)
        return len(self._vertices[vid][0])+len(self._vertices[vid][1])
    nb_edges.__doc__=IEdgeListGraph.nb_edges.__doc__

    def in_edges(self, vid):
        if vid not in self :
            raise InvalidVertex(vid)
        for eid in self._vertices[vid][0] : yield eid
    in_edges.__doc__=IEdgeListGraph.in_edges.__doc__

    def out_edges(self, vid):
        if vid not in self :
            raise InvalidVertex(vid)
        for eid in self._vertices[vid][1] : yield eid
    out_edges.__doc__=IEdgeListGraph.out_edges.__doc__

    def nb_in_edges(self, vid):
        if vid not in self :
            raise InvalidVertex(vid)
        return len(self._vertices[vid][0])
    nb_in_edges.__doc__=IEdgeListGraph.nb_in_edges.__doc__

    def nb_out_edges(self, vid):
        if vid not in self :
            raise InvalidVertex(vid)
        return len(self._vertices[vid][1])
    nb_out_edges.__doc__=IEdgeListGraph.nb_out_edges.__doc__

    # ##########################################################
    #
    # Mutable Vertex Graph concept
    #
    # ##########################################################
    def add_vertex(self, vid=None):
        return self._vertices.add((set(),set()),vid )
    add_vertex.__doc__=IMutableVertexGraph.add_vertex.__doc__

    def remove_vertex(self, vid):
        if vid not in self :
            raise InvalidVertex(vid)
        link_in,link_out=self._vertices[vid]
        for edge in list(link_in) : self.remove_edge(edge)
        for edge in list(link_out) : self.remove_edge(edge)
        del self._vertices[vid]
    remove_vertex.__doc__=IMutableVertexGraph.remove_vertex.__doc__

    def clear(self):
        self._edges.clear()
        self._vertices.clear()
    clear.__doc__=IMutableVertexGraph.clear.__doc__

    # ##########################################################
    #
    # Mutable Edge Graph concept
    #
    # ##########################################################
    def add_edge(self, sid, tid, eid=None):
        if sid not in self :
            raise InvalidVertex(sid)
        if tid not in self :
            raise InvalidVertex(tid)
        eid=self._edges.add((sid,tid),eid)
        self._vertices[sid][1].add(eid)
        self._vertices[tid][0].add(eid)
        return eid
    add_edge.__doc__=IMutableEdgeGraph.add_edge.__doc__

    def remove_edge(self,eid):
        if not self.has_edge(eid) :
            raise InvalidEdge(eid)
        sid,tid=self._edges[eid]
        self._vertices[sid][1].remove(eid)
        self._vertices[tid][0].remove(eid)
        del self._edges[eid]
    remove_edge.__doc__=IMutableEdgeGraph.remove_edge.__doc__

    def clear_edges(self):
        self._edges.clear()
        for vid,(in_edges,out_edges) in self._vertices.iteritems() :
            in_edges.clear()
            out_edges.clear()
    clear_edges.__doc__=IMutableEdgeGraph.clear_edges.__doc__
    
    # ##########################################################
    #
    # Extend Graph concept
    #
    # ##########################################################
    def extend(self, graph):
        #vertex adding
        trans_vid={}
        for vid in list(graph.vertices()) :
            trans_vid[vid]=self.add_vertex()

        #edge adding
        trans_eid={}
        for eid in list(graph.edges()) :
            sid=trans_vid[graph.source(eid)]
            tid=trans_vid[graph.target(eid)]
            trans_eid[eid]=self.add_edge(sid,tid)

        return trans_vid,trans_eid
    extend.__doc__=IExtendGraph.extend.__doc__
    
    def sub_graph(self, vids):
        """
        """
        from copy import deepcopy
        vids = set(vids)
        
        result = deepcopy(self)
        result._vertices.clear()
        result._edges.clear()
        
        for key, edges in self._vertices.items():
            if key in vids:
                inedges, outedges = edges
                sortedinedges = set([eid for eid in inedges if self.source(eid) in vids])
                sortedoutedges = set([eid for eid in outedges if self.target(eid) in vids])
                result._vertices.add((sortedinedges,sortedoutedges), key)
                for eid in sortedoutedges:
                    result._edges.add(self._edges[eid], eid)
        
        return result
