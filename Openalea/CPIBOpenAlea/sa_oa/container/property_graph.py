# -*- python -*-
#
#       OpenAlea.Core
#
#       Copyright 2006-2009 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>
#                       Fred Theveny <frederic.theveny@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://sa_oa.gforge.inria.fr
#
################################################################################
"""This module provide a set of concepts to add properties to graph elements"""

__license__ = "Cecill-C"
__revision__ = " $Id: property_graph.py 17312 2014-08-05 11:18:52Z jlegra02 $ "

from interface.property_graph import IPropertyGraph, PropertyError
from graph import Graph, InvalidVertex, InvalidEdge
import numpy as np
from heapq import heappop, heappush
import warnings

VertexProperty, EdgeProperty, GraphProperty = range(3)
VertexIdType, EdgeIdType, ValueType = range(3)

class PropertyGraph(IPropertyGraph, Graph):
    """
    Simple implementation of IPropertyGraph using
    dict as properties and two dictionaries to
    maintain these properties
    """
    metavidtypepropertyname = "valueproperty_as_vid"
    metaeidtypepropertyname = "valueproperty_as_eid"
    
    def __init__(self, graph=None, **kwds):
        self._vertex_property = {}
        self._edge_property = {}
        self._graph_property = {}
        Graph.__init__(self, graph, **kwds)

    def vertex_property_names(self):
        """todo"""
        return self._vertex_property.iterkeys()
    vertex_property_names.__doc__ = IPropertyGraph.vertex_property_names.__doc__

    def vertex_properties(self):
        """todo"""
        return self._vertex_property
    # vertex_properties.__doc__ = IPropertyGraph.vertex_properties.__doc__

    def vertex_property(self, property_name, vids = None):
        """todo"""
        try:
            if vids is not None:
                return dict([(k,v) for k,v in self._vertex_property[property_name].iteritems() if k in vids])
            else:
                return self._vertex_property[property_name]
        except KeyError:
            raise PropertyError("property %s is undefined on vertices"
                                % property_name)
    vertex_property.__doc__=IPropertyGraph.vertex_property.__doc__

    def edge_property_names(self):
        """todo"""
        return self._edge_property.iterkeys()
    edge_property_names.__doc__ = IPropertyGraph.edge_property_names.__doc__

    def edge_properties(self):
        """todo"""
        return self._edge_property
    #  edge_properties.__doc__ = IPropertyGraph. edge_properties.__doc__

    def edge_property(self, property_name):
        """todo"""
        try:
            return self._edge_property[property_name]
        except KeyError:
            raise PropertyError("property %s is undefined on edges"
                                % property_name)
    edge_property.__doc__ = IPropertyGraph.edge_property.__doc__

    def graph_property(self, property_name):
        """todo"""
        try:
            return self._graph_property[property_name]
        except KeyError:
            raise PropertyError("property %s is undefined on graph"
                                % property_name)
    graph_property.__doc__ = IPropertyGraph.graph_property.__doc__
                                
    def graph_properties(self):
        return self._graph_property
    
    
    def graph_property_names(self):
        """todo"""
        return self._graph_property.iterkeys()

    def add_vertex_property(self, property_name, values = None):
        """todo"""
        if property_name in self._vertex_property:
            raise PropertyError("property %s is already defined on vertices"
                                % property_name)
        if values is None: values = {}                                
        self._vertex_property[property_name] = values
    add_vertex_property.__doc__ = IPropertyGraph.add_vertex_property.__doc__

    def extend_vertex_property(self, property_name, values ):
        """todo AND TO CHECK AND TEST !!"""
        if not isinstance(values, dict):
            raise TypeError("Values %s is not a type 'dict'" % values)                                
        if property_name not in self._vertex_property:
            print PropertyError("Property %s is not defined on vertices"
                                % property_name)
            print "Creating vertex property %s" % property_name
            self._vertex_property[property_name] = {}
        
        for k,v in values.iteritems():
            if k in self.vertices():
                if not self._vertex_property[property_name].has_key(k):
                    self._vertex_property[property_name][k] = v
                else:
                    print "Vertex id {} already has a value for vertex property {}".format(k, property_name)
            else:
                print "Vertex id {} doesn't exist in the graph !!".format(k)

    def remove_vertex_property(self, property_name):
        """todo"""
        try:
            del self._graph_property['units'][property_name]
        except:
            pass
        try:
            del self._vertex_property[property_name]
        except KeyError:
            raise PropertyError("property %s is undefined on vertices"
                                % property_name)
    remove_vertex_property.__doc__ = IPropertyGraph.remove_vertex_property.__doc__

    def add_edge_property(self, property_name, values =  None):
        """todo"""
        if property_name in self._edge_property:
            raise PropertyError("property %s is already defined on edges"
                                % property_name)
        if values is None: values = {}                                
        self._edge_property[property_name] = values
    add_edge_property.__doc__ = IPropertyGraph.add_edge_property.__doc__

    def remove_edge_property(self, property_name):
        """todo"""
        try:
            del self._graph_property['units'][property_name]
        except:
            pass
        try:
            del self._edge_property[property_name]
        except KeyError:
            raise PropertyError("property %s is undefined on edges"
                                % property_name)
    remove_edge_property.__doc__ = IPropertyGraph.remove_edge_property.__doc__

    def add_graph_property(self, property_name, values = None):
        """todo"""
        if property_name in self._graph_property:
            raise PropertyError("property %s is already defined on graph"
                                % property_name)
        if values is None: values = {}
        self._graph_property[property_name] = values
    
    def extend_graph_property(self, property_name, values):
        """todo"""
        assert values is not None
        if property_name not in self._graph_property:
            raise PropertyError("property %s is not defined on graph"
                                % property_name)

        if self.graph_property(property_name) is not None or self.graph_property(property_name) != []:
            assert isinstance(values, type(self.graph_property(property_name)))

        if isinstance(self.graph_property(property_name), list):
            self._graph_property[property_name].extend(values)
        elif isinstance(self.graph_property(property_name), dict):
            self._graph_property[property_name].update( dict([(k,v) for k,v in values.iteritems() if k not in self.graph_property(property_name).keys()]) )
        else:
            print "Unable to extend 'graph_property' with this type of data: {}".format(type(values))

    def remove_graph_property(self, property_name):
        """todo"""
        try:
            del self._graph_property[property_name]
            del self._graph_property['units'][property_name]
        except KeyError:
            raise PropertyError("property %s is undefined on graph"
                                % property_name)

    def remove_vertex(self, vid):
        """todo"""
        for prop in self._vertex_property.itervalues():
            prop.pop(vid, None)
        Graph.remove_vertex(self, vid)
    remove_vertex.__doc__ = Graph.remove_vertex.__doc__

    def clear(self):
        """todo"""
        for prop in self._vertex_property.itervalues():
            prop.clear()
        for prop in self._edge_property.itervalues():
            prop.clear()
        for prop in self._graph_property.itervalues():
            prop.clear()
        Graph.clear(self)
    clear.__doc__ = Graph.clear.__doc__

    def remove_edge(self, eid):
        """todo"""
        for prop in self._edge_property.itervalues():
            prop.pop(eid, None)
        Graph.remove_edge(self, eid)
    remove_edge.__doc__ = Graph.remove_edge.__doc__

    def clear_edges(self):
        """todo"""
        for prop in self._edge_property.itervalues():
            prop.clear()
        Graph.clear_edges(self)
    clear_edges.__doc__ = Graph.clear_edges.__doc__

    @staticmethod
    def _translate_property(values, trans_vid, trans_eid, key_translation = ValueType, value_translation = ValueType):
        # translation function
        from copy import deepcopy
        
        id_value = lambda value: value
        
        trans_vid = deepcopy(trans_vid)
        trans_vid[None] = None
        def translate_vid(vid): 
            if isinstance(vid, list): return [trans_vid[i]  for i in vid]
            if isinstance(vid, tuple): return tuple([trans_vid[i]  for i in vid])
            return trans_vid[vid]
        
        trans_eid = deepcopy(trans_eid)
        trans_eid[None] = None
        def translate_eid(eid):
            if isinstance(eid, list): return [trans_eid[i]  for i in eid]
            if isinstance(eid, tuple): return tuple([trans_eid[i]  for i in eid])
            return trans_eid[eid]

        translator = { ValueType : id_value, VertexIdType :  translate_vid, EdgeIdType : translate_eid }

        key_translator = translator[key_translation]
        value_translator = translator[value_translation]

        # translate vid and value
        return dict([(key_translator(vid),value_translator(val)) for vid, val in values.iteritems()])


    def _relabel_and_add_vertex_edge_properties(self,graph, trans_vid, trans_eid):
        
        # update properties on vertices
        for prop_name in graph.vertex_property_names():
            if prop_name not in self._vertex_property:
                self.add_vertex_property(prop_name)
            value_translator = graph.get_property_value_type(prop_name,VertexProperty)

            # import property into self. translate vid and value
            self.vertex_property(prop_name).update(self._translate_property(graph.vertex_property(prop_name), trans_vid, trans_eid, VertexIdType, value_translator))

        # update properties on edges
        for prop_name in graph.edge_property_names():
            if prop_name not in self._edge_property:
                self.add_edge_property(prop_name)
            
            # Check what type of translation is required for value of the property
            value_translator = graph.get_property_value_type(prop_name,EdgeProperty)
            
            # import property into self. translate vid and value
            self.edge_property(prop_name).update(self._translate_property(graph.edge_property(prop_name), trans_vid, trans_eid, EdgeIdType, value_translator))

    def translate_graph_property(self, prop_name, trans_vid, trans_eid):
        """ Translate a graph property according to meta info """
        old_prop = self.graph_property(prop_name)
        print type(old_prop)
        key_translator = self.get_graph_property_key_type(prop_name)
        value_translator = self.get_property_value_type(prop_name, GraphProperty)
        print 'translate_graph_property',prop_name, key_translator, value_translator
        
        return self._translate_property(old_prop, trans_vid, trans_eid, key_translator, value_translator)
    
    
    def extend(self, graph):
        """todo"""
        
        # add and translate the vertex and edge ids of the second graph
        trans_vid, trans_eid = Graph.extend(self,graph)
        
        # relabel the edge and vertex property
        self.__relabel_and_add_vertex_edge_properties(graph, trans_vid, trans_eid)
        
        # update properties on graph
        #gproperties = self.graph_property()
        newgproperties = {}
        for pname, prop in graph.graph_property_names():
            newgproperty = graph.translate_graph_property(pname,trans_vid, trans_eid)
            newgproperties[pname] = newgproperty

            prop.update(newgproperties)

        return trans_vid, trans_eid
        
    extend.__doc__ = Graph.extend.__doc__


    def set_graph_property_value_to_vid_type(self, propertyname, property_type = VertexProperty):
        """ Give meta info on property value type. Associate it to Vertex Id type """
        if not self._graph_property.has_key(self.metavidtypepropertyname):
            self.add_graph_property(self.metavidtypepropertyname,([],[],[],[]))
        prop = self.graph_property(self.metavidtypepropertyname)[property_type]
        prop.append(propertyname)

    def set_graph_property_value_to_eid_type(self, propertyname,  property_type = EdgeProperty):
        """ Give meta info on property value type. Associate it to Edge Id type """
        if not self._graph_property.has_key(self.metaeidtypepropertyname):
            self.add_graph_property(self.metaeidtypepropertyname,([],[],[],[]))
        prop = self.graph_property(self.metaeidtypepropertyname)[property_type]
        prop.append(propertyname)

    def set_graph_property_key_to_vid_type(self, propertyname):
        """ Give meta info on graph property key type. Associate it to Vertex Id type"""
        self.set_graph_property_value_to_vid_type(propertyname, 3)

    def set_graph_property_key_to_eid_type(self, propertyname):
        """ Give meta info on graph property key type. Associate it to Edge Id type """
        self.set_graph_property_value_to_eid_type(propertyname, 3)

    def get_property_value_type(self, propertyname, property_type = VertexProperty):
        """ Return meta info on property value type. """
        try:
            prop = self.graph_property(self.metavidtypepropertyname)[property_type]
            if propertyname in prop : return VertexIdType
        except:
            pass
        try:
            prop = self.graph_property(self.metaeidtypepropertyname)[property_type]
            if propertyname in prop : return EdgeIdType
        except:
            return ValueType

    def get_graph_property_key_type(self, propertyname):
        """ Return meta info on graph property key type. """
        return self.get_property_value_type(propertyname, 3)


    def __to_set(self, s):
        if not isinstance(s, set):
            if isinstance(s, list):
                s=set(s)
            else:
                s=set([s])
        return s

    def in_neighbors(self, vid, edge_type=None):
        """ Return the in vertices of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider (can be a set)

        :Returns:
        - `neighbors_list` : the set of parent vertices of the vertex vid
        """
        
        if vid not in self :
            raise InvalidVertex(vid)
        
        if edge_type==None:
            neighbors_list=set([self.source(eid) for eid in self._vertices[vid][0] ])
        else:
            edge_type=self.__to_set(edge_type) 
            edge_type_property = self._edge_property['edge_type']
            neighbors_list=set([self.source(eid) for eid in self._vertices[vid][0] if edge_type_property[eid] in edge_type])
        return neighbors_list

    def iter_in_neighbors(self, vid, edge_type=None):
        """ Return the in vertices of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider (can be a set)

        :Returns:
        - `iterator` : an iterator on the set of parent vertices of the vertex vid
        """
        return iter(self.in_neighbors(vid, edge_type))

    def out_neighbors(self, vid, edge_type=None):
        """ Return the out vertices of the vertex vid

        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider (can be a set)

        :Returns:
        - `neighbors_list` : the set of child vertices of the vertex vid
        """
        if vid not in self :
            raise InvalidVertex(vid)

        if edge_type==None:
            neighbors_list=set([self.target(eid) for eid in self._vertices[vid][1] ])
        else:
            edge_type=self.__to_set(edge_type) 
            edge_type_property = self._edge_property['edge_type']
            neighbors_list=set([self.target(eid) for eid in self._vertices[vid][1] if edge_type_property[eid] in edge_type])
        return neighbors_list


    def iter_out_neighbors(self, vid, edge_type=None):
        """ Return the out vertices of the vertex vid

        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider (can be a set)

        :Returns:
        - `iterator` : an iterator on the set of child vertices of the vertex vid
        """
        return iter(self.out_neighbors(vid, edge_type))

    def neighbors(self, vid, edge_type=None):
        """ Return the neighbors vertices of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider (can be a set)

        :Returns:
        - `neighbors_list` : the set of neighobrs vertices of the vertex vid
        """
        return self.in_neighbors(vid, edge_type) | self.out_neighbors(vid, edge_type)

    def iter_neighbors(self, vid, edge_type=None):
        """ Return the neighbors vertices of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider (can be a set)

        :Returns:
        - `iterartor` : iterator on the set of neighobrs vertices of the vertex vid
        """
        return iter(self.neighbors(vid, edge_type))
  

    def in_edges(self, vid, edge_type=None):
        """ Return in edges of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider (can be a set)

        :Returns:
        - `edge_list` : the set of the in edges of the vertex vid
        """
        if vid not in self :
            raise InvalidVertex(vid)

        if not edge_type:
            edge_list=set([eid for eid in self._vertices[vid][0]])
        else:
            edge_type=self.__to_set(edge_type)
            edge_type_property = self._edge_property['edge_type']
            edge_list=set([eid for eid in self._vertices[vid][0] if edge_type_property[eid] in edge_type])
        return  edge_list
        
    def iter_in_edges(self, vid, edge_type=None):
        """ Return in edges of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider (can be a set)

        :Returns:
        - `iterator` : an iterator on the set of the in edges of the vertex vid
        """  
        return iter(self.in_edges(vid, edge_type))


    def out_edges(self, vid, edge_type=None):
        """ Return out edges of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider (can be a set)

        :Returns:
        - `edge_list` : the set of the out edges of the vertex vid
        """
        if vid not in self :
            raise InvalidVertex(vid)
        
        if edge_type==None:
            edge_list=set([eid for eid in self._vertices[vid][1]])
        else:
            edge_type=self.__to_set(edge_type)
            edge_type_property = self._edge_property['edge_type']
            edge_list=set([eid for eid in self._vertices[vid][1] if edge_type_property[eid] in edge_type])
        return  edge_list

    def iter_out_edges(self, vid, edge_type=None):
        """ Return in edges of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider

        :Returns:
        - `iterator` : an iterator on the set of the in edges of the vertex vid
        """  
        return iter(self.out_edges(vid, edge_type))


    def edges(self, vid=None, edge_type=None):
        """ Return edges of the vertex vid
        If vid=None, return all edges of the graph
        
        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider

        :Returns:
        - `edge_list` : the set of the edges of the vertex vid
        """
        if vid==None:
            return set(self._edges.keys())       
        return self.out_edges(vid, edge_type) | self.in_edges(vid, edge_type)

    def iter_edges(self, vid, edge_type=None):
        """ Return in edges of the vertex vid
        If vid=None, return all edges of the graph
        
        :Parameters:
        - `vid` : a vertex id
        - `edges_type` : type of edges we want to consider

        :Returns:
        - `iterator` : an iterator on the set of the edges of the vertex vid
        """  
        return iter(self.edges(vid, edge_type))

    def neighborhood(self, vid, max_distance=1, edge_type=None):
        """ Return the neighborhood of the vertex vid at distance max_distance (the disc, not the circle)
        
        :Parameters:
        - `vid` : vertex id

        :Returns:
        - `neighbors_list` : the set of the vertices at distance below max_distance of the vertex vid (including vid)
        """
        dist=self.topological_distance(vid, edge_type=edge_type, max_depth=max_distance, full_dict=False)
        return set(dist.keys())

    def iter_neighborhood(self, vid, n, edge_type=None):
        """ Return the neighborhood of the vertex vid at distance n (the disc, not the circle)
        
        :Parameters:
        - `vids` : a set of vertex id

        :Returns:
        - `iterator` : an iterator on the set of the vertices at distance n of the vertex vid
        """
        return iter(self.neighborhood(vid, n, edge_type))    


    def topological_distance(self, vid, edge_type = None, edge_dist = lambda x,y : 1, max_depth=float('inf'), full_dict=True, return_inf = True):
        """ Return the distances of each vertices from the vertex `vid` according a cost function
        
        :Parameters:
        - `vid` (int) - a vertex id
        - `edges_type` (str) - type of edges we want to consider e.g. 's' or 't'
        - `edge_dist` (function) - the cost function
        - `max_depth` (float) - the maximum depth that we want to reach
        - `full_dict` (bool) - if True this function will return the entire dictionary (with inf values)
        - `return_inf` (bool) - if True (default) return 'inf' values, else 'nan'.

        :Returns:
        - `dist_dict` : a dictionary of the distances, key : vid, value : distance
        """
        dist={}
        reduced_dist={}
        reduced_dist[vid]=0
        Q=[]

        infinite_distance = float('inf')
        for k in self._vertices.iterkeys():
            dist[k] = infinite_distance
            heappush(Q, (dist[k], k))

        dist[vid]=0
        heappush(Q, (dist[vid], vid))
        treated=set()
        modif=True
        n = self.nb_vertices()
        while (len(treated)!=n and modif):
            modif = False
            actualVid = heappop(Q)
            while actualVid[1] in treated and actualVid[0] == infinite_distance:
                actualVid = heappop(Q)

            if actualVid[0] != infinite_distance:
                actualVid = actualVid[1]
                treated.add(actualVid)
            
                for neighb in self.iter_neighbors(actualVid, edge_type):
                    dist_neighb = dist[neighb]
                    if (((dist_neighb == infinite_distance) or (dist_neighb > dist[actualVid] + edge_dist(neighb, actualVid)))
                          and (dist[actualVid] + edge_dist(neighb, actualVid) < max_depth+1 ) ):
                        
                        dist_neighb = dist[actualVid] + edge_dist(neighb, actualVid)
                        dist[neighb] = dist_neighb
                        reduced_dist[neighb] = dist_neighb
                        heappush(Q, (dist_neighb, neighb))
                    modif = True
        #~ return (reduced_dist, dist)[full_dict], Q
        if not return_inf:
            dist = dict( [(k,(v if v != infinite_distance else np.nan)) for k,v in dist.iteritems()] )

        return (reduced_dist, dist)[full_dict]


    def adjacency_matrix(self, edge_type = None, edge_dist = 1, no_edge_val = 0, oriented = True, reflexive = True, reflexive_value = 0):
        """
        Return the adjacency matrix of the graph.
        :Parameters:
        - `edge_type` : type of edges we want to consider
        - `edge_dist` : cost ot cost function to apply between two edges, default : 1
        - `no_edge_val` : cost to put if there is no edge between two vertices, default : 0
        - `oriented` : if True, the graph is considered oriented and we always add an edge j -> i if i -> j exists
        - `reflexive` : if True, the graph is considered reflexive and we will put the cost or the cost_function `reflexive_value` on the diagonal of the adjacency_matrix, default : 0

        :Return:
        - `numpy.array` : a NxN matrix
        """
        if not isinstance(edge_dist, type(lambda m: 1)):
            val_edge_dist = edge_dist
            edge_dist = lambda g, x, y : val_edge_dist
        if not isinstance(reflexive_value, type(lambda m: 1)):
            val_reflexive_value = reflexive_value
            reflexive_value = lambda g, x, y : val_reflexive_value
        
        n = self.nb_vertices()
        adjacency_matrix = np.array(n*[n*[no_edge_val]])
        for edge in self.edges(edge_type=edge_type):
            v1, v2 = self.edge_vertices(edge)
            adjacency_matrix[v1, v2] = edge_dist(self, v1, v2)
            if not oriented:
                adjacency_matrix[v2, v1] = edge_dist(self, v2, v1)
        if reflexive:
            for i in range(n) : adjacency_matrix[i, i] = reflexive_value(self, i, i)
        return adjacency_matrix

    def floyd_warshall(self, edge_type = None, edge_dist = 1, oriented = False):
        adjacency_matrix = self.adjacency_matrix(edge_type, edge_dist, float('inf'), oriented)
        n = self.nb_vertices()
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    adjacency_matrix[i, j]=min(adjacency_matrix[i, j], 
                                               adjacency_matrix[i, k] + adjacency_matrix[k, j])
        return adjacency_matrix
        

    def _add_vertex_to_region(self, vids, region_name):
        """
        add a set of vertices to a region
        """
        if not "regions" in self._vertex_property:
            warnings.warn("Property 'regions' is not defined on vertex. Adding it!")
            self._vertex_property["regions"] = {}

        for vid in vids:
            if self._vertex_property["regions"].has_key(vid):
                self._vertex_property["regions"][vid].append(region_name)
            else:
                self._vertex_property["regions"][vid]=[region_name]

            self._graph_property[region_name].append(vid)

    def _remove_vertex_from_region(self, vids, region_name):
        """
        remove a set of vertices to a region
        """
        for vid in vids:
            self._vertex_property["regions"][vid].remove(region_name)
            if self._vertex_property["regions"][vid]==[]:
                self._vertex_property["regions"].pop(vid)
                
            self._graph_property[region_name].remove(vid)

    def add_vertex_to_region(self, vids, region_name):
        """
        add a set of vertices to a region
        """
        if not region_name in self._graph_property:
            #~ raise PropertyError("property %s is not defined on graph" % region_name)
            warnings.warn("Property %s is not defined for vertices on the graph, adding it..." % region_name)
            self._graph_property[region_name] = vids
        
        self._add_vertex_to_region(self.__to_set(vids), region_name)

    def remove_vertex_from_region(self, vids, region_name):
        """
        remove a set of vertices to a region
        """
        if not region_name in self._graph_property:
            raise PropertyError("property %s is not defined on graph"
                                % region_name)
        self._remove_vertex_from_region(self.__to_set(vids), region_name)

    def add_region_from_func(self, func, region_name):
        """ Create a region of vertices according a function
        
        :Parameters:
        - `func` : the function to make the region (might return True or False)
        - `region_name` : the name of the region
        
        """
        if region_name in self._graph_property:
            raise PropertyError("property %s is already defined on graph"
                                % region_name)
        self._graph_property[region_name]=[]
        if not "regions" in self._vertex_property.keys():
            self.add_vertex_property("regions")
        for vid in self._vertices.keys():
            if func(self, vid):
                self._add_vertex_to_region(set([vid]), region_name)

    def add_regions_from_dict(self, dict_regions, region_names):
        """
        If one already posses a dict indicating for a list of vertex which region they belong to, it can be given to the graph directly.
        
        :Parameters:
        - `dict_regions` (dict) - *keys = ids (SpatialImage); *values = intergers indicating the region(s)
        - `region_name` (list) - a list containing the name of the region(s)
        """
        list_regions = np.unique(dict_regions.values())
        if len(region_names) != len(list_regions):
            warnings.warn("You didn't provided the same number of regions and region names.")
            pass
        
        if not "regions" in self._vertex_property.keys():
            self.add_vertex_property("regions")
        
        for region, region_name in enumerate(region_names):
            if region_name in self._graph_property:
                raise PropertyError("property %s is already defined on graph"
                                    % region_name)
            self._graph_property[region_name]=[]
            for vid in dict_regions:
                if dict_regions[vid] == list_regions[region]:
                    self._add_vertex_to_region(set([vid]), region_name)

    def iter_region(self, region_name):
        if not region_name in self._graph_property:
            raise PropertyError("property %s is not defined on graph"
                                % region_name)
        return iter(self._graph_property[region_name])

    def remove_region(self, region_name):
        """ Remove a region 
        
        :Parameters:
        - `region_name` : the name of the region
        
        """
        if not region_name in self._graph_property:
            raise PropertyError("property %s is not defined on graph"
                                % region_name)

        for vid in self.iter_region(region_name):
            self._vertex_property["regions"][vid].remove(region_name)
            if self._vertex_property["regions"][vid]==[]:
                self._vertex_property["regions"].pop(vid)

        return self._graph_property.pop(region_name)

    def is_connected_region(self, region_name, edge_type=None):
        """
        Return True if a region is connected
        """
        if not region_name in self._graph_property:
            raise PropertyError("property %s is not defined on graph"
                                % region_name)
        region_sub_graph=Graph.sub_graph(self, self._graph_property[region_name])
        distances=region_sub_graph.topological_distance(region_sub_graph._vertices.keys()[0], edge_type=edge_type)
        return not float('inf') in distances.values()


    def to_networkx(self):
        import networkx as nx
        """ Return a NetworkX Graph from a graph.

        :Parameters: 
         - `g` - TemporalPropertyGraph (property graphs temporaly linked)

        :Returns: 
         - A NetworkX graph.
        """

        g = self

        graph = nx.Graph()
        graph.add_nodes_from(g.vertices())
        graph.add_edges_from(( (g.source(eid), g.target(eid)) for eid in g.edges()))

        # Add graph, vertex and edge properties
        for k, v in g.graph_properties().iteritems():
            graph.graph[k] = v

        vp = g._vertex_property
        for prop in vp:
            for vid, value in vp[prop].iteritems():
                graph.node[vid][prop] = value
        
        ep = g._edge_property
        for eid in g.edges():
            graph.edge[g.source(eid)][g.target(eid)]['eid'] = eid

        for prop in ep:
            for eid, value in ep[prop].iteritems():
                graph.edge[g.source(eid)][g.target(eid)][prop] = value

        return graph 

    def from_networkx(self, graph):
        """ Return a Graph from a NetworkX Directed graph.
        :Parameters: 
            - `graph` : A NetworkX graph.

        :Returns: 
            - `g`: a :class:`~sa_oa.container.interface.Graph`.
        """
        self.clear()
        g = self

        if not graph.is_directed():
            graph = graph.to_directed()

        vp = self._vertex_property

        for vid in graph.nodes_iter():
            g.add_vertex(vid)
            d = graph.node[vid]
            for k, v in d.iteritems():
                vp.setdefault(k,{})[vid] = v

        ep = self._edge_property
        for source, target in graph.edges_iter():
            d = graph[source][target]
            eid = d.get('eid')
            eid = g.add_edge(source, target, eid)
            for k, v in d.iteritems():
                if k != 'eid':
                    ep.setdefault(k,{})[eid] = v

        gp = self._graph_property
        gp.update(graph.graph)

        return g

