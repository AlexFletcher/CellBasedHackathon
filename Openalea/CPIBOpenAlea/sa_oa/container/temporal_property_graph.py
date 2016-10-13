# -*- python -*-
#
#       OpenAlea.Container
#
#       Copyright 2011 - 2013 INRIA - CIRAD - INRA
#
#       File author(s): Jonathan Legrand <jonathan.legrand@ens-lyon.fr>
#                       Christophe Pradal
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://sa_oa.gforge.inria.fr
#
################################################################################
"""This module provide a class that extends the PropertyGraph with different types of edges"""

__license__ = "Cecill-C"
__revision__ = " $Id: temporal_property_graph.py 17387 2014-08-31 18:28:14Z jlegra02 $ "

from property_graph import *
import warnings
import numpy as np

from sa_oa.container.temporal_graph_analysis import translate_ids_Image2Graph, translate_keys_Graph2Image, exist_all_relative_at_rank

class TemporalPropertyGraph(PropertyGraph):
    """
    Simple implementation of PropertyGraph using
    dict as properties and two dictionaries to
    maintain these properties
    """
    STRUCTURAL = 's'
    TEMPORAL = 't'

    def __init__(self, graph=None, **kwds):
        PropertyGraph.__init__(self, graph, idgenerator='max',**kwds)
        self.add_edge_property('edge_type')

        # list of dict
        # each dict define the mapping between the new and the old vid.
        # old label define both graph index and local id.
        self.add_edge_property('old_label')
        self.add_vertex_property('old_label')
        self.add_vertex_property('index')
        self._old_to_new_ids = []
        self.nb_time_points = 0


    def extend(self, graphs, mappings, time_steps = None):
        """ Extend the structure with graphs and mappings.
        Each graph contains structural edges. 
        Mapping define dynamic edges between two graphs.

        :Parameters:
            - graphs - a list of PropertyGraph
            - mappings - a list defining the dynamic or temporal edges between two graphs.
            - time_steps - time corresponding to each graph

        :warning:: len(graphs) == len(mappings)-1
        """

        assert len(graphs) == len(mappings)+1

        self.append(graphs[0])
        #~ self.add_graph_property('nb_time_points')
        #~ len(mappings) = self.graph_property('nb_time_points')
        for g, m in zip(graphs[1:],mappings):
            self.append(g,m)
            self.nb_time_points += 1

        if time_steps is not None:
            assert len(graphs) == len(time_steps)
            self.add_graph_property('time_steps',time_steps)

        return self._old_to_new_ids


    def append(self, graph, mapping=None):
        """

        """
        if mapping:
            assert len(self._old_to_new_ids) >= 1, 'You have to create temporal edges between two graphs. Add a graph first without mapping'

        current_index = len(self._old_to_new_ids)

        edge_types = self.edge_property('edge_type')
        old_edge_labels = self.edge_property('old_label')
        old_vertex_labels = self.vertex_property('old_label')
        indices = self.vertex_property('index')

        # add and translate the vertex and edge ids of the second graph
        relabel_ids = Graph.extend(self,graph)
        old_to_new_vids, old_to_new_eids = relabel_ids
        # relabel the edge and vertex property
        self._relabel_and_add_vertex_edge_properties(graph, old_to_new_vids, old_to_new_eids)

        # update properties on graph
        temporalgproperties = self.graph_properties()

        # while on a property graph, graph_property are just dict of dict, 
        # on a temporal property graph, graph_property are dict of list of dict
        # to keep the different values for each time point.

        for gname in graph.graph_property_names():
            if gname in [self.metavidtypepropertyname,self.metavidtypepropertyname]:
                temporalgproperties[gname] = graph.graph_property(gname)
            else:
                newgproperty = graph.translate_graph_property(gname, old_to_new_vids, old_to_new_eids)
                temporalgproperties[gname] = temporalgproperties.get(gname,[])+[newgproperty]

        self._old_to_new_ids.append(relabel_ids)

        # set edge_type property for structural edges
        for old_eid, eid in old_to_new_eids.iteritems():
            edge_types[eid] = self.STRUCTURAL
            old_edge_labels[eid] = old_eid

        for old_vid, vid in old_to_new_vids.iteritems():
            old_vertex_labels[vid] = old_vid
            indices[vid] = current_index

        def use_sub_lineage(mother, daughters, on_ids_source, on_ids_target):
            found_sub_lineage=False; tmp_daughters = []
            for d in daughters:
                if iterable(d):
                    found_sub_lineage=True; tmp_d = []
                    if "sub_lineage" not in self.graph_properties():
                        self.add_graph_property("sub_lineage")
                    for sub_d in d:
                        if iterable(sub_d):
                            use_sub_lineage(mother, sub_d, on_ids_source, on_ids_target)
                        else:
                            tmp_d.append(on_ids_target[0][sub_d])
                    tmp_daughters.append(tmp_d)
                else:
                    tmp_daughters.append(on_ids_target[0][d])
            if found_sub_lineage:
                self.graph_property("sub_lineage").update({on_ids_source[0][mother]:tmp_daughters})

        def flatten(l):
            import collections
            for el in l:
                if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
                    for sub in flatten(el):
                        yield sub
                else:
                    yield el

        if mapping:
            unused_lineage = {}
            on_ids_source, on_ids_target = self._old_to_new_ids[-2:]
            for k, l in mapping.iteritems():
                l_flat = list(flatten(l)) # flatten the lineage if not
                # Check if the mother cell and ALL daugthers are present in their respective graph : WE DON'T WANT TO CREATE A PARTIAL LINEAGE !!!!
                if on_ids_source[0].has_key(k) and ( sum([on_ids_target[0].has_key(v) for v in l_flat]) == len(l_flat) ):
                    use_sub_lineage(k, l, on_ids_source, on_ids_target)
                    for v in l_flat:
                        eid = self.add_edge(on_ids_source[0][k], on_ids_target[0][v])
                        edge_types[eid] = self.TEMPORAL
                else:
                    unused_lineage.update({k:l})
            if unused_lineage != {}:
                print "Un-used lineage info to avoid partial lineage, t{} to t{}: {}".format(current_index-1,current_index,unused_lineage)

        return relabel_ids

    def clear(self):
        PropertyGraph.clear(self)
        self._old_to_new_ids = []

    def __to_set(self, s):
        if not isinstance(s, set):
            if isinstance(s, list):
                s=set(s)
            else:
                s=set([s])
        return s

    def children(self, vid):
        """ Return children of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id

        :Returns:
        - `children_list` : the set of the children of the vertex vid
        """
        return self.out_neighbors(vid, 't')

    def iter_children(self, vid):
        """ Return children of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id

        :Returns:
        - `iterator` : an iterator on the set of the children of the vertex vid
        """
        return iter(self.children(vid))

    def has_children(self, vid):
        """
        Return True if the vid `vid` has a child or children.
        """
        if self.children(vid) != set():
            return True
        else:
            return False

    def parent(self, vid):
        """ Return parents of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id

        :Returns:
        - `parents_list` : the set of the parents of the vertex vid
        """
        return self.in_neighbors(vid, 't')

    def iter_parent(self, vid):
        """ Return parents of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id

        :Returns:
        - `iterator` : the set of the children of the vertex vid
        """
        return iter(self.parent(vid))

    def has_parent(self, vid):
        """
        Return True if the vid `vid` has a parent.
        """
        if self.parent(vid) != set():
            return True
        else:
            return False

    def sibling(self, vid):
        """ Return sibling of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id

        :Returns:
        - `sibling_list` : the set of sibling of the vertex vid
        """
        if self.parent(vid):
            return self.children(self.parent(vid).pop())-set([vid])
        else:
            return None

    def iter_sibling(self, vid):
        """ Return of the vertex vid
        
        :Parameters:
        - `vid` : a vertex id

        :Returns:
        - `iterator` : an iterator on the set of sibling of the vertex vid
        """
        return iter(self.sibling(vid))    

    def descendants(self, vids, n = None):
        """ Return the 0, 1, ..., nth descendants of the vertex vid
        
        :Parameters:
        - `vids` : a set of vertex id

        :Returns:
        - `descendant_list` : the set of the 0, 1, ..., nth descendant of the vertex vid
        """
        edge_type='t'
        neighbs=set()
        vids=self.__to_set(vids)
        if n==0 :
            return vids
        elif n==1 :
            for vid in vids:
                neighbs |= (self.out_neighbors(vid, edge_type) | set([vid]))
            return neighbs
        else :
            if n is None :
                n = self.nb_time_points                
            for vid in vids :
                neighbs |= (self.descendants(self.out_neighbors(vid, edge_type), n-1) | set([vid]))
                if list(neighbs)==self._vertices.keys():
                    return neighbs
        return neighbs

    def iter_descendant(self, vids, n = None):
        """ Return the 0, 1, ..., nth descendants of the vertex vid
        
        :Parameters:
        - `vids` : a set of vertex id

        :Returns:
        - `iterator` : an iterator on the set of the 0, 1, ..., nth descendants of the vertex vid
        """
        return iter(self.descendant(vids, n))

    def ancestors(self, vids, n = None):
        """Return the 0, 1, ..., nth ancestors of the vertex vid
        
        :Parameters:
        - `vids` : a set of vertex id
        
        :Returns:
        - `anestors_list` : the set of the 0, 1, ..., nth ancestors of the vertex vid
        """
        edge_type='t'
        neighbs=set()
        vids=self.__to_set(vids)
        if n==0 :
            return vids
        elif n==1 :
            for vid in vids:
                neighbs |= (self.in_neighbors(vid, edge_type) | set([vid]))
            return neighbs
        else :
            if n is None:
                n = self.nb_time_points
            for vid in vids :
                neighbs |= (self.ancestors(self.in_neighbors(vid, edge_type), n-1) | set([vid]))
                if list(neighbs)==self._vertices.keys():
                    return neighbs
        return neighbs

    def iter_ancestors(self, vids, n):
        """ Return the 0, 1, ..., nth ancestors of the vertex vid
        
        :Parameters:
        - `vids` : a set of vertex id

        :Returns:
        - `iterator` : an iterator on the set of the 0, 1, ..., nth ancestors of the vertex vid
        """
        return iter(self.ancestors(vids, n))

    def lineaged_vertex(self, fully_lineaged=True):
        """
        Return ids of lineaged vertices.
        One can ask for strict lineage, i.e. only vertices temporally linked from the beginning (`self.vertex_property('index')`==`0`) to the end (`self.nb_time_points`).
        
        :Parameter:
         - `fully_lineaged` (bool) : if True, return vertices temporally linked from the beginning to the end, otherwise return vertices having at least a parent or a child(ren).
        """
        if fully_lineaged:
            last_tp_ids_lineaged_from_0 = [k for k in self.vertices() if exist_all_relative_at_rank(self, k, -self.nb_time_points)]
            first_tp_ids_lineaged_from_0 = [k for k in self.ancestors(last_tp_ids_lineaged_from_0, self.nb_time_points) if self.vertex_property('index')[k]==0]
            first_tp_ids_fully_lineaged = [k for k in first_tp_ids_lineaged_from_0 if exist_all_relative_at_rank(self, k, self.nb_time_points)]
            return list(np.sort(list(self.descendants(first_tp_ids_fully_lineaged, self.nb_time_points))))
        else:
            return [k for k in self.vertices() if (self.has_children(k) or self.has_parent(k))]

    def vertex_at_time(self, time_point, lineaged=False, fully_lineaged=False, as_parent=False, as_children=False):
        """
        Return vertices ids corresponding to a given time point in the TPG.
        Optional parameters can be used to filter the list of vertices ids returned.
        
        :Parameters:
         - `time_point` (int) : integer giving which time point should be considered
         - `lineaged` (bool) : if True, return vertices having at least a parent or a child(ren)
         - `fully_lineaged` (bool) : if True, return vertices linked from the beginning to the end of the graph
         - `as_parent` (bool) : if True, return vertices lineaged as parents
         - `as_children` (bool) : if True, return vertices lineaged as children
        """
        if lineaged and (not as_parent and not as_children):
            as_parent = as_children = True
        if as_parent or as_children:
            lineaged = True

        if fully_lineaged:
            return [k for k in self.lineaged_vertex(fully_lineaged=True) if self.vertex_property('index')[k]==time_point]
        if lineaged:
            return [ k for k in self.vertices() if (self.vertex_property('index')[k]==time_point) and
             ( (self.has_children(k) if as_parent else False ) or ( self.has_parent(k) if as_children else False ) ) ]
        else:
            return [k for k in self.vertices() if self.vertex_property('index')[k]==time_point]

    def vertex_property_with_image_labels(self, vertex_property, graph_labels2keep = None, image_labels2keep = None, time_point = None):
        """
        Return a dictionary extracted from the graph.vertex_property(`vertex_property`) with relabelled keys into "images labels" thanks to the dictionary graph.vertex_property('old_labels').
        
        :Parameters:
         - `vertex_property` : can be an existing 'graph.vertex_property' or a string refering to an existing graph.vertex_property to extract.
         - `image_labels2keep` (int|list) : a list of "image type" labels (i.e. from SpatialImages) to return in the `relabelled_dictionnary`
         - `graph_labels2keep` (int|list) : a list of "graph type" ids (i.e. from PropertyGraphs) to return in the `relabelled_dictionnary`

        :Returns:
         - *key_ vertex/cell label, *values_ `vertex_property`
        
        :Examples:
         graph.vertex_property_with_image_labels( graph.vertex_property('volume') )
         graph.vertex_property_with_image_labels( 'volume' )
         graph.vertex_property_with_image_labels( 'volume' , SpatialImageAnalysis.L1() )
        
        """
        if isinstance(vertex_property,str):
            if (graph_labels2keep is not None):
                return translate_keys_Graph2Image(self, self.vertex_property(vertex_property,graph_labels2keep))
            else:
                vertex_property = self.vertex_property(vertex_property)

        if (graph_labels2keep is not None):
            return translate_keys_Graph2Image(self, dict( [(k,v) for k,v in vertex_property.iteritems() if k in graph_labels2keep] )) 

        if (image_labels2keep is not None):
            return translate_keys_Graph2Image(self, dict( [(k,v) for k,v in vertex_property.iteritems() if k in translate_ids_Image2Graph(self,image_labels2keep,time_point)] )) 


    def region_vids(self, region_name):
        """
        Return a list of vids (TPG vertex id type) that belong to the region `region_name` according to graph.
        :Parameters:
         - `region_name` (str) : the nama of a previously defined region via `self.add_vertex_to_region()`
        """
        if 'regions' in list(self.vertex_properties()):
            return sorted(list(set([k for k,v in self.vertex_property('regions').iteritems() for r in v if r==region_name])))
        else:
            print "No property 'regions' added to the graph yet!"
            return None


def iterable(obj):
    try :
        iter(obj)
        return True
    except TypeError,te:
        return False
