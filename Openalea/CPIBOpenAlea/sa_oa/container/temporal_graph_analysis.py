# -*- python -
#
#       OpenAlea.Container
#
#       Copyright 2012 INRIA - CIRAD - INRA
#
#       File author(s):  Jonathan Legrand <jonathan.legrand@ens-lyon.fr>
#                        Frederic Boudon <frederic.boudon@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://sa_oa.gforge.inria.fr
#
################################################################################
"""This module helps to analyse TemporalPropertyGraph from Spatial Images."""

import warnings, types, numpy as np, copy, math
import matplotlib.pyplot as plt
from numpy.linalg import svd
#from sa_oa.image.algo.analysis import return_list_of_vectors


def regions_from_lineage(graph, starting_tp=0):
    """
    Generate a dict where vertex ids (vids) -from a starting time point- will receive the id of their ancestor at that starting point.
    """
    rlineage = {}
    for vid in graph.lineaged_vertex(True):
        tp = graph.vertex_property('index')[vid]
        if tp >= starting_tp:
            ancestors = list(graph.ancestors(vid,n=tp))
            ancestor_starting_tp = [vid2 for vid2 in ancestors if graph.vertex_property('index')[vid2]==starting_tp]
            assert len(ancestor_starting_tp)==1
            rlineage[vid] = ancestor_starting_tp[0]

    return rlineage


def lineage_colors(graph, list_cell, alea_range=None, distance_between_colors=5):
    """
    Generate a dict where keys -from a given list `list_cell`- receive a random integer from the list as value.
    """
    import random
    import math
    if alea_range==None:
        alea_range = [0,255]
    if isinstance(alea_range,list) and (len(alea_range)==2):
        lineage_color = {}
        for n,k in enumerate(list_cell):
            print n,'/',len(list_cell)
            if n==0:
                lineage_color[k] = random.randint(alea_range[0], alea_range[1])
            else:
                next_color = random.randint(alea_range[0], alea_range[1])
                nei = graph.neighbors(k)
                colored_nei = list(nei & set(lineage_color.keys()))
                if colored_nei != []:
                    seq = list(set(np.arange(alea_range[1]-alea_range[0]))-set([lineage_color[c_n] for c_n in colored_nei])-set(np.arange(next_color,next_color+distance_between_colors))-set(np.arange(next_color-distance_between_colors,next_color)))
                    tmp_distance_between_colors = distance_between_colors
                    while seq == []:
                        tmp_distance_between_colors -= 1
                        if tmp_distance_between_colors <= 1:
                            warnings.warn('The range was not wide enought!')
                            return None
                        seq = list(set(np.arange(alea_range[0], alea_range[1]))-set([lineage_color[c_n] for c_n in colored_nei])-set(np.arange(next_color,next_color+tmp_distance_between_colors))-set(np.arange(next_color-tmp_distance_between_colors,next_color)))
                    lineage_color[k] = random.choice(seq)
                else:
                    lineage_color[k] = next_color

        return lineage_color

def add_graph_vertex_property_from_dictionary(graph, name, dictionary, unit=None):
    """
    Add a vertex property with name 'name' to the graph build from an image.
    The values of the property are given as by a dictionary where keys are TemporalPropertyGraph vertex labels.
    """
    if name in graph.vertex_properties():
        if (unit is not None) and (not graph._graph_property["units"].has_key(name) or graph._graph_property["units"](name) is None):
            graph._graph_property["units"].update({name:unit})
            print '{} unit upgraded in graph to {}.'.format(name, unit)
        raise ValueError('Existing vertex property {}'.format(name))

    graph.add_vertex_property(name)
    graph.vertex_property(name).update( dictionary )
    if unit is not None:
        graph._graph_property["units"].update({name:unit})


def keys_to_mother_id(graph, data, rank = 1):
    """
    Translate a dict with daughter ids as key, to a dict with mother ids as keys.
    :Parameters:
     - 'graph' (TPG): the TPG to be used for translation
     - 'data' (dict): the dictionary to translate
    """
    translated_dict = {}
    for vid in data.keys():
        if exist_relative_at_rank(graph,vid,-rank):
            mother = list(graph.ancestors(vid,rank)-graph.ancestors(vid,rank-1))[0]
            if not translated_dict.has_key(mother):
                translated_dict[mother] = data[vid]

    return translated_dict


def keys_to_daughter_id(graph, data, rank = 1):
    """
    Translate a dict with mother ids as key, to a dict with daughter ids as keys.
    :Parameters:
     - 'graph' (TPG): the TPG to be used for translation
     - 'data' (dict): the dictionary to translate
    """
    translated_dict = {}
    for vid in data.keys():
        if exist_relative_at_rank(graph,vid,rank):
            daughter = list(graph.descendants(vid,rank)-graph.descendants(vid,rank-1))
            for d in daughter:
                translated_dict[d] = data[vid]

    return translated_dict


def translate_ids_Graph2Image(graph, id_list):
    """
    Return a list which contains SpatialImage ids type translated from the TPG ids type `id_list`.

    :Parameters:
     - `graph` (TPG): the TemporalPropertyGraph containing the translation informations
     - `id_list` (list) - graphs ids type
    """
    if isinstance(id_list,int):
        return graph.vertex_property('old_label')[id_list]

    if isinstance(id_list,set):
        id_list = list(id_list)
    if not isinstance(id_list,list):
        raise ValueError('This is not an "int" or a "list" type variable.')

    return [graph.vertex_property('old_label')[k] for k in id_list]

def translate_ids_Image2Graph(graph, id_list, time_point):
    """
    Return a list which contains TPG ids type translated from the SpatialImage ids type `id_list`.

    :Parameters:
     - `graph` (TPG): the TemporalPropertyGraph containing the translation informations
     - `id_list` (list) - SpatialImage ids type
     - `time_point` (int) - index of the SpatialImage in the TemporalPropertyGraph

    :WARNING:
        `time_point` numbers starts at '0'
    """
    if isinstance(id_list,set):
        id_list = list(id_list)
    if (not isinstance(id_list,list)) and (not isinstance(id_list,int)):
        raise ValueError('This is not an "int" or a "list" type variable.')

    graph_labels_at_time_point = graph.vertex_at_time(time_point)
    Image2Graph_labels_at_time_point = dict( (v,k) for k,v in graph.vertex_property('old_label').iteritems() if k in graph_labels_at_time_point )

    if isinstance(id_list,int):
        if id_list in Image2Graph_labels_at_time_point:
            return Image2Graph_labels_at_time_point[id_list]
        else:
            print 'Label {0} is not in the graph at t{1}'.format(id_list,time_point+1)
            return []

    Image2Graph_labels, Image_labels_not_found = [], []
    for k in id_list:
        if k in Image2Graph_labels_at_time_point.keys():
            Image2Graph_labels.append(Image2Graph_labels_at_time_point[k])
        else:
            Image_labels_not_found.append(k)

    if Image_labels_not_found != []:
        warnings.warn("The cell ids"+str(Image_labels_not_found)+"were not in the graph!")

    return Image2Graph_labels

def translate_keys_Graph2Image(graph, dictionary, time_point=None):
    """
    Return a dictionary which keys are SpatialImage ids type .
    Initial keys are graph ids type and need to be translated into SpatialImage ids type.

    :Parameters:
    - `dictionary` (dict) - keys are SpatialImage ids type;
    - `time_point` (int) - index of the SpatialImage in the TemporalPropertyGraph;

    :WARNING:
        `time_point` numbers starts at '0'
    """
    if not isinstance(dictionary,dict):
        raise ValueError('This is not a "dict" type variable.')

    translated_dict = {}
    if time_point is None:
        for k in dictionary:
            if graph.vertex_property('old_label').has_key(k):
                translated_dict[graph.vertex_property('old_label')[k]] = dictionary[k]
            else:
                raise KeyError("The dictionary you want to translate contain label from more than one time point, found redundant keys!")

        return translated_dict
    else:
        return dict( (graph.vertex_property('old_label')[k], v) for k, v in dictionary.iteritems() if graph.vertex_property('index')[k] == time_point )

def translate_keys_Image2Graph(graph, dictionary, time_point):
    """
    Return a dictionary which keys are graph ids type .
    Initial keys are SpatialImage ids type and need to be translated into graph ids type.

    :Parameters:
    - `dictionary` (dict) - keys are graph ids type;
    - `time_point` (int) - index of the SpatialImage in the TemporalPropertyGraph;

    :WARNING:
        `time_point` numbers starts at '0'
    """
    if not isinstance(dictionary,dict):
        raise ValueError('This is not a "dict" type variable.')

    return dict( (k,dictionary[v]) for k,v in graph.vertex_property('old_label').iteritems() if (graph.vertex_property('index')[k] == time_point) and (v in dictionary) )


def __normalized_parameters(func):
    def wrapped_function(graph, vertex_property, vids = None, rank = 1 , verbose = False):
        """
        :Parameters:
        - 'graph' : a TPG.
        - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
        - 'vids' : by default a vertex id or a list of vertex ids. If 'vids=None' the mean absolute deviation will be computed for all ids present in the graph provided.
        - 'rank' : neighborhood at distance 'rank' will be used.

        :Return:
        - a single value if vids is an interger, or a dictionnary of *keys=vids and *values= "result of applyed fucntion `func`"
        """
        # -- If a name is given, we use vertex_property stored in the graph with this name.
        if isinstance(vertex_property,str):
            vertex_property = graph.vertex_property(vertex_property)
        # -- If an instancemethod is given, we use create a dictionary for the vids base ont the method.
        if isinstance(vertex_property,types.MethodType):
            vertex_property = dict([(vid,vertex_property(vid)) for vid in graph.vertices()])

        # -- If no vids provided we compute the function for all keys present in the vertex_property
        if vids==None:
            vids = vertex_property.keys()

        # -- Now execute the called 'func':
        if isinstance(vids, int):
            # - for single id, compute single result
            return func(graph, vertex_property, vids, rank)
        else:
            # - for set of ids, we compute a dictionary of resulting values.
            l={}
            for k in vids:
                if verbose and k%10==0: print k,'/',len(vids)
                try:
                    l[k] = func(graph, vertex_property, k, rank, edge_type='s')
                except:
                    print 'Error computing {} value for vid {}...'.format(func, k)
            return l

    return  wrapped_function


@__normalized_parameters
def laplacian(graph, vertex_property, vid, rank, edge_type):
    """
    Sub-function computing the laplacian between ONE vertex ('vid') and its neighbors at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' : a vertex id.
    - 'rank' : neighborhood at distance 'rank' will be used.

    :Return:
    - a single value = laplacian between vertex 'vid' and its neighbors at rank 'rank'.
    """
    if rank == 1:
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)
        vid_neighborhood.remove(vid)
    else: # if rank > 1, we want to compute the change only over the cell at `rank` and not for all cells between rank 1 and `rank`.
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)-graph.neighborhood(vid,rank-1, edge_type)

    nb_neighborhood = len(vid_neighborhood)

    result = 0
    ivalue = vertex_property[vid]
    k=0
    if nb_neighborhood != 0 : # if ==0 it's mean that there is no neighbors for the vertex vid.
        for i in vid_neighborhood:
            if i in vertex_property.keys():
                result = result + vertex_property[i]
                k+=1
        if k!=0:
            return ivalue - (result / float(k))

@__normalized_parameters
def mean_abs_dev(graph, vertex_property, vid, rank, edge_type):
    """
    Sub-function computing the mean sum of absolute difference between ONE vertex ('vid') and its neighbors at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' : a vertex id.
    - 'rank' : neighborhood at distance 'rank' will be used.

    :Return:
    - a single value = the mean absolute deviation between vertex 'vid' and its neighbors at rank 'rank'.
    """
    if rank == 1:
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)
        vid_neighborhood.remove(vid)
    else: # if rank > 1, we want to compute the change only over the cell at `rank` and not for all cells between rank 1 and `rank`.
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)-graph.neighborhood(vid,rank-1, edge_type)

    nb_neighborhood = len(vid_neighborhood)

    result = 0
    ivalue = vertex_property[vid]
    k=0
    if nb_neighborhood != 0 : # if ==0 it's mean that there is no neighbors for the vertex vid.
        for i in vid_neighborhood:
            if i in vertex_property.keys():
                result = result + abs(ivalue - vertex_property[i])
                k+=1
        if k!=0:
            return result / float(k)


@__normalized_parameters
def mean_neigh(graph, vertex_property, vid, rank, edge_type):
    #"""
    #Sub-function computing the laplacian between ONE vertex ('vid') and its neighbors at rank 'rank'.
#
    #:Parameters:
    #- 'graph' : a TPG.
    #- 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    #- 'vid' : a vertex id.
    #- 'rank' : neighborhood at distance 'rank' will be used.
#
    #:Return:
    #- a single value = laplacian between vertex 'vid' and its neighbors at rank 'rank'.
    #"""
    if rank == 1:
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)
        vid_neighborhood.remove(vid)
    else: 
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)-graph.neighborhood(vid,rank-1, edge_type)

    nb_neighborhood = len(vid_neighborhood)

    result = 0
    ivalue = vertex_property[vid]
    k=0
    if nb_neighborhood != 0 : 
        for i in vid_neighborhood:
            if i in vertex_property.keys():
                result = result + vertex_property[i]
                k+=1
        if k!=0:
            return (ivalue + result) / float(k+1)


@__normalized_parameters
def change(graph, vertex_property, vid, rank, edge_type):
    """
    Sub-function computing the difference between ONE vertex ('vid') and its neighbors at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' : a vertex id.
    - 'rank' : neighborhood at distance 'rank' will be used.

    :Return:
    - a single value = laplacian between vertex 'vid' and its neighbors at rank 'rank'.
    """
    if rank == 1:
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)
        vid_neighborhood.remove(vid)
    else: # if rank > 1, we want to compute the change only over the cell at `rank` and not for all cells between rank 1 and `rank`.
        vid_neighborhood = graph.neighborhood(vid,rank, edge_type)-graph.neighborhood(vid,rank-1, edge_type)

    nb_neighborhood = len(vid_neighborhood)
    result = 0
    ivalue = vertex_property[vid]
    k=0
    if nb_neighborhood != 0 : # if ==0 it's mean that there is no neighbors for the vertex vid.
        for i in vid_neighborhood:
            if i in vertex_property.keys():
                result = result + vertex_property[i]
                k+=1
        if k!=0:
            return result/float(k) - ivalue


def __normalized_temporal_parameters(func):
    def wrapped_function(graph, vertex_property, vids = None, rank = 1, labels_at_t_n = False, check_exist_all_relative_at_rank = True, rank_lineage_check = None, verbose = False):
        """
        :Parameters:
        - 'graph' : a TPG.
        - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
        - 'vids' : by default a vertex id or a list of vertex ids. If 'vids=None' the function `func` will be computed for all ids present in the graph provided.
        - 'rank' : temporal neighborhood at distance 'rank' will be used.
        - 'rank_lineage_check' : usefull if you want to check the lineage for a different rank than the temporal neighborhood rank.

        :Example:
        VG12 = g.translate_keys_Graph2Image(relative_temporal_change(g, 'volume', rank = 1, rank_lineage_check = 4), 0 )
        VG15 = g.translate_keys_Graph2Image(relative_temporal_change(g, 'volume', rank = 4, rank_lineage_check = 4), 0 )
        This make sure that the same lineage is used for volumetric growth computation between t1-t2 and t1-t5.

        :Return:
        - a single value if vids is an interger, or a dictionnary of *keys=vids and *values= value computed by `func`
        """
        # -- If a name is given, we use vertex_property stored in the graph with this name.
        if isinstance(vertex_property,str):
            vertex_property = graph.vertex_property(vertex_property)

        # -- If no vids provided we compute the function for all keys present in the vertex_property
        no_warn = False
        if vids==None:
            vids = vertex_property.keys()
            no_warn = True
        if isinstance(vids,int):
            vids=[vids] # for single id, compute single result

        if rank<0:
            if not labels_at_t_n:
                rank = -rank
            else:
                raise ValueError("Rank should be positive if you want labels @t_n.")
        if rank_lineage_check == None:
            rank_lineage_check = rank

        try:
            graph.graph_property('time_steps')
        except:
            warnings.warn("You did not defined the `time_steps` when creating the TemporalPropertyGraph.")
            raise ValueError("Use graph.add_graph_property('time_steps',time_steps) to add it.")

        # -- For a list of ids, we create a dictionary of resulting values from temporal function `func`.
        temporal_func={}
        for n,vid in enumerate(vids):
            if verbose and n%10==0: print n,'/',len(vids)
            # -- We compute the `time_interval` each time in case `vids` have differents indexes:
            index_1 = graph.vertex_property('index')[vid]
            try:
                time_interval = graph.graph_property('time_steps')[index_1+rank]-graph.graph_property('time_steps')[index_1]
            except:
                continue #will go to the next `vid` if its 'index_1+rank' is out of range!
            if check_exist_all_relative_at_rank:
                if exist_all_relative_at_rank(graph, vid, rank_lineage_check): # Check if ALL descendants up to `rank_lineage_check` exists !
                    temporal_func[vid] = func(graph, vertex_property, vid, rank, time_interval)
            else:
                if exist_relative_at_rank(graph, vid, rank): # Check if there is at least one descendants for `vid` at `rank`!
                    temporal_func[vid] = func(graph, vertex_property, vid, rank, time_interval)

        # -- If there was any problems with some vis, we print it before returning the results:
        omit = list( set(vids) - set(temporal_func.keys()) )
        if omit != [] and not no_warn:
            if len(omit)<=20 or verbose:
                print "Some of the `vids` have been omitted: ", omit
            else:
                print "Some of the `vids` have been omitted."

        # -- Now we return the results of temporal differentiation function:
        if labels_at_t_n:
            return temporal_func
        else:
            return keys_to_daughter_id(graph, temporal_func, rank)

    return  wrapped_function


@__normalized_temporal_parameters
def temporal_rate(graph, vertex_property, vid, rank, time_interval):
    """
    Sub-function computing the temporal rate of change between ONE vertex ('vid') and its descendants at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' (str) : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' (int|list) : a vertex id.
    - 'rank' (int) : neighborhood at distance 'rank' will be used.

    :Return:
    - a single value = temporal change between vertex 'vid' and its neighbors at rank 'rank'.
    """
    if rank == 1:
        vid_descendants = graph.children(vid)
        vid_parent = vid
    # - If rank > 1, we want to compute the change only over the cell at `rank` and not for all cells between rank 1 and `rank`.
    if rank > 1:
        vid_descendants = graph.descendants(vid,rank)-graph.descendants(vid,rank-1)
        vid_parent = vid

    try:
        parent_value = vertex_property[vid_parent]
    except KeyError as e:
        print KeyError("temporal_rate "+"error 1 "+str(vid_parent))
        return np.nan
    try:
        descendants_value = sum([vertex_property[id_descendant] for id_descendant in vid_descendants])
    except KeyError as e:
        print KeyError("temporal_rate "+"error 2 "+str(e))
        return np.nan

    return (descendants_value / parent_value) * 1. / float(time_interval)


@__normalized_temporal_parameters
def log_temporal_rate(graph, vertex_property, vid, rank, time_interval):
    """
    Sub-function computing the log temporal rate of change between ONE vertex ('vid') and its descendants at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' (str) : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' (int|list) : a vertex id.
    - 'rank' (int) : neighborhood at distance 'rank' will be used.

    :Return:
    - a single value = temporal change between vertex 'vid' and its neighbors at rank 'rank'.
    """
    if rank == 1:
        vid_descendants = graph.children(vid)
        vid_parent = vid
    # - If rank > 1, we want to compute the change only over the cell at `rank` and not for all cells between rank 1 and `rank`.
    if rank > 1:
        vid_descendants = graph.descendants(vid,rank)-graph.descendants(vid,rank-1)
        vid_parent = vid

    try:
        parent_value = vertex_property[vid_parent]
    except KeyError as e:
        print KeyError("log_temporal_rate "+"error 1 "+str(vid_parent))
        return np.nan
    try:
        descendants_value = sum([vertex_property[id_descendant] for id_descendant in vid_descendants])
    except KeyError as e:
        print KeyError("log_temporal_rate "+"error 2 "+str(e))
        return np.nan

    return np.log2(descendants_value / parent_value) * 1. / float(time_interval)


@__normalized_temporal_parameters
def temporal_change(graph, vertex_property, vid, rank, time_interval):
    """
    Sub-function computing the temporal change between ONE vertex ('vid') and its descendants at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' : a vertex id.
    - 'rank' : neighborhood at distance 'rank' will be used.

    :Return:
    - a single value = temporal change between vertex 'vid' and its neighbors at rank 'rank'.
    """
    if rank == 1:
        vid_descendants = graph.children(vid)
        vid_parent = vid
    # - If rank > 1, we want to compute the change only over the cell at `rank` and not for all cells between rank 1 and `rank`.
    if rank > 1:
        vid_descendants = graph.descendants(vid,rank)-graph.descendants(vid,rank-1)
        vid_parent = vid

    try:
        parent_value = vertex_property[vid_parent]
    except KeyError as e:
        print KeyError("temporal_rate "+"error 1 "+str(vid_parent))
        return np.nan
    try:
        descendants_value = sum([vertex_property[id_descendant] for id_descendant in vid_descendants])
    except KeyError as e:
        print KeyError("temporal_rate "+"error 2 "+str(e))
        return np.nan

    return (descendants_value - parent_value) / float(time_interval)


@__normalized_temporal_parameters
def relative_temporal_change(graph, vertex_property, vid, rank, time_interval):
    """
    Sub-function computing the relative temporal change between ONE vertex ('vid') and its descendants at rank 'rank'.

    :Parameters:
    - 'graph' : a TPG.
    - 'vertex_property' : the dictionnary TPG.vertex_property('property-of-interest'), or the string 'property-of-interest'.
    - 'vid' : a vertex id.
    - 'rank' : neighborhood at distance 'rank' will be used.

    :Return:
    - a single value = relative temporal change between vertex 'vid' and its neighbors at rank 'rank'.
    """
    return temporal_change(graph, vertex_property, vid, rank, time_interval).values()[0] / float(vertex_property[vid])


def shape_anisotropy_2D(graph):
    """
    Sub-function computing the shape anisotropy of one cell using.

    :Parameters:
     - 'graph' (TGP) - a TPG.
    :Return:
     - shape_anisotropy (dict) .
    """
    assert len(graph.vertex_property('inertia_values').values()[0])==2
    return dict([ (vid, (inertia[0]-inertia[1])/(inertia[0]+inertia[1])) for vid, inertia in graph.vertex_property('inertia_values').iteritems() ])


def fractional_anisotropy(eigenvalues):
    """
    Compute fractional anisotropy of a tensor considered to represent a diffusion ellipsoid.
    $$    \text{FA} = \sqrt{\frac{3}{2}} \frac{\sqrt{(\lambda_1 - \hat{\lambda})^2 + (\lambda_2 - \hat{\lambda})^2 + (\lambda_3 - \hat{\lambda})^2}}{\sqrt{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}$$
    with the trace $\hat{\lambda} = (\lambda_1 + \lambda_2 + \lambda_3)/3$
    :!Can not be negative!:
    """
    assert len(eigenvalues)==3
    l1, l2, l3 = eigenvalues
    if (l1+l2+l3) == 0: return None
    l = (l1+l2+l3)/3.
    return math.sqrt(3./2.) * (math.sqrt( (l1-l)**2+(l2-l)**2+(l3-l)**2 )/math.sqrt(l1**2+l2**2+l3**2))


def shape_anisotropy_3D(graph):
    """
    Sub-function computing the shape anisotropy of one cell using the Fractional anisotropy scalar

    :Parameters:
     - 'graph' (TGP) - a TPG.
    :Return:
     - shape_anisotropy = temporal division rate between vertex 'vid' and its descendants at rank 'rank'.
    """
    return dict([ (vid, fractional_anisotropy(inertia)) for vid, inertia in graph.vertex_property('inertia_values').iteritems() ])


def reduced_shape_anisotropy_3D(graph):
    """
    Sub-function computing the shape anisotropy of one cell using the Fractional anisotropy scalar

    :Parameters:
     - 'graph' (TGP) - a TPG.
    :Return:
     - shape_anisotropy = temporal division rate between vertex 'vid' and its descendants at rank 'rank'.
    """
    return dict([ (vid, fractional_anisotropy(inertia)) for vid, inertia in graph.vertex_property('reduced_inertia_values').iteritems() ])


def epidermis_wall_gaussian_curvature(graph):
    """
    Use the graph vertex property `epidermis_wall_principal_curvature_value` to compute the Gaussian curvature.
    The principal curvature values saved there are based only on wall voxels
    Gaussian curvature is the product of principal curvatures 'k1*k2'.
    """
    return dict([ (vid, curv_values[0] * curv_values[1]) for vid, curv_values in graph.vertex_property('epidermis_wall_principal_curvature_values').iteritems()])


def epidermis_local_gaussian_curvature(graph, radius):
    """
    Use the graph vertex property `epidermis_wall_principal_curvature_value` to compute the Gaussian curvature.
    The principal curvature values saved there are based on voxels present in a within a certain radius around the wall median.
    Gaussian curvature is the product of principal curvatures 'k1*k2'.
    """
    #assert radius in graph.graph_property('radius_local_principal_curvature_estimation')
    return dict([ (vid, curv_values[0] * curv_values[1]) for vid, curv_values in graph.vertex_property('epidermis_local_principal_curvature_values_r{}'.format(radius)).iteritems()])


def epidermis_local_curvature_ratio(graph, radius):
    """
    Use the graph vertex property `epidermis_wall_principal_curvature_value` to compute the Gaussian curvature.
    The principal curvature values saved there are based on voxels present in a within a certain radius around the wall median.
    Gaussian curvature is the product of principal curvatures 'k1*k2'.
    """
    #assert radius in graph.graph_property('radius_local_principal_curvature_estimation')
    return dict([ (vid, curv_values[0] / curv_values[1]) for vid, curv_values in graph.vertex_property('epidermis_local_principal_curvature_values_r{}'.format(radius)).iteritems()])


def division_rate(graph, rank=1, labels_at_t_n = False):
    """
    Division rate
    :Parameters:
     - 'graph' (TGP) - a TPG.
     - 'rank' (int) - children at distance 'rank' will be used.
     - 'labels_at_t_n' (bool) - specify if the division rate values returned should be associated with parent ids or children ids.

    :Return:
     - div_rate = temporal division rate between vertex 'vid' and its descendants at rank 'rank'.
    """
    if (rank > 1) and (labels_at_t_n is False):
        raise ValueError("The translation function `translate_keys2daughters_ids` doesn't work for rank != 1.")

    div_rate = {}
    for vid in graph.vertices():
        index_1 = graph.vertex_property('index')[vid]
        try:
            time_interval = graph.graph_property('time_steps')[index_1+rank]-graph.graph_property('time_steps')[index_1]
        except:
            continue #will go to the next `vid` if its 'index_1+rank' is out of range!

        descendants=graph.descendants(vid,rank)-set([vid])
        if descendants == set([]):
            div_rate[vid] = 1 / float(time_interval)
        else:
            div_rate[vid] = len(descendants) / float(time_interval)

    if labels_at_t_n:
        return div_rate
    else:
        return translate_keys2daughters_ids(graph, div_rate)


def exist_relative_at_rank(graph, vid, rank):
    """
    Check if there is a relative (descendant or ancestor) of the 'vid' at 'rank'.

    :Parameters:
     - 'graph' (TPG): the TPG to be used for translation
     - 'vid' (int): the initial point to look-out for rank existence.
     - 'rank' (int): the rank to test.
    """
    if rank == 0 :
        return True
    if (rank > 0) :
        try :
            return graph.descendants(vid,rank)-graph.descendants(vid,rank-1) != set()
        except:
            return False
    if (rank < 0) :
        try :
            return graph.ancestors(vid,abs(rank))-graph.ancestors(vid,abs(rank)-1) != set()
        except:
            return False


def exist_all_relative_at_rank(graph, vid, rank):
    """
    Check if lineage is complete over several ranks.
    i.e. every decendants cells from `vid` have a lineage up to rank `rank`.

    :Parameters:
     - `graph` TPG to browse;
     - 'vid' : a vertex id.
     - 'rank' : neighborhood at distance 'rank' will be used.
    """
    if rank == 0 or not exist_relative_at_rank(graph, vid, rank):
        return False

    if rank == 1 or rank == -1:
        return exist_relative_at_rank(graph, vid, rank)

    if (rank > 0):
        descendants_at_rank = {}
        descendants_at_rank[1] = graph.children(vid)
        for r in xrange(2,rank+1):
            for v in descendants_at_rank[r-1]:
                if graph.children(v) == set():
                    return False
                if descendants_at_rank.has_key(r):
                    descendants_at_rank[r].update(graph.children(v))
                else:
                    descendants_at_rank[r] = graph.children(v)
        return True

    if (rank < 0):
        rank = -rank
        descendants_at_rank = {}
        descendants_at_rank[1] = graph.parent(vid)
        for r in xrange(2,rank+1):
            for v in descendants_at_rank[r-1]:
                if graph.parent(v) == set():
                    return False
                if descendants_at_rank.has_key(r):
                    descendants_at_rank[r].update(graph.parent(v))
                else:
                    descendants_at_rank[r] = graph.parent(v)
        return True


def time_point_property(graph, time_point, vertex_property, lineaged=False, fully_lineaged=False, as_parent=False, as_children=False):
    """
    Allow to extract a property 'vertex_property' from the temporal graph for one time-point.

    :Parameters:
    - `graph` (TemporalPropertyGraph) - Spatio-temporal graph to browse;
    - `time_point` (int) - define the time-point to consider;
    - `vertex_property` (str) - name of the vertex property to extract;
    :Return:
    - dictionnary of vertex property extracted from the time-point 'time_point';
    """
    # if a name is given, we use vertex_property stored in the graph with this name.
    if isinstance(vertex_property,str):
        vertex_property = graph.vertex_property(vertex_property)

    if time_point not in graph.vertex_property('index').values():
        raise ValueError(str(time_point)+"not in"+str(graph))

    vids_at_time = graph.vertex_at_time(time_point, lineaged, fully_lineaged, as_parent, as_children)

    return dict([(i,vertex_property[i]) for i in vids_at_time if vertex_property.has_key(i)])


def time_point_property_by_regions(graph, time_point, vertex_property, lineaged=False, fully_lineaged=False, as_parent=False, as_children=False):
    """
    Allow to extract a property 'vertex_property' from the temporal graph for one time-point, sorted by regions.
    Return a dict of dict, first level of keys are region(s) name(s) and second layer are vertex ids and the values of their associated property.

    :Parameters:
    - `graph` (TemporalPropertyGraph) - Spatio-temporal graph to browse;
    - `time_point` (int) - define the time-point to consider;
    - `vertex_property` (str) - name of the vertex property to extract;
    :Return:
    - dictionary of regions which values are a dictionnary of vertex property extracted from the time-point;
    """
    extracted_property = time_point_property(graph, time_point, vertex_property, lineaged, fully_lineaged, as_parent, as_children)

    regions_names = list(np.unique([v[0] for k,v in graph.vertex_property('regions').iteritems() if graph.vertex_property('index')[k]==time_point]))

    property_by_regions = {}
    for region_name in regions_names:
        property_by_regions[region_name] = dict( (k,v) for k,v in extracted_property.iteritems() if graph.vertex_property('regions').has_key(k) and (graph.vertex_property('regions')[k][0]==region_name) )

    return property_by_regions


def extend_margins(mini, maxi, percent=0.04):
    """
    Extend the mini and maxi value by a percentage of their difference.
    Used for display purpose when defining axes range so mini and maxi values are well separated form the surrounding box.
    """
    assert mini < maxi
    diff = maxi-mini
    pc = diff*percent

    if mini-pc<0 and mini>0:
        mini = 0
    else:
        mini = mini-pc
    if maxi+pc>0 and maxi<0:
        maxi = 0
    else:
        maxi = maxi+pc

    return mini, maxi


def histogram_property_by_time_points(graph, vertex_property, time_points=None, **kwargs):
    """
    Display an histogram or barplot of the provided `vertex_property` by `time_points`.
    """
    # Handle initial data to work with:
    property_name = None
    if isinstance(vertex_property,str):
        assert vertex_property in list(graph.vertex_properties())
        property_name = vertex_property
        vertex_property = graph.vertex_property(vertex_property)
    if time_points is None:
        time_points = range(graph.nb_time_points+1)

    # kwargs associated to graph properties:
    ppt_kwargs = {}
    if 'lineaged' in kwargs: ppt_kwargs.update({'lineaged':kwargs['lineaged']})
    if 'as_parent' in kwargs: ppt_kwargs.update({'as_parent':kwargs['as_parent']})
    if 'as_children' in kwargs: ppt_kwargs.update({'as_children':kwargs['as_children']})
    vids, data = [], []
    for tp in time_points:
        ppty_dict = time_point_property(graph, tp, property_name, **ppt_kwargs)
        vids.append(np.array(ppty_dict.keys()))
        data.append(np.array(ppty_dict.values()))

    # Handle kwargs passed to 'plt.hist':
    h_kwargs= {}
    if 'bins' in kwargs: h_kwargs.update({'bins':kwargs['bins']})
    if 'range' in kwargs: h_kwargs.update({'range':kwargs['range']})
    if 'normed' in kwargs: h_kwargs.update({'normed':kwargs['normed']})
    if 'histtype' in kwargs: h_kwargs.update({'histtype':kwargs['histtype']})
    # Other kwargs:
    xlim = kwargs['xlim'] if 'xlim' in kwargs else None
    outliers = kwargs['outliers'] if 'outliers' in kwargs else None# if a type 'dict', should be *key=vid:*values=True/False; if a type 'int', will serve as threshold for MAD estimator
    thres = kwargs['threshold'] if 'threshold' in kwargs else np.nan
    sidebyside = kwargs['sidebyside'] if 'sidebyside' in kwargs else False

    # Handling "outliers" infos:
    if isinstance(outliers, int):
        #~ print 'Outliers threshold provided... selecting outliers by MAD estimator !'
        from sa_oa.container.graph_clusterer import mad_based_outlier
        thres = copy.copy(outliers)
        all_outliers = mad_based_outlier(vertex_property,thres)
        outliers = []
        for tp in time_points:
            tp_ids = [i for i in graph.vertex_at_time(tp) if i in all_outliers.keys()]
            outliers.append([vertex_property[k] for k in tp_ids if all_outliers[k]])
    elif isinstance(outliers, dict):
        #~ print 'True/False dictionary of outliers provided !'
        outliers_dict = copy.copy(outliers)
        outliers = []
        for tp in time_points:
            outliers.append( [data[tp][n] for n,vid in enumerate(vids[tp]) if outliers_dict[vid]] )
    else:
        #~ print 'No outliers infos !'
        pass

    # Define margins of outliers:
    if outliers is not None:
        p_data, n_data = [i for k in outliers if k!=[] for i in k if i>0], [i for k in outliers if k!=[] for i in k if i<0]
        if p_data != []:
            min_out_p, max_out_p = min(p_data), max(p_data)
        else:
            min_out_p, max_out_p = None, None
        if n_data != []:
            min_out_n, max_out_n = min(n_data), max(n_data)
        else:
            min_out_n, max_out_n = None, None

    # -- Initialisation of the matplotlib figure:
    fsize = [14,4] if sidebyside else [14,8]
    fig = plt.figure(figsize=fsize, dpi=80)

    # - First make an histogram of the data:
    # Create a subfigure and the histogram:
    histo = fig.add_subplot(121) if sidebyside else fig.add_subplot(211)
    n, bins, patches = histo.hist(data, label = ["time point #{}".format(tp) for tp in time_points], rwidth=1., **h_kwargs)
    # Add x, y axis labels:
    if kwargs.has_key('xlabel'):
        histo.set_xlabel(kwargs['xlabel'])
    elif property_name is not None:
        if graph.graph_property("units").has_key(property_name):
            histo.set_xlabel(property_name+" ("+graph.graph_property("units")[property_name]+")",family='freesans')
        else:
            histo.set_xlabel(property_name)
    if h_kwargs.has_key('normed') and h_kwargs['normed']:
        histo.set_ylabel("Relative Frequency")
    elif h_kwargs.has_key('normed') and not h_kwargs['normed']:
        histo.set_ylabel("Frequency")
    # Add outliers positions if available:
    try:
        histo.vlines([min_out_p, max_out_p, min_out_n, max_out_n], ymin=0, ymax=np.max(n)/2.5, color='r', linestyles='dashed', label='outliers', hold=True)
    except:
        pass
    if xlim is not None: histo.set_xlim(xlim[0],xlim[1])
    # Add a legend:
    histo.legend(loc=0, framealpha=0.7, fontsize='small')

    # - Then make a boxplot of the data:
    # Create a subfigure and the boxplot:
    bp = fig.add_subplot(122) if sidebyside else fig.add_subplot(212)
    bp.boxplot(data, vert=0, positions=time_points)
    # Add x, y axis labels:
    bp.set_ylabel('Time points')
    if kwargs.has_key('xlabel'):
        bp.set_xlabel(kwargs['xlabel'])
    elif property_name is not None:
        if graph.graph_property("units").has_key(property_name):
            bp.set_xlabel(property_name+" ("+graph.graph_property("units")[property_name]+")",family='freesans')
        else:
            bp.set_xlabel(property_name)
    # Add outliers positions if available:
    if isinstance(outliers, list) and len(outliers)==len(time_points):
        for tp, outs in enumerate(outliers):
            plt.plot(outs, np.zeros_like(outs)+tp, 'ro', scalex=False, label='outliers' if tp==0 else None)
    elif outliers is not None and min_out is not None:
        bp.vlines([min_out_p, max_out_p, min_out_n, max_out_n], ymin=min(time_points), ymax=max(time_points), color='r', linestyles='dashed', label='outliers', hold=True)
    else:
        pass
    # Add outliers detection method infos if availables:
    if outliers is not None:
        MAD_kwargs = dict(y=0.01, x=0.99, ha='right', va='bottom', fontsize='small')
        bp.axes.set_title('MAD estimator threshold={}'.format(thres), **MAD_kwargs)
    if xlim is not None: bp.set_xlim(xlim[0],xlim[1])
    # Add a legend:
    bp.legend(loc=1, framealpha=0.7, fontsize='small', numpoints=1)

    # Add a global title to the figure:
    if kwargs.has_key('title'):
        plt.suptitle(kwargs['title'])
    elif property_name is not None:
        plt.suptitle("Histogram and Boxplot of {} property".format(property_name))
    # Make a thight layout around the figures:
    plt.tight_layout()

    return n, bins, patches

def boxplot_property_by_time_points_and_regions(graph, vertex_property, regions, remove_outliers=None, **kwargs):
    """
    Display an histogram or barplot of the provided `vertex_property` by `time_points`.
    :Parameters:
     - `graph` (temporal_property_graph) - represent the spatio-temporal relations between clusters
     - `vertex_property` (str|dict) - string matching a vertex_property known to the graph or a dictionary
     - `regions` (dict) - dictionary sorting vertex_ids (vids) into groups, *keys=vids: *values=group_ids
     - `remove_outliers` (int) - if not None, will detect outliers according to given threshold and remove them from boxplot display 
    
    """
    # Handle initial data to work with:
    property_name = None
    if isinstance(vertex_property,str):
        assert vertex_property in list(graph.vertex_properties())
        property_name = vertex_property
        vertex_property = graph.vertex_property(vertex_property)

    temporal = kwargs['temporal'] if 'temporal' in kwargs else False
    time_points = range(graph.nb_time_points) if temporal else range(graph.nb_time_points+1) 
    if temporal:
        vertex_property = keys_to_mother_id(graph, vertex_property)

    if remove_outliers is not None:
        from basics import mad_based_outlier
        outliers = mad_based_outlier(vertex_property,remove_outliers)
    else:
        outliers = dict([(k,False) for k in vertex_property])

    fully_lineaged_vtx = graph.lineaged_vertex(fully_lineaged=True)
    tp_index_vtx = graph.vertex_property('index')
    data = {}; mini, maxi = np.inf, -np.inf
    for tp in time_points:
        ppty_dict = dict([(k,v) for k,v in vertex_property.iteritems() if (tp_index_vtx[k]==tp)and(k in fully_lineaged_vtx)and(not outliers[k])])
        mini = np.min([mini,np.min(ppty_dict.values())])
        maxi = np.max([maxi,np.max(ppty_dict.values())])
        data[tp] = {}
        for vid, region in regions.iteritems():
            if ppty_dict.has_key(vid):
                if data[tp].has_key(region):
                    data[tp][region].append(ppty_dict[vid])
                else:
                    data[tp][region]= [ppty_dict[vid]]

    regions_ids = np.unique(regions.values())
    N_regions = len(regions_ids)
    N_box = len(time_points)-1

    bp_mini, bp_maxi = extend_margins(mini, maxi, 0.04)
    # kwargs associated to graph properties:
    bp_kwargs = {}
    bp_kwargs['range'] = kwargs['range'] if 'range' in kwargs else [bp_mini, bp_maxi]

    # -- Initialisation of the matplotlib figure:
    fsize = [5*N_box,7]
    fig = plt.figure(figsize=fsize, dpi=80)
    bp = {}
    # - Then make a boxplot of the data:
    for tp in time_points[1:]:
        # Create a subfigure and the boxplot:
        bp[tp] = fig.add_subplot(1,N_box,tp)
        tp_data = [data[tp][region] for region in regions_ids]
        bp[tp].boxplot(tp_data, vert=1, positions=range(N_regions))
        plt.plot(range(N_regions), [np.mean(data[0][r]) if data[0].has_key(r) else [None] for r in regions_ids] , 'rx', scalex=False, label='Initial value')
        if tp >=2:
            plt.plot(range(N_regions), [np.mean(data[tp-1][r]) for r in regions_ids], 'gx', scalex=False, label='t_n-1 mean')
        bp[tp].set_ylim(tuple(bp_kwargs['range']))
        # Add x, y axis labels:
        bp[tp].set_xlabel('regions')
        if kwargs.has_key('xlabel'):
            bp[tp].set_ylabel(kwargs['xlabel'])
        elif property_name is not None:
            if graph.graph_property("units").has_key(property_name):
                bp[tp].set_ylabel(property_name+" ("+graph.graph_property("units")[property_name]+")",family='freesans')
            else:
                bp[tp].set_ylabel(property_name)

        if temporal:
            plt.title('t{} -> t{}'.format(tp, tp+1))
        else:
            plt.title('Time point {}'.format(tp))
    # Add a legend:
    leg = plt.legend(loc=2, framealpha=0.7, fontsize='small')
    # Add outliers detection method infos if availables:
    if remove_outliers is not None:
        MAD_kwargs = dict(y=mini, x=0.1, ha='left', va='top', fontsize='small')
        bp[time_points[1]].text(s='MAD estimator threshold={}'.format(remove_outliers), **MAD_kwargs)
    # Add a global title to the figure:
    if kwargs.has_key('title'):
        plt.suptitle(kwargs['title'])
    elif property_name is not None:
        plt.suptitle("Boxplot of {} property by time points and regions{}".format(property_name, ' (without outliers)' if remove_outliers is not None else ''))
    # Make a thight layout around the figures:
    plt.tight_layout()

    return mini,maxi,'Done!'



def translate_keys2daughters_ids(graph, dictionary):
    """
    Translate keys of a dictionary to daughter ids according to the graph (Temporal Property Graph).
    """
    dictionary_daughters={}
    no_descendants=[]
    for vid in dictionary:
        vid_descendants = graph.descendants(vid ,1)-graph.descendants(vid,0)
        if vid_descendants != set():
            for id_descendant in vid_descendants:
                dictionary_daughters[id_descendant]=dictionary[vid]
        elif graph.vertex_property('index')[vid] < graph.nb_time_points-1: #if `vid``belong to the last time point, it's perfectly normal that there is no descendants
            no_descendants.append(vid)

    if no_descendants!=[]:
        warnings.warn("No daughter found for those vertex:"+str(no_descendants))

    return dictionary_daughters


def translate_list2daughters_ids(graph, ids_list):
    """
    Translate keys of a dictionary to daughter ids according to the graph (Temporal Property Graph).
    """
    list_daughters=[]
    no_descendants=[]
    for vid in ids_list:
        vid_descendants = graph.descendants(vid ,1)-graph.descendants(vid,0)
        if vid_descendants != set():
            for id_descendant in vid_descendants:
                list_daughters.append(id_descendant)
        elif graph.vertex_property('index')[vid] < graph.nb_time_points-1: #if `vid``belong to the last time point, it's perfectly normal that there is no descendants
            no_descendants.append(vid)

    if no_descendants!=[]:
        warnings.warn("No daughter found for those vertex:"+str(no_descendants))

    return list_daughters


def weighted_mean( values, weights ):
    """
    Function computing a weighted mean of `values` according to `weights`.
         - `values` (list)
         - `weights` (list)
    Return:
         - wm = sum_i( value_i * weight_i/sum(weights) )
    """
    assert len(values)==len(weights)

    if isinstance(values[0],int) or isinstance(values[0],np.ndarray):
        return sum( [v * w / float(sum(weights)) for v,w in zip(values, weights)] )
    if isinstance(values[0],list) or isinstance(values[0],tuple):
        return sum( [np.array(v) * w / float(sum(weights)) for v,w in zip(values, weights)] )


def __strain_parameters(func):
    def wrapped_function(graph, vids = None, labels_at_t_n = True, use_projected_anticlinal_wall = False, verbose = False):
        """
        :Parameters:
         - `graph` (TPG)
         - `vids` (int|list) - list of (graph) vertex ids @ t_n. If :None: it will be computed for all possible vertex
         - `use_projected_anticlinal_wall` (bool) - if True, use the medians of projected anticlinal wall ('projected_anticlinal_wall_median') to compute a strain in 2D.
        """
        # Check if cells landmarks have been recovered ans stored in the graph structure:
        if use_projected_anticlinal_wall:
            assert 'projected_anticlinal_wall_median' in graph.edge_property_names()
            assert 'epidermis_wall_median' in graph.vertex_property_names()
            assert 'epidermis_surface' in graph.vertex_property_names()
        else:
            assert 'wall_median' in graph.edge_property_names()
            assert 'unlabelled_wall_median' in graph.vertex_property_names()
            assert 'epidermis_wall_median' in graph.vertex_property_names()

        # -- If the vid is not associated through time it's not possible to compute the strain.
        if vids is None:
            vids = list(graph.vertices())
        if isinstance(vids,int):
            vids=[vids]

        # -- If the vid is not associated through time it's not possible to compute the strain.
        tmp = copy.copy(vids)
        for vid in vids:
            if graph.descendants(vid, 1) == set([vid]):
                tmp.remove(vid)

        vids = tmp
        if use_projected_anticlinal_wall:
            spatial_vertices_edge = dict( [(tuple(sorted(v)),k) for k,v in graph._edges.iteritems() if k in graph.edges(edge_type='s')] )
            # Next line could be changed: we could create a smaller dict by looking in something smaller than all graphs edges!!!
            wall_median_2 = dict( ( tuple(sorted(graph.edge_vertices(eid))), graph.edge_property('projected_anticlinal_wall_median')[eid]) for eid in graph.edges(edge_type='s') if graph.edge_property('projected_anticlinal_wall_median').has_key(eid) )
            stretch_mat = {}

            for n,vid in enumerate(vids):
                if verbose and n%10==0: print n, " / ", len(vids)
                missing_data = False
                spatial_edges = list(graph.edges( vid, 's' ))
                # We create the dictionary of medians associate over time: {(label_1,label_2):[x,y,z]} with label_1<label_2 and label_1&label_2 in the first layer !
                wall_median = dict( ( tuple(sorted(graph.edge_vertices(eid))), graph.edge_property('projected_anticlinal_wall_median')[eid] ) for eid in spatial_edges if graph.edge_property('projected_anticlinal_wall_median').has_key(eid) )

                # -- Now we want to associate landmarks over time:
                xyz_t1, xyz_t2 = [], []
                descendants_vid = list(graph.descendants(vid,1)-set([vid]))
                neighbors_t1 = [nei for nei in list(graph.neighbors(vid,'s')) if wall_median.has_key(tuple(sorted([nei,vid])))]
                # - Using barycenter of the epidermis wall :
                xyz_t1.append(np.array(graph.vertex_property('epidermis_wall_median')[vid]))
                if len(descendants_vid)>1:
                    xyz_t2.append(weighted_mean( [graph.vertex_property('epidermis_wall_median')[desc_id_cp] for desc_id_cp in descendants_vid], [graph.vertex_property('epidermis_surface')[desc_id_cp] for desc_id_cp in descendants_vid if graph.vertex_property('epidermis_surface').has_key(desc_id_cp)] ))
                else:
                    xyz_t2.append(np.array(graph.vertex_property('epidermis_wall_median')[descendants_vid[0]]))

                # - Now we need to find the time correspondance between each median:
                for neighbor in neighbors_t1:
                    descendants_neighbor = [desc_nei for desc_nei in graph.descendants(neighbor,1)-set([neighbor])]
                    nei = list(set(descendants_vid) | set(descendants_neighbor))
                    # - In order to find the decendants that share the 'same' wall we have to find those from each (two) groups that have a topological distance of 1:
                    topological_distance = dict( ((vtx_id, graph.topological_distance(vtx_id, 's', full_dict=False))) for vtx_id in nei )
                    neighbors_descendants = []
                    for nei_1 in descendants_vid:
                        for nei_2 in descendants_neighbor:
                            if wall_median_2.has_key(tuple(sorted((nei_1,nei_2)))) and topological_distance[nei_1][nei_2] == 1:
                                neighbors_descendants.append(tuple(sorted((nei_1,nei_2))))
                    # - If correspondance is found we retreive landmarks coordiantes
                    if neighbors_descendants != [] and len([wall_median_2[desc_id_cp] for desc_id_cp in neighbors_descendants])==len(neighbors_descendants):
                        id_couple = (min(vid,neighbor),max(vid,neighbor))
                        xyz_t1.append(np.array(wall_median[id_couple]))
                        # - We compute a weighted average for the t_n+1 landmark if the cell has divided: x,y,z = sum_i{ wall_surface_i * [x_i,y_i,z_i] }
                        xyz_t2.append(weighted_mean( [wall_median_2[desc_id_cp] for desc_id_cp in neighbors_descendants], [graph.edge_property('wall_surface')[spatial_vertices_edge[desc_id_cp]] for desc_id_cp in neighbors_descendants if graph.edge_property('wall_surface').has_key(spatial_vertices_edge[desc_id_cp])] ))

                if not missing_data:
                    assert len(xyz_t1) == len(xyz_t2)
                    N = len(xyz_t1)
                    stretch_mat[vid] = func(graph, xyz_t1, xyz_t2)
        else:
            ############################################################
            # If graph.vertex_property('unlabelled_wall_median')[vid] is None, it means we should not compute a 3D strain for this cell !!!
            # It's beacause unlabelled walls were not contiguous, so it could not be used to compute a median !!
            ############################################################
            spatial_vertices_edge = dict( ((min(v[0],v[1]),max(v[0],v[1])),k) for k,v in graph._edges.iteritems() if k in graph.edges(edge_type='s'))
            # Next line need to be changed: we should create a smaller dict by looking in something smaller than all graphs edges!!!
            wall_median_2 = dict( ( (min(graph.edge_vertices(eid)[0],graph.edge_vertices(eid)[1]),max(graph.edge_vertices(eid)[0],graph.edge_vertices(eid)[1])) ,graph.edge_property('wall_median')[eid]) for eid in graph.edges(edge_type='s') if graph.edge_property('wall_median').has_key(eid) )
            stretch_mat = {}
            for n,vid in enumerate(vids):
                if verbose and n%10==0: print n, " / ", len(vids)
                missing_data=False
                spatial_out_edges = list(graph.edges( vid, 's' ))
                wall_median = dict( ((min(graph.edge_vertices(eid)[0],graph.edge_vertices(eid)[1]),max(graph.edge_vertices(eid)[0],graph.edge_vertices(eid)[1])),graph.edge_property('wall_median')[eid]) for eid in spatial_out_edges if graph.edge_property('wall_median').has_key(eid) )

                # - If not every 's'type edges from a `vid` have a median and there is no 'unlabelled_wall_median' associated, it means we don't have all the info for strain computation:
                if (len(wall_median) == len(spatial_out_edges)) or (graph.vertex_property('unlabelled_wall_median').has_key(vid)):
                    xyz_t1, xyz_t2 = [], []
                    descendants_vid = graph.descendants(vid,1)-set([vid])
                    neighbors_t1 = list(graph.neighbors(vid,'s'))
                    # -- Now we need to find the time correspondance between each median:
                    for neighbor in neighbors_t1:
                        descendants_neighbor = graph.descendants(neighbor,1)-set([neighbor])
                        nei = list(descendants_vid | descendants_neighbor)
                        # -- In order to find the decendants that share the 'same' wall we have to find those from each (two) groups that have a topological distance of 1:
                        topological_distance = dict( ((vtx_id, graph.topological_distance(vtx_id, 's', full_dict=False))) for vtx_id in nei )
                        neighbors_descendants = []
                        for nei_1 in descendants_vid:
                            for nei_2 in descendants_neighbor:
                                if topological_distance[nei_1][nei_2] == 1:
                                    if nei_1 < nei_2: neighbors_descendants.append((nei_1,nei_2))
                                    if nei_1 > nei_2: neighbors_descendants.append((nei_2,nei_1))

                        if neighbors_descendants != [] and len([wall_median_2[desc_id_cp] for desc_id_cp in neighbors_descendants])==len(neighbors_descendants):
                            id_couple = (min(vid,neighbor),max(vid,neighbor))
                            xyz_t1.append(np.array(wall_median[id_couple]))
                            xyz_t2.append(weighted_mean( [wall_median_2[desc_id_cp] for desc_id_cp in neighbors_descendants], [graph.edge_property('wall_surface')[spatial_vertices_edge[desc_id_cp]] for desc_id_cp in neighbors_descendants] ))
                        else:
                            missing_data=True

                    # -- We need to make sure `unlabelled_wall_median` is both @ t_n and at least for one daughter cell @ t_n+1
                    if not missing_data and graph.vertex_property('unlabelled_wall_median').has_key(vid) and (sum([graph.vertex_property('unlabelled_wall_median').has_key(id2) for id2 in descendants_vid])>=1):
                        xyz_t1.append(graph.vertex_property('unlabelled_wall_median')[vid])
                        xyz_t2.append( weighted_mean( [graph.vertex_property('unlabelled_wall_median')[desc_vid] for desc_vid in descendants_vid], [graph.vertex_property('unlabelled_wall_surface')[desc_id] for desc_id in descendants_vid] ))

                    # -- We need to make sure `epidermis_wall_median` is both @ t_n and at least for one daughter cell @ t_n+1
                    if not missing_data and graph.vertex_property('epidermis_wall_median').has_key(vid) and (sum([graph.vertex_property('epidermis_wall_median').has_key(id2) for id2 in descendants_vid])>=1):
                        xyz_t1.append(graph.vertex_property('epidermis_wall_median')[vid])
                        xyz_t2.append( weighted_mean( [graph.vertex_property('epidermis_wall_median')[desc_vid] for desc_vid in descendants_vid], [graph.vertex_property('epidermis_surface')[desc_id] for desc_id in descendants_vid] ))

                if not missing_data:
                    assert len(xyz_t1) == len(xyz_t2)
                    N = len(xyz_t1)
                    stretch_mat[vid] = func(graph, xyz_t1, xyz_t2)

        # -- Now we return the results of temporal differentiation function:
        if labels_at_t_n:
            return stretch_mat
        else:
            stretch_mat_daughters={}
            print "You have asked for labels @ t_n+1"
            for vid in stretch_mat:
                vid_descendants = graph.descendants(vid ,1)-graph.descendants(vid,0)
                for id_descendant in vid_descendants:
                    stretch_mat_daughters[id_descendant]=stretch_mat[vid]
            return stretch_mat_daughters

    return  wrapped_function

def __strain_parameters2(func):
    def wrapped_function(graph, vids = None, labels_at_t_n = True, use_projected_anticlinal_wall = False, verbose = False):
        """
        :Parameters:
         - `graph` (TPG)
         - `vids` (int|list) - list of (graph) vertex ids @ t_n. If :None: it will be computed for all possible vertex
         - `use_projected_anticlinal_wall` (bool) - if True, use the medians of projected anticlinal wall ('projected_anticlinal_wall_median') to compute a strain in 2D.
        """
        # Check if cells landmarks have been recovered ans stored in the graph structure:
        if use_projected_anticlinal_wall:
            assert 'surfacic_3D_landmarks' in graph.edge_property_names()
            assert 'epidermis_wall_median' in graph.vertex_property_names()
            assert 'daughters_fused_epidermis_wall_median' in graph.vertex_property_names()
            assert 'L1' in graph.vertex_property_names()
        else:
            assert '3D_landmarks' in graph.edge_property_names()
            assert 'epidermis_wall_median' in graph.vertex_property_names()
            assert 'unlabelled_wall_median' in graph.vertex_property_names()
            assert 'daughters_fused_epidermis_wall_median' in graph.vertex_property_names()
            assert 'daughters_fused_unlabelled_wall_median' in graph.vertex_property_names()
            #~ assert 'daughters_fused_wall_median' in graph.vertex_property_names()
            unlabelled_data = False
            try: graph.vertex_property('unlabelled_wall_median')
            except: unlabelled_data = True

        # -- If the vid is not associated through time it's not possible to compute the strain.
        if vids is None:
            vids = list(set(graph.lineaged_vertex(fully_lineaged=False))-set(graph.vertex_at_time(graph.nb_time_points,fully_lineaged=False)))
        if isinstance(vids,int):
            assert vids in graph.lineaged_vertex(fully_lineaged=False)
            vids = [vids]

        N = len(vids); percent=0
        missing_rank2_proj_mat, missing_daughters_fused_rank2_proj_mat, missing_epidermis_wall_median = [], [], []
        stretch_mat, score = {}, {}
        for n,vid in enumerate(vids):
            if verbose and n*100/float(N)>=percent: print "{}%...".format(percent),; percent += 10
            if verbose and n+1==N: print "100%"
            spatial_edges = list(graph.edges( vid, 's' ))
            # -- We recover the landmarks associated to each vertex:
            landmarks_t1, landmarks_t2 = [], []
            # - We use the epidermis wall median as an extra landmark if the vertex is in the L1:
            ep_wm = graph.vertex_property('epidermis_wall_median')
            daughters_fused_ep_wm = graph.vertex_property('daughters_fused_epidermis_wall_median')
            if graph.vertex_property('L1').has_key(vid):
                if ep_wm.has_key(vid) and daughters_fused_ep_wm.has_key(vid):
                    landmarks_t1.append(ep_wm[vid])
                    landmarks_t2.append(daughters_fused_ep_wm[vid])
                else:
                    missing_epidermis_wall_median.append(vid)
            if use_projected_anticlinal_wall:
                ppt = 'surfacic_3D_landmarks'
            else:
                ppt = '3D_landmarks'
                unlab_wm = graph.vertex_property('unlabelled_wall_median')
                fused_unlab_wm = graph.vertex_property('daughters_fused_unlabelled_wall_median')
                if unlabelled_data and unlab_wm.has_key(vid) and fused_unlab_wm.has_key(vid):
                    landmarks_t1.append(unlab_wm[vid])
                    landmarks_t2.append(fused_unlab_wm[vid])
 
            # - We create two matrix of landmarks positions (before and after deformation) to compute the strain:
            nb_missing_data = 0
            for eid in spatial_edges:
                # - If the spatial edge we are looking at is pointing to a target outside the L1, we do not want to use it for surfacic 3D:
                target = graph.edge_vertices(eid)[0] if graph.edge_vertices(eid)[0] != vid else graph.edge_vertices(eid)[1]
                if use_projected_anticlinal_wall and target not in graph.vertex_property('L1'):
                    continue # no need to worry, these are not the droids you're looking for !
                # - Now we add medians between used vertex as landmarks:
                if graph.edge_property(ppt).has_key(eid):
                    landmarks_t1.append(graph.edge_property(ppt)[eid][0])
                    landmarks_t2.append(graph.edge_property(ppt)[eid][1])
                else:
                    nb_missing_data+=1

            # - Make a projection into the rank-2 subspace:
            if use_projected_anticlinal_wall and graph._vertex_property.has_key('epidermis_rank-2_projection_matrix'):
                try:
                    H2_t1 = graph.vertex_property('epidermis_rank-2_projection_matrix')[vid]
                    landmarks_t1 = np.array([np.dot(H2_t1,pts) for pts in landmarks_t1])
                except:
                    missing_rank2_proj_mat.append(vid)
                try:
                    H2_t2 = graph.vertex_property('daughters_fused_epidermis_rank-2_projection_matrix')[vid]
                    landmarks_t2 = np.array([np.dot(H2_t2,pts) for pts in landmarks_t2])
                except:
                    missing_daughters_fused_rank2_proj_mat.append(vid)

            #~ if nb_missing_data != 0:
                #~ warnings.warn("Missing {} landmark{} for the t_n vertex {} at time {}".format(nb_missing_data, "s" if nb_missing_data>=2 else "", vid, graph.vertex_property('index')[vid]))
            if nb_missing_data == 0 and len(landmarks_t1)>=4: 
                assert len(landmarks_t1) == len(landmarks_t2)
                # - Convert voxel based metric into real-worlds units:
                vid_index = graph.vertex_property('index')[vid]
                res_1, res_2 = np.array(graph.graph_property("images_voxelsize")[vid_index]), np.array(graph.graph_property("images_voxelsize")[vid_index+1])
                landmarks_t1 = landmarks_t1*res_1
                landmarks_t2 = landmarks_t2*res_2
                stretch_mat[vid], score[vid] = func(graph, landmarks_t1, landmarks_t2)

        if missing_epidermis_wall_median != []:
            print 'Could not use the epidermis wall median as an extra landmark for vids: {}'.format(missing_epidermis_wall_median)
        if missing_rank2_proj_mat != []:
            print "Missing epidermis_rank-2_projection_matrix for vid: {}".format(vid)
        if missing_daughters_fused_rank2_proj_mat != []:
            print "Missing daughters_fused_epidermis_rank-2_projection_matrix for vid: {}".format(vid)

        # -- Now we return the results of temporal differentiation function:
        if labels_at_t_n:
            return stretch_mat, score
        else:
            return translate_keys2daughters_ids(graph, stretch_mat), translate_keys2daughters_ids(graph, score)

    return  wrapped_function

#~ @__strain_parameters
@__strain_parameters2
def stretch_matrix(graph, xyz_t1, xyz_t2):
    """
    Compute the stretch / deformation matrix.
    """
    from sklearn import linear_model
    # - Compute the centroids:
    c_t1 = np.mean(xyz_t1,0)
    c_t2 = np.mean(xyz_t2,0)
    # - Compute the centered matrix:
    centered_coord_t1=np.array(xyz_t1-c_t1)
    centered_coord_t2=np.array(xyz_t2-c_t2)
    # - A is the affine transformation matrix between the centered vertices position of two time points:
    regr = linear_model.Ridge (alpha = .01, fit_intercept=False)
    regr.fit(centered_coord_t1,centered_coord_t2)

    return regr.coef_, r2_scoring(centered_coord_t2, np.dot(regr.coef_,centered_coord_t1.T).T)

def stretch_main_orientations(graph, stretch_mat=None, **kwargs):
    """
    Return the stretch main directions and the associated values before deformation (default) or after.
    """
    try: vids = kwargs['vids']
    except: vids = None
    try: labels_at_t_n = kwargs['labels_at_t_n']
    except: labels_at_t_n = True
    try: use_projected_anticlinal_wall = kwargs['use_projected_anticlinal_wall']
    except: use_projected_anticlinal_wall = False
    if stretch_mat is None:
        print 'Computing the strecht matrix...'
        stretch_mat, score = stretch_matrix(graph, vids, True, use_projected_anticlinal_wall, True)

    directions = {}; values={}
    for vid in stretch_mat:
        ##  Singular Value Decomposition (SVD) of A.
        R,D_A,Q=svd(stretch_mat[vid])
        # Compute Strain Rates :
	# removed next 4 lines: NLM
        #if not labels_at_t_n:
        #    directions[vid] = return_list_of_vectors(Q, by_row=1)
        #else:
        #    directions[vid] = return_list_of_vectors(R, by_row=0)
        values[vid] = D_A

    if labels_at_t_n:
        return directions, values
    else:
        return translate_keys2daughters_ids(graph,directions), translate_keys2daughters_ids(graph,values)

def strain_rates(graph, stretch_mat=None, **kwargs):
    """
    Return the strain rate: sr[c][i] = np.log(D_A[i])/deltaT, for i = [0,1] if 2D, i = [0,1,2] if 3D.
    Dimensionnality is imposed by the one of the strain matrix.
    """
    try: vids = kwargs['vids']
    except: vids = None
    try: labels_at_t_n = kwargs['labels_at_t_n']
    except: labels_at_t_n = True
    try: use_projected_anticlinal_wall = kwargs['use_projected_anticlinal_wall']
    except: use_projected_anticlinal_wall = False
    if stretch_mat is None:
        print 'Computing the strecht matrix...'
        stretch_mat, score = stretch_matrix(graph, vids, True, use_projected_anticlinal_wall, True)

    sr = {}
    for vid in stretch_mat:
        ##  Singular Value Decomposition (SVD) of A.
        R,D_A,Q=svd(stretch_mat[vid])
        # Compute Strain Rates :
        sr[vid] = np.log(D_A)/float(time_interval(graph, vid, rank =1))

    if labels_at_t_n:
        return sr
    else:
        return translate_keys2daughters_ids(graph,sr)

def expansion_anisotropy(graph, stretch_mat=None, **kwargs):
    """
    Compute the expansion anisotropy in 2D: ea[c] = np.log(D_A[0]/D_A[1])/np.log(D_A[0]*D_A[1]).
    """
    try: vids = kwargs['vids']
    except: vids = None
    try: labels_at_t_n = kwargs['labels_at_t_n']
    except: labels_at_t_n = True
    try: use_projected_anticlinal_wall = kwargs['use_projected_anticlinal_wall']
    except: use_projected_anticlinal_wall = False
    if stretch_mat is None:
        print 'Computing the strecht matrix...'
        stretch_mat, score = stretch_matrix(graph, vids, True, use_projected_anticlinal_wall, True)

    ea = {}
    for vid in stretch_mat:
        ##  Singular Value Decomposition (SVD) of A.
        R,D_A,Q=svd(stretch_mat[vid])
        # Compute Strain Rates and Areal Strain Rate:
        ea[vid] = np.log(D_A[0]/D_A[1])/np.log(D_A[0]*D_A[1])

    if labels_at_t_n:
        return ea
    else:
        return translate_keys2daughters_ids(graph,ea)

def areal_strain_rates(graph, stretch_mat=None, **kwargs):
    """
    Compute the areal strain rate: asr[c] = sum_i(np.log(D_A[i])/deltaT), for i = [0,1].
    """
    try: vids = kwargs['vids']
    except: vids = None
    try: labels_at_t_n = kwargs['labels_at_t_n']
    except: labels_at_t_n = True
    try: use_projected_anticlinal_wall = kwargs['use_projected_anticlinal_wall']
    except: use_projected_anticlinal_wall = False
    if stretch_mat is None:
        print 'Computing the strecht matrix...'
        stretch_mat, score = stretch_matrix(graph, vids, True, use_projected_anticlinal_wall, True)

    asr = {}
    for vid in stretch_mat:
        ##  Singular Value Decomposition (SVD) of A.
        R,D_A,Q=svd(stretch_mat[vid])
        # Compute Strain Rates and Areal Strain Rate:
        asr[vid] = sum(np.log(D_A[0:2])/float(time_interval(graph, vid, rank =1)))

    if labels_at_t_n:
        return asr
    else:
        return translate_keys2daughters_ids(graph,asr)

def volumetric_strain_rates(graph, stretch_mat=None, **kwargs):
    """
    Compute the volumetric strain rate: asr[c] = sum_i(np.log(D_A[i])/deltaT), for i = [0,1,2].
    """
    try: vids = kwargs['vids']
    except: vids = None
    try: labels_at_t_n = kwargs['labels_at_t_n']
    except: labels_at_t_n = True
    try: use_projected_anticlinal_wall = kwargs['use_projected_anticlinal_wall']
    except: use_projected_anticlinal_wall = False
    if stretch_mat is None:
        print 'Computing the strecht matrix...'
        stretch_mat, score = stretch_matrix(graph, vids, True, use_projected_anticlinal_wall, True)

    vsr = {}
    for vid in stretch_mat:
        ##  Singular Value Decomposition (SVD) of A.
        R,D_A,Q=svd(stretch_mat[vid])
        # Compute Strain Rates and Areal Strain Rate:
        vsr[vid] = sum(np.log(D_A)/float(time_interval(graph, vid, rank =1)))

    if labels_at_t_n:
        return vsr
    else:
        return translate_keys2daughters_ids(graph,vsr)

def anisotropy_ratios(graph, stretch_mat=None, **kwargs):
    """
    Anisotropy ratio :
     R, strain_values, Q = svd( stretch_mat ), then
    -if 2D:
        return strain_values[0]/strain_values[1]
    - if 3D:
        return [ sv[0]/sv[1], sv[1]/sv[2], sv[0]/sv[2] ]

    Dimensionnality is imposed by the one of the strain matrix `stretch_mat`.
    """
    try: vids = kwargs['vids']
    except: vids = None
    try: labels_at_t_n = kwargs['labels_at_t_n']
    except: labels_at_t_n = True
    try: use_projected_anticlinal_wall = kwargs['use_projected_anticlinal_wall']
    except: use_projected_anticlinal_wall = False
    if stretch_mat is None:
        print 'Computing the strecht matrix...'
        stretch_mat, score = stretch_matrix(graph, vids, True, use_projected_anticlinal_wall, True)

    anisotropy_ratio = {}
    for vid in stretch_mat:
        ##  Singular Value Decomposition (SVD) of A.
        R,D_A,Q = svd(stretch_mat[vid])
        if len(D_A) == 3:
            anisotropy_ratio[vid] = [ D_A[0]/D_A[1], D_A[1]/D_A[2], D_A[0]/D_A[2] ]
        else:
            anisotropy_ratio[vid] = D_A[0]/D_A[1]

    if labels_at_t_n:
        return anisotropy_ratio
    else:
        return translate_keys2daughters_ids(graph,anisotropy_ratio)


def time_interval(graph,vid,rank=1):
    """
    Compute the time interval for the vexterx id `vid` thanks to data saved in graph.graph_property('time_steps').
    """
    index_1 = graph.vertex_property('index')[vid]
    return (graph.graph_property('time_steps')[index_1+rank]-graph.graph_property('time_steps')[index_1])


def sibling_volume_ratio(graph):
    """
    """
    svr={}
    used_vtx = []
    for vtx in graph.vertices():
        sibling = graph.sibling(vtx)
        if sibling is not None and len(sibling)==1 and list(graph.sibling(vtx))[0] not in used_vtx:
            used_vtx.append(vtx)
            sibling = list(graph.sibling(vtx))[0]
            ratio = graph.vertex_property('volume')[vtx]/graph.vertex_property('volume')[sibling]
            if ratio >=1:
                svr[sibling,vtx]=1./ratio
            else:
                svr[vtx,sibling]=ratio
    return svr


def subspace_projection(point_set, subspace_rank = 2, centering=True, verbose=True):
    """
    Project a point set o coordinate into a subspace.
    
    :Parameters:
     - point_set (np.array): list of coordinates of shape (n_point, init_dim).
     - dimension_reduction (int) : the dimension reduction to apply
    """
    point_set = np.array(point_set)
    nb_coord = point_set.shape[0]
    init_dim = point_set.shape[1]
    assert init_dim > subspace_rank

    if centering:
        # - Compute the centered matrix:
        centered_point_set=point_set-point_set.mean(axis=0)
    else:
        centered_point_set = point_set
        if point_set.mean(axis=0) != np.zeros([init_dim, init_dim]):
            warnings.warn("The provided point set is not centered!")
    
    # -- Compute the Singular Value Decomposition (SVD) of centered coordinates:
    U,D,V = svd(centered_point_set, full_matrices=False)
    V = V.T

    # -- Compute the projection matrix:
    H = np.dot(V[:,0:subspace_rank], V[:,0:subspace_rank].T)
    # -- Projection of the coordinate into the defined subspace by the previously computed eigenvector:
    reduced_point_set = np.array([np.dot(H, centered_point_set[k]) for k in range(nb_coord)])
    assert reduced_point_set.shape[1] == init_dim

    return reduced_point_set


def r2_scoring( Y, Y_pred ):
    """
    Function computing an r^2-like scoring in 3D.
    """
    Y = np.array(Y); Y_pred = np.array(Y_pred)
    assert Y.shape == Y_pred.shape
    if Y.shape[0]==3 and Y.shape[0]>3:
        Y = Y.T; Y_pred = Y_pred.T
    n_coord = Y.shape[0]

    from sa_oa.image.algo.analysis import distance
    SCR = sum( [distance( Y[n], Y_pred[n] )**2 for n in range(n_coord)] )
    SCT = sum( [distance( Y[n], np.mean(Y,0) )**2 for n in range(n_coord)] )
    return 1-SCR/SCT


def triplot(graphs_list, values2plot, labels_list=None, values_name="", normed=False):
    """
    TO DO
    """
    import numpy as np
    if labels_list==None:
        labels_list=[]
        for g in graphs_list:
            labels_list.append(g.vertex_property('label'))

    values=[]
    abs_dev_values=[]
    laplacian_values=[]
    #-- if 'values2plot' is a string, it must be a property found in all graphs in the 'graph_list'.
    if type(values2plot)==type(str('str')):
        for g in graphs_list:
            if values2plot not in g.vertex_property_names():
                #sys.exit(1)
                raise ValueError(values2plot)
            else:
                if (values_name==""):
                    values_name=values2plot
                values.append(g.vertex_property(values2plot).values())
                abs_dev_values.append(dev_abs(g,values2plot,True))
                laplacian_values.append(laplacian(g,values2plot,True))

    import matplotlib.pyplot as plt
    fig = plt.figure()
    fig.subplots_adjust( wspace=0.13, left=0.05, right=0.95, top=0.95)
    main=fig.add_subplot(1,2,1)
    main.hist(values, bins=20,normed=normed,
        label=( ('t1, n='+str(len(values[0]))+', mean='+str(np.round(np.mean(values[0]), 2))) ,
        ('t2, n='+str(len(values[1]))+', mean='+str(np.round(np.mean(values[1]), 2))) ,
        ('t3, n='+str(len(values[2]))+', mean='+str(np.round(np.mean(values[2]), 2))) ), histtype='bar' )
    plt.title("L1 cells' "+values_name)
    if values_name=='volume':
        plt.xlabel('Volumes'+ r' ($\mu m^3$)')
    else:
        plt.xlabel(values_name)
    if normed:
            plt.ylabel('Frequency')
    else:
        plt.ylabel('Number of observations')
    plt.legend()

    dev=fig.add_subplot(2,2,2)
    dev.hist(abs_dev_values, bins=20,normed=normed,
        label=( ('t1, n='+str(len(abs_dev_values[0]))+', mean='+str(np.round(np.mean(abs_dev_values[0]), 2))) ,
        ('t2, n='+str(len(abs_dev_values[1]))+', mean='+str(np.round(np.mean(abs_dev_values[1]), 2))) ,
        ('t3, n='+str(len(abs_dev_values[2]))+', mean='+str(np.round(np.mean(abs_dev_values[2]), 2))) ), histtype='bar' )
    plt.title("L1 cells' absolute deviance from neighbors in "+values_name)
    plt.xlabel('Deviance from neighbors in volumes'+ r' ($\mu m^3$)')
    if normed:
            plt.ylabel('Frequency')
    else:
        plt.ylabel('Number of observations')
    plt.legend()

    lap=fig.add_subplot(2,2,4)
    lap.hist(laplacian_values, bins=20,normed=normed,
        label=( ('t1, n='+str(len(laplacian_values[0]))+', mean='+str(np.round(np.mean(laplacian_values[0]), 2))) ,
        ('t2, n='+str(len(laplacian_values[1]))+', mean='+str(np.round(np.mean(laplacian_values[1]), 2))) ,
        ('t3, n='+str(len(laplacian_values[2]))+', mean='+str(np.round(np.mean(laplacian_values[2]), 2))) ), histtype='bar' )
    plt.title("L1 cells' laplacian from neighbors in "+values_name)
    plt.xlabel('Laplacian from neighbors in volumes'+ r' ($\mu m^3$)')
    if normed:
            plt.ylabel('Frequency')
    else:
        plt.ylabel('Number of observations')
    plt.legend()
    plt.show()


#~ def shape_anisotropy_2D(graph, vids=None, add2vertex_property = True):
    #~ """
    #~ Compute shape anisotropy in 2D based on the two largest inertia axis length.
    #~ !!!!! SHOULD BE COMPARED WITH NORMAL VECTOR FROM CURVATURE COMPUTATION !!!!!
    #~ """
    #~ if vids == None:
        #~ vids = graph.vertex_property('inertia_axis').keys()
    #~ else:
        #~ for vid in vids:
            #~ if not graph.vertex_property('inertia_axis').has_key(vid):
                #~ warnings.warn("Inertia axis hasn't been computed for vid #"+str(vid))
                #~ vids.remove(vid)
    #~
    #~ shape_anisotropy_2D = {}
    #~ for vid in vids:
        #~ axis_len = graph.vertex_property('inertia_axis')[vid][1]
        #~ shape_anisotropy_2D[vid] = float(axis_len[0]-axis_len[1])/(axis_len[0]+axis_len[1])
    #~
    #~ if add2vertex_property:
        #~ graph.add_vertex_property("2D shape anisotropy",shape_anisotropy_2D)
    #~
    #~ return shape_anisotropy_2D
#~
#~
