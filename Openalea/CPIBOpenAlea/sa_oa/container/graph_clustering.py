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
"""This module helps to use clustering and standardization methods on graphs."""

import warnings
import numpy as np
from numpy import ndarray
import copy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d.axes3d import Axes3D
from scipy.sparse import csr_matrix
from sa_oa.container.temporal_graph_analysis import exist_relative_at_rank

def distance_matrix_from_vector(data, variable_types, no_dist_index = []):
    """
    Function creating a distance matrix based on a vector (list) of values.
    Each values are attached to an individual.

    :Parameters:
     - `data` (list) - vector/list of value
     - `variable_types` (str) - type of variable

    :Returns:
     - `dist_mat` (np.array) - distance matrix
    """
    N = len(data)
    dist_mat = np.zeros( shape = [N,N], dtype=float )

    if variable_types == "Numeric":
        for i in xrange(N):
            for j in xrange(i+1,N):# we skip when i=j because in that case the distance is 0.
                if i in no_dist_index or j in no_dist_index or data[i] is None or data[j] is None:
                    dist_mat[i,j] = dist_mat[j,i] = None
                else:
                    dist_mat[i,j] = dist_mat[j,i] = abs(data[i]-data[j])

    if variable_types == "Ordinal":
        rank = data.argsort() # In case of ordinal variables, observed values are replaced by ranked values.
        for i in xrange(N):
            for j in xrange(i+1,N): # we skip when i=j because in that case the distance is 0.
                if i in no_dist_index or j in no_dist_index or data[i] is None or data[j] is None:
                    dist_mat[i,j] = dist_mat[j,i] = None
                else:
                    dist_mat[i,j] = dist_mat[j,i] = abs(rank[i]-rank[j])

    return dist_mat


def standardisation(data, norm, variable_types = None, verbose = False):
    """
    :Parameters:
     - `mat` (np.array) - distance matrix
     - `norm` (str) - "L1" or "L2", select which standarisation metric to apply to the data;
     - `variable_types` (str) - "Numeric" or "Ordinal" or "Interval"
    
    :Returns:
     - `standard_mat` (np.array) - standardized distance matrix
    """
    
    if (norm != 'L1') and (norm != 'L2'):
        raise ValueError("Undefined standardisation metric")

    # -- Identifying case where numpy.array are vectors:
    if isinstance(data,ndarray) and (((data.shape[0] == 1) and (data.shape[1] > 1)) or ((data.shape[1] == 1) and (data.shape[0] > 1))):
        data.tolist()

    # -- Creating the distance matrix if not provided in `data`:
    if isinstance(data,list):
        if isinstance(data[0],list) and len(data)==len(data[0]):
            data = np.array(data)
        elif not isinstance(data[0],list):
            if verbose: print "You provided a vector of variable, we compute the distance matrix first !"
            distance_matrix = distance_matrix_from_vector(data, variable_types)
        else:
            raise ValueError("Can not convert the provided data.")

    if isinstance(data,ndarray) and (data.shape[0]==data.shape[1]) and ((data.shape[0]!=1)and(data.shape[1]!=1)):
        if verbose: print "You provided a distance matrix !"
        distance_matrix = data

    # -- Now we can start the standardisation:
    nan_index = np.isnan(distance_matrix)
    if True in nan_index:
        nb_missing_values = len(np.where(nan_index is True)[0])
        distance_matrix = np.nan_to_num(distance_matrix)
        #~ warnings.warn("It seems there are {0} missing values to take into account!".format(np.sqrt(nb_missing_values)))
    else:
        nb_missing_values = 0.

    N = distance_matrix.shape[0]
    if norm == "L1":
        absd = np.nansum(np.nansum(abs(distance_matrix))) / (N*(N-1)-nb_missing_values)
        return distance_matrix / absd

    if norm == "L2":
        sd = np.nansum(np.nansum(distance_matrix**2)) / (N*(N-1)-nb_missing_values)
        return distance_matrix / sd


def standardisation_scikit(data, norm, variable_types):
    """
    :Parameters:
     - `mat` (np.array) - distance matrix
     - `norm` (str) - "L1" or "L2", select which standarisation metric to apply to the data;
     - `variable_types` (str) - "Numeric" or "Ordinal" or "Interval"
    
    :Returns:
     - `standard_mat` (np.array) - standardized distance matrix
    """
    if (norm != 'L1') and (norm != 'L2'):
        raise ValueError("Undefined standardisation metric")

    X = distance_matrix_from_vector(data, variable_types)

    from sklearn import preprocessing

    if norm == "L1":
        return preprocessing.normalize(X, norm='l1')

    if norm == "L2":
        return preprocessing.normalize(X, norm='l2')


def weighted_distance_matrix(distance_matrix_list, weights_list):
    """
    Compute a distance matrix from a list of distance matrix and a list of weight:
        D_w = sum_i(w_i * D_i); where sum_i(w_i) = 1.
    """
    if sum(weights_list) != 1:
        raise ValueError("The weights do not sum to 1 !")

    if len(distance_matrix_list) != len(weights_list):
        raise ValueError("You did not provided the same number of distance matrix and weigths !")

    # weighted_matrix initialisation:
    weighted_matrix = distance_matrix_list[0].copy()
    weighted_matrix.fill(0)
    for n, weigth in enumerate(weights_list):
        weighted_matrix += weigth * distance_matrix_list[n]

    return weighted_matrix


def csr_matrix_from_graph(graph, vids2keep):
    """
    Create a sparse matrix representing a connectivity matrix recording the topological information of the graph.
    Defines for each vertex the neighbouring vertex following a given structure of the data.

    :Parameters:
     - `graph` (Graph | PropertyGraph | TemporalPropertyGraph) - graph from which to extract connectivity.
     - `vids2keep` (list) - list of vertex ids to build the sparse matrix from (columns and rows will be ordered according to this list).
    """    
    N = len(vids2keep)
    data,row,col = [],[],[]

    for edge in graph.edges(edge_type='s'):
        s,t = graph.edge_vertices(edge)
        if (s in vids2keep) and (t in vids2keep):
            row.extend([vids2keep.index(s),vids2keep.index(t)])
            col.extend([vids2keep.index(t),vids2keep.index(s)])
            data.extend([1,1])

    return csr_matrix((data,(row,col)), shape=(N,N))


def weighted_global_distance_matrix(graph, variable_list=[], variable_weights=[], variable_types=[], topo_weight=0., temporal_list=[], temporal_weights=[], temporal_types=[], standardisation_method = "L1", only_lineaged_vertices = True, vids = None, rank = 1 ):
    """
    :Parameters:
     - `graph` (Graph | PropertyGraph | TemporalPropertyGraph) - graph from which to extract spatial, spatio-temporal and topological/euclidean distance variables.
     - `variable_list` (list) - list of `vertex_property_names` related to spatial information (ex. volume).
     - `variable_weights` (list) - list of weights related to spatial information (ex. volume). If only an integer is given for several variables (in `variable_list`), we divide it by the number of variables in `variable_list`.
     - `variable_types` (list) - list of variable types. Can be "Ordinal" or "Numeric".
     - `topo_weight` (int) - weight related to topological/euclidian distance.
     - `temporal_list` (list) - list of `vertex_property_names` related to spatio-temporal information (ex. volumetric growth).
     - `temporal_weights` (list) - list of weights related to spatio-temporal information (ex. volumetric growth). If only an integer is given for several variables (in `temporal_list`), we divide it by the number of variables in `temporal_list`.

    :NOTE:
     - We use :ABSOLUTE: temporal change (from sa_oa.container.temporal_graph_analysis import temporal_change)!
     - We assign temporal change to t_n+1.
     - One can provide a vector (type list) of spatial or/and spatio-temporal data => NEED to be detected before creating the `vtx_list` and would force to provide this list of vertex!!! (or we could use a dictionary)  NOT DONE YET !!!
    """
    # -- Taking care of "stupid" cases:
    if topo_weight == []: topo_weight = 0.
    if variable_weights == 0.: variable_weights = []
    if temporal_weights == 0.: temporal_weights = []
    if float(topo_weight) == 1.:
        warnings.warn("No topological/euclidean distance between each time point !!")
    if isinstance(topo_weight,list) and len(topo_weight) == 1:
        topo_weight = float(topo_weight[0])
    elif not isinstance(topo_weight,float):
        raise ValueError("Check your value for `topo_weight`, should be an float in [0.,1.].")

    # -- Making sure all spatial variable asked for in `variable_list` are present in the `graph`:
    if isinstance(variable_list,str):
        assert variable_list in graph.vertex_property_names()
    if isinstance(variable_list,list):
        for variable_name in variable_list:
            if isinstance(variable_name,str):
                assert variable_name in graph.vertex_property_names()

    # -- Handling multiple types of inputs for vertex related variables (usually spatial properties):
    if isinstance(variable_list,str) or isinstance(variable_list,dict): variable_list = [variable_list]
    if isinstance(variable_types,str): variable_types = [variable_types]
    if isinstance(variable_weights,int) or isinstance(variable_weights,float): 
        variable_weights = [float(variable_weights)]
    nb_variables = len(variable_list)
    if nb_variables != 0:
        assert len(variable_list)==len(variable_types)
        # - If only an integer is given for several variables (in `variable_list`), we divide it by the number of variables in `variable_list`:
        if len(variable_weights)==1 and len(variable_weights) != nb_variables:
            variable_weights = [variable_weights[0] / float(nb_variables) for i in xrange(nb_variables)]
        if len(variable_weights) != nb_variables:
            raise ValueError("The list `variable_weights` and `variable_list` should be of the same length or len(variable_weights)==1")

    # -- Handling multiple types of inputs for spatio-temporal weight, variables and types:
    if isinstance(temporal_list,str) or isinstance(temporal_list,dict): temporal_list = [temporal_list]
    if isinstance(temporal_types,str): temporal_types = [temporal_types]
    nb_temporal_variables = len(temporal_list)
    if isinstance(temporal_weights,int) or isinstance(temporal_weights,float):
        temporal_weights = [float(temporal_weights)]
    if nb_temporal_variables != 0:
        assert len(temporal_list)==len(temporal_types)
        # - If only an integer is given for several temporal variables (in `temporal_list`), we divide it by the number of temporal variable in `temporal_list`:
        if len(temporal_weights)==1 and len(temporal_weights) != nb_temporal_variables:
            temporal_weights = [temporal_weights[0] / float(nb_temporal_variables) for i in xrange(nb_temporal_variables)]
        if len(temporal_weights) != nb_temporal_variables:
            raise ValueError("The list `temporal_weights` and `temporal_list` should be of the same length or len(temporal_weights)==1")

    assert sum(variable_weights)+sum(temporal_weights)+topo_weight==1.

    # -- Creating the list of vertices:
    # -- We keep vertex ids only if they are temporally linked in the graph at `rank` or -`rank`
    if vids is None:
        vtx_list = [vid for vid in graph.vertices() if exist_relative_at_rank(graph, vid, rank) or exist_relative_at_rank(graph, vid, -rank)]
    else:
        vtx_list = [vid for vid in vids if exist_relative_at_rank(graph, vid, rank) or exist_relative_at_rank(graph, vid, -rank)]
    N = len(vtx_list)
    index = dict( (vid,graph.vertex_property('index')[vid]) for vid in vtx_list )

    # -- We compute the standardized distance matrix related to spatial variables:
    if nb_variables != 0:
        print("Computing the standardized distance matrix related to spatial variables...")
        variable_standard_distance_matrix = {}
        for n, variable_name in enumerate(variable_list):
            if isinstance(variable_name,str):
                variable_vector = [graph.vertex_property(variable_name)[vid] if graph.vertex_property(variable_name).has_key(vid) else None  for vid in vtx_list]# we need to do that if we want to have all matrix ordered the same way
            if isinstance(variable_name,dict):
                variable_vector = [variable_name[vid] if variable_name.has_key(vid) else None for vid in vtx_list]# we need to do that if we want to have all matrix ordered the same way
            variable_standard_distance_matrix[n] = standardisation(variable_vector, standardisation_method, variable_types[n])
    if nb_variables != 0 and variable_weights[0] == 1.: # there is no data to 're-norm'...
        return vtx_list, variable_standard_distance_matrix[0]

    # -- We compute the standardized distance matrix related to temporal variables:
    if nb_temporal_variables != 0:
        print("Computing the standardized distance matrix related to temporal variables...")
        # If we want to work with temporally differentiated variables, we will have to filter vertex without parent (since we assign spatio-temporal variable @ t_n+1):
        temporal_standard_distance_matrix = {}
        for n, temporal_name in enumerate(temporal_list):
            if isinstance(temporal_name,dict):
                dict_temporal = temporal_name
            elif isinstance(temporal_name,str):
                from sa_oa.container.temporal_graph_analysis import temporal_change
                dict_temporal = temporal_change(graph, temporal_name, vtx_list, rank, labels_at_t_n = False)
            else:
                raise ValueError("Unrecognized type of data.")
            # - Now we create the 'vector' of data sorted by vertices to create the distance matrix:
            temporal_distance_list = [dict_temporal[vid] if dict_temporal.has_key(vid) else None for vid in vtx_list]# we need to do that if we want to have all matrix ordered the same way
            temporal_standard_distance_matrix[n] = standardisation(temporal_distance_list, standardisation_method, temporal_types[n])

    if sum(temporal_weights) == 1.:
        warnings.warn("You have asked only for a pairwise distance matrix based on temporally differentiated variables affected @t_n+1. There will be no distances for the first time-point!")
    if nb_temporal_variables != 0 and temporal_weights[0] == 1.: # there is no data to 're-norm'...
        return vtx_list, temporal_standard_distance_matrix[0]

    # -- We compute the standardized topological distance matrix:
    if topo_weight != 0.:
        print("Computing the standardized topological distance matrix...")
        import time
        t = time.time()
        # - Extraction of the topological distance between all vertex of a time point:
        topo_dist = dict( [(vid, graph.topological_distance(vid, 's',return_inf=False)) for vid in vtx_list] )
        print "Time to compute the topological distances: {0}s".format(time.time() - t)
        # - Transformation into a distance matrix
        mat_topo_dist_rank = np.array([[(topo_dist[i][j] if i!=j else 0) for i in vtx_list] for j in vtx_list])

        # - Topological distance standardisation
        mat_topo_dist_rank_standard = standardisation(mat_topo_dist_rank, "L1")
    else:
        mat_topo_dist_rank_standard = np.zeros( [N,N], dtype=int ) # `mat_topo_dist_rank_standard``should exist, but the values will not be used !!
    if topo_weight == 1.:
        return vtx_list, mat_topo_dist_rank_standard
        #~ return vtx_list, mat_topo_dist_rank

    def renorm(line, column, mat_topo, var_mat, temp_mat, w_topo, w_var, w_temp):
        w_renorm_topo = 0.
        if w_topo != 0.:
            if np.isnan(mat_topo[i,j]):
                w_renorm_topo = w_topo
                w_topo = 0.

        w_renorm_var = 0.
        for n in var_mat:
            if np.isnan(var_mat[n][i,j]):
                w_renorm_var += w_var[n]
                w_var[n] = 0.

        w_renorm_temp = 0.
        for n in temp_mat:
            if np.isnan(temp_mat[n][i,j]):
                w_renorm_temp += w_temp[n]
                w_temp[n] = 0.

        renorm = (1.-(w_renorm_topo+w_renorm_var+w_renorm_temp))
        if renorm != 0.:
            return w_topo/renorm, np.array(w_var)/renorm if w_var!=[] else [], np.array(w_temp)/renorm if w_temp!=[] else []
        else:
            return w_topo, w_var, w_temp

    print("Creating the global weighted pairwise standard distance matrix...")
    # - Replacing nan by zeros for computation.
    mat_topo = np.nan_to_num(mat_topo_dist_rank_standard)
    if nb_variables != 0.:
        var_mat = [np.nan_to_num(variable_standard_distance_matrix[n]) for n in xrange(len(variable_standard_distance_matrix))]
    else:
        var_mat, variable_standard_distance_matrix = [], []
    if nb_temporal_variables != 0.:
        temp_mat = [np.nan_to_num(temporal_standard_distance_matrix[n]) for n in xrange(len(temporal_standard_distance_matrix))]
    else:
        temp_mat, temporal_standard_distance_matrix = [], []

    # Finally making the global weighted pairwise standard distance matrix:
    D = np.zeros( shape = [N,N], dtype=float )
    w_mat = np.zeros( shape = [N,N], dtype=float )
    for i in xrange(D.shape[0]):
        for j in xrange(D.shape[1]):
            if i>j: #D[i,j]=D[j,i] and if i==j, D[i,j]=D[j,i]=0
                # - Computing weight according to missing values.
                w_topo, w_var, w_temp = renorm(i,j,mat_topo_dist_rank_standard, variable_standard_distance_matrix, temporal_standard_distance_matrix, copy.copy(topo_weight), copy.copy(variable_weights), copy.copy(temporal_weights))
                # - Pairwise weighted standard distance matrix
                D[i,j] = D[j,i] = w_topo * mat_topo[i,j] + sum([w_var[n]*var_mat[n][i,j] for n in xrange(len(w_var))]) + sum([w_temp[n]*temp_mat[n][i,j] for n in xrange(len(w_temp))])
    
    return vtx_list, D


def cluster_distance_matrix( distance_matrix, clustering ):
    """
    Function computing distance between clusters.
    For $\ell \eq q$  :
    \[ D(q,\ell) = \dfrac{ \sum_{i,j \in q; i \neq j} D(i,j) }{(N_{q}-1)N_{q}} , \]
    For $\ell \neq q$  :
    \[ D(q,\ell) = \dfrac{ \sum_{i \in q} \sum_{j \in \ell} D(i,j) }{N_{q} N_{\ell} } , \]
    where $D(i,j)$ is the distance matrix, $N_{q}$ and $N_{\ell}$ are the number of elements found in clusters $q$ and $\ell$.

    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    clusters_ids = list(set(clustering))
    nb_clusters = len(clusters_ids)
    nb_ids_by_clusters = [len(np.where(clustering == q)[0]) for q in clusters_ids]

    D = np.zeros( shape = [nb_clusters,nb_clusters], dtype = float )
    for n,q in enumerate(clusters_ids):
        for m,l in enumerate(clusters_ids):
            if n==m:
                index_q = np.where(clustering == q)[0]
                D[n,m] = sum( [distance_matrix[i,j] for i in index_q for j in index_q if i!=j] ) / ( (nb_ids_by_clusters[n]-1) * nb_ids_by_clusters[n])
            if n>m:
                index_q = np.where(clustering == q)[0]
                index_l = np.where(clustering == l)[0]
                D[n,m] = D[m,n]= sum( [distance_matrix[i,j] for i in index_q for j in index_l] ) / (nb_ids_by_clusters[n] * nb_ids_by_clusters[m])

    return D


def within_cluster_distance(distance_matrix, clustering):
    """
    Function computing within cluster distance.
    $$ D(q) = \dfrac{ \sum_{i,j \in q; i \neq j} D(i,j) }{(N_{q}-1)N_{q}} ,$$
    where $D(i,j)$ is the distance matrix, $N$ is the total number of elements and $N_{q}$ is the number of elements found in clusters $q$.

    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    clusters_ids = list(set(clustering))
    nb_clusters = len(clusters_ids)
    nb_ids_by_clusters = [len(np.where(clustering == q)[0]) for q in clusters_ids]

    D_within = {}
    for n,q in enumerate(clusters_ids):
        index_q = np.where(clustering == q)[0]
        D_within[q] = 2. * sum( [distance_matrix[i,j] for i in index_q for j in index_q if j>i] ) / ( (nb_ids_by_clusters[n]-1) * nb_ids_by_clusters[n])

    if nb_clusters == 1:
        return D_within.values()
    else:
        return D_within


def between_cluster_distance(distance_matrix, clustering):
    """
    Function computing within cluster distance.
    $$ D(q) = \dfrac{ \sum_{i \in q} \sum_{j \not\in q} D(i,j) }{ (N - N_q) N_q }, $$
    where $D(i,j)$ is the distance matrix, $N$ is the total number of elements and $N_{q}$ is the number of elements found in clusters $q$.

    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    N = len(clustering)
    clusters_ids = list(set(clustering))
    nb_clusters = len(clusters_ids)
    nb_ids_by_clusters = [len(np.where(clustering == q)[0]) for q in clusters_ids]

    if 1 in nb_ids_by_clusters:
        raise ValueError("A cluster contain only one element!")

    D_between = {}
    for n,q in enumerate(clusters_ids):
        index_q = np.where(clustering == q)[0]
        index_not_q = list(set(xrange(len(clustering)))-set(index_q))
        D_between[q] = sum( [distance_matrix[i,j] for i in index_q for j in index_not_q] ) / ( (N-nb_ids_by_clusters[n]) * nb_ids_by_clusters[n])

    if nb_clusters == 1:
        return D_between.values()
    else:
        return D_between


def cluster_diameters(distance_matrix, clustering):
    """
    Function computing within cluster diameter, i.e. the max distance between two vertex from the same cluster.
    $$ \max_{i,j \in q} D(j,i) ,$$
    where $D(i,j)$ is the distance matrix and $q$ a cluster, .
    
    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    clusters_ids = list(set(clustering))
    nb_clusters = len(clusters_ids)

    diameters = {}
    for q in clusters_ids:
        index_q = np.where(clustering == q)[0]
        diameters[q] = max([distance_matrix[i,j] for i in index_q for j in index_q])

    if nb_clusters == 1:
        return diameters.values()
    else:
        return diameters


def clusters_separation(distance_matrix, clustering):
    """
    Function computing within cluster diameter, i.e. the min distance between two vertex from two diferent clusters.
    $$ \min_{i \in q, j \not\in q} D(j,i) ,$$
    where $D(i,j)$ is the distance matrix and $q$ a cluster.
    
    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    clusters_ids = list(set(clustering))
    nb_clusters = len(clusters_ids)

    separation = {}
    for q in clusters_ids:
        index_q = np.where(clustering == q)[0]
        index_not_q = np.where(clustering != q)[0]
        separation[q] = min([distance_matrix[i,j] for i in index_q for j in index_not_q ])

    if nb_clusters == 1:
        return separation.values()
    else:
        return separation


def global_cluster_distance(distance_matrix, clustering):
    """
    Function computing global cluster distances, i.e. return the sum of within_cluster_distance and between_cluster_distance.

    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    w = within_cluster_distance(distance_matrix, clustering)
    b = between_cluster_distance(distance_matrix, clustering)

    N = len(clustering)
    clusters_ids = list(set(clustering))
    nb_ids_by_clusters = [len(np.where(clustering == q)[0]) for q in clusters_ids]

    gcd_w = sum( [ (nb_ids_by_clusters[q]*(nb_ids_by_clusters[q]-1))/float(sum([nb_ids_by_clusters[l]*(nb_ids_by_clusters[l]-1) for l in clusters_ids if l != q])) * w[q] for q in clusters_ids] )
    gcd_b = sum( [(N-nb_ids_by_clusters[q])*nb_ids_by_clusters[q]/float(sum([(N-nb_ids_by_clusters[l])*nb_ids_by_clusters[l] for l in clusters_ids if l != q])) * b[q] for q in clusters_ids] )
    return gcd_w, gcd_b


def CH_estimator(k,w,b,N):
    """
    Index of Calinski and Harabasz (1974).
    $$ \text{CH}(k) = \dfrac{B(k) / (k-1)}{W(k) / (n-k)} $$
    where $B(k)$ and $W(k)$ are the between- and within-cluster sums of squares, with $k$ clusters.
    The idea is to maximize $\text{CH}(k)$ over the number of clusters $k$. $\text{CH}(1)$ is not defined.

    :Parameters:
     - k: number of clusters;
     - w: global WITHIN cluster distances;
     - b: global BETWEEN cluster distances;
     - N: population size.
    """
    return (b/(k-1))/(w/(N-k))

def Hartigan_estimator(k,w_k,N):
    """
    Index of Hartigan (1975).
    Hartigan (1975) proposed the statistic:
    $$ \text{H}(k) = \left\{ \dfrac{W(k)}{W(k+1)} - 1 \right\} / (n-k-1)$$
    The idea is to start with $k=1$ and to add a cluster as long as $\text{H}(k)$ is sufficiently large.\\
    One can use an approximate $F$-distribution cut-off; instead Hartigan suggested that a cluster be added if $H(k) > 10$. Hence the estimated number of clusters is the smallest $k \geqslant 1$ such that $\text{H}(k) \leqslant 10$. This estimate is defined for $k=1$ and can potentially discriminate between one \textit{versus} more than one cluster.

    :Parameters:
     - k: number of clusters;
     - w_k: DICTIONARY of global WITHIN cluster distances (must contain k and k+1);
     - N: population size.
    """
    return (w_k[k]/w_k[k+1]-1)/(N-k-1)


def clustering_estimators(distance_matrix, clustering_method, clustering_range=[1,10]):
    """
    Compute various estimators based on clustering results.
    
    :Parameters:
     - distance_matrix (np.array): distance matrix to be used for clustering;
     - clustering_method (str): clustering methods to be applyed, must be "Ward" or "Spectral"
     - clustering_range (list): range of clusters to consider.
    """
    from sklearn import metrics
    from sklearn.cluster import spectral_clustering, Ward
    clustering, w, N = {}, {}, {}
    CH, sil = {}, {}

    for k in xrange(clustering_range[0], clustering_range[1]+1):
        if clustering_method == "Ward":
            ward = Ward(n_clusters=k).fit(distance_matrix)
            clustering[k] = ward.labels_
        if clustering_method == "Spectral":
            beta = 1
            similarity = np.exp(-beta * distance_matrix / distance_matrix.std())
            clustering[k] = spectral_clustering( similarity, n_clusters=k )

        N[k] = len(clustering[k])
        if k!=1:
            w[k], b = global_cluster_distance(distance_matrix, clustering[k])
            CH[k] = CH_estimator(k,w[k],b,N[k])
            sil[k] = metrics.silhouette_score(distance_matrix, clustering[k], metric='euclidean')
        else:
            w[k] = within_cluster_distance(distance_matrix, clustering[k])

    Hartigan = dict( [(k, Hartigan_estimator(k,w,N[k]) ) for k in xrange(clustering_range[0], clustering_range[1])] )

    return clustering, CH, Hartigan, sil


def vertex2clusters_distance(distance_matrix, clustering):
    """
    Compute the mean distance between a vertex and those from each group.
    $$ D(i,q) = \dfrac{ \sum_{i \neq j} D(i,j) }{ N_q } \: , \quad \forall i \in [1,N] \:, \: q \in [1,Q],$$
    where $D(i,j)$ is the distance matrix and $q$ a cluster.
    
    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    N = len(clustering)
    clusters_ids = list(set(clustering))
    nb_clusters = len(clusters_ids)
    nb_ids_by_clusters = [len(np.where(clustering == q)[0]) for q in clusters_ids]

    # -- Compute clusters index once and for all:
    index_q = {}
    for n,q in enumerate(clusters_ids):
        index_q[q] = np.where(clustering == q)[0]

    vertex_cluster_distance = {}
    for i in xrange(N):
        tmp = np.zeros( shape=[nb_clusters], dtype=float )
        for n,q in enumerate(clusters_ids):
            tmp[n] = sum([distance_matrix[i,j] for j in index_q[q] if i!=j])/(nb_ids_by_clusters[n]-1)

        vertex_cluster_distance[i] = tmp

    return vertex_cluster_distance


def vertex_distance2cluster_center(distance_matrix, clustering):
    """
    Compute the distance between a vertex and the center of its group.

    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    N = len(clustering)

    D_iq = vertex2clusters_distance(distance_matrix, clustering)
    vertex_distance2center = {}
    for i in xrange(N):
        vertex_distance2center[i] = D_iq[i][clustering[i]]

    return vertex_distance2center


def plot_vertex_distance2cluster_center(distance_matrix, clustering):
    """
    Plot the distance between a vertex and the center of its group.

    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    vtx2center = vertex_distance2cluster_center(distance_matrix, clustering)

    N = len(clustering)
    clusters_ids = list(set(clustering))
    nb_ids_by_clusters = [len(np.where(clustering == q)[0]) for q in clusters_ids]

    # -- Compute clusters index once and for all:
    index_q = {}
    fig = plt.figure()
    for n,q in enumerate(clusters_ids):
        index_q[q] = np.where(clustering == q)[0]
        vector = [vtx2center[i] for i in index_q[q]]
        vector.sort()
        plt.plot( vector, label = "Cluster "+str(clusters_ids[n]) )
        plt.axis([0,max(nb_ids_by_clusters), min(vtx2center.values()), max(vtx2center.values())])

    plt.legend(ncol=3)


def time_point_by_clusters(graph, distance_matrix, clustering):
    """
    Return the representativity of each time points (index+1) by clusters.

    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    N = len(clustering)
    clusters_ids = list(set(clustering))
    ids_by_clusters = dict( (q, np.where(clustering == q)[0]) for q in clusters_ids )
    nb_ids_by_clusters = dict( (q, len(np.where(clustering == q)[0])) for q in clusters_ids )
    index_time_points = list(set(graph.vertex_property('index').values()))

    percent = {}
    for q in clusters_ids:
        for t in index_time_points:
            nb = len([k for k in ids_by_clusters[q] if k in graph.vertex_ids_at_time(t)])
            percent[(q,t)] = [nb/float(nb_ids_by_clusters[q]), nb]

    return dict( ((q,t+1), percent[(q,t)]) for q in clusters_ids for t in index_time_points if percent[(q,t)][0] != 0)


def forward_projection_match(graph, distance_matrix, clustering):
    """
    Compute the temporal evolution of ids of each clusters.
    """
    N = len(clustering)
    clusters_ids = list(set(clustering))
    ids_by_clusters = dict( (q, np.where(clustering == q)[0]) for q in clusters_ids )
    nb_ids_by_clusters = dict( (q, len(np.where(clustering == q)[0])) for q in clusters_ids )
    index_time_points = list(set(graph.vertex_property('index').values()))

    
    forward_projection_cluster = {}
    for t in index_time_points[:-1]:
        for vid in graph.vertex_ids_at_time(t):
            if graph.has_children(vid):
                children = graph.children(vid)
                for vid_children in children:
                    forward_projection_cluster[(t,clustering[vid])] = 1

    return None

def isolation(distance_matrix, clustering, index = None):
    """
    Isolation, a vertex is isolated if the max distance to a member of its group if superior to the min distance to a member of another group.

    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    N = len(clustering)
    clusters_ids = list(set(clustering))
    nb_clusters = len(clusters_ids)

    if index is None:
        index = xrange(N)

    index_q = {}
    for q in clusters_ids:
        index_q[q] = np.where(clustering == q)[0]
        index_not_q = np.where(clustering != q)[0]

    isolated = {}
    for i in index:
        if max([distance_matrix[i,j] for j in index_q]) >  min([distance_matrix[i,j] for j in index_not_q if distance_matrix[i,j]!=0]):
            isolated[i] = True
        else:
            isolated[i] = False

    if len(index) == 1:
        return isolated.values()
    else:
        return isolated



# An implementation of the gap statistic algorithm from Tibshirani, Walther, and Hastie's "Estimating the number of clusters in a data set via the gap statistic".
#~ library(plyr)
#~ library(ggplot2)

#~ # Given a matrix `data`, where rows are observations and columns are individual dimensions, compute and plot the gap statistic (according to a uniform reference distribution).
#~ def gap_statistic(data, min_num_clusters = 1, max_num_clusters = 10, num_reference_bootstraps = 10):
    #~ num_clusters = xrange(min_num_clusters,max_num_clusters)
    #~ actual_dispersions = maply(num_clusters, function(n) dispersion(data, n))
    #~ ref_dispersions = maply(num_clusters, function(n) reference_dispersion(data, n, num_reference_bootstraps))
    #~ mean_ref_dispersions = ref_dispersions[ , 1]
    #~ stddev_ref_dispersions = ref_dispersions[ , 2]
    #~ gaps = mean_ref_dispersions - actual_dispersions
#~ 
    #~ print(plot_gap_statistic(gaps, stddev_ref_dispersions, num_clusters))
#~ 
    #~ print("The estimated number of clusters is ", num_clusters[which.max(gaps)], ".", sep = ""))
#~ 
    #~ list(gaps = gaps, gap_stddevs = stddev_ref_dispersions)
#~ 
#~ 
#~ # Plot the gaps along with error bars.
#~ def plot_gap_statistic(gaps, stddevs, num_clusters) {
    #~ plot(num_clusters, gaps, xlab = "# clusters", ylab = "gap", geom = "line", main = "Estimating the number of clusters via the gap statistic") + geom_errorbar(aes(num_clusters, ymin = gaps - stddevs, ymax = gaps + stddevs), size = 0.3, width = 0.2, colour = "darkblue")
    #~ plt.plot(num_clusters, gaps, xlab = "# clusters", ylab = "gap")
    #~ geom_errorbar(aes(num_clusters, ymin = gaps - stddevs, ymax = gaps + stddevs), size = 0.3, width = 0.2, colour = "darkblue")
    #~ plt.title("Estimating the number of clusters via the gap statistic")
#~ 
#~ # Calculate log(sum_i(within-cluster_i sum of squares around cluster_i mean)).
#~ def dispersion(data, num_clusters):
    #~ # R's k-means algorithm doesn't work when there is only one cluster.
    #~ if (num_clusters == 1):
        #~ cluster_mean = np.mean(data, 0) #column wise mean
        #~ distances_from_mean = aaply((data - cluster_mean)^2, 1, sum)
        #~ log(sum(distances_from_mean))
    #~ else:
        #~ # Run the k-means algorithm `nstart` times. Each run uses at most `iter.max` iterations.
        #~ k = kmeans(data, centers = num_clusters, nstart = 10, iter.max = 50)
        #~ # Take the sum, over each cluster, of the within-cluster sum of squares around the cluster mean. Then take the log. This is `W_k` in TWH's notation.
        #~ log(sum(k$withinss))
#~ 
#~ 
#~ 
#~ # For an appropriate reference distribution (in this case, uniform points in the same range as `data`), simulate the mean and standard deviation of the dispersion.
#~ def reference_dispersion(data, num_clusters, num_reference_bootstraps):
    #~ dispersions = maply(1:num_reference_bootstraps, function(i) dispersion(generate_uniform_points(data), num_clusters))
    #~ mean_dispersion = mean(dispersions)
    #~ stddev_dispersion = sd(dispersions) / sqrt(1 + 1 / num_reference_bootstraps) # the extra factor accounts for simulation error
    #~ c(mean_dispersion, stddev_dispersion)
#~ 
#~ 
#~ # Generate uniform points within the range of `data`.
#~ def generate_uniform_points(data):
    #~ # Find the min/max values in each dimension, so that we can generate uniform numbers in these ranges.
    #~ mins = aaply(data, 2, min)
    #~ maxs = apply(data, 2, max)
#~ 
    #~ num_datapoints = nrow(data)
    #~ # For each dimension, generate `num_datapoints` points uniformly in the min/max range.
    #~ uniform_pts = maply(1:length(mins), function(dim) runif(num_datapoints, min = mins[dim], max = maxs[dim]))
    #~ uniform_pts = t(uniform_pts)



def MDS_graph(dataset, dimension=3, clustering=None):
    """
    Compute and display a Multi Dimensional Scaling of the `dataset`.

    :Parameters:
     - `dataset` (np.array) - set of similarities/dissimilarities to embedd.
     - `dimension` () - number of dimension in which to immerse the provided dataset.
     - `clustering` (list) - list of int giving a label to each point (should be ordered the same way than dataset)
     - `dissimilarity` (string) - Which dissimilarity measure to use. Supported are `euclidean` and `precomputed`.
    """
    if dimension != 2 and dimension!=3:
        print "There will be no representation of the embedded coordinates."

    if clustering is None:
        D = dataset # BE CAREFUL: D is not defined.
        clustering = xrange(len(D))

    from sklearn import manifold
    xy = manifold.MDS(n_components=dimension, metric=True, dissimilarity='precomputed').fit_transform(dataset)
    fig = plt.figure()
    if dimension == 2:
        plt.scatter(xy[:,0],xy[:,1], c=clustering, cmap=cm.jet)
    elif dimension == 3:
        ax = Axes3D(fig)
        ax.scatter(xy[:,0],xy[:,1],xy[:,2], c=clustering, marker='s', cmap=cm.jet)
        ax_min = np.min([xy[:,0],xy[:,1],xy[:,2]]); ax_max=np.max([xy[:,0],xy[:,1],xy[:,2]])
        ax.set_xlim3d(ax_min, ax_max)
        ax.set_ylim3d(ax_min, ax_max)
        ax.set_zlim3d(ax_min, ax_max)
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
    else:
        return xy

    plt.title("Multi Dimensional Scaling")
    plt.show()


def plot_cluster_distances(cluster_distances):
    """
    Display a heat-map of cluster distances with matplotlib.
    """
    fig = plt.figure()
    plt.imshow(cluster_distances, cmap=cm.jet, interpolation='nearest')

    numrows, numcols = cluster_distances.shape
    cd = cluster_distances
    def format_coord(x, y):
        col = int(x+0.5)
        row = int(y+0.5)
        if col>=0 and col<numcols and row>=0 and row<numrows:
            z = cd[row,col]
            return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, z)
        else:
            return 'x=%1.4f, y=%1.4f'%(x, y)

    plt.format_coord = format_coord
    plt.title("Cluster distances heat-map")
    plt.colorbar()
    plt.show()


def distance_matrix_histogram(distance_matrix, frequency = True, bins = None):
    """
    Plot the histogram of the upper part of the distance matrix.
    """
    N = distance_matrix.shape[0]
    if bins == None:
        bins = np.floor_divide(N,3)

    fig = plt.figure()
    plt.hist([distance_matrix[i,j] for i in xrange(N) for j in xrange(N) if j>i and not np.isnan(distance_matrix[i,j])], bins=bins, normed = False if frequency else True)
    plt.title("Distance matrix Histogram")
    plt.xlabel("Pairwise distance")
    plt.ylabel("Frequency" if frequency else "Relative Frequency")
    

def sorted_k_dist_graph(dist_matrix, k =None, plot=True, filled=False, rstride=20, cstride=20):
    """
    :Parameters:
     - `dist_matrix` (np.array) - pairwise distance matrix.
     - `k` (None|int|list/tuple) - k nearest neighbors to look at. On can search for a range of k.
     - `plot` (bool) - specify if we shoud plot the sorted k-dist graph.

    """
    distance_matrix = copy.copy(dist_matrix)
    if not isinstance(distance_matrix,ndarray):
        distance_matrix = np.array(distance_matrix)

    assert distance_matrix.shape[0] == distance_matrix.shape[1]

    distance_matrix.sort(axis=1)
    N = distance_matrix.shape[0]

    sorted_k_dist = {}

    if isinstance(k,int):
        k_dist = distance_matrix[:,k]
        sorted_k_dist[k] = sorted(k_dist)
        dimension = 2
    elif k is None:
        k_range = [1,N]
        dimension = 3
    elif len(k)==2:
        start,stop = min(k),max(k)
        assert start > 0
        assert stop <= N
        k_range = [start,stop]
        dimension = 3

    if dimension == 3:
        for k in xrange(k_range[0],k_range[1]):
            k_dist = distance_matrix[:,k]
            sorted_k_dist[k] = sorted(k_dist)

    if plot:
        fig = plt.figure()
        if dimension == 2:
            plt.scatter(x = xrange(N), y = sorted_k_dist[k])
            plt.title("Sorted k-dist graph, for k={0}".format(k))
            plt.xlabel('Distance sorted ranks')
            plt.ylabel('k-distance')
        if dimension == 3:
            dataset=np.array([sorted_k_dist[k] for k in sorted_k_dist])
            ax = fig.add_subplot(1,1,1, projection='3d')
            xx, yy = np.mgrid[:dataset.shape[0],:dataset.shape[1]]
            if filled:
                surf = ax.plot_surface(xx, yy, dataset, rstride = rstride, cstride = cstride, cmap=cm.jet, linewidth=0)
            else:
                surf = ax.plot_wireframe(xx, yy, dataset, rstride = rstride, cstride = cstride)
            Z_min = np.min(dataset)
            Z_max = np.max(dataset)
            ax.set_zlim3d(Z_min, Z_max)
            ax.set_xlabel('k-neighbor')
            ax.set_ylabel('Distance sorted ranks')
            ax.set_zlabel('k-distance')
            plt.title( "Sorted k-dist graph, for k=[{0},{1}]".format(k_range[0],k_range[1]) )
            if filled:
                fig.colorbar(surf)

        plt.show()

    return sorted_k_dist

