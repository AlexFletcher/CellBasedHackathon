# -*- coding: utf-8 -*-
# -*- python -*-
#
#       OpenAlea.Container
#
#       Copyright 2008-2009 INRIA - CIRAD - INRA
#
#       File author(s): Christophe Pradal <christophe.pradal.at.cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://sa_oa.gforge.inria.fr
#
###############################################################################

'''
This module provides an implementation of the graph interface using networkx backend.
Backend implementation provide a way to reuse existing algorithms implemented
in different graph libraries, and to compare algorithms and graph implementation.
'''

__docformat__ = "restructuredtext"
__license__ = "Cecill-C"
__revision__ = " $Id: networkx.py 11356 2011-11-16 17:44:53Z pradal $ "

import networkx as nx
from sa_oa.container import Graph

def to_networkx(g):
    """ Return a NetworkX Graph from a graph.

    :Parameters: 
        - `g`: a graph implementing :class:`sa_oa.container.interface.IEdgeListGraph` interface.

    :Returns: 
        - A NetworkX graph.

    """
    graph = nx.DiGraph()
    graph.add_node_from(g.vertices())
    graph.add_edge_from(( (g.source(eid), g.target(eid)) for eid in g.edges()))
    return graph

def from_networkx(graph, klass=Graph):
    """ Return a Graph from a NetworkX Directed graph.

    :Parameters: 
        - `graph` : A NetworkX graph.
        - `klass` : 

    :Returns: 
        - `g`: a :class:`~sa_oa.container.interface.Graph`.

    """
    if not graph.directed:
        graph = graph.to_directed()

    g = klass()
    for vid in graph.nodes_iter():
        g.add_vertex(vid)

    for source, target in graph.edges_iter():
        g.add_edge(source, target)

    return g
