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

'''Display a graph using graphviz.

Display a graph with the graphviz package by using the pydot module.
'''

__docformat__ = "restructuredtext"
__license__ = "Cecill-C"
__revision__ = " $Id: graphviz.py 7865 2010-02-08 18:04:39Z cokelaer $ "

import pydot
from sa_oa.container.traversal.tree import pre_order

def graph_to_pydot(g, **kwds):
    """ Return a pydot Graph from a graph.

    :Parameters: 

        - `g`: a graph implementing :class:`sa_oa.container.interface.graph.IEdgeListGraph` interface.
        -`prog`: default is dot. value must one of :

            * dot for hierarchical layout of directed graph
            * neato and fdp for spring model layouts
            * twopi for radial layout
            * circo for circular layout

    :Returns: 

        - A pydot graph.

    :Examples:

        >>> dotg = graph_to_pydot(g)
        >>> dotg.write_svg('toto.svg', prog='circo')

    """
    edges = set((str(g.source(eid)), str(g.target(eid))) for eid in g.edges())
    graph = pydot.graph_from_edges(edges, directed=True)

    return graph

def tree_to_pydot(t, root=None, **kwds):
    """ Return a pydot directed Graph from a tree.

    :Parameters: 
        - `t`: a tree implementing :class:`sa_oa.container.interface.tree.ITree` interface.
    :Returns: 
        - A pydot graph.

    :Examples:

        >>> dotg = tree_to_pydot(t)
        >>> dotg.write_svg('toto.svg', prog='dot')

    """
    if root is None:
        root = t.root
    edges = set((str(t.parent(vid)), str(vid)) for vid in pre_order(t, root) if vid != root)
    graph = pydot.graph_from_edges(edges, directed=True)

    return graph


