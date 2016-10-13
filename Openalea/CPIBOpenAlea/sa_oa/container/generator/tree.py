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
This module provides different traversal for tree data structure
implemented :class:`sa_oa.container.interface.tree` interface.
'''

__docformat__ = "restructuredtext"
__license__ = "Cecill-C"
__revision__ = " $Id: tree.py 7865 2010-02-08 18:04:39Z cokelaer $ "

from random import randint
from collections import deque

def regular_tree(tree, vtx_id, nb_children=3, nb_vertices=20):
    """ Build a regular tree where each vertex has the same number of children.

    :Parameters:
        - `tree` : an existing tree implementing the :class:sa_oa.container.interface.tree.ITree interface
        - `vtx_id` : a vertex id from which the generated tree will be added.
        - `nb_children` : the number of children added at each node.
        - `nb_vertices` : the number of vertices that will be added to the initial tree

    :Returns:
        - `tree`: the reference to the inital tree augmented by a sub_tree containing `nb_vertices`.

    :Example:

    >>> from sa_oa.container.tree import PropertyTree
    >>> from sa_oa.container.generator.tree import regular
    >>> tree = PropertyTree()
    >>> tree = regular(tree, tree.root, nb_children=3, nb_vertices=20)
    >>> print len(tree)
    """

    vid = vtx_id
    l=[vid]
    while nb_vertices > 0:
        n = min(nb_children, nb_vertices)
        vid = l.pop(0)
        for i in range(n):
            v = tree.add_child(vid)
            nb_vertices -= 1
            l.append(v)
    return tree

def random_tree(tree, root, nb_children=3, nb_vertices=20):
    """ Build a random tree where each vertex has a random number of children between 1 and `nb_children`.

    :Parameters:
        - `tree` : an existing tree implementing the :class:sa_oa.container.interface.tree.ITree interface
        - `root` : a vertex id from which the generated tree will be added.
        - `nb_children` : the number of children added at each node.
        - `nb_vertices` : the number of vertices that will be added to the initial tree

    :Returns:
        - `vid`: the last vertex id added to the inital tree.

    :Example:

    >>> from sa_oa.container.tree import PropertyTree
    >>> from sa_oa.container.generator.tree import random_tree
    >>> tree = Tree()
    >>> vid = random_tree(tree, tree.root, nb_children=3, nb_vertices=20)
    >>> print len(tree)
    """

    vid = root
    l = deque([vid])
    while nb_vertices > 0:
        n = min(randint(1,nb_children), nb_vertices)
        vid = l.popleft()
        for i in range(n):
            edge_type = '+'
            if i == n/2:
                edge_type='<'

            v=tree.add_child(vid, edge_type=edge_type)
            nb_vertices -= 1
            l.append(v)
    return l[-1]
