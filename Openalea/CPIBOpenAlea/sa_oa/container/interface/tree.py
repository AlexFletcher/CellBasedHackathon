# -*- python -*-
# -*- coding: utf-8 -*-
#
#       Graph : graph package
#
#       Copyright  or Copr. 2006 INRIA - CIRAD - INRA
#
#       File author(s): Christophe Pradal <christophe pradal at cirad fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://sa_oa.gforge.inria.fr/
#
################################################################################

"""
This module provide a set of tree concepts to form a tree interface.
"""

__license__ = "Cecill-C"
__revision__ = " $Id: tree.py 7865 2010-02-08 18:04:39Z cokelaer $ "

class ITree(object):
    """Rooted Tree interface.

    depth(vid), depth() and sub_tree(vid) can be extenal algorithms.
    """

    def parent(self, vtx_id):
        """
        Return the parent of `vtx_id`.

        :param vtx_id: The vertex identifier.
        :returns: vertex identifier
        """
        raise NotImplementedError

    def children(self, vtx_id):
        """
        Return a vertex iterator

        :param vtx_id: The vertex identifier.
        :returns: iter of vertex identifier
        """
        raise NotImplementedError


    def nb_children(self, vtx_id):
        """
        Return the number of children

        :param vtx_id: The vertex identifier.
        :rtype: int
        """
        raise NotImplementedError


    def siblings(self, vtx_id):
        """
        Return an iterator of vtx_id siblings.
        vtx_id is not include in siblings.

        :param vtx_id: The vertex identifier.
        :returns: iter of vertex identifier
        """
        raise NotImplementedError


    def nb_siblings(self, vtx_id):
        """
        Return the number of siblings

        :returns: int
        """
        raise NotImplementedError


    def is_leaf(self, vtx_id):
        """
        Test if `vtx_id` is a leaf.

        :returns: bool
        """
        raise NotImplementedError



class IOrderedTree(object):
    """
    An ordered tree is a rooted tree where an order relation is
    defined between chidren.
    """

    def first_child(self, vid):
        """
        Return the first child of vid

        :param vid: The vertex identifier.
        :returns: vid
        """
        raise NotImplementedError


    def last_child(self, vid):
        """
        Return the last child of vid

        :param vid: The vertex identifier.
        :returns: vid
        """
        raise NotImplementedError


    def previous_sibling(self, vid):
        """
        Return previous sibling

        :param vid: The vertex identifier.

        :returns: vid
        """
        raise NotImplementedError


    def next_sibling(self, vid):
        """
        Return next sibling

        :param vid: The vertex identifier.

        :returns: vid
        """
        raise NotImplementedError


class IMutableTree(object):
    """
    A mutable rooted tree is a tree with specific methods to add vertex.
    The substitute method is defined outside the interface.
    substitute(self,vid,tree)
    """

    def add_child(self, parent, child ):
        """
        Add a child at the end of children

        :param parent: The parent identifier.
        :param child: The child identifier.
        """
        raise NotImplementedError


    def insert_sibling(self, vid1, vid2):
        """
        Insert vid2 before vid1.

        :param vid1: a vertex identifier
        :param vid2: the vertex to insert
        """
        raise NotImplementedError

    def insert_parent(self, vtx_id, parent_id):
        """
        Insert parent_id between vtx_id and its actual parent.

        :param vtx_id: a vertex identifier
        :param parent_id: a vertex identifier
        """
        raise NotImplementedError


class IEditableTree(object):
    """
    An editable tree is a mutable tree where you can add
    or retrieve sub trees.
    """

    def sub_tree(self, vtx_id):
        """
        Return a reference of the tree rooted on `vtx_id` in O(1).

        :returns: Editable Tree
        """
        raise NotImplementedError

    def insert_sibling_tree(self, vid, tree ):
        """
        Insert a tree before the vid.
        vid and the root of the tree are siblings.
        Complexity have to be O(1) if tree comes from the actual tree
        ( tree= sef.sub_tree() )

        :param vid: vertex identifier
        :param tree: a rooted tree
        """
        raise NotImplementedError

    def add_child_tree(self, parent, tree):
        """
        Add a tree after the children of the parent vertex.
        Complexity have to be O(1) if tree == sub_tree()

        :param parent: vertex identifier
        :param tree: a rooted tree
        """
        raise NotImplementedError

    def remove_tree(self, vtx_id):
        """
        Remove the sub tree rooted on `vtx_id`.

        :returns: bool
        """
        raise NotImplementedError
