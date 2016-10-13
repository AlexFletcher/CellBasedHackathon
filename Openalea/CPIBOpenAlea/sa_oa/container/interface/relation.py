# -*- python -*-
# -*- coding: utf-8 -*-
#
#       IRelation : container package
#
#       Copyright  or Copr. 2006 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       VPlants WebSite : https://gforge.inria.fr/projects/vplants/
#

__doc__="""
This module provide a topological relation interface
"""

__license__= "Cecill-C"
__revision__=" $Id: relation.py 7865 2010-02-08 18:04:39Z cokelaer $ "

class RelationError (Exception) :
    """
    base class for all exception in a relation
    """

class InvalidLeft (RelationError,KeyError) :
    """
    exception raised when a wrong id for left element is provided
    """

class InvalidRight (RelationError,KeyError) :
    """
    exception raised when a wrong id for right element is provided
    """

class InvalidLink (RelationError,KeyError) :
    """
    exception raised when a link between two elements does not exist
    """

class IRelation (object) :
    """
    interface definition of a topological relation
    a relation links elements of two separate sets E1 (left) and E2 (right)
    """
    def is_valid (self) :
        """
        test wether the relation fulfill all relations properties
        """
        raise NotImplementedError

    def has_left (self, elmid) :
        """
        return true if the element
        specified by its id is in the relation
        """
        raise NotImplementedError

    def has_right (self, elmid) :
        """
        return true if the element
        specified by its id is in the relation
        """
        raise NotImplementedError

    def has_link (self, lid) :
        """
        return True if the link
        specified by its id is in the relation
        """
        raise NotImplementedError

class ILeftListRelation (object) :
    """
    relation view as a collection of left elements
    """
    def left_elements (self) :
        """
        iterator on all left elements
        of this relation
        """
        raise NotImplementedError

    def nb_left_elements (self) :
        """
        number of elements in the left set
        of this relation
        """
        raise NotImplementedError

class IRightListRelation (object) :
    """
    relation view as a collection of right elements
    """
    def right_elements (self) :
        """
        iterator on all right elements
        of this relation
        """
        raise NotImplementedError

    def nb_right_elements (self) :
        """
        number of elements in the right set
        of this relation
        """
        raise NotImplementedError

class ILinkRelation (object) :
    """
    explicit links
    """
    def links (self) :
        """
        iterator on all links in the relation
        return : iter of lid
        """
        raise NotImplementedError

    def from_left (self, elmid) :
        """
        iterator on all links connected
        to an element of left set
        return : iter of lid
        """
        raise NotImplementedError

    def nb_links_from_left (self, elmid) :
        """
        number of links attached to this left elm
        return : int
        """
        raise NotImplementedError

    def from_right (self, elmid) :
        """
        iterator on all links connected
        to an element of right set
        return : iter of lid
        """
        raise NotImplementedError

    def nb_links_from_right (self, elmid) :
        """
        number of links attached to this right elm
        return : int
        """
        raise NotImplementedError

    def left (self, lid) :
        """
        left element corresponding to the source of the link
        return : cid
        """
        raise NotImplementedError

    def right (self, lid) :
        """
        right element corresponding to the target of the link
        return : pid
        """
        raise NotImplementedError

class IMutableRelation (object) :
    """
    interface for editing methods on relation
    """
    def add_left_element (self, elmid=None) :
        """
        add a new left element connected to nothing
        if elmid is None, create a free id
        return used elmid
        """
        raise NotImplementedError

    def remove_left_element (self, elmid) :
        """
        remove an element and all the references to attached
        right elements
        """
        raise NotImplementedError

    def add_right_element (self, elmid=None) :
        """
        add a new right element connected to nothing
        if elmid is None, create a free id
        return used elmid
        """
        raise NotImplementedError

    def remove_right_element (self, elmid) :
        """
        remove an element and all the references to attached
        left elements
        """
        raise NotImplementedError

    def add_link (self, left_elmid, right_elmid, lid=None) :
        """
        add a link between a left element
        and a right element that must already exist
        return the link id used
        """
        raise NotImplementedError

    def remove_link (self, lid) :
        """
        remove a link between two elements
        """
        raise NotImplementedError
