# -*- python -*-
# -*- coding: utf-8 -*-
#
#       Relation : container package
#
#       Copyright or  or Copr. 2006 INRIA - CIRAD - INRA
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
This module provide a simple pure python implementation
for a IRelation interface
"""

__license__= "Cecill-C"
__revision__=" $Id: relation.py 10381 2011-04-04 11:52:57Z pradal $ "

from interface.relation import (InvalidLeft,InvalidRight,InvalidLink,
                                IRelation,ILeftListRelation,IRightListRelation,
                                ILinkRelation,IMutableRelation)
from utils import IdDict

class StrInvalidLeft (InvalidLeft) :
    """
    exception raised when a wrong left id is provided
    """
    def __init__ (self, elmid) :
        InvalidLeft.__init__(self,"left element %d does not exist" % elmid)

class StrInvalidRight (InvalidRight) :
    """
    exception raised when a wrong right id is provided
    """
    def __init__ (self, elmid) :
        InvalidRight.__init__(self,"right element %d does not exist" % elmid)

class StrInvalidLink (InvalidLink) :
    """
    exception raised when a link between two elements does not exist
    """
    def __init__ (self, lid) :
        InvalidLink.__init__(self,"link %d does not exist" % lid)

class Relation (IRelation,ILeftListRelation,IRightListRelation,ILinkRelation,IMutableRelation) :
    """
    implementation of a relation
    """
    def __init__ (self, idgenerator = "set") :
        """
        constructor of an empty relation
        """
        self._left_links=IdDict(idgenerator = idgenerator)
        self._right_links=IdDict(idgenerator = idgenerator)
        self._link_extremities=IdDict(idgenerator = idgenerator)

    ########################################################################
    #
    #               Relation concept
    #
    ########################################################################
    def is_valid (self) :
        return True
    is_valid.__doc__=IRelation.is_valid.__doc__

    def has_left (self, elmid) :
        return elmid in self._left_links
    has_left.__doc__=IRelation.has_left.__doc__

    def has_right (self, elmid) :
        return elmid in self._right_links
    has_right.__doc__=IRelation.has_right.__doc__

    def has_link (self, lid) :
        return lid in self._link_extremities
    has_link.__doc__=IRelation.has_link.__doc__
    ########################################################################
    #
    #               Left list concept
    #
    ########################################################################
    def left_elements (self) :
        return self._left_links.iterkeys()
    left_elements.__doc__=ILeftListRelation.left_elements.__doc__

    def nb_left_elements (self) :
        return len(self._left_links)
    nb_left_elements.__doc__=ILeftListRelation.nb_left_elements.__doc__
    ########################################################################
    #
    #               Right list concept
    #
    ########################################################################
    def right_elements (self) :
        return self._right_links.iterkeys()
    right_elements.__doc__=IRightListRelation.right_elements.__doc__

    def nb_right_elements (self) :
        return len(self._right_links)
    nb_right_elements.__doc__=IRightListRelation.nb_right_elements.__doc__
    ########################################################################
    #
    #               Link Relation concept
    #
    ########################################################################
    def links (self) :
        return self._link_extremities.iterkeys()
    links.__doc__=ILinkRelation.links.__doc__

    def from_left (self, elmid) :
        try :
            return iter(self._left_links[elmid])
        except KeyError :
            raise StrInvalidLeft(elmid)
    from_left.__doc__=ILinkRelation.from_left.__doc__

    def nb_links_from_left (self, elmid) :
        try :
            return len(self._left_links[elmid])
        except KeyError :
            raise StrInvalidLeft(elmid)
    nb_links_from_left.__doc__=ILinkRelation.nb_links_from_left.__doc__

    def from_right (self, elmid) :
        try :
            return iter(self._right_links[elmid])
        except KeyError :
            raise StrInvalidRight(elmid)
    from_right.__doc__=ILinkRelation.from_right.__doc__

    def nb_links_from_right (self, elmid) :
        try :
            return len(self._right_links[elmid])
        except KeyError :
            raise StrInvalidRight(elmid)
    nb_links_from_right.__doc__=ILinkRelation.nb_links_from_right.__doc__

    def left (self, lid) :
        try :
            return self._link_extremities[lid][0]
        except KeyError :
            raise StrInvalidLink(lid)
    left.__doc__=ILinkRelation.left.__doc__

    def right (self, lid) :
        try :
            return self._link_extremities[lid][1]
        except KeyError :
            raise StrInvalidLink(lid)
    right.__doc__=ILinkRelation.right.__doc__
    ########################################################################
    #
    #               Mutable relation concept
    #
    ########################################################################
    def add_left_element (self, elmid=None) :
        return self._left_links.add(set(),elmid)
    add_left_element.__doc__=IMutableRelation.add_left_element.__doc__

    def remove_left_element (self, elmid) :
        for lid in list(self.from_left(elmid)) :
            self.remove_link(lid)
        del self._left_links[elmid]
    remove_left_element.__doc__=IMutableRelation.remove_left_element.__doc__

    def add_right_element (self, elmid=None) :
        return self._right_links.add(set(),elmid)
    add_right_element.__doc__=IMutableRelation.add_right_element.__doc__

    def remove_right_element (self, elmid) :
        for lid in list(self.from_right(elmid)) :
            self.remove_link(lid)
        del self._right_links[elmid]
    remove_right_element.__doc__=IMutableRelation.remove_right_element.__doc__

    def add_link (self, left_elmid, right_elmid, lid=None) :
        if not self.has_left(left_elmid) :
            raise StrInvalidLeft(left_elmid)
        if not self.has_right(right_elmid) :
            raise StrInvalidRight(right_elmid)
        lid=self._link_extremities.add( (left_elmid,right_elmid),lid )
        self._left_links[left_elmid].add(lid)
        self._right_links[right_elmid].add(lid)
        return lid
    add_link.__doc__=IMutableRelation.add_link.__doc__

    def remove_link (self, lid) :
        self._left_links[self.left(lid)].remove(lid)
        self._right_links[self.right(lid)].remove(lid)
        del self._link_extremities[lid]
    remove_link.__doc__=IMutableRelation.remove_link.__doc__


