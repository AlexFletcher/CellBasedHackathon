# -*- python -*-
# -*- coding: utf-8 -*-
#
#       ITopoMesh : container package
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
This module provide a topological mesh interface
"""

__license__= "Cecill-C"
__revision__=" $Id: topomesh.py 13343 2012-12-07 04:08:32Z revesansparole $ "

class TopomeshError (Exception) :
    """
    base class for all exception in a topomesh
    """

class InvalidWisp (TopomeshError,KeyError) :
    """
    exception raised when a wrong wisp id is provided
    """

class InvalidDegree (TopomeshError,ValueError) :
    """
    exception raised when a the degree is outside of bounds
    """

class ITopomesh (object) :
    """
    interface definition of a topological mesh
    a mesh is formed of elements called wisps
    separated by elements of degree-1
    """
    def degree (self) :
        """
        maximum degree (or scale) of elements of the mesh
        """
        raise NotImplementedError

    def is_valid (self) :
        """
        test wether the mesh fulfill all mesh properties
        """
        raise NotImplementedError

    def has_wisp (self, degree, wid) :
        """
        return true if the wisp
        specified by its id is in mesh
        """
        raise NotImplementedError

    def borders (self, degree, wid, degree_offset=1) :
        """
        iterator on all border of this wisp
        """
        raise NotImplementedError

    def nb_borders (self, degree, wid, degree_offset=1) :
        """
        number of border of this wisp
        """
        raise NotImplementedError

    def regions (self, degree, wid) :
        """
        iterator on all regions this wisp separate
        """
        raise NotImplementedError

    def nb_regions (self, degree, wid) :
        """
        number of regions this wisp separate
        """
        raise NotImplementedError

class IWispListMesh (object) :
    """
    mesh view as a collection of wisps
    """
    def wisps (self, degree) :
        """
        iterator on all wisps of a given degree
        """
        raise NotImplementedError

    def nb_wisps (self, degree) :
        """
        number of wisps of the given degree
        """
        raise NotImplementedError

class INeighborhoodMesh (object) :
    """
    implicit neighborhood between wisps at the same scale
    """
    def border_neighbors (self, degree, wid) :
        """
        iterator on all wisps at the same degree
        that share a border with this wisp
        """
        raise NotImplementedError

    def nb_border_neighbors (self, degree, wid) :
        """
        number of border_neighbors of this wisp
        """
        raise NotImplementedError

    def region_neighbors (self, degree, wid) :
        """
        iterator on all wisps at the same degree
        that share a region with this wisp
        """
        raise NotImplementedError

    def nb_region_neighbors (self, degree, wid) :
        """
        number of region_neighbors of this wisp
        """
        raise NotImplementedError

class IMutableMesh (object) :
    """
    interface for editing methods on mesh
    """
    def add_wisp (self, degree, wid=None) :
        """
        add a new wisp connected to nothing
        if wid is None, create a free id
        return used wid
        """
        raise NotImplementedError

    def remove_wisp (self, degree, wid) :
        """
        remove wisp from the mesh
        remove all attached links
        """
        raise NotImplementedError

    def link (self, degree, wid, border_id) :
        """
        link a wisp of degree degree with id wid
        with another wisp of degree degree-1 with id border_id
        """
        raise NotImplementedError

    def unlink (self, degree, wid, border_id) :
        """
        remove links between a region and its border
        """
        raise NotImplementedError
