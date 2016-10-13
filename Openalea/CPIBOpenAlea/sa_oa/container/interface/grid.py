# -*- python -*-
# -*- coding: utf-8 -*-
#
#       Grid : grid package
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
This module provide a grid interface
"""

__license__= "Cecill-C"
__revision__=" $Id: grid.py 7865 2010-02-08 18:04:39Z cokelaer $ "

class IGrid (object) :
    """
    interface definition of simple N dimensional grids
    with finite number of case per dimension
    """

    def dim (self) :
        """
        dmension of the grid
        number of coordinates
        :rtype: int
        """
        raise NotImplementedError

    def shape (self) :
        """
        return the shape of the grid,
        number of cases per dimension

        :rtype: iter of int
        """
        raise NotImplementedError

class ICaseListGrid (object) :
    """
    grid is seen as a collection of cases
    regularly positionned in space
    """
    def __len__ (self) :
        """
        number of cases in the grid

        :rtype: int
        """
        raise NotImplementedError

    def __iter__ (self) :
        """
        iterator on case indexes

        :rtype: iter of int
        """
        raise NotImplementedError

    def index (self, coord ) :
        """
        compute the index of a case from his position
        inverse function of `coordinates`

        :param coord: position in each dimension
        :type coord: tuple of int
        :rtype: int
        """
        raise NotImplementedError

    def coordinates (self, ind) :
        """
        compute the position in each dimension from the index of the case
        inverse function of `index`

        :param ind: index of the case
        :type ind: int
        :rtype: tuple of int
        """
        raise NotImplementedError
