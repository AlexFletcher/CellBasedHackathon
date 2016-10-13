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
This module provide a simple pure python implementation
for a grid interface
"""

__license__= "Cecill-C"
__revision__=" $Id: grid.py 7865 2010-02-08 18:04:39Z cokelaer $ "

from interface.grid import IGrid,ICaseListGrid

class Grid (IGrid,ICaseListGrid) :
    """
    interface definition of simple N dimensional grids
    with finite number of case per dimension
    """

    def __init__ (self, shape) :
        """
        constructor of a finite grid
        :param shape: number of case in each dimension
        :type shape: iter of int
        """
        self._shape=[int(s) for s in shape]
        offset=[1]
        for i,incr in enumerate(self._shape[:-1]) :
            offset.append(offset[i]*incr)
        self._offset=offset

    # ##########################################################
    #
    #               Grid concept
    #
    # ##########################################################
    def dim (self) :
        return len(self._shape)
    dim.__doc__=IGrid.dim.__doc__

    def shape (self) :
        return iter(self._shape)
    shape.__doc__=IGrid.shape.__doc__

    # ##########################################################
    #
    #               Case list concept
    #
    # ##########################################################
    def __len__ (self) :
        s=1
        for incr in self._shape : s*=incr
        return s
    __len__.__doc__=ICaseListGrid.__len__.__doc__

    def __iter__ (self) :
        return iter(xrange(self.__len__()))
    __iter__.__doc__=ICaseListGrid.__iter__.__doc__

    def index (self, coord ) :
        return sum([coord[i]*offset for i,offset in enumerate(self._offset)])
    index.__doc__=ICaseListGrid.index.__doc__

    def coordinates (self, ind) :
        if not (0<=ind<len(self)) :
            raise IndexError("index out of range index: %d max : %d" % (ind,len(self)))
        reste=ind
        coord=[]
        for i in xrange(self.dim()-1,-1,-1) :
            coord.append(reste/self._offset[i])
            reste=reste%self._offset[i]
        coord.reverse()
        return coord
    coordinates.__doc__=ICaseListGrid.coordinates.__doc__
