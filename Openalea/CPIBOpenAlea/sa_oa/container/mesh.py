# -*- python -*-
# -*- coding: utf-8 -*-
#
#       Topomesh : container package
#
#       Copyright or  or Copr. 2006 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <revesansparole@gmail.com>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       VPlants WebSite : https://gforge.inria.fr/projects/vplants/
#

__doc__="""
This module provide a simple pure python implementation
for a mesh, i.e. a topomesh associated with point positions
"""

__license__= "Cecill-C"
__revision__=" $Id: mesh.py 14917 2013-09-27 12:28:55Z pradal $ "

__all__ = ["Mesh"]


from sa_oa.container.topomesh import Topomesh, InvalidDart, TopomeshError

class Mesh (Topomesh) :
    """Geometrical mesh

    Associate position to elements of a topomesh
    """

    def __init__ (self, id_gen = 'max') :
        """Constructor of an empty mesh
        
        :Parameters:
         - `idgenerator` (str) - type of id gen used in the
                                 mesh
        """
        Topomesh.__init__(self, id_gen)
        
        self._pos = {}
    
    def clear (self) :
        Topomesh.clear(self)
        self._pos.clear()
    
    def position (self, did) :
        """Return position in space of a given dart
        
        if dart degree is greater than 0, a centroid
        position is returned
        
        :Parameters:
         - `did` (did) - id of dart
        
        :Returns: Vec or None if position is undefined
        """
        if self.degree(did) == 0 :
            return self._pos[did]
        else : #compute centroid
            pos = tuple(self.position(bid) \
                        for bid in self.borders(did, self.degree(did) ) )
            if len(pos) == 0 or None in pos :
                return None
            else :
                return reduce(lambda x, y: x + y, pos) / len(pos)
    
    def set_position (self, did, pos) :
        if self.degree(did) > 0 :
            raise TopomeshError("degree of dart must be 0")
        
        self._pos[did] = pos
    
    def apply_geom_transfo (self, transfo) :
        """Apply a transfo on all points positions
        
        :Parameters:
         - `transfo` (func) - function that take a pos as
                              argument and return new pos
        """
        for pid, vec in self._pos.items() :
            self._pos[pid] = transfo(vec)
    
    def add_dart (self, degree, did = None) :
        did = Topomesh.add_dart(self, degree, did)
        if degree == 0 :
            self._pos[did] = None
        
        return did
    
    def remove_dart (self, did) :
        Topomesh.remove_dart(self, did)
        self._pos.pop(did, None)
    
    def local_view (self, did) :
        """Return a local view on a part of a mesh
        
        This local view contains the given dart and all its
        borders until borders of degree 0.
        
        :Returns: (MeshView)
        """
        return MeshView(self, did)


class MeshView (object) :
    """A local view on a mesh
    """
    
    def __init__ (self, mesh, ref_did) :
        self._mesh = mesh
        self._ref_did = ref_did

        self._local_dids = set([ref_did])
        front = [ref_did]
        for i in range(mesh.degree(ref_did) ) :
            bids = set()
            for did in front :
                bids |= set(mesh.borders(did) )
            
            front = bids
            self._local_dids |= front
    
    #copy mesh
    def instance (self) :
        """Return a Mesh object, copy of the information
        this view is referencing
        """
        m = Mesh()
        
        #add darts
        for did in self._local_dids :
            deg = self._mesh.degree(did)
            m.add_dart(deg, did)
            if deg == 0 :
                m.set_position(did, self._mesh.position(did) )#TODO copy?
        
        #add links
        for did in self._local_dids :
            for bid in self.borders(did) :
                m.link(did, bid)
        
        #return
        return m
    
    #view function
    def degree (self, did) :
        return self._mesh.degree(did)
    
    def has_dart (self, did) :
        return did in self._local_dids
    
    def borders (self, did, offset = 1) :
        if did not in self._local_dids :
            raise InvalidDart(did)
        
        #no need for fancy testing since all borders
        #are in the view by definition
        #contrary to :func:regions below
        return self._mesh.borders(did, offset)
    
    def nb_borders (self, did) :
        if did not in self._local_dids :
            raise InvalidDart(did)
        
        return self._mesh.nb_borders(did)
    
    def regions (self, did, offset = 1) :
        if did not in self._local_dids :
            raise InvalidDart(did)
        
        for rid in self._mesh.regions(did, offset) :
            if rid in self._local_dids :
                yield rid
    
    def nb_regions (self, did) :
        return len(tuple(self.regions(did) ) )
    
    def _darts (self, degree) :
        for did in self._local_dids :
            if self._mesh.degree(did) == degree :
                yield did
    
    def darts (self, degree = None) :
        if degree is None :
            return iter(self._local_dids)
        else :
            return self._darts(degree)
    
    def nb_darts (self, degree = None) :
        if degree is None :
            return len(self._local_dids)
        else :
            nb = 0
            for did in self._local_dids :
                if self._mesh.degree(did) == degree :
                    nb += 1
            
            return nb
    
    def border_neighbors (self, did) :
        for bid in self._mesh.border_neighbors(did) :
            if bid in self._local_dids :
                yield bid
    
    def nb_border_neighbors (self, did) :
        return len(tuple(self.border_neighbors(did) ) )
    
    def region_neighbors (self, did) :
        for rid in self._mesh.region_neighbors(did) :
            if rid in self._local_dids :
                yield rid
    
    def nb_region_neighbors (self, did) :
        return len(tuple(self.region_neighbors(did) ) )
    
    def position (self, did) :
        if did not in self._local_dids :
            raise InvalidDart(did)
        
        return self._mesh.position(did)








