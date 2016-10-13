# -*- python -*-
# -*- coding: utf-8 -*-
#
#       Topomesh : container package
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

"""
This module provide a python definition of a property
"""

__license__= "Cecill-C"
__revision__=" $Id: data_prop.py 14917 2013-09-27 12:28:55Z pradal $ "

class Quantity (object) :
    """Associate some semantic and unit to a value
    """
    
    def __init__ (self, value, unit = "", type = "", description = "") :
        """Constructor
        
        value: the actual value of the quantity
               may be a single value or a map of values
        unit: string description of the physical unit
        type: string description of the type of object
              (e.g.: 'float','int','Vector3')
        description: a string describing the quantity
        """
        self._value = value
        self._unit = unit
        self._type = type
        self._description = description
    
    ##########################################################
    #
    #        accessors
    #
    ##########################################################
    def value (self) :
        """Retrieve the value of this quantity.
        """
        return self._value
    
    def set_value (self, value) :
        """Set the value of this quantity.
        """
        self._value = value
    
    def unit (self) :
        """Retrieve unit of this quantity.
        """
        return self._unit
    
    def set_unit (self, unit) :
        """Set unit of this quantity.
        """
        self._unit = unit
    
    def type (self) :
        """Retrieve type of data in this quantity.
        """
        return self._type
    
    def set_type (self, type) :
        """Set type of data in this quantity.
        """
        self._type = type
    
    def description (self) :
        """Return the description associated with this quantity.
        """
        return self._description
    
    def set_description (self, description) :
        """Set the description associated with this quantity.
        """
        self._description = description
    
    ##########################################################
    #
    #        number interface
    #
    ##########################################################
    def __iadd__ (self, value) :
        """Add a value to this quantity.
        """
        if isinstance(self._value,dict) :
            for k,v in self._value.iteritems() :
                self._value[k] += value
        else :
            self._value += value
        return self
    
    def __isub__ (self, value) :
        """Substract a value from this quantity.
        """
        if isinstance(self._value,dict) :
            for k,v in self._value.iteritems() :
                self._value[k] -= value
        else :
            self._value -= value
        return self
    
    def __imul__ (self, value) :
        """Multiply this quantity by a value.
        """
        if isinstance(self._value,dict) :
            for k,v in self._value.iteritems() :
                self._value[k] *= value
        else :
            self._value *= value
        return self
    
    def __idiv__ (self, value) :
        """Divide this quantity by a value.
        """
        if isinstance(self._value,dict) :
            for k,v in self._value.iteritems() :
                self._value[k] /= value
        else :
            self._value /= value
        return self
    ##########################################################
    #
    #        dict interface for dict quantities
    #
    ##########################################################
    def clear (self) :
        self._value.clear()
    
    def __len__ (self) :
        return len(self._value)
    
    def __iter__ (self) :
        return iter(self._value)
    
    def itervalues (self) :
        return self._value.itervalues()
    
    def iteritems (self) :
        return self._value.iteritems()
    
    def __getitem__ (self, key) :
        return self._value[key]
    
    def __setitem__ (self, key, val) :
        self._value[key] = val
    
    def __delitem__ (self, key) :
        del self._value[key]
    
    def get (self, *args, **kwds) :
        return self._value.get(*args,**kwds)
    
    def pop (self, *args, **kwds) :
        return self._value.pop(*args,**kwds)
    
    def update (self, *args, **kwds) :
        self._value.update(*args,**kwds)

class DataProp (dict) :
    """
    Implementation of property with unit, type and description
    """
    
    def __init__ (self, *args, **kwds) :
        self._name = kwds.pop("name","")
        self._unit = kwds.pop("unit","")
        self._type = kwds.pop("type","")
        self._description = kwds.pop("description","")
        dict.__init__(self,*args,**kwds)
    
    ##########################################################
    #
    #        accessors
    #
    ##########################################################
    def name (self) :
        """Retrieve the name of this property.
        """
        return self._name
    
    def unit (self) :
        """Retrieve unit of this property.
        """
        return self._unit
    
    def type (self) :
        """Retrieve type of data in this property.
        """
        return self._type
    
    def description (self) :
        """Return the description associated with this property.
        """
        return self._description
