# -*- python -*-
#
#       tissueshape: function used to deal with tissue geometry
#
#       Copyright 2006 INRIA - CIRAD - INRA  
#
#       File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
# 
#       OpenAlea WebSite : http://sa_oa.gforge.inria.fr
#

__doc__ = """
This module defines functions to convert a topomesh into geometrical pgl objects
"""

__license__ = "Cecill-C"
__revision__ = " $Id: $ "

from numpy import array,add
from sa_oa.container import Quantity

def tovec (pos, unit = "m", centered = False) :
	"""Format a property into geometrical vectors
	
	inverse of :func:`totup`
	
	:Parameters:
	 - `pos` (dict of (pid|tuple of float) ) - geometrical
	         position of points in space
	 - `unit` (str) - spatial unit
	 - `centered` (bool) - if True, translate positions
	    so that the barycenter of pos is (0,0,0)
	
	:Returns Type: :class:sa_oa.container.Quantity
	"""
	#test for quantity, will be deprecated in the future I hope
	if isinstance(pos,Quantity) :
		prop = pos.value()
		unit = pos.unit()
		descr = pos.description()
	else :
		prop = pos
		descr = ""
	
	#change for array
	pos = Quantity(dict( (pid,array(tup) ) for pid,tup in prop.iteritems() ),
	               unit,
	               "array",
	               descr)
	
	if centered :
		bary = reduce(add,pos.itervalues() ) / len(pos)
		pos.set_value(dict( (pid,vec - bary) for pid,vec in pos.iteritems() ) )
	
	return pos

def totup (pos) :
	"""Format a geometrical vector quantity into tuples
	
	inverse of :func:`tovec`
	
	:Parameters:
	 - `pos` (Quantity of (pid|array) ) - geometrical
	         position of points in space
	
	:Returns Type: :class:sa_oa.container.Quantity
	"""
	#test for quantity, will be deprecated in the future I hope
	if isinstance(pos,Quantity) :
		prop = pos.value()
		unit = pos.unit()
		descr = pos.description()
	else :
		prop = pos
		unit = ""
		descr = ""
	
	#change for tuples
	vec = prop.itervalues().next()
	typ = vec.dtype.name
	
	try :
		dim = len(vec)
		if 'int' in typ :
			return Quantity(dict( (pid,tuple(int(v) for v in vec) ) \
			                       for pid,vec in prop.iteritems() ),
			                unit,
			                "tuple",
			                descr)
		elif 'float' in typ :
			return Quantity(dict( (pid,tuple(float(v) for v in vec) ) \
			                       for pid,vec in prop.iteritems() ),
			                unit,
			                "tuple",
			                descr)
		else :
			raise UserWarning("unhandled type of elements")
	except TypeError :
		if 'int' in typ :
			return Quantity(dict( (pid,int(val) ) \
			                       for pid,val in prop.iteritems() ),
			                unit,
			                "int",
			                descr)
		elif 'float' in typ :
			return Quantity(dict( (pid,float(val) ) \
			                       for pid,val in prop.iteritems() ),
			                unit,
			                "float",
			                descr)
		else :
			raise UserWarning("unhandled type of elements")

