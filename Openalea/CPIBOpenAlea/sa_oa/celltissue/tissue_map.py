# -*- python -*-
#
#       celltissue: main tissue object and functions to use it
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

__doc__="""
This module defines map on a tissue that extends classical dictionaries
"""

__license__= "Cecill-C"
__revision__=" $Id: $ "

class TUniformMap (dict) :
	"""
	return a unique value for each elements of the tissue
	"""
	def __init__ (self, value, tissue, elm_type=None) :
		self._value=value
		self._tissue=tissue
		self._elm_type=elm_type
	#################################################################
	#
	#		numeric concept
	#
	#################################################################
	def __iadd__ (self, val) :
		self._value+=val
		return self
	
	def __isub__ (self, val) :
		self._value-=val
		return self
	
	def __imul__ (self, val) :
		self._value*=val
		return self
	
	def __idiv__ (self, val) :
		self._value/=val
		return self
	################################################################
	#
	#		dict concept
	#
	################################################################
	def __getitem__ (self, key) :
		return self._value
	
	def __setitem__ (self, key, val) :
		self._value=val
	
	def __iter__ (self) :
		return self._tissue.elements(self._elm_type)
	
	def iterkeys (self) :
		return iter(self)
	
	def itervalues (self) :
		val=self._value
		for wid in  self :
			yield val
	
	def iteritems (self) :
		val=self._value
		for wid in self :
			yield wid,val

