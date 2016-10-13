from interface import IPropertyMap,INumericProperty

class DensityMap (IPropertyMap) :
	"""
	classical map plus some methods that depends on elm volume
	"""
	def __init__ (self, property_map) :
		self._prop=property_map
	#########################################################################
	#
	#		IPropertyMap
	#
	#########################################################################
	def associated_map (self) :
		return self._prop
	
	def __getitem__ (self, elm_id) :
		return self._prop[elm_id]
	
	def __setitem__ (self, elm_id, val) :
		self._prop[elm_id]=val
	
	def __iter__ (self) :
		return iter(self._prop)
	
	def iterkeys (self) :
		return self._prop.iterkeys()
	
	def itervalues (self) :
		return self._prop.itervalues()
	
	def iteritems (self) :
		return self._prop.iteritems()
	
	#########################################################################
	#
	#		DensityMap
	#
	#########################################################################
	def volume (self, elm_id) :
		"""
		return the volume of element elm_id
		"""
		return 1.
	
	def as_quantity (self) :
		return QuantityAdaptor(self)
	########################################################################
	#
	#		usefull functions to access quantities
	#		instead of densities
	#
	########################################################################
	def quantity (self, elm_id) :
		"""
		return the amount of substance in the element
		:rtype: float
		"""
		return self[elm_id]*self.volume(elm_id)
	
	def set_quantity (self, elm_id, amount) :
		"""
		store amount of substance in the element
		"""
		self[elm_id]=amount/self.volume(elm_id)
	
	def add_quantity (self, elm_id, amount) :
		"""
		add an amount of substance in the element
		"""
		self[elm_id]+=amount/self.volume(elm_id)
	
	def quantities (self) :
		"""
		iterator on all stored quantities
		:rtype: iter of (elm_id,float)
		"""
		for elm_id in self :
			yield elm_id,self.quantity(elm_id)
	
	def iterquantities (self) :
		"""
		iterator on all stored quantities
		:rtype: iter of float
		"""
		for elm_id in self :
			yield self.quantity(elm_id)

class QuantityAdaptor (IPropertyMap) :
	"""
	adaptor to look at a density property
	that return quantities by default
	"""
	def __init__ (self, density_property) :
		self._prop=density_property
	def __getitem__ (self, elm_id) :
		return self._prop.quantity(elm_id)
	def __setitem__ (self, elm_id, amount) :
		return self._prop.set_quantity(elm_id,amount)
	def __iter__ (self) :
		return iter(self._prop)
	def iterkeys (self) :
		return self._prop.iterkeys()
	def itervalues (self) :
		return self._prop.iterquantities()
	def iteritems (self) :
		return self._prop.quantities()

