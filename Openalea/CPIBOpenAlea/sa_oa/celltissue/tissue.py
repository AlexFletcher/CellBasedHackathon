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
This module defines the main Tissue object
"""

__license__= "Cecill-C"
__revision__=" $Id: $ "

from sa_oa.container import IdDict
from relation import GraphRelation,RelationRelation,MeshRelation

relation_type_def={GraphRelation.name:GraphRelation,
				RelationRelation.name:RelationRelation,
				MeshRelation.name:MeshRelation}

class TissueError (Exception) :
	"""
	base class of all tissue exceptions
	"""

class InvalidElement (TissueError,KeyError) :
	"""
	raised when a wrong element id has been provided
	"""

class InvalidElementType (TissueError,KeyError) :
	"""
	raised when a wrong element type id has been provided
	"""

class InvalidRelation (TissueError,KeyError) :
	"""
	raised when a wrong relation id has been provided
	"""

class Tissue (object) :
	"""
	base class for all tissues
	a tissue is a collection of elements
	and a set of links between these elements
	"""
	def __init__ (self) :
		self._type_name = IdDict()
		self._element_type = IdDict()
		self._relation_object = IdDict()
		self._type_relation = {}
	#################################################
	#
	#		ITissue
	#
	#################################################
	def _elements (self, type_id) :
		for elm_id,elm_type in self._element_type.iteritems() :
			if elm_type == type_id :
				yield elm_id
	
	def elements (self, type_id=None) :
		"""
		iterator on all elements of tissue of type elm_type_id
		if elm_type_id is None iterator on all elements
		return : iter of elm_id
		"""
		if type_id is None :
			return iter(self._element_type)
		else :
			return self._elements(type_id)
	
	def nb_elements (self, type_id=None) :
		"""
		number of elements of the given type
		total number of elements in the tissue if None
		return : int
		"""
		if type_id is None :
			return len(self._element_type)
		else :
			nb = 0
			for elm_id,elm_type in self._element_type.iteritems() :
				if elm_type == type_id :
					nb += 1
			return nb
	
	def type (self, elm_id) :
		"""
		return the id of the provided element
		return : elm_type_id
		"""
		try :
			return self._element_type[elm_id]
		except KeyError :
			raise InvalidElement(elm_id)
	
	def type_name (self, type_id) :
		"""
		return the provided name for this type
		return : str
		"""
		try :
			return self._type_name[type_id]
		except KeyError :
			raise InvalidElementType(type_id)
	
	def set_type_name (self, type_id, name) :
		"""
		modify the name of a given type
		"""
		try :
			self._type_name[type_id] = name
		except KeyError :
			raise InvalidElementType(type_id)
	
	def types (self) :
		"""
		iterator on all element types in this tissue
		return : iter of elm_type_id
		"""
		return iter(self._type_name)
	
	def relation (self, relation_id) :
		"""
		return a reference on the relation object
		specified by its id
		"""
		try :
			return self._relation_object[relation_id]
		except KeyError :
			raise InvalidRelation(relation_id)
	
	def relations (self, type_id=None) :
		"""
		iterator on all relations in which the type
		is involved
		if type_id is None, iterate on all relations in the tissue
		return : iter of relation_id
		"""
		if type_id is None :
			return iter(self._relation_object)
		else :
			try :
				return iter(self._type_relation[type_id])
			except KeyError :
				raise InvalidElementType(type_id)

	def nb_relations (self, type_id=None) :
		"""
		number of relation the given type is used for
		all relations maintained by the tissue if None
		return : int
		"""
		if type_id is None :
			return len(self._relation_object)
		else  :
			try :
				return len(self._type_relation[type_id])
			except KeyError :
				raise InvalidElementType(type_id)
	###################################################
	#
	#		IMutableTissue
	#
	###################################################
	def add_type (self, type_name="elm", type_id=None) :
		"""
		add a new type of elements in the tissue
		using the provided id if not None
		return the id used for this type
		return : type_id
		"""
		type_id = self._type_name.add(type_name,type_id)
		self._type_relation[type_id] = set()
		return type_id
	
	def add_relation (self, relation_type="graph", type_used=[], relation_id=None, id_gen = "set") :
		"""
		add a relation between different type of elements,
		if the provided type is None, create a dummy type instead
		return : id used for this relation
		"""
		type_used = list(type_used)
		for ind,type_id in enumerate(type_used) :
			if type_id is None :
				type_used[ind] = self.add_type("unknown")
		try :
			Relation = relation_type_def[relation_type]
		except KeyError :
			raise InvalidRelation("%s unkown, available relation are : %s" % (str(relation_type),str(relation_type_def.keys())))
		relation = Relation(self,tuple(type_used),id_gen)
		relation_id = self._relation_object.add(relation,relation_id)
		for type_id in type_used :
			self._type_relation[type_id].add(relation_id)
			#update of the relation with elements
			#already in the tissue
			for elm_id in self.elements(type_id) :
				relation.add_element(type_id,elm_id)
		return relation_id
	
	def add_element (self, type_id, elm_id=None, safe=True) :
		"""
		add a new element of the given type in the tissue
		return id used for this element
		return : elm_id
		"""
		try :
			elm_id = self._element_type.add(type_id,elm_id)
		except KeyError :
			if safe :
				raise InvalidElement(elm_id)
			else :
				return elm_id
		for relation_id in self.relations(type_id) :
			self.relation(relation_id).add_element(type_id,elm_id)
		return elm_id
	
	def remove_element (self, elm_id) :
		"""
		remove the given element and all links
		in all relations that points to this element
		"""
		type_id = self.type(elm_id)
		#remove this element from all
		#relation it is involved in
		for relation_id in self.relations(type_id) :
			self.relation(relation_id).remove_element(type_id,elm_id)
		#remove the element
		del self._element_type[elm_id]
	
	def clear (self) :
		"""
		clear all elements in the tissue
		"""
		self._type_name.clear()
		self._element_type.clear()
		self._relation_object.clear()
		self._type_relation.clear()

