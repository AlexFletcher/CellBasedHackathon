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
This module defines main types of relations
"""

__license__= "Cecill-C"
__revision__=" $Id: $ "

from sa_oa.container import IdDict,Graph,Relation,Topomesh

class TissueRelation (object) :
	"""
	abstract object to express relations between two elements
	"""
	
	name = "trel" #shortcut used to add this relation in a tissue
	
	def __init__ (self, tissue, type_used) :
		self._tissue = tissue
	
	def tissue (self) :
		"""Tissue in which this relation exists.
		"""
		return self._tissue
	
	def involved_elements (self) :
		"""Iterator on elements types
		involved int his relation.
		
		Order has no meaning.
		"""
		return iter([])
	
	#######################################################
	#
	#		Mutable
	#
	#######################################################
	def add_element (self, type_id, elm_id) :
		"""
		add a new element
		using the provided id
		"""
		pass
	
	def remove_element (self, type_id, elm_id) :
		"""
		remove an element and all links connected to it
		"""
		pass
	#######################################################
	#
	#		Dict translation
	#
	#######################################################
	def description (self) :
		"""
		return a tuple of arguments used as a description
		to reconstruct this relation
		"""
		pass
	
	def reduce (self) :
		"""
		return a representation of internal links as basic python elements
		return a tuple ( list of relations | list of properties )
		"""
		return ([],[])
	
	def fill (self, relations, properties) :
		"""
		construct internal links from basic python elements
		"""
		pass

class GraphRelation (TissueRelation,Graph) :
	"""
	specification of a relation object for a graph
	"""
	
	name = "graph"
	
	def __init__ (self, tissue, type_used, id_gen) :
		TissueRelation.__init__(self,tissue,type_used)
		Graph.__init__(self,idgenerator = id_gen)
		self._vertex_type = type_used[0]
		self._edge_type = type_used[1]
	
	def involved_elements (self) :
		yield self._vertex_type
		yield self._edge_type
	
	#######################################################
	#
	#		Mutable
	#
	#######################################################
	def add_element (self, elm_type, elm_id) :
		if elm_type==self._vertex_type :
			if not self.has_vertex(elm_id) :
				Graph.add_vertex(self,elm_id)
	
	def remove_element (self, elm_type, elm_id) :
		if elm_type == self._edge_type :
			if self.has_edge(elm_id) :
				Graph.remove_edge(self,elm_id)
		elif elm_type == self._vertex_type :
			Graph.remove_vertex(self,elm_id)
	#######################################################
	#
	#		Mutable Graph
	#
	#######################################################
	def add_vertex (self, vid=None) :
		return self.tissue().add_element(self._vertex_type,vid)
	
	def remove_vertex (self, vid) :
		self.tissue().remove_element(vid)
	
	def add_edge (self, sid, tid, eid=None) :
		eid = self.tissue().add_element(self._edge_type,eid,False)
		Graph.add_edge(self,sid,tid,eid)
		return eid
	
	def remove_edge (self, eid) :
		Graph.remove_edge(self,eid)
	#######################################################
	#
	#		Dict translation
	#
	#######################################################
	def description (self) :
		return (self.name,self._vertex_type,self._edge_type)
	
	def reduce (self) :
		rel = dict( (eid,(self.source(eid),self.target(eid))) for eid in self.edges() )
		return ([rel],[])
	
	def fill (self, relations, properties) :
		for eid,(sid,tid) in relations[0].iteritems() :
			Graph.add_edge(self,sid,tid,eid)

class RelationRelation (TissueRelation,Relation) :
	"""
	specification of a relation object for a relation
	"""
	
	name = "relation"
	
	def __init__ (self, tissue, type_used, id_gen) :
		TissueRelation.__init__(self,tissue,type_used)
		Relation.__init__(self,idgenerator = id_gen)
		self._left_type = type_used[0]
		self._right_type = type_used[1]
		self._link_type = type_used[2]
	
	def involved_elements (self) :
		yield self._left_type
		yield self._right_type
		yield self._link_type
	
	#######################################################
	#
	#		Mutable
	#
	#######################################################
	def add_element (self, elm_type, elm_id) :
		if elm_type == self._left_type :
			if not self.has_left(elm_id) :
				Relation.add_left_element(self,elm_id)
		elif elm_type == self._right_type :
			if not self.has_right(elm_id) :
				Relation.add_right_element(self,elm_id)
	
	def remove_element (self, elm_type, elm_id) :
		if elm_type == self._link_type :
			if self.has_link(elm_id) :
				Relation.remove_link(self,elm_id)
		elif elm_type == self._left_type :
			Relation.remove_left_element(self,elm_id)
		elif elm_type == self._right_type :
			Relation.remove_right_element(self,elm_id)
	#######################################################
	#
	#		Mutable Graph
	#
	#######################################################
	def add_left_element (self, elmid=None) :
		return self.tissue().add_element(self._left_type,elmid)
	
	def remove_left_element (self, elmid) :
		self.tissue().remove_element(elmid)
	
	def add_right_element (self, elmid=None) :
		return self.tissue().add_element(self._right_type,elmid)
	
	def remove_right_element (self, elmid) :
		self.tissue().remove_element(elmid)
	
	def add_link (self, left_elmid, right_elmid, lid=None) :
		lid = self.tissue().add_element(self._link_type,lid,False)
		Relation.add_link(self,left_elmid,right_elmid,lid)
		return lid
	
	def remove_link (self, lid) :
		Relation.remove_link(self,lid)
	#######################################################
	#
	#		Dict translation
	#
	#######################################################
	def description (self) :
		return (self.name,self._left_type,self._right_type,self._link_type)
	
	def reduce (self) :
		rel = dict( (lid,(self.left(lid),self.right(lid))) for lid in self.links() )
		return ([rel],[])
	
	def fill (self, relations, properties) :
		for lid,(left_elmid,right_elmid) in relations[0].iteritems() :
			Relation.add_link(self,left_elmid,right_elmid,lid)

class MeshRelation (TissueRelation,Topomesh) :
	"""
	specification of a relation object for a topomesh
	"""
	
	name = "mesh"
	
	def __init__ (self, tissue, type_used, id_gen) :
		TissueRelation.__init__(self,tissue,type_used)
		Topomesh.__init__(self,len(type_used)-1,idgenerator = id_gen)
		self._degree_type = type_used
	
	def involved_elements (self) :
		return iter(self._degree_type)
	
	#######################################################
	#
	#		Mutable
	#
	#######################################################
	def add_element (self, elm_type, elm_id) :
		for degree,typ in enumerate(self._degree_type) :
			if (typ==elm_type) and (not self.has_wisp(degree,elm_id)) :
				Topomesh.add_wisp(self,degree,elm_id)
	
	def remove_element (self, elm_type, elm_id) :
		for degree,typ in enumerate(self._degree_type) :
			if typ == elm_type :
				Topomesh.remove_wisp(self,degree,elm_id)
	##########################################################
	#
	#		Mutable Mesh
	#
	##########################################################
	def add_wisp (self, degree, wid=None) :
		return self.tissue().add_element(self._degree_type[degree],wid)
	
	def remove_wisp (self, degree, wid) :
		self.tissue().remove_element(wid)
	#######################################################
	#
	#		Dict translation
	#
	#######################################################
	def description (self) :
		return (self.name,)+tuple(self._degree_type)
	
	def reduce (self) :
		relations = []
		for d in xrange(1,self.degree()+1) :
			rel = {}
			key = 0
			for rid in self.wisps(d) :
				for bid in self.borders(d,rid) :
					rel[key] = (rid,bid)
					key += 1
			relations.append(rel)
		return (relations,[])
	
	def fill (self, relations, properties) :
		for d,rel in enumerate(relations) :
			for lid,(rid,bid) in rel.iteritems() :
				Topomesh.link(self,d+1,rid,bid)


