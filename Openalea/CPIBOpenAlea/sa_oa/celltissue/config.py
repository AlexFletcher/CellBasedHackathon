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
This module defines a Config object to store informations in a nice way
"""

__license__= "Cecill-C"
__revision__=" $Id: $ "

from xml.dom.minidom import parseString,Document

class ConfigItem (object) :
	"""
	a container that contains a value
	"""
	def __init__ (self, name, value, unit = "", info = "") :
		self.name = name
		self.value = value
		self.unit = unit
		self.__doc__ = info

class Config (object) :
	"""
	a simple class to manage config files
	allow to access all properties defined into a python
	file using the dot operator
	"""
	def __init__ (self, name, items = [], info = "") :
		self.name = name
		self.elms = list(items)
		self.__doc__ = info
	
	def add_section (self, section) :
		self.elms.append(section)
		return len(self.elms)-1
	
	def add_item (self, item, section = None) :
		if section is None :#add item into the list
			self.elms.append(item)
		else :				#assume the section refers to a section
							#and add element to it
			self.elms[section].add_item(item)
		#add item into the local dict to access
		#it using the dot operator
		setattr(self,item.name,item.value)
	
	def __setattr__ (self, name, value) :
		object.__setattr__(self,name,value)
		if name not in ("name","elms","__doc__") :
			#update value in the tree too
			for elm in self.elms :
				if isinstance(elm,ConfigItem) :
					if elm.name == name :
						elm.value = value
						return
				else :
					elm.__setattr__(name,value)
	
	def __getitem__ (self, name) :
		"""
		access using dict interface
		"""
		if name not in ("name","elms","__doc__") :
			return getattr(self,name)
	
	def __iter__ (self) :
		"""
		emulation of list interface
		"""
		return iter(self.elms)

class ConfigFormat (object) :
	"""
	a simple class to help creating config files
	"""
	def __init__ (self, context, config = None) :
		self.context = context
		if config is None :
			self.cfg = Config("main")
		else :
			self.cfg = config
		self.current_section = None
	
	def add_section (self, title, info = "") :
		section = Config(title,[],info)
		self.current_section = self.cfg.add_section(section)
	
	def add_item (self, item) :
		self.cfg.add_item(item,self.current_section)
	
	def add (self, name, unit = "", info = "") :
		item = ConfigItem(name,self.context[name],unit,info)
		self.add_item(item)
	
	def config (self) :
		return self.cfg

class ConfigFile (object) :
	"""
	read and write config files in xml
	"""
	def __init__ (self, config = None) :
		self.cfg = config
		self.current_section = None
	
	def load_node (self, node) :
		"""
		load a specific node
		"""
		if node.nodeType == Document.ELEMENT_NODE :
			name = node.nodeName
			if name == "section" :
				title = node.attributes["a_name"].value
				info = node.attributes["b_info"].value
				#save current section
				section_mem = self.current_section
				#add a new section
				self.current_section = self.cfg.add_section(Config(title,[],info))
				#fill section with all elements
				for child in node.childNodes :
					self.load_node(child)
				#pop the previous section id
				self.current_section = section_mem
			elif name == "item" :
				title = node.attributes["a_name"].value
				unit = node.attributes["b_unit"].value
				info = node.attributes["c_info"].value
				value_node, = (child for child in node.childNodes if child.nodeType == Document.ELEMENT_NODE)
				pyval = value_node.attributes["value"].value
				val = eval(pyval)
				#create item
				item = ConfigItem(title,val,unit,info)
				self.cfg.add_item(item,self.current_section)
			else :
				raise UserWarning("xml element %s not recognized" % name)
	
	def load (self, data) :
		"""
		load a config object from its string representation
		"""
		doc = parseString(data)
		main, = (child for child in doc.childNodes if child.nodeType == Document.ELEMENT_NODE)
		for node in main.childNodes :
			self.load_node(node)
	
	def read (self, data) :
		"""
		return a config from the given xml representation
		"""
		self.cfg = Config("main")
		self.load(data)
		return self.cfg
	
	def save_node (self, obj, xmlparent) :
		"""
		create and return an xmlnode
		that contains all information in item
		"""
		doc = xmlparent.ownerDocument
		if isinstance(obj,Config) :
			#section
			node = doc.createElement("section")
			node.setAttribute("a_name",obj.name)
			node.setAttribute("b_info",obj.__doc__)
			#list of items
			for item in obj :
				self.save_node(item,node)
			#add node to parent
			xmlparent.appendChild(node)
		elif isinstance(obj,ConfigItem) :
			#item
			node = doc.createElement("item")
			node.setAttribute("a_name",obj.name)
			node.setAttribute("b_unit",obj.unit)
			node.setAttribute("c_info",obj.__doc__)
			valnode = doc.createElement("value")
			if type(obj.value) == str :
				strval = "'%s'" % str(obj.value)
			else :
				strval = str(obj.value)
			valnode.setAttribute("value",strval)
			node.appendChild(valnode)
			#add node to parent
			xmlparent.appendChild(node)
		else :
			raise UserWarning ("item of type %s undefined" % str(type(obj)))
	
	def save (self) :
		"""
		return a string representation of the config object
		"""
		doc = Document()
		node = doc.createElement(self.cfg.name)
		for elm in self.cfg :
			self.save_node(elm,node)
		doc.appendChild(node)
		return doc.toprettyxml()

