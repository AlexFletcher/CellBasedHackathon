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
This module defines functions to serialize a tissue
"""

__license__= "Cecill-C"
__revision__=" $Id: $ "

from os.path import splitext
from zipfile import ZipFile,ZipInfo,ZIP_DEFLATED
import pickle
from time import localtime
from tissue import Tissue
from config import ConfigFile

class TissueFile (object) :
	"""
	object used to serialize a tissue
	"""
	def __init__ (self, filename, mode = 'r') :
		"""
		create a tissue file that will be stored
		in a file caled filename
		"""
		self._filename = filename
		self._mode = mode
		self._elms = {}
		if mode in ('r','a') :
			f = ZipFile(filename,'r')
			for name in f.namelist() :
				info = f.getinfo(name)
				data = f.read(name)
				self._elms[name] = (info,data)
			f.close()
	
	def close (self) :
		"""
		write everything inside the file
		and close it
		"""
		if self._mode in ('w','a') :
			f = ZipFile(self._filename,'w')
			for name,(info,data) in self._elms.iteritems() :
				f.writestr(info,data)
			f.close()
	
	def filenames (self) :
		"""
		return name of all files currently in the archive
		"""
		return self._elms.iterkeys()
	
	def read_file (self, filename) :
		"""
		read a file in the archive
		"""
		info,data = self._elms[filename]
		return data,info.comment
	
	def write_file (self, data, filename, description = "") :
		"""
		write a file into the archive
		"""
		info = ZipInfo(filename)
		info.comment = description
		info.date_time = localtime()[:6]
		info.external_attr = 0644 << 16L
		info.compress_type = ZIP_DEFLATED
		self._elms[filename] = (info,data)
	
	def remove_file (self, filename) :
		"""
		remove a file from the tissue archive
		"""
		del self._elms[filename]
	
	def chmod (self, filename, mod = 0644) :
		"""
		change file permission
		"""
		info = self._elms[filename][0]
		info.external_attr = mod << 16L
	
	def set_time (self, filename, datetime = None) :
		if datetime is None :
			datetime = localtime()[:6]
		info = self._elms[filename][0]
		info.date_time = datetime
	
	def set_compression (self, filename, comp = ZIP_DEFLATED) :
		info = self._elms[filename][0]
		info.compress_type = comp
	
	########################################
	#
	#	read methods
	#
	########################################
	def read_dict (self, filename) :
		data,descr = self.read_file(filename)
		d = pickle.loads(data)
		return d,descr
	
	def read_property (self, property_name) :
		return self.read_dict("%s.tip" % property_name)
	
	def properties (self) :
		"""
		iterator on all property names inside this tissue
		"""
		for filename in self.filenames() :
			propname,ext = splitext(filename)
			if ext == ".tip" :
				yield propname

	def read_tissue (self) :
		t=Tissue()
		#add types
		types,descr=self.read_dict("_types.tis")
		for type_id,name in types.iteritems() :
			t.add_type(name,type_id)
		#add elements
		elements,descr=self.read_dict("_elements.tis")
		for elm_id,type_id in elements.iteritems() :
			t.add_element(type_id,elm_id)
		#add relations
		relation_descr,d = self.read_dict("_relations.tis")
		for relation_id,(relation_args,link_files,prop_files) in relation_descr.iteritems() :
			relation = t.relation(t.add_relation(relation_args[0],relation_args[1:],relation_id))
			#fill relation with elements
			relations = []
			for filename in link_files :
				link,descr = self.read_dict(filename)
				relations.append(link)
			properties = []
			for filename in prop_files :
				prop,descr = self.read_dict(filename)
				properties.append(prop)
			relation.fill(relations,properties)
		return t,descr
	
	def read (self, name=None) :
		if name is None :
			return self.read_tissue()
		else :
			return self.read_property(name)
	
	def configs (self) :
		"""
		iterator on all config names inside this tissue
		"""
		for filename in self.filenames() :
			propname,ext = splitext(filename)
			if ext == ".xml" :
				yield propname

	def read_config (self, config_name) :
		filename="%s.xml" % config_name
		data,descr = self.read_file(filename)
		cfg_file = ConfigFile()
		return cfg_file.read(data)
	
	def external_files (self) :
		"""Iterator on all external files.
		
		i.e. not a tissue, not a property,
		not a config file.
		"""
		for filename in self.filenames() :
			propname,ext = splitext(filename)
			if ext not in (".tis",".tip",".xml") :
				yield filename
	
	########################################
	#
	#	write methods
	#
	########################################
	def write_dict (self, d, filename, description="") :
		data = pickle.dumps(d)
		self.write_file(data,filename,description)
	
	def write_property (self, prop, property_name, description="") :
		#prop = dict(prop.iteritems() )
		self.write_dict(prop,"%s.tip" % property_name,description)
	
	def write_tissue (self, tissue) :
		types = dict( (type_id,tissue.type_name(type_id)) for type_id in tissue.types() )
		self.write_dict(types,"_types.tis","name of element types")
		elements = dict( (elm_id,tissue.type(elm_id)) for elm_id in tissue.elements() )
		self.write_dict(elements,"_elements.tis","elements and their types")
		relation_descr = {}
		for relation_id in tissue.relations() :
			relation = tissue.relation(relation_id)
			relations,properties = relation.reduce()
			link_files = []
			for ind,rel in enumerate(relations) :
				filename = "_rel%d_link%d.tis" % (relation_id,ind)
				self.write_dict(rel,filename,"see file _relations.tis for a description of this relation")
				link_files.append(filename)
			prop_files = []
			for ind,prop in enumerate(properties) :
				filename = "_rel%d_prop%d.tis" % (relation_id,ind)
				self.write_dict(prop,filename,"see file _relations.tis for a description of this relation")
				prop_files.append(filename)
			relation_descr[relation_id]=(relation.description(),link_files,prop_files)
		self.write_dict(relation_descr,"_relations.tis","description of the relations between elements")
	
	def write (self, obj, name=None, description="") :
		if name is None :
			self.write_tissue(obj)
		else :
			self.write_property(obj,name,description)
	
	def write_config (self, config, config_name) :
		filename="%s.xml" % config_name
		cfg_file = ConfigFile(config)
		self.write_file(cfg_file.save(),filename,"config file")		

def topen (filename, mode='r') :
	return TissueFile(filename,mode)


