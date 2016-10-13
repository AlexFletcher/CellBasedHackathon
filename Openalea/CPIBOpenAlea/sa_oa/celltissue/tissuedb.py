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
"""
Defines a database to store tissue and properties
"""

__license__= "Cecill-C"
__revision__=" $Id: $ "

from serial import topen

class TissueDB (object) :
	"""A simple container to maintain a tissue with a set
	of properties.
	"""
	
	def __init__ (self) :
		self._tissue = None
		self._property = {}
		self._descr = {}
		self._config = {}
		self._external = {}
	
	############################################
	#
	#	simple accessors
	#
	############################################
	def tissue (self) :
		"""Access to raw tissue structure.
		"""
		return self._tissue
	
	def set_tissue (self, tissue) :
		"""Set the tissue in this db.
		"""
		self._tissue = tissue
	
	def get_property (self, propname) :
		"""Access to a given property.
		"""
		return self._property[propname]
	
	def set_property (self, propname, prop) :
		"""Set the given property.
		"""
		self._property[propname] = prop
	
	def description (self, propname) :
		"""Retrieve the description associated
		to a given property.
		"""
		return self._descr[propname]
	
	def set_description (self, propname, descr) :
		"""Set the description associated
		to a given property.
		"""
		self._descr[propname] = descr
	
	def get_config (self, cfgname) :
		"""Access to a given config.
		"""
		return self._config[cfgname]
	
	def set_config (self, cfgname, cfg) :
		"""Store a cfg in this db.
		"""
		self._config[cfgname] = cfg
	
	def get_external_data (self, filename) :
		"""Access to data in external format.
		"""
		return self._external[filename]
	
	def set_external_data (self, filename, data) :
		"""Store any data in the tissue.
		"""
		self._external[filename] = data
	############################################
	#
	#	topology accessors
	#
	############################################
	def get_topology (self, toponame, cfg = "config") :
		"""Access to a given topology.
		
		toponame is the name of the topology as defined
		in the configuration file given by cfg
		"""
		return self._tissue.relation(self._config[cfg][toponame])
	############################################
	#
	#	inout
	#
	############################################
	def read (self, filename) :
		"""Fill this DB with info in filename.
		"""
		f = topen(filename,'r')
		#tissue
		self._tissue,descr = f.read()
		#properties
		self._property.clear()
		for propname in f.properties() :
			self._property[propname],self._descr[propname] = f.read(propname)
		#configs
		self._config.clear()
		for cfgname in f.configs() :
			self._config[cfgname] = f.read_config(cfgname)
		#external data
		self._external.clear()
		for filename in f.external_files() :
			self._external[filename],descr = f.read_file(filename)
		#close
		f.close()
	
	def write (self, filename) :
		"""Write this DB into a file.
		"""
		f = topen(filename,'w')
		#tissue
		f.write(self._tissue)
		#properties
		for propname,prop in self._property.iteritems() :
			f.write(prop,propname,self.description(propname) )
		#configs
		for cfgname,cfg in self._config.iteritems() :
			f.write_config(cfg,cfgname)
		#external data
		for filename,data in self._external.iteritems() :
			f.write_file(data,filename)
		#close
		f.close()

	def properties (self):
		return list(self._property.iterkeys())
