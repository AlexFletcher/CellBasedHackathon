from os import path
from PyQt4.QtCore import Qt,SIGNAL,QObject
from PyQt4.QtGui import QAction,QActionGroup,QToolBar,\
                        QMenu,QIcon,QFileDialog,\
                        QColorDialog,QPixmap,QColor,\
                        QStatusBar,QLineEdit
from PyQGLViewer import Camera
from sa_oa.pglviewer.constants import VIEW,ROTATION_ANCHOR
from sa_oa.pglviewer.elm_gui import ElmGUI
import sa_oa.pglviewer.pyshell.pyshell_rc

#######################################################
#
#	view
#
#######################################################


class ActionSaveDB (QAction) :
	"""Take a snapshot of the current frame.
	"""
	def __init__ (self, parent,db) :
		QAction.__init__(self,"savedb",parent)
		self.setIcon(QIcon(":image/snapshot") )
		self.setShortcut("Ctrl+D")
		self.connect(self,SIGNAL("triggered(bool)"),self.savedb)
		self.db=db
	
	def set_view (self, view) :
		self._view = view
	
	def savedb (self) :
		filename = 'test_savedb.zip'
		self.db.write(filename)


class ActionRenorm (QAction) :
	"""Take a snapshot of the current frame.
	"""
	def __init__ (self, parent,db,cv) :
		QAction.__init__(self,"renormalise colour map",parent)
		self.setIcon(QIcon(":image/snapshot") )
		self.setShortcut("Ctrl+R")
		self.connect(self,SIGNAL("triggered(bool)"),self.renorm)
		self.db=db
		self.cv=cv
	
	def set_view (self, view) :
		self._view = view
	
	def renorm (self) :

		#self.cv[0].change_property(self.db.get_property('cytokinin'))
		self.cv[0].renormalise()		
		self.cv[0].redraw()


class ActionChangeProp (QAction) :
	
	def __init__ (self, parent,db,cv,propname,vwindow) :
		QAction.__init__(self,propname,parent)
		self._propname = propname
		self.setCheckable(True)

		self.connect(self,SIGNAL("triggered(bool)"),self.change_prop)
		self._db=db
		self._cv=cv
		if propname==cv[0]._prop_name:
			self.setChecked(True)
		self._vwindow=vwindow
	def set_view (self, view) :
		self._view = view
	
	def change_prop (self) :

		self._cv[0].change_property(self._db.get_property(self._propname),self._propname)
		self._cv[0].renormalise()
		self._vwindow.setWindowTitle(self._propname)		
		self._cv[0].redraw()
		

		
#######################################################
#
#	GUI
#
#######################################################
class ViewDBGUI (ElmGUI) :
	"""Some standard GUI for a 3D view.
	"""
	def __init__ (self,db,cv,vwindow) :
		ElmGUI.__init__(self)
		self._db=db
		self._cv=cv
		self._vwindow=vwindow
	def setup_ui (self) :
		if ElmGUI.setup_ui(self) :
			self._action_bar = QToolBar("viewDB")
			self._menu = QMenu("DB")
			#view
			self._action_savedb = ActionSaveDB(self._action_bar,self._db)
			self._action_bar.addAction(self._action_savedb)
			self._menu.addAction(self._action_savedb)
			self._action_renorm = ActionRenorm(self._action_bar,self._db,self._cv)
			self._action_bar.addAction(self._action_renorm)
			self._menu.addAction(self._action_renorm)
			
			#cell types menu
			self._menu.addSeparator()
			self._menu_types = self._menu.addMenu("types")
			

			def change_view_types (active_types) :
				for i in self._cv:
					i.cache_geometry_subset(self._db.get_property('cell_type'),self.active_types)
					i.redraw()
			
			def check_types(checked):
				self.active_types=[]
				for ctype,ac in self.acdict.iteritems():
					if ac.isChecked():
						self.active_types.append(ctype)
				change_view_types(self.active_types)
						
			
			ctypes=self._db.get_property('cell_type_strings')
			self.acdict={}
			for code,ctype in ctypes.iteritems() :
				ac = self._menu_types.addAction(ctype)
				ac.setCheckable(True)
				ac.setChecked(True)
				QObject.connect(ac,
				                SIGNAL("toggled(bool)"),
				                check_types)
				self.acdict[code]=ac
			
			self._action_bar.addAction(self._menu_types.menuAction() )


			#property menu
			self._menu_props = self._menu.addMenu("property")
			self._actions_view_prop = []
			self._prop_group = QActionGroup(self._action_bar)
			species_desc=self._db.get_property('species_desc')
			for propname,ntype in species_desc.iteritems():
				if ntype=='CELL' or ntype=='cell':
					action = ActionChangeProp(self._action_bar,self._db,self._cv,propname,self._vwindow)
					self._menu_props.addAction(action)
					self._prop_group.addAction(action)
					self._actions_view_prop.append(action)
			self._action_bar.addAction(self._menu_props.menuAction() )

			
			#v index menu
			self._menu.addSeparator()
			self._menu_vindex = self._menu.addMenu("v index")
			

			def change_vindex_types (vtypes) :
				for i in self._cv:
					i.cache_geometry_subset(self._db.get_property('vindex'),self.v_types)
					i.redraw()
			
			def check_vtypes(checked):
				self.v_types=[]
				for ctype,ac in self.vindict.iteritems():
					if ac.isChecked():
						self.v_types.append(ctype)
				change_vindex_types(self.v_types)
						
			
			vtypes=range(min(list(self._db.get_property('vindex').itervalues())),max(list(self._db.get_property('vindex').itervalues()))+1)
			#vtypes=range(0,2)
			
			self.vindict={}
			for code in vtypes:
				ac = self._menu_vindex.addAction(str(code))
				ac.setCheckable(True)
				ac.setChecked(True)
				QObject.connect(ac,
				                SIGNAL("toggled(bool)"),
				                check_vtypes)
				self.vindict[code]=ac
			
			self._action_bar.addAction(self._menu_vindex.menuAction() )
			

			#cursor

			
			self._xcoord = QLineEdit("")
			self._xcoord.setMaximumWidth(50)
			self._xcoord.setAlignment(Qt.AlignRight)
			
			self._ycoord = QLineEdit("")
			self._ycoord.setMaximumWidth(50)
			self._ycoord.setAlignment(Qt.AlignRight)
			
			self._zcoord = QLineEdit("")
			self._zcoord.setMaximumWidth(50)
			self._zcoord.setAlignment(Qt.AlignRight)
			
			QObject.connect(self._xcoord,
			                SIGNAL("editingFinished()"),
			                self.change_cursor_position)
			QObject.connect(self._ycoord,
			                SIGNAL("editingFinished()"),
			                self.change_cursor_position)
			QObject.connect(self._zcoord,
			                SIGNAL("editingFinished()"),
			                self.change_cursor_position)



			return True
		else :
			return False
	
	def install (self, main_window) :
		ElmGUI.install(self,main_window)
		#install
		main_window.menuBar().addMenu(self._menu)
		self.add_action_bar(main_window,self._action_bar)
		#register view
		self._action_savedb.set_view(main_window.view() )
		self._action_renorm.set_view(main_window.view() )


		self.cursor_position_changed(main_window.view().cursor().position() )
		QObject.connect(main_window.view().cursor(),
		                SIGNAL("set_position"),
		                self.cursor_position_changed)
		
		#status bar
		self.add_status_widget(main_window,self._xcoord,True)
		self.add_status_widget(main_window,self._ycoord,True)
		self.add_status_widget(main_window,self._zcoord,True)


	def uninstall (self, main_window) :
		ElmGUI.uninstall(self,main_window)
		main_window.menuBar().removeAction(self._menu.menuAction() )
		self.remove_bar(main_window,self._action_bar)
		#discard view
		self._action_savedb.set_view(None)
		self._action_renorm.set_view(None)

		QObject.disconnect(main_window.view().cursor(),
		                   SIGNAL("set_position"),
		                   self.cursor_position_changed)
		
		#status bar
		self.remove_status_widget(main_window,self._xcoord)
		self.remove_status_widget(main_window,self._ycoord)
		self.remove_status_widget(main_window,self._zcoord)



	def change_cursor_position (self) :
		"""Change the cursor position
		according to information in statusbar
		"""
		main_window = self.installed()
		if main_window is not None :
			x = float(str(self._xcoord.text() ) )
			y = float(str(self._ycoord.text() ) )
			z = float(str(self._zcoord.text() ) )
			main_window.view().cursor().set_position( (x,y,z) )
	
	def cursor_position_changed (self, pos) :
		"""Display the new cursor position
		"""
		self._xcoord.setText("%.2f" % pos.x)
		self._ycoord.setText("%.2f" % pos.y)
		self._zcoord.setText("%.2f" % pos.z)

