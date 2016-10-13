
from sa_vp.plantgl.scenegraph import Material
from sa_oa.tissueview import MeshView2D,MeshView
from sa_oa.pglviewer import  LoopView, LoopGUI, View3DGUI,\
                                TemplateGUI, Viewer, \
                                ViewerGUI, ColorScaleGUI
from model_output.view_db_gui import ViewDBGUI
from model_output.CustomScalarPropView import CustomScalarPropView
from model_output.CustomScalarPropView_edge import CustomScalarPropView_edge

from model_output.cellWallPropView import CellWallPropView2D
from model_output.edgePropView import EdgePropView2D

from sa_vp.plantgl.ext.color import JetMap
from PyQt4.QtCore import Qt, SIGNAL, SLOT, QObject


from model_utils.db_utilities import get_mesh

class SimViewer(QObject):
	def __init__(self, db):
		QObject.__init__(self)

		self.db=db
		self.cv=[]
		self.wv=[]

		mesh=get_mesh(db)
		pos=db.get_property('position')

		species_names = list(db.get_property('species_desc').iterkeys())
		species_names.sort(key=str.lower)
		self.view_names=[species_names[0]]#,'PIN7loc']
		self.cmaps=[JetMap(0,0.01, outside_values=True),JetMap(0.0,1.0, outside_values=True)]
			    
		sv=CustomScalarPropView(db,self.view_names[0],self.cmaps[0])
		sv._triangulation_method='topo'
		sv.cache_geometry()
		self.cv.append(sv)
		if mesh.degree()==3:        
			self.wv.append(MeshView(mesh, pos, 1, 0, Material((0,0,0))))
		if mesh.degree()==2:
			self.wv.append(MeshView2D(mesh, pos, 1, 0, Material((0,0,0)), 0.1))
		"""
		sv=CustomScalarPropView_edge(db,self.view_names[1],self.cmaps[1])
		sv._triangulation_method='topo'
		sv.cache_geometry()
		self.cv.append(sv)
		if mesh.degree()==3:        
			self.wv.append(MeshView(mesh, pos, 1, 0, Material((0,0,0))))
		if mesh.degree()==2:
			self.wv.append(MeshView2D(mesh, pos, 1, 0, Material((0,0,0)), 0.1))
		"""
	def redraw(self):		
		for i in self.cv:
			i.redraw()
		for i in self.wv:
			i.redraw()

			
	def create_windows(self, qapp, sch):
		loop = LoopView(sch)
		gui = TemplateGUI("Viewer")
		gui.setup_ui()

		v=[]
		
		v.append(Viewer())
		v[0].setWindowTitle(self.view_names[0])
		v[0].add_world(self.cv[0])
		if len(self.wv)>0:
			v[0].add_world(self.wv[0])
		v[0].add_gui(ViewerGUI(vars() ) )
	
		v[0].add_gui(gui)
		v[0].add_gui(LoopGUI(loop) )
		#v[0].add_gui(View3DGUI() )
		v[0].add_gui(ViewDBGUI(self.db,self.cv,v[0]) )
		v[0].show()
		mesh = get_mesh(self.db)
		v[0].view().set_dimension(mesh.degree())
		
		for i in range(1,len(self.view_names)):
			v.append(Viewer())
			v[i].setWindowTitle(self.view_names[i])
			v[i].add_world(self.cv[i])
			if len(self.wv)>0:
				v[i].add_world(self.wv[i])
			v[i].show()
			v[i].view().set_dimension(mesh.degree())
			v[0].synchronize(v[i])
		self.v=v


  
