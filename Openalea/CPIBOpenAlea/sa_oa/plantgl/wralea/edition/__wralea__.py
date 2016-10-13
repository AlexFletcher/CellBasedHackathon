from sa_oa.core import *
from pgl_interface import *
from pgl_interface_widget import *

__name__ = "vplants.plantgl.edition"

__version__ = '0.0.2'
__license__ = 'CECILL-V2'
__authors__ = 'F. Boudon'
__institutes__ = 'INRIA/CIRAD'
__description__ = 'PlantGL Edition nodes.'
__url__ =  'http://sa_oa.gforge.inria.fr'

__all__ = ['curve2d', 'nurbs']

curve2d = Factory( name= "Curve2D", 
                  description= "Display and edition of a curve 2D", 
                  category = "Visualisation, Edition", 
                  nodemodule = "pgl_edition_node",
                  nodeclass = "curve2D",
                  widgetmodule = "pgl_interface_widget",
                  widgetclass = "Curve2DWidget",
                  inputs=(dict(name="curve", interface=None,),dict(name="curve", interface=None,)),
                  outputs=(dict(name="curve", interface=None,),),
                  lazy = False
                  )

nurbs= Factory( name= "NurbsPatch", 
                  description= "Display and edition of a Nurbs Patch", 
                  category = "Visualisation, Edition", 
                  nodemodule = "pgl_edition_node",
                  nodeclass = "nurbs",
                  widgetmodule = "pgl_interface_widget",
                  widgetclass = "NurbsPatchWidget",
                  inputs=(dict(name="nurbs", interface=None,),),
                  outputs=(dict(name="nurbs",),),
                  lazy = False
                  )





