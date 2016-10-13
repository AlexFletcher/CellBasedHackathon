from sa_oa.core import *


__name__ = "vplants.plantgl.visualization"
__alias__ = ["PlantGL.Visualization"]

__version__ = '0.0.2'
__license__ = 'CECILL-V2'
__authors__ = 'F. Boudon'
__institutes__ = 'INRIA/CIRAD'
__description__ = 'PlantGL Object Generator nodes.'
__url__ =  'http://sa_oa.gforge.inria.fr'

__all__ = ['plot3d']
        
plot3d = Factory( name= "plot3D", 
                  description= "Viewer Display", 
                  category = "Visualisation, plot", 
                  nodemodule = "sa_oa.plantgl.wralea.visualization.viewernode",
                  nodeclass = "Plot3D",
                  #lazy = False
                  )


