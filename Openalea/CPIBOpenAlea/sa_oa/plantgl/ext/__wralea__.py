# -*- python -*-
#
#       OpenAlea.PlantGL.ext Library: PlantGL external Library module
#
#       Copyright 2006 INRIA - CIRAD - INRA  
#
#       File author(s): BOUDON Frederic <frederic.boudon@cirad.fr>
#                       Da SILVA David <david.da_silva@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
# 
#       OpenAlea WebSite : http://sa_oa.gforge.inria.fr
#


from sa_oa.core import *

__name__ = "vplants.plantgl.dresser"
__alias__ = ["PlantGL.Dresser"]

__version__ = '0.0.1'
__license__ = 'CECILL-C'
__authors__ = 'F. Boudon and D. Da Silva'
__institutes__ = 'INRIA/CIRAD'
__description__ = 'Spatial distribution module.'
__url__ = 'http://www.scipy.org'
    
__all__ = ["dresser"]    

dresser = Factory( name="PGL Dresser",
                   description="Add Pgl geometry an object",
                   category="scene.PGL",
                   nodemodule="dresser",
                   nodeclass="dresser",
                   )
                  

