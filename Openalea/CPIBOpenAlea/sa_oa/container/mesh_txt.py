# -*- python -*-
# -*- coding: latin-1 -*-
#
#       Mesh : container package
#
#       Copyright or © or Copr. 2006 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <revesansparole@gmail.com>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       VPlants WebSite : https://gforge.inria.fr/projects/vplants/
#

__doc__="""
This module provide a simple way to serialize a mesh in a txt file
"""

__license__= "Cecill-C"
__revision__=" $Id: mesh_txt.py 14917 2013-09-27 12:28:55Z pradal $ "

__all__ = ["mesh_to_txt"
         , "mesh_from_txt"]


def mesh_to_txt (mesh) :
    """Construct a txt representation of the given mesh
    
    :Returns: (str)
    """
    txt = ""
    
    #write points
    txt += "BEGIN mesh points (pid: x, y, z, ...)\n"
    for pid in mesh.darts(0) :
        txt += "%d: " % pid
        txt += ", ".join([str(v) for v in mesh.position(pid)] )
        txt += "\n"
    
    txt += "END mesh points\n"
    
    #find borders
    borders = [(mesh.degree(did), did, tuple(mesh.borders(did) ) ) \
               for did in mesh.darts() if mesh.degree(did) > 0]
    borders.sort()
    
    #write borders
    txt += "BEGIN mesh borders (deg - did: bid1, ..., bidn)\n"
    for deg, did, bids in borders :
        txt += "%d - %d: " % (deg, did)
        txt += ", ".join(["%d" % bid for bid in bids])
        txt += "\n"
    
    txt += "END mesh borders\n"
    
    #return
    return txt


def mesh_from_txt (mesh, txt) :
    """Clear the mesh and fill it with description contained
    in txt
    """
    mesh.clear()
    
    lines = [line.strip() for line in txt.splitlines() \
             if len(line.strip() ) > 0]
    
    #find points
    line = lines.pop(0)
    while "BEGIN mesh points" not in line :
        if len(lines) == 0 :
            raise UserWarning("Not a valid mesh file")
        line = lines.pop(0)
    
    line = lines.pop(0)
    while "END mesh points" not in line :
        if len(lines) == 0 :
            raise UserWarning("Not a valid mesh file")
        
        pid_str, pos_str = line.split(": ")
        pid = mesh.add_dart(0, int(pid_str) )
        mesh.set_position(pid, tuple(float(s) for s in pos_str.split(", ") ) )
        
        line = lines.pop(0)

    #find border description
    line = lines.pop(0)
    while "BEGIN mesh borders" not in line :
        if len(lines) == 0 :
            raise UserWarning("Not a valid mesh file")
        line = lines.pop(0)

    line = lines.pop(0)
    while "END mesh borders" not in line :
        if len(lines) == 0 :
            raise UserWarning("Not a valid mesh file")
        
        did_str, bids_str = line.split(": ")
        deg, did = (int(s) for s in did_str.split(" - ") )
        did = mesh.add_dart(deg, did)
        for bid in [int(s) for s in bids_str.split(", ")] :
            mesh.link(did, bid)
        
        line = lines.pop(0)

















