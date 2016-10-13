# -*- python -*-
# -*- coding: latin-1 -*-
#
#       Topomesh : container package
#
#       Copyright or © or Copr. 2006 INRIA - CIRAD - INRA
#
#       File author(s): Jerome Chopard <jerome.chopard@sophia.inria.fr>
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
__revision__=" $Id: topomesh_txt.py 14917 2013-09-27 12:28:55Z pradal $ "

from data_prop import DataProp
from topomesh import Topomesh
from utils.utils_txt import write_description,read_description

def topomesh_to_txt (f, mesh, description, props) :
    """Write the txt representation of a topomesh.
    
    f: an open file
    mesh: a topomesh object
    description: a textual description
    props: a list of property for each degree
    """
    #create base node
    f.write("BEGIN topomesh degree %d\n" % mesh.degree())
    f.write("\n")
    #description
    write_description(f,description)
    f.write("\n")
    
    #write properties
    f.write("BEGIN properties description\n")
    for deg,deg_props in enumerate(props) :
        f.write("deg%d" % deg)
        for name,prop in deg_props :
            f.write("\t%s(%s,%s)" % (name,prop.type(),prop.unit() ) )
        f.write("\n")
    f.write("END properties description\n")
    f.write("\n")

    #write wisps
    for deg in xrange(mesh.degree() + 1) :
        if deg < len(props) :
            deg_props = props[deg]
        else :
            deg_props = []
        f.write("BEGIN wisp degree %d\n" % deg)
        for wid in mesh.wisps(deg) :
            f.write("id %d" % wid)
            for name,prop in deg_props :
                f.write("\t%s" % prop[wid])
            f.write("\n")
        f.write("END wisp degree %d\n" % deg)

    #write links between wisps
    f.write("BEGIN decomposition\n")
    for deg in xrange(1,mesh.degree() + 1) :
        for wid in mesh.wisps(deg) :
            for bid in mesh.borders(deg,wid) :
                f.write("link degree %d wid %d bid %d\n" % (deg,wid,bid))
    f.write("END decomposition\n")

    f.write("END topomesh\n")

    #return
    return f

def txt_to_topomesh (f, method = "set") :
    """
    retrieve the topomesh structure from txt stream
    returns topomesh,description
    """
    #create topomesh
    line = ""
    while "BEGIN topomesh" not in line :
        line = f.readline()
    deg = int(line.split(" ")[3])
    mesh = Topomesh(deg,method)
    #read description
    descr = read_description(f)
    
    #read properties
    props = [[] for i in xrange(mesh.degree() + 1)]
    line = ""
    while "BEGIN properties description" not in line :
        line = f.readline()
    line = f.readline()
    while "END properties description" not in line :
        gr = line[:-1].split("\t")
        prop_deg = int(gr[0][3:])
        for prop_descr in gr[1:] :
            name,prop = prop_descr.split("(")
            typ,unit = prop[:-1].split(",")
            props[prop_deg].append( (name,DataProp(type = typ,unit = unit) ) )
        
        line = f.readline()

    #wisps
    for i in xrange(mesh.degree() + 1) :
        line = ""
        while "BEGIN wisp" not in line :
            line = f.readline()
        deg = int(line.split(" ")[3])
        line = f.readline().rstrip()
        while "END wisp" not in line :
            gr = line.split()
            wid = int(gr[1])
            mesh.add_wisp(deg,wid)
            #props
            for ind,val in enumerate(gr[2:]) :
                if props[deg][ind][1].type() == 'str' :
                    props[deg][ind][1][wid] = val
                else :
                    props[deg][ind][1][wid] = eval("%s(%s)" % (props[deg][ind][1].type(),val) )
            line = f.readline().rstrip()

    #links
    line = ""
    while "BEGIN decomposition" not in line :
        line = f.readline()
    line = f.readline()
    while "END decomposition" not in line :
        gr = line.split(" ")
        mesh.link(int(gr[2]),int(gr[4]),int(gr[6]))
        line = f.readline()

    #return
    return mesh,descr,props

def write_topomesh (filename, mesh, description, props = []) :
    """
    write a topomesh in a file using a txt representation
    """
    f = open(filename,'w')
    topomesh_to_txt(f,mesh,description,props)
    f.close()

def read_topomesh (filename, method = "set") :
    """
    read a topomesh stored in a file as a txt representation
    """
    f = open(filename,'r')
    mesh,descr,props = txt_to_topomesh(f,method)
    f.close()
    return mesh,descr,props

