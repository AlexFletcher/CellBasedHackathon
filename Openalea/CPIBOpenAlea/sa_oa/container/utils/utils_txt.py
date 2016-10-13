# -*- python -*-
# -*- coding: latin-1 -*-
#
#       utils : container package
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
This module provide a simple way to serialize elements in a txt file
"""

__license__= "Cecill-C"
__revision__=" $Id: utils_txt.py 7865 2010-02-08 18:04:39Z cokelaer $ "

def write_description (f, description) :
    """Add a description paragraph to a given stream
    if description is not empty
    """
    if description == "" :
        return
    f.write("BEGIN description\n")
    f.write(description)
    f.write("\n")
    f.write("END description\n")

def read_description (f) :
    """Read the description in a stream.
    """
    line = ""
    while "BEGIN description" not in line :
        line = f.readline()
    descr = ""
    line = f.readline()
    while "END description" not in line :
        descr += line
        line = f.readline()
    return descr[:-1]

def property_to_xml (doc, prop, description,
                                        key_type = "", key_description = "",
                                        value_type = "", value_description = "") :
    """
    return an xml description of this dictionary
    """
    #create base node
    prop_node = doc.createElement("property")
    #description
    write_description(doc,prop_node,description)

    #keys
    key_node = doc.createElement("keys")
    key_node.setAttribute("type",key_type)
    key_node.setAttribute("description",key_description)
    prop_node.appendChild(key_node)

    #values
    value_node = doc.createElement("values")
    value_node.setAttribute("type",value_type)
    value_node.setAttribute("description",value_description)
    prop_node.appendChild(value_node)

    #datas
    datas_node = doc.createElement("datas")
    for key,value in prop.iteritems() :
        node = doc.createElement("data")
        node.setAttribute("key",str(key))
        node.setAttribute("val",str(value))
        datas_node.appendChild(node)
    prop_node.appendChild(datas_node)

    #return
    return prop_node

def xml_to_property (property_node) :
    """
    retrieve the property structure from an xml node
    returns property,(description,(key_type,key_descr),(value_type,value_descr))
    """
    #create property
    prop = {}
    #read description
    descr = read_description(property_node)

    #keys
    key_node, = (node for node in property_node.childNodes if node.nodeName == "keys")
    key_type = key_node.getAttribute("type")
    key_descr = key_node.getAttribute("description")

    #values
    value_node, = (node for node in property_node.childNodes if node.nodeName == "values")
    value_type = value_node.getAttribute("type")
    value_descr = value_node.getAttribute("description")

    #datas
    datas_node, = (node for node in property_node.childNodes if node.nodeName == "datas")
    for node in datas_node.childNodes :
        if node.nodeName == "data" :
            key = eval(node.getAttribute("key"))
            value = eval(node.getAttribute("val"))
            prop[key] = value

    #return
    return prop,(descr,(key_type,key_descr),(value_type,value_descr))
