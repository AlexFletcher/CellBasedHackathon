# -*- python -*-
#
#       tissueshape: function used to deal with tissue geometry
#
#       Copyright 2006 INRIA - CIRAD - INRA  
#
#       File author(s): Jerome Chopard <revesansparole@gmail.com>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
# 
#       OpenAlea WebSite : http://sa_oa.gforge.inria.fr
#

__doc__ = """
This module defines functions to create predefined meshes
with regular geometry
"""

__license__ = "Cecill-C"
__revision__ = " $Id: mesh_shapes.py 14917 2013-09-27 12:28:55Z pradal $ "

__all__ = ["line"
         , "polygon"
         , "regular_grid"]

from math import sqrt,sin,cos,pi
from grid import Grid
from mesh import Mesh

def line (nb_segments) :
    """Create a 1D mesh along a single line
    made of nb segments
    """
    #create mesh
    m = Mesh()
    
    #create points
    pids = []
    for i in range(nb_segments + 1) :
        pid = m.add_dart(0)
        m.set_position(pid, float(i) )
        pids.append(pid)
    
    #create segments
    for i in range(nb_segments) :
        did = m.add_dart(1)
        m.link(did, pids[i])
        m.link(did, pids[i + 1])
    
    #return
    return m


def polygon (nb_sides) :
    """Create a regular polygonal 2D mesh
    """
    #create mesh
    m = Mesh()
    
    #create points
    pids = []
    for i in range(nb_sides) :
        pid = m.add_dart(0)
        vec = (cos(2 * pi * i / nb_sides), sin(2 * pi * i / nb_sides) )
        m.set_position(pid, vec)
        pids.append(pid)
    
    #create segments
    eids = []
    for i in range(nb_sides) :
        eid = m.add_dart(1)
        m.link(eid, pids[i])
        m.link(eid, pids[(i + 1) % nb_sides])
        eids.append(eid)
    
    #create face
    fid = m.add_dart(2)
    for eid in eids :
        m.link(fid, eid)
    
    #return
    return m


def _regular_grid2D (shape) :
    cell_grid = Grid(shape)
    point_grid = Grid( (i + 1 for i in shape) )
    imax,jmax = shape
    
    #create mesh
    m = Mesh()
    
    #create cells
    cell = [m.add_dart(2) for ind in cell_grid]
    
    #create points
    point = [m.add_dart(0) for ind in point_grid]
    for i,pid in enumerate(point) :
        m.set_position(pid, point_grid.coordinates(i) )
    
    #create vertical walls
    for i in xrange(imax + 1) :
        for j in xrange(jmax) :
            wid = m.add_dart(1)
            m.link(wid, point[point_grid.index( (i,j) )])
            m.link(wid, point[point_grid.index( (i,j + 1) )])
            if i > 0 :
                m.link(cell[cell_grid.index( (i - 1,j) )], wid)
            if i < imax :
                m.link(cell[cell_grid.index( (i,j) )], wid)
    
    #create horizontal walls
    for i in xrange(imax) :
        for j in xrange(jmax + 1) :
            wid = m.add_dart(1)
            m.link(wid, point[point_grid.index( (i,j) )])
            m.link(wid, point[point_grid.index( (i + 1,j) )])
            if j > 0 :
                m.link(cell[cell_grid.index( (i,j - 1) )], wid)
            if j < jmax :
                m.link(cell[cell_grid.index( (i,j) )], wid)
    
    #return
    return m

def _regular_grid3D (shape) :
    cell_grid = Grid(shape)
    point_grid = Grid( (i + 1 for i in shape) )
    imax,jmax,kmax = shape
    
    #create mesh
    m = Mesh()
    
    #create cells
    cell = [m.add_dart(3) for ind in cell_grid]
    
    #create points
    point = [m.add_dart(0) for ind in point_grid]
    for i,pid in enumerate(point) :
        m.set_position(pid, point_grid.coordinates(i) )
    
    #create edges
    point_to_edge = {}
    #between i and i + 1
    for i in xrange(imax) :
        for j in xrange(jmax + 1) :
            for k in xrange(kmax + 1) :
                eid = m.add_dart(1)
                pid1 = point[point_grid.index( (i,j,k) )]
                pid2 = point[point_grid.index( (i + 1,j,k) )]
                m.link(eid, pid1)
                m.link(eid, pid2)
                key = (min(pid1, pid2), max(pid1, pid2) )
                point_to_edge[key] = eid
    #between j and j + 1
    for i in xrange(imax + 1) :
        for j in xrange(jmax) :
            for k in xrange(kmax + 1) :
                eid = m.add_dart(1)
                pid1 = point[point_grid.index( (i,j,k) )]
                pid2 = point[point_grid.index( (i,j + 1,k) )]
                m.link(eid, pid1)
                m.link(eid, pid2)
                key = (min(pid1, pid2), max(pid1,pid2) )
                point_to_edge[key] = eid
    #between k and k + 1
    for i in xrange(imax + 1) :
        for j in xrange(jmax + 1) :
            for k in xrange(kmax) :
                eid = m.add_dart(1)
                pid1 = point[point_grid.index( (i,j,k) )]
                pid2 = point[point_grid.index( (i,j,k + 1) )]
                m.link(eid, pid1)
                m.link(eid, pid2)
                key = (min(pid1, pid2), max(pid1,pid2) )
                point_to_edge[key] = eid
    
    #create walls
    #in the (i,j) plane, k constant
    for k in xrange(kmax + 1) :
        for i in xrange(imax) :
            for j in xrange(jmax) :
                wid = m.add_dart(2)
                pids = (point[point_grid.index( (i,j,k) )],
                        point[point_grid.index( (i + 1,j,k) )],
                        point[point_grid.index( (i + 1,j + 1,k) )],
                        point[point_grid.index( (i,j + 1,k) )] )
                for ind in xrange(4) :
                    pid1 = pids[ind]
                    pid2 = pids[(ind + 1) % 4]
                    eid = point_to_edge[(min(pid1, pid2), max(pid1, pid2) )]
                    m.link(wid, eid)
                if k > 0 :
                    m.link(cell[cell_grid.index( (i,j,k - 1) )], wid)
                if k < kmax :
                    m.link(cell[cell_grid.index( (i,j,k) )], wid)
    
    #in the (j,k) plane, i constant
    for i in xrange(imax + 1) :
        for j in xrange(jmax) :
            for k in xrange(kmax) :
                wid = m.add_dart(2)
                pids = (point[point_grid.index( (i,j,k) )],
                        point[point_grid.index( (i,j + 1,k) )],
                        point[point_grid.index( (i,j + 1,k + 1) )],
                        point[point_grid.index( (i,j,k + 1) )] )
                for ind in xrange(4) :
                    pid1 = pids[ind]
                    pid2 = pids[(ind + 1) % 4]
                    eid = point_to_edge[(min(pid1, pid2), max(pid1, pid2) )]
                    m.link(wid, eid)
                if i > 0 :
                    m.link(cell[cell_grid.index( (i - 1,j,k) )], wid)
                if i < imax :
                    m.link(cell[cell_grid.index( (i,j,k) )], wid)
    
    #in the (k,i) plane, j constant
    for j in xrange(jmax + 1) :
        for k in xrange(kmax) :
            for i in xrange(imax) :
                wid = m.add_dart(2)
                pids = (point[point_grid.index( (i,j,k) )],
                        point[point_grid.index( (i,j,k + 1) )],
                        point[point_grid.index( (i + 1,j,k + 1) )],
                        point[point_grid.index( (i + 1,j,k) )] )
                for ind in xrange(4) :
                    pid1 = pids[ind]
                    pid2 = pids[(ind + 1) % 4]
                    eid = point_to_edge[(min(pid1, pid2), max(pid1, pid2) )]
                    m.link(wid, eid)
                if j > 0 :
                    m.link(cell[cell_grid.index( (i,j - 1,k) )], wid)
                if j < jmax :
                    m.link(cell[cell_grid.index( (i,j,k) )], wid)
    
    #return
    return m

def regular_grid (shape) :
    """Create a mesh with a regular grid shape
    
    :Parameters:
     - `shape` (int or tuple of int) - nb of cells in each dimension
    
    :Returns: (Mesh)
    """
    if type(shape) == int :
        return line(shape)
    elif len(shape) == 1 :
        return line(shape[0])
    elif len(shape) == 2 :
        return _regular_grid2D(shape)
    elif len(shape) == 3 :
        return _regular_grid3D(shape)
    else :
        raise NotImplementedError("nD grid still in the box, %s" % str(shape) )


def _hexagonal_grid2D (shape, shape_geom = "hexa") :#TODO
    raise NotImplementedError
    '''
    cell_grid = Grid(shape)
    imax,jmax = shape
    assert imax > 1
    assert jmax > 1
    point_grid = Grid( (2 * imax + 1, 3 + (jmax - 1) ) )
    
    #create tissue
    t = Tissue()
    CELL = t.add_type("cell")
    WALL = t.add_type("wall")
    POINT = t.add_type("point")
    mesh_id = t.add_relation("mesh",(POINT,WALL,CELL) )
    mesh = t.relation(mesh_id)
    
    #create cells
    cell = [mesh.add_wisp(2) for ind in cell_grid]
    
    #create points
    point = [mesh.add_wisp(0) for ind in point_grid]
    
    #create horizontal walls
    for i in xrange(imax) :
        for j in xrange(jmax + 1) :
            wid = mesh.add_wisp(1)
            if j % 2 == 0 :
                mesh.link(1,wid,point[point_grid.index( (2 * i,j) )])
                mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j) )])
            else :
                mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j) )])
                mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j) )])
            
            if j > 1 :
                mesh.link(2,cell[cell_grid.index( (i,j - 2) )],wid)
            if j < jmax :
                mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
    
    j = jmax
    if j % 2 != 0 :
        for i in xrange(imax) :
            wid = mesh.add_wisp(1)
            mesh.link(1,wid,point[point_grid.index( (2 * i,j + 1) )])
            mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
            
            mesh.link(2,cell[cell_grid.index( (i,j - 1) )],wid)
    else :
        for i in xrange(imax) :
            wid = mesh.add_wisp(1)
            mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
            mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 1) )])        
            
            mesh.link(2,cell[cell_grid.index( (i,j - 1) )],wid)
    
    #create vertical walls
    for i in xrange(imax) :
        for j in xrange(jmax) :
            #bottom left wall
            wid = mesh.add_wisp(1)
            if j % 2 == 0 :
                mesh.link(1,wid,point[point_grid.index( (2 * i,j) )])
                mesh.link(1,wid,point[point_grid.index( (2 * i,j + 1) )])
            else :
                mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j) )])
                mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
            
            mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
            if j > 0 :
                if j % 2 == 0 :
                    if i > 0 :
                        mesh.link(2,cell[cell_grid.index( (i - 1,j - 1) )],wid)
                else :
                    mesh.link(2,cell[cell_grid.index( (i,j - 1) )],wid)
            #bottom right wall
            wid = mesh.add_wisp(1)
            if j % 2 == 0 :
                mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j) )])
                mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
            else :
                mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j) )])
                mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 1) )])
            
            mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
            if j > 0 :
                if j % 2 == 0 :
                    mesh.link(2,cell[cell_grid.index( (i,j - 1) )],wid)
                else :
                    if i < (imax - 1) :
                        mesh.link(2,cell[cell_grid.index( (i + 1,j - 1) )],wid)
    
    #top left, first column
    i = 0
    for j in (2 * k for k in xrange( (jmax + 1) / 2) ) :
        wid = mesh.add_wisp(1)
        mesh.link(1,wid,point[point_grid.index( (i,j + 1) )])
        mesh.link(1,wid,point[point_grid.index( (i,j + 2) )])
        mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
        
    #top right, last column
    i = imax - 1
    for j in (2 * k + 1 for k in xrange(jmax / 2) ) :
        wid = mesh.add_wisp(1)
        mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 1) )])
        mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 2) )])
        mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
    
    #top left and right for top line
    j = jmax - 1
    if j % 2 == 0 :
        #top left
        for i in xrange(1,imax) :
            wid = mesh.add_wisp(1)
            mesh.link(1,wid,point[point_grid.index( (2 * i,j + 1) )])
            mesh.link(1,wid,point[point_grid.index( (2 * i,j + 2) )])
            mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
        #top right
        for i in xrange(imax) :
            wid = mesh.add_wisp(1)
            mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
            mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 2) )])
            mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
    else :
        #top left
        for i in xrange(imax) :
            wid = mesh.add_wisp(1)
            mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 1) )])
            mesh.link(1,wid,point[point_grid.index( (2 * i + 1,j + 2) )])
            mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
        #top right
        for i in xrange(imax - 1) :
            wid = mesh.add_wisp(1)
            mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 1) )])
            mesh.link(1,wid,point[point_grid.index( (2 * i + 2,j + 2) )])
            mesh.link(2,cell[cell_grid.index( (i,j) )],wid)
        
    #create property position of points
    if shape_geom == 'box' :
        pos = dict( (pid,point_grid.coordinates(i) ) for i,pid in enumerate(point) )
    elif shape_geom == 'hexa' :
        pos = {}
        for ind,pid in enumerate(point) :
            i,j = point_grid.coordinates(ind)
            if j % 2 == 0 :
                x = 0.5 + (i / 2) * 3 + (i % 2)
                y = j * sqrt(3) / 2.
            else :
                x = - 1 + ( (i + 1) / 2) * 3 + ( (i + 1) % 2)
                y = j * sqrt(3) / 2.
            pos[pid] = (x,y)
    else :
        raise NotImplementedError("shape_geom: %s not recognized. only 'box' or 'hexa'" % str(shape_geom) )
    
    #create config
    cfg = Config("topology")
    cfg.add_item(ConfigItem("CELL",CELL) )
    cfg.add_item(ConfigItem("WALL",WALL) )
    cfg.add_item(ConfigItem("POINT",POINT) )
    cfg.add_item(ConfigItem("mesh_id",mesh_id) )
    
    #return
    return t,cfg,pos
    '''

def hexagonal_grid (shape, shape_geom = 'hexa') :#TODO
    """Create a tissue with a regular hexagonal grid topology
    
    The shape of each cell may either be a box or an hexagon
    
    :Parameters:
     - `shape` (int or tuple of int) - nb of cells in each dimension
     - `shape_geom` (str) : either 'hexa' or 'box'
    
    :Return: a tissuedb containing :
                - a mesh
                - a config
                - a pos property
    
    :Returns Type: TissueDB
    """
    raise NotImplementedError
    '''
    if type(shape) == int or len(shape) == 1 :
        t,cfg,pos = _regular_grid1D(shape)
    elif len(shape) == 2 :
        imax,jmax = shape
        if imax == 1 and jmax == 1 :
            t,cfg,pos = _single_cell2D(6)
        else :
            t,cfg,pos = _hexagonal_grid2D(shape,shape_geom)
    
    #create tissueDB
    db = TissueDB()
    db.set_tissue(t)
    db.set_property("position",Quantity(pos,"None") )
    db.set_description("position","geometrical position of points")
    db.set_config("config",cfg)
    
    return db
    '''





