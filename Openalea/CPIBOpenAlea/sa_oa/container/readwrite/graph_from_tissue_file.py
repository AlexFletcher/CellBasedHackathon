from sa_oa.image.algo.graph_from_image import *
import numpy as np

################################################################################
def graph_from_tissue(mesh, position):
    """

    """
    dim = mesh.degree()
    if not 1 < dim < 4 :
        print "bad dimension (expect 2 or 3)"
        return
        
    
    from sa_oa.tissueshape import centroid
    labels = [i for i in mesh.wisps(dim)]
    neighborhood=dict([(i,[j for j in mesh.border_neighbors(dim,i)]) for i in mesh.wisps(dim)])
    graph, label2vertex, edges = generate_graph_topology(labels, neighborhood)
    barycenters = None
    barycenters = dict( [(i, centroid(mesh, position, dim, i)) for i in labels] )
    add_vertex_property_from_dictionary(graph,'barycenter',barycenters)
    volumes = None
    
    from sa_oa.tissueshape import face_surface_2D, cell_volume
    if dim == 2:
        volumes = dict( [(i, face_surface_2D(mesh, position, i)) for i in labels] )
    else :
        volumes = dict( [(i, cell_volume(mesh, position, i)) for i in labels] )
    add_vertex_property_from_dictionary(graph,'volume',volumes)
    
    
    wall_vertices_coordinates = dict( [((min([t for t in mesh.regions(dim-1,w)]),
                                         max([t for t in mesh.regions(dim-1,w)])),
                                        [position[pt] for pt in [b for b in mesh.borders(dim-1,w,offset=dim-1)]]) 
                                       for w in mesh.wisps(dim-1) if len([i for i in mesh.regions(dim-1, w)])==2 ] )
    add_edge_property_from_label_property(graph,'wall_vertices_coordinates',wall_vertices_coordinates)
    
    
    from numpy.linalg import norm
    barycenter_distance = None
    barycenter_distance = dict( [((min([t for t in mesh.regions(dim-1,w)]),max([t for t in mesh.regions(dim-1,w)])), 
                                  norm(barycenters[[t for t in mesh.regions(dim-1,w)][0]] - barycenters[[t for t in mesh.regions(dim-1,w)][1]])) 
                                 for w in mesh.wisps(dim-1) if len([i for i in mesh.regions(dim-1, w)])==2] )
    add_edge_property_from_label_property(graph,'barycenter_distance',barycenter_distance)
    
    contact_size = None
    if dim == 2 :
        from sa_oa.tissueshape import edge_length
        contact_size = dict( [((min([t for t in mesh.regions(dim-1,w)]),max([t for t in mesh.regions(dim-1,w)])),
                               edge_length(mesh,position, w)) 
                              for w in mesh.wisps(dim-1) if len([i for i in mesh.regions(dim-1, w)])==2] )
    else :
        from sa_oa.tissueshape import face_surface_3D
        contact_size = dict( [((min([t for t in mesh.regions(dim-1,w)]),max([t for t in mesh.regions(dim-1,w)])), 
                               face_surface_3D(mesh,position, w)) 
                              for w in mesh.wisps(dim-1) if len([i for i in mesh.regions(dim-1, w)])==2] )
    add_edge_property_from_label_property(graph,'contact_size',contact_size)
    
    
    return graph



################################################################################
def graph_from_svg(filename):
    """
    """
    from vplants.plantgl.math import Vector2
    from sa_oa.container import Topomesh
    from sa_oa.svgdraw import open_svg, SVGSphere, SVGConnector
    from sa_oa.tissueshape import edge_loop_around
    
    mesh = Topomesh(2)
    f = open_svg(filename,'r')
    sc = f.read()
    f.close()
    
    zone_pos = {}
    zone_svg_id = {}
    lay = sc.get_layer("zones")
    for elm in lay.elements() :
        if isinstance(elm,SVGSphere) :
            cid = mesh.add_wisp(2)
            zone_svg_id[elm.id()] = cid
            zone_pos[cid] = Vector2(*sc.natural_pos(*elm.scene_pos(elm.center() ) ) )
    
    vertex_pos = {}
    vertex_svg_id = {}
    lay = sc.get_layer("walls")
    for elm in lay.elements() :
        if isinstance(elm,SVGSphere) :
            pid = mesh.add_wisp(0)
            vertex_svg_id[elm.id()] = pid
            vertex_pos[pid] = Vector2(*sc.natural_pos(*elm.scene_pos(elm.center() ) ) )
    
    lay = sc.get_layer("walls")
    for elm in lay.elements() :
        if isinstance(elm,SVGConnector) :
            eid = mesh.add_wisp(1)
            mesh.link(1,eid,vertex_svg_id[elm.source()])
            mesh.link(1,eid,vertex_svg_id[elm.target()])
    for cid,ref_point in zone_pos.iteritems() :
        for eid in edge_loop_around(mesh,vertex_pos,ref_point) :
            mesh.link(2,cid,eid)
    
    position = dict([ (v,np.array(p)) for v,p in vertex_pos.iteritems() ])
    
    return graph_from_tissue(mesh, position)


################################################################################
def graph_from_TissueDB(filename, meshname="mesh_id", posname="position"):
    """
    """
    from sa_oa.celltissue import TissueDB
    db = TissueDB()
    db.read(filename)
    mesh = db.get_topology(meshname)
    pos = db.get_property(posname)
    position = dict([(i,np.array(j)) for i,j in pos.iteritems()])
    
    return graph_from_tissue(mesh, position)


#################################################
#def graph_from_image(image, 
#                     labels = None, 
#                     background = 1, 
#                     default_properties = default_properties,
#                     default_real_property = True,
#                     bbox_as_real = False,
#                     remove_stack_margins_cells = True):
#    """ 
#        Construct a PropertyGraph from a SpatialImage (or equivalent) representing a segmented image.

#        :Parameters:
#         - `labels` (list) - sequence of label numbers of the objects to be measured.
#            If labels is None, all labels are used.
#         - `background` (int) - label representing background.
#         - `default_properties` (list) - the list of name of properties to create. It should be in default_properties.
#         - `default_real_property` (bool) - If default_real_property = True, property is in real-world units else in voxels.
#         - `bbox_as_real` (bool) - If bbox_as_real = True, bounding boxes are in real-world units else in voxels.

#        :rtype: PropertyGraph

#        :Examples:

#        >>> import numpy as np
#        >>> image = np.array([[1, 2, 7, 7, 1, 1],
#                          [1, 6, 5, 7, 3, 3],
#                          [2, 2, 1, 7, 3, 3],
#                          [1, 1, 1, 4, 1, 1]])

#        >>> from sa_oa.image.algo.graph_from_image import graph_from_image
#        >>> graph = graph_from_image(image)

#    """

#    analysis = SpatialImageAnalysis(image)
#    if remove_stack_margins_cells:
#        analysis.remove_margins_cells()

#    if labels is None: 
#        filter_label = False
#        labels = list(analysis.labels())
#        if background in labels : del labels[labels.index(background)]
#        neigborhood = analysis.neighbors(labels)
#    else:
#        filter_label = True
#        if isinstance(labels,int) : labels = [labels]
#        # -- We don't want to have the "outer cell" (background) and "removed cells" (0) in the graph structure.
#        if 0 in labels: labels.remove(0)
#        if background in labels: labels.remove(background)
#        neigborhood = analysis.neighbors(labels)

#    labelset = set(labels)

#    graph, label2vertex, edges = generate_graph_topology(labels, neigborhood)

#    if 'boundingbox' in default_properties : 
#        add_vertex_property_from_label_and_value(graph,'boundingbox',labels,analysis.boundingbox(labels,real=bbox_as_real),mlabel2vertex=label2vertex)

#    if 'volume' in default_properties : 
#        add_vertex_property_from_dictionary(graph,'volume',analysis.volume(labels,real=default_real_property),mlabel2vertex=label2vertex)

#    barycenters = None
#    if 'barycenter' in default_properties :
#        barycenters = analysis.center_of_mass(labels,real=default_real_property)
#        add_vertex_property_from_label_and_value(graph,'barycenter',labels,barycenters,mlabel2vertex=label2vertex)

#    background_neighbors = set(analysis.neighbors(background))
#    background_neighbors.intersection_update(labelset)
#    if 'L1' in default_properties :         
#        add_vertex_property_from_label_and_value(graph,'L1',labels,[(l in background_neighbors) for l in labels],mlabel2vertex=label2vertex)

#    if 'border' in default_properties : 
#        border_cells = analysis.border_cells()
#        try: border_cells.remove(background)
#        except: pass
#        border_cells = set(border_cells)
#        add_vertex_property_from_label_and_value(graph,'border',labels,[(l in border_cells) for l in labels],mlabel2vertex=label2vertex)

#    if 'inertia_axis' in default_properties : 
#        inertia_axis, inertia_values = analysis.inertia_axis(labels,barycenters)
#        add_vertex_property_from_label_and_value(graph,'inertia_axis',labels,zip(inertia_axis,inertia_values),mlabel2vertex=label2vertex)

#    if 'wall_surface' in default_properties : 
#        filtered_edges = {}
#        for source,targets in neigborhood.iteritems():
#            if source in labelset :
#                filtered_edges[source] = [ target for target in targets if source < target and target in labelset ]
#        wall_surfaces = analysis.wall_surfaces(filtered_edges,real=default_real_property)
#        add_edge_property_from_label_property(graph,'wall_surface',wall_surfaces,mlabelpair2edge=edges)

#    if 'epidermis_surface' in default_properties :
#        def not_background(indices):
#            a,b = indices
#            if a == background: 
#                if b == background: raise ValueError(indices)
#                else : return b
#            elif b == background: return a
#            else: raise ValueError(indices)
#        epidermis_surfaces = analysis.cell_wall_surface(background,list(background_neighbors) ,real=default_real_property)
#        epidermis_surfaces = dict([(not_background(indices),value) for indices,value in epidermis_surfaces.iteritems()])
#        add_vertex_property_from_label_property(graph,'epidermis_surface',epidermis_surfaces,mlabel2vertex=label2vertex)

#    return graph       
