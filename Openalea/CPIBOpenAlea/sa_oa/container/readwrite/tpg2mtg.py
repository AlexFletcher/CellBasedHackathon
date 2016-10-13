# -*- python -
#
#       OpenAlea.Container
#
#       Copyright 2012 INRIA - CIRAD - INRA
#
#       File author(s):  Jonathan Legrand <jonathan.legrand@ens-lyon.fr>
#                        Frederic Boudon <frederic.boudon@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://sa_oa.gforge.inria.fr
#
################################################################################
"""This module helps to analyse TemporalPropertyGraph from Spatial Images."""

def tpg2mtg(graph, properties=['index'], root_ids=None, mtg = None, binary_sorter = None):
    """
    :Parameters:
     - `graph`(TPG) - a tempora property graph,
     - `properties`(list) - list of properties to export in MTG
     - `root_ids`(list) - list of root ids (@t0) to build the MTG
    """
    if root_ids is None:
        root_ids = [i for i in graph.vertices() if len(graph.parents(i)) == 0 and len(graph.children(i)) > 0]
            
    from sa_oa.mtg import MTG
    
    if mtg is None:
        mtg = MTG()
    
    mtg.add_property('tpg_cid')
    for p in properties:
        mtg.add_property(p)
    
    meristem_id = mtg.add_component(mtg.root, label='M')

    def translate_property(mtg, graph, mcid, tcid):
        mtg.property('tpg_cid')[mcid] = tcid
        for p in properties:
            try:
                mtg.property(p)[mcid] = graph.vertex_property(p)[tcid]
            except:
                continue
    
    def iterable(obj):
        try :
            iter(obj)
            return True
        except TypeError,te:
            return False

    def determine_parent(sortedchcids, mtg, parentmtgcid, parents):
        for chcid in sortedchcids:
            if iterable(chcid):
                # creating a virtual cell to have binary lineage
                mcid = mtg.add_child(parentmtgcid, label='C')
                for lchcid in chcid:
                    if iterable(lchcid) :
                        # if we have a new level of recursion, we apply the same algo.
                        determine_parent(lchcid, mtg, mcid, parents)
                    else:
                        parents[lchcid] = mcid
            else:
                parents[chcid] = parentmtgcid
        
    tcid2mcid = {}
    
    for trcid in root_ids:
        mcid = mtg.add_component(meristem_id, label='C')
        translate_property(mtg, graph, mcid, trcid)
        tcid2mcid[trcid] = mcid
        
        tcidstack = [trcid]
        
        while len(tcidstack) > 0:
            parentcid = tcidstack.pop(0)
            parentmtgcid = tcid2mcid[parentcid]
            chtcids = graph.children(parentcid)
            parents = {}
            
            # in this step, we check if we have a binary lineage. if not we create virtual cells.
            if len(chtcids) > 2 and binary_sorter:
                sortedchcids = binary_sorter(graph, chtcids)
                determine_parent(sortedchcids, mtg, parentmtgcid, parents)
            else:
                for chtcid in chtcids:
                    parents[chtcid] = parentmtgcid
            
            for chtcid in chtcids:
                mcid = mtg.add_child(parents[chtcid], label='C')
                translate_property(mtg, graph, mcid, chtcid)
                tcid2mcid[chtcid] = mcid
                tcidstack.append(chtcid)
    
    return mtg

def expert_sorter(graph, chtcids):
    """
    """
    assert 'sub_lineage' in graph.graph_properties()
    parent = list(graph.parent(list(chtcids)[0]))[0]
    return graph.graph_property('sub_lineage')[parent]


def auto_binary_sorter(graph, chtcids):
    """
    """
    # !!!!!!!!!!!!!!!!! the latest binary sibling should have the closest ppt !!!!!!!!!!!!!!!!!
    if len(chtcids) == 3:
        if sum([graph.vertex_property('L1').has_key(cid) for cid in chtcids])==3:
            return epidermis_surface_sorter(graph, chtcids)
        else:
            return volume_sorter(graph, chtcids)
    elif len(chtcids) == 4:
        return topological_sorter(graph, chtcids)
    else:
        warnings.warn("Not able to restore binary graph for more than 4 children...")
        return alea_sorter(graph, chtcids)


def volume_sorter(graph, chtcids):
    """
    Create the most probable binary lineage for a set of 3 siblings cells in a tissue.
    It is based on a comparison of the cells volume where we hypothesize a low difference in volumetric growth in between them.
    :Example:
     vol_A = 100 and vol_B = vol_C = 50, then the most probable lineage is:
     ABC -> A & BC -> A & B & C , in recursive list writting: [A,[B,C]]

    :Parameters:
     - `graph` (TPG) - a tempora property graph,
     - `chtcids` (list) - a list of sibling ids, ex. [A,B,C]

    :Returns:
     - `siblings` (list) - a recursive list indicating the deducted binary lineage, ex. [A,[B,C]]
    """
    assert len(chtcids) == 3
    volumes = dict([(vid, graph.vertex_property('volume')[vid]) for vid in chtcids])
    
    N = len(chtcids); other_ratios=[]
    ref_ratio = 0
    for n,vid in enumerate(chtcids):
        other_siblings = list(set(chtcids)-set([vid]))
        ratio = volumes[vid] / sum([volumes[cid] for cid in other_siblings])
        other_ratios.append(ratio)
        if ratio > ref_ratio:
            siblings = [vid,other_siblings]
            ref_ratio = ratio
            ref = n
        
    other_ratios.pop(ref)
    
    if ref_ratio < 0.7:
        print "Warning: the volume_sorter might be wrong in creating the most probable binary lineage for {}".format(siblings)
        print "Selected volume ratio = {}, against: {}".format(round(ref_ratio,3), [round(ratio,3) for ratio in other_ratios])
        if sum([graph.vertex_property('L1').has_key(cid) for cid in chtcids])==3:
            ep_sorted = epidermis_surface_sorter(graph, chtcids)
            if ep_sorted == siblings:
                print "Epidermis_surface_sorter return the same result..."
            else:
                print "Epidermis_surface_sorter does not return the same result: {}".format(ep_sorted)
            
    return siblings


def epidermis_surface_sorter(graph, chtcids):
    """
    Create the most probable binary lineage for a set of 3 siblings cells in a tissue.
    It is based on a comparison of the cells volume where we hypothesize a low difference in volumetric growth in between them.
    :Example:
     vol_A = 100 and vol_B = vol_C = 50, then the most probable lineage is:
     ABC -> A & BC -> A & B & C , in recursive list writting: [A,[B,C]]

    :Parameters:
     - `graph` (TPG) - a tempora property graph,
     - `chtcids` (list) - a list of sibling ids, ex. [A,B,C]

    :Returns:
     - `siblings` (list) - a recursive list indicating the deducted binary lineage, ex. [A,[B,C]]
    """
    assert len(chtcids) == 3
    assert sum([graph.vertex_property('L1').has_key(cid) for cid in chtcids])==3
    ep_surf = dict([(vid, graph.vertex_property('epidermis_surface')[vid]) for vid in chtcids])

    N = len(chtcids); other_ratios=[]
    ref_ratio = 0
    for n,vid in enumerate(chtcids):
        other_siblings = list(set(chtcids)-set([vid]))
        ratio = ep_surf[vid] / sum([ep_surf[cid] for cid in other_siblings])
        other_ratios.append(ratio)
        if ratio > ref_ratio:
            siblings = [vid,other_siblings]
            ref_ratio = ratio
            ref = n
        
    other_ratios.pop(ref)
    
    if ref_ratio < 0.7:
        print "Warning: the epidermis_surface_sorter might be wrong in creating the most probable binary lineage for {}".format(siblings)
        print "Selected surface ratio = {}, against: {}".format(round(ref_ratio,3), [round(ratio,3) for ratio in other_ratios])
    
    return siblings


def topological_sorter(graph, chtcids):
    """
    Create virtual intermediary cells to create a binary graph (one vertex can only give 2 new vertices)
    """
    reduced_neighborhood = dict([(vid, set(graph.neighbors(vid,'s'))&set(chtcids) ) for vid in chtcids])
    extern_nei = dict([(vid, set(graph.neighbors(vid,'s'))-set(chtcids) ) for vid in chtcids])

    not_adjacent = []
    for vid in chtcids:
        other_siblings = list(set(chtcids)-set([vid]))
        for cid in other_siblings:
            if vid not in reduced_neighborhood[cid]:
                not_adjacent.append([vid,cid])

    adjacent = []
    for vid in chtcids:
        other_siblings = list(set(chtcids)-set([vid]))
        for cid in other_siblings:
            # - Trying to find a common neighbor outside the set of siblings:
            if vid < cid and extern_nei[vid] & extern_nei[cid] != set([]):
                adjacent.append([vid,cid])
                if [vid,cid] in not_adjacent:
                    print "Warning: the topological_sorter might be wrong in creating the most probable binary lineage for {}".format(adjacent)

    return adjacent
