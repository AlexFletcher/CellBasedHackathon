import pickle, gzip
import numpy as np
from sa_oa.image.serial.basics import imread
from sa_oa.image.algo.analysis import SpatialImageAnalysis, DICT
from sa_oa.image.algo.graph_from_image import graph_from_image
from vplants.tissue_analysis import LienTissuTXT
from sa_oa.container import TemporalPropertyGraph

"""
basic script with functions to save TemporalPropertyGraph as tab delimited table with properties called for as column (named)
"""

def open_pkz_tpg(filename,path=""):
    f = gzip.open(path+filename,'r')
    graph=pickle.load(f)
    f.close()
    return graph


def write_txt(filename, vertex_properties, path=""):
    s="vid"+"\t"+"index"+"\t"+"cell_spi_id"+"\t"+"time_step"+"\t"
    s+="\t".join([p for p in properties])
    s+="\n"
    print s
    graph=open_pkz_tpg(filename)
    vids=sorted(list(graph.vertices()))
    for vid in vids:
        index = graph.vertex_property('index')[vid]
        cell_spi_id = graph.vertex_property('old_label')[vid]
        time_step = graph.graph_property('time_step')[index ]
        s+=str(vid)+"\t"+str(index)+"\t"+str(cell_spi_id)+"\t"+str(time_step)+"\t"
        for p in vertex_properties:
            print "Saving ",; "{}... ".format(p)
            if c in graph.vertex_property(p):
                s+=str(graph.vertex_property(p)[vid])+"\t"
            else:
                s+="\t"
        s=s[:-1]+"\n"

    f=open(path+filename+".txt","w")
    f.write(s)
    f.close()


def plot_mesures(filename, vertex_properties):
    #regarder le nombre de paramètres correspondant à telle catégorie de nombre de cellules
    graph=open_pkz_tpg(fleur)
    x=[]
    for p in vertex_properties:
        data=len([i for i in graph.vertex_property(p)])
        x.append(data)
    pylab.plot(np.sort(x))

### vertex_properties EXAMPLES
#vertex_properties=[ 'index',
 #'label',
 #'L1',
 #'volume_mean_neigh',
 #'Epidermis strain anisotropy',
 #'volume_growth',
 #'volume_growth_mean_neigh',
 #'Cell gaussian curvature','Epidermis strain anisotropy_mean_neigh',
 #'Cell gaussian curvature_mean_neigh',
 #'shape_anisotropy_mean_neigh',
 #'Division rate',
 #'Division rate_mean_neigh',
 #'volume',
 #'shape_anisotropy']
 
#vertex_properties=[ 'index',
 #'label',
 #'L1',
 #'volume_laplacian',
 #'Epidermis strain anisotropy',
 #'volume_growth',
 #'volume_growth_laplacian',
 #'Cell gaussian curvature','Epidermis strain anisotropy_laplacian',
 #'Cell gaussian curvature_laplacian',
 #'shape_anisotropy_laplacian',
 #'Division rate',
 #'Division rate_laplacian',
 #'volume',
 #'shape_anisotropy']
 
#vertex_properties=[ 'index',
# 'label',
# 'L1',
# "AHP6_signal",
# 'volume_mean_abs_dev',
# 'Epidermis strain anisotropy',
# 'volume_growth',
# 'volume_growth_mean_abs_dev',
# 'Cell gaussian curvature','Epidermis strain anisotropy_mean_abs_dev',
# 'Cell gaussian curvature_mean_abs_dev',
# 'shape_anisotropy_mean_abs_dev',
# 'Division rate',
# 'Division rate_mean_abs_dev',
# 'volume',
# 'shape_anisotropy']
 
#vertex_properties=[
 #'index',
 #'label',
 #'L1',
 #'Division rate',
 #'volume',
 #'Epidermis strain anisotropy',
 #'volume_growth',
 #'Cell gaussian curvature']
