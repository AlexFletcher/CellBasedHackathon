import sa_oa.plantgl.math as mt
import sa_oa.plantgl.scenegraph as sg
import sa_oa.plantgl.algo as alg

def rgb2intensity(c) : return c.getRGBAverage()

class AscCodec (sg.SceneCodec):
    """ Ascii Point File Format 

    """

    def __init__(self):
        """
        Initialisation of the codec info
        """
        sg.SceneCodec.__init__(self,"ASC",sg.SceneCodec.Mode.ReadWrite)

    def formats(self):
        """ return formats """
        return [ sg.SceneFormat("Asc Codec",["asc","pts","xyz","pwn"],"The Ascii point file format") ]
        # pts format :
        #         first line : nb of points 
        #          then x y z w r g b
        # pwn format :
        #         first line : nb of points 
        #          then x y z
        # asc format :
        #          lines : x y z [r g b]
        # xyz format :
        #          lines : x y z
    def read(self,fname):
        """ read an ascii point file """
        import warnings
        pts = sg.Point3Array([])
        col = sg.Color4Array([])
        isptsfile = ('.pts' in fname)
        ispwnfile = ('.pwn' in fname)
        f = file(fname,"r")
        if isptsfile or ispwnfile:
            f.readline() # skip line number
        i = 0
        for line in f.readlines():
            if line[0] == '#': continue
            values = line.split()
            try:
                pts.append(mt.Vector3(float(values[0]),float(values[1]),float(values[2])))
                if len(values) > 3:
                    if not isptsfile : col.append(sg.Color4(int(values[3]),int(values[4]),int(values[5]),0))
                    else : col.append(sg.Color4(int(values[4]),int(values[5]),int(values[6]),0))
            except Exception,e:
                if isptsfile and len(values) == 1:
                    warnings.warn("Skip line "+str(i+isptsfile)+" in file '"+fname+"'.")
                else:
                    if len(values) != 0:
                        warnings.warn("Error in file '"+fname+"' at line "+str(i+isptsfile))
                        raise e
            i+=1
        f.close()
        center = pts.getCenter()
        if len(col) == len(pts) :
            pointset = sg.PointSet(pts,col)
        else:
            pointset = sg.PointSet(pts)
        pointset = sg.Translated(-center,pointset)
        sc = sg.Scene()
        sc += pointset
        return sc

    def write(self,fname,scene):
        """ write an ascii point file """
        print("Write "+fname)
        d = alg.Discretizer()
        f = file(fname,'w')
        isptsfile = ('.pts' in fname)
        ispwnfile = ('.pwn' in fname)
        isxyz = ('.xyz' in fname)
        if isptsfile or ispwnfile:
            nbpoints =  0
            for i in scene:
                if i.apply(d):
                    p = d.discretization
                    if isinstance(p,sg.PointSet) :
                        nbpoints += len(p.pointList)
            f.write(str(nbpoints)+'\n')
        for i in scene:
            if i.apply(d):
                p = d.discretization
                if isinstance(p,sg.PointSet) :
                    hasColor = not p.colorList is  None and len(p.colorList) > 0
                    col = i.appearance.ambient
                    for i,pt in enumerate(p.pointList):
                        f.write(str(pt.x)+' '+str(pt.y)+' '+str(pt.z)+' ')                        
                        if not isxyz and not ispwnfile:
                            if hasColor:                          
                                col = p.colorList[i]
                            if isptsfile:
                                f.write(str(rgb2intensity(col)))
                            f.write(str(col.red)+' '+str(col.green)+' '+str(col.blue)+'\n')
                        else: f.write('\n')
        f.close()
    

codec = AscCodec()
sg.SceneFactory.get().registerCodec(codec)
