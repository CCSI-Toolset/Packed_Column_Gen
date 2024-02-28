# ***************************************************************************
# *                                                                         *
# *   Copyright: See License.md file                                        *
# *                                                                         *
# *   Authors: Yash Girish Shah, Grigorios Panagakos                        *
# *                                                                         *
# ***************************************************************************

from .preamble import *

def getWorkbenchFolder():

    import os.path
    from os import path
    
    recommended_folders = (
    "Mod/packedColumnGen" ,
    "Mod/FreeCAD-packedColumnGen" 
    )
    
    basedir = str(FreeCAD.getUserAppDataDir())
    folder = ""
    
    for tryfolder in recommended_folders:
            if path.exists(basedir + tryfolder):
                    folder = basedir + tryfolder
                    return folder
    
                    
    return ""
        
def datum_settings(S):
    S.AttachmentOffset = App.Placement(App.Vector(0.0000000000, 0.0000000000, 0.0000000000),  App.Rotation(0.0000000000, 0.0000000000, 0.0000000000))
    S.MapReversed = False
    S.MapPathParameter = 0.000000

def convertDistanceUnits(l):
    return "%.6f mm"%(l)

def convertAngleUnits(a):
    return "%.6f deg"%(a)

def setDatum_distance(base,l):
    base.setDatum(base.ConstraintCount-1,App.Units.Quantity(convertDistanceUnits(l)))

def setDatum_angle(base,l):
    base.setDatum(base.ConstraintCount-1,App.Units.Quantity(convertAngleUnits(l)))

def line2vec(base,line_i = None):
    if line_i is None: line_i = base.GeometryCount - 1 
    return base.getPoint(line_i,2) - base.getPoint(line_i,1)

def setLength(base,l,line_i = None):
    if line_i is None: line_i = base.GeometryCount - 1 
    base.addConstraint(Sketcher.Constraint("Distance",line_i,l))
    setDatum_distance(base,l)

def horizontal_constraint(base,line_i = None):
    if line_i is None: line_i = base.GeometryCount - 1 
    base.addConstraint(Sketcher.Constraint("Horizontal",line_i))

def vertical_constraint(base,line_i = None):
    if line_i is None: line_i = base.GeometryCount - 1 
    base.addConstraint(Sketcher.Constraint("Vertical",line_i))

def rotateVector(vec,theta):
    v = vec[0] + 1j*vec[1]
    v *= np.exp(1j*theta)
    return App.Vector(np.real(v),np.imag(v),0.)

def connectedLineGen(base, i = 1, vec = None, line_i = None):
    if vec is None: vec = App.Vector(1.,0.,0.)
    if line_i is None: line_i = base.GeometryCount - 1 
    base.addGeometry(Part.LineSegment(base.getPoint(line_i,i),base.getPoint(line_i,i) + vec),False)
    line_j = base.GeometryCount - 1
    base.addConstraint(Sketcher.Constraint("Coincident",line_j,1,line_i,i))
    return base.GeometryCount - 1

def connectedLine1(base,vec = None, line_i = None):
    return connectedLineGen(base, i = 1, vec = vec, line_i = line_i)

def connectedLine2(base,vec = None, line_i = None):
    return connectedLineGen(base, i = 2, vec = vec, line_i = line_i)

def deg2rad(angle):
    return angle*np.pi/180.

def connecterLine(base,line_i,point_i,line_j,point_j):
    base.addGeometry(Part.LineSegment(base.getPoint(line_i,point_i),base.getPoint(line_j,point_j)),False)
    ii = base.GeometryCount - 1
    base.addConstraint(Sketcher.Constraint("Coincident",ii,1,line_i,point_i))
    base.addConstraint(Sketcher.Constraint("Coincident",ii,2,line_j,point_j))
    return ii

def setAngle(base,line_i,line_j,theta):
    # Sets counter-clockwise angle through distances. The angle is always the one obtained by rotating the line_i in the counter-clockwise direction till line_j is encountered.
    # NOTE: setAngle should always be run AFTER setting lengths of the lines involved.

    isect = None
    # identify common point and store in isect (i'nter'sect)
    for ii,jj in zip([1,1,2,2],[1,2,1,2]):
        vtmp = base.getPoint(line_i,ii) - base.getPoint(line_j,jj)
        if vtmp.Length < 1e-14: isect = (ii,jj)

    if isect is None:
        print("Cannot find intersection")
        sys.exit()

    #print("isect = ",isect)
    i0 = isect[0]; lst = [1,2]; lst.remove(i0); i1 = lst[0]
    j0 = isect[1]; lst = [1,2]; lst.remove(j0); j1 = lst[0]
    #print("outer points:",i1,j1)

    vi = base.getPoint(line_i,i1) - base.getPoint(line_i,i0)
    vj = base.getPoint(line_j,j1) - base.getPoint(line_j,j0)
    vsub = vi - vj

    isign = 1.
    if np.cos(theta) < 0.:
        isign = -1.

    if np.abs(np.cos(theta)) < 1e-16 :
        # Lines are orthogonal
        base.addConstraint(Sketcher.Constraint("Perpendicular",line_i,line_j))
    else:
        icons = connectedLineGen(base,i = i0,vec = isign*vi.normalize(),line_i = line_i)
        base.toggleConstruction(icons)
        base.addConstraint(Sketcher.Constraint("Parallel",icons,line_i))
        setLength(base,vj.Length*np.abs(np.cos(theta)),line_i = icons)

        vtmp = base.getPoint(line_j,j1) - base.getPoint(icons,2)

        jcons = connectedLineGen(base,i = 2,vec = vtmp,line_i = icons)
        base.toggleConstruction(jcons)
        setLength(base,vj.Length*np.abs(np.sin(theta)),line_i = jcons)
        #base.addConstraint(Sketcher.Constraint("Perpendicular",jcons,icons))
        base.addConstraint(Sketcher.Constraint("Coincident",jcons,2,line_j,j1))

def getLineNum(base):
    return base.GeometryCount - 1

def cell_outside_circle(bounding_box,vec,rcut):
    vec = np.array([vec[0],vec[1]])
    vlist = []
    for i,j in zip(["xmin","xmin","xmax","xmax"],["ymin","ymax","ymin","ymax"]):
        vlist.append(np.array([bounding_box[i]*1.,bounding_box[j]*1.]) + vec)

    outside = True
    for vi in vlist:
        if np.sum(vi**2.) < (rcut*2.)**2.:
            outside = False
    return outside

def within(x,x1,x2):
    lst = [x1,x2]
    x1 = np.min(lst)
    x2 = np.max(lst)
    bw = False
    if (x < x2) & (x > x1): bw = True
    return bw


def cell_outside_circle2(bounding_box,vec,rcut):
    vec = np.array([vec[0],vec[1]])
    xmin = bounding_box["xmin"] + vec[0]
    xmax = bounding_box["xmax"] + vec[0]
    ymin = bounding_box["ymin"] + vec[1]
    ymax = bounding_box["ymax"] + vec[1]
    vlist = []
    for i,j in zip(["xmin","xmin","xmax","xmax"],["ymin","ymax","ymin","ymax"]):
        vlist.append(np.array([bounding_box[i]*1.,bounding_box[j]*1.]) + vec)

    outside = True
    # Case 1: one of vlist is inside the circle
    #for vi in vlist:
    #    if np.sum(vi**2.) < (rcut)**2.:
    #        outside = False
    #        print("Case1")
    #        return outside

    # Case 2: Bbox is outside but one of the Bbox edges intersects the circle
    for xi in [xmin,xmax]:
        if within(xi,-rcut,rcut):
            yi = np.sqrt(rcut**2. - xi**2.)
            if within(ymin,-yi,yi) | within(ymax,-yi,yi):
                outside = False
                print("Case2")
                return outside

    for yi in [ymin,ymax]:
        if within(yi,-rcut,rcut):
            xi = np.sqrt(rcut**2. - yi**2.)
            if within(xmin,-xi,xi) | within(xmax,-xi,xi):
                outside = False
                print("Case2")
                return outside

    # Case 3: None of the Bbox edges intersect but the Bbox contains the circle
    if within(-rcut,xmin,xmax) & within(rcut,xmin,xmax):
        if within(-rcut,ymin,ymax) & within(rcut,ymin,ymax):
            print("Case3")
            outside = False
            return outside

    return outside

def cell_outside_circle1(original_vertices,vec,rcut):
    vlist = [vi + vec for vi in original_vertices]
    rmin = np.sqrt(np.min([vi[0]*vi[0] + vi[1]*vi[1] for vi in vlist]))
    outside = True
    if rmin <= rcut*1.1: outside = False
    return outside

def getFaceNormal(face):
    vlist = [i.Point for i in face.Vertexes]
    # Compute centroid at o
    o = App.Vector(0.,0.,0.)
    for v in vlist:
        o += v
    o /= len(vlist)
    if len(vlist) > 2:
        v1 = vlist[0] - o
        v2 = vlist[1] - o
        v = v1.cross(v2).normalize()
        return v

def checkFaceNormal(face,normal):
    match = False
    if (normal.Length != 0):
        normal = normal.normalize()
        v = getFaceNormal(face)
        if np.abs(np.abs(v.dot(normal)) - v.Length) < 1e-14:
            match = True
    return match

def getFaceListMatch(bdy,normal):
    faces = bdy.Shape.Faces
    flist = [i for i,f in enumerate(faces) if checkFaceNormal(f,normal)]
    return flist

def getLtransfYScale(pol_pat,lproj):
    flist = getFaceListMatch(pol_pat,App.Vector(0.,0.,1.))
    Amax = np.max([pol_pat.Shape.Faces[i].Area for i in flist])
    flist = [i for i in flist if pol_pat.Shape.Faces[i].Area/Amax > 0.001]
    zlist = np.array([pol_pat.Shape.Faces[i].Vertexes[0].Point[2] for i in flist])
    i1 = int(np.argmin(np.abs(zlist[1:] - zlist[0])) + 1)
    f1 = pol_pat.Shape.Faces[flist[0]]
    f2 = pol_pat.Shape.Faces[flist[i1]]
    vlist1 = [i.Point for i in f1.Vertexes]
    vlist2 = [i.Point for i in f2.Vertexes]
    d = -1e10
    for i,vi in enumerate(vlist1):
        for j,vj in enumerate(vlist2):
            dij = np.abs((vi - vj).dot(lproj))
            if dij > d:
                d = dij
                tup = (i,j)
    #print("Faces = %d,%d; tup = %d,%d "%(flist[0]+1,flist[i1]+1,tup[0],tup[1]))
    return d

def getNbodies(top,bot,rcut,tx,ty):
    Nx = []
    lx = tx[0]
    ly = tx[1]
    Nx.append((rcut - top["xmin"])/lx)
    Nx.append((top["xmax"] + rcut)/lx)
    Nx.append((rcut - bot["xmin"])/lx)
    Nx.append((bot["xmax"] + rcut)/lx)

    Nx.append((rcut - top["ymin"])/ly)
    Nx.append((top["ymax"] + rcut)/ly)
    Nx.append((rcut - bot["ymin"])/ly)
    Nx.append((bot["ymax"] + rcut)/ly)

    Ny = []
    lx = ty[0]
    ly = ty[1]
    Ny.append((rcut - top["xmin"])/lx)
    Ny.append((top["xmax"] + rcut)/lx)
    Ny.append((rcut - bot["xmin"])/lx)
    Ny.append((bot["xmax"] + rcut)/lx)
    Ny.append((rcut - top["ymin"])/ly)
    Ny.append((top["ymax"] + rcut)/ly)
    Ny.append((rcut - bot["ymin"])/ly)
    Ny.append((bot["ymax"] + rcut)/ly)

    Nx = int(np.ceil(np.max(Nx))) + 1
    Ny = int(np.ceil(np.max(Ny))) + 1
    return Nx,Ny



