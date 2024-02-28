# -*- coding: utf-8 -*-
# Column code v1: Created Feb 02, 2023
# Code for intensified and non-intensified columns

import FreeCAD
import PartDesign
import PartDesignGui
import Sketcher
import Part
import numpy as np
import sys
#import File
import ImportGui
import Mesh

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
    base.addConstraint(constraint("Distance",line_i,l))
    setDatum_distance(base,l)

def horizontal_constraint(base,line_i = None):
    if line_i is None: line_i = base.GeometryCount - 1 
    base.addConstraint(constraint("Horizontal",line_i))

def vertical_constraint(base,line_i = None):
    if line_i is None: line_i = base.GeometryCount - 1 
    base.addConstraint(constraint("Vertical",line_i))

def rotateVector(vec,theta):
    v = vec[0] + 1j*vec[1]
    v *= np.exp(1j*theta)
    return vector(np.real(v),np.imag(v),0.)

def connectedLineGen(base, i = 1, vec = None, line_i = None):
    if vec is None: vec = vector(1.,0.,0.)
    if line_i is None: line_i = base.GeometryCount - 1 
    base.addGeometry(line(base.getPoint(line_i,i),base.getPoint(line_i,i) + vec),False)
    line_j = base.GeometryCount - 1
    base.addConstraint(constraint("Coincident",line_j,1,line_i,i))
    return base.GeometryCount - 1

def connectedLine1(base,vec = None, line_i = None):
    return connectedLineGen(base, i = 1, vec = vec, line_i = line_i)

def connectedLine2(base,vec = None, line_i = None):
    return connectedLineGen(base, i = 2, vec = vec, line_i = line_i)

def deg2rad(angle):
    return angle*np.pi/180.

def connecterLine(base,line_i,point_i,line_j,point_j):
    base.addGeometry(line(base.getPoint(line_i,point_i),base.getPoint(line_j,point_j)),False)
    ii = base.GeometryCount - 1
    base.addConstraint(constraint("Coincident",ii,1,line_i,point_i))
    base.addConstraint(constraint("Coincident",ii,2,line_j,point_j))
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
        base.addConstraint(constraint("Perpendicular",line_i,line_j))
    else:
        icons = connectedLineGen(base,i = i0,vec = isign*vi.normalize(),line_i = line_i)
        base.toggleConstruction(icons)
        base.addConstraint(constraint("Parallel",icons,line_i))
        setLength(base,vj.Length*np.abs(np.cos(theta)),line_i = icons)

        vtmp = base.getPoint(line_j,j1) - base.getPoint(icons,2)

        jcons = connectedLineGen(base,i = 2,vec = vtmp,line_i = icons)
        base.toggleConstruction(jcons)
        setLength(base,vj.Length*np.abs(np.sin(theta)),line_i = jcons)
        #base.addConstraint(constraint("Perpendicular",jcons,icons))
        base.addConstraint(constraint("Coincident",jcons,2,line_j,j1))

def getLineNum(base):
    return base.GeometryCount - 1

def cell_outside_circle(bounding_box,vec):
    vec = np.array([vec[0],vec[1]])
    vlist = []
    for i,j in zip(["xmin","xmin","xmax","xmax"],["ymin","ymax","ymin","ymax"]):
        vlist.append(np.array(bounding_box[i]*1.,bounding_box[j]*1.) + vec)

    outside = True
    for vi in vlist:
        if np.sum(vi**2.) < (rcut*2)**2.:
            outside = False
    return outside

def checkFaceNormal(face,normal):
    match = False
    if (normal.Length != 0):
        normal = normal.normalize()
        vlist = [i.Point for i in face.Vertexes]
        o = vector(0.,0.,0.)
        for v in vlist:
            o = o + v
        o = o/len(vlist)
        if len(vlist) > 2:
            v1 = vlist[0] - o
            v2 = vlist[1] - o
            v = v1.cross(v2)
            if np.abs(np.abs(v.dot(normal)) - v.Length) < 1e-14:
                match = True
    return match

def getFaceListMatch(bdy,normal):
    faces = bdy.Shape.Faces
    flist = [i for i,f in enumerate(faces) if checkFaceNormal(f,normal)]
    return flist

def getLtransfYScale(pol_pat,lproj):
    flist = getFaceListMatch(pol_pat,vector(0.,0.,1.))
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
    print("Faces = %d,%d; tup = %d,%d "%(flist[0]+1,flist[i1]+1,tup[0],tup[1]))
    return d
            

        

# ========= User Parameters ===============
export_dir  = u"D:/OneDrive/NETL/FreeCAD_demo/3Dparametrization_PCCC/"
exportName = "B90G60" #"M250Y"
dName       = "unit_cell"
H           = 14.75 #M327Y equivalent
#H           = 19.58 # M250Y equivalent exact
#alpha_angle = 45       # 45. #Alpha_angle is a dependent parameter
beta_angle  = 60 #67.06 #89.      # 90.
gamma_angle = 60.      # 45.
t1          = 0.890955*0.33 # 0.890955
t2          = 0.18*0.33     # 0.18
h           = 0.891155*0.33      # 0.891155
rcut        = 38.1 #38.1 # 1.5 inch
rcol        = rcut #+ 6.*t1
NbodiesX = 7
NbodiesY = 7
NbodiesZ = 4
tcol     = t1*1. # Column thickness
hcol     = 100 # Column height
nInlet   = 4
rInlet   = 2
hInlet   = 15 - 8 # Inlet dripping points depth
column_Z_offset = 17.5 #+ 10
column_Z_offset_bottom = 20
packing_height = hcol - column_Z_offset - column_Z_offset_bottom - hInlet  #210
#packing_height = 14.78*NbodiesZ
print("packing height = ",packing_height)
lgap = 0
cut = True
alpha_angle = 90. - beta_angle/2. # <-- Parallelity
intensified_device = False
derive_Nz_from_packing_height = True #True
axis_overlap = 4e-3 # 4e-5
modify_hcol = False
apply_fusion = False
packing_gap_toleranceZ = -0.055 #-0.055
# =========================================

t3          = t1/2.*np.tan(alpha_angle*np.pi/180.)  #0.45     # 0.72

try:
    App.closeDocument(dName)
except:
    print("Opening new document")

App.newDocument(dName)
dmt = App.activeDocument()
dmt.addObject("PartDesign::Body","Body")


E = App.ActiveDocument
F = Gui.ActiveDocument
E.recompute()

#Gui.runCommand("PartDesign_NewSketch",0)
Gui.activateWorkbench("PartDesignWorkbench")

G = Gui.getDocument(dName)
A = App.getDocument(dName)

#G.ActiveView.setActiveObject("pdbody",A.getObject("Body"),"")

body = A.getObject("Body")
body.newObject("Sketcher::SketchObject","Base")


base = A.getObject("Base")
base.Support = (A.getObject("XY_Plane"),[""])
base.MapMode = "FlatFace"
E.recompute()

#G.setEdit(body,0,"Base")
vector = App.Vector
constraint = Sketcher.Constraint
line = Part.LineSegment
hvec = vector(1.,0.,0.)
vvec = vector(0.,1.,0.)

lines = []

# Line -- 0
base.addGeometry(line(vector(0.,0.,0.),vector(-1.,1.,0)),False)
lines.append(getLineNum(base))
base.addConstraint(constraint("Coincident",lines[0],1,-1,1))

s1 = np.abs((2*t1 + h)/np.cos(deg2rad(alpha_angle)))
print(convertDistanceUnits(s1))
setLength(base,s1)
E.recompute()

# Line -- 1
lines.append(connectedLine2(base,vec = vector(-1,0,0)*H))
horizontal_constraint(base)
setLength(base,H)
setAngle(base,lines[1],lines[0],deg2rad(alpha_angle + 90.))
E.recompute()

# Line -- 2
x1 = np.abs(base.getPoint(lines[1],2)[0] + (2*t1+h)/np.tan(deg2rad(beta_angle/2.)))
base.addGeometry(line(vector(0.,0.,0.),vector(-1.,0.,0.)*x1),False)
lines.append(getLineNum(base))
base.addConstraint(constraint("Coincident",lines[2],1,lines[0],1))
base.addConstraint(constraint("Parallel",lines[2],lines[1]))
E.recompute()


# Line -- 3
lines.append(connecterLine(base,lines[1],2,lines[2],2))
setLength(base,(2*t1+h)/np.sin(deg2rad(beta_angle/2.)))
base.toggleConstruction(lines[3])
E.recompute()


# Line -- 4
dirvec = line2vec(base,lines[3])
dirvec = dirvec.normalize()
l1 = np.abs(t1/np.sin(deg2rad(beta_angle/2.)))
print(convertDistanceUnits(l1))
pt1 = base.getPoint(lines[1],2) + dirvec*l1

base.addGeometry(line(pt1 ,pt1 + hvec),False)
lines.append(getLineNum(base))
base.addConstraint(constraint("PointOnObject",lines[4],1,lines[3]))
base.addConstraint(constraint("Parallel",lines[4],lines[1]))
if not intensified_device:
    base.toggleConstruction(lines[4]) # Toggled construction to remove cooling channels 
E.recompute()

# Line -- 5
pt2 = base.getPoint(lines[2],2) - dirvec*l1
base.addGeometry(line(pt2,pt2 + hvec),False)
lines.append(getLineNum(base))
base.addConstraint(constraint("PointOnObject",lines[5],1,lines[3]))
base.addConstraint(constraint("Parallel",lines[5],lines[4]))
if not intensified_device:
    base.toggleConstruction(lines[5]) # Toggled construction to remove cooling channels 
E.recompute()


# Line -- 6
lines.append(connecterLine(base,lines[4],2,lines[5],2))
base.addConstraint(constraint('Parallel',lines[6],lines[0])) 
if not intensified_device:
    base.toggleConstruction(lines[6]) # Toggled construction to remove cooling channels 
E.recompute()


# Line -- 7
lines.append(connecterLine(base,lines[1],2,lines[4],1))
setLength(base,l1)
base.toggleConstruction(lines[7])
E.recompute()

# Line -- 8
lines.append(connecterLine(base,lines[5],1,lines[2],2))
base.toggleConstruction(lines[8])
base.addConstraint(constraint("Equal",lines[7],lines[8]))
E.recompute()


# Line -- 9
lines.append(connectedLine2(base,hvec,lines[4]))
base.addConstraint(constraint("PointOnObject",lines[9],2,lines[0]))
base.addConstraint(constraint("Perpendicular",lines[0],lines[9]))
base.toggleConstruction(lines[-1])
setLength(base,t2)

E.recompute()


# Line -- 10
dirvec = -1.*line2vec(base,lines[1])
dirvec = rotateVector(dirvec,-deg2rad(beta_angle)).normalize()
lines.append(connectedLine2(base,dirvec*H,lines[1]))
base.addConstraint(constraint("Equal",lines[-1],lines[1]))
setAngle(base,lines[-1],lines[1],deg2rad(beta_angle))
E.recompute()

# Line -- 11
lines.append(connectedLine2(base,dirvec,lines[2]))
base.addConstraint(constraint('Equal',lines[-1],lines[2])) 
base.addConstraint(constraint('Parallel',lines[-1],lines[10])) 
E.recompute()

# Line -- 12
lines.append(connecterLine(base,lines[10],2,lines[11],2))
E.recompute()


# Line -- 13
lth = line2vec(base,lines[4]).Length
lines.append(connectedLine1(base,dirvec*lth,lines[4]))
base.addConstraint(constraint("Equal",lines[-1],lines[4]))
base.addConstraint(constraint('Parallel',lines[-1],lines[10])) 
#setAngle(base,lines[-1],lines[4],deg2rad(beta_angle))
if not intensified_device:
    base.toggleConstruction(lines[-1]) # Toggled construction to remove cooling channels 
E.recompute()

# Line -- 14
lth = line2vec(base,lines[5]).Length
lines.append(connectedLine1(base,dirvec*lth,lines[5]))
base.addConstraint(constraint("Equal",lines[-1],lines[5]))
base.addConstraint(constraint('Parallel',lines[-1],lines[10])) 
#setAngle(base,lines[-1],lines[5],deg2rad(beta_angle))
if not intensified_device:
    base.toggleConstruction(lines[-1]) # Toggled construction to remove cooling channels 
E.recompute()

# Line -- 15
lines.append(connecterLine(base,lines[13],2,lines[14],2))
if not intensified_device:
    base.toggleConstruction(lines[-1]) # Toggled construction to remove cooling channels 
E.recompute()
E.recompute()
E.recompute()

lsweep = (base.getPoint(lines[2],1) - base.getPoint(lines[11],2)).Length
print("l1 = ",convertDistanceUnits(l1))
L = l1/np.abs(np.sin(deg2rad(gamma_angle)))

if derive_Nz_from_packing_height:
    ii = int(packing_height/lsweep + 1. - 1e-5)
    if ii != NbodiesZ:
        if modify_hcol: hcol += np.max([0.,(ii - NbodiesZ)*lsweep])
        NbodiesZ = ii
        print("Modifying NbodiesZ. NbodiesZ = %d, hcol = %.3g"%(NbodiesZ,hcol))

# Datum line -- Diagonal

body.newObject("PartDesign::Line","Diagonal")
diagonal_datum = A.getObject("Diagonal")
#datum_settings(diagonal_datum)
diagonal_datum.Support = [(base,"Vertex1"),(base,"Vertex5")]
diagonal_datum.MapMode = 'TwoPointLine'
diagonal_datum.Visibility = False
E.recompute()

body.newObject("PartDesign::Plane","DiagPlane")
diagPlane = A.getObject("DiagPlane")
#datum_settings(diagPlane)
diagPlane.Support = [(diagonal_datum,"")]
diagPlane.MapMode = 'ObjectYZ'
diagPlane.Visibility = False
E.recompute()


body.newObject("Sketcher::SketchObject","sweepProfile")
sweepProfile = A.getObject("sweepProfile")
sweepProfile.Support = (diagPlane,[""])
sweepProfile.MapMode = "FlatFace"
E.recompute()


# Line -- 0
sweepProfile.addGeometry(line(vector(0.,0.,0.),vector(-1.,-1.,0.)*L),False)
sweepProfile.addConstraint(constraint("Coincident",0,1,-1,1))
lsweep1 = packing_height*1.
setLength(sweepProfile,lsweep1/np.cos(deg2rad(gamma_angle)))
E.recompute()


# Line -- 1 (Construction)
connectedLine2(sweepProfile,vector(0.,lsweep1,0.),0)
sweepProfile.toggleConstruction(1)
sweepProfile.addConstraint(constraint("PointOnObject",1,2,-1))
sweepProfile.addConstraint(constraint("Perpendicular",-1,1))
setLength(sweepProfile,lsweep1)
E.recompute()



# 3D Sweep
#Gui.runCommand("PartDesign_AdditivePipe",0)
body.newObject("PartDesign::AdditivePipe","HalfCell")
halfCell = A.getObject("HalfCell")
#G.setEdit(body,0,"HalfCell")

halfCell.Profile = base
halfCell.Spine = sweepProfile
base.Visibility = False
sweepProfile.Visibility = False
E.recompute()


body.newObject("Sketcher::SketchObject","oSk")
oSk = A.getObject("oSk")
oSk.Support = (halfCell,["Face8"])
oSk.MapMode = "FlatFace"
E.recompute()


# Line -- 0
oSk.addGeometry(line(vector(0.,0.,0.),vector(-1.,0.,0.)),False)
oSk.addConstraint(constraint("Coincident",0,1,-1,1))
oSk.addConstraint(constraint("PointOnObject",0,2,-1))
oSk.addConstraint(constraint("Distance",0,t3))
oSk.setDatum(2,App.Units.Quantity(convertDistanceUnits(t3)))
E.recompute()
oSk.Visibility = False


body.newObject("PartDesign::Plane","YZof")
YZof = A.getObject("YZof")
#datum_settings(YZof)
YZof.Support = [(oSk,"Vertex2"),(A.getObject("Y_Axis"),""),(A.getObject("Z_Axis"),"")]
YZof.MapMode = 'OXY'
E.recompute()



body.newObject('Sketcher::SketchObject','dum')
dum = A.getObject("dum")
#dum.Support = (YZof,['']) # latest mod <--- Just Checking
dum.Support = (A.getObject("YZ_Plane"),[""])
dum.MapMode = 'FlatFace'
E.recompute()

dum.addGeometry(Part.LineSegment(vector(0.000000,0.000000,0),vector(0.,-1.,0)),False)
#dum.addConstraint(constraint("Coincident",0,1,-1,1))
dum.addConstraint(constraint("PointOnObject",0,1,-1)) # <--
dum.addConstraint(constraint("Vertical",0))
setLength(dum,1.)
#dum.addConstraint(constraint("DistanceX",-1,1,0,1,t1/2.)) #<--
dum.addConstraint(constraint("DistanceX",-1,1,0,1,axis_overlap)) #<--
#dum.setDatum(3,App.Units.Quantity(convertDistanceUnits(t1/2.))) #<--
dum.setDatum(3,App.Units.Quantity(convertDistanceUnits(axis_overlap))) #<--
E.recompute()

YZof.Visibility = False

body.newObject('PartDesign::ShapeBinder','Rotation_axis')
sbd = A.getObject("Rotation_axis")
sbd.Support = [(dum,"Edge1")]
E.recompute()

dum.Visibility = False

body.newObject('PartDesign::PolarPattern','PolarPattern')
pol_pat = A.getObject("PolarPattern")
E.recompute()
pol_pat.Originals = [halfCell,]
pol_pat.Axis = (sbd,['Edge1'])
pol_pat.Angle = 360
pol_pat.Occurrences = 2
body.Tip = pol_pat
E.recompute()

sbd.Visibility = False

bounding_box = {}
vlist = [i.Point for i in pol_pat.Shape.Vertexes]
xlist = np.array([i[0] for i in vlist])
ylist = np.array([i[1] for i in vlist])
zlist = np.array([i[2] for i in vlist])
bounding_box["xmin"] = xlist.min()
bounding_box["ymin"] = ylist.min()
bounding_box["zmin"] = zlist.min()
bounding_box["xmax"] = xlist.max()
bounding_box["ymax"] = ylist.max()
bounding_box["zmax"] = zlist.max()
print(bounding_box)


# Lid covers

try:
    shp = Part.getShape(body,"",needSubElement=True,refine=True)
except:
    shp = Part.getShape(body,"",needSubElement=True,refine=False)

lTransfX = (base.getPoint(lines[2],1) - base.getPoint(lines[11],2))
lTransfX *= (lTransfX.Length + lgap)/lTransfX.Length
lTransfY = lTransfX.cross(vector(0,0,-1))
lTransfY /= lTransfY.Length
lTransfZ = sweepProfile.getPoint(0,2) - sweepProfile.getPoint(0,1)
lTransfZ = (lTransfZ[0] + packing_gap_toleranceZ)*vector(0,0,1) + lTransfX*np.abs(lTransfZ[1])/lTransfX.Length

#if intensified_device:
#    ll = np.abs((pol_pat.Shape.Vertexes[2].Point - pol_pat.Shape.Vertexes[53].Point).dot(lTransfY)) + lgap # Projection of the point to point distance on the vertical translation vector direction
#else:
#    ll = np.abs((pol_pat.Shape.Vertexes[2].Point - pol_pat.Shape.Vertexes[35].Point).dot(lTransfY)) + lgap # Projection of the point to point distance on the vertical translation vector direction
ll = getLtransfYScale(pol_pat,lTransfY)
lTransfY *= ll

# Test translate

#A.addObject("Part::Feature","test").Shape = shp
#A.getObject("test").Placement = App.Placement(lTransfX*0 + lTransfY*1+lTransfZ*0,App.Rotation(vector(0,0,1),0),vector(0,0,0))
#sys.exit()


# Cut Packing in a cylinder
A.addObject('PartDesign::Body','BodyCut')
body_cut = A.BodyCut

body_cut.newObject("Sketcher::SketchObject","BaseCut")
base_cut = A.getObject("BaseCut")
base_cut.Support = (A.getObject("XY_Plane001"),[""])
base_cut.MapMode = "FlatFace"
base_cut.Visibility = False
E.recompute()


# Line -- 0 (Circle)
base_cut.addGeometry(Part.Circle(vector(0.,0.,0.),vector(0.,0.,1.),rcut),False)
base_cut.addConstraint(constraint("Coincident",0,3,-1,1))
base_cut.addConstraint(constraint("Radius",0,rcut))
base_cut.setDatum(1,App.Units.Quantity(convertDistanceUnits(rcut)))
E.recompute()

body_cut.newObject('PartDesign::Pad',"cut_pad")
cut_pad = A.getObject("cut_pad")
cut_pad.Profile = base_cut
cut_pad.Length = lsweep*NbodiesZ
if derive_Nz_from_packing_height:
    cut_pad.Length = packing_height
cut_pad.Reversed = True
E.recompute()

shp_cut = Part.getShape(body_cut,"",needSubElement=False,refine=False)

commons = []
labs0 = [i.Label for i in A.findObjects() if "Common" in i.Label]
A.addObject("Part::MultiCommon","Common")
labs = [i.Label for i in A.findObjects() if "Common" in i.Label]
labs = [i for i in labs if i not in labs0]
commons.append(A.getObject(labs[0]))
commons[-1].Shapes = [body_cut,body]

E.recompute()
vertical_faces_lst = [getFaceListMatch(commons[-1],vector(0.,0.,1.))]
coords = [(0,0,0)]


f = open(export_dir+u"log.txt","w")
f.writelines(["started \n"])

bd = []
bd_cut = []
matrix = {}
#for n in range(0,NbodiesZ):
for n in [0]:
    for k in range(-NbodiesY,NbodiesY):
        for j in range(-NbodiesX-n,NbodiesX-n):
            if np.abs(j) + np.abs(k) + np.abs(n) == 0:
                continue
            if cell_outside_circle(bounding_box,lTransfX*j + lTransfY*k+lTransfZ*n):
                continue
            bdname = "bd_%d_%d_%d"%(j,k,n)
            labs0 = [i.Label for i in A.findObjects() if "Body" in i.Label]
            A.addObject("Part::Feature","Body").Shape = shp
            labs = [i.Label for i in A.findObjects() if "Body" in i.Label]
            labs = [i for i in labs if i not in labs0]
            matrix[(j,k)] = len(bd)
            bd.append(A.getObject(labs[0]))
            bd[-1].Label = bdname
            bd[-1].Placement = App.Placement(lTransfX*j + lTransfY*k+lTransfZ*n,App.Rotation(vector(0,0,1),0),vector(0,0,0))

            if cut:
                bdcutname = "bdcut_%d_%d"%(j,k)
                labs0 = [i.Label for i in A.findObjects() if "Body" in i.Label]
                A.addObject("Part::Feature","Body").Shape = shp_cut
                labs = [i.Label for i in A.findObjects() if "Body" in i.Label]
                labs = [i for i in labs if i not in labs0]
                bd_cut.append(A.getObject(labs[0]))
                bd_cut[-1].Label = bdcutname

                labs0 = [i.Label for i in A.findObjects() if "Common" in i.Label]
                A.addObject("Part::MultiCommon","Common")
                labs = [i.Label for i in A.findObjects() if "Common" in i.Label]
                labs = [i for i in labs if i not in labs0]
                commons.append(A.getObject(labs[0]))
                commons[-1].Shapes = [bd_cut[-1],bd[-1]]
                vertical_faces_lst.append(getFaceListMatch(commons[-1],vector(0.,0.,1.)))
                coords.append((j,k,n))
            E.recompute()
            f.writelines(["j = %d, k = %d, n = %d/%d \n"%(j,k,n+1,NbodiesZ)])

E.recompute()

f.writelines(["Removing cols outside the cylinder \n"])
ilist = []
for i in range(len(commons)):
    if (commons[i].Shape.Volume < 1e-16):
        #labs = (commons[i].Shapes[0].Label,commons[i].Shapes[1].Label)
        labs = (commons[i].Shapes[0].Name,commons[i].Shapes[1].Name)
        #print(commons[i].Label,labs)
        A.removeObject(commons[i].Label)
        A.recompute()
        A.removeObject(labs[0])
        A.recompute()
        A.removeObject(labs[1])
        A.recompute()
        ilist.append(i)

f.writelines(["Removing complete \n"])

tmp = commons + []
commons = [tmp[i] for i in range(len(commons)) if i not in ilist]

tmp = vertical_faces_lst + []
vertical_faces_lst = [tmp[i] for i in range(len(vertical_faces_lst)) if i not in ilist]

tmp = coords + []
coords = [tmp[i] for i in range(len(coords)) if i not in ilist]

Area_topBot_each = []
Area_topBot_only = 0.
for i,c in enumerate(commons):
    flist = vertical_faces_lst[i]
    Asum = np.sum([c.Shape.Faces[j].Area for j in flist])
    Area_topBot_each.append(Asum)
    ni = coords[i][2]
    if (ni == 0) | (ni == NbodiesZ-1):
        Area_topBot_only += Asum

Area_topBot_all = np.sum(Area_topBot_each)

vol_bdcut = body_cut.Shape.Volume
print("vol_bdcut = %.5g cu.mm"%(vol_bdcut)) 
print("lsweep = %.3g"%(lsweep))
A_packing = np.sum(np.array([i.Shape.Area for i in commons]))
A_packing_true = A_packing - Area_topBot_all + Area_topBot_only
vol_packing = np.sum(np.array([i.Shape.Volume for i in commons]))
print("A_packing = %.5g sq.mm"%(A_packing)) 
print("Area_topBot_all = %.5g sq.mm"%(Area_topBot_all)) 
print("A_packing_true = %.5g sq.mm"%(A_packing_true)) 
print("vol_packing = %.5g cu.mm"%(vol_packing)) 
print("A_packing_normalized = %.5g sq.mm/cu.mm"%(A_packing_true/vol_bdcut)) 
print("A_packing_normalized1 = %.5g sq.mm/cu.mm"%(A_packing_true/(vol_bdcut-vol_packing))) 
vol_unitCell = (bounding_box["xmax"] - bounding_box["xmin"])*(bounding_box["ymax"] - bounding_box["ymin"])*(bounding_box["zmax"]-bounding_box["zmin"])/2.
print("A_packing_unitCell = %.5g sq.mm/cu.mm"%(commons[0].Shape.Area/vol_unitCell))

f.writelines(["A_packing = %.5g sq.mm"%(A_packing) + "\n"]) 
f.writelines(["Area_topBot_all = %.5g sq.mm"%(Area_topBot_all) + "\n"]) 
f.writelines(["A_packing_true = %.5g sq.mm"%(A_packing_true) + "\n"]) 
f.writelines(["vol_packing = %.5g cu.mm"%(vol_packing) + "\n"]) 
f.writelines(["vol_bdcut = %.5g cu.mm"%(vol_bdcut) + "\n"]) 
f.writelines(["A_packing_normalized = %.5g sq.mm/cu.mm"%(A_packing_true/vol_bdcut) + "\n"]) 
f.writelines(["A_packing_normalized1 = %.5g sq.mm/cu.mm"%(A_packing_true/(vol_bdcut-vol_packing)) + "\n"]) 
E.recompute()


# Create column walls 
f.writelines(["Creating column outer shell \n"])
A.addObject('PartDesign::Body','Column')
column = A.Column

column.newObject("PartDesign::Plane","colPadPlane")
colPadPlane = A.getObject("colPadPlane")

if derive_Nz_from_packing_height:
    colPadPlane.AttachmentOffset = App.Placement(vector(0.,0.,1.)*(packing_height - hcol/2.) + vector(0.,0.,column_Z_offset_bottom),  App.Rotation(0.,0.,0.))
else:
    #colPadPlane.AttachmentOffset = App.Placement(vector(0.,0.,0.5)*lsweep*NbodiesZ - vector(0.,0.,column_Z_offset),  App.Rotation(0.,0.,0.))
    colPadPlane.AttachmentOffset = App.Placement(vector(0.,0.,1.)*(lsweep*NbodiesZ - hcol/2.) + vector(0.,0.,column_Z_offset_bottom),  App.Rotation(0.,0.,0.))
colPadPlane.Support = [A.getObject("XY_Plane002"),""]
colPadPlane.MapMode = "FlatFace"
colPadPlane.MapReversed = True
colPadPlane.Visibility = False
E.recompute()

column.newObject("Sketcher::SketchObject","column_profile")
s = A.getObject("column_profile")
s.Support = (colPadPlane,[""])
s.MapMode = "FlatFace"
s.Visibility = False
E.recompute()



# Line -- 0 (Circle)
s.addGeometry(Part.Circle(vector(0.,0.,0.),vector(0.,0.,1.),rcol),False)
s.addConstraint(constraint("Coincident",0,3,-1,1))
s.addConstraint(constraint("Radius",0,rcol))
s.setDatum(1,App.Units.Quantity(convertDistanceUnits(rcol)))
E.recompute()

# Line -- 1 (Circle)
s.addGeometry(Part.Circle(vector(0.,0.,0.),vector(0.,0.,1.),rcol+tcol),False)
s.addConstraint(constraint("Coincident",1,3,-1,1))
s.addConstraint(constraint("Radius",1,rcol+tcol))
s.setDatum(3,App.Units.Quantity(convertDistanceUnits(rcol+tcol)))
E.recompute()

column.newObject('PartDesign::Pad',"col_pad")
col_pad = A.getObject("col_pad")
col_pad.Profile = s
col_pad.Length = hcol
col_pad.Midplane = 1
E.recompute()

column.newObject("Sketcher::SketchObject","lid_sketch")
s = A.getObject("lid_sketch")
s.Support = (col_pad,["Face3"])
s.MapMode = "FlatFace"
s.Visibility = False

# Line -- 0 (Circle)
s.addGeometry(Part.Circle(vector(0.,0.,0.),vector(0.,0.,1.),rcol),False)
s.addConstraint(constraint("Coincident",0,3,-1,1))
s.addConstraint(constraint("Radius",0,rcol))
s.setDatum(1,App.Units.Quantity(convertDistanceUnits(rcol)))
E.recompute()

column.newObject("PartDesign::Pad","lid")
lid = A.getObject("lid")
lid.Profile = s
lid.Length = tcol
lid.Reversed = True

#dummyPlane_inlet.newObject("PartDesign::Plane","dummyPlane_inlet")
#dummyPlane_inlet = A.getObject("dummyPlane_inlet")
#dummyPlane_inlet.Support = [(lid,["Face4"]),""]
#dummyPlane_inlet.MapMode = "FlatFace"
#dummyPlane_inlet.MapReversed = True
#dummyPlane_inlet.Visibility = False
#E.recompute()



xr = rcol/(nInlet+1)
if (xr - 2.*rInlet) <= tcol:
    print("CANNOT PROCEED FURTHER. Reduce either the inlet radius of dripping points or the number of dripping points.")
    sys.exit()
k = 0
hname = "inletSketch_%d"%(k)
column.newObject("Sketcher::SketchObject",hname)
s = A.getObject(hname)
s.Support = (lid,["Face4"])
s.MapMode = "FlatFace"
s.Visibility = False
E.recompute()

for i in range(-nInlet-3,nInlet+4):
    for j in range(-nInlet-3,nInlet+4):
        #s.Visibility = False
        vec1 = np.array([xr*i,xr*j])
        val = 1. + rInlet/(np.sqrt(np.sum(np.power(vec1,2.))) + 1e-10)
        vec1 *= val
        val = np.sum(np.power(vec1,2.)) - (rcol - tcol)**2.
        if (val >= 0):
            continue
        s.addGeometry(Part.Circle(vector(xr*i,xr*j,0.),vector(0.,0.,1.),rInlet),False)
        s.addConstraint(constraint("DistanceX",k,3,xr*i))
        s.addConstraint(constraint("DistanceY",k,3,xr*j))
        s.addConstraint(constraint("Radius",k,rInlet))
        #print("hname = ",hname,xr*i,xr*j)
        s.setDatum(2+3*k,App.Units.Quantity(convertDistanceUnits(rInlet)))
        E.recompute()

        k += 1


hname = "gOutlet_pocket"
column.newObject("PartDesign::Pocket",hname)
gpkt = A.getObject(hname)
gpkt.Profile = s
gpkt.Length = tcol

k = 0
hname = "inletSketch_1"
column.newObject("Sketcher::SketchObject",hname)
s = A.getObject(hname)
s.Support = (lid,["Face6"])
s.MapMode = "FlatFace"
s.AttachmentOffset = App.Placement(App.Vector(0.0000000000, 0.0000000000, -tcol/2.),  App.Rotation(0.0000000000, 0.0000000000, 0.0000000000))
s.Visibility = False
E.recompute()

for i in range(-nInlet-3,nInlet+4):
    for j in range(-nInlet-3,nInlet+4):
        #s.Visibility = False
        vec1 = np.array([xr*i,xr*j])
        val = 1. + rInlet/(np.sqrt(np.sum(np.power(vec1,2.))) + 1e-10)
        vec1 *= val
        val = np.sum(np.power(vec1,2.)) - (rcol - tcol)**2.
        if (val >= 0):
            continue
        s.addGeometry(Part.Circle(vector(xr*i,xr*j,0.),vector(0.,0.,1.),rInlet),False)
        s.addConstraint(constraint("DistanceX",k,3,xr*i))
        s.addConstraint(constraint("DistanceY",k,3,xr*j))
        s.addConstraint(constraint("Radius",k,rInlet))
        #print("hname = ",hname,xr*i,xr*j)
        s.setDatum(2+3*k,App.Units.Quantity(convertDistanceUnits(rInlet)))
        E.recompute()

        k += 1

hname = "inlet"
#column.newObject("PartDesign::Pocket",hname)
column.newObject("PartDesign::Pad",hname)
hole = A.getObject(hname)
hole.Profile = s
hole.Length = tcol/2. + hInlet
hole.Reversed = False
#inlet = A.getObject(hname)

column.ViewObject.Transparency = 70
f.writelines(["Column outer shell complete \n"])


E.recompute()
E.recompute()

f.writelines(["Saving File \n"])
Gui.SendMsgToActiveView("ViewFit")
A.saveAs(export_dir + exportName + u".FCStd")
export_objects = [column] + commons
#ImportGui.export(export_objects,export_dir + exportName + u".step")
f.writelines(["Saving Complete \n"])

f.close()

commonLabsList = [i.Label for i in commons]
common_indices = np.array([0] + [int(i.replace("Common","")) for i in commonLabsList[1:]],dtype = np.int)
print(commonLabsList)
print(common_indices)

xlist = np.sort(np.unique([i[1] for i in coords]))
arrlist = "A,B,C,D,E,F,G,H,I,J,K,L,M,N".split(",")

j = 0
fusions = []
if apply_fusion:
    for xi,x in enumerate(xlist):
        clist = [commons[i] for i in range(len(commons)) if coords[i][1] == x]

        if len(clist) < 2:
            fusions.append(clist[0])
            continue

        fusionName = "fusion%d"%(j)
        A.addObject("Part::Fuse",fusionName)
        fusionObject = A.getObject(fusionName)
        fusionObject.Base = clist[0]
        fusionObject.Tool = clist[1]
        fusions = fusions + [fusionObject]
        fusionObject = None
        j += 1
        E.recompute()

        if len(clist) > 2:
            for kk in range(2,len(clist)):
                fusionName = "fusion%d"%(j)
                A.addObject("Part::Fuse",fusionName)
                fusionObject = A.getObject(fusionName)
                fusionObject.Base = clist[kk]
                fusionObject.Tool = fusions[-1]
                fusions = fusions + [fusionObject]
                fusionObject = None
                j += 1
                E.recompute()

        fusions[-1].Label = "Fusion" + arrlist[xi]

sys.exit()


A.addObject("Part::Compound","Compound")
A.Compound.Links = commons
E.recompute()


Gui.SendMsgToActiveView("ViewFit")

export_objects = [A.getObject("Compound"),column]

#IGES export
#Part.export(export_objects,export_dir + dName + u".iges")

# STEP export
ImportGui.export(export_objects,export_dir + exportName + u".step")

# STL export
#Mesh.export(export_objects,export_dir + dName + u".stl")
