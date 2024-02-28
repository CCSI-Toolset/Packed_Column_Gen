# -*- coding: utf-8 -*-

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
    a1 = a*180./np.pi
    return "%.6f deg"%(a1)

def setDatum_distance(base,l):
    base.setDatum(base.ConstraintCount-1,App.Units.Quantity(convertDistanceUnits(l)))

def setDatum_angle(base,l):
    base.setDatum(base.ConstraintCount-1,App.Units.Quantity(convertAngleUnits(l)))

def addline_pack(base,nprev,theta_val,Lpacking = None,attach_point = 2):
    vvec = vector(0.,1.,0.)
    l1 = base.getPoint(nprev,1) - base.getPoint(nprev,2)
    if attach_point == 1:
        l1 = -1.*l1
    if Lpacking is None: Lpacking = l1.Length
    l1 /= l1.Length
    z1 = (l1.x + 1j*(l1.y))*np.exp(1j*theta_val)*Lpacking
    base.addGeometry(line(base.getPoint(nprev,attach_point),base.getPoint(nprev,attach_point) - vvec),False)
    ncurr = base.GeometryCount - 1
    base.addConstraint(constraint("Coincident",ncurr,1,nprev,attach_point))

    if z1.real > 0:
        base.addConstraint(constraint("DistanceX",nprev,attach_point,ncurr,2,np.abs(z1.real)))
    if z1.real < 0:
        base.addConstraint(constraint("DistanceX",ncurr,2,nprev,attach_point,np.abs(z1.real)))
    setDatum_distance(base,np.abs(z1.real))

    if z1.imag > 0:
        base.addConstraint(constraint("DistanceY",nprev,attach_point,ncurr,2,np.abs(z1.imag)))
    if z1.imag < 0:
        base.addConstraint(constraint("DistanceY",ncurr,2,nprev,attach_point,np.abs(z1.imag)))
    setDatum_distance(base,np.abs(z1.imag))
    return ncurr

def offsetline_pack(base,noffset,nNext,nattach):
    l1 = base.getPoint(noffset,2) - base.getPoint(noffset,1)
    base.addGeometry(line(base.getPoint(nattach,2),base.getPoint(nattach,2) + l1),False)
    ncurr = base.GeometryCount - 1
    base.addConstraint(constraint("Coincident",ncurr,1,nattach,2))
    base.addConstraint(constraint("Parallel",ncurr,noffset))
    base.addGeometry(line(base.getPoint(ncurr,2),base.getPoint(ncurr,2) + vector(0,-1,0)),False)
    base.addConstraint(constraint("Coincident",ncurr+1,1,ncurr,2))
    if nNext > 0:
        base.addConstraint(constraint("PointOnObject",ncurr+1,2,nNext))
        base.addConstraint(constraint("Perpendicular",ncurr+1,nNext))
        base.addConstraint(constraint("Distance",ncurr+1,tcol))
        setDatum_distance(base,tcol)
        base.toggleConstruction(ncurr+1)
    else:
        base.addConstraint(constraint("Coincident",ncurr+1,2,noffset,2))
        base.addConstraint(constraint("Horizontal",ncurr+1))
    return ncurr


def get_NpackY(theta1,column_height,dpack,Ldrip,Lpack):
    Ny = np.int((column_height - Ldrip - 2.*dpack)/(Lpack*np.sin(theta1*np.pi/180.)))
    return Ny

def get_NpackX(theta1,tcol,pack_gap,column_radius,Lpack):
    Lh = 2.*column_radius
    Lunit = 2.*Lpack*np.cos(theta1*np.pi/180.) + 2.*tcol/np.sin(theta1*np.pi/180.) + 2.*pack_gap
    Nx = np.int((Lh - 1.*pack_gap )/Lunit)
    return Nx

def get_dpack_wall(theta1,tcol,pack_gap,column_radius,Lpack):
    Lh = 2.*column_radius
    Nx = get_NpackX(theta1,tcol,pack_gap,column_radius,Lpack)
    Lunit = 2.*Lpack*np.cos(theta1*np.pi/180.) + 2.*tcol/np.sin(theta1*np.pi/180.) + pack_gap
    dpack_wall = (Lh - Nx*Lunit)/(Nx+1)
    return dpack_wall

def Lrange(theta1,tcol,pack_gap,column_radius,dpack_wall = 0., Nx = 2):
    Lh = 2.*column_radius
    Lunit = 2.*tcol/np.sin(theta1*np.pi/180.) + pack_gap + dpack_wall
    Lmax = (Lh - Nx*Lunit - dpack_wall)/(2.*Nx*np.cos(theta1*np.pi/180.))
    return Lmax

# ========= User Parameters ===============
#export_dir  = u"C:/Users/shahya/Desktop/FreeCAD/2D/"
export_dir  = u"D:/OneDrive/NETL/FreeCAD_demo/2D/step/"
dName = "doc1"
column_radius = 38.1
column_height = 100
Ldrip = 3
Wdrip = 4
Ndrip = 9
Wgas = (2.*column_radius - Ndrip*Wdrip)/(Ndrip + 1)
tcol = 0.89
Lpack = 12 #13.5 #13.5*np.sqrt(2.)
dpack = 17
dpack_wall = Wgas*1.
pack_gap = 2.*tcol
#NpackY = 7
#theta1 = 30.
#theta2 = 2.*theta1
# =========================================

theta_mat = [30., 45., 60.]
L_mat = [10., 13., 14.8]
pack_mat  = [2.*tcol,3.*tcol,4.*tcol]

permutations = []
for k in range(len(pack_mat)):
    for j in range(len(L_mat)):
        for i in range(len(theta_mat)):
            permutations.append((i,j,k))

fname = export_dir + "log"
f = open(fname,"w")
f.writelines(["# N, theta1, Lpack, pack_gap, NpackY, dpack_wall\n"])
f.close()
for i,ii in enumerate(permutations):
    dName       = "design%d"%(i)
    theta1 = theta_mat[ii[0]] #30.
    theta2 = 2.*theta1
    Lpack = L_mat[ii[1]]
    pack_gap = pack_mat[ii[2]]
    NpackY = get_NpackY(theta1,column_height,dpack,Ldrip,Lpack)
    NpackX = get_NpackX(theta1,tcol,pack_gap,column_radius,Lpack)
    dpack_wall = get_dpack_wall(theta1,tcol,pack_gap,column_radius,Lpack)

    print("i = %d, Lmax = %.3g, dpack_wall = %.3g, NpackY = %d, NpackX = %d"%(i,Lrange(theta1,tcol,pack_gap,column_radius,dpack_wall = pack_gap),dpack_wall,NpackY, NpackX))

    f = open(fname,"a")
    f.writelines(["%d, %.3g, %.3g, %.3g, %d, %.3g \n"%(i, theta1, Lpack, pack_gap, NpackY, dpack_wall)])
    f.close()

    #ilist = [None,9]
    #if i not in ilist:
        #continue

    theta1_deg = theta1*1.
    theta2_deg = theta2*1.

    theta1 *= np.pi/180.
    theta2 *= np.pi/180.

    if Wgas < 0.:
        print("CANNOT PROCEED FURTHER. REDUCE THE NUMBER OF DRIPPING POINTS")
        print("EXITING... \n\n\n")
        sys.exit()

    if theta2 < theta1:
        print("CANNOT PROCEED FURTHER. THETA2 SHOULD AT LEAST BE %.3g degrees"%(theta1_deg))
        print("EXITING... \n\n\n")
        sys.exit()

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


    # Line -- 0
    base.addGeometry(line(vector(0.,0.,0.),vector(-1.*column_radius,0.,0)),False)
    base.addConstraint(constraint("Coincident",0,1,-1,1))
    base.addConstraint(constraint("Distance",0,column_radius))
    setDatum_distance(base,column_radius)
    base.addConstraint(constraint("Horizontal",0))

    # Line -- 1
    base.addGeometry(line(vector(0.,0.,0.),vector(column_radius,0.,0)),False)
    base.addConstraint(constraint("Coincident",1,1,-1,1))
    base.addConstraint(constraint("Horizontal",1))
    base.addConstraint(constraint("Equal",1,0))

    # Line -- 2
    base.addGeometry(line(base.getPoint(0,2),base.getPoint(0,2) + vector(0.,1.,0.)*column_height),False)
    base.addConstraint(constraint("Coincident",2,1,0,2))
    base.addConstraint(constraint("Vertical",2))
    base.addConstraint(constraint("Distance",2,column_height))
    setDatum_distance(base,column_height)

    # Line -- 3
    base.addGeometry(line(base.getPoint(1,2),base.getPoint(1,2) + vector(0.,1.,0.)*column_height),False)
    base.addConstraint(constraint("Coincident",3,1,1,2))
    base.addConstraint(constraint("Vertical",3))
    base.addConstraint(constraint("Equal",3,2))

    # Line -- 4
    base.addGeometry(line(base.getPoint(2,2),base.getPoint(2,2) + vector(1.,0.,0.)*Wgas),False)
    base.addConstraint(constraint("Coincident",4,1,2,2))
    base.addConstraint(constraint("Horizontal",4))
    base.addConstraint(constraint("Distance",4,Wgas))
    setDatum_distance(base,Wgas)

    # Line -- 5
    base.addGeometry(line(base.getPoint(4,2),base.getPoint(4,2) + vector(0.,-1.,0.)*Ldrip),False)
    base.addConstraint(constraint("Coincident",5,1,4,2))
    base.addConstraint(constraint("Vertical",5))
    base.addConstraint(constraint("Distance",5,Ldrip))
    setDatum_distance(base,Ldrip)

    # Line -- 6
    base.addGeometry(line(base.getPoint(5,2),base.getPoint(5,2) + vector(1.,0.,0.)*Wdrip),False)
    base.addConstraint(constraint("Coincident",6,1,5,2))
    base.addConstraint(constraint("Horizontal",6))
    base.addConstraint(constraint("Distance",6,Wdrip))
    setDatum_distance(base,Wdrip)

    num1 = 4 # Wgas
    num2 = 5 # Ldrip
    num3 = 6 # Wdrip
    hvec = vector(1.,0.,0.)
    vvec = vector(0.,1.,0.)

    # Line -- 7
    base.addGeometry(line(base.getPoint(6,2),base.getPoint(6,2) + vvec*Ldrip),False)
    base.addConstraint(constraint("Coincident",7,1,6,2))
    base.addConstraint(constraint("Vertical",7))
    base.addConstraint(constraint("Equal",7,num2))

    # Line -- 8
    base.addGeometry(line(base.getPoint(7,2),base.getPoint(7,2) + hvec*Wgas),False)
    base.addConstraint(constraint("Coincident",8,1,7,2))
    base.addConstraint(constraint("Horizontal",8))
    base.addConstraint(constraint("Equal",8,num1))

    n = base.GeometryCount - 1
    for i in range(1,Ndrip):
        n += 1
        base.addGeometry(line(base.getPoint(n-1,2),base.getPoint(n-1,2) - vvec*Ldrip),False)
        base.addConstraint(constraint("Coincident",n,1,n-1,2))
        base.addConstraint(constraint("Vertical",n))
        base.addConstraint(constraint("Equal",n,num2))
        
        n += 1
        base.addGeometry(line(base.getPoint(n-1,2),base.getPoint(n-1,2) + hvec*Wdrip),False)
        base.addConstraint(constraint("Coincident",n,1,n-1,2))
        base.addConstraint(constraint("Horizontal",n))
        base.addConstraint(constraint("Equal",n,num3))

        n += 1
        base.addGeometry(line(base.getPoint(n-1,2),base.getPoint(n-1,2) + vvec*Ldrip),False)
        base.addConstraint(constraint("Coincident",n,1,n-1,2))
        base.addConstraint(constraint("Vertical",n))
        base.addConstraint(constraint("Equal",n,num2))

        n += 1
        base.addGeometry(line(base.getPoint(n-1,2),base.getPoint(n-1,2) + hvec*Wgas),False)
        base.addConstraint(constraint("Coincident",n,1,n-1,2))
        if i != (Ndrip - 1):
            base.addConstraint(constraint("Horizontal",n))
            base.addConstraint(constraint("Equal",n,num1))
        else:
            base.addConstraint(constraint("Coincident",n,2,3,2))
            
    list_start  = []
    list_pack   = []
    list_offset = []

    # Line -- n+1
    Nstart = n + 1
    pt1 = base.getPoint(6,1) - vvec*dpack + hvec*(Lpack*np.cos(theta1))

    base.addGeometry(line(pt1,pt1 + hvec),False)
    base.addConstraint(constraint("Horizontal",n+1))
    base.addConstraint(constraint("DistanceY",n+1,1,6,1,dpack))
    setDatum_distance(base,dpack)
    s1 = Lpack*np.cos(theta1) - Wgas + dpack_wall
    base.addConstraint(constraint("DistanceX",6,1,n+1,1,s1))
    setDatum_distance(base,s1)
    s1 += tcol/np.sin(theta1)
    base.addConstraint(constraint("DistanceX",6,1,n+1,2,s1))
    setDatum_distance(base,s1)

    # Line -- n+2 
    Npack_list = [n+2]
    base.addGeometry(line(base.getPoint(n+1,1),base.getPoint(n+1,1) - vvec),False)
    agl = np.pi - theta1
    base.addConstraint(constraint("Coincident",n+2,1,n+1,1))
    base.addConstraint(constraint("Angle",n+2,1,n+1,1,agl))
    setDatum_angle(base,agl)
    base.addConstraint(constraint("Distance",n+2,Lpack))
    setDatum_distance(base,Lpack)

    agl = theta2*1.
    for i in range(1,NpackY):
        agl *= -1.
        Npack_list.append(addline_pack(base,Npack_list[-1],agl))

    Noffset_list = [offsetline_pack(base,Npack_list[0], Npack_list[1],Nstart)]
    for i in range(1,len(Npack_list)):
        if i == len(Npack_list)-1:
            Noffset_list.append(offsetline_pack(base,Npack_list[i],-1,Noffset_list[-1]))
        else:
            Noffset_list.append(offsetline_pack(base,Npack_list[i], Npack_list[i+1],Noffset_list[-1]))

    list_start.append(Nstart)
    list_pack.append(Npack_list)
    list_offset.append(Noffset_list)

    Nstart1 = base.GeometryCount

    l1 = base.getPoint(list_start[-1],2) - base.getPoint(list_start[-1],1)
    pt1 = base.getPoint(list_start[-1],2) + hvec*pack_gap
    base.addGeometry(line(pt1,pt1 + hvec),False)
    base.addConstraint(constraint("Horizontal",Nstart1))
    base.addConstraint(constraint("DistanceX",list_start[-1],2,Nstart1,1,pack_gap))
    setDatum_distance(base,pack_gap)
    base.addConstraint(constraint("DistanceX",Nstart1,1,Nstart1,2,l1.Length))
    setDatum_distance(base,l1.Length)
    base.addConstraint(constraint("DistanceY",Nstart1,1,6,1,dpack))
    setDatum_distance(base,dpack)

    Npack_list1 = [addline_pack(base,Nstart1,-theta1,Lpacking = Lpack,attach_point = 1)]
    agl = -1.*theta2
    for i in range(1,NpackY):
        agl *= -1.
        Npack_list1.append(addline_pack(base,Npack_list1[-1],agl))


    Noffset_list1 = [offsetline_pack(base,Npack_list1[0], Npack_list1[1],Nstart1)]
    for i in range(1,len(Npack_list1)):
        if i == len(Npack_list1)-1:
            Noffset_list1.append(offsetline_pack(base,Npack_list1[i],-1,Noffset_list1[-1]))
        else:
            Noffset_list1.append(offsetline_pack(base,Npack_list1[i], Npack_list1[i+1],Noffset_list1[-1]))

    list_start.append(Nstart1)
    list_pack.append(Npack_list1)
    list_offset.append(Noffset_list1)
    #x_max = -1e10
    #for i in range(len(Npack_list1)):
    #    x1 = np.max([base.getPoint(Npack_list1[i],1).x,base.getPoint(Npack_list1[i],2).x])
    #    x1 = np.max([x1,base.getPoint(Noffset_list1[i],1).x,base.getPoint(Noffset_list1[i],2).x])
    #    if x1 > x_max:
    #        x_max = x1*1.

    #if (x_max > 0.):
    #    print("CANNOT PROCEED FURTHER. INFEASIBLE GEOMETRY")
    #    print("EXITING... \n\n\n")
    #    sys.exit()

    for jj in range(1,NpackX):
        print("jj = ",jj)

        Nstart2 = base.GeometryCount

        l1 = Lpack*np.cos(theta1) + dpack_wall + Lpack*np.cos(theta1)
        pt1 = base.getPoint(list_start[-1],2) + hvec*l1*2.
        base.addGeometry(line(pt1,pt1+hvec),False)
        base.addConstraint(constraint("Horizontal",Nstart2))
        base.addConstraint(constraint("DistanceX",list_start[-1],2,Nstart2,1,l1))
        setDatum_distance(base,l1)
        base.addConstraint(constraint("DistanceY",Nstart2,1,6,1,dpack))
        setDatum_distance(base,dpack)
        l1 = base.getPoint(list_start[-1],2) - base.getPoint(list_start[-1],1)
        base.addConstraint(constraint("DistanceX",Nstart2,1,Nstart2,2,l1.Length))
        setDatum_distance(base,l1.Length)

        Npack_list2 = [addline_pack(base,Nstart2,np.pi + theta1,Lpacking = Lpack,attach_point = 1)]
        agl = 1.*theta2
        for i in range(1,NpackY):
            agl *= -1.
            Npack_list2.append(addline_pack(base,Npack_list2[-1],agl))


        Noffset_list2 = [offsetline_pack(base,Npack_list2[0], Npack_list2[1],Nstart2)]
        for i in range(1,len(Npack_list2)):
            if i == len(Npack_list2)-1:
                Noffset_list2.append(offsetline_pack(base,Npack_list2[i],-1,Noffset_list2[-1]))
            else:
                Noffset_list2.append(offsetline_pack(base,Npack_list2[i], Npack_list2[i+1],Noffset_list2[-1]))

        list_start.append(Nstart2)
        list_pack.append(Npack_list2)
        list_offset.append(Noffset_list2)

        Nstart2 = base.GeometryCount

        l1 = base.getPoint(list_start[-1],2) - base.getPoint(list_start[-1],1)
        pt1 = base.getPoint(list_start[-1],2) + hvec*pack_gap
        base.addGeometry(line(pt1,pt1 + hvec),False)
        base.addConstraint(constraint("Horizontal",Nstart2))
        base.addConstraint(constraint("DistanceX",list_start[-1],2,Nstart2,1,pack_gap))
        setDatum_distance(base,pack_gap)
        base.addConstraint(constraint("DistanceX",Nstart2,1,Nstart2,2,l1.Length))
        setDatum_distance(base,l1.Length)
        base.addConstraint(constraint("DistanceY",Nstart2,1,6,1,dpack))
        setDatum_distance(base,dpack)

        Npack_list2 = [addline_pack(base,Nstart2,-theta1,Lpacking = Lpack,attach_point = 1)]
        agl = -1.*theta2
        for i in range(1,NpackY):
            agl *= -1.
            Npack_list2.append(addline_pack(base,Npack_list2[-1],agl))


        Noffset_list2 = [offsetline_pack(base,Npack_list2[0], Npack_list2[1],Nstart2)]
        for i in range(1,len(Npack_list2)):
            if i == len(Npack_list2)-1:
                Noffset_list2.append(offsetline_pack(base,Npack_list2[i],-1,Noffset_list2[-1]))
            else:
                Noffset_list2.append(offsetline_pack(base,Npack_list2[i], Npack_list2[i+1],Noffset_list2[-1]))

        list_start.append(Nstart2)
        list_pack.append(Npack_list2)
        list_offset.append(Noffset_list2)


    body.newObject('PartDesign::Pad',"Pad")
    pad = A.getObject("Pad")
    pad.Profile = base
    pad.Length = 1.
    pad.Reversed = True
    E.recompute()

    base.Visibility = False

    export_objects = [pad]
    # STEP export
    ImportGui.export(export_objects,export_dir + dName + u".step")

    #sys.exit()
    App.closeDocument(dName)

