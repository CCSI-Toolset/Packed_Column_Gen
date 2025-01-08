# ***************************************************************************
# *                                                                         *
# *   Copyright: See License.md file                                        *
# *                                                                         *
# *   Authors: Yash Girish Shah, Grigorios Panagakos                        *
# *                                                                         *
# ***************************************************************************

from .preamble import *
from .commonFunctions import *

def addline_pack(base,nprev,theta_val,Lpacking = None,attach_point = 2):
    vvec = App.Vector(0.,1.,0.)
    l1 = base.getPoint(nprev,1) - base.getPoint(nprev,2)
    if attach_point == 1:
        l1 = -1.*l1
    if Lpacking is None: Lpacking = l1.Length
    l1 /= l1.Length
    z1 = (l1.x + 1j*(l1.y))*np.exp(1j*theta_val)*Lpacking
    base.addGeometry(Part.LineSegment(base.getPoint(nprev,attach_point),base.getPoint(nprev,attach_point) - vvec),False)
    ncurr = base.GeometryCount - 1
    base.addConstraint(Sketcher.Constraint("Coincident",ncurr,1,nprev,attach_point))

    if z1.real > 0:
        base.addConstraint(Sketcher.Constraint("DistanceX",nprev,attach_point,ncurr,2,np.abs(z1.real)))
    if z1.real < 0:
        base.addConstraint(Sketcher.Constraint("DistanceX",ncurr,2,nprev,attach_point,np.abs(z1.real)))
    setDatum_distance(base,np.abs(z1.real))

    if z1.imag > 0:
        base.addConstraint(Sketcher.Constraint("DistanceY",nprev,attach_point,ncurr,2,np.abs(z1.imag)))
    if z1.imag < 0:
        base.addConstraint(Sketcher.Constraint("DistanceY",ncurr,2,nprev,attach_point,np.abs(z1.imag)))
    setDatum_distance(base,np.abs(z1.imag))
    return ncurr

def offsetline_pack(base,noffset,nNext,nattach,tcol):
    l1 = base.getPoint(noffset,2) - base.getPoint(noffset,1)
    base.addGeometry(Part.LineSegment(base.getPoint(nattach,2),base.getPoint(nattach,2) + l1),False)
    ncurr = base.GeometryCount - 1
    base.addConstraint(Sketcher.Constraint("Coincident",ncurr,1,nattach,2))
    base.addConstraint(Sketcher.Constraint("Parallel",ncurr,noffset))
    base.addGeometry(Part.LineSegment(base.getPoint(ncurr,2),base.getPoint(ncurr,2) + App.Vector(0,-1,0)),False)
    base.addConstraint(Sketcher.Constraint("Coincident",ncurr+1,1,ncurr,2))
    if nNext > 0:
        base.addConstraint(Sketcher.Constraint("PointOnObject",ncurr+1,2,nNext))
        base.addConstraint(Sketcher.Constraint("Perpendicular",ncurr+1,nNext))
        base.addConstraint(Sketcher.Constraint("Distance",ncurr+1,tcol))
        setDatum_distance(base,tcol)
        base.toggleConstruction(ncurr+1)
    else:
        base.addConstraint(Sketcher.Constraint("Coincident",ncurr+1,2,noffset,2))
        base.addConstraint(Sketcher.Constraint("Horizontal",ncurr+1))
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
