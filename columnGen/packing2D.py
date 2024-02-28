# ***************************************************************************
# *                                                                         *
# *   PackedColumnGen Copyright (c) 2012 - 2023, by the software owners:    *
# *   Oak Ridge Institute for Science and Education (ORISE), TRIAD          *
# *   National Security, LLC., Lawrence Livermore National Security, LLC.,  *
# *   The Regents of the University of California, through Lawrence         *
# *   Berkeley National Laboratory, Battelle Memorial Institute, Pacific    *
# *   Northwest Division through Pacific Northwest National Laboratory,     *
# *   Carnegie Mellon University, West Virginia University, Boston          *
# *   University, the Trustees of Princeton University, The University of   *
# *   Texas at Austin, URS Energy & Construction, Inc., et al.              *
# *   All rights reserved.                                                  *
# *                                                                         *
# *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS   * 
# *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     * 
# *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS     * 
# *   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE        * 
# *   COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT    * 
# *   , INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            * 
# *   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR    * 
# *   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)    * 
# *   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN             * 
# *   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR          * 
# *   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,        * 
# *   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    * 
# *                                                                         *
# *   Authors: Yash Girish Shah, Grigorios Panagakos                        *
# *                                                                         *
# ***************************************************************************

from .preamble import *
from .commonFunctions import *
from .Functions2D import *

class PackedColumn2D:


    def __init__(self,obj\
            ,column_radius = 38.1
            ,column_height = 100
            ,Ldrip = 3 \
            ,Wdrip = 4 \
            ,Ndrip = 9 \
            ,tcol = 0.89 \
            ,Lpack = 12 \
            ,dpack = 17 \
            ,pack_gap = 1.78\
            ,theta1 = 45\
            ): 

        self.Object = obj
        obj.Proxy = self
        self.Type_ = "Column2D"

        self.propertiesDict = {}
        self.removedPropertiesDict = {}
        self.MainShp = None
        self.objects = []
        self.resetNeeded = False
        self.columnCreated = False

        alen = "App::PropertyLength"
        aang = "App::PropertyAngle"
        abool = "App::PropertyBool"
        aenum = "App::PropertyEnumeration"
        aint = "App::PropertyInteger"
        aflt = "App::PropertyFloat"

        subsect = "Column Dimensions"
        self._addProperty(alen,"ColumnWidth",subsect,"Width of the column.").ColumnWidth = 2.*column_radius
        self._addProperty(alen,"ColumnHeight",subsect,"Height of the column. ").ColumnHeight = column_height

        subsect = "Packing Parameters"
        self._addProperty(alen,"CellLength",subsect,"Length of the repetitive elementary cell").CellLength = Lpack
        self._addProperty(alen,"PackingPlacement",subsect,"Distance between the solvent inlet (dripping points) and the top face of the packing.").PackingPlacement = dpack
        self._addProperty(alen,"PackingSpacing",subsect,"Spacing between adjacent 2D packing sheets. (along X)").PackingSpacing = pack_gap
        self._addProperty(alen,"Thickness",subsect,"Thickness of the packing material").Thickness = tcol
        self._addProperty(aang,"PackingAngle",subsect,"Angle of the packing walls from the horizontal.").PackingAngle = theta1

        subsect = "Dripping Points"
        self._addProperty(alen,"DripWidth", subsect,"Width of the dripping points.").DripWidth = Wdrip
        self._addProperty(alen,"Depth", subsect,"Depth of the dripping points.").Depth = Ldrip
        self._addProperty(aint,"NumberOfDrippingPoints", subsect,"Number of dripping points along the column width.").NumberOfDrippingPoints = Ndrip

        subsect = "General Properties"
        self._addProperty(aenum,"GeometryConstructionType",subsect,"Specify whether the geometry to be constructed is 3D or 2D").GeometryConstructionType = ["3D","2D"]
        self._addProperty(alen,"ColumnThickness",subsect,"Thickness of the column for 3D construction.").ColumnThickness = tcol
        obj.GeometryConstructionType = "2D"
        if (obj.GeometryConstructionType == "2D"):
            self._removeProperty("ColumnThickness")

    def _addProperty(self,propertyType,propertyName,propertySection,propertyDescription):
        if propertyName not in self.propertiesDict:
            self.propertiesDict[propertyName] = [propertyType,propertyName,propertySection,propertyDescription]
            self.propertiesDict[propertyName].append(self.Object.addProperty(propertyType,propertyName,propertySection,propertyDescription))
        return self.propertiesDict[propertyName][-1]

    def _removeProperty(self,item):
        obj = self.Object
        if hasattr(obj,item):
            self.removedPropertiesDict[item] = obj.__getattribute__(item)
            obj.removeProperty(item)
            return True
        else:
            return False

    def _addRemovedProperty(self,item):
        obj = self.Object
        if (not hasattr(obj,item)) & (item in self.removedPropertiesDict):
            obj.addProperty(self.propertiesDict[item][0],self.propertiesDict[item][1],self.propertiesDict[item][2],self.propertiesDict[item][3])
            obj.__setattr__(item,self.removedPropertiesDict.pop(item))
            return True
        else:
            return False

    def propertiesList(self):
        return list(self.propertiesDict.keys())

    def onChanged(self,obj,prop):

        if prop == "GeometryConstructionType":
            if obj.GeometryConstructionType == "2D":
                self._removeProperty("ColumnThickness")

            if obj.GeometryConstructionType == "3D":
                self._addRemovedProperty("ColumnThickness")

        if hasattr(self,"propertiesDict"):
            if prop in self.propertiesList():
                self.resetNeeded = True

    def createColumn(self):
        obj = self.Object
        column_radius = float(obj.ColumnWidth)/2.
        column_height = float(obj.ColumnHeight)
        Ldrip = float(obj.Depth)
        Wdrip = float(obj.DripWidth)
        Ndrip = int(obj.NumberOfDrippingPoints)
        Wgas = (2.*column_radius - Ndrip*Wdrip)/(Ndrip + 1)
        tcol = float(obj.Thickness)
        Lpack = float(obj.CellLength)
        dpack = float(obj.PackingPlacement)
        pack_gap = float(obj.PackingSpacing)
        theta1 = float(obj.PackingAngle)
        output_type = obj.GeometryConstructionType
        if output_type == "3D":
            col_thickness = float(obj.ColumnThickness)



        theta2 = 2.*theta1
        NpackY = get_NpackY(theta1,column_height,dpack,Ldrip,Lpack)
        NpackX = get_NpackX(theta1,tcol,pack_gap,column_radius,Lpack)
        dpack_wall = get_dpack_wall(theta1,tcol,pack_gap,column_radius,Lpack)


        #print("i = %d, Lmax = %.3g, dpack_wall = %.3g, NpackY = %d, NpackX = %d"%(i,Lrange(theta1,tcol,pack_gap,column_radius,dpack_wall = pack_gap),dpack_wall,NpackY, NpackX))

        theta1_deg = theta1*1.
        theta2_deg = theta2*1.

        theta1 *= np.pi/180.
        theta2 *= np.pi/180.

        # Perform checks
        if Wgas <= 0.:
            QtGui.QMessageBox.information(None,"2D Column Building Error", "CANNOT PROCEED FURTHER as there are too many dripping points than the column width can accommodate. REDUCE THE NUMBER OF DRIPPING POINTS. Building Failed...")
            return

        if NpackY <= 1:
            QtGui.QMessageBox.information(None,"2D Column Building Error", "Cannot proceed further as the column is not tall enough to accommodate a single unit of packing. Reduce cell length or increase column height. Building Failed...")
            return

        if NpackX == 0:
            QtGui.QMessageBox.information(None,"2D Column Building Error", "Cannot proceed further as the column is not wide enough to accommodate a single unit of packing. Reduce cell length or increase column width. Building Failed...")
            return


        if theta2 < theta1:
            print("CANNOT PROCEED FURTHER. THETA2 SHOULD AT LEAST BE %.3g degrees"%(theta1_deg))
            print("EXITING... \n\n\n")
            sys.exit()



        Gui = FreeCADGui
        E = App.ActiveDocument
        F = Gui.ActiveDocument
        E.recompute()

        dName = E.Name

        G = Gui.getDocument(dName)
        A = App.getDocument(dName)

        body = A.addObject("PartDesign::Body","Body")

        base = body.newObject("Sketcher::SketchObject","Base")
        base.Support = (A.getObject("XY_Plane"),[""])
        base.MapMode = "FlatFace"
        E.recompute()

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
        if Ndrip > 1:
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
        # Old code
        #s1 = Lpack*np.cos(theta1) - Wgas + dpack_wall
        #base.addConstraint(constraint("DistanceX",6,1,n+1,1,s1))
        #setDatum_distance(base,s1)
        #s1 += tcol/np.sin(theta1)
        #base.addConstraint(constraint("DistanceX",6,1,n+1,2,s1))
        #setDatum_distance(base,s1)

        # New code
        s1 = dpack_wall + Lpack*np.cos(theta1)
        base.addConstraint(constraint("DistanceX",2,1,n+1,1,s1))
        setDatum_distance(base,s1)
        s1 += tcol/np.sin(theta1)
        base.addConstraint(constraint("DistanceX",2,1,n+1,2,s1))
        setDatum_distance(base,s1)


        # Line -- n+2 
        Npack_list = [n+2]
        base.addGeometry(line(base.getPoint(n+1,1),base.getPoint(n+1,1) - vvec),False)
        agl = np.pi - theta1
        agl_deg = 180. - theta1_deg
        base.addConstraint(constraint("Coincident",n+2,1,n+1,1))
        base.addConstraint(constraint("Angle",n+2,1,n+1,1,agl))
        setDatum_angle(base,agl_deg)
        base.addConstraint(constraint("Distance",n+2,Lpack))
        setDatum_distance(base,Lpack)

        agl = theta2*1.
        for i in range(1,NpackY):
            agl *= -1.
            Npack_list.append(addline_pack(base,Npack_list[-1],agl))

        Noffset_list = [offsetline_pack(base,Npack_list[0], Npack_list[1],Nstart,tcol)]
        for i in range(1,len(Npack_list)):
            if i == len(Npack_list)-1:
                Noffset_list.append(offsetline_pack(base,Npack_list[i],-1,Noffset_list[-1],tcol))
            else:
                Noffset_list.append(offsetline_pack(base,Npack_list[i], Npack_list[i+1],Noffset_list[-1],tcol))

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


        Noffset_list1 = [offsetline_pack(base,Npack_list1[0], Npack_list1[1],Nstart1,tcol)]
        for i in range(1,len(Npack_list1)):
            if i == len(Npack_list1)-1:
                Noffset_list1.append(offsetline_pack(base,Npack_list1[i],-1,Noffset_list1[-1],tcol))
            else:
                Noffset_list1.append(offsetline_pack(base,Npack_list1[i], Npack_list1[i+1],Noffset_list1[-1],tcol))

        list_start.append(Nstart1)
        list_pack.append(Npack_list1)
        list_offset.append(Noffset_list1)



        for jj in range(1,NpackX):

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


            Noffset_list2 = [offsetline_pack(base,Npack_list2[0], Npack_list2[1],Nstart2,tcol)]
            for i in range(1,len(Npack_list2)):
                if i == len(Npack_list2)-1:
                    Noffset_list2.append(offsetline_pack(base,Npack_list2[i],-1,Noffset_list2[-1],tcol))
                else:
                    Noffset_list2.append(offsetline_pack(base,Npack_list2[i], Npack_list2[i+1],Noffset_list2[-1],tcol))

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


            Noffset_list2 = [offsetline_pack(base,Npack_list2[0], Npack_list2[1],Nstart2,tcol)]
            for i in range(1,len(Npack_list2)):
                if i == len(Npack_list2)-1:
                    Noffset_list2.append(offsetline_pack(base,Npack_list2[i],-1,Noffset_list2[-1],tcol))
                else:
                    Noffset_list2.append(offsetline_pack(base,Npack_list2[i], Npack_list2[i+1],Noffset_list2[-1],tcol))

            list_start.append(Nstart2)
            list_pack.append(Npack_list2)
            list_offset.append(Noffset_list2)


        pad = body.newObject('PartDesign::Pad',"Pad")
        pad.Profile = base
        if output_type == "3D":
            pad.Length = col_thickness
        else:
            pad.Length = 1.
        pad.Reversed = True
        E.recompute()

        base.Visibility = False
        E.recompute()

        self.Children = [] 

        if output_type == "2D":
            flist = pad.Shape.Faces
            normalZ = vector(0.,0.,1.)
            flist_indices = [i for i,f in enumerate(flist) if np.abs(getFaceNormal(f).dot(normalZ)) == 1.]
            flist_indices = [i for i in flist_indices if flist[i].CenterOfMass[2] == 0.]
            if len(flist_indices) == 0:
                QtGui.QMessageBox.information(None, "Column Building Error", "Error occured in building column. Building Failed...")
                return
            i = flist_indices[0]

            faceShp = Part.getShape(body,"Pad.Face%d"%(i+1),needSubElement = True,refine = True)
            self.MainShp = faceShp
        else:
            self.MainShp = body.Shape

        obj.Shape = self.MainShp

        body.removeObjectsFromDocument()
        A.removeObject(body.Label)
        A.recompute()

        self.columnCreated = True
        self.resetNeeded = False

    def execute (self,obj):
        #if self.MainShp is not None:
            #obj.Shape = self.MainShp
        pass

    def reset(self):
        Gui = FreeCADGui
        E = App.ActiveDocument
        F = Gui.ActiveDocument
        dName = E.Name
        G = Gui.getDocument(dName)
        A = App.getDocument(dName)
        
        self.resetNeeded = False
        self.columnCreated = False


class ViewProviderColumn2D:
    
    obj_name = "Packed Column"
    
    def __init__(self, vobj, obj_name):
        self.obj_name = obj_name
        self.Object = vobj.Object # Original Object
        vobj.Proxy = self

    def attach(self, vobj):
        #vobj.addExtension("Gui::ViewProviderOriginGroupExtensionPython")
        self.ViewObject = vobj
        return

    def updateData(self, fp, prop):
        return

    def getDisplayModes(self,obj):
        return "As Is"
        
    def getDefaultDisplayMode(self):
        return "As Is"

    def setDisplayMode(self,mode):
        return "As Is"

    def onChanged(self, vobj, prop):
        pass
        
    def getIcon(self):
        return getWorkbenchFolder() + "/Resources/Icons/columnTree.svg" # + (self.obj_name).lower() + '.svg'"
        
    def __getstate__(self):
        return None

    def __setstate__(self,state):
        return None

    def onDelete(self,feature,subelements):
        try:
            self.Object.Proxy.reset()
        except:
            msg = "Test: This object is being deleted"
        #App.Console.PrintMessage(msg)
        return True


    def claimChildren(self):
        objs = []
        if hasattr(self.Object.Proxy,"Children"):
            objs.extend(self.Object.Proxy.Children)
        if hasattr(self.Object.Proxy,"Base"):
            objs.append(self.Object.Proxy.Base)
        if hasattr(self.Object,"Base"):
            objs.append(self.Object.Base)
        if hasattr(self.Object,"Objects"):
            objs.extend(self.Object.Objects)
        if hasattr(self.Object,"Components"):
            objs.extend(self.Object.Components)
        if hasattr(self.Object,"Children"):
            objs.extend(self.Object.Children)

        return objs
