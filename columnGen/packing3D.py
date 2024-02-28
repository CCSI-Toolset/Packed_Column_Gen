# ***************************************************************************
# *                                                                         *
# *   Copyright: See License.md file                                        *
# *                                                                         *
# *   Authors: Yash Girish Shah, Grigorios Panagakos                        *
# *                                                                         *
# ***************************************************************************

from .preamble import *
from .commonFunctions import *


class PackedColumn:
 
    def __init__(self, obj, H = 14.75, t1 = 0.890955, t2          = 0.18, h = 0.891155 \
        ,beta_angle  = 90. \
        ,gamma_angle = 45. \
        ,rcut        = 38.1 \
        ,hcol     = 100 \
        ,NbodiesX = 7 \
        ,NbodiesY = 7 \
        ,NbodiesZ = 4 \
        ,nInlet   = 4 \
        ,rInlet   = 2 \
        ,hInlet   = 15 - 8 \
        ,intensified_device = False \
        ,packing_height = 55.5 \
        ,tcol       = 0.890955*1./3. \
        ,column_Z_offset = 17.5 \
        ):

        self.Object = obj
        obj.Proxy = self
        self.Type_ = "Column3D"

        self.propertiesDict = {}
        self.removedPropertiesDict = {}
        self._addProperty("App::PropertyLength","CellLength","Packing Dimensions","Length of the repetitive elementary cell").CellLength = H
        self._addProperty("App::PropertyLength","Thickness1", "Packing Dimensions","Thickness of the packing material").Thickness1 = t1

        self._addProperty("App::PropertyEnumeration","HeightSpecificationMethod", "Packing Dimensions","Method for defining packing height. \n (1) Define by height. \n (2) Define by number of packing layers.").HeightSpecificationMethod = ["Packing Height", "Layers"]
        obj.HeightSpecificationMethod = "Packing Height"
        self._addProperty("App::PropertyLength","PackingHeight", "Packing Dimensions","Packing height. Applicable only if packing height is specified as the method of determining height.").PackingHeight = packing_height
        self._addProperty("App::PropertyInteger","PackingLayers", "Packing Dimensions","Number of layers of the elementary cells along packing height. Applicable only if layers is specified as the method of determining height.").PackingLayers = NbodiesZ
        self._removeProperty("PackingLayers")
        #self._PackingLayers = obj.PackingLayers; obj.removeProperty("PackingLayers")
        self._addProperty("App::PropertyLength","PackingPlacement", "Packing Dimensions","Distance between the solvent inlet (dripping points) and the top face of the packing.").PackingPlacement = column_Z_offset

        #self._addProperty("App::PropertyInteger","AlongX", "Number of packing elementary cells","Number of packing elementary cell units along the X coordinate").AlongX = NbodiesX
        #self._addProperty("App::PropertyInteger","AlongY", "Number of packing elementary cells","Number of packing elementary cell units along the Y coordinate").AlongY = NbodiesY

        self._addProperty("App::PropertyLength","DripRadius", "Dripping Points","Radius of the dripping points.").DripRadius = rInlet
        self._addProperty("App::PropertyLength","Depth", "Dripping Points","Depth of the dripping points.").Depth = hInlet
        self._addProperty("App::PropertyInteger","NumberAlongRadius", "Dripping Points","Number of dripping points along the column radius.").NumberAlongRadius = nInlet

        self._addProperty("App::PropertyAngle","CorrugationAngle", "Packing Angles","Corrugation angle (beta) of the packing. ").CorrugationAngle = beta_angle
        self._addProperty("App::PropertyAngle","SweepAngle", "Packing Angles","Sweep angle (gamma) of the packing. ").SweepAngle = gamma_angle

        self._addProperty("App::PropertyLength","ColumnRadius", "Column Dimensions","Inner radius of the column body. ").ColumnRadius = rcut
        self._addProperty("App::PropertyLength","ColumnHeight", "Column Dimensions","Height of the column. ").ColumnHeight = hcol
        self._addProperty("App::PropertyLength","ColumnWallThickness", "Column Dimensions","Thickness of the column wall. ").ColumnWallThickness = tcol

        self._addProperty("App::PropertyBool","EnableIntensifiedDevice", "Intensified Device Settings","Specify whether to construct intensified packing. ").EnableIntensifiedDevice = intensified_device
        self._addProperty("App::PropertyLength","CoolingChannelWidth", "Intensified Device Settings","Width of the cooling channels (only applicable if creating intensified packing).").CoolingChannelWidth = h
        self._addProperty("App::PropertyLength","Thickness2", "Intensified Device Settings","Thickness of the packing material at cell edges. (only applicable if creating intensified packing)").Thickness2 = t2
        if not intensified_device:
            self._removeProperty("CoolingChannelWidth")
            self._removeProperty("Thickness2")

        self.MainShp = None
        self.objects = []
        self.resetNeeded = False
        self.columnCreated = False
        #self.createColumn()

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

        if prop == "HeightSpecificationMethod":

            if obj.HeightSpecificationMethod == "Packing Height":
                self._removeProperty("PackingLayers")
                self._addRemovedProperty("PackingHeight")

            if obj.HeightSpecificationMethod == "Layers":
                self._removeProperty("PackingHeight")
                self._addRemovedProperty("PackingLayers")

        if prop == "EnableIntensifiedDevice":

            for item in "CoolingChannelWidth,Thickness2".split(","):
                if obj.EnableIntensifiedDevice:
                    self._addRemovedProperty(item)
                else:
                    self._removeProperty(item)


        if hasattr(self,"propertiesDict"):
            if prop in self.propertiesList():
                self.resetNeeded = True

    def createColumn(self):
        obj = self.Object
        H                  = float(obj.CellLength) #14.75 #M327Y equivalent
        beta_angle         = float(obj.CorrugationAngle) #67.06 #89.      # 90.
        gamma_angle        = float(obj.SweepAngle)      # 45.
        t1                 = float(obj.Thickness1) #0.890955*0.33 # 0.890955
        rcol               = float(obj.ColumnRadius) #38.1 # 1.5 inch
        rcut               = rcol #+ 6.*t1
        #NbodiesX           = int(obj.AlongX)
        #NbodiesY           = int(obj.AlongY)
        tcol               = float(obj.ColumnWallThickness) # Column thickness
        hcol               = float(obj.ColumnHeight) # Column height
        nInlet             = int(obj.NumberAlongRadius)
        rInlet             = float(obj.DripRadius)
        hInlet             = float(obj.Depth) # Inlet dripping points depth
        column_Z_offset    = float(obj.PackingPlacement)
        intensified_device = obj.EnableIntensifiedDevice

        if obj.HeightSpecificationMethod == "Packing Height":
            derive_Nz_from_packing_height = True
            packing_height  = float(obj.PackingHeight)
        else:
            derive_Nz_from_packing_height = False
            NbodiesZ        = int(obj.PackingLayers)

        if intensified_device:
            t2              = float(obj.Thickness2) #0.18*0.33     # 0.18
            h               = float(obj.CoolingChannelWidth) #0.891155*0.33      # 0.891155
        else:
            t1 *= 1./3.
            h  = t1*1.
            t2 = t1/5.


        cut = True
        alpha_angle = 90. - beta_angle/2. # <-- Parallelity
        axis_overlap = 4e-3 # 4e-5
        modify_hcol = False
        apply_fusion = False
        apply_compound = True


        t3  = t1/2.*np.tan(alpha_angle*np.pi/180.)  #0.45     # 0.72
              

        # Perform checks
        xr = rcol/(nInlet+1)
        if (xr - 2.*rInlet) <= tcol:
            QtGui.QMessageBox.information(None, "Column Building Error", "Too many dripping points. Building Failed. Reduce either the inlet radius or the number of dripping points.")
            return

        if t1 <= 0:
            QtGui.QMessageBox.information(None, "Column Building Error", "A non-zero thickness is required to build the column. Buliding Failed...")
            return

        if intensified_device:
            if (h <= 0) | (t2 <= 0):
                QtGui.QMessageBox.information(None, "Column Building Error", "A non-zero thickness and cooling channel width is required to build the column. Buliding Failed...")
                return

        if (beta_angle < 30.) | (beta_angle > 90.) \
           | (gamma_angle < 30.) | (gamma_angle > 60.):
               QtGui.QMessageBox.information(None,"Warning","The range of angles tested are between 30 and 90 degrees for the corrugation angle, and between 30 and 60 degrees for the sweep angle. Proceed with caution.")




        Gui = FreeCADGui
        E = App.ActiveDocument
        F = Gui.ActiveDocument
        E.recompute()

        #Gui.activateWorkbench("PartDesignWorkbench")

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
        hvec = vector(1.,0.,0.)
        vvec = vector(0.,1.,0.)

        lines = []

        # Line -- 0
        base.addGeometry(line(vector(0.,0.,0.),vector(-1.,1.,0)),False)
        lines.append(getLineNum(base))
        base.addConstraint(Sketcher.Constraint("Coincident",lines[0],1,-1,1))

        s1 = np.abs((2*t1 + h)/np.cos(deg2rad(alpha_angle)))
        #print(convertDistanceUnits(s1))
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
        #print(convertDistanceUnits(l1))
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
        #print("l1 = ",convertDistanceUnits(l1))
        L = l1/np.abs(np.sin(deg2rad(gamma_angle)))

        if not derive_Nz_from_packing_height:
            packing_height = lsweep*NbodiesZ

        column_Z_offset_bottom = hcol - column_Z_offset - packing_height - hInlet - tcol


        # Datum line -- Diagonal

        diagonal_datum = body.newObject("PartDesign::Line","Diagonal")
        diagonal_datum.Support = [(base,"Vertex1"),(base,"Vertex5")]
        diagonal_datum.MapMode = 'TwoPointLine'
        diagonal_datum.Visibility = False
        E.recompute()

        diagPlane = body.newObject("PartDesign::Plane","DiagPlane")
        diagPlane.Support = [(diagonal_datum,"")]
        diagPlane.MapMode = 'ObjectYZ'
        diagPlane.Visibility = False
        E.recompute()


        sweepProfile = body.newObject("Sketcher::SketchObject","sweepProfile")
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
        halfCell = body.newObject("PartDesign::AdditivePipe","HalfCell")
        halfCell.Profile = base
        halfCell.Spine = sweepProfile
        base.Visibility = False
        sweepProfile.Visibility = False
        E.recompute()



        oSk = body.newObject("Sketcher::SketchObject","oSk")
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


        YZof = body.newObject("PartDesign::Plane","YZof")
        YZof.Support = [(oSk,"Vertex2"),(A.getObject("Y_Axis"),""),(A.getObject("Z_Axis"),"")]
        YZof.MapMode = 'OXY'
        E.recompute()



        dum = body.newObject('Sketcher::SketchObject','dum')
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

        sbd = body.newObject('PartDesign::ShapeBinder','Rotation_axis')
        sbd.Support = [(dum,"Edge1")]
        E.recompute()

        dum.Visibility = False

        pol_pat = body.newObject('PartDesign::PolarPattern','PolarPattern')
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
        self.vlist = vlist
        xlist = np.array([i[0] for i in vlist])
        ylist = np.array([i[1] for i in vlist])
        zlist = np.array([i[2] for i in vlist])
        bounding_box["xmin"] = xlist.min()
        bounding_box["ymin"] = ylist.min()
        bounding_box["zmin"] = zlist.min()
        bounding_box["xmax"] = xlist.max()
        bounding_box["ymax"] = ylist.max()
        bounding_box["zmax"] = zlist.max()
        #print(bounding_box)
        topBB = {}
        vtop = [i for i in vlist if np.abs(i[2] - bounding_box["zmax"]) < 1e-5]
        topBB["xmin"] = np.min([i[0] for i in vtop])
        topBB["xmax"] = np.max([i[0] for i in vtop])
        topBB["ymin"] = np.min([i[1] for i in vtop])
        topBB["ymax"] = np.max([i[1] for i in vtop])
        #print(topBB)

        botBB = {}
        vbot = [i for i in vlist if np.abs(i[2] - bounding_box["zmin"]) < 1e-5]
        botBB["xmin"] = np.min([i[0] for i in vbot])
        botBB["xmax"] = np.max([i[0] for i in vbot])
        botBB["ymin"] = np.min([i[1] for i in vbot])
        botBB["ymax"] = np.max([i[1] for i in vbot])
        #print(botBB)

        org_bounding_box = {}
        bb = pol_pat.Shape.BoundBox
        org_bounding_box["xmin"] = bb.XMin
        org_bounding_box["ymin"] = bb.YMin
        org_bounding_box["zmin"] = bb.ZMin
        org_bounding_box["xmax"] = bb.XMax
        org_bounding_box["ymax"] = bb.YMax
        org_bounding_box["zmax"] = bb.ZMax


        # Lid covers

        try:
            shp = Part.getShape(body,"",needSubElement=True,refine=True)
        except:
            shp = Part.getShape(body,"",needSubElement=True,refine=False)

        lTransfX = (base.getPoint(lines[2],1) - base.getPoint(lines[11],2))
        lTransfX *= (lTransfX.Length)/lTransfX.Length
        lTransfY = lTransfX.cross(vector(0,0,-1))
        lTransfY /= lTransfY.Length

        ll = getLtransfYScale(pol_pat,lTransfY)
        lTransfY *= ll

        self.lx = lTransfX
        self.ly = lTransfY

        NbodiesX,NbodiesY = getNbodies(topBB,botBB,rcut,lTransfX,lTransfY)
        self.NbodiesX = NbodiesX
        self.NbodiesY = NbodiesY
        #print("Computed Nbodies:", tstX,tstY)

        # Cut Packing in a cylinder
        body_cut = A.addObject('PartDesign::Body','BodyCut')

        base_cut = body_cut.newObject("Sketcher::SketchObject","BaseCut")
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

        cut_pad = body_cut.newObject('PartDesign::Pad',"cut_pad")
        cut_pad.Profile = base_cut
        cut_pad.Length = packing_height
        cut_pad.Reversed = True
        E.recompute()

        shp_cut = Part.getShape(body_cut,"",needSubElement=False,refine=False)

        commons = [A.addObject("Part::MultiCommon","Common")]
        commons[-1].Shapes = [body_cut,body]
        commons[-1].Refine = True
        E.recompute()

        coords = [(0,0)]
        vertical_faces_lst = [getFaceListMatch(commons[-1],vector(0.,0.,1.))]



        bd = [body]
        bd_cut = [body_cut]
        matrix = {}
        matrix[(0,0)] = 0
        for k in range(-NbodiesY,NbodiesY):
            for j in range(-NbodiesX,NbodiesX):
                if np.abs(j) + np.abs(k) == 0:
                    continue
                if cell_outside_circle(bounding_box,lTransfX*j + lTransfY*k,rcut):
                    continue

                #if cell_outside_circle2(org_bounding_box,lTransfX*j + lTransfY*k,rcut):
                    #continue
                #if cell_outside_circle1(vlist,lTransfX*j + lTransfY*k,rcut):
                    #continue

                bdname = "bd_%d_%d"%(j,k)
                matrix[(j,k)] = len(bd)
                bd.append(A.addObject("Part::Feature","Body"))
                bd[-1].Shape = shp
                bd[-1].Label = bdname
                bd[-1].Placement = App.Placement(lTransfX*j + lTransfY*k,App.Rotation(vector(0,0,1),0),vector(0,0,0))

                if cut:
                    bdcutname = "bdcut_%d_%d"%(j,k)
                    bd_cut.append(A.addObject("Part::Feature","Body"))
                    bd_cut[-1].Shape = shp_cut
                    bd_cut[-1].Label = bdcutname

                    # Check if volume of the intersection is non-zero and do not proceed if it is zero.
                    vol = bd[-1].Shape.common(bd_cut[-1].Shape).Volume
                    if vol < 1e-16:
                        A.removeObject(bd_cut.pop().Name)
                        A.removeObject(bd.pop().Name)
                        matrix.pop((j,k))
                        continue


                    commons.append(A.addObject("Part::MultiCommon","Common"))
                    commons[-1].Shapes = [bd_cut[-1],bd[-1]]
                    commons[-1].Refine = True
                    vertical_faces_lst.append(getFaceListMatch(commons[-1],vector(0.,0.,1.)))
                    coords.append((j,k))
                #if commons[-1].Label == "Common001":
                #    App.Console.PrintMessage("Column Building Complete! %s \n"%(commons[-1].Label))
                #    App.Console.PrintMessage("j = %d, k = %d \n"%(j,k))
                #    print("Shape BB:",bd[-1].Shape.BoundBox)
                #    print("Min Radii:",np.min([np.sqrt(item.Point[0]**2. + item.Point[1]**2.) for item in bd[-1].Shape.Vertexes]))
                E.recompute()

        E.recompute()


        ilist = []
        for i in range(len(commons)):
            if (commons[i].Shape.Volume < 1e-16):
                labs = (commons[i].Shapes[0].Name,commons[i].Shapes[1].Name)
                A.removeObject(commons[i].Label)
                A.recompute()
                A.removeObject(labs[0])
                A.recompute()
                A.removeObject(labs[1])
                A.recompute()
                ilist.append(i)

        tmp = commons + []
        commons = [tmp[i] for i in range(len(commons)) if i not in ilist]

        tmp = vertical_faces_lst + []
        vertical_faces_lst = [tmp[i] for i in range(len(vertical_faces_lst)) if i not in ilist]

        tmp = coords + []
        coords = [tmp[i] for i in range(len(coords)) if i not in ilist]
        
        column = A.addObject('PartDesign::Body','ColumnAlias')

        colPadPlane = column.newObject("PartDesign::Plane","colPadPlane")

        colPadPlane.AttachmentOffset = App.Placement(vector(0.,0.,1.)*(packing_height - hcol/2.) + vector(0.,0.,column_Z_offset_bottom),  App.Rotation(0.,0.,0.))

        colPadPlane.Support = [A.getObject("XY_Plane002"),""]
        colPadPlane.MapMode = "FlatFace"
        colPadPlane.MapReversed = True
        colPadPlane.Visibility = False
        E.recompute()

        column_profile = column.newObject("Sketcher::SketchObject","column_profile")
        #s = A.getObject("column_profile")
        s = column_profile
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

        col_pad = column.newObject('PartDesign::Pad',"col_pad")
        #col_pad = A.getObject("col_pad")
        col_pad.Profile = s
        col_pad.Length = hcol
        col_pad.Midplane = 1
        E.recompute()

        s = column.newObject("Sketcher::SketchObject","lid_sketch")
        #s = A.getObject("lid_sketch")
        s.Support = (col_pad,["Face3"])
        s.MapMode = "FlatFace"
        s.Visibility = False

        # Line -- 0 (Circle)
        s.addGeometry(Part.Circle(vector(0.,0.,0.),vector(0.,0.,1.),rcol),False)
        s.addConstraint(constraint("Coincident",0,3,-1,1))
        s.addConstraint(constraint("Radius",0,rcol))
        s.setDatum(1,App.Units.Quantity(convertDistanceUnits(rcol)))
        E.recompute()

        lid = column.newObject("PartDesign::Pad","lid")
        #lid = A.getObject("lid")
        lid.Profile = s
        lid.Length = tcol
        lid.Reversed = True

        lid_sketch = s


        k = 0
        hname = "inletSketch_%d"%(k)
        s = column.newObject("Sketcher::SketchObject",hname)
        inlet_sketch = s
        #s = A.getObject(hname)
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
        gpkt = column.newObject("PartDesign::Pocket",hname)
        #gpkt = A.getObject(hname)
        gpkt.Profile = s
        gpkt.Length = tcol

        k = 0
        hname = "inletSketch_1"
        s = column.newObject("Sketcher::SketchObject",hname)
        inletsketch1 = s
        #s = A.getObject(hname)
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
        hole = column.newObject("PartDesign::Pad",hname)
        #hole = A.getObject(hname)
        hole.Profile = s
        hole.Length = tcol/2. + hInlet
        hole.Reversed = False
        #inlet = A.getObject(hname)

        column.ViewObject.Transparency = 70

        E.recompute()
        E.recompute()
        Gui.SendMsgToActiveView("ViewFit")

        commonLabsList = [i.Label for i in commons]
        common_indices = np.array([0] + [int(i.replace("Common","")) for i in commonLabsList[1:]],dtype = np.int)

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
                fusionObject = A.addObject("Part::Fuse",fusionName)
                #fusionObject = A.getObject(fusionName)
                fusionObject.Base = clist[0]
                fusionObject.Tool = clist[1]
                fusions = fusions + [fusionObject]
                fusionObject = None
                j += 1
                E.recompute()

                if len(clist) > 2:
                    for kk in range(2,len(clist)):
                        fusionName = "fusion%d"%(j)
                        fusionObject = A.addObject("Part::Fuse",fusionName)
                        #fusionObject = A.getObject(fusionName)
                        fusionObject.Base = clist[kk]
                        fusionObject.Tool = fusions[-1]
                        fusions = fusions + [fusionObject]
                        fusionObject = None
                        j += 1
                        E.recompute()

                fusions[-1].Label = "Fusion" + arrlist[xi]

        if apply_compound:
            packing = A.addObject("Part::Compound","Packing")
            packing.Links = commons #+ [column]
            E.recompute()

        # create column alias
        column_final = A.addObject("Part::Feature","Column")
        column_final.Shape = Part.getShape(column,"",needSubElement=True,refine=True)
        column_final.ViewObject.Transparency = 70


        self.MainShp = packing.Shape
        #self.objects = commons + [column,packing]
        self.objects = [column_final]
        #self.Object.Children = [column,packing]
        self.Children = [column_final] #[column,packing]
        obj.Shape = self.MainShp



        # Cleanup 
        A.removeObject(packing.Label)
        A.recompute()
        column.removeObjectsFromDocument()
        A.removeObject(column.Label)
        A.recompute()

        for O in commons:
            common_child_objects = [item for item in O.Shapes]
            labs = [item.Name for item in O.Shapes]
            A.removeObject(O.Label)
            A.recompute()
            for ishape in common_child_objects:
                ilab = ishape.Name
                if "PartDesign::" in ishape.TypeId:
                    ishape.removeObjectsFromDocument()
                A.removeObject(ilab)
        A.recompute()


        self.columnCreated = True
        self.resetNeeded = False


    #def attach(self,obj):
        #obj.addExtension("App::OriginGroupExtensionPython")
        #obj.Origin = FreeCAD.ActiveDocument.addObject("App::Origin","Origin")

    
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

        # Break compounds first
        ilist = []
        for i,item in enumerate(self.objects):
            if item.TypeId == "Part::Compound":
                ilist.append(i)
                A.removeObject(item.Label)
        A.recompute()
        tmp = self.objects + []
        self.objects = [tmp[i] for i in range(len(tmp)) if i not in ilist]

        #for O in self.objects:
        while len(self.objects) > 0:
            O = self.objects.pop()
            if O.TypeId == "Part::MultiCommon":
                common_child_objects = [item for item in O.Shapes]
                labs = [item.Name for item in O.Shapes]
                A.removeObject(O.Label)
                A.recompute()
                for ishape in common_child_objects:
                    ilab = ishape.Name
                    if "PartDesign::" in ishape.TypeId:
                        ishape.removeObjectsFromDocument()
                    A.removeObject(ilab)
                A.recompute()
            elif "PartDesign::" in O.TypeId:
                O.removeObjectsFromDocument()
                A.removeObject(O.Label)
                A.recompute()
            else:
                A.removeObject(O.Label)
                A.recompute()

        self.resetNeeded = False
        self.columnCreated = False

        
class ViewProviderColumn:
    
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
