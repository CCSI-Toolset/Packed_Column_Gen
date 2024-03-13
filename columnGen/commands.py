# ***************************************************************************
# *                                                                         *
# *   Copyright: See License.md file                                        *
# *                                                                         *
# *   Authors: Yash Girish Shah, Grigorios Panagakos                        *
# *                                                                         *
# ***************************************************************************

from .preamble import *
from .commonFunctions import getWorkbenchFolder as getWorkbenchFolder
import columnGen.packing3D as packing3D
import columnGen.packing2D as packing2D

class PackedColCommand2D:
    
    def GetResources(self):
        return {'Pixmap'  : getWorkbenchFolder() + '/Resources/Icons/settings_icon2D.svg',
                'Accel' : "Shift+O", 
                'MenuText': "2D Packed Column",
                'ToolTip' : "Define Parameters for a 2D Packed Column"}

    def Activated(self):
        obj=FreeCAD.ActiveDocument.addObject("Part::FeaturePython","Packed Column")
        packing2D.PackedColumn2D(obj)
        packing2D.ViewProviderColumn2D(obj.ViewObject, "Packed Column") 
        FreeCAD.ActiveDocument.recompute()
        FreeCADGui.SendMsgToActiveView("ViewFit")
        return

    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
               return False
        else:
               return True

FreeCADGui.addCommand('PackedColumn2D',PackedColCommand2D())

class PackedColCommand:
    
    def GetResources(self):
        return {'Pixmap'  : getWorkbenchFolder() + '/Resources/Icons/settings_icon.svg',
                'Accel' : "Shift+P", 
                'MenuText': "Packed Column",
                'ToolTip' : "Define Parameters for a Packed Column"}

    def Activated(self):
        obj=FreeCAD.ActiveDocument.addObject("Part::FeaturePython","Packed Column")
        packing3D.PackedColumn(obj)
        packing3D.ViewProviderColumn(obj.ViewObject, "Packed Column") 
        FreeCAD.ActiveDocument.recompute()
        FreeCADGui.SendMsgToActiveView("ViewFit")
        return

        
    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
               return False
        else:
               return True

FreeCADGui.addCommand('PackedColumn',PackedColCommand())

       


class GenerateCommand:
    
    def GetResources(self):
        return {'Pixmap'  : getWorkbenchFolder() + '/Resources/Icons/generate.svg',
                'Accel' : "Shift+G", 
                'MenuText': "Generate",
                'ToolTip' : "Generate the Packed Column"}

    def Activated(self):
        selectionList = FreeCADGui.Selection.getSelection()
        for item in selectionList:
            if item.TypeId == "Part::FeaturePython":
                proxy = item.Proxy

                if proxy.Type_ in ["Column2D","Column3D"]:

                    if proxy.resetNeeded:
                        proxy.reset()
                        FreeCAD.ActiveDocument.recompute()

                    if not proxy.columnCreated:
                        proxy.createColumn()
                        FreeCAD.ActiveDocument.recompute()

                    item.touch()
                    FreeCAD.ActiveDocument.recompute()

        FreeCADGui.SendMsgToActiveView("ViewFit")
        return

        
    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
               return False
        else:
               return True

FreeCADGui.addCommand('Generate',GenerateCommand())


class ResetCommand:
    
    def GetResources(self):
        return {'Pixmap'  : getWorkbenchFolder() + '/Resources/Icons/hexahedron.svg',
                'Accel' : "Shift+R", 
                'MenuText': "Reset",
                'ToolTip' : "Reset the Packed Column"}

    def Activated(self):
        Gui = FreeCADGui
        selectionList = Gui.Selection.getSelection()
        for item in selectionList:
            if item.TypeId == "Part::FeaturePython":
                proxy = item.Proxy
                proxy.reset()
                FreeCAD.ActiveDocument.recompute()

        FreeCAD.ActiveDocument.recompute()
        FreeCADGui.SendMsgToActiveView("ViewFit")
        return

        
    def IsActive(self):
        if FreeCAD.ActiveDocument == None:
               return False
        else:
               return True

FreeCADGui.addCommand('Reset',ResetCommand())

