# ***************************************************************************
# *                                                                         *
# *   Copyright: See License.md file                                        *
# *                                                                         *
# *   Authors: Yash Girish Shah, Grigorios Panagakos                        *
# *                                                                         *
# ***************************************************************************



class packedColGenWorkbench (Workbench):

    MenuText = "Packed Columns"
    ToolTip = "A workbench for generating packed columns with structured packing"
    
    def getWorkbenchFolder(self):

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
    

    def __init__(self):
        
        self.__class__.Icon = self.getWorkbenchFolder() + "/Resources/Icons/column_icon.svg"

    def Initialize(self):
        """This function is executed when FreeCAD starts"""
        import FreeCAD
        import columnGen 
        self.list = ["PackedColumn","PackedColumn2D","Generate"]
        self.appendToolbar("ColumnGenerator",self.list) 
        self.appendMenu("ColumnGenerator",self.list) 
        #self.appendMenu(["An existing Menu","My submenu"],self.list) 

    def Activated(self):
        return

    def Deactivated(self):
        return

    def ContextMenu(self, recipient):
        self.appendContextMenu("ColumnGenerator",self.list) 

    def GetClassName(self): 
        return "Gui::PythonWorkbench"
       
Gui.addWorkbench(packedColGenWorkbench())
