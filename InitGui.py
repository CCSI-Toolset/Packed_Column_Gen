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
