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

import FreeCAD
import FreeCADGui
import FreeCAD as App
from FreeCAD import Base
import Part
import PartDesign
import PartDesignGui
import math
import sys
import Sketcher
import numpy as np
#import File
import ImportGui
import Mesh
from PySide import QtGui
from PySide import QtCore

