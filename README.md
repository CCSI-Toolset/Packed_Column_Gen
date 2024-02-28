# PackedColumnGen
PackedColumnGen is a tool for building CAD models of absorber columns with customized structured packing designs. The tool generates these designs from a user-defined set of input parameters. This tool is designed as [Workbench](https://wiki.freecad.org/Workbenches) for [FreeCAD](https://www.freecad.org/), a free and open-source CAD software with a user-friendly interface. The tool currently supports the following features:

* Parametric construction of three-dimensional columns with structured packing.
* Addition of cooling channels to create process intensified structured packing.
* Built-in functions that extend support for scripted packing generation.

Refer to the [user guide](./doc/README.md) to learn more and get started with example cases.

## Installation

### Prerequisites
FreeCAD version 0.19 or newer needs to be installed before running the setup. The latest version can be downloaded [here](https://www.freecad.org/).

### Installation using the Setup file

The tool source code is packaged into the `Setup.exe` executable file, which extracts to the default FreeCAD modules directory. On Windows, the module directory is `%APPDATA%\FreeCAD\Mod\`, which is usually `C:\Users\username\Appdata\Roaming\FreeCAD\Mod\`. Restart FreeCAD after extraction for changes to take effect and the module should be accessible as `Packed Columns` workbench in FreeCAD as shown below.

![Workbench demonstration](./doc/Images/Workbench.png)

### Manual installation
1. Find the FreeCAD workbench directory by typing `App.getUserAppDataDir()` in FreeCAD [Python console](https://wiki.freecad.org/Python_console).
2. Navigate to the workbench directory and copy or extract the contents of this repository to a new folder within this directory.
3. **IMPORTANT** Rename the folder to "packedColumnGen" (case sensitive).
4. Restart FreeCAD for changes to take effect. The workbench `Packed Columns` should now in the list of workbenches in FreeCAD.

## Authors
 > [Yash Girish Shah](mailto:yashgirish.shah@netl.doe.gov)   
 > [Grigorios Panagakos](mailto:gpanagak@andrew.cmu.edu)
