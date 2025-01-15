# User Guide

## Getting started

## Example 1: Default structured packing 

###  1. Create a new document `File` &rarr; `New` and select the `Packed Columns`	Workbench as shown below.

![Workbench figure](./Images/Workbench.png)

### 2. In the menu bar, select `ColumnGenerator` &rarr; `Packed Column` to create an instance of packed column.

![Packed Column Selection](./Images/Packed_column_menu.png)

### A `Packed_Column` object appears in the tree view.

![Packed Column treeview](./Images/packed_column_treeview.png)

### 3. Select the packed column instance and review the column parameters in the property editor. Make any changes if needed.

![Property editor](./Images/property_editor.png)

### 4. With the packed column selected, click on `ColumnGenerator` &rarr; `Generate` to generate the column.

![Generate](./Images/menu_generate.png)

### 5. This constructs the default column with default structured packing. For a detailed description of the column parameters refer to the `Parameter Description` section below.

![Default Column](./Images/default_column.png)

## Example 2: Intensified structured packing

### 1. Repeat the steps 1-3 in Example 1. Select the `Enable Intensified Device` option as `true`. This populates two new parameters, `Cooling Channel Width` and `Thickness2` in the Property editor. Review these parameters and make modifications if needed.

![Intensified Settings](./Images/intensified_settings.png)


### 2. With the packed column selected, click on `ColumnGenerator` &rarr; `Generate` to generate the column. 

![Generate](./Images/menu_generate.png)

### 3. This constructs the default column with the default intensified structured packing. 

![Intensified Column](./Images/intensified_column.png)

## Parameter Description

### Column Dimensions

| Parameter                | Description                      |
| ---------------------    | -------------------------------- |
| Column Height            | Height of the column.            |
| Column Radius            | Inner radius of the column body. |
| Column Wall Thickness    | Thickness of the column wall.    |

### Dripping Points

![Dripping point parameters](./Images/DrippingPointsParameters.png)

### Packing Angles

![Packing Angles](./Images/packingAngles.png)

### Packing Dimensions

![Packing Dimensions](./Images/packingDimensions.png)

## Authors
 > [Yash Girish Shah](mailto:yashgirish.shah@netl.doe.gov)   
 > [Grigorios Panagakos](mailto:gpanagak@andrew.cmu.edu)
