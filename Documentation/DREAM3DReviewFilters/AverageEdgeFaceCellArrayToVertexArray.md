# Average Edge/Face/Cell Array to Vetex Array #

## Group (Subgroup) ##

DREAM3D Review (Statistics)

## Description ##

This **Filter** averages the selected **Edge**, **Face**, or **Cell** array onto the **Vertices** of a mesh-like **Geometry**.  Mesh-like **Geometries** are those that require explicit definition of connectivity between **Vertices** and **Elements**, specifically **Edge**, **Triangle**, **Quadrilateral**, and **Tetrahedral** **Geometries**.  The value stored on each **Vertex** is the _arithmetic mean_ of the values for all **Elements** (i.e., the **Edges**, **Faces**, or **Cells**) that connect to the **Vertex**.

The selected **Edge**, **Face**, or **Cell** array may be scalar or multicomponent; if the array is multicomponent, then each component value is averaged individually.  Therefore, the output averaged **Vertex** array will have the same component dimensions as the input array.  The input array may also be of any primitive type; however, the output averaged array will always be of type _float_.

## Parameters ##

None

## Required Geometry ###

Any mesh-like **Geometry**

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Edge/Face/Cell Attribute Array** | None | Any | Any | The **Edge**, **Face**, or **Cell** array to average |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Vertex Attribute Array** | None | float | Any | The averaged **Vertex** array |

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users

