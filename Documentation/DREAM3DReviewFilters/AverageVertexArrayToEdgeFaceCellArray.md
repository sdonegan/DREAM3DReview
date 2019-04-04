# Average Vertex Array to Edge/Face/Cell Array #

## Group (Subgroup) ##

DREAM3D Review (Statistics)

## Description ##

This **Filter** averages the selected **Vertex** array onto the **Edges**, **Faces**, or **Cells** of a mesh-like **Geometry**.  Mesh-like **Geometries** are those that require explicit definition of connectivity between **Vertices** and **Elements**, specifically **Edge**, **Triangle**, **Quadrilateral**, and **Tetrahedral** **Geometries**.  The **Element** type on which the averaged array will be stored depends on the **Geometry**: **Edge** **Geometries** result in an **Edge** array, **Triangle** and **Quadrilateral** **Geometries** result in a **Face** array, and **Tetrahedral** **Geometries** result in a **Cell** array.  

The value stored on each **Edge**, **Face**, or **Cell** is the _arithmetic mean_ of the values for all **Vertices** that belong to that **Edge**, **Face**, or **Cell**.  The user may also optionally select for the mean to be _weighted_; in this case, the average values are weighted by the distance of the **Vertex** from the **Edge**, **Face**, or **Cell** _centroid_ such that **Vertices** that are closer to the centroid contribute more strongly to the average value.

The selected **Vertex** array may be scalar or multicomponent; if the array is multicomponent, then each component value is averaged individually.  Therefore, the output averaged **Edge**, **Face**, or **Cell** array will have the same component dimensions as the input array.  The input array may also be of any primitive type; however, the output averaged array will always be of type _float_.

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Perform Weighted Average | bool | Whether the weight the average by the distance of the **Vertices** from the **Edge**, **Face**, or **Cell** centroids |

## Required Geometry ###

Any mesh-like **Geometry**

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Vertex Attribute Array** | None | Any | Any | The **Vertex** array to average |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Edge/Face/Cell Attribute Array** | None | float | Any | The averaged **Edge**, **Face**, or **Cell** array |

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users

