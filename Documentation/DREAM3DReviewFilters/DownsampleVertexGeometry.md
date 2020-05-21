Downsample Vertex Geometry
=============

## Group (Subgroup) ##
Sampling (Geometry)

## Description ##
This **Filter** downsamples (i.e., reduces the number of points in) an input **Vertex Geometry**.  Three approaches for downsampling are provided by the **Filter**:

- Remove every Nth point:
    - Each Nth point, in the order the points are stored, is removed from the **Vertex Geometry**.  The parameter N is set by the user.  For example, if N=2, half of the points in the **Vertex Geometry** would be removed.
- Remove a fixed fraction of random points:
    - A user-defined fraction of points are removed at random from the **Vertex Geometry**.  For example, if the user selects a fraction of 0.2, 20% of the points would be removed from the **Vertex Geometry** at random.
- Downsample the geometry on a grid:
    - The user defines the resolution (i.e., voxel spacing) of a structure rectilinear grid.  This sampling grid is overlaid on the **Vertex Geoemtry**.  All points that fall in a given voxel are averaged together, producing a single point at each voxel whose (x,y,z) coordinates are the mean of the coordinates of all points in that voxel.  Additionally, any **Attribute Arrays** that exist on these points are averaged together within each sampling voxel. 

Note that the **Vertex Geometry** is modified in place.

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Downsample Type | Enumeration | The downsampling approach to use, either *Remove Every Nth Point*, *Remove a Fixed Random Fraction of Points*, or *Downsample the Geometry on a Grid* |
| Decimation Frequency | int | The cadence at which to remove points, if *Remove Every Nth Point* is chosen |
| Fraction to Remove | float | The fraction of points to remove, if *Remove a Fixed Random Fraction of Points* is chosen |
| Grid Resolution | float 3x | The resolution of the downsampling grid, if *Downsample the Geometry on a Grid* is chosen |

## Required Geometry ###

Vertex

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Data Container** | None | N/A | N/A | **Data Container** holding the **Vertex Geometry** to downsample |

## Created Objects ##

None

## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users