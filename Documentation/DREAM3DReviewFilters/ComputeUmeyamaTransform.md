# Compute Umeyama Transform #

## Group (Subgroup) ##

DREAM3D Review (Registration)

## Description ##

This **Filter** computes the Umeyama transform between two unstructured **Geometries**.  An "unstructured" **Geometry** is any geometry that requires explicit definition of **Vertex** positions.  Specifically, **Vertex**, **Edge**, **Triangle**, **Quadrilateral**, and **Tetrahedral** **Geometries** are considered unstructured **Geometries**.  The Umeyama transform estimates the affine transformation that maximizes the _similarity_ (specifically, the covariance) between two n-dimensional point sets.  For two unstructured **Geometries**, this has the effect of estimating, in closed form, the affine transformation that best brings the **Vertices** of the moving **Geometry** into the reference frame of the **Vertices** of the fixed **Geometry**.  The algorithm may optionally incorporate scaling into the computation of the affine transformation matrix.

The moving and fixed unstructured **Geometries** do not need to be the same type, since only the **Vertex** positions are considered; however, the two unstructured **Geometries** need to share the same _number_ of **Vertices**.  This **Filter** will create an **Attribute Matrix** located in the **Data Container** of the moving **Geometry** in which to store the transformation matrix.  Note that the stored transformation matrix will be in _row major_ order.

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Use Scaling | bool | Whether to include scaling into the transformation matrix |

## Required Geometry ###

Any unstructured **Geometry**

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| Moving **Data Container** | None | N/A | N/A | The **Data Container** holding the moving **Geometry** |
| Fixed **Data Container** | None | N/A | N/A | The **Data Container** holding the fixed **Geometry** |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Generic Attribute Matrix** | TransformationData | Generic | N/A | The **Attribute Matrix** used to store the computed transformation matrix  |
| **Attribute Array** | TransformationMatrix | float | (2, 4) | The computed transformation matrix |

## References ## 

[1] Least-squares estimation of transformation parameters between two point patterns, S. Umeyama, PAMI 1991, DOI: 10.1109/34.88573

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
