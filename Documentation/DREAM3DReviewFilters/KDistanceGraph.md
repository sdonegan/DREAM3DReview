# K Distance Graph #

## Group (Subgroup) ##

DREAM3D Review (Clustering)

## Description ##

This **Filter** computes the distance between each point and its k<sup>th</sup> nearest neighbor.  For example, if \f$ k = 1 \f$, this **Filter** will store the distance bewteen each point and its closest nearest neighbor (i.e., the distance that is smallest among all pair-wise distances).  The user may select from a number of options to use as the distance metric.  When sorted smallest-to-largest, the k distance array forms a graph that is useful for estimating parameters in some clustering algorithms, such as [DBSCAN](@ref dbscan).  The user may opt to use a mask array to ignore points in the distance computation; these points will contain a distance value of 0 in the output array.

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| K<sup>th</sup> Nearest Neighbor | int32_t | Which nearest neighbor for which to compute the distance |
| Distance Metric | Enumeration | The metric used to determine the distances between points |
| Use Mask | bool | Whether to use a boolean mask array to ignore certain points flagged as _false_ from the algorithm |

## Required Geometry ###

None

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| Any **Attribute Array** | None | Any| Any | The **Attribute Array** for which the compute the k<sup>th</sup> nearest neighbor distances |
| **Attribute Array** | Mask | bool | (1) | Specifies if the point is to be counted in the algorithm, if _Use Mask_ is checked |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Attribute Array** | KDistance | double | (1) | Distance to the k<sup>th</sup> nearest neighbor for each point  |

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users

