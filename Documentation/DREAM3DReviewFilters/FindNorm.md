# Find Norm #

## Group (Subgroup) ##

DREAM3D Review (Statistics)

## Description ##

This **Filter** computes the p<sup>th</sup> norm of an **Attribute Array**.  Specifically, for each tuple of the array, the following is computed:

\f[ \left\| \mathbf{x} \right\| _p := \bigg( \sum_{i=1}^n \left| x_i \right| ^p \bigg) ^{1/p} \f]   

where \f$ n \f$ is the number of components for the **Attribute Array**.  When \f$ p = 2 \f$, this results in the _Euclidean norm_; when \f$ p = 1 \f$, this results in the _Manhattan norm_ (also called the _taxicab norm_).  The p-space value may be any real number greater than or equal to zero.  When \f$ 0 \leq p < 1 \f$, the result may not strictly be a _norm_ in the exact sense.  Additionally, when \f$ p = 0 \f$, the result is simply the number of components for the **Attribute Array**.     

_Note:_ If the input array is a scalar array, the output array will contain the same values as the input array, but in 32-bit floating point precision.
 
## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| p-Space Value | float | p-Value used for computing the norm |

## Required Geometry ###

Any

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| Any **Attribute Array** | None | Any | Any | Input **Attribute Array** for computing the norm |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Attribute Array** | Norm | float | (1) | Norm of the input **Attribute Array** |

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users
