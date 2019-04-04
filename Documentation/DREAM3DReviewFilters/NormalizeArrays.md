# Normalize Attribute Arrays #

## Group (Subgroup) ##

DREAM3D Review (Statistics)

## Description ##

This **Filter** performs a _normalization_ operation on any number of scalar input **Attribute Arrays**.  The user may select to either rescale the arrays to a specified range or standardize the arrays to have mean of 0 and unit variance.  The rescaling operation applies to following to each tuple of the input **Attribute Arrays**:

\f[ X' = a + \frac{\left(X-X_{\min}\right)\left(b-a\right)}{X_{\max} - X_{\min}} \f] 

where \f$ X \f$ is the original array value, \f$ X_{min} \f$ is the original array minimum, \f$ X_{max} \f$ is the original array maximum, \f$ X' \f$ is new array value, and \f$ [a, b] \f$ is the new range.  The standardization operation applies the following to each tuple:

\f[ X' = \frac{X - \mu}{\sigma} \f] 

where \f$ X \f$ is the original array value, \f$ \mu \f$ is the original array mean, \f$ \sigma \f$ is the original array standard deviation, and \f$ X' \f$ is the new array value.

The output **Attribute Arrays** will all be in 64-bit floating point precision (doubles).  The user may opt to use a mask to ignore certain points; where the mask is _false_, the point will not be included when computing \f$ X_{min} \f$, \f$ X_{max} \f$, \f$ \mu \f$, or \f$ \sigma \f$.  If a mask is used, the user must also provide a default value to initialize ignored points in the output **Attribute Arrays**.

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Operation Type | Enumeration | Whether to rescale to a range standardize the selected **Attribute Arrays** |
| Postfix | string | Postfix to attach to the output **Attribute Arrays** |
| Use Mask | bool | Whether to use a boolean mask array to ignore certain points flagged as _false_ from the algorithm |
| Default Masked Value | double | Default value to use for masked points in the output **Attribute Arrays**, if _Use Mask_ is checked |

## Required Geometry ###

None

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| Any **Attribute Arrays** | None | Any | Any | **Attribute Arrays** to rescale or standardize |
| **Attribute Array** | Mask | bool | (1) | Specifies if the point is to be counted in the algorithm, if _Use Mask_ is checked |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Attribute Arrays** | Original name + postfix | double | (1) | Rescaled or standardize **Attribute Arrays** |

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
