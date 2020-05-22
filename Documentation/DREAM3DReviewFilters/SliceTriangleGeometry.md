# Slice Triangle Geometry #

## Group (Subgroup) ##

Sampling (Geometry)

## Description ##

This **Filter** slices an input **Triangle Geometry**, producing an **Edge Geometry**.  The user can control the direction along which the slicing occurs, the range over which to slice (either the entire range of the geoemtry or a specified subregion), and the spacing bewteen slices.  The total area and perimieter of each slice is also computed and stored as an attribute on each created slice.

Additionally, if the input **Triangle Geometry** is labeled with an identifier array (such as different regions or features), the user may select this array and the resulting edges will inherit these identifiers.


## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Slice Direction (ijk) | float (3x) | Direction on which to slice the **Triangle Geometry** |
| Slice Range | Enumeration | Type of slice range to use, either *Full Range* or *User Defined Range* |
| Slice Spacing | float | Spacing between slices |
| Have Region Ids | bool | Whether to supply an id array that propagates to the created edges |

## Required Geometry ###

Triangle

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Data Container** | None| N/A | N/A | **Data Container** that contains the **Triangle Geometry** to be sliced |
| **Triangle Attribute Array** | None | int32_t | (1) | Optional identifier array, if *Have Region Ids* is selected |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Data Container** | SliceDataContainer | N/A | N/A | **Data Container** for the resulting **Edge Geoemtry** |
| **Attribute Matrix** | EdgeData | Edge | N/A | **Attribute Matrix** to store information about the created edges |
| **Edge Attribute Array** | SliceIds | int32_t | (1) | Identifies the slice to which each edge belongs |
| **Edge Attribute Array** | RegionIds | int32_t | (1) | Identifies the region from which each edge came from in the original **Triangle Geoemtry**, if *Have Region Ids* is selected |
| **Attribute Matrix** | SliceData | Edge Feature | N/A | **Attribute Matrix** to store information about the created edges |
| **Feature Attribute Array** | SliceAreas | Feature | (1) | The total area (i.e., summed area of each enclosed polygon) of a given slice |
| **Feature Attribute Array** | SlicePerimeters | Feature | (1) | The total perimeter (i.e., summed edge length) of a given slice |

## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users