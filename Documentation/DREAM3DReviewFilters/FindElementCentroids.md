# Find Element Centroids #

## Group (Subgroup) ##

DREAM3D Review (Geometry)

## Description ##

This **Filter** computes the centroids of all **Elements** in the supplied **Geometry**.  Currently, **Edge**, **Triangle**, **Quadrilateral**, and **Tetrahedral** **Geometries** are supported.  The user must select where to store the computed centroids array; the **Geometry** used will belong to the **Data Container** in the path for the selected for the centroids.

**Element** centroids represent fundamental topological information about the **Geometry**, and thus will be stored within the **Geometry** object when computed.  This **Filter** will additionally store the centroids as an **Attribute Array** within the data structure.  The array represents the (x, y, z) coordinates of the centroid.  This **Filter** will enforce that the centroids are stored in an **Attribute Matrix** that matches the fundamental **Element** type of the selected **Geometry**:

| Geometry Type | Destination Attribute Matrix |
|----------|-----------|
| Edge | Edge |
| Triangle | Face |
| Quadrilateral | Face |
| Tetrahedral | Cell | 

This **Filter** may optionally create a new **Data Container** that holds a **Vertex Geometry** whose point positions are the computed centroids.  A **Vertex Attribute Matrix** will be created for this new **Geometry** to store **Vertex**-level information. 

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Create Vertex Geometry from Centroids  | bool | Whether to create a **Vertex Geometry** from the computed centroids |
| Vertex Data Container | string | Name of the **Data Container** to hold the new centroids **Vertex Geometry**, if _Create Vertex Geometry from Centroids_ is checked |
| Vertex Attribute Matrix | string | Name of the **Vertex Attribute Matrix** for the new **Vertex Data Container**, if _Create Vertex Geometry from Centroids_ is checked  |

## Required Geometry ###

**Edge**, **Triangle**, **Quadrilateral**, or **Tetrahedral**

## Required Objects ##

None

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Attribute Array** | Centroids | float | (3) | Centroid coordinates of each **Element** in the supplied **Geometry** |

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users
