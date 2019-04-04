# Extract Internal Surfaces from Triangle Geometry #

## Group (Subgroup) ##

DREAM3D Review (Geometry)

## Description ##

This **Filter** extracts any **Triangles** from the supplied **Triangle Geometry** that contain any _internal nodes_, then uses these extracted **Triangles** to create a new **Data Container** with the reduced **Triangle Geometry**.  This operation is the same as removing all **Triangles** that only lie of the outer surface of the supplied **Triangle Geometry**.  The user must supply a **Vertex Attribute Array** that defines the type for each node of the **Triangle Geometry**.  Node types may take the following values:

| Id Value | Node Type |
|----------|-----------|
| 2 | Normal **Vertex** |
| 3 | Triple Line |
| 4 | Quadruple Point |
| 12 | Normal **Vertex** on the outer surface |
| 13 | Triple Line on the outer surface |
| 14 | Quadruple Point on the outer surface |

This **Filter** has the effect of removing any **Triangles** that only contain **Vertices** whose node Id values are 12, 13, or 14.  In general, this _node type_ array is created when the [original surface mesh is created](@ref quicksurfacemesh).   

It is unknown until run time how the **Geometry** will be changed by removing certain **Vertices** and **Triangles**.  Therefore, this **Filter** requires that a new **Data Container** be created to contain the new internal **Triangle Geometry**.  This new **Data Container** will contain copies of any **Feature** or **Ensemble** **Attribute Matrices** from the original **Data Container**.  Additionally, all **Vertex** and **Face** data will be copied, with tuples _removed_ for any **Vertices** or **Faces** removed by the **Filter**.  The user must supply a name for the new **Data Container**, but all other copied objects (**Attribute Matrices** and **Attribute Arrays**) will retain the same names as the original source.

## Parameters ##

None

## Required Geometry ###

Triangle

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Data Container** | TriangleDataContainer | N/A | N/A | **Data Container** holding the **Triangle Geometry** from which to extract internal surfaces |
| **Vertex Attribute Array** | NodeTypes | int8_t | (1) | Specifies the type of node in the **Geometry** |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| Internal Surfaces **Data Container** | InternalTrianglesDataContainer | N/A | N/A | **Data Container** holding the new **Triangle Geometry** and any copied **Attribute Matrices** and **Attribute Arrays** |

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
