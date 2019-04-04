# Principal Component Analysis #

## Group (Subgroup) ##

DREAM3D Review (Dimensionality Reduction)

## Description ##

This **Filter** performs principal component analysis on a data space defined by the user-selected scalar **Attribute Arrays**.  Specifically, the dimensionality of the data space analyzed is equal to the number of selected arrays, where each array represents one dimension.  The arrays may be of any primitive type, and are not required to share the same primitive type.  This **Filter** uses a matrix approach to compute the principal components; the user may opt to use either the _correlation_ or _covariance_ matrix. The algorithm proceeds as follows:

1. Form a \f$ m \f$ x \f$ n \f$ matrix, \f$ \mathbf{X} \f$, where the columns are the input **Attribute Arrays** (\f$ n \f$ is the number of selected arrays (columns) and \f$ m \f$ is the number of tuples for the arrays (rows))
  * When forming \f$ \mathbf{X} \f$, the values for each of the arrays are cast to double precision
2. If the _correlation_ approach is being used, standardize each column of \f$ \mathbf{X} \f$ such that it has mean of 0 and unit variance; this approach may be useful if each **Attribute Array** has different scale
3. Center \f$ \mathbf{X} \f$ by subtracting the column-wise mean from each row value; call this centered matrix \f$ \mathbf{B} \f$
4. Compute the covariance/correlation matrix, \f$ \mathbf{C} \f$, as follows:

\f[
\mathbf{C} = { 1 \over {m-1} } \mathbf{B}^{*} \otimes \mathbf{B}
\f]

  where \f$ * \f$ is the conjugate transpose; if \f$ m = 1 \f$, then the prefactor denominator is set to 1

5. Perform the eigen decomposition of \f$ \mathbf{C} \f$ to find the eigenvalues and eigenvectors; the eigenvalues are stored in ascending order, and the eigenvectors correspond to the same order as the eigenvalues.

The computed eigenvalues and eigenvectors are stored in a new **Attribute Matrix**.  The number of eigenvalues/eigenvectors computed is \f$ n \f$; the dimensionality of the eigenvectors is also \f$ n \f$.

The user may opt to project the data space to a lower dimensionality using the computed eigenvectors.  A lower dimensionality, \f$ d \f$, must be specified; the \f$ d \f$ eigenvectors that have the highest eigenvalues are used to project the data space.  The \f$ d \f$ eigenvectors form a \f$ d \f$ x \f$ n \f$ matrix that is used to post-multiply each row of the centered matrix \f$ \mathbf{B} \f$.  The result is an **Attribute Array** that is \f$ m \f$ tuples long with \f$ d \f$ dimensions.  It may be useful to visualize this space as a point cloud; this can be accomplished by [creating a Vertex Geometry](@ref creategeometry) using the projected array as coordinates.  Note that a **Vertex Geometry** requires three coordinates for point positions, so if \f$ d < 3 \f$, additional components must be added to the projected data space.  It is possible to [create an array of all zeros](@ref createdataarray) and then [append it to the projected array](@ref combineattributearrays) to get the correct dimensionality.  Also note that **Vertex Geometry** coordinates must be 32-bit floating point values, but the projected space created by this **Filter** will be 64-bit floating point; the different precision can be created by [converting the primitive type](@ref convertdata) of the projected array.  Finally, the data values from the original arrays may be visualized within this projected space by [moving the original Attribute Arrays](@ref movedata) into the **Vertex Attribute Matrix** of the new **Vertex Geometry**.

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Matrix Approach | Enumeration | Whether to use the correlation or covariance matrix approach |
| Project Data Space | bool | Whether to project the input data space to a lower dimensionality based on the principal component eigenvectors |
| Number of Dimensions for Projection | int32_t | The number of dimensions on which to project the data space, if _Project Data Space_ is checked |

## Required Geometry ###

None

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| Any **Attribute Array**  | None | Any | (1) | **Attribute Arrays** that form the data space for the principal component analysis |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Generic Attribute Matrix** | PrincipalComponentAnalysis | Generic | N/A | **Attribute Matrix** in which to store the results of the eigen analysis |
| **Attribute Array**  | PrincipalComponentEigenvalues | double | (1) | Eigenvalues from the principal component analysis |
| **Attribute Array**  | PrincipalComponentEigenvectors | double | (n) | Eigenvectors from the principal component analysis |
| **Attribute Array**  | ProjectedDataSpace | double | (_Number of Dimensions for Projection_) | Projected data space, if _Project Data Space_ is checked |

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users
