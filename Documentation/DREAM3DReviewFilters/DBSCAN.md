# DBSCAN #

## Group (Subgroup) ##

DREAM3D Review (Clustering)

## Description ##

This **Filter** applies the DBSCAN (density-based spatial clustering of applications with noise) algorithm to an **Attribute Array**.  DBSCAN is a _clustering algorithm_ that assigns to each point of the **Attribute Array** a _cluster Id_; points that have the same cluster Id are grouped together more densely (in the sense that the _distance_ between them is small) in the data space (i.e., points that have many nearest neighbors will belong to the same cluster).  The user may select from a number of options to use as the distance metric.  Points that are in sparse regions of the data space are considered "outliers"; these points will belong to cluster Id 0.  Additionally, the user may opt to use a mask to ignore certain points; where the mask is _false_, the points will be categorized as outliers and placed in cluster 0.  The algorithm requires two parameters: a _neighborbood region_, called epsilon; and the minimum number of points needed to form a cluster.  The algorithm, in pseudocode, proceeds as follows:

    for each point p in dataset
    {
      cluster = 0
      if p has been visited
      {
        continue to next point
      }
      mark p as visited
      neighbor_points = all points within epsilon distance from p
      if the number of neighbor_points < minimum number of points
      {
        mark p is outlier (cluster Id = 0)
      }
      else
      {
        cluster++
        expand_cluster(p, neighbor_points, cluster, epsilon, minimum number of points)
      }
    }

The `expand_cluster()` function, in pseudocode, is as follows:

    expand_cluster(p, neighbor_points, cluster, epsilon, minimum number of points)
    {
      add p to cluster
      for each point p_prime in neighbor_points
      {
        if p_prime has not been visited
        {
          mark p_prime as visited
          neighbor_points_prime = all points within epsilon distance from p_prime
          if the number of neighbor_points_prime >= minimum number of points
          {
            adjoin neighbor_points_prime to neighbor_points
          }
        }
        if p_prime is not a member of any cluster
        {
          add p_prime to cluster
        }
      }
    }

An advantage of DBSCAN over other clustering approaches (e.g., [k means](@ref kmeans)) is that the number of clusters is not defined _a priori_.  Additionally, DBSCAN is capable of finding arbitrarily shaped, nonlinear clusters, and is robust to noise.  However, the choice of epsilon and the minimum number of points affects the quality of the clustering.  In general, a reasonable rule of thumb for choosing the minimum number of points is that it should be, at least, greater than or equal to the dimensionality of the data set plus 1 (i.e., the number of components of the **Attribute Array** plus 1).  The epsilon parameter may be estimated using a _k distance graph_, which can be computed using [this Filter](@ref kdistancegraph).  When computing the k distance graph, set the k nearest neighbors value equal to the minimum number of points intended for DBSCAN.  A reasonable choice of epsilon will be where the graph shows a strong bend.  If using this approach to help estimate epsilon, remember to use the same distance metric in both **Filters**!  An alternative method to choosing the two parameters for DBSCAN is to rely on _domain knowledge_ for the data, considering things like what neighbor distances between points make sense for a given metric.  
    
A clustering algorithm can be considered a kind of segmentation; this implementation of DBSCAN does not rely on the **Geometry** on which the data lie, only the _topology_ of the space that the array itself forms.  Therefore, this **Filter** has the effect of creating either **Features** or **Ensembles** depending on the kind of array passed to it for clustering.  If an **Element** array (e.g., voxel-level **Cell** data) is passed to the **Filter**, then **Features** are created (in the previous example, a **Cell Feature Attribute Matrix** will be created).  If a **Feature** array is passed to the **Filter**, then an **Ensemble Attribute Matrix** is created.  The following table shows what type of **Attribute Matrix** is created based on what sort of array is used for clustering:

| Attribute Matrix Source             | Attribute Matrix Created |
|------------------|--------------------|
| Generic | Generic |
| Vertex | Vertex Feature | 
| Edge | Edge Feature |
| Face | Face Feature | 
| Cell | Cell Feature| 
| Vertex Feature | Vertex Ensemble |
| Edge Feature | Edge Ensemble |
| Face Feature | Face Ensemble |
| Cell Feature | Cell Ensemble|
| Vertex Ensemble | Vertex Ensemble |
| Edge Ensemble | Edge Ensemble |
| Face Ensemble | Face Ensemble |
| Cell Ensemble | Cell Ensemble|

## Parameters ##

| Name | Type | Description |
|------|------|-------------|
| Epsilon | float | The epsilon-neighborbood around each point is queried |
| Minimum Number of Points | int32_t | The minimum number of points needed to form a _dense region_ (i.e., the minimum number of points needed to be called a cluster) |
| Distance Metric | Enumeration | The metric used to determine the distances between points |
| Use Mask | bool | Whether to use a boolean mask array to ignore certain points flagged as _false_ from the algorithm |

## Required Geometry ###

None

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| Any **Attribute Array** | None | Any| Any | The **Attribute Array** to cluster |
| **Attrubute Array** | Mask | bool | (1) | Specifies if the point is to be counted in the algorithm, if _Use Mask_ is checked |

## Created Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Attribute Matrix** | ClusterData | Feature/Ensemble | N/A | The **Attribute Matrix** in which to store information associated with the created clusters |
| **Attribute Array** | ClusterIds | int32_t | (1) | Specifies to which cluster each point belongs  |

## References ## 

[1] A density-based algorithm for discovering clusters in large spatial databases with noise, M. Ester, H.P. Kriegel, J. Sander, and X. Xu, Proceedings of the Second International Conference on Knowledge Discovery and Data Mining, pp. 226-231, 1996.

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
https://groups.google.com/forum/?hl=en#!forum/dream3d-users
