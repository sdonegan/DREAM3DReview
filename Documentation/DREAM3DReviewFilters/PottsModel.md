# Potts Model #

## Group (Subgroup) ##


DREAM3DReview (Coarsening)

## Description ##

This **Filter** simulates grain growth using the Potts model.  The Potts model is a generalization of the Ising model to \f$ S \f$ states, or _spins_, on a regular lattice.  This version of the Potts model functions in the context of **Feature** Ids on an **Image Geometry**, and thus will have the effect of coarsening **Features**; additionally, both 2D and 3D **Image Geometries** may be utilized.  The **Feature** Ids are coarsened _in place_.  The present implementation uses a Monte Carlo approach; the user may enter the desired number of Monte Carlo iterations to perform.  For each Monte Carlo iteration, the algorithm proceeds as follows:

1. Pick a **Feature** Id (_spin_) at random
2. Compute the change in the Hamiltonian (\f$ \Delta E \f$) associated with flipping the chosen _spin_ to another _spin_ (_candidate_)
3. If the energy change computed in step 2 is negative, then flip the _spin_ to the _candidate_
4. If the energy change is positive, then flip the spin if the following is true:

\f[ r < e^{\frac{-\Delta E}{kT}} \f]  
  
  where \f$ r \f$ is a random number on the interval \f$ [0, 1) \f$, \f$ \Delta E \f$ is the energy change computed in step 2, \f$ k \f$ is Boltzmann's constant (1.380 x 10<sup>-23</sup> J/K), and \f$ T \f$ is a user-defined temperature value (in Kelvin).
  
5. Increment the time step by \f$ 1 / N \f$, where \f$ N \f$ is the total number of _lattice sites_ (i.e., the number of **Cells** in the **Image Geometry**, or the number of _trues_ in the mask, if _Use Mask_ is checked); this ensures that, on average, a spin flip is attempted once per lattice site for each iteration

The Hamiltonian used for this Potts model implementation is _isotropic_; thus, all boundaries are treated similarly.  The major driving force that encourages flipping is the local neighborhood around each spin.  Specifically, the energy change computed in step 2 is equivalent to the following:

\f[ \Delta E = \frac{1}{2} \sum_{n} \delta(s_{n}, s_{c}) - \delta(s_{n}, s_{i})  \f]

where \f$ s_{i} \f$ is the local spin attempting to flip, \f$ n \f$ are the _valid_ neighbors of \f$ s_{i} \f$, \f$ s_{c} \f$ is the candidate spin to which \f$ s_{i} \f$ may flip, and \f$ \delta \f$ is the Kronecker delta function.  The neighborhood used for each \f$ s_{i} \f$ is _fully connected_; in 2D, this equates to 8 neighbors, and in 3D this equates to 26 neighbors.  If _Periodic Boundaries_ is checked, then **Cells** on the surface will consider **Cells** on the opposite surface as valid neighbors.  If _Periodic Boundaries_ is not checked, the surface **Cells** will not consider **Cells** on the opposite surface as valid neighbors.

The user may specify a mask to ignore certain points from the simulation; lattice sites where the mask is _false_ will never be selected to potentially flip, may not be selected as candidate spins for other lattice sites, and will not be considered valid neighbors when computing the energy change for the Hamiltonian.  Masked points therefore act as sites where boundary motion is pinned; this phenomenon is known as _Zener pinning_.  Note that spins with the Id value 0 (**Feature** Id = 0) also share this behavior (**Feature** Id 0 will act the same as if a mask value of _false_ is at that position, even if no mask is being used).

This implementation of the Potts model uses several techniques to speed up the overall computation.  First, after step 1, if the selected site's neighbors all have the same spin as the local site (i.e., the local site is completely within a grain), no spin flip will be attempted.  Additionally, when selecting candidate spins, only spins that are among neighboring sites may be selected.  

## Parameters ##

| Name | Type | Description |
|------|------|------|
| Iterations | int32_t | Number of Monte Carlo time steps |
| Temperature | double | Temperature value to use when computing \f$ kT \f$, in Kelvin |
| Periodic Boundaries | bool | Whether to enforce periodic boundary conditions when computing neighbors |
| Use Mask | bool | Whether to use a boolean mask array to ignore certain points flagged as _false_ from the algorithm |

## Required Geometry ##

Image

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|-------------|---------|-----|
| **Cell Attribute Array** | FeatureIds | int32_t | (1) | Specifies to which **Cell** each **Feature** belongs |
| **Cell Attrubute Array** | Mask | bool | (1) | Specifies if the point is to be counted in the algorithm, if _Use Mask_ is checked |

## Created Objects ##

None

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this plugin.

## DREAM3D Mailing Lists ##

If you need more help with a filter, please consider asking your question on the DREAM3D Users mailing list:
