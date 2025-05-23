---
layout: page
title: SCALE DG regional and global model (SCALE-DG)
---

SCALE-DG is a three-dimensional (3D) regional and global atmospheric model
in which both of the horizontal and the vertical discretizations are based on the discontinuous Galerkin (DG) method.

<!-- Now, only the dynamical process has been implemented in addition to simple turbulent models.
In the near feature, using the [SCALE library](https://scale.riken.jp) physical processes such as cloud microphysics and radiation will be supported. -->

## Brief Description

### Features

* Support both 3D regional and global simulations only in SCALE-DG
  * This feature is provided by a model framework provided by our DG library
* Execute parallel simulations in supercomputers  
  * Process and thread parallelizations are performed by MPI library (we use 2D domain compositions) and OpenMP, respectively.
* Generate initial data for idealized experiments
* Treat NetCDF files as input and output data

### Dynamical core
* Governing equation: 3D fully compressible non-hydrostatic equations
* Grid system: hexahedral finite element and horizontally curvilinear coordinate system
  * For the global mode, cubed sphere coordinates are used
  * Topography is treated using a terrain-following coordinate
* Spatial discretization: nodal discontinuous Galerkin method (e.g., Hesthave and Warburton, 2007)
  * Polynomial order associated with numerical accuracy can be arbitrarily chosen.
  * Numerical flux: Rusanov flux
  * Stabilization mechanisms: Inherent numerical dissipation with numerical fluxes and explicit modal filtering
* Temporal discretization: various type of Runge-Kutta (RK) schemes
  * Full explicit (horizontally and vertically explicit; HEVE)
    * Classical 4s4o RK scheme (where 4s4o means 4th-order and 4 stage)
    * Strong-stability preserving (SSP) RK schemes: 3s3o, 4s3o, 5s3o (Higueras and Roldan, 2018), 10s4o (Ketcheson, 2008)
  * Horizontally explicit and vertically implicit (HEVI)
    * Implicit and explicit RK schemes: ARK232 (Giraldo et al., 2013), ARK324 (Kennedy and Carpenter, 2003)
* Tracer transport
  * Preserving the non-negativity is ensured by a limiter (Light and Durran, 2016) and SSP RK schemes in conjunction with fundamental spatial discretization by nodal DGM.

### Physical processes

* Turbulence process
  * Smagorinsky (1963) and Lilly (1962)-type sub-grid scale model corrected by Brown et al. (1994) and Scotti et al. (1993)

* Cloud microphysics (using the [SCALE library](https://scale.riken.jp))
  * 3-class 1 moment bulk scheme (Kessler 1969) 

* Surface flux scheme
  * Simplified bulk formulation of momentum, heat, and latent flux with constant coefficients for idealized experiments

## Documents

The description of SCALE-DG is available at [Document page]({{ '/documents/' | relative_url }}).

## Example of numerical experiments

Simulation results of standard test cases such as baroclinic wave and Held Suarez experiments 
are shown in [Gallery page]({{ '/gallery/' | relative_url }}).
