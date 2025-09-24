---
layout: page
title: Gallery / Simulation with regional domain using SCALE-DG (This page is under construction)
---

{% assign files = "test_advect2d.f90,advect2d.conf,Makefile,visualize.py" | split: "," %}

## Mountain wave

This page shows a result of quasi 2D mountain wave using the discontinuous Galerkin method (DGM). 
An example of the setting files is given in rootdir/model/atmosphere3d/test/case/mountain_wave/.

### 1. Description of dynamical core

Please see XX.

### 2. Experimental setup

The experiment setup for a linear and hydrostatic regime is based on Case 6 in Giraldo and Restelli (2008). 
A bell-shape mountain with characteristic radius 10 km are located x=120 km. Initially, we set a uniform horizontal wind of 20 m/s. A periodic boundary condition is imposed at horizontal boundaries. Near model top and lateral boundaries, sponge layers are set to avoid the reflection of waves. 


In next section, we show results in the case of using the following model parameters 

- Spatial resolution: NeX=36, NeY=1, NeZ=5, p=7 (LGL nodes)
- Element-wise 16th-order exponential filter is applied for all modes and the factor of filter strength, $\alpha$, is 0.05; for example, the highest mode is dumped by $\exp{(-\alpha)}$.  
- $\Delta t=1.6$ [sec] (with a full explicit ten-stage and fourth-order RK scheme)


### 3. Result

<div class="container">
  <div class="item">
    This animation shows a result of mountain wave in linear and hydrostatic regime. The upper and lower panels represent vertical section (at $y=0$ m) of horizontal and vertical winds, respectively.
  </div>
  <div class="item">  
    <div class="youtube">    
      <iframe  width="448" height="600" src="https://www.youtube.com/embed/{{ site.data.gallery.mountain_wave_linear_hydrostatic_movie_id }}?rel=0" frameborder="0" allowfullscreen></iframe>
    </div>
  </div>
</div>

### 4. Reference

- Giraldo, F. X., and Restelli, M. 2008: A study of spectral element and discontinuous Galerkin methods for the Navier–Stokes equations in nonhydrostatic mesoscale atmospheric modeling: Equation sets and test cases. J. Comput. Phys., 227(8), 3849-3877.

- Hesthaven, J. S., and T. Warburton, 2007: Nodal discontinuous Galerkin methods: algorithms, analysis, and applications, Springer Science & Business Media
