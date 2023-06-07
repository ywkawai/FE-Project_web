---
layout: page
title: Gallery / Simulation with regional domain using SCALE-DG (This page is under construction)
---

{% assign files = "test_advect2d.f90,advect2d.conf,Makefile,visualize.py" | split: "," %}

## Baroclinic wave

This page shows a simulation result of idealized numerical experiment of baroclinic wave 
using the discontinuous Galerkin method (DGM). 
An example of the setting files is given in rootdir/model/atmosphere3d/test/case/baroclinic_wave/.

### 1. Description of dynamical core

Please see XX.

### 2. Experimental setup

The experiment setup is based on Ullrich et al. (2015). 
We consider a three-dimensional channel domain defined as $D = \{(x,y,z)| 0 ≤ x ≤ L_x,0 ≤ y ≤ L_y,0 ≤ z ≤ z_T \}$, where $L_x = 40000$ km, $L_y = 6000$ km, and $z_T = 30$ km. 
For the meridional, top, and bottom boundaries, no-flux boundary condition is applied. The zonal boundaries are periodic. 
As for rotation effect, we assume the beta-plane approximation. 

The basic state is a steady-state geostrophically balanced flow, and the analytic expressions are derived by Ullrich et al. (2015). 
To trigger the baroclinic instability, 
we add a perturbation into the basic state of zonal flow. 

In next section, we show results in the case of using the following model parameters 

- Spatial resolution: NeX=80, NeY=12, NeZ=12, p=7 (LGL nodes)
- Explicit diffusion with 4th-order differential operator with the decay coefficient 4x10^15[ m4/s] 
- Element-wise 16th-order exponential filter is applied for all modes and the factor of filter strength, $\alpha$, is 1.0; for example, the highest mode is dumped by $\exp{(-\alpha)}$.  
- $\Delta t=30$ [sec] (with HEVI method and temporal scheme of ARK324 )


### 3. Result

<div class="container">
  <div class="item">
    In this animation, upper panel shows the horizontal distributions of surface pressure (tone) and surface temperature (contour). The rectangle region bounded by red line is extended in lower panel. 
  </div>
  <div class="item">  
    <div class="youtube">
      <iframe  width="448" height="252" src="https://www.youtube.com/embed/{{ site.data.gallery.barocwavetest_movie_id }}?rel=0" frameborder="0" allowfullscreen></iframe>
    </div>
  </div>
</div>

### 4. Reference

- Hesthaven, J. S., and T. Warburton, 2007: Nodal discontinuous Galerkin methods: algorithms, analysis, and applications, Springer Science & Business Media

- Ullrich, P. A., K. A. Reed, and C. Jablonowski. "Analytical initial conditions and an analysis of baroclinic instability waves in f‐and β‐plane 3D channel models." Quarterly Journal of the Royal Meteorological Society 141.693 (2015): 2972-2988.