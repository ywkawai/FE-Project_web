---
layout: page
title: Gallery / Simulation with regional domain using SCALE-DG (This page is under construction)
---

{% assign files = "test_advect2d.f90,advect2d.conf,Makefile,visualize.py" | split: "," %}

## Rising warm bubble

This page shows a result of large eddy simulation of idealized planetary boundary layer turbulence using the discontinuous Galerkin method (DGM). 
An example of the setting files is given in rootdir/model/atmosphere3d/test/case/rising_thermal_bubble/.

### 1. Description of dynamical core

Please see XX.

### 2. Experimental setup

We follow an experimental setup for Robert smooth bubble described in section 3.2 in Giraldo and Restelli (2008), although the numerical experiment is originally presented in Robert (1993). 

In next section, we show results in the case of using the following model parameters 

- Spatial resolution: NeX=20, NeZ=30, p=10 (LGL nodes)
- Element-wise 32th-order exponential filter is applied for all modes and the factor of filter strength, $\alpha$, is 12.0; for example, the highest mode is dumped by $\exp{(-\alpha)}$.  
- $\Delta t=0.001$ [sec] (with a full explicit SSP third-order RK scheme)


### 3. Result

<div class="container">
  <div class="item">
    In this animation, the tone and contour represent the disturbance field of potential temperature, and the vectors represent wind fields, respectively. 
  </div>
  <div class="item">  
    <div class="youtube">    
      <iframe  width="340" height="605" src="https://www.youtube.com/embed/{{ site.data.gallery.rising_warm_bubble_movie_id }}?rel=0" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>
    </div>
  </div>
</div>

### 4. Reference

- Giraldo, F. X., & Restelli, M. (2008): A study of spectral element and discontinuous Galerkin methods for the Navier–Stokes equations in nonhydrostatic mesoscale atmospheric modeling: Equation sets and test cases. Journal of Computational Physics, 227(8), 3849-3877.

- Hesthaven, J. S., and T. Warburton, 2007: Nodal discontinuous Galerkin methods: algorithms, analysis, and applications, Springer Science & Business Media

- Robert, A. (1993): Bubble convection experiments with a semi-implicit formulation of the Euler equations, Journal of the Atmospheric
Sciences, 50, 1865–1873. 