---
layout: page
title: Gallery (under construction)
---

<div class="container">
  <div class="item">
    This page summarizes links to pages where we demonstrate results of various numerical experiments using our DG models. In each page, the experiment setup and model configuration are also provided in addition to the simulation results.
  </div>
  <div class="item">  
    <div class="youtube">
      <iframe  width="448" height="252" src="https://www.youtube.com/embed/{{ site.data.gallery.barocwavetest_movie_id }}?rel=0" frameborder="0" allowfullscreen></iframe>
    </div>
  </div>
</div>

### Simple examples and source codes

* Linear advection equation
  *  <a href="{{ '/gallery/simple/linear_advection_1d/' | relative_url }}">1D</a>, <a href="{{ '/gallery/simple/linear_advection_2d/' | relative_url }}">2D plane</a>, 3D plane, 2D global, 3D global
* Linear advection-diffusion equation
  * <a href="{{ '/gallery/simple/linear_adv_diffusion_1d/' | relative_url }}">1D</a>
* Euler equation

### Simulations with regional domains using SCALE-DG

* Tracer advection
* Sound wave
* Inertia gravity wave
* <a href="{{ '/gallery/regional/rising_warm_bubble/' | relative_url }}">Rising warm bubble</a>
* Density current
* Kelvin Helmholtz wave
* Mountain wave
* <a href="{{ '/gallery/regional/planetary_boundary_layer_tb/' | relative_url }}">Planetary boundary layer turbulence</a>
* <a href="{{ '/gallery/regional/baroclinic_instability/' | relative_url }}">Baroclinic wave</a>

### Simulations with global domains using SCALE-DG

* Tracer advection
* Sound wave
* Inertia gravity wave
* Mountain wave
* Equatorial wave
* Baroclinic wave
* Held Suarez test

### Simulations using global shallow water DG model

* Standard test cases proposed by Williamson (1992)
  * case 1: advection of cosine bell
  * case 2: steady state of nonlinear zonal geostrophic flow
  * case 5: zonal flow over an isolated mountain
  * case 6: Rossby-Haurwitz wave
* Barotropic instability
* Cross polar flow
