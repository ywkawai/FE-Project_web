---
layout: page
title: Gallery / Simulation with regional domain using SCALE-DG (This page is under construction)
---

{% assign files = "test_advect2d.f90,advect2d.conf,Makefile,visualize.py" | split: "," %}

## Planetary boundary layer turbulence

This page shows a result of large eddy simulation of idealized planetary boundary layer turbulence using the discontinuous Galerkin method (DGM). 
An example of the setting files is given in rootdir/model/atmosphere3d/test/case/pbl_turbulence/.

### 1. Description of dynamical core

Please see XX.

### 2. Experimental setup

The experiment setup is based on Nishizawa et al. (2015). 
To drive thermal convection, a constant heat flux with 200 W/m2 is injected at the surface.
To parametrize the sub-grid scale effect of turbulence, 
a Smagorinsky-Lilly type turbulent model (Smagorinsky, 1963; Lilly, 1962; Brown et al., 1994) is used. 

In next section, we show results in the case of using the following model parameters 

- Spatial resolution: NeX=48, NeY=48, NeZ=10, p=7 (LGL nodes)
- Element-wise 16th-order exponential filter is applied for all modes and the factor of filter strength, $\alpha$, is 0.05; for example, the highest mode is dumped by $\exp{(-\alpha)}$.  
- $\Delta t=0.003$ [sec] (with a full explicit ten-stage and fourth-order RK scheme)


### 3. Result

<div class="container">
  <div class="item">
    In this animation, the horizontal cross section (at $z=500$ m) and the vertical section (at $y=0$ m) of vertical wind are presented in upper and lower panels, respectively. 
  </div>
  <div class="item">  
    <div class="youtube">    
      <iframe  width="448" height="600" src="https://www.youtube.com/embed/{{ site.data.gallery.pbl_turbulence_movie_id }}?rel=0" frameborder="0" allowfullscreen></iframe>
    </div>
  </div>
</div>

### 4. Reference

- Brown, A. R., S. H. Derbyshire, and P. J. Mason, 1994: Large-eddy simulation of stable atmo- spheric boundary layers with a revised stochastic subgrid model. Quarterly Journal of the Royal Meteorological Society, 120 (520), 1485–1512.

- Hesthaven, J. S., and T. Warburton, 2007: Nodal discontinuous Galerkin methods: algorithms, analysis, and applications, Springer Science & Business Media

- Nishizawa, S., H. Yashiro, Y. Sato, Y. Miyamoto, and H. Tomita, 2015: Influence of grid aspect ratio on planetary boundary layer turbulence in large-eddy simulations. Geoscientific Model Development, 8 (10), 3393–3419.

- Smagorinsky, J., 1963: GENERAL CIRCULATION EXPERIMENTS WITH THE PRIMITIVE EQUATIONS: I. THE BASIC EXPERIMENT. Monthly Weather Review, 91 (3), 99–164.

- Lilly, D. K., 1962: On the numerical simulation of buoyant convection. Tellus, 14 (2), 148–172.