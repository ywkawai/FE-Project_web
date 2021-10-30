---
layout: page
title: Install Document
---

## Install Depending Softwares

### Required Softwares

* Fortran compiler supporting Fortran 2008
* MPI library
* LAPACK
* NetCDF
* [SCALE library](https://scale.riken.jp)

### Optional Softwares

* visualization
  * [matplotlib](https://matplotlib.org)
  * [gpview (GPhys)](http://ruby.gfd-dennou.org/products/gphys/)

## Install SCALE
See [Users Guide of SCALE](https://scale.riken.jp/documents/index.html#users-guide).

## Build FE-Project

1. preparation

  - set SCALE_FE_SYS environmental variable (see the sysdef directory)

  `% export SCALE_FE_SYS=MacOSX-gnu-ompi`   (for example)

  - set a directory in which SCALE library is contained

  `% export SCALE="~/workspace/scale-5.4.4/"`   (for example)

  - If you use the develop version of SCALE library, set a variable as

  `% export SCALE_DEVELOP=T`

  - set a directory in which NetCDF library is contained (if necessary).

  `% export NETCDF="/ap/netcdf4-fortran/4.7.3/"`   (for example)

2. build the library in the directory of FElib

 `% cd rootdir/FElib/src/`

 `% make`

## Compile and run simple sample programs

 For example, in the case of sample/advect1d, 
 
 `% cd rootdir/sample/advect1d/`

 `% make`

## Compile and run atmospheric models

 If you want to build a three-dimensional nonhydrostatic atmospheric model, 
 and conduct a idealized test case, such as density current, using it, 
 
 `% cd rootdir/model/atm_nonhydro3d/test/case/density_current`

 `% make`

 `% make run`

 In the directory of 'visualize', some python scripts with matplotlib are prepared for visualizing simulation results.
