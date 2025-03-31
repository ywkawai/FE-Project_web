---
layout: page
title: Gallery / Simple example
---

{% assign files = "test_advdiff1d.f90,advdiff1d.conf,Makefile,visualize.py" | split: "," %}

## One-dimensional Linear Advection-Diffusion Problem

This page show an example 
to numerically solve a one-dimensional linear advection equation 
with the discontinuous Galerkin method (DGM) using FElib. 
The equation solved is as 
\\[ \dfrac{\partial u}{\partial t} + \dfrac{\partial (c u)}{\partial x} = \dfrac{\partial^2 (D u)}{\partial x^2}   \;\;\; {\rm in} \;\; \Omega \times T \tag{1} \\]
where
u=u(x,t) is the scalar quantity, $c$ is a constant velocity, $D$ is a diffusivity, and 
$t \in (0,T]$ for $T < \infty$, and $\Omega$ is $0 \le x \le L$. 
We consider an initial condition as
\\[ u(x,0) = \sin(k_x x) \;\;\; x \in \Omega \\]
and periodic boundary condition is imposed at $x=0, L$. 
The exact solution is given as
\\[ u(x,t) = \exp\left(-k_x^2 D t\right) \sin\left[ k_x (x-ct) \right]. \tag{2} \\]
In this page, we set $L=1$, $k_x=2\pi$, $c=0.5$, and $D=0.05$.

### 1. Discretization

$\Omega$ is divided by non-overlapping finite element $\Omega_e$ 
and in each element we introduce a local coordinate as $\xi=2(x-x_e)/h_e$
where $x_e$ and $h_e$ are the center position and the width of the element.

The numerical solution within $\Omega_e$ is represented by
\\[ u(\xi,t)|\_{\Omega_e} \sim u^e(\xi,t)=\sum_{j=0}^p U^e_j(\xi,t) l_j(\xi) \tag{3} \\] 
where $l_j(\xi)$ is the Lagrange polynomials associated with the j's Lagrange-Gauss-Lobatto (LGL) node. 
To treat the diffusion term, Eq.(1) is rewritten as
\\[ \dfrac{\partial u}{\partial t} + \dfrac{\partial (c u)}{\partial x} = \dfrac{\partial (D q)}{\partial x}, \;\;\; q = \dfrac{\partial u}{\partial x} \\]
where $q$ is an auxiliary variable represented by the polynomial expansion same as Eq. (3). 
Applying the Galerkin projection into the equations, 
the following equations are obtained
\\[ \dfrac{h_e}{2} \int_{-1}^1 \dfrac{\partial u}{\partial t} l_j d\xi +  \int_{-1}^1 \dfrac{\partial (c u - D q)}{\partial \xi} l_j d\xi  = 0, \;\;\;  \dfrac{h_e}{2} \int_{-1}^1 q l_j d\xi = \int_{-1}^1 \dfrac{\partial u}{\partial \xi} l_j d\xi \tag{4}  \\]
After performing integration by parts to the terms with the differential twice, 
we have the strong form of the semi-discretized equation as
\\[ \dfrac{h_e}{2} \int_{-1}^1 \dfrac{\partial u^e}{\partial t} l_j d\xi +  \int_{-1}^1 \dfrac{\partial (c u^e - q^e)}{\partial \xi} l_j d\xi  + [ (\widehat{cu} + \widehat{q}) - (cu^e - q^e)]^{1}_{-1} = 0  \tag{5a} \\]
and 
\\[ \dfrac{h_e}{2} \int^1\_{-1}  q^e l_j d\xi = \int^1\_{-1} \dfrac{\partial u^e}{\partial \xi} l_j d\xi + [ \widehat{u} - u^e ]^{1}\_{-1}  \tag{5b} \\]
where
$\widehat{*}$ indicates a numerical fluxes at an element boundary. 
For the advection term, we use the full upwind flux as 
\\[ \widehat{cu}(u^l,u^r) = \dfrac{c}{2}(u^l + u^r) - \dfrac{|c|}{2}(u^l - u^r) \tag{6} \\] 
where $u^l$ and $u^r$ are the values of $u$ at left and the right sides of the element boundary. 
For the diffusive term, we ues the alternating flux used in the local DGM as
\\[ \widehat{q} = q^l, \;\;\; \widehat{u} = u^r \\]
Equivalently, if we introduce the component of normal vector $n_x$ whose sign is defined such that the outward flux is positive, 
Eq. (6) can be written as 
\\[ \widehat{cu}(u^+,u^-) n_x  = \dfrac{c}{2}(u^+ + u^-)n_x - \beta \dfrac{|c|}{2}(u^+ - u^-) \tag{7} \\]
where $u^+$ and $u^-$ are the value of $u$ at the boundary of neighbor and own elements, respectively. 
Substituting Eq. (3) to Eq. (5) gives the semi-discretized equation based on the nodal DGM 
in the matrix-vector form as
\\[ \dfrac{h_e}{2} M \dfrac{d \vec{U^e}}{dt} = - S^T (c \vec{U^e} - D \vec{Q}^e) - B (\widehat{cu} + D\widehat{q} - c\vec{U^e} - D \vec{Q}^e), \\]
\\[ \dfrac{h_e}{2} M \vec{Q}^e = S^T \vec{U}^e + B [\widehat{u} - \vec{U}^e] \\]
where
$\vec{U^e}=(U^e_0,\cdots,U^e_p)$, $\vec{Q^e}=(Q^e_0,\cdots,Q^e_p)$, 
$M$ is the mass matrix, $S$ is the stiffness matrix, and $B={\rm diag}(-1,\cdots,1)$ is the matrix associated with boundary integration operator. 
If we multiply the both sides of the above equation by $2/h_e M^{-1}$, 
the differential form is obtained as follows
\\[ \dfrac{d \vec{U^e}}{dt} = - e D (c \vec{U^e} - D \vec{Q}^e) -  f L (\widehat{cu} + D\widehat{q} - c\vec{U}^e - D \vec{Q}^e), \;\;\;  \vec{Q}^e = e D \vec{U}^e + f L(\widehat{u} - \vec{U}^e) \tag{8} \\]
where
$D=M^{-1} S^T$ is the differential matrix, 
$L=M^{-1} B$ is the lifting matrix, and $e=f=2/h_e$.  

The temporal discretization for the spatially semi-discretized equation in Eq. (8) is based on the method of lines 
which treats ordinary differential equations $d \vec{U^e}/dt = F(\vec{U^e})$.
For the full explicit DGM, 
the temporal discretization is often performed by strong Stability Preserving Runge-Kutta (SSP RK) 
and here we adopt a four stage and third order schemes.


### 2. Implementation

An implementation of the discretized equation described in the above section is shown below.
This source code is essentially the same as that in rootdir/sample/advdiff1d/test_advdiff1d.f90. 

{::options parse_block_html="true" /}

<details><summary markdown="span">Source code: test_advdiff1d.f90</summary>
```Fortran
#include "scaleFElib.h"
program linear_advdiff_eq
  use scale_precision
  use scale_io
  use scale_prc
  use scale_sparsemat
  use scale_element_base
  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d
  use scale_meshfield_base
  use scale_meshfieldcomm_base
  use scale_meshfieldcomm_1d
  use scale_file_history_meshfield
  use scale_file_history 
  use scale_time_manager  
  use scale_timeint_rk
  implicit none
  !-----------------------------------------------------------------------------

  real(RP), parameter :: dom_xmin = 0.0_RP, dom_xmax = 1.0_RP
  real(RP), parameter :: c = 0.5_RP 
  real(RP), parameter :: D = 0.05_RP 
  real(RP) :: beta = 1.0_RP

  type(LineElement)  :: refElem
  type(SparseMat) :: Dx, Lift
  type(MeshLineDom1D), target :: mesh
  type(LocalMesh1D), pointer :: lcmesh

  type(MeshField1D), target :: u, q

  type(MeshFieldComm1D) :: fields_comm
  type(MeshFieldComm1D) :: auxvars_comm
  type(MeshFieldContainer), save :: field_list(1)
  type(MeshFieldContainer), save :: auxvars_list(1)

  type(timeint_rk), allocatable :: tinteg_lc(:)
  integer :: nowstep
  integer :: rkstage
  integer :: tintbuf_ind
  integer, parameter :: RKVAR_U = 1

  integer :: n
  integer :: ke
  integer :: HST_ID
  !-----------------------------------------------------------------------------

  call init()
  call set_initcond()

  do nowstep=1, TIME_NSTEP

    do rkstage=1, tinteg_lc(1)%nstage
    
      !* Exchange halo data (prognostic variables)
      call fields_comm%Put(field_list, 1)
      call fields_comm%Exchange()
      call fields_comm%Get(field_list, 1)

      !* Update auxiliary variables
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        call cal_aux( q%local(n)%val, &
          u%local(n)%val, lcmesh, lcmesh%refElem1D )
      end do

      !* Exchange halo data (auxiliary variables)
      call auxvars_comm%Put(auxvars_list, 1)
      call auxvars_comm%Exchange()
      call auxvars_comm%Get(auxvars_list, 1)

      !* Update prognostic variables
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        tintbuf_ind = tinteg_lc(n)%tend_buf_indmap(rkstage)      

        call cal_prg_tend( tinteg_lc(n)%tend_buf2D_ex(:,:,RKVAR_U,tintbuf_ind), &
          u%local(n)%val, q%local(n)%val, lcmesh, lcmesh%refElem1D              )
 
        call tinteg_lc(n)%Advance( rkstage, u%local(n)%val, RKVAR_U,            &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
      end do
    end do

    !* Advance time
    call TIME_manager_advance()    
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWSUBSEC, TIME_NOWSTEP )

    !* Output
    call FILE_HISTORY_meshfield_put(HST_ID, u)
    call FILE_HISTORY_meshfield_write()
  end do

  call final()
  LOG_INFO("check",*) "-----------", nowstep

contains
subroutine cal_prg_tend( dudt, u_, q_, lmesh, elem)
    use scale_sparsemat, only: sparsemat_matmul
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: dudt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)

    integer :: ke
    !------------------------------------------------------------------------

    call cal_del_flux_prg( del_flux,           & ! (out)
      u_, q_, lmesh%normal_fn(:,:,1),          & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem )    ! (in)

    !$omp parallel do private(Fx, LiftDelFlx)
    do ke = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, c * u_(:,ke) - D * q_(:,ke), Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke), LiftDelFlx)

      dudt(:,ke) = - (  lmesh%Escale(:,ke,1,1) * Fx(:) &
                     + LiftDelFlx(:) )
    end do

    return
  end subroutine cal_prg_tend

  subroutine cal_del_flux_prg( del_flux, u_, q_, nx, vmapM, vmapP, lmesh, elem )    
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: alpha
    !------------------------------------------------------------------------

    !$omp parallel do private(iM, iP, alpha)
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      alpha = 0.5_RP * abs( c )
      del_flux(i) = 0.5_RP * (                    &
          ( c * u_(iP) - c * u_(iM) ) * nx(i)     &
        - alpha * ( u_(iP) - u_(iM) )             & 
        - D * ( 1.0_RP - nx(i) )                  &
           * ( q_(iP) - q_(iM) ) * nx(i)          )
    end do

    return
  end subroutine cal_del_flux_prg

  subroutine cal_aux( q_, u_, lmesh, elem)
    use scale_sparsemat, only: sparsemat_matmul
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: u_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)

    integer :: ke
    !------------------------------------------------------------------------

    call cal_bnd_flux_aux( del_flux,           & ! (out)
      u_, lmesh%normal_fn(:,:,1),              & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem )    ! (in)

    !-----
    !$omp parallel do private(Fx, LiftDelFlx)
    do ke = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, u_(:,ke), Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke), LiftDelFlx)

      q_(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) &
               + LiftDelFlx(:)
    end do

    return
  end subroutine cal_aux

  subroutine cal_bnd_flux_aux( del_flux, u_, nx, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: delVar
    !------------------------------------------------------------------------

    !$omp parallel do private(i, iM, iP, delVar)
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)
      delVar = 0.5_RP * ( u_(iP) - u_(iM) )
      del_flux(i) = ( 1.0_RP + nx(i) ) * delVar * nx(i)
    end do

    return
  end subroutine cal_bnd_flux_aux  

  subroutine set_initcond()
    use scale_const, only: PI => CONST_PI
    implicit none
    !-------------------------------------------------

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do ke=lcmesh%NeS, lcmesh%NeE
        u%local(n)%val(:,ke) = sin( 2.0_RP * PI * lcmesh%pos_en(:,ke,1) / ( dom_xmax - dom_xmin ) )
      end do
    end do

    call FILE_HISTORY_meshfield_put(HST_ID, u)
    call FILE_HISTORY_meshfield_write()   

    return
  end subroutine set_initcond

  subroutine init()
    use scale_calendar
    implicit none  

    integer :: NeGX      = 8
    integer :: PolyOrder = 4
    character(len=H_MID) :: TINTEG_SCHEME_TYPE  = 'ERK_SSP_3s3o'

    namelist /PARAM_ADVDIFF1D/ &
      NeGX, PolyOrder,        &
      TINTEG_SCHEME_TYPE,     &
      beta
        
    integer :: comm, myrank, nprocs
    logical :: ismaster      
    integer :: ierr
    
    !---------------------------------

    call PRC_MPIstart( comm )
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    call IO_setup( "advdiff1d.conf", allow_noconf = .true. )
    call IO_LOG_setup( myrank, ismaster )  
    
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ADVDIFF1D,iostat=ierr)
    if( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_ADVDIFF1D. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ADVDIFF1D)

    
    call CALENDAR_setup
    call TIME_manager_Init

    !--

    call refElem%Init( PolyOrder, .false. )
    call Dx%Init( refElem%Dx1 )
    call Lift%Init( refElem%Lift )

    call mesh%Init( NeGX, dom_xmin, dom_xmax, refElem, 1 )
    call mesh%Generate()

    call u%Init( "u", "1", mesh )
    call q%Init( "q", "1", mesh )
    call fields_comm%Init( 1, 0, mesh )
    call auxvars_comm%Init( 1, 0, mesh )
    field_list(1)%field1d => u
    auxvars_list(1)%field1d => q

    call FILE_HISTORY_meshfield_setup( mesh )
    call FILE_HISTORY_reg( u%varname, "u", u%unit, HST_ID, dim_type='X')

    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call tinteg_lc(n)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1,      &
                              2, (/ lcmesh%refElem%Np, lcmesh%NeA /)  )
    end do

    return
  end subroutine init

  subroutine final()
    implicit none
    !---------------------------------

    call FILE_HISTORY_meshfield_finalize()
    
    call u%Final(); call q%Final()
    call fields_comm%Final()
    call auxvars_comm%Final()
    call mesh%Final()
    call Dx%Final(); call Lift%Final()    
    call refElem%Final()

    call TIME_manager_Final()
    call PRC_MPIfinish()

    return
  end subroutine final

end program linear_advdiff_eq
```
</details>
<br/>

{::options parse_block_html="false" /}

### 3. Build and Execution

- Build FElib following [install document]({{ '/documents/install.html' | relative_url }}). 
   Note that the environment variables set here are also used in subsequent steps.  

- Make a new directory in rootdir/sample/, and copy the following four files: 
{% for ver in files %}
  * [{{ ver }}]({{'/gallery/simple/linear_adv_diffusion_1d/' | relative_url}}/{{ver}}) 
{% endfor %}

- Execute `% make` to compile test_advdiff1d.f90. 
  After this step, a binary file (test_advdiff1d) will be obtained.

- Execute `% make run` to run the program with advdiff1d.conf. 
  If this step succeeds, history.pe000000.nc will be output.

- Execute `% make vis` to visualize the simulation result. 
  By this command, the python script will generate advdiff1D.mp4.


Simulation paramters in advdiff1d.conf are configured as folllows:

- As for the spatial discretization, the number of element is 16 and the polynominal order ($p$) is 3.

- The integration period is $T=3$ and the timestep is $\Delta t=0.0005$. 

### 4. Result

<div class="container">
  <div class="item">
   Left animation shows the simulation result for the parameters mentioned above.
   Initial profile are simply moved with the advected speed $c$ while 
   decreasing the amplitude with the e-folding time of $(Dk_x^2)^{-1} \sim 0.5$.
  </div>
  <div class="item">  
    <div class="youtube">
      <iframe  width="448" height="252" src="https://www.youtube.com/embed/{{ site.data.gallery.1d_linear_advdiff_movie_id }}?rel=0" frameborder="0" allowfullscreen></iframe>
    </div>
  </div>
</div>

Although the details of numerical behavior associated with high-order DGM are not described here (see excercise described below), 
we can verify that the $L^2$ error norm for $p=3$ and the linear advection-diffusion equation with $C^\infty$ initial profile 
decreases with $O(h_e^4)$ (i.e., the fourth-order accuracy) 
when the element size decreases. 

### 5. Exercise

- Change initial condition

- Investigate convergence rate of numerical solution

### 6. Reference

- Hesthaven, J. S., and T. Warburton, 2007: Nodal discontinuous Galerkin methods: algorithms, analysis, and applications, Springer Science & Business Media
