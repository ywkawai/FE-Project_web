---
layout: page
title: Gallery / Simple example
---

{% assign files = "test_advect1d.f90,advect1d.conf,Makefile,visualize.py" | split: "," %}

## One-dimensional Linear Advection Problem

This page show an example 
to numerically solve a one-dimensional linear advection equation 
with the discontinuous Galerkin method (DGM) using FElib. 
The equation solved is as 
\\[ \dfrac{\partial u}{\partial t} + \dfrac{\partial (c u)}{\partial x} = 0  \;\;\; {\rm in} \;\; \Omega \times T \tag{1} \\]
where
u=u(x,t) is the scalar quantity, $c$ is a constant velocity, 
$t \in (0,T]$ for $T < \infty$, and $\Omega$ is $0 \le x \le L$. 
The initial condition is
\\[ u(x,0) = u_0(x) \;\; x \in \Omega \\]
and periodic boundary condition is imposed at $x=0, L$. 
The exact solution is given as
\\[ u(x,t) = u_0(x-ct). \tag{2} \\]
In this page, we assume that $L=1$, $c=1$, and $u_0(x)=\sin(4 \pi x/L)$.

### 1. Discretization

$\Omega$ is divided by non-overlapping finite element $\Omega_e$ 
and in each element we introduce a local coordinate as $\xi=2(x-x_e)/h_e$
where $x_e$ and $h_e$ are the center position and the width of the element.

In a nodal DGM, 
the numerical solution within $\Omega_e$ is represented by
\\[ u(\xi,t)|\_{\Omega_e} \sim u^e(\xi,t)=\sum_{j=0}^p U^e_j(\xi,t) l_j(\xi) \tag{3} \\] 
where $l_j(\xi)$ is the Lagrange polynomials associated with the j's Lagrange-Gauss-Lobatto (LGL) node. 
When we minimize L2 norm of a residual $\partial u^e/\partial t + \partial (c u^e)/\partial x$  against the tendency of each $U^e_j(\xi,t)$, 
the following equation is obtained
\\[ \dfrac{h_e}{2} \int_{-1}^1 \dfrac{\partial u^e}{\partial t} l_j d\xi +  \int_{-1}^1 \dfrac{\partial c u^e}{\partial \xi} l_j d\xi  = 0. \tag{4}  \\]
After performing integration by parts to the second term in the left of Eq. (4) twice, 
we have the strong form of the semi-discretized equation as
\\[ \dfrac{h_e}{2} \int_{-1}^1 \dfrac{\partial u^e}{\partial t} l_j d\xi +  \int_{-1}^1 \dfrac{\partial c u^e}{\partial \xi} l_j d\xi  + [ (\widehat{cu} - cq^e)]^{1}_{-1} = 0   \;\;\;\; (j=0,\cdots,p) \tag{5} \\]
where
$\widehat{Uq}$ is the numerical flux a a element boundary 
and the form is defined as
\\[ \widehat{cu}(u^l,u^r) = \dfrac{c}{2}(u^l + u^r) - \beta \dfrac{|c|}{2}(u^l - u^r) \tag{6} \\]
where $u^l$ and $u^r$ are the value of $u$ at left and the right sides of the element boundary, respectively.
The numerical flux for $\beta=0$ is the central flux,  while for $\beta=1$ is the full upwind flux. 
Equivalently, if we introduce the component of normal vector $n_x$ whose sign is defined such that the outward flux is positive, 
Eq. (6) can be written as 
\\[ \widehat{cu}(u^+,u^-) n_x  = \dfrac{c}{2}(u^+ + u^-)n_x - \beta \dfrac{|c|}{2}(u^+ - u^-) \tag{7} \\]
where $u^+$ and $u^-$ are the value of $u$ at the boundary of neighbor and own elements, respectively. 
Substituting Eq. (3) to Eq. (5) gives the semi-discretized equation based on the nodal DGM 
in the matrix-vector form as
\\[ \dfrac{h_e}{2} M \dfrac{d \vec{U^e}}{dt} = - S^T (c \vec{U^e}) - B (\widehat{cu} - c\vec{U^e})\\]
where
$\vec{U^e}=(U^e_0,\cdots,U^e_p)$, 
$M$ is the mass matrix, $S$ is the stiffness matrix, and $B={\rm diag}(-1,\cdots,1)$ is the matrix associated with boundary integration operator. 
If we multiply the both sides of the above equation by $2/h_e M^{-1}$, 
the differential form is obtained as follows
\\[ \dfrac{d \vec{U^e}}{dt} = - e D (c \vec{U^e}) -  f L (\widehat{cu} - c\vec{U^e}) \tag{8} \\]
where
$D=M^{-1} S^T$ is the differential matrix, 
$L=M^{-1} B$ is the lifting matrix, and $e=f=2/h_e$.  

The temporal discretization for the spatially semi-discretized equation in Eq. (8) is based on the method of lines 
which treats ordinary differential equations $d \vec{U^e}/dt = F(\vec{U^e})$.
For the full explicit DGM, 
the temporal discretization is often performed by strong Stability Preserving Runge-Kutta (SSP RK) 
and here we adopt a three stage and third order schemes (Shu 1988) defined as 
\\[ \vec{U^e}^{(1)} = \vec{U^e}^n + \Delta t  F(\vec{U^e}),  \\] 
\\[ \vec{U^e}^{(2)} = \dfrac{3}{4} \vec{U^e}^n + \dfrac{1}{4} \vec{U^e}^{(1)} + \dfrac{1}{4} \Delta t  F(\vec{U^e}^{(1)}),  \\] 
\\[ \vec{U^e}^{n+1} = \dfrac{1}{3} \vec{U^e}^n + \dfrac{2}{3} \vec{U^e}^{(2)} + \dfrac{2}{3} \Delta t  F(\vec{U^e}^{(2)})  \\]
where $\Delta t$ is the time step. 
Note that, for $p+1$ th-order and $p+1$ stage Runge-Kutta DG method with $p$th-order polynomials, 
$\Delta t$ should be less than $h_e/[c(2p+1)]$ to ensure the numerical stability.


### 2. Implementation

An implementation of the discretized equation described in the above section is shown below.
This source code is essentially the same as that in rootdir/sample/advect1d/test_advect1d.f90. 

{::options parse_block_html="true" /}

<details><summary markdown="span">Source code: test_advect1d.f90</summary>
```Fortran
#include "scalelib.h"
program linear_adv_eq
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
  real(RP), parameter :: c = 1.0_RP 
  real(RP) :: beta = 1.0_RP

  type(LineElement)  :: refElem
  type(SparseMat) :: Dx, Lift
  type(MeshLineDom1D), target :: mesh
  type(LocalMesh1D), pointer :: lcmesh

  type(MeshField1D), target :: u

  type(MeshFieldComm1D) :: fields_comm
  type(MeshFieldContainer) :: field_list(1)

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
    
      !* Exchange halo data
      call fields_comm%Put(field_list, 1)
      call fields_comm%Exchange()
      call fields_comm%Get(field_list, 1)

      !* Update prognostic variables
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        tintbuf_ind = tinteg_lc(n)%tend_buf_indmap(rkstage)      

        call cal_tend( tinteg_lc(n)%tend_buf2D_ex(:,:,RKVAR_U,tintbuf_ind), &
          u%local(n)%val, lcmesh, lcmesh%refElem1D                          )
 
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

contains
  subroutine cal_tend( dudt, u_, lmesh, elem )
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: dudt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
    !------------------------------------------------------------------------

    call cal_del_flux_dyn( del_flux,           & ! (out)
      u_, lmesh%normal_fn(:,:,1),              & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem )    ! (in)

    do ke = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, c * u_(:,ke), Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke), LiftDelFlx)
      dudt(:,ke) = - ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx )
    end do

    return
  end subroutine cal_tend

  subroutine cal_del_flux_dyn( del_flux, u_, nx, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA) 
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    !------------------------------------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      del_flux(i) = 0.5_RP * ( c * ( u_(iP) - u_(iM) ) * nx(i)     &
                             - beta * abs(c) * ( u_(iP) - u_(iM) ) )
    end do

    return
  end subroutine cal_del_flux_dyn

  subroutine set_initcond()
    use scale_const, only: PI => CONST_PI
    implicit none
    !-------------------------------------------------

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do ke=lcmesh%NeS, lcmesh%NeE
        u%local(n)%val(:,ke) = sin( 4.0_RP * PI * lcmesh%pos_en(:,ke,1) / ( dom_xmax - dom_xmin ) )
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

    namelist /PARAM_ADVECT1D/ &
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
    
    call IO_setup( "advect1d.conf", allow_noconf = .true. )
    call IO_LOG_setup( myrank, ismaster )  
    
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ADVECT1D,iostat=ierr)
    if( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_TEST. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ADVECT1D)

    
    call CALENDAR_setup
    call TIME_manager_Init

    !--

    call refElem%Init( PolyOrder, .false. )
    call Dx%Init( refElem%Dx1 )
    call Lift%Init( refElem%Lift )

    call mesh%Init( NeGX, dom_xmin, dom_xmax, refElem, 1 )
    call mesh%Generate()

    call u%Init( "u", "1", mesh )
    call fields_comm%Init( 1, 0, mesh )
    field_list(1)%field1d => u

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
    
    call u%Final()
    call fields_comm%Final()
    call mesh%Final()
    call Dx%Final(); call Lift%Final()    
    call refElem%Final()

    call TIME_manager_Final()
    call PRC_MPIfinish()

    return
  end subroutine final

end program linear_adv_eq
```
</details>
<br/>

{::options parse_block_html="false" /}

### 3. Build and Execution

- Build FElib following [install document]({{ '/documents/install.html' | relative_url }}). 
   Note that the environment variables set here are also used in subsequent steps.  

- Make a new directory in rootdir/sample/, and copy the following four files: 
{% for ver in files %}
  * [{{ ver }}]({{'/gallery/simple/linear_advection_1d/' | relative_url}}/{{ver}}) 
{% endfor %}

- Execute `% make` to compile test_advect1d.f90. 
  After this step, a binary file (test_advect1d) will be obtained.

- Execute `% make run` to run the program with advect1d.conf. 
  If this step succeeds, history.pe000000.nc will be output.

- Execute `% make vis` to visualize the simulation result. 
  By this command, the python script will generate advect1D.mp4.


Simulation paramters in advect1d.conf are configured as folllows:

- As for the spatial discretization, the number of element is 16 and the polynominal order ($p$) is 3.

- The integration period is $T=5$ and the timestep is $\Delta t=0.005$. 

### 4. Result

<div class="container">
  <div class="item">
   Left animation shows the simulation result for the parameters mentioned above.
   The numerical solution accurately maintins the sine curve of the inital condition. 
   Initial profile are simply moved with about the advected speed (U=1)  
   and the change in the amplitude of wave is quite small. 
  </div>
  <div class="item">  
    <div class="youtube">
      <iframe  width="448" height="252" src="https://www.youtube.com/embed/{{ site.data.gallery.1d_linear_advection_movie_id }}?rel=0" frameborder="0" allowfullscreen></iframe>
    </div>
  </div>
</div>

Although the details of numerical behavior associated with high-order DGM are not described here (see excercise described below), 
we can verify the $L^2$ error norm for $p=3$ and the linear advection equation with $C^\infty$ initial profile 
decreases with $O(h_e^4)$ (i.e., the fourth-order accuracy) 
when the element size decreases. 
This convergence rate is mainly related to the interpolation.

### 5. Exercise

- Change initial condition

  In this example, as an initial condition, 
  we set a sine wave with $C^\infty$. 
  Thus, when the polynomial order increases or the element size decreases, 
  the numerical solution rapidly converges.
  If we give a cosine bell or top-hat function which cannot be infinitely differentiated, 
  the numerical behaviors observed in the advected profile will be dramatically changed. 
  Confirm these behaviros by modifiying lines of source code in which the initial condition is set. 

- Investigate convergence rate of numerical solution

  Because the exact solution of the linear advection problem can be obtained at any time by Eq. (2), 
  we can easily investigate the numerical errors with the discretization. 
  For the detail of the source code, please see  rootdir/sample/advection1d/test_advect1d.f90. 

- Confirm effect of disspation with full upwinding numerical flux

  In section 4, 
  we show the results for the case of the full upwinding numerical flux ($\beta=1$ in Eq. (6)). 
  This numerical flux appropriately provides numerical disspation in the range of short wavelengths.
  On the other hand, if we set $\beta=0$ in Eq.

- Change temporal scheme

  FElib provides various type of SSP RK schemes. 
  Although a three-stage and third-order SSP RK scheme was applied in this example, 
  other scheme is also available by setting the parameter TINTEG_SCHEME_TYPE in advect1d.conf.   
  When the temporal scheme is replaced, confirm how the numerical stability change or the temporal error decreases.

### 6. Reference

- Hesthaven, J. S., and T. Warburton, 2007: Nodal discontinuous Galerkin methods: algorithms, analysis, and applications, Springer Science & Business Media
