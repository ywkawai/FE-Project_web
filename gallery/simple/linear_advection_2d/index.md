---
layout: page
title: Gallery / Simple example
---

{% assign files = "test_advect2d.f90,advect2d.conf,Makefile,visualize.py" | split: "," %}

## Two-dimensional Linear Advection Problem

This page shows an example 
to numerically solve a two-dimensional linear advection equation 
with the discontinuous Galerkin method (DGM) using FElib. 
The equation solved is as 
$$
\begin{align} 
\dfrac{\partial u}{\partial t} + \dfrac{\partial (c_x u)}{\partial x} +  \dfrac{\partial (c_y u)}{\partial x} = 0  \;\;\; {\rm in} \;\; \Omega \times T 
\end{align}
$$
where
u=u(x,y,t) is the scalar quantity,  
$t \in (0,T]$ for $T < \infty$, $\Omega$ is the computational domain with $0 \le x \le L_x$ and $0 \le y \le L_y$, and the components of advective velocity in x-direction and y-direction are denoted by $c_x$ and $c_y$, respectively. 
The initial condition is
\\[ u(x,y,0) = u_0(x,y) \;\; x \in \Omega \\]
and periodic boundary condition is imposed at $x=0, L_x$ for the x-direction and $y=0, L_y$ for the y-direction. 

In this page, we assume that $L_x=L_y=1$.
an exponential profile $u_0(x,y) = \exp{[-(x^2+y^2)/R^2}]$ where $R=1/10$ is used as a initial condition.
As for the advective velocity, we consider the following two cases:
i) a uniform flow with constant velocity with $c_x=c_y=1/\sqrt{2}$
ii) the swirling flow with $c_x=\sin^2{(\pi x)} \sin{(2\pi y)} \cos{(2\pi t)}$ 
and $c_y=-\sin^2{(\pi y)} \sin{(2\pi x)} \cos{(2\pi t)}$.
Note that, for the both cases, 
the profile at $t=n$ (where $n=1,2,\dots$) should be the same as the initial one if the equation would be exactly solved.

### 1. Discretization

The computational domain is divided by non-overlapping quadrilateral element $\Omega_e$ 
whose sides have a width of $h_e$.
In each element, we introduce a local coordinate as $\xi=2(x-x_e)/h_e$ and $\zeta=2(y-y_e)/h_e$.
where $x_e$ and $y_e$ are the center position in the x- and y- directions, respectively.
When the quadrilateral elements are used, 
we can represent the numerical solution within $\Omega_e$ with tensor-product manner as
$$
\begin{align}
 u(\xi,\zeta,t)|\_{\Omega_e} \sim u^e(\xi,t)=\sum_{j1=0}^p \sum_{j2=0}^p U^e_{j1,j2}(\xi,\zeta,t) l_{j1}(\xi) l_{j2}(\zeta)
\end{align}
$$ 
where $l_{j1}(\xi)$ and $l_{j2}(\zeta)$ are the Lagrange polynomials associated with the $j$'s and Lagrange-Gauss-Lobatto (LGL) node. 
In the tensor-product formulation, 
the discretization strategy for the one-dimensional linear advection problem can be straightforwardly extended to the multi-dimensional one.  
The corresponding strong form of the semi-discretized equation with nodal DGM is written as
$$
\begin{align}
 \dfrac{h_e^2}{4} \int_{-1}^1 \int_{-1}^1 \dfrac{\partial u^e}{\partial t} l_{j1} l_{j2} d\xi d\zeta &+  \int_{-1}^1 \int_{-1}^1 \left[ \dfrac{\partial c_x u^e}{\partial \xi} + \dfrac{\partial c_y u^e}{\partial \xi} \right] l_{j1} l_{j2}  d\xi d\eta  \nonumber \\
 &+ \left[ \int_{-1}^1 ( \widehat{c_x u} - c_x u^e ) d\eta \right]^{\xi=1}_{\xi=-1}  + \left[  \int_{-1}^1  ( \widehat{c_y u} - c_y u^e ) d\xi \right]^{\zeta=1}_{\zeta=-1} 
  = 0  
 \;\;\;\; (j=0,\cdots,p)
\end{align}
$$
where 
the numerical flux at element boundary in $x-$ or $y-$ direction is defined as 
$$
\begin{align} 
 \widehat{cu}(u^+,u^-) n  = \dfrac{c}{2}(u^+ + u^-)n - \beta \dfrac{|c|}{2}(u^+ - u^-)
\end{align}
$$
where 
$u^+$ and $u^-$ are the value of $u$ at the boundary of neighbor and own elements, respectively.
For the x- and y- directions, $n=n_x, c=c_x$  and $n=n_y, c=c_y$ for the y-direction 
(where $n_x$ and $n_y$ are the component of normal vector at element boundary), respectively. 

Substituting Eq. (2) to Eq. (3) gives the semi-discretized equation based on the nodal DGM 
in the matrix-vector form as
$$
\begin{align} 
 \dfrac{h_e^2}{4} M \dfrac{d \vec{U^e}}{dt} 
 =& -\left[ S_x^T (c_x \vec{U^e}) + S_y^T (c_y \vec{U^e}) \right] \nonumber \\
  &- B_x (\widehat{c_x u} - c_x \vec{U^e}) - B_y (\widehat{c_y u} - c_y \vec{U^e})
\end{align}
$$
where
$\vec{U^e}=(U^e_{0,0},\cdots,U^e_{p,p})$, 
$M$ is the mass matrix, 
$S_x$ and $S_y$ are the stiffness matrix in $x-$ and $y-$ direction, 
and $S_{x,1D}$ and $S_{y,1D}$  are the matrix associated with boundary integration in $x-$ and $y-$ direction. 
If we multiply the both sides of the above equation by $4/h_e^2 M^{-1}$, 
the differential form is obtained as follows
$$
\begin{align}
  \dfrac{d \vec{U^e}}{dt} = - e [ D_x (c_x \vec{U^e}) + D_y (c_y \vec{U^e}) ] -  f L (\widehat{cu} - c\vec{U^e})
\end{align}
$$
where
$D=M^{-1} S^T$ is the differential matrix, 
$L=M^{-1} B$ is the lifting matrix, and $e=f=2/h_e$.  

As in the one-dimensional linear advection problem, 
the temporal discretization for the spatially semi-discretized equation in Eq. (8) is based on the method of lines 
and we adopt a strong Stability Preserving Runge-Kutta (SSP RK) method with three stages and third order accuracy. 


### 2. Implementation

An implementation of the discretized equation described in the above section is shown below.
This source code is essentially the same as that in rootdir/sample/advect2d/test_advect2d.f90. 

{::options parse_block_html="true" /}

<details><summary markdown="span">Source code: test_advect2d.f90</summary>
```Fortran
#include "scalelib.h"
program linear_adv_eq
  use scale_precision
  use scale_io
  use scale_prc
  use scale_sparsemat
  use scale_element_base
  use scale_element_quadrilateral
  use scale_localmesh_2d
  use scale_mesh_rectdom2d
  use scale_meshfield_base
  use scale_meshfieldcomm_base
  use scale_meshfieldcomm_rectdom2d
  use scale_file_history_meshfield
  use scale_file_history 
  use scale_time_manager  
  use scale_timeint_rk
  implicit none
  !-----------------------------------------------------------------------------

  real(RP), parameter :: dom_xmin = 0.0_RP, dom_xmax = 1.0_RP
  real(RP), parameter :: dom_ymin = 0.0_RP, dom_ymax = 1.0_RP  
  real(RP), parameter :: c_x = 1.0_RP, c_y = 1.0_RP
  real(RP) :: beta = 1.0_RP

  type(QuadrilateralElement)  :: refElem
  type(SparseMat) :: Dx, Dy, Lift
  type(MeshRectDom2D), target :: mesh
  type(LocalMesh2D), pointer :: lcmesh

  type(MeshField2D), target :: u

  type(MeshFieldCommRectDom2D) :: fields_comm
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
          u%local(n)%val, lcmesh, lcmesh%refElem2D                          )
 
        call tinteg_lc(n)%Advance( rkstage, u%local(n)%val, RKVAR_U,              &
                                   1, lcmesh%refElem2D%Np, lcmesh%NeS, lcmesh%NeE )
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
  write(*,*) "Final"

contains
  subroutine cal_tend( dudt, u_, lmesh, elem )
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: dudt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
    !------------------------------------------------------------------------

    call cal_elembnd_flux( del_flux,                      & ! (out)
      u_, lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem               ) ! (in)

    do ke = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, c_x * u_(:,ke), Fx)
      call sparsemat_matmul(Dy, c_y * u_(:,ke), Fy)      
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke), LiftDelFlx)
      dudt(:,ke) = - ( lmesh%Escale(:,ke,1,1) * Fx(:) &
                     + lmesh%Escale(:,ke,2,2) * Fy(:) &
                     + LiftDelFlx(:) )
    end do

    return
  end subroutine cal_tend

  subroutine cal_elembnd_flux( del_flux, u_, nx, ny, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA) 
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)    
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: vel
    !------------------------------------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      vel = c_x * nx(i) + c_y * ny(i)
      del_flux(i) = 0.5_RP * ( vel * ( u_(iP) - u_(iM) )             &
                             - beta * abs(vel) * ( u_(iP) - u_(iM) ) )
    end do

    return
  end subroutine cal_elembnd_flux

  subroutine set_initcond()
    use scale_const, only: PI => CONST_PI
    implicit none
    !-------------------------------------------------

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do ke=lcmesh%NeS, lcmesh%NeE
        u%local(n)%val(:,ke) = sin( 4.0_RP * PI * lcmesh%pos_en(:,ke,1) / ( dom_xmax - dom_xmin ) ) &
                             * sin( 4.0_RP * PI * lcmesh%pos_en(:,ke,2) / ( dom_ymax - dom_ymin ) )
      end do
    end do

    call FILE_HISTORY_meshfield_put(HST_ID, u)
    call FILE_HISTORY_meshfield_write()   

    return
  end subroutine set_initcond

  subroutine init()
    use scale_calendar
    implicit none  

    integer :: NeGX      = 8, NeGY      = 8
    integer :: NprcX     = 1, NprcY     = 1
    integer :: PolyOrder = 4
    character(len=H_MID) :: TINTEG_SCHEME_TYPE  = 'ERK_SSP_3s3o'

    namelist /PARAM_ADVECT2D/ &
      NprcX, NeGX, NprcY, NeGY, &
      PolyOrder,                &
      TINTEG_SCHEME_TYPE,       &
      beta
        
    integer :: comm, myrank, nprocs
    logical :: ismaster      
    integer :: ierr
    !---------------------------------

    call PRC_MPIstart( comm )
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    call IO_setup( "advect2d.conf", allow_noconf = .true. )
    call IO_LOG_setup( myrank, ismaster )  
    
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ADVECT2D,iostat=ierr)
    if( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_TEST. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ADVECT2D)

    
    call CALENDAR_setup
    call TIME_manager_Init

    !--

    call refElem%Init( PolyOrder, .false. )
    call Dx%Init( refElem%Dx1 )
    call Dy%Init( refElem%Dx2 )    
    call Lift%Init( refElem%Lift )

    call mesh%Init( NeGX, NeGY, dom_xmin, dom_xmax, dom_ymin, dom_ymax, &
      .true., .true., refElem, 1, NprcX, NprcY )
    call mesh%Generate()

    call u%Init( "u", "1", mesh )
    call fields_comm%Init( 1, 0, mesh )
    field_list(1)%field2d => u

    call FILE_HISTORY_meshfield_setup( mesh2D_=mesh )
    call FILE_HISTORY_reg( u%varname, "u", u%unit, HST_ID, dim_type='XY')

    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call tinteg_lc(n)%Init( TINTEG_SCHEME_TYPE, TIME_DTSEC, 1,       &
                              2, (/ lcmesh%refElem2D%Np, lcmesh%NeA /) )
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
    call Dx%Final(); call Dy%Final(); call Lift%Final()    
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
  * [{{ ver }}]({{'/gallery/simple/linear_advection_2d/' | relative_url}}/{{ver}}) 
{% endfor %}

- Execute `% make` to compile test_advect2d.f90. 
  After this step, a binary file (test_advect1d) will be obtained.

- Execute `% make run` to run the program with advect2d.conf. 
  If this step succeeds, history.pe000000.nc will be output.

- Execute `% make vis` to visualize the simulation result. 
  By this command, the python script will generate advect2D.mp4.


Simulation paramters in advect2d.conf are configured as folllows:

- As for the spatial discretization, the number of element is $16\times 16$ and the polynominal order ($p$) is 3.

- The integration period is $T=5$ and the timestep is $\Delta t=0.005$. 

### 4. Result

<div class="container">
  <div class="item">
   Left animation shows the simulation result for the parameters mentioned above.
   The numerical solution accurately maintins the sine curve of the inital condition. 
   Initial profile are simply moved with about the advected speed ($c_x=c_y=1$)  
   and the change in the amplitude of wave is quite small. 
  </div>
  <div class="item">  
    <div class="youtube">
      <iframe  width="448" height="252" src="https://www.youtube.com/embed/{{ site.data.gallery.2d_linear_advection_movie_id }}?rel=0" frameborder="0" allowfullscreen></iframe>
    </div>
  </div>
</div>

Although the details of numerical behavior associated with high-order DGM are not described here (see excercise described below), 
we can verify the $L^2$ error norm for $p=3$ and the linear advection equation with $C^\infty$ initial profile 
decreases with $O(h_e^4)$ (i.e., the fourth-order accuracy) for significantly small $\Delta t$ 
when the element size decreases. 
This convergence rate is mainly related to the interpolation.

### 5. Exercise

- Change initial condition
- Investigate convergence rate of numerical solution
- Confirm effect of disspation with full upwinding numerical flux
- Change temporal scheme

### 6. Reference

- Hesthaven, J. S., and T. Warburton, 2007: Nodal discontinuous Galerkin methods: algorithms, analysis, and applications, Springer Science & Business Media
