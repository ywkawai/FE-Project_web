!-------------------------------------------------------------------------------
!> FE-Project example program
!!
!! @par Description
!!          one-dimensional linear advection problem
!!
!! @author Yuta Kawai
!!
!<
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
  type(MeshFieldContainer), save :: field_list(1)

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