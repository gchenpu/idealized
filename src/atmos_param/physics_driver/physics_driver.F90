                           module physics_driver_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="">
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!     Provides high level interfaces for calling the entire
!     FMS atmospheric physics package.
!
!    physics_driver_mod accesses the model's physics modules and
!    obtains tendencies and boundary fluxes due to the physical
!    processes that drive atmospheric time tendencies and supply 
!    boundary forcing to the surface models.
! </OVERVIEW>
! <DESCRIPTION>
!     This version of physics_driver_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Radiation, Rayleigh damping, gravity wave drag, vertical diffusion of
!     momentum and tracers, and the downward pass of vertical diffusion for
!     temperature and specific humidity are performed in the down routine.
!     The up routine finishes the vertical diffusion and computes moisture
!     related terms (convection,large-scale condensation, and precipitation).
! </DESCRIPTION>
! <DIAGFIELDS>
! </DIAGFIELDS>
! <DATASET NAME="physics_driver.res">
! native format restart file
! </DATASET>
!
! <DATASET NAME="physics_driver.res.nc">
! netcdf format restart file
! </DATASET>


! <INFO>

!   <REFERENCE>            </REFERENCE>
!   <COMPILER NAME="">     </COMPILER>
!   <PRECOMP FLAG="">      </PRECOMP>
!   <LOADER FLAG="">       </LOADER>
!   <TESTPROGRAM NAME="">  </TESTPROGRAM>
!   <BUG>                  </BUG>
!   <NOTE> 
!   </NOTE>
!   <FUTURE> Deal with conservation of total energy?              </FUTURE>

! </INFO>
!   shared modules:

use time_manager_mod,        only: time_type, get_time, operator (-), &
                                   time_manager_init, operator(*)
use field_manager_mod,       only: field_manager_init, MODEL_ATMOS
use tracer_manager_mod,      only: tracer_manager_init, &
                                   get_number_tracers, &
                                   get_tracer_names, &
                                   get_tracer_index, NO_TRACER
use block_control_mod,       only: block_control_type
use atmos_tracer_driver_mod, only: atmos_tracer_driver_init,    &
                                   atmos_tracer_driver_time_vary, &
                                   atmos_tracer_driver_endts, &
                                   atmos_tracer_driver,  &
                                   atmos_tracer_driver_end
use mpp_mod,                 only: input_nml_file
use fms_mod,                 only: mpp_clock_id, mpp_clock_begin,   &
                                   mpp_clock_end, CLOCK_MODULE_DRIVER, &
                                   fms_init,  &
                                   open_namelist_file, stdlog, stdout,  &
                                   write_version_number, field_size, &
                                   file_exist, error_mesg, FATAL,   &
                                   WARNING, NOTE, check_nml_error, &
                                   close_file, mpp_pe, mpp_root_pe, &
                                   mpp_error, mpp_chksum, string
use fms_io_mod,              only: restore_state, &
                                   register_restart_field, restart_file_type, &
                                   save_restart, get_mosaic_tile_file
use constants_mod,           only: RDGAS

use diag_manager_mod,        only: register_diag_field, send_data

!    shared radiation package modules:

use aerosol_types_mod,       only: aerosol_type, aerosol_time_vary_type

use physics_radiation_exch_mod, only: exchange_control_type, &
                                      clouds_from_moist_type, &
                                      clouds_from_moist_block_type, &
                                      cosp_from_rad_type, &
                                      cosp_from_rad_control_type, &
                                      cosp_from_rad_block_type, &
                                      radiation_flux_control_type, & 
                                      radiation_flux_block_type, & 
                                      alloc_clouds_from_moist_type, &
                                      alloc_cloud_scheme_data_type

use physics_types_mod,       only: alloc_physics_tendency_type, &
                                   physics_tendency_type, & 
                                   physics_tendency_block_type, &
                                   physics_type, & 
                                   physics_control_type, & 
                                   physics_input_block_type, &
                                   dealloc_physics_tendency_type
  
use aerosol_mod,             only: aerosol_init, aerosol_driver, &
                                   aerosol_time_vary, &
                                   aerosol_endts, &
                                   aerosol_dealloc, aerosol_end
!    component modules:

use cosp_driver_mod,         only: cosp_driver_init, cosp_driver, &
                                   cosp_driver_end, cosp_driver_time_vary, &
                                   cosp_driver_endts
use  moist_processes_mod,    only: moist_processes,    &
                                   moist_processes_init,  &
                                   set_cosp_precip_sources, &
                                   moist_processes_time_vary, &
                                   moist_processes_endts, &
                                   moist_processes_end

use vert_turb_driver_mod,    only: vert_turb_driver,  &
                                   vert_turb_driver_init,  &
                                   vert_turb_driver_end, &
                                   vert_turb_driver_restart

use vert_diff_driver_mod,    only: vert_diff_driver_down,  &
                                   vert_diff_driver_up,    &
                                   vert_diff_driver_init,  &
                                   vert_diff_driver_end,   &
                                   surf_diff_type
 
use damping_driver_mod,      only: damping_driver,      &
                                   damping_driver_init, &
                                   damping_driver_time_vary,  &
                                   damping_driver_endts, &
                                   damping_driver_end,  &
                                   damping_driver_restart

use grey_radiation_mod,       only: grey_radiation_init, grey_radiation, &
                                    grey_radiation_end

use clubb_driver_mod,         only: clubb_init, clubb, clubb_end
use monin_obukhov_mod,        only: monin_obukhov_init
#ifdef SCM
! Option to add SCM radiative tendencies from forcing to lw_tendency
! and radturbten

use scm_forc_mod,            only: use_scm_rad, add_scm_tdtlw, add_scm_tdtsw

#endif

!-----------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    physics_driver_mod accesses the model's physics modules and
!    obtains tendencies and boundary fluxes due to the physical
!    processes that drive atmospheric time tendencies and supply 
!    boundary forcing to the surface models.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public  physics_driver_init, physics_driver_down,   &
        physics_driver_down_time_vary, physics_driver_up_time_vary, &
        physics_driver_down_endts, physics_driver_up_endts, &
        physics_driver_up, physics_driver_end, &
        do_moist_in_phys_up, get_diff_t, &
        get_radturbten, zero_radturbten, physics_driver_restart, &
        cosp_driver_init, set_cosp_precip_sources

private          &

!  called from physics_driver_down:
         check_args, &

!  called from check_args:
         check_dim

interface check_dim
     module procedure check_dim_2d, check_dim_3d, check_dim_4d
end interface

!---------------------------------------------------------------------
!------- namelist ------

logical :: do_moist_processes = .true.
                               ! call moist_processes routines
real    :: tau_diff = 3600.    ! time scale for smoothing diffusion 
                               ! coefficients
integer :: do_clubb = 0        ! activate clubb parameterization ?
logical :: do_cosp = .false.   ! activate COSP simulator ?
logical :: do_modis_yim = .true. ! activate simple modis simulator ?
logical :: do_radiation = .true.
                               ! calculating radiative fluxes and
                               ! heating rates?
logical :: do_grey_radiation = .false. ! do grey radiation scheme?
real    :: R1 = 0.25           ! rif:(09/10/09) In Grey radiation we are computing just the total   
real    :: R2 = 0.25           ! SW radiation. We need to divide it into 4 components
real    :: R3 = 0.25           ! to go through the Coupler and Ice modules.
real    :: R4 = 0.25           ! 	Sum[R(i)*SW] = SW  


real    :: diff_min = 1.e-3    ! minimum value of a diffusion 
                               ! coefficient beneath which the
                               ! coefficient is reset to zero
logical :: diffusion_smooth = .true.
                               ! diffusion coefficients should be 
                               ! smoothed in time?
logical :: donner_meso_is_largescale = .true.
                               ! donner meso clouds are treated as 
                               ! largescale (rather than convective)
                               ! as far as the COSP simulator is 
                               ! concerned ?
logical :: allow_cosp_precip_wo_clouds = .true.
                               ! COSP will see {ls, cv} precip in grid-
                               ! boxes w/o {ls, cv} clouds ?
logical :: override_aerosols_cloud = .false.
                               ! use offline aerosols for cloud calculation
                               ! (via data_override in aerosol_driver)?

! <NAMELIST NAME="physics_driver_nml">
!  <DATA NAME="do_radiation" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!calculating radiative fluxes and
! heating rates?
!  </DATA>
!  <DATA NAME="do_moist_processes" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!call moist_processes routines
!  </DATA>
!  <DATA NAME="tau_diff" UNITS="" TYPE="real" DIM="" DEFAULT="3600.">
!time scale for smoothing diffusion 
! coefficients
!  </DATA>
!  <DATA NAME="diff_min" UNITS="" TYPE="real" DIM="" DEFAULT="1.e-3">
!minimum value of a diffusion 
! coefficient beneath which the
! coefficient is reset to zero
!  </DATA>
!  <DATA NAME="diffusion_smooth" UNITS="" TYPE="logical" DIM="" DEFAULT=".true.">
!diffusion coefficients should be 
! smoothed in time?
!  </DATA>
!  <DATA NAME="override_aerosols_cloud" UNITS="" TYPE="logical" DIM="" DEFAULT=".false.">
!use offline aerosols for cloud calculation
! (via data_override in aerosol_driver)?
!  </DATA>
! </NAMELIST>
!

! ---> h1g, 2012-08-28, add option of applying surface fluxes in host-model
!                       by default .true. (that is, applying surface fluxes in host-model)
logical :: l_host_applies_sfc_fluxes = .true.
! <--- h1g, 2012-08-28

namelist / physics_driver_nml / do_radiation, &
                                do_clubb, & 
                                do_cosp, &
                                do_modis_yim, &
                                donner_meso_is_largescale, &
                                allow_cosp_precip_wo_clouds, &
                                do_moist_processes, tau_diff,      &
                                diff_min, diffusion_smooth, &
                                do_grey_radiation, R1, R2, R3, R4,  &
                                override_aerosols_cloud, &
                                l_host_applies_sfc_fluxes

!---------------------------------------------------------------------
!------- public data ------
! <DATA NAME="surf_diff_type" UNITS="" TYPE="surf_diff_type" DIM="" DEFAULT="">
! Defined in vert_diff_driver_mod, republished here. See vert_diff_mod for details.
! </DATA>

public  surf_diff_type   ! defined in  vert_diff_driver_mod, republished
                         ! here

!---------------------------------------------------------------------
!------- private data ------

!--------------------------------------------------------------------
! list of restart versions readable by this module:
!
! version 1: initial implementation 1/2003, contains diffusion coef-
!            ficient contribution from cu_mo_trans_mod. This variable
!            is generated in physics_driver_up (moist_processes) and
!            used on the next step in vert_diff_down, necessitating
!            its storage.
!
! version 2: adds pbltop as generated in vert_turb_driver_mod. This 
!            variable is then used on the next timestep by topo_drag
!            (called from damping_driver_mod), necessitating its 
!            storage.
!
! version 3: adds the diffusion coefficients which are passed to 
!            vert_diff_driver.  These diffusion are saved should
!            smoothing of vertical diffusion coefficients be turned
!            on.
!
! version 4: adds a logical variable, convect, which indicates whether
!            or not the grid column is convecting. This diagnostic is
!            needed by the entrain_module in vert_turb_driver.
!
! version 5: adds radturbten when strat_cloud_mod is active, adds 
!            lw_tendency when edt_mod or entrain_mod is active.
!
! version 6: adds donner cell and meso cloud variables when donner_deep
!            is activated.

! version 7: adds shallow convection cloud variables when uw_conv
!            is activated.

! version 8: adds lsc cloud props for radiation. only readable when in
!            netcdf mode.


!---------------------------------------------------------------------
integer, dimension(8) :: restart_versions = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)

!--------------------------------------------------------------------
!    the following allocatable arrays are either used to hold physics 
!    data between timesteps when required, or hold physics data between
!    physics_down and physics_up.
!  
!    diff_cu_mo     contains contribution to difusion coefficient
!                   coming from cu_mo_trans_mod (called from 
!                   moist_processes in physics_driver_up) and then used 
!                   as input on the next time step to vert_diff_down 
!                   called in physics_driver_down.
!    diff_t         vertical diffusion coefficient for temperature
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    diff_m         vertical diffusion coefficient for momentum
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    radturbten     the sum of the radiational and turbulent heating,
!                   generated in both physics_driver_down (radiation)
!                   and physics_driver_up (turbulence) and then used
!                   in moist_processes
!    pbltop         top of boundary layer obtained from vert_turb_driver
!                   and then used on the next timestep in topo_drag_mod
!                   called from damping_driver_down        
!    convect        flag indicating whether convection is occurring in
!                   a grid column. generated in physics_driver_up and
!                   then used in vert_turb_driver called from 
!                   physics_driver_down on the next step.
!----------------------------------------------------------------------
real,    dimension(:,:,:), allocatable :: diff_cu_mo, diff_t, diff_m
real,    dimension(:,:,:), allocatable :: radturbten
real,    dimension(:,:)  , allocatable :: pbltop, cush, cbmf, hmint, cgust !miz 
real,    dimension(:,:)  , allocatable :: tke, pblhto, rkmo, taudpo        !miz
integer, dimension(:,:,:), allocatable :: exist_shconv, exist_dpconv       !miz

logical, dimension(:,:)  , allocatable :: convect

real,    dimension(:,:,:), allocatable ::       &
                           temp_last, q_last

real,    dimension(:,:,:), allocatable ::  fl_lsrain, fl_lssnow, &
                                           fl_lsgrpl, &
                                           fl_donmca_rain, fl_donmca_snow,&
                                           fl_ccrain, fl_ccsnow

real   ,    dimension(:,:)  , allocatable :: tsurf_save

real,    dimension(:,:,:), allocatable ::  diff_t_clubb
   
!--- for netcdf restart
type(restart_file_type), pointer, save :: Phy_restart => NULL()
type(restart_file_type), pointer, save :: Til_restart => NULL()
logical                                :: in_different_file = .false.
integer                                :: vers
integer                                :: now_doing_strat = 0
integer                                :: now_doing_entrain = 0
integer                                :: now_doing_edt = 0
real, allocatable                      :: r_convect(:,:)

type(aerosol_time_vary_type) :: Aerosol_cld

!---------------------------------------------------------------------
!    internal timing clock variables:
!---------------------------------------------------------------------
integer :: damping_clock, turb_clock,   &
           tracer_clock, diff_up_clock, diff_down_clock, &
           moist_processes_clock, cosp_clock

!--------------------------------------------------------------------
!    miscellaneous control variables:
!---------------------------------------------------------------------
logical   :: do_check_args = .true.   ! argument dimensions should 
                                      ! be checked ?
logical   :: module_is_initialized = .false.
                                      ! module has been initialized ?
logical   :: doing_edt                ! edt_mod has been activated ?
logical   :: doing_entrain            ! entrain_mod has been activated ?
logical   :: doing_strat              ! stratiform clouds has been activated ?
logical   :: doing_donner             ! donner_deep_mod has been 
                                      ! activated ?
logical   :: doing_uw_conv            ! uw_conv shallow cu mod has been 
                                      ! activated ?
logical   :: doing_liq_num = .false.  ! Prognostic cloud droplet number has 
                                      ! been activated?
integer   :: nt                       ! total no. of tracers
integer   :: ntp                      ! total no. of prognostic tracers
integer   :: ncol                     ! number of stochastic columns
 
integer   :: num_uw_tracers


!---------------------------------------------------------------------
!---------------------------------------------------------------------

character(len=4)     :: mod_name = 'phys'
character(len=32)    :: tracer_units, tracer_name
  character(len=128) :: diaglname
real                 :: missing_value = -999.
logical              :: include_donmca_in_cosp

integer                            :: id_tdt_phys,         &
                                      id_tdt_phys_vdif_dn, &
                                      id_tdt_phys_vdif_up, &
                                      id_tdt_phys_turb,    &
                                      id_tdt_phys_moist

integer, dimension(:), allocatable :: id_tracer_phys,         &
                                      id_tracer_phys_vdif_dn, &
                                      id_tracer_phys_vdif_up, &
                                      id_tracer_phys_turb,    &
                                      id_tracer_phys_moist

type (clouds_from_moist_block_type) :: Restart


                            contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="physics_driver_init">
!  <OVERVIEW>
!    physics_driver_init is the constructor for physics_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_init is the constructor for physics_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_init (Time, lonb, latb, axes, pref, &
!                             trs, Surf_diff, phalf )
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="pref" TYPE="real">
!   reference prssure profiles
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   array of model latitudes at cell corners [radians]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   array of model longitudes at cell corners [radians]
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   axis indices, (/x,y,pf,ph/)
!                (returned from diag axis manager)
!  </IN>
!  <INOUT NAME="trs" TYPE="real">
!   atmospheric tracer fields
!  </INOUT>
!  <INOUT NAME="Surf_diff" TYPE="surf_diff_type">
!   surface diffusion derived type
!  </INOUT>
!  <IN NAME="phalf" TYPE="real">
!   pressure at model interface levels
!  </IN>
! <ERROR MSG="physics_driver_init must be called first" STATUS="FATAL">
! </ERROR>
! </SUBROUTINE>
!
subroutine physics_driver_init (Time, lonb, latb, lon, lat, axes, &
                                Surf_diff, Exch_ctrl, Atm_block, Moist_clouds, &
                                Physics, Physics_tendency, diffm, difft)

!---------------------------------------------------------------------
!    physics_driver_init is the constructor for physics_driver_mod.
!---------------------------------------------------------------------

type(time_type),              intent(in)    :: Time
real,    dimension(:,:),      intent(in)    :: lonb, latb
real,    dimension(:,:),      intent(in)    :: lon, lat
integer, dimension(4),        intent(in)    :: axes
type(surf_diff_type),         intent(inout) :: Surf_diff
type (exchange_control_type), intent(inout) :: Exch_ctrl
type (block_control_type),    intent(in)    :: Atm_block
type(clouds_from_moist_type), intent(inout) :: Moist_clouds(:)
type(physics_tendency_type),  intent(inout) :: Physics_tendency
type(physics_type),           intent(inout) :: Physics
real,    dimension(:,:,:),    intent(out),  optional :: diffm, difft

!---------------------------------------------------------------------
!  intent(in) variables:
!     Time       current time (time_type)
!     lonb       longitude of the grid box corners [ radians ]
!     latb       latitude of the grid box corners [ radians ]
!     axes       axis indices, (/x,y,pf,ph/)
!                (returned from diag axis manager)
!
!   intent(inout) variables:
!     Surf_diff  surface diffusion derived type variable
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      integer :: nb, ibs, ibe, jbs, jbe
      real, dimension (size(lonb,1)-1, size(latb,2)-1) :: sgsmtn
      integer          ::  id, jd, kd, n, k, nc
      integer          ::  ierr, io, unit, logunit, outunit
      integer          ::  ndum
      character(len=16)::  cosp_precip_sources_modified

      integer          ::  moist_processes_init_clock, damping_init_clock, &
                           turb_init_clock, diff_init_clock, &
                           aerosol_init_clock, &
                           grey_radiation_init_clock , &
                           tracer_init_clock
      real, dimension(:,:,:),   allocatable :: phalf
      real, dimension(:,:,:,:), allocatable :: trs
!---------------------------------------------------------------------
!  local variables:
!
!       sgsmtn        sgs orography obtained from mg_drag_mod;
!                     appears to not be currently used
!       aerosol_names names associated with the activated aerosols
!                     that will be seen by the radiation package
!       aerosol_family_names
!              names associated with the activated aerosol
!              families that will be seen by the radiation package
!       id,jd,kd      model dimensions on the processor  
!       ierr          error code
!       io            io status returned from an io call
!       unit          unit number used for an i/ operation

!---------------------------------------------------------------------
!    if routine has already been executed, return.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that the modules used by this module that are not called 
!    later in this subroutine have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call tracer_manager_init
      call field_manager_init (ndum)
 
!--------------------------------------------------------------------
!    read namelist.
!--------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=physics_driver_nml, iostat=io)
      ierr = check_nml_error(io,"physics_driver_nml")
#else
      if ( file_exist('input.nml')) then
        unit = open_namelist_file ()
        ierr=1; do while (ierr /= 0)
        read  (unit, nml=physics_driver_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'physics_driver_nml')
        enddo
10      call close_file (unit)
      endif
#endif

      if(do_radiation .and. do_grey_radiation) & 
        call error_mesg('physics_driver_init','do_radiation and do_grey_radiation cannot both be .true.',FATAL)

!--------------------------------------------------------------------
!    write version number and namelist to log file.
!--------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
               write(logunit, nml=physics_driver_nml)
 
!---------------------------------------------------------------------
!    define the model dimensions on the local processor.
!---------------------------------------------------------------------
      id = size(lonb,1)-1 
      jd = size(latb,2)-1 
      kd = Atm_block%npz
      call get_number_tracers (MODEL_ATMOS, num_tracers=nt, &
                               num_prog=ntp)

!---------------------------------------------------------------------
      cosp_clock       =       &
                mpp_clock_id( '   Physics_down: COSP',    &
                   grain=CLOCK_MODULE_DRIVER )
      damping_clock         =     &
                mpp_clock_id( '   Physics_down: Damping',    &
                  grain=CLOCK_MODULE_DRIVER )
      turb_clock            =      &
                mpp_clock_id( '   Physics_down: Vert. Turb.', &
                  grain=CLOCK_MODULE_DRIVER )
      tracer_clock          =      &
                mpp_clock_id( '   Physics_down: Tracer',    &
                 grain=CLOCK_MODULE_DRIVER )
      diff_down_clock       =     &
                mpp_clock_id( '   Physics_down: Vert. Diff.',   &
                 grain=CLOCK_MODULE_DRIVER )
      diff_up_clock         =     &
                mpp_clock_id( '   Physics_up: Vert. Diff.',     &
                grain=CLOCK_MODULE_DRIVER )
      moist_processes_clock =      &
                mpp_clock_id( '   Physics_up: Moist Processes', &
                grain=CLOCK_MODULE_DRIVER )

      moist_processes_init_clock =      &
        mpp_clock_id( '   Physics_driver_init: Moist Processes: Initialization', &
                grain=CLOCK_MODULE_DRIVER )
      damping_init_clock         =     &
        mpp_clock_id( '   Physics_driver_init: Damping: Initialization',    &
                  grain=CLOCK_MODULE_DRIVER )
      turb_init_clock            =      &
        mpp_clock_id( '   Physics_driver_init: Vert. Turb.: Initialization', &
                  grain=CLOCK_MODULE_DRIVER )
      diff_init_clock       =     &
        mpp_clock_id( '   Physics_driver_init: Vert. Diff.: Initialization',   &
                 grain=CLOCK_MODULE_DRIVER )
      aerosol_init_clock       =       &
        mpp_clock_id( '   Physics_driver_init: Aerosol: Initialization', &
                       grain=CLOCK_MODULE_DRIVER )
      grey_radiation_init_clock       =       &
        mpp_clock_id( '   Physics_driver_init: Grey Radiation: Initialization', &
                       grain=CLOCK_MODULE_DRIVER )
      tracer_init_clock          =      &
        mpp_clock_id( '   Physics_driver_init: Tracer: Initialization',    &
                 grain=CLOCK_MODULE_DRIVER )

!-----------------------------------------------------------------------
!--- allocate Physics_tendency to hold the physics tendencies
!-----------------------------------------------------------------------
      call alloc_physics_tendency_type (Physics_tendency, Atm_block)

!--- define trs and p_half on the full domain 
    allocate (trs(id,jd,kd,nt), phalf(id,jd,kd+1))
    do nb = 1, Atm_block%nblks
      ibs = Atm_block%ibs(nb)-Atm_block%isc+1
      ibe = Atm_block%ibe(nb)-Atm_block%isc+1
      jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
      jbe = Atm_block%jbe(nb)-Atm_block%jsc+1
      trs(ibs:ibe,jbs:jbe,:,1:ntp)    = Physics%block(nb)%q
      trs(ibs:ibe,jbs:jbe,:,ntp+1:nt) = Physics%block(nb)%tmp_4d
      phalf(ibs:ibe,jbs:jbe,:)        = Physics%block(nb)%p_half
!--- the 'temp' variable inside of Physics is no longer needed - deallocate it
      deallocate(Physics%block(nb)%tmp_4d)
    enddo

!-----------------------------------------------------------------------
!---------- initialize physics -------
    if (do_moist_processes) then
      call mpp_clock_begin ( moist_processes_init_clock )
      call moist_processes_init (id, jd, kd, lonb, latb, lon, lat, phalf, &
                                 Physics%glbl_qty%pref(:,1),&
                                 axes, Time, doing_donner,  &
                                 doing_uw_conv,  &
                                 num_uw_tracers, doing_strat, &
                                 do_clubb_in=do_clubb, &
                                 do_cosp_in=do_cosp, &
                                 donner_meso_is_largescale_in= &
                                         donner_meso_is_largescale, &
                                 include_donmca_in_cosp_out = &
                                         include_donmca_in_cosp)

      call mpp_clock_end ( moist_processes_init_clock )
    else
      diff_cu_mo = 0.0
      convect = .false.
    endif
     
!-----------------------------------------------------------------------
!    initialize damping_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( damping_init_clock )
      call damping_driver_init (lonb, latb, Physics%glbl_qty%pref(:,1), axes, Time, &
                                sgsmtn)
      call mpp_clock_end ( damping_init_clock )

!-----------------------------------------------------------------------
!    initialize vert_turb_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( turb_init_clock )
      call vert_turb_driver_init (lonb, latb, id, jd, kd, axes, Time, &
                                  doing_edt, doing_entrain, do_clubb)
      call mpp_clock_end ( turb_init_clock )

!-----------------------------------------------------------------------
!    initialize vert_diff_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( diff_init_clock )
      call vert_diff_driver_init (Surf_diff, id, jd, kd, axes, Time, do_clubb )
      call mpp_clock_end ( diff_init_clock )

      if (do_moist_processes) then
!-----------------------------------------------------------------------
!    initialize aerosol_mod     
!-----------------------------------------------------------------------
        call mpp_clock_begin ( aerosol_init_clock )
        call aerosol_init (lonb, latb, Aerosol_cld)
        call mpp_clock_end ( aerosol_init_clock )
      endif ! do_moist_processes

      if(do_grey_radiation) then
         call mpp_clock_begin ( grey_radiation_init_clock )
         call grey_radiation_init(axes, Time) 
         call mpp_clock_end ( grey_radiation_init_clock )
      endif
        
!-----------------------------------------------------------------------
!    initialize atmos_tracer_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( tracer_init_clock )
      call atmos_tracer_driver_init (lonb, latb, trs, axes, Time, phalf)
      call mpp_clock_end ( tracer_init_clock )

!---------------------------------------------------------------------
!    allocate space for the module variables.
!---------------------------------------------------------------------
      allocate ( diff_t     (id, jd, kd) ) ; diff_t = 0.0
      allocate ( diff_m     (id, jd, kd) ) ; diff_m = 0.0
      allocate ( diff_cu_mo (id, jd, kd) ) ; diff_cu_mo = 0.0
      allocate ( pbltop     (id, jd) )     ; pbltop     = -999.0
      allocate ( cush       (id, jd) )     ; cush=-1. !miz
      allocate ( cbmf       (id, jd) )     ; cbmf=0.0 !miz
      allocate ( hmint      (id, jd) )     ; hmint=0. !miz
      allocate ( cgust      (id, jd) )     ; cgust=0.0 !miz
      allocate ( tke        (id, jd) )     ; tke  =0.0 !miz
      allocate ( pblhto     (id, jd) )     ; pblhto=0.0 !miz
      allocate ( rkmo       (id, jd) )     ; rkmo=15.0 !miz
      allocate ( taudpo     (id, jd) )     ; taudpo=28800.    !miz
      allocate ( exist_shconv(id, jd,48) ) ; exist_shconv = 0 !miz
      allocate ( exist_dpconv(id, jd,48) ) ; exist_dpconv = 0 !miz
      allocate ( convect    (id, jd) )     ; convect = .false.
      allocate ( radturbten (id, jd, kd))  ; radturbten = 0.0
      allocate ( r_convect  (id, jd) )     ; r_convect   = 0.0
       
      if (do_clubb > 0 ) then
        allocate ( diff_t_clubb(id, jd, kd) );      diff_t_clubb      = 0.0
      end if

!--------------------------------------------------------------------
!    these variables needed to preserve rain fluxes, q and T from end 
!    of one step for use in COSP simulator on next step.
!--------------------------------------------------------------------
      allocate (fl_lsrain  (id, jd, kd))
      allocate (fl_lssnow  (id, jd, kd))
      allocate (fl_lsgrpl  (id, jd, kd))
      allocate (fl_ccrain  (id, jd, kd))
      allocate (fl_ccsnow  (id, jd, kd))
      allocate (fl_donmca_snow  (id, jd, kd))
      allocate (fl_donmca_rain  (id, jd, kd))
      allocate ( temp_last (id, jd, kd))
      allocate ( q_last    (id, jd, kd))
      fl_lsrain = 0.
      fl_lssnow = 0.
      fl_lsgrpl = 0.
      fl_ccrain = 0.
      fl_ccsnow = 0.
      fl_donmca_rain = 0.
      fl_donmca_snow = 0.
      temp_last = 0.
      q_last    = 0.

      if (do_cosp .or. do_modis_yim) then
        allocate ( tsurf_save (id, jd))
        tsurf_save = 0.
      endif

!-----------------------------------------------------------------------
!--- store some control variables for radiation and cosp
!-----------------------------------------------------------------------
      Exch_ctrl%doing_strat               = doing_strat
      Exch_ctrl%doing_donner              = doing_donner
      Exch_ctrl%doing_uw_conv             = doing_uw_conv
      Exch_ctrl%donner_meso_is_largescale = donner_meso_is_largescale
      Exch_ctrl%do_cosp                   = do_cosp
      Exch_ctrl%do_modis_yim              = do_modis_yim

      ! count the number of cloud schemes
      Exch_ctrl%ncld = 0
      if (Exch_ctrl%doing_strat)   Exch_ctrl%ncld = Exch_ctrl%ncld + 1
      if (Exch_ctrl%doing_donner)  Exch_ctrl%ncld = Exch_ctrl%ncld + 2
      if (Exch_ctrl%doing_uw_conv) Exch_ctrl%ncld = Exch_ctrl%ncld + 1
!-----------------------------------------------------------------------
!    allocate derived-type that stores cloud properties
!    return from moist processes
!-----------------------------------------------------------------------
     call error_mesg('physics_driver_mod', 'number of cloud schemes found = '//trim(string(Exch_ctrl%ncld)), NOTE)

     call alloc_clouds_from_moist_type(Moist_clouds, Exch_ctrl, Atm_block)

!--------------------------------------------------------------------
!    call physics_driver_read_restart to obtain initial values for the module
!    variables. Also register restart fields to be ready for intermediate 
!    restart.
!--------------------------------------------------------------------
      allocate(Restart%Cloud_data(Exch_ctrl%ncld))

      do nc = 1, Exch_ctrl%ncld
        ! restart values allocated on the full domain
        call alloc_cloud_scheme_data_type( Moist_clouds(1)%block(1)%Cloud_data(nc)%scheme_name, &
                                           id, jd, kd, Restart%Cloud_data(nc))
      enddo

      call physics_driver_register_restart (Restart)
      if(file_exist('INPUT/physics_driver.res.nc')) then
         call restore_state(Phy_restart)
         if(in_different_file) call restore_state(Til_restart)
      endif
!---------------------------------------------------------------------
!    a flag indicating columns in which convection is occurring is
!    present beginning with v4. if not present, set it to .false.
!---------------------------------------------------------------------
      convect = .false.
      where(r_convect .GT. 0.) 
         convect = .true.
      end where
         
100 FORMAT("CHECKSUM::",A32," = ",Z20)
      outunit = stdout()
      write(outunit,*) 'BEGIN CHECKSUM(physics_driver_init):: '
      write(outunit,100) 'diff_cu_mo             ', mpp_chksum(diff_cu_mo            )
      write(outunit,100) 'pbltop                 ', mpp_chksum(pbltop                )
      write(outunit,100) 'cush                   ', mpp_chksum(cush                  )
      write(outunit,100) 'cbmf                   ', mpp_chksum(cbmf                  )
      write(outunit,100) 'hmint                  ', mpp_chksum(hmint                 )
      write(outunit,100) 'cgust                  ', mpp_chksum(cgust                 )
      write(outunit,100) 'tke                    ', mpp_chksum(tke                   )
      write(outunit,100) 'pblhto                 ', mpp_chksum(pblhto                )
      write(outunit,100) 'rkmo                   ', mpp_chksum(rkmo                  )
      write(outunit,100) 'taudpo                 ', mpp_chksum(taudpo                )
      write(outunit,100) 'exist_shconv           ', mpp_chksum(exist_shconv          )
      write(outunit,100) 'exist_dpconv           ', mpp_chksum(exist_dpconv          )
      write(outunit,100) 'diff_t                 ', mpp_chksum(diff_t                )
      write(outunit,100) 'diff_m                 ', mpp_chksum(diff_m                )
      write(outunit,100) 'r_convect              ', mpp_chksum(r_convect             )
  if ( doing_strat ) then
      write(outunit,100) 'radturbten             ', mpp_chksum(radturbten            )
  endif
      do nc = 1, size(Restart%Cloud_data,1)
        ! NOTE: the order of the checksums in stdout will be different
        if ( trim(Restart%Cloud_data(nc)%scheme_name).eq.'donner_cell' ) then
          write(outunit,100) 'cell_cld_frac          ', mpp_chksum(Restart%Cloud_data(nc)%cloud_area )
          write(outunit,100) 'cell_liq_amt           ', mpp_chksum(Restart%Cloud_data(nc)%liquid_amt )
          write(outunit,100) 'cell_liq_size          ', mpp_chksum(Restart%Cloud_data(nc)%liquid_size)
          write(outunit,100) 'cell_ice_amt           ', mpp_chksum(Restart%Cloud_data(nc)%ice_amt    )
          write(outunit,100) 'cell_ice_size          ', mpp_chksum(Restart%Cloud_data(nc)%ice_size   )
        endif
        if ( trim(Restart%Cloud_data(nc)%scheme_name).eq.'donner_meso' ) then
          write(outunit,100) 'meso_cld_frac          ', mpp_chksum(Restart%Cloud_data(nc)%cloud_area )
          write(outunit,100) 'meso_liq_amt           ', mpp_chksum(Restart%Cloud_data(nc)%liquid_amt )
          write(outunit,100) 'meso_liq_size          ', mpp_chksum(Restart%Cloud_data(nc)%liquid_size)
          write(outunit,100) 'meso_ice_amt           ', mpp_chksum(Restart%Cloud_data(nc)%ice_amt    )
          write(outunit,100) 'meso_ice_size          ', mpp_chksum(Restart%Cloud_data(nc)%ice_size   )
          write(outunit,100) 'nsum_out               ', mpp_chksum(Restart%Cloud_data(nc)%nsum_out   )
        endif
        if ( trim(Restart%Cloud_data(nc)%scheme_name).eq.'uw_conv' ) then
          write(outunit,100) 'shallow_cloud_area     ', mpp_chksum(Restart%Cloud_data(nc)%cloud_area    )
          write(outunit,100) 'shallow_liquid         ', mpp_chksum(Restart%Cloud_data(nc)%liquid_amt    )
          write(outunit,100) 'shallow_ice            ', mpp_chksum(Restart%Cloud_data(nc)%ice_amt       )
          write(outunit,100) 'shallow_droplet_number ', mpp_chksum(Restart%Cloud_data(nc)%droplet_number)
          write(outunit,100) 'shallow_ice_number     ', mpp_chksum(Restart%Cloud_data(nc)%ice_number    )
        endif
        if ( trim(Restart%Cloud_data(nc)%scheme_name).eq.'strat_cloud' ) then
          write(outunit,100) 'lsc_cloud_area         ', mpp_chksum(Restart%Cloud_data(nc)%cloud_area    )
          write(outunit,100) 'lsc_liquid             ', mpp_chksum(Restart%Cloud_data(nc)%liquid_amt    )
          write(outunit,100) 'lsc_ice                ', mpp_chksum(Restart%Cloud_data(nc)%ice_amt       )
          write(outunit,100) 'lsc_droplet_number     ', mpp_chksum(Restart%Cloud_data(nc)%droplet_number)
          write(outunit,100) 'lsc_ice_number         ', mpp_chksum(Restart%Cloud_data(nc)%ice_number    )
          write(outunit,100) 'lsc_snow               ', mpp_chksum(Restart%Cloud_data(nc)%snow          )
          write(outunit,100) 'lsc_rain               ', mpp_chksum(Restart%Cloud_data(nc)%rain          )
          write(outunit,100) 'lsc_snow_size          ', mpp_chksum(Restart%Cloud_data(nc)%snow_size     )
          write(outunit,100) 'lsc_rain_size          ', mpp_chksum(Restart%Cloud_data(nc)%rain_size     )
        endif
      enddo ! nc

      do nb = 1, Atm_block%nblks
        ibs = Atm_block%ibs(nb)-Atm_block%isc+1
        ibe = Atm_block%ibe(nb)-Atm_block%isc+1
        jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
        jbe = Atm_block%jbe(nb)-Atm_block%jsc+1

        !-- copy cloud data from restart
        do nc = 1, size(Restart%Cloud_data,1)

          ! common to all cloud schemes
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%cloud_area     = Restart%Cloud_data(nc)%cloud_area    (ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_amt     = Restart%Cloud_data(nc)%liquid_amt    (ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_amt        = Restart%Cloud_data(nc)%ice_amt       (ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%droplet_number = Restart%Cloud_data(nc)%droplet_number(ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name    = Restart%Cloud_data(nc)%scheme_name

          ! properties specific to large-scale/stratiform clouds
          if (trim(Restart%Cloud_data(nc)%scheme_name) .eq. 'strat_cloud') then
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number = Restart%Cloud_data(nc)%ice_number (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain       = Restart%Cloud_data(nc)%rain       (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow       = Restart%Cloud_data(nc)%snow       (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain_size  = Restart%Cloud_data(nc)%rain_size  (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow_size  = Restart%Cloud_data(nc)%snow_size  (ibs:ibe,jbs:jbe,:)
          endif
  
          ! properties specific to donner deep clouds (both cell and meso)
          if (trim(Restart%Cloud_data(nc)%scheme_name) .eq. 'donner_cell' .or. &
              trim(Restart%Cloud_data(nc)%scheme_name) .eq. 'donner_meso') then
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_size = Restart%Cloud_data(nc)%liquid_size (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_size    = Restart%Cloud_data(nc)%ice_size    (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%nsum_out    = Restart%Cloud_data(nc)%nsum_out    (ibs:ibe,jbs:jbe)
          endif

          ! properties specific to uw shallow convective clouds
          if (trim(Restart%Cloud_data(nc)%scheme_name) .eq. 'uw_conv') then
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number = Restart%Cloud_data(nc)%ice_number (ibs:ibe,jbs:jbe,:)
          endif

        enddo

        !--- return trs to the blocked data structure
        Physics%block(nb)%q = trs(ibs:ibe,jbs:jbe,:,1:ntp)
        Physics_tendency%block(nb)%qdiag = trs(ibs:ibe,jbs:jbe,:,ntp+1:nt)
      enddo
      deallocate (trs, phalf)

      vers = restart_versions(size(restart_versions(:)))

!---------------------------------------------------------------------
!    if desired, define variables to return diff_m and diff_t.
!---------------------------------------------------------------------
      if (present(difft)) then
        difft = diff_t
      endif
      if (present(diffm)) then
        diffm = diff_m
      endif

!---------------------------------------------------------------------
!    initialize module diagnostics
!---------------------------------------------------------------------

      id_tdt_phys_vdif_dn = register_diag_field ( mod_name,    &
         'tdt_phys_vdif_dn', axes(1:3), Time,                  &
         'temperature tendency from physics driver vdif down', &
         'K/s', missing_value=missing_value)

      id_tdt_phys_vdif_up = register_diag_field ( mod_name,    &
         'tdt_phys_vdif_up', axes(1:3), Time,                  &
         'temperature tendency from physics driver vdif up',   &
         'K/s', missing_value=missing_value)

      id_tdt_phys_turb = register_diag_field ( mod_name,       &
         'tdt_phys_turb', axes(1:3), Time,                     &
         'temperature tendency from physics driver vdif turb', &
         'K/s', missing_value=missing_value)

      id_tdt_phys_moist = register_diag_field ( mod_name,            &
         'tdt_phys_moist', axes(1:3), Time,                          &
         'temperature tendency from physics driver moist processes', &
         'K/s', missing_value=missing_value)

      id_tdt_phys = register_diag_field ( mod_name,            &
         'tdt_phys', axes(1:3), Time,                          &
         'temperature tendency from physics ', &
         'K/s', missing_value=missing_value)

      allocate (id_tracer_phys(ntp))
      allocate (id_tracer_phys_vdif_dn(ntp))
      allocate (id_tracer_phys_vdif_up(ntp))
      allocate (id_tracer_phys_turb(ntp))
      allocate (id_tracer_phys_moist(ntp))

      do n = 1,ntp

        call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                               units = tracer_units)
        
        diaglname = trim(tracer_name)//  &
                    ' tendency from physics'
        id_tracer_phys(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver vdif down'
        id_tracer_phys_vdif_dn(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_vdif_dn',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver vdif up'
        id_tracer_phys_vdif_up(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_vdif_up',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver vert turb'
        id_tracer_phys_turb(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_turb',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver moist processes'
        id_tracer_phys_moist(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_moist',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

      end do

      call monin_obukhov_init
!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!-----------------------------------------------------------------------

 end subroutine physics_driver_init


!######################################################################
! <SUBROUTINE NAME="physics_driver_down_time_vary">
!  <OVERVIEW>
!    physics_driver_time_vary makes sure that all time-dependent, spacially-
!    independent calculations are completed before entering window or thread
!    loops. Resultant fields are usually saved as module variables in the
!    module where needed.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_time_vary makes sure that all time-dependent, spacially-
!    independent calculations are completed before entering window or thread
!    loops. Resultant fields are usually saved as module variables in the
!    module where needed.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_down_time_vary (Time, Time_next)
!
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   time of next time step
!  </IN>
! </SUBROUTINE>
!

subroutine physics_driver_down_time_vary (Time, Time_next, dt)

!---------------------------------------------------------------------
!    physics_driver_down_time_vary makes sure that all time-dependent, 
!    spacially-independent calculations are completed before entering window 
!    or thread loops. Resultant fields are usually saved as module variables in 
!    the module where needed.
!-----------------------------------------------------------------------

type(time_type),         intent(in)             :: Time, Time_next
real,                    intent(in)             :: dt

type(time_type) :: Time_last
!---------------------------------------------------------------------      
      call damping_driver_time_vary (dt)
      call atmos_tracer_driver_time_vary (Time)

!-------------------------------------------------------------------------      

end subroutine physics_driver_down_time_vary



!######################################################################

subroutine physics_driver_down_endts(is,js)

integer, intent(in)  :: is,js

      call damping_driver_endts
      call atmos_tracer_driver_endts

!--------------------------------------------------------------------
!    set a flag to indicate that this check was done and need not be
!    done again.
!--------------------------------------------------------------------
      do_check_args = .false.

end subroutine physics_driver_down_endts

!###################################################################

subroutine physics_driver_up_time_vary (Time, Time_next, dt)

!---------------------------------------------------------------------
!    physics_driver_up_time_vary makes sure that all time-dependent, 
!    spacially-independent calculations are completed before entering 
!    window or thread loops. Resultant fields are usually saved as 
!    module variables in the module where needed.
!-----------------------------------------------------------------------

type(time_type),         intent(in)             :: Time
type(time_type),         intent(in)             :: Time_next
real,                    intent(in)             :: dt

    if (do_moist_processes) then
      call aerosol_time_vary (Time, Aerosol_cld)
      call moist_processes_time_vary (dt)
    endif
    if (do_cosp) call cosp_driver_time_vary (Time_next)

!----------------------------------------------------------------------      

end subroutine physics_driver_up_time_vary


!######################################################################

subroutine physics_driver_up_endts (is,js)

integer, intent(in)  :: is,js

    if (do_cosp) call cosp_driver_endts
    if (do_moist_processes) then
      call moist_processes_endts (is,js)
      call aerosol_endts (Aerosol_cld)
    endif

end subroutine physics_driver_up_endts


!######################################################################


!######################################################################
! <SUBROUTINE NAME="physics_driver_down">
!  <OVERVIEW>
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.    
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_down (is, ie, js, je,                       &
!                                Time_prev, Time, Time_next,           &
!                                lat, lon, area,                       &
!                                p_half, p_full, z_half, z_full,       &
!                                u, v, t, q, r, um, vm, tm, qm, rm,    &
!                                frac_land, rough_mom,                 &
!                                albedo,    t_surf_rad, t_ref, q_ref,  &
!                                u_star,    b_star, q_star,            &
!                                dtau_du,  dtau_dv,  tau_x,  tau_y,    &
!                                udt, vdt, tdt, qdt, rdt,              &
!                                flux_sw,  flux_lw,  coszen,  gust,    &
!                                Surf_diff
!  </TEMPLATE>
!  <IN NAME="Time_prev" TYPE="time_type">
!   previous time, for variable um, vm, tm, qm, rm
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   next time, used for diagnostics
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <INOUT NAME="rd" TYPE="real">
!   multiple 3d diagnostic tracer fields 
!  </INOUT>
!  <IN NAME="frac_land" TYPE="real">
!   fraction of land coverage in a model grid point
!  </IN>
!  <IN NAME="rough_mom" TYPE="real">
!   boundary layer roughness
!  </IN>
!  <IN NAME="albedo" TYPE="real">
!   surface albedo
!  </IN>
!  <IN NAME="t_surf_rad" TYPE="real">
!   surface radiative temperature
!  </IN>
!  <IN NAME="u_star" TYPE="real">
!   boundary layer wind speed (frictional speed)
!  </IN>
!  <IN NAME="b_star" TYPE="real">
!   ???
!  </IN>
!  <IN NAME="q_star" TYPE="real">
!   boundary layer specific humidity
!  </IN>
!  <IN NAME="dtau_du" TYPE="real">
!   derivative of zonal surface stress w.r.t zonal wind speed
!  </IN>
!  <IN NAME="dtau_dv" TYPE="real">
!   derivative of meridional surface stress w.r.t meridional wind speed
!  </IN>
!  <INOUT NAME="tau_x" TYPE="real">
!   boundary layer meridional component of wind shear
!  </INOUT>
!  <INOUT NAME="tau_y" TYPE="real">
!   boundary layer zonal component of wind shear
!  </INOUT>
!  <INOUT NAME="udt" TYPE="real">
!   zonal wind tendency
!  </INOUT>
!  <INOUT NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </INOUT>
!  <INOUT NAME="tdt" TYPE="real">
!   temperature tendency
!  </INOUT>
!  <INOUT NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </INOUT>
!  <INOUT NAME="rdt" TYPE="real">
!   multiple tracer tendencies
!  </INOUT>
!  <OUT NAME="flux_sw" TYPE="real">
!   Shortwave flux from radiation package
!  </OUT>
!  <OUT NAME="flux_lw" TYPE="real">
!   Longwave flux from radiation package
!  </OUT>
!  <OUT NAME="coszen" TYPE="real">
!   cosine of zenith angle
!  </OUT>
!  <OUT NAME="gust" TYPE="real">
!  </OUT>
!  <INOUT NAME="Surf_diff" TYPE="surface_diffusion_type">
!   Surface diffusion 
!  </INOUT>
!
! </SUBROUTINE>
!
subroutine physics_driver_down (is, ie, js, je, npz,              &
                                Time_prev, Time, Time_next,       &
                                lat, lon, area,                   &
                                Physics_input_block,              &
                                frac_land, rough_mom,             &
                                frac_open_sea,                    &
                                albedo,                           &
                                t_surf_rad, t_ref, q_ref,         &
                                u_star,    b_star, q_star,        &
                                dtau_du, dtau_dv,  tau_x,  tau_y, &
                                Physics_tendency_block,           &
                                Surf_diff,                        &
                                gust,                             &
                                Rad_flux_control,                 &
                                Rad_flux_block,                   &
                                diffm, difft  )

!---------------------------------------------------------------------
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.
!-----------------------------------------------------------------------

integer,                 intent(in)             :: is, ie, js, je, npz
type(time_type),         intent(in)             :: Time_prev, Time, Time_next
real,dimension(:,:),     intent(in)             :: lat, lon, area
type(physics_input_block_type), intent(in)      :: Physics_input_block
real,dimension(:,:),     intent(in)             :: frac_land,   &
                                                   rough_mom, &
                                                   albedo, t_surf_rad, &
                                                   t_ref, q_ref, &  ! cjg: PBL depth mods
                                                   u_star, b_star,    &
                                                   q_star, dtau_du,   &
                                                   dtau_dv, frac_open_sea
real,dimension(:,:),     intent(inout)          :: tau_x,  tau_y
type(physics_tendency_block_type), intent(inout):: Physics_tendency_block
real,dimension(:,:),     intent(out)            :: gust
type(surf_diff_type),    intent(inout)          :: Surf_diff
type(radiation_flux_control_type),  intent(in)  :: Rad_flux_control
type(radiation_flux_block_type),    intent(in)  :: Rad_flux_block
real,  dimension(:,:,:), intent(out)  ,optional :: diffm, difft 

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Time_prev      previous time, for variables um,vm,tm,qm,rm 
!                     (time_type)
!      Time           current time, for variables u,v,t,q,r  (time_type)
!      Time_next      next time, used for diagnostics   (time_type)
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous 
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      frac_land
!      rough_mom
!      albedo
!      albedo_vis_dir surface visible direct albedo [ dimensionless ]
!      albedo_nir_dir surface nir direct albedo [ dimensionless ]
!      albedo_vis_dif surface visible diffuse albedo [ dimensionless ]
!      albedo_nir_dif surface nir diffuse albedo [ dimensionless ]
!      t_surf_rad
!      u_star
!      b_star
!      q_star
!      dtau_du
!      dtau_dv
!
!  intent(inout) variables:
!
!      tau_x
!      tau_y
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency 
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!      rd             multiple 3d diagnostic tracer fields 
!                     [ unit / unit / sec ]
!      Surf_diff      surface_diffusion_type variable
!
!   intent(out) variables:
!
!      flux_sw
!      flux_sw_dir            net shortwave surface flux (down-up) [ w / m^2 ]
!      flux_sw_dif            net shortwave surface flux (down-up) [ w / m^2 ]
!      flux_sw_down_vis_dir   downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_down_vis_dif   downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_down_total_dir total downward shortwave surface flux [ w / m^2 ]
!      flux_sw_down_total_dif total downward shortwave surface flux [ w / m^2 ]
!      flux_sw_vis            net downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_vis_dir        net downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_sw_vis_dif        net downward shortwave surface flux in visible spectrum [ w / m^2 ]
!      flux_lw
!      coszen
!      gust
!
!   intent(in), optional variables:
!
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    local variables:

      real, dimension(ie-is+1,je-js+1,npz) :: diff_t_vert, &
                                                        diff_m_vert
      real, dimension(ie-is+1,je-js+1,npz) :: tdt_rad, &
                                                        tdt_lw
      real, dimension(ie-is+1,je-js+1) :: z_pbl
      integer          ::    sec, day, n, nextinct
      real             ::    dt, alpha, dt2
      logical          ::    used

!---> h1g, 2015-08-11
      real, dimension(ie-is+1,je-js+1) :: tke_avg
!<--- h1g, 2015-08-11

!---------------------------------------------------------------------
!   local variables:
!
!      diff_t_vert     vertical diffusion coefficient for temperature
!                      calculated on the current step
!      diff_m_vert     vertical diffusion coefficient for momentum   
!                      calculated on the current step
!      z_pbl           height of planetary boundary layer
!      sec, day        second and day components of the time_type 
!                      variable
!      dt              model physics time step [ seconds ]
!      alpha           ratio of physics time step to diffusion-smoothing
!                      time scale
!
!---------------------------------------------------------------------
      real, dimension(:,:,:,:), pointer :: r, rm
      real, dimension(:,:,:), pointer :: p_full, p_half, z_full, z_half
      real, dimension(:,:,:), pointer :: udt, vdt, tdt
      real, dimension(:,:,:,:), pointer :: rdt, rdiag
      real, dimension(:,:,:), pointer :: u, v, t, um, vm, tm 

      u => Physics_input_block%u
      v => Physics_input_block%v
      t => Physics_input_block%t
      r => Physics_input_block%q
      if (associated(Physics_input_block%um)) then
        um => Physics_input_block%um
        vm => Physics_input_block%vm
        tm => Physics_input_block%tm
        rm => Physics_input_block%qm
      else
        um => Physics_input_block%u
        vm => Physics_input_block%v
        tm => Physics_input_block%t
        rm => Physics_input_block%q
      endif
      p_full => Physics_input_block%p_full
      p_half => Physics_input_block%p_half
      z_full => Physics_input_block%z_full
      z_half => Physics_input_block%z_half
      udt => Physics_tendency_block%u_dt
      vdt => Physics_tendency_block%v_dt
      tdt => Physics_tendency_block%t_dt
      rdt => Physics_tendency_block%q_dt
      rdiag => Physics_tendency_block%qdiag

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
                         'module has not been initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    if COSP is activated, save the surface (skin) temperature for
!    its use.
!---------------------------------------------------------------------
      if (do_cosp) then
        tsurf_save(is:ie,js:je) = t_surf_rad
      endif

!---------------------------------------------------------------------
!    check the size of the input arguments. this is only done on the
!    first call to physics_driver_down.
!---------------------------------------------------------------------
      if (do_check_args) call check_args  &
                   (lat, lon, area, p_half, p_full, z_half, z_full, &
                    u, v, t, r(:,:,:,1), r, um, vm, tm, r(:,:,:,1), rm, &
                    udt, vdt, tdt, rdt(:,:,:,1), rdt, rdiag)

!---------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
      call get_time (Time_next - Time_prev, sec, day)
      dt = real(sec + day*86400)

!rab      if(do_grey_radiation) then !rif:(09/10/09) 
!rab        call grey_radiation(is, js, Time, Time_next, lat, lon, phalfgrey, albedo, t_surf_rad, t, tdt, flux_sw, flux_lw)
!rab        coszen = 1.0
!rab        flux_sw_dir     = R1*flux_sw
!rab        flux_sw_dif     = R2*flux_sw
!rab        flux_sw_vis_dir = R3*flux_sw
!rab        flux_sw_vis_dif = R4*flux_sw
!rab      endif

      if (do_radiation) then
        radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) + Rad_flux_block%tdt_rad(:,:,:)
	surf_diff%tdt_rad(is:ie,js:je,:)=Rad_flux_block%tdt_rad(:,:,:) !miz
      endif
#ifdef SCM
! Option to add SCM radiative tendencies from forcing to Rad_flux_block%tdt_lw
! and radturbten

      if (use_scm_rad) then
        call add_scm_tdtlw( Rad_flux_block%tdt_lw )
        call add_scm_tdtlw( radturbten (is:ie,js:je,:) )
        call add_scm_tdtsw( radturbten (is:ie,js:je,:) )
      endif

#endif

!----------------------------------------------------------------------
!    call damping_driver to calculate the various model dampings that
!    are desired. 
!----------------------------------------------------------------------
      z_pbl(:,:) = pbltop(is:ie,js:je) 
      call mpp_clock_begin ( damping_clock )
      call damping_driver (is, js, lat, Time_next, dt, area,        &
                           p_full, p_half, z_full, z_half,          &
                           um, vm, tm, rm(:,:,:,1), rm(:,:,:,1:ntp),&
                           udt, vdt, tdt, rdt(:,:,:,1), rdt, z_pbl)
     call mpp_clock_end ( damping_clock )

!---------------------------------------------------------------------
!    call vert_turb_driver to calculate diffusion coefficients. save
!    the planetary boundary layer height on return.
!---------------------------------------------------------------------

      if (id_tdt_phys_turb > 0) then
        used = send_data ( id_tdt_phys_turb, -2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_turb(n) > 0) then
          used = send_data ( id_tracer_phys_turb(n), -2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

      call mpp_clock_begin ( turb_clock )
      call vert_turb_driver (is, js, Time, Time_next, dt,            &
                             Rad_flux_block%tdt_lw, frac_land,  &
                             p_half, p_full, z_half, z_full,         &
                             t_ref, q_ref,                           &  ! cjg: PBL depth mods
                             u_star, b_star, q_star, rough_mom,      &
                             lat, convect(is:ie,js:je),              &
                             u, v, t, r(:,:,:,1), r, um, vm,                  &
                             tm, rm(:,:,:,1), rm, rdiag,                      &
                             udt, vdt, tdt, rdt(:,:,:,1), rdt,                &
                             diff_t_vert, diff_m_vert, gust, z_pbl, tke_avg = tke_avg)
     call mpp_clock_end ( turb_clock )
     pbltop(is:ie,js:je) = z_pbl(:,:)
     tke   (is:ie,js:je) = tke_avg(:,:)

      if (id_tdt_phys_turb > 0) then
        used = send_data ( id_tdt_phys_turb, +2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_turb(n) > 0) then
          used = send_data ( id_tracer_phys_turb(n), +2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

!-----------------------------------------------------------------------
!    process any tracer fields.
!-----------------------------------------------------------------------

      nextinct = get_tracer_index(MODEL_ATMOS,'Extinction')
      if (Rad_flux_control%do_rad .and. nextinct /= NO_TRACER) then
        rdiag(:,:,:,nextinct) = Rad_flux_block%extinction(:,:,:)
      endif

      call mpp_clock_begin ( tracer_clock )
      call atmos_tracer_driver (is, ie, js, je, Time, lon, lat,  &
                                area, z_pbl, rough_mom,         &
                                frac_open_sea, frac_land, &
                                p_half, p_full,  &
                                u, v, t, r(:,:,:,1), r, rm, rdt, rdiag, dt, &
                                u_star, b_star, q_star, z_half, z_full, &
                                t_surf_rad, albedo, Time_next, &
                                Rad_flux_block%flux_sw_down_vis_dir, &
                                Rad_flux_block%flux_sw_down_vis_dif)
      call mpp_clock_end ( tracer_clock )

!-----------------------------------------------------------------------
!    optionally use an implicit calculation of the vertical diffusion 
!    coefficients.
!
!    the vertical diffusion coefficients are solved using an implicit
!    solution to the following equation:
!
!    dK/dt   = - ( K - K_cur) / tau_diff
!
!    where K         = diffusion coefficient
!          K_cur     = diffusion coefficient diagnosed from current 
!                      time steps' state
!          tau_diff  = time scale for adjustment
!
!    in the code below alpha = dt / tau_diff
!---------------------------------------------------------------------
      if (diffusion_smooth) then
        call get_time (Time_next - Time, sec, day)
        dt2 = real(sec + day*86400)
        alpha = dt2/tau_diff
        diff_m(is:ie,js:je,:) = (diff_m(is:ie,js:je,:) +       &
                                 alpha*(diff_m_vert(:,:,:) +  &
                                 diff_cu_mo(is:ie,js:je,:)) )/&
                                 (1. + alpha)
        where (diff_m(is:ie,js:je,:) < diff_min)
          diff_m(is:ie,js:je,:) = 0.0
        end where
        diff_t(is:ie,js:je,:) = (diff_t(is:ie,js:je,:) +      &
                                 alpha*diff_t_vert(:,:,:) )/  &
                                 (1. + alpha)
        where (diff_t(is:ie,js:je,:) < diff_min)
          diff_t(is:ie,js:je,:) = 0.0
        end where
      else
        diff_t(is:ie,js:je,:) = diff_t_vert
        diff_m(is:ie,js:je,:) = diff_m_vert + diff_cu_mo(is:ie, js:je,:)
      end if

!-----------------------------------------------------------------------
!    call vert_diff_driver_down to calculate the first pass atmos-
!    pheric vertical diffusion.
!-----------------------------------------------------------------------

      if (id_tdt_phys_vdif_dn > 0) then
        used = send_data ( id_tdt_phys_vdif_dn, -2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_dn(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_dn(n), -2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

      call mpp_clock_begin ( diff_down_clock )
      radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) - tdt(:,:,:)
      if (do_clubb > 0) then
        call vert_diff_driver_down (is, js, Time_next, dt, p_half,   &
                                    p_full, z_full,   &
                                    diff_m(is:ie,js:je,:),         &
                                    diff_t(is:ie,js:je,:),         &
                                    u ,v ,t ,r(:,:,:,1) ,r(:,:,:,1:ntp), &
                                    dtau_du, dtau_dv, tau_x, tau_y,  &
                                    udt, vdt, tdt, rdt(:,:,:,1), rdt,       &
                                    Surf_diff,                     &
                                    diff_t_clubb=diff_t_clubb(is:ie,js:je,:))
      else
        call vert_diff_driver_down (is, js, Time_next, dt, p_half,   &
                                    p_full, z_full,   &
                                    diff_m(is:ie,js:je,:),         &
                                    diff_t(is:ie,js:je,:),         &
                                    u ,v ,t ,r(:,:,:,1) ,r(:,:,:,1:ntp), &
                                    dtau_du, dtau_dv, tau_x, tau_y,  &
                                    udt, vdt, tdt, rdt(:,:,:,1), rdt,        &
                                    Surf_diff)
      endif

      if (id_tdt_phys_vdif_dn > 0) then
        used = send_data ( id_tdt_phys_vdif_dn, +2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_dn(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_dn(n), +2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

!---------------------------------------------------------------------
!    if desired, return diff_m and diff_t to calling routine.
!-----------------------------------------------------------------------
      if (present(difft)) then
        difft = diff_t(is:ie,js:je,:)
      endif
      if (present(diffm)) then
        diffm = diff_m(is:ie,js:je,:)
      endif

     call mpp_clock_end ( diff_down_clock )

      u => null()
      v => null()
      t => null()
      r => null()
      um => null()
      vm => null()
      tm => null()
      rm => null()
      p_full => null()
      p_half => null()
      z_full => null()
      z_half => null()
      udt => null()
      vdt => null()
      tdt => null()
      rdt => null()
      rdiag => null()

 end subroutine physics_driver_down



!#######################################################################
! <SUBROUTINE NAME="physics_driver_up">
!  <OVERVIEW>
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_up (is, ie, js, je,                    &
!                               Time_prev, Time, Time_next,        &
!                               lat, lon, area,                    &
!                               p_half, p_full, z_half, z_full,    & 
!                               omega,                             &
!                               u, v, t, q, r, um, vm, tm, qm, rm, &
!                               frac_land,                         &
!                               udt, vdt, tdt, qdt, rdt,           &
!                               Surf_diff,                         &
!                               lprec,   fprec, gust  )            &
!  </TEMPLATE>
!  <IN NAME="Time_prev" TYPE="time_type">
!   previous time, for variable um, vm, tm, qm, rm
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   next time, used for diagnostics
!  </IN>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="omega" TYPE="real">
!   Veritical pressure tendency
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <IN NAME="frac_land" TYPE="real">
!   fraction of land coverage in a model grid point
!  </IN>
!  <INOUT NAME="udt" TYPE="real">
!   zonal wind tendency
!  </INOUT>
!  <INOUT NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </INOUT>
!  <INOUT NAME="tdt" TYPE="real">
!   temperature tendency
!  </INOUT>
!  <INOUT NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </INOUT>
!  <INOUT NAME="rdt" TYPE="real">
!   multiple tracer tendencies
!  </INOUT>
!  <OUT NAME="lprec" TYPE="real">
!  </OUT>
!  <OUT NAME="fprec" TYPE="real">
!  </OUT>
!  <OUT NAME="gust" TYPE="real">
!  </OUT>
!  <INOUT NAME="Surf_diff" TYPE="surface_diffusion_type">
!   Surface diffusion 
!  </INOUT>
! </SUBROUTINE>
!
 subroutine physics_driver_up (is, ie, js, je, npz,        &
                               Time_prev, Time, Time_next, &
                               lat, lon, area,             &
                               Physics_control,            &
                               Physics_input_block,        &
                               frac_land,                  &
                               u_star, b_star, q_star,     &
                               shflx, lhflx,               &!miz
                               Physics_tendency_block,     &
                               Moist_clouds_block,         &
                               Cosp_control, Cosp_block,   &
                               Exch_ctrl,  Surf_diff,      &
                               lprec, fprec, gust)

!----------------------------------------------------------------------
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!---------------------------------------------------------------------

integer,                intent(in)                :: is, ie, js, je, npz
type(time_type),        intent(in)                :: Time_prev, Time, Time_next
real,dimension(:,:),    intent(in)                :: lat, lon, area
type(physics_control_type), intent(in)            :: Physics_control
type(physics_input_block_type), intent(inout)     :: Physics_input_block
real,dimension(:,:),    intent(in)                :: frac_land
real,dimension(:,:),    intent(in)                :: u_star, b_star, q_star, shflx, lhflx!miz
type(physics_tendency_block_type), intent(inout)  :: Physics_tendency_block
type(clouds_from_moist_block_type), intent(inout) :: Moist_clouds_block
type(cosp_from_rad_control_type),   intent(in)    :: Cosp_control
type(cosp_from_rad_block_type),     intent(in)    :: Cosp_block
type (exchange_control_type), intent(in)          :: Exch_ctrl
type(surf_diff_type),   intent(inout)             :: Surf_diff
real,dimension(:,:),    intent(out)               :: lprec, fprec
real,dimension(:,:),    intent(inout)             :: gust

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      Time_prev      previous time, for variables um,vm,tm,qm,rm 
!                     (time_type)
!      Time           current time, for variables u,v,t,q,r  (time_type)
!      Time_next      next time, used for diagnostics   (time_type)
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      omega
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous 
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      frac_land
!      rough_mom
!      albedo
!      t_surf_rad
!      u_star
!      b_star
!      q_star
!      dtau_du
!      dtau_dv
!
!  intent(inout) variables:
!
!      tau_x
!      tau_y
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency 
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!      Surf_diff      surface_diffusion_type variable
!      gust
!
!   intent(out) variables:
!
!      lprec     
!      fprec       
!
!   intent(in), optional variables:
!
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!   local variables:

      real, dimension(ie-is+1, je-js+1, npz) :: diff_cu_mo_loc
      real, dimension(ie-is+1, je-js+1)            :: gust_cv
      real, dimension(ie-is+1, je-js+1)            :: land_mask
      integer :: sec, day
      real    :: dt
      real, dimension(ie-is+1, je-js+1) :: u_sfc, v_sfc
      real, dimension(ie-is+1, je-js+1, npz+1) :: pflux
      real, dimension(ie-is+1, je-js+1, npz)   ::  &
                             tca, cca, rhoi, lsliq, lsice, ccliq,  &
                             ccice, reff_lsclliq, reff_lsclice, &
                             reff_ccclliq, reff_ccclice, &
                             reff_lsprliq, reff_lsprice, &
                             reff_ccprliq, reff_ccprice, &
                             fl_lsrain_loc, fl_lssnow_loc,  &
                             fl_lsgrpl_loc, &
                             fl_donmca_rain_loc, fl_donmca_snow_loc, &
                             fl_ccrain_loc, fl_ccsnow_loc, mr_ozone_loc
      real, dimension(ie-is+1, je-js+1, npz, ncol) ::  &
                             stoch_mr_liq, stoch_mr_ice, &
                             stoch_size_liq, stoch_size_frz
      type(aerosol_type) :: Aerosol
      integer :: i, j , k, n
      integer :: nls, ncc
      real    :: alphb
      integer :: flag_ls, flag_cc
      integer :: kmax
      logical :: used
      logical :: hydrostatic, phys_hydrostatic
      integer :: istrat, icell, imeso, ishallow

! save the temperature and moisture tendencies from sensible and latent heat fluxes
      real, dimension(ie-is+1, je-js+1) :: tdt_shf,  qdt_lhf
   
!---------------------------------------------------------------------
!   local variables:
!
!        diff_cu_mo_loc   diffusion coefficient contribution due to 
!                         cumulus momentum transport
!        gust_cv
!        sec, day         second and day components of the time_type 
!                         variable
!        dt               physics time step [ seconds ]
!      Aerosol         aerosol_type variable describing the aerosol
!                      fields to be seen by the moist processes routines
!
!---------------------------------------------------------------------
      real, dimension(:,:,:),   pointer :: u, v, w, t, um, vm, tm, omega
      real, dimension(:,:,:,:), pointer :: r, rm
      real, dimension(:,:,:),   pointer :: p_full, p_half, z_full, z_half
      real, dimension(:,:,:),   pointer :: udt, vdt, tdt
      real, dimension(:,:,:,:), pointer :: rdt, rdiag

      u => Physics_input_block%u
      v => Physics_input_block%v
      w => Physics_input_block%w
      t => Physics_input_block%t
      r => Physics_input_block%q
      if (associated(Physics_input_block%um)) then
        um => Physics_input_block%um
        vm => Physics_input_block%vm
        tm => Physics_input_block%tm
        rm => Physics_input_block%qm
      else
        um => Physics_input_block%u
        vm => Physics_input_block%v
        tm => Physics_input_block%t
        rm => Physics_input_block%q
      endif
      omega => Physics_input_block%omega
      p_full => Physics_input_block%p_full
      p_half => Physics_input_block%p_half
      z_full => Physics_input_block%z_full
      z_half => Physics_input_block%z_half
      udt => Physics_tendency_block%u_dt
      vdt => Physics_tendency_block%v_dt
      tdt => Physics_tendency_block%t_dt
      rdt => Physics_tendency_block%q_dt
      rdiag => Physics_tendency_block%qdiag
      hydrostatic = Physics_control%hydrostatic
      phys_hydrostatic = Physics_control%phys_hydrostatic

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
             'module has not been initialized', FATAL)
      endif

!----------------------------------------------------------------------
!    define number of model layers.
!----------------------------------------------------------------------
      kmax = size(u,3)

!----------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
      call get_time (Time_next-Time_prev, sec, day)
      dt = real(sec+day*86400)

!------------------------------------------------------------------
!    call vert_diff_driver_up to complete the vertical diffusion
!    calculation.
!------------------------------------------------------------------

      if (id_tdt_phys_vdif_up > 0) then
        used = send_data ( id_tdt_phys_vdif_up, -2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_up(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_up(n), -2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

! ---> h1g, 2012-08-28, save temperature and moisture tendencies due to surface fluxes at lowest-level
      if( .not. l_host_applies_sfc_fluxes ) then
          tdt_shf(:,:) = tdt(:, :, kmax)
          qdt_lhf(:,:) = rdt(:, :, kmax, 1)
      endif

      call mpp_clock_begin ( diff_up_clock )
      call vert_diff_driver_up (is, js, Time_next, dt, p_half,   &
                                Surf_diff, tdt, rdt(:,:,:,1), rdt )

! ---> h1g, 2012-08-28, save temperature and moisture tendencies due to surface fluxes at lowest-level
      if( .not. l_host_applies_sfc_fluxes ) then
          tdt_shf(:,:) = tdt(:, :, kmax) - tdt_shf(:,:)
          qdt_lhf(:,:) = rdt(:, :, kmax, 1) - qdt_lhf(:,:)

          tdt(:, :, kmax) = tdt(:, :, kmax) - tdt_shf(:,:)
          rdt(:, :, kmax, 1) = rdt(:, :, kmax, 1) - qdt_lhf(:,:)
      endif

      radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) + tdt(:,:,:)
      call mpp_clock_end ( diff_up_clock )

      if (id_tdt_phys_vdif_up > 0) then
        used = send_data ( id_tdt_phys_vdif_up, +2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_up(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_up(n), +2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

!-----------------------------------------------------------------------
!    if the fms integration path is being followed, call moist processes
!    to compute moist physics, including convection and processes 
!    involving condenstion.
!-----------------------------------------------------------------------
      if (do_moist_processes) then

        if (id_tdt_phys_moist > 0) then
          used = send_data ( id_tdt_phys_moist, -2.0*tdt(:,:,:), &
                             Time_next, is, js, 1)
        endif

        do n=1,ntp
          if (id_tracer_phys_moist(n) > 0) then
            used = send_data ( id_tracer_phys_moist(n), -2.0*rdt(:,:,:,n), &
                               Time_next, is, js, 1)
          endif
        end do

        call mpp_clock_begin ( moist_processes_clock )

!-----------------------------------------------------------------------
! to avoid a call to Aerosol when using do_grey_radiation (rif, 09/02/09)
        if (.not. do_grey_radiation .and. do_moist_processes) then
        ! get aerosol mass concentrations
          pflux(:,:,1) = 0.0e+00
          do i=2,size(p_full,3)
            pflux(:,:,i) = 0.5e+00*(p_full(:,:,i-1) + p_full(:,:,i))
          end do
          pflux(:,:,size(p_full,3)+1) = p_full(:,:,size(p_full,3)) 
          call aerosol_driver (is, js, Time, r, p_half, pflux, &
                               Aerosol_cld, Aerosol, override_aerosols_cloud)
        end if
!--------------------------------------------------------------------
!    on steps on which the cosp simulator is called, move the values
!    of precip flux saved on the previous step so they will not be 
!    overwritten on the upcoming call to moist_processes.
!--------------------------------------------------------------------
        if (do_cosp) then
          if (Cosp_control%step_to_call_cosp) then
            fl_lsrain_loc(:,:,:) = fl_lsrain(is:ie,js:je,:)
            fl_lssnow_loc(:,:,:) = fl_lssnow(is:ie,js:je,:)
            fl_lsgrpl_loc(:,:,:) = fl_lsgrpl(is:ie,js:je,:)
            fl_ccrain_loc(:,:,:) = fl_ccrain(is:ie,js:je,:)
            fl_ccsnow_loc(:,:,:) = fl_ccsnow(is:ie,js:je,:)
            fl_donmca_rain_loc(:,:,:) = fl_donmca_rain(is:ie,js:je,:)
            fl_donmca_snow_loc(:,:,:) = fl_donmca_snow(is:ie,js:je,:)
            mr_ozone_loc(:,:,:)  = Cosp_block%mr_ozone(:,:,:)
          endif
        endif

       ! NOTE: moist_processes is always called with strat
       ! do_strat_will always have to be activated
       ! (or at least Cloud_data allocated for it)
       istrat   = Moist_clouds_block%index_strat
       icell    = Moist_clouds_block%index_donner_cell
       imeso    = Moist_clouds_block%index_donner_meso
       ishallow = Moist_clouds_block%index_uw_conv

       if (doing_donner .and. doing_uw_conv) then
         call moist_processes (is, ie, js, je, Time_next, dt, &
           frac_land, p_half, p_full, z_half, z_full, omega,    &
           diff_t(is:ie,js:je,:), radturbten(is:ie,js:je,:),    &
           surf_diff%tdt_rad(is:ie,js:je,:), &
           surf_diff%tdt_dyn(is:ie,js:je,:), &
           surf_diff%qdt_dyn(is:ie,js:je,:), &
           surf_diff%dgz_dyn(is:ie,js:je,:), &
           surf_diff%ddp_dyn(is:ie,js:je,:), &
	   hmint(is:ie,js:je),               &
           cush(is:ie,js:je), cbmf(is:ie,js:je), cgust(is:ie,js:je),  &!miz
	   tke(is:ie,js:je),  pblhto(is:ie,js:je), rkmo(is:ie,js:je), &!miz
	   taudpo(is:ie,js:je),                                       &!miz
           exist_shconv(is:ie,js:je,:), exist_dpconv(is:ie,js:je,:),  &!miz
           pbltop(is:ie,js:je), u_star, b_star, q_star, shflx, lhflx, &!miz
           t,   r(:,:,:,1),   r,    u,  v,  w,             &
           tm,  rm(:,:,:,1),  rm,   um, vm,                &
           tdt, rdt(:,:,:,1), rdt, rdiag, udt, vdt,    &
           diff_cu_mo_loc, convect(is:ie,js:je), lprec, &
           fprec, fl_lsrain(is:ie,js:je,:), fl_lssnow(is:ie,js:je,:),  &
           fl_ccrain(is:ie,js:je,:), fl_ccsnow(is:ie,js:je,:),    &
           fl_donmca_rain(is:ie,js:je,:), fl_donmca_snow(is:ie,js:je,:), &
           gust_cv, area, lon, lat, &
           Moist_clouds_block%Cloud_data(istrat)%cloud_area,  &
           Moist_clouds_block%Cloud_data(istrat)%liquid_amt,  &
           Moist_clouds_block%Cloud_data(istrat)%ice_amt,     &
           Moist_clouds_block%Cloud_data(istrat)%droplet_number,   &
           Moist_clouds_block%Cloud_data(istrat)%ice_number,  &
      ! snow, rain
           Moist_clouds_block%Cloud_data(istrat)%snow,  &
           Moist_clouds_block%Cloud_data(istrat)%rain,  & 
           Moist_clouds_block%Cloud_data(istrat)%snow_size,  &
           Moist_clouds_block%Cloud_data(istrat)%rain_size,  &
           diff_t_clubb =diff_t_clubb,                       &
           tdt_shf = tdt_shf,                                  &
           qdt_lhf = qdt_lhf,                                  &
           Aerosol=Aerosol, &
           shallow_cloud_area    =Moist_clouds_block%Cloud_data(ishallow)%cloud_area, &
           shallow_liquid        =Moist_clouds_block%Cloud_data(ishallow)%liquid_amt, &
           shallow_ice           =Moist_clouds_block%Cloud_data(ishallow)%ice_amt,    &
           shallow_droplet_number=Moist_clouds_block%Cloud_data(ishallow)%droplet_number, &
           shallow_ice_number    =Moist_clouds_block%Cloud_data(ishallow)%ice_number, &
           cell_cld_frac      =Moist_clouds_block%Cloud_data(icell)%cloud_area,  &
           cell_liq_amt       =Moist_clouds_block%Cloud_data(icell)%liquid_amt,  &
           cell_liq_size      =Moist_clouds_block%Cloud_data(icell)%liquid_size, &
           cell_ice_amt       =Moist_clouds_block%Cloud_data(icell)%ice_amt,     &
           cell_ice_size      =Moist_clouds_block%Cloud_data(icell)%ice_size,    &
           cell_droplet_number=Moist_clouds_block%Cloud_data(icell)%droplet_number, &
           meso_cld_frac      =Moist_clouds_block%Cloud_data(imeso)%cloud_area,  &
           meso_liq_amt       =Moist_clouds_block%Cloud_data(imeso)%liquid_amt,  &
           meso_liq_size      =Moist_clouds_block%Cloud_data(imeso)%liquid_size, &
           meso_ice_amt       =Moist_clouds_block%Cloud_data(imeso)%ice_amt,     &
           meso_ice_size      =Moist_clouds_block%Cloud_data(imeso)%ice_size,    &
           meso_droplet_number=Moist_clouds_block%Cloud_data(imeso)%droplet_number, &
           nsum_out           =Moist_clouds_block%Cloud_data(imeso)%nsum_out,   &
           hydrostatic=hydrostatic, phys_hydrostatic=phys_hydrostatic  )
       else if (doing_donner) then
        call moist_processes (is, ie, js, je, Time_next, dt, frac_land, &
                           p_half, p_full, z_half, z_full, omega,    &
                           diff_t(is:ie,js:je,:),                    &
                           radturbten(is:ie,js:je,:),                &
                           surf_diff%tdt_rad(is:ie,js:je,:), & !miz
			   surf_diff%tdt_dyn(is:ie,js:je,:), & !miz
			   surf_diff%qdt_dyn(is:ie,js:je,:), & !miz
			   surf_diff%dgz_dyn(is:ie,js:je,:), & !miz
			   surf_diff%ddp_dyn(is:ie,js:je,:), & !miz
			   hmint(is:ie,js:je),               & !miz
                           cush           (is:ie,js:je),             &
                           cbmf           (is:ie,js:je),             &
                           cgust          (is:ie,js:je),             &
	                   tke            (is:ie,js:je),             &
                           pblhto         (is:ie,js:je),             &!miz
                           rkmo           (is:ie,js:je),             &!miz
                           taudpo         (is:ie,js:je),             &!miz
                           exist_shconv   (is:ie,js:je,:),           &
                           exist_dpconv   (is:ie,js:je,:),           &
                            pbltop(is:ie,js:je),         &!miz
                            u_star, b_star, q_star, shflx, lhflx,    &!miz
                           t,   r(:,:,:,1),   r,    u,  v,  w,         &
                           tm,  rm(:,:,:,1),  rm,   um, vm,            &
                           tdt, rdt(:,:,:,1), rdt, rdiag, udt, vdt,    &
                           diff_cu_mo_loc , convect(is:ie,js:je), lprec, fprec, &
            fl_lsrain(is:ie,js:je,:), fl_lssnow(is:ie,js:je,:),  &
            fl_ccrain(is:ie,js:je,:), fl_ccsnow(is:ie,js:je,:),    &
            fl_donmca_rain(is:ie,js:je,:), fl_donmca_snow(is:ie,js:je,:), &
                           gust_cv, area, lon, lat,  &
           Moist_clouds_block%Cloud_data(istrat)%cloud_area,  &
           Moist_clouds_block%Cloud_data(istrat)%liquid_amt,  &
           Moist_clouds_block%Cloud_data(istrat)%ice_amt,     &
           Moist_clouds_block%Cloud_data(istrat)%droplet_number,   &
           Moist_clouds_block%Cloud_data(istrat)%ice_number,  &
      ! snow, rain
           Moist_clouds_block%Cloud_data(istrat)%snow,  &
           Moist_clouds_block%Cloud_data(istrat)%rain,  & 
           Moist_clouds_block%Cloud_data(istrat)%snow_size,  &
           Moist_clouds_block%Cloud_data(istrat)%rain_size,  &
                           diff_t_clubb   =diff_t_clubb,                                           &
                           tdt_shf = tdt_shf,                                                      &
                           qdt_lhf = qdt_lhf,                                                      &
                           Aerosol=Aerosol,       &
           cell_cld_frac      =Moist_clouds_block%Cloud_data(icell)%cloud_area,  &
           cell_liq_amt       =Moist_clouds_block%Cloud_data(icell)%liquid_amt,  &
           cell_liq_size      =Moist_clouds_block%Cloud_data(icell)%liquid_size, &
           cell_ice_amt       =Moist_clouds_block%Cloud_data(icell)%ice_amt,     &
           cell_ice_size      =Moist_clouds_block%Cloud_data(icell)%ice_size,    &
           cell_droplet_number=Moist_clouds_block%Cloud_data(icell)%droplet_number, &
           meso_cld_frac      =Moist_clouds_block%Cloud_data(imeso)%cloud_area,  &
           meso_liq_amt       =Moist_clouds_block%Cloud_data(imeso)%liquid_amt,  &
           meso_liq_size      =Moist_clouds_block%Cloud_data(imeso)%liquid_size, &
           meso_ice_amt       =Moist_clouds_block%Cloud_data(imeso)%ice_amt,     &
           meso_ice_size      =Moist_clouds_block%Cloud_data(imeso)%ice_size,    &
           meso_droplet_number=Moist_clouds_block%Cloud_data(imeso)%droplet_number, &
           nsum_out           =Moist_clouds_block%Cloud_data(imeso)%nsum_out,   &
           hydrostatic=hydrostatic, phys_hydrostatic=phys_hydrostatic  )
                          
       else if (doing_uw_conv) then
        call moist_processes (is, ie, js, je, Time_next, dt, frac_land,         &
                            p_half, p_full, z_half, z_full, omega,    &
                            diff_t(is:ie,js:je,:),                    &
                            radturbten(is:ie,js:je,:),                &
                            surf_diff%tdt_rad(is:ie,js:je,:),         &!miz
                            surf_diff%tdt_dyn(is:ie,js:je,:),         &!miz
                            surf_diff%qdt_dyn(is:ie,js:je,:),         &!miz
                            surf_diff%dgz_dyn(is:ie,js:je,:),         &!miz
                            surf_diff%ddp_dyn(is:ie,js:je,:),         &!miz
			    hmint(is:ie,js:je),                       &!miz
                            cush           (is:ie,js:je),             &!
                            cbmf           (is:ie,js:je),             &!
                            cgust          (is:ie,js:je),             &!
	                    tke            (is:ie,js:je),             &
                            pblhto         (is:ie,js:je),             &!miz
                            rkmo           (is:ie,js:je),             &!miz
                            taudpo         (is:ie,js:je),             &!miz
                            exist_shconv   (is:ie,js:je,:),           &!
                            exist_dpconv   (is:ie,js:je,:),           &!
                            pbltop(is:ie,js:je),         &!miz
                            u_star, b_star, q_star, shflx, lhflx,     &!miz
                            t,   r(:,:,:,1),   r,    u,  v,  w,         &
                            tm,  rm(:,:,:,1),  rm,   um, vm,            &
                            tdt, rdt(:,:,:,1), rdt, rdiag, udt, vdt,    &
                            diff_cu_mo_loc , convect(is:ie,js:je), lprec, fprec,&
               fl_lsrain(is:ie,js:je,:), fl_lssnow(is:ie,js:je,:),  &
               fl_ccrain(is:ie,js:je,:), fl_ccsnow(is:ie,js:je,:),    &
           fl_donmca_rain(is:ie,js:je,:), fl_donmca_snow(is:ie,js:je,:), &
                            gust_cv, area, lon, lat,   &
           Moist_clouds_block%Cloud_data(istrat)%cloud_area,  &
           Moist_clouds_block%Cloud_data(istrat)%liquid_amt,  &
           Moist_clouds_block%Cloud_data(istrat)%ice_amt,     &
           Moist_clouds_block%Cloud_data(istrat)%droplet_number,   &
           Moist_clouds_block%Cloud_data(istrat)%ice_number,  &
      ! snow, rain
           Moist_clouds_block%Cloud_data(istrat)%snow,  &
           Moist_clouds_block%Cloud_data(istrat)%rain,  & 
           Moist_clouds_block%Cloud_data(istrat)%snow_size,  &
           Moist_clouds_block%Cloud_data(istrat)%rain_size,  &
                           diff_t_clubb   =diff_t_clubb,                       &
                           tdt_shf = tdt_shf,                                  &
                           qdt_lhf = qdt_lhf,                                  &
                           Aerosol=Aerosol, &
           shallow_cloud_area    =Moist_clouds_block%Cloud_data(ishallow)%cloud_area, &
           shallow_liquid        =Moist_clouds_block%Cloud_data(ishallow)%liquid_amt, &
           shallow_ice           =Moist_clouds_block%Cloud_data(ishallow)%ice_amt,    &
           shallow_droplet_number=Moist_clouds_block%Cloud_data(ishallow)%droplet_number, &
           shallow_ice_number    =Moist_clouds_block%Cloud_data(ishallow)%ice_number, &
                           hydrostatic=hydrostatic, phys_hydrostatic=phys_hydrostatic  )
       else
        call moist_processes (is, ie, js, je, Time_next, dt, frac_land, &
                           p_half, p_full, z_half, z_full, omega,    &
                           diff_t(is:ie,js:je,:),                    &
                           radturbten(is:ie,js:je,:),                &
                           surf_diff%tdt_rad(is:ie,js:je,:),         &!miz
                           surf_diff%tdt_dyn(is:ie,js:je,:), 	     &!miz
                           surf_diff%qdt_dyn(is:ie,js:je,:), 	     &!miz
                           surf_diff%dgz_dyn(is:ie,js:je,:), 	     &!miz
                           surf_diff%ddp_dyn(is:ie,js:je,:), 	     &!miz
			   hmint(is:ie,js:je),                       &!miz
                           cush           (is:ie,js:je),             &!
                           cbmf           (is:ie,js:je),             &!
                           cgust          (is:ie,js:je),             &!
	                   tke            (is:ie,js:je),             &
                           pblhto         (is:ie,js:je),             &!miz
                           rkmo           (is:ie,js:je),             &!miz
                           taudpo         (is:ie,js:je),             &!miz
                           exist_shconv   (is:ie,js:je,:),           &!
                           exist_dpconv   (is:ie,js:je,:),           &!
                            pbltop(is:ie,js:je),         &!miz
                            u_star, b_star, q_star, shflx, lhflx,    &!miz
                            t,   r(:,:,:,1),   r,    u,  v,  w,         &
                            tm,  rm(:,:,:,1),  rm,   um, vm,            &
                            tdt, rdt(:,:,:,1), rdt, rdiag, udt, vdt,    &
                            diff_cu_mo_loc , convect(is:ie,js:je), lprec, fprec, &
               fl_lsrain(is:ie,js:je,:), fl_lssnow(is:ie,js:je,:),  &
               fl_ccrain(is:ie,js:je,:), fl_ccsnow(is:ie,js:je,:),    &
           fl_donmca_rain(is:ie,js:je,:), fl_donmca_snow(is:ie,js:je,:), &
                           gust_cv, area, lon, lat,   &
           Moist_clouds_block%Cloud_data(istrat)%cloud_area,  &
           Moist_clouds_block%Cloud_data(istrat)%liquid_amt,  &
           Moist_clouds_block%Cloud_data(istrat)%ice_amt,     &
           Moist_clouds_block%Cloud_data(istrat)%droplet_number,   &
           Moist_clouds_block%Cloud_data(istrat)%ice_number,  &
      ! snow, rain
           Moist_clouds_block%Cloud_data(istrat)%snow,  &
           Moist_clouds_block%Cloud_data(istrat)%rain,  & 
           Moist_clouds_block%Cloud_data(istrat)%snow_size,  &
           Moist_clouds_block%Cloud_data(istrat)%rain_size,  &
                           diff_t_clubb   =diff_t_clubb,                                           &
                           tdt_shf = tdt_shf,                                                      &
                           qdt_lhf = qdt_lhf,                                                      &
                           Aerosol=Aerosol, &
                           hydrostatic=hydrostatic, phys_hydrostatic=phys_hydrostatic  )
        endif
        call mpp_clock_end ( moist_processes_clock )
        diff_cu_mo(is:ie, js:je,:) = diff_cu_mo_loc(:,:,:)
        radturbten(is:ie,js:je,:) = 0.0

!---------------------------------------------------------------------
!    add the convective gustiness effect to that previously obtained 
!    from non-convective parameterizations.
!---------------------------------------------------------------------
        gust = sqrt( gust*gust + gust_cv*gust_cv)

        if (id_tdt_phys_moist > 0) then
          used = send_data ( id_tdt_phys_moist, +2.0*tdt(:,:,:), &
                             Time_next, is, js, 1)
        endif
        do n=1,ntp
          if (id_tracer_phys_moist(n) > 0) then
            used = send_data ( id_tracer_phys_moist(n), +2.0*rdt(:,:,:,n), &
                               Time_next, is, js, 1)
          endif
        end do
        if (id_tdt_phys > 0) then
           used = send_data ( id_tdt_phys, tdt(:,:,:), &
                              Time_next, is, js, 1)
        endif
        do n=1,ntp
          if (id_tracer_phys(n) > 0) then
            used = send_data ( id_tracer_phys(n), rdt(:,:,:,n), &
                               Time_next, is, js, 1)
          endif
        end do

      endif ! do_moist_processes

      if (do_moist_processes) then  
        call aerosol_dealloc (Aerosol)
      endif
      
      if (do_cosp) then
        call mpp_clock_begin ( cosp_clock )
        if (Cosp_control%step_to_call_cosp) then

!---------------------------------------------------------------------
!    on the first step of a job segment, the values of t,q and precip 
!    flux will not be available at the proper time level. in this case
!    denoted by temp-_last = 0.0, use values from the current step for 
!    t, q and precip flux.
!---------------------------------------------------------------------
          alphb = SUM(temp_last(is:ie,js:je,:))
          if (alphb == 0.) then
            temp_last(is:ie,js:je,:) = t(:,:,:) + dt*tdt(:,:,:)
            q_last(is:ie,js:je,:) = r(:,:,:,1) + dt*rdt(:,:,:,1)
            fl_lsrain_loc(:,:,:) = fl_lsrain(is:ie,js:je,:)
            fl_lssnow_loc(:,:,:) = fl_lssnow(is:ie,js:je,:)
            fl_lsgrpl_loc(:,:,:) = fl_lsgrpl(is:ie,js:je,:)
            fl_ccrain_loc(:,:,:) = fl_ccrain(is:ie,js:je,:)
            fl_ccsnow_loc(:,:,:) = fl_ccsnow(is:ie,js:je,:)
            fl_donmca_rain_loc(:,:,:) = fl_donmca_rain(is:ie,js:je,:)
            fl_donmca_snow_loc(:,:,:) = fl_donmca_snow(is:ie,js:je,:)
          endif

!----------------------------------------------------------------------
!    define the total and convective cloud fractions in each grid box as
!    the average over the stochastic columns.
!----------------------------------------------------------------------
          tca = 0.
          cca = 0.
          do n=1,ncol                      
            where (Cosp_block%stoch_cloud_type(:,:,:,n) > 0.) 
              tca(:,:,:)  = tca(:,:,:) +  1.0
            end where
            where (Cosp_block%stoch_cloud_type(:,:,:,n) == 2.) 
              cca(:,:,:)  = cca(:,:,:) +  1.0
            end where
          end do
          tca = tca/ float(ncol)                
          cca = cca/ float(ncol)

!--------------------------------------------------------------------
!    define the atmospheric density to use in converting concentrations
!    to mixing ratios.
!--------------------------------------------------------------------
          do k=1, size(Cosp_block%stoch_cloud_type,3)
            do j=1, size(t,2)
              do i=1, size(t,1)
                rhoi(i,j,k) =  RDGAS*temp_last(i+is-1,j+js-1,k)/ &
                                                          p_full(i,j,k) 
              end do
            end do
          end do

!--------------------------------------------------------------------
!   convert the condensate concentrations in each stochastic column to 
!   mixing ratios. 
!--------------------------------------------------------------------
          do n=1,ncol                       
            do k=1, size(Cosp_block%stoch_cloud_type,3)
              do j=1, size(t,2)
                do i=1, size(t,1)
                  stoch_mr_liq(i,j,k,n) = 1.0e-03*  &
                         Cosp_block%stoch_conc_drop(i,j,k,n)*rhoi(i,j,k)
                  stoch_mr_ice(i,j,k,n) = 1.0e-03*  &
                         Cosp_block%stoch_conc_ice (i,j,k,n)*rhoi(i,j,k)
                  stoch_size_liq(i,j,k,n) = 1.0e-06*  &
                         Cosp_block%stoch_size_drop(i,j,k,n)
                  stoch_size_frz(i,j,k,n) = 1.0e-06*  &
                         Cosp_block%stoch_size_ice (i,j,k,n)
                end do
              end do
            end do
          end do
          stoch_mr_liq = stoch_mr_liq/(1.0-stoch_mr_liq)
          stoch_mr_ice = stoch_mr_ice/(1.0-stoch_mr_ice)

!---------------------------------------------------------------------
!    define the grid box mean largescale and convective condensate 
!    mixing ratios and sizes.
!---------------------------------------------------------------------
          lsliq = 0.
          lsice = 0.
          ccliq = 0.
          ccice = 0.
          reff_lsclliq = 0.
          reff_lsclice = 0.
          reff_ccclliq = 0.
          reff_ccclice = 0.
          reff_lsprliq = 0.
          reff_lsprice = 0.
          reff_ccprliq = 0.
          reff_ccprice = 0.
          do k=1, size(Cosp_block%stoch_cloud_type,3)
            do j=1, size(t,2)
              do i=1, size(t,1)
                nls = 0
                ncc = 0
                do n=1,ncol                       
                  if (Cosp_block%stoch_cloud_type(i,j,k,n) == 1.) then
                    nls = nls + 1
                    lsliq(i,j,k) = lsliq(i,j,k) +  &
                                     Cosp_block%stoch_conc_drop(i,j,k,n)
                    lsice(i,j,k) = lsice(i,j,k) +   &
                                     Cosp_block%stoch_conc_ice (i,j,k,n)
                    reff_lsclliq(i,j,k) = reff_lsclliq(i,j,k) +  &
                                     Cosp_block%stoch_size_drop(i,j,k,n)
                    reff_lsclice(i,j,k) = reff_lsclice(i,j,k) +  &
                                     Cosp_block%stoch_size_ice (i,j,k,n)
                  else if (Cosp_block%stoch_cloud_type(i,j,k,n) == 2.) then
                    ncc = ncc + 1
                    ccliq(i,j,k) = ccliq(i,j,k) +  &
                                     Cosp_block%stoch_conc_drop(i,j,k,n)
                    ccice(i,j,k) = ccice(i,j,k) +  &
                                     Cosp_block%stoch_conc_ice (i,j,k,n)
                    reff_ccclliq(i,j,k) = reff_ccclliq(i,j,k) +  &
                                     Cosp_block%stoch_size_drop(i,j,k,n)
                    reff_ccclice(i,j,k) = reff_ccclice(i,j,k) +  &
                                     Cosp_block%stoch_size_ice (i,j,k,n)
                  endif
                end do
                if (nls > 0) then
                  lsliq(i,j,k) = 1.0e-03*lsliq(i,j,k)/float(nls)
                  lsice(i,j,k) = 1.0e-03*lsice(i,j,k)/float(nls)
                  reff_lsclliq(i,j,k) = 1.0e-06*  &
                                      reff_lsclliq (i,j,k)/float(nls)
                  reff_lsclice(i,j,k) = 1.0e-06*  &
                                      reff_lsclice (i,j,k)/float(nls)
                endif
                if (ncc > 0) then
                  ccliq(i,j,k) = 1.0e-03*ccliq(i,j,k)/float(ncc)
                  ccice(i,j,k) = 1.0e-03*ccice(i,j,k)/float(ncc)
                  reff_ccclliq(i,j,k) = 1.0e-06*  &
                                    reff_ccclliq (i,j,k) /float(ncc)
                  reff_ccclice(i,j,k) = 1.0e-06*  &
                                    reff_ccclice (i,j,k) /float(ncc)
                endif
                ccliq(i,j,k) = ccliq(i,j,k)*rhoi(i,j,k)/ &
                                                  (1.0-ccliq(i,j,k))
                ccice(i,j,k) = ccice(i,j,k)*rhoi(i,j,k)/  &
                                                  (1.0-ccice(i,j,k))
                lsliq(i,j,k) = lsliq(i,j,k)*rhoi(i,j,k)/  &
                                                  (1.0-lsliq(i,j,k))
                lsice(i,j,k) = lsice(i,j,k)*rhoi(i,j,k)/  &
                                                  (1.0-lsice(i,j,k))
              end do
            end do
          end do
     
!---------------------------------------------------------------------
!   define land_mask array. set it to 1 over land, 0 over ocean; define
!   based on frac_land > 0.5 being land.
!---------------------------------------------------------------------
          where (frac_land > 0.50)
            land_mask(:,:) =  1.0
          elsewhere
            land_mask(:,:) =  0.0
          end where

          if (allow_cosp_precip_wo_clouds) then
          else
!--------------------------------------------------------------------
!    allow ls precip only in columns containing ls cloud. allow 
!    convective precip only in columns with convective cloud,
!--------------------------------------------------------------------
          do j=1, size(t,2)
            do i=1, size(t,1)
              flag_ls = 0
              flag_cc = 0
              do k=1, size(Cosp_block%stoch_cloud_type,3)
                do n=1,ncol                        
                  if (Cosp_block%stoch_cloud_type(i,j,k,n) == 1.) then
                    flag_ls = 1
                    exit
                  else if(Cosp_block%stoch_cloud_type(i,j,k,n) == 2.) then
                    flag_cc = 1
                    exit
                  endif
                end do
                if (flag_ls == 1 .and. flag_cc == 1) exit
              end do
              if (flag_ls == 0) then
                fl_lsrain_loc(i,j,:) = 0.
                fl_lssnow_loc(i,j,:) = 0.
                fl_lsgrpl_loc(i,j,:) = 0.
              endif 
              if (flag_cc == 0) then
                fl_ccrain_loc(i,j,:) = 0.
                fl_ccsnow_loc(i,j,:) = 0.
              endif 
            end do
          end do
          endif

          if (include_donmca_in_cosp) then
            fl_ccrain_loc = fl_ccrain_loc + fl_donmca_rain_loc
            fl_ccsnow_loc = fl_ccsnow_loc + fl_donmca_snow_loc
          endif

!---------------------------------------------------------------------
!    pass in the large-scale graupel flux, lowest-level u and v wind
!    components.
!---------------------------------------------------------------------
          fl_lsgrpl = 0.
          u_sfc = u(:,:,kmax)
          v_sfc = v(:,:,kmax)

!---------------------------------------------------------------------
!    call the cosp simulator to produce the desired outputs.
!---------------------------------------------------------------------
          call cosp_driver (lat*180./ACOS(-1.0), lon*180./ACOS(-1.0),  &
                            Cosp_block%daytime, &
                            p_half, &
                            p_full, z_half,  &
                            z_full, u_sfc, v_sfc, mr_ozone_loc, &
                            temp_last(is:ie,js:je,:),  &
                            q_last(is:ie,js:je,:), tca, cca, lsliq,  &
                            lsice, ccliq, ccice,  &
                           fl_lsrain_loc, fl_lssnow_loc, fl_lsgrpl_loc,&
                            fl_ccrain_loc, fl_ccsnow_loc,&
                            0.5*reff_lsclliq, 0.5*reff_lsclice,  &
                            reff_lsprliq, reff_lsprice,  &
                            0.5*reff_ccclliq, 0.5*reff_ccclice, &
                            reff_ccprliq, reff_ccprice, &
                            tsurf_save(is:ie, js:je), land_mask, &
                            Time_next, is, js, &
                stoch_mr_liq_in =stoch_mr_liq,  &
                stoch_mr_ice_in =stoch_mr_ice,  &
                stoch_size_liq_in =0.5*stoch_size_liq, &
                stoch_size_frz_in = 0.5*stoch_size_frz,  &
                tau_stoch_in = Cosp_block%tau_stoch,&
                lwem_stoch_in = Cosp_block%lwem_stoch, &
                stoch_cloud_type_in = Cosp_block%stoch_cloud_type)
        endif ! (Cosp_control%step_to_call_cosp)
        call mpp_clock_end ( cosp_clock )
      endif ! (do_cosp)

!--------------------------------------------------------------------
! save t and q from end of step for use with next call to COSP
!--------------------------------------------------------------------
      temp_last(is:ie,js:je,:) = t(:,:,:) + tdt(:,:,:)*dt
      q_last   (is:ie,js:je,:) = r(:,:,:,1) + rdt(:,:,:,1)*dt
       
!-----------------------------------------------------------------------

      u => null()
      v => null()
      w => null()
      t => null()
      r => null()
      um => null()
      vm => null()
      tm => null()
      rm => null()
      omega => null()
      p_full => null()
      p_half => null()
      z_full => null()
      z_half => null()
      udt => null()
      vdt => null()
      tdt => null()
      rdt => null()
      rdiag => null()

 end subroutine physics_driver_up


!#######################################################################
! <SUBROUTINE NAME="physics_driver_end">
!  <OVERVIEW>
!   physics_driver_end is the destructor for physics_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_end is the destructor for physics_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_end (Time)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
! </SUBROUTINE>
!
subroutine physics_driver_end (Time, Physics, Moist_clouds, Physics_tendency, Atm_block)

!---------------------------------------------------------------------
!    physics_driver_end is the destructor for physics_driver_mod.
!---------------------------------------------------------------------

type(time_type), intent(in) :: Time
type(physics_type), intent(in) :: Physics
type(clouds_from_moist_type), intent(inout) :: Moist_clouds(:)
type(physics_tendency_type),  intent(inout) :: Physics_tendency
type(block_control_type), intent(in) :: Atm_block

!--------------------------------------------------------------------
!   intent(in) variables:
! 
!      Time      current time [ time_type(days, seconds) ]
!
!--------------------------------------------------------------------
integer :: n, nb, nc, ibs, ibe, jbs, jbe
integer :: moist_processes_term_clock, damping_term_clock, turb_term_clock, &
           diff_term_clock, aerosol_term_clock, clubb_term_clock, &
           grey_radiation_term_clock, tracer_term_clock, cosp_term_clock

      clubb_term_clock =      &
        mpp_clock_id( '   Phys_driver_term: clubb: Termination', &
                grain=CLOCK_MODULE_DRIVER )
      moist_processes_term_clock =      &
        mpp_clock_id( '   Phys_driver_term: MP: Termination', &
                grain=CLOCK_MODULE_DRIVER )
      damping_term_clock         =     &
        mpp_clock_id( '   Phys_driver_term: Damping: Termination',    &
                  grain=CLOCK_MODULE_DRIVER )
      turb_term_clock            =      &
        mpp_clock_id( '   Phys_driver_term: Vert. Turb.: Termination', &
                  grain=CLOCK_MODULE_DRIVER )
      diff_term_clock       =     &
        mpp_clock_id( '   Phys_driver_term: Vert. Diff.: Termination',   &
                 grain=CLOCK_MODULE_DRIVER )
      cosp_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: COSP: Termination', &
                       grain=CLOCK_MODULE_DRIVER )
      if (do_moist_processes) &
      aerosol_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: Aerosol: Termination', &
                       grain=CLOCK_MODULE_DRIVER )
      grey_radiation_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: Grey Radiation: Termination', &
                       grain=CLOCK_MODULE_DRIVER )
      tracer_term_clock          =      &
        mpp_clock_id( '   Phys_driver_term: Tracer: Termination',    &
                 grain=CLOCK_MODULE_DRIVER )
!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
              'module has not been initialized', FATAL)
      endif

      do nb = 1, Atm_block%nblks
        ibs = Atm_block%ibs(nb)-Atm_block%isc+1
        ibe = Atm_block%ibe(nb)-Atm_block%isc+1
        jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
        jbe = Atm_block%jbe(nb)-Atm_block%jsc+1

        do nc = 1, size(Moist_clouds(1)%block(nb)%Cloud_data,1)

          ! common to all cloud schemes
          Restart%Cloud_data(nc)%cloud_area    (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%cloud_area
          Restart%Cloud_data(nc)%liquid_amt    (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_amt
          Restart%Cloud_data(nc)%ice_amt       (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_amt
          Restart%Cloud_data(nc)%droplet_number(ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%droplet_number
          Restart%Cloud_data(nc)%scheme_name                       = Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name

          ! properties specific to large-scale/stratiform clouds
          if (trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'strat_cloud') then
            Restart%Cloud_data(nc)%ice_number(ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number
            Restart%Cloud_data(nc)%rain      (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain
            Restart%Cloud_data(nc)%snow      (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow
            Restart%Cloud_data(nc)%rain_size (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain_size
            Restart%Cloud_data(nc)%snow_size (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow_size
          endif
 
          ! properties specific to donner deep clouds (both cell and meso)
          if (trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'donner_cell' .or. &
              trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'donner_meso') then
            Restart%Cloud_data(nc)%liquid_size(ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_size
            Restart%Cloud_data(nc)%ice_size   (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_size
            Restart%Cloud_data(nc)%nsum_out   (ibs:ibe,jbs:jbe)   = Moist_clouds(1)%block(nb)%Cloud_data(nc)%nsum_out
          endif

          ! properties specific to uw shallow convective clouds
          if (trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'uw_conv') then
            Restart%Cloud_data(nc)%ice_number(ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number
          endif

        enddo
      enddo

      call physics_driver_netcdf

!--------------------------------------------------------------------
!    call the destructor routines for those modules who were initial-
!    ized from this module.
!--------------------------------------------------------------------
      call mpp_clock_begin ( turb_term_clock )
      call vert_turb_driver_end
      call mpp_clock_end ( turb_term_clock )
      call mpp_clock_begin ( diff_term_clock )
      call vert_diff_driver_end
      call mpp_clock_end ( diff_term_clock )

      if (do_clubb > 0) then
         deallocate ( diff_t_clubb )
      end if

!--------------------------------------------------------------------
!    terminate radiation routines and data
!--------------------------------------------------------------------
      if (do_moist_processes) then
        call mpp_clock_begin ( aerosol_term_clock )
        call aerosol_end (Aerosol_cld)
        call mpp_clock_end ( aerosol_term_clock )
      endif

      call mpp_clock_begin ( grey_radiation_term_clock )
      if(do_grey_radiation) call grey_radiation_end 
      call mpp_clock_end ( grey_radiation_term_clock )

      if (do_moist_processes) then  
        call mpp_clock_begin ( moist_processes_term_clock )
        call moist_processes_end (clubb_term_clock)
        call mpp_clock_end ( moist_processes_term_clock )
      endif

      call mpp_clock_begin ( tracer_term_clock )
      call atmos_tracer_driver_end
      call mpp_clock_end ( tracer_term_clock )
      call mpp_clock_begin ( damping_term_clock )
      call damping_driver_end
      call mpp_clock_end ( damping_term_clock )
      if (do_cosp) then
        call mpp_clock_begin ( cosp_term_clock )
        call cosp_driver_end
        call mpp_clock_end ( cosp_term_clock )
      endif

!---------------------------------------------------------------------
!    deallocate the module variables.
!---------------------------------------------------------------------
      deallocate (diff_cu_mo, diff_t, diff_m, pbltop,  hmint, cush, cbmf, cgust, tke, pblhto,&
                  rkmo, taudpo, exist_shconv, exist_dpconv, convect, radturbten, r_convect)
      deallocate (fl_lsrain, fl_lssnow, fl_lsgrpl, fl_ccrain, fl_ccsnow, &
                  fl_donmca_snow, fl_donmca_rain, &
                  temp_last, q_last)
 
      deallocate (id_tracer_phys_vdif_dn)
      deallocate (id_tracer_phys_vdif_up)
      deallocate (id_tracer_phys_turb)
      deallocate (id_tracer_phys_moist)

      if (do_cosp .or. do_modis_yim) then
         deallocate (tsurf_save)
      endif

      call dealloc_physics_tendency_type (Physics_tendency)

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.


!-----------------------------------------------------------------------

 end subroutine physics_driver_end

!#######################################################################
! <SUBROUTINE NAME="physics_driver_restart">
!
! <DESCRIPTION>
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine physics_driver_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp


  if (mpp_pe() == mpp_root_pe() ) then
     call error_mesg('physics_driver_mod', 'Writing netCDF formatted restart file: RESTART/physics_driver.res.nc', NOTE)
  endif
  call physics_driver_netcdf(timestamp)
  call vert_turb_driver_restart(timestamp)

!    call moist_processes_restart(timestamp)
  call damping_driver_restart(timestamp)

end subroutine physics_driver_restart
! </SUBROUTINE> NAME="physics_driver_restart"

! <SUBROUTINE NAME="physics_driver_netcdf">
!
! <DESCRIPTION>
! Write out restart file for physics driver.
! This routine is needed so that physics_driver_restart and physics_driver_end
! can call a routine which will not result in multiple copies of restart files 
! being written by the destructor routines.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
! </DESCRIPTION>
!
subroutine physics_driver_netcdf(timestamp)
  character(len=*), intent(in), optional :: timestamp

    r_convect = 0.
    where(convect)
       r_convect = 1.0
    end where
    call save_restart(Phy_restart, timestamp)
    if(in_different_file) call save_restart(Til_restart, timestamp)

end subroutine physics_driver_netcdf
! </SUBROUTINE> NAME="physics_driver_netcdf"

!#######################################################################
! <FUNCTION NAME="do_moist_in_phys_up">
!  <OVERVIEW>
!    do_moist_in_phys_up returns the value of do_moist_processes
!  </OVERVIEW>
!  <DESCRIPTION>
!    do_moist_in_phys_up returns the value of do_moist_processes
!  </DESCRIPTION>
!  <TEMPLATE>
!   logical = do_moist_in_phys_up()
!  </TEMPLATE>
! </FUNCTION>
!
function do_moist_in_phys_up()

!--------------------------------------------------------------------
!    do_moist_in_phys_up returns the value of do_moist_processes
!----------------------------------------------------------------------

logical :: do_moist_in_phys_up

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('do_moist_in_phys_up',  &
              'module has not been initialized', FATAL)
      endif
 
!-------------------------------------------------------------------
!    define output variable.
!-------------------------------------------------------------------
      do_moist_in_phys_up = do_moist_processes

 
end function do_moist_in_phys_up

!#####################################################################
! <FUNCTION NAME="get_diff_t">
!  <OVERVIEW>
!    returns the values of array diff_t
!  </OVERVIEW>
!  <DESCRIPTION>
!    returns the values of array diff_t
!  </DESCRIPTION>
!  <TEMPLATE>
!   diff_t(:,:,:) = get_diff_t()
!  </TEMPLATE>
! </FUNCTION>
!
!#####################################################################
function get_diff_t() result(diff_t_out)
real, dimension(size(diff_t,1),size(diff_t,2),size(diff_t,3)) :: diff_t_out

  if ( .not. module_is_initialized) then
    call error_mesg ('get_diff_t','module has not been initialized', FATAL)
  endif

  diff_t_out = diff_t

end function get_diff_t

!#####################################################################
! <FUNCTION NAME="get_radturbten">
!  <OVERVIEW>
!    returns the values of array radturbten
!  </OVERVIEW>
!  <DESCRIPTION>
!    returns the values of array radturbten
!  </DESCRIPTION>
!  <TEMPLATE>
!   radturbten(:,:,:) = get_radturbten()
!  </TEMPLATE>
! </FUNCTION>
!
!#####################################################################
function get_radturbten() result(radturbten_out)
real, dimension(size(radturbten,1),size(radturbten,2),size(radturbten,3)) :: radturbten_out

  if ( .not. module_is_initialized) then
    call error_mesg ('get_radturbten','module has not been initialized', FATAL)
  endif

  radturbten_out = radturbten

end function get_radturbten
!#####################################################################
! <SUBROUTINE NAME="zero_radturbten">
!  <OVERVIEW>
!    sets all values of array radturbten to zero
!  </OVERVIEW>
!  <DESCRIPTION>
!    sets all values of array radturbten to zero
!  </DESCRIPTION>
!  <TEMPLATE>
!   call zero_radturbten()
!  </TEMPLATE>
! </SUBROUTINE>
!
!#####################################################################
subroutine zero_radturbten()

  if ( .not. module_is_initialized) then
    call error_mesg ('zero_radturbten','module has not been initialized', FATAL)
  endif

  radturbten = 0.0

end subroutine zero_radturbten



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
               
     
!#####################################################################
! <SUBROUTINE NAME="physics_driver_register_restart">
!  <OVERVIEW>
!    physics_driver_register_restart will register restart field when do_netcdf file 
!    is true. 
!  </OVERVIEW>
subroutine physics_driver_register_restart (Restart)
  type(clouds_from_moist_block_type), intent(inout), target :: Restart
  character(len=64) :: fname, fname2
  integer           :: id_restart
  integer           :: nc
  logical           :: reproduce_ulm_restart = .true.
  integer           :: index_strat

  if (do_moist_processes) then  
    if(doing_strat) then 
       now_doing_strat = 1
    else
       now_doing_strat = 0
    endif

    if(doing_edt) then 
       now_doing_edt = 1
    else
       now_doing_edt = 0
    endif

    if(doing_entrain) then 
       now_doing_entrain = 1
    else
       now_doing_entrain = 0
    endif
  endif

  fname = 'physics_driver.res.nc'
  call get_mosaic_tile_file(fname, fname2, .false. ) 
  allocate(Phy_restart)
  if(trim(fname2) == trim(fname)) then
     Til_restart => Phy_restart
     in_different_file = .false.
  else
     in_different_file = .true.
     allocate(Til_restart)
  endif

  id_restart = register_restart_field(Phy_restart, fname, 'vers',          vers,              no_domain=.true.)
  id_restart = register_restart_field(Phy_restart, fname, 'doing_strat',   now_doing_strat,   no_domain=.true.)
  id_restart = register_restart_field(Phy_restart, fname, 'doing_edt',     now_doing_edt,     no_domain=.true.)
  id_restart = register_restart_field(Phy_restart, fname, 'doing_entrain', now_doing_entrain, no_domain=.true.)

  id_restart = register_restart_field(Til_restart, fname, 'diff_cu_mo', diff_cu_mo)
  id_restart = register_restart_field(Til_restart, fname, 'pbltop',     pbltop)
  id_restart = register_restart_field(Til_restart, fname, 'cush',       cush, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'cbmf',       cbmf, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'hmint',      hmint, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'cgust',      cgust, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'tke',        tke, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'pblhto',     pblhto, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'rkmo',       rkmo, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'taudpo',     taudpo, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'exist_shconv', exist_shconv, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'exist_dpconv', exist_dpconv, mandatory = .false.)
  id_restart = register_restart_field(Til_restart, fname, 'diff_t',     diff_t)
  id_restart = register_restart_field(Til_restart, fname, 'diff_m',     diff_m)
  id_restart = register_restart_field(Til_restart, fname, 'convect',    r_convect) 
  if (do_clubb > 0) then
    id_restart = register_restart_field(Til_restart, fname, 'diff_t_clubb', diff_t_clubb, mandatory = .false.)
  end if
  if (doing_strat) then
    id_restart = register_restart_field(Til_restart, fname, 'radturbten',       radturbten)
  endif

  index_strat = 0
  do nc = 1, size(Restart%Cloud_data,1)
    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'strat_cloud' .and. .not. reproduce_ulm_restart) then
      id_restart = register_restart_field(Til_restart, fname, 'lsc_cloud_area',     Restart%Cloud_data(nc)%cloud_area,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_liquid',         Restart%Cloud_data(nc)%liquid_amt,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_ice',            Restart%Cloud_data(nc)%ice_amt,        mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_droplet_number', Restart%Cloud_data(nc)%droplet_number, mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_ice_number',     Restart%Cloud_data(nc)%ice_number,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_snow',           Restart%Cloud_data(nc)%snow,           mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_rain',           Restart%Cloud_data(nc)%rain,           mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_snow_size',      Restart%Cloud_data(nc)%snow_size,      mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_rain_size',      Restart%Cloud_data(nc)%rain_size,      mandatory = .false.)
    endif
    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'strat_cloud' .and. reproduce_ulm_restart) index_strat = nc

    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'donner_cell') then
      id_restart = register_restart_field(Til_restart, fname, 'cell_cloud_frac',  Restart%Cloud_data(nc)%cloud_area,  mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'cell_liquid_amt',  Restart%Cloud_data(nc)%liquid_amt,  mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'cell_liquid_size', Restart%Cloud_data(nc)%liquid_size, mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'cell_ice_amt',     Restart%Cloud_data(nc)%ice_amt,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'cell_ice_size',    Restart%Cloud_data(nc)%ice_size,    mandatory = .false.)
    endif

    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'donner_meso') then
      id_restart = register_restart_field(Til_restart, fname, 'meso_cloud_frac',  Restart%Cloud_data(nc)%cloud_area,  mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'meso_liquid_amt',  Restart%Cloud_data(nc)%liquid_amt,  mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'meso_liquid_size', Restart%Cloud_data(nc)%liquid_size, mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'meso_ice_amt',     Restart%Cloud_data(nc)%ice_amt,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'meso_ice_size',    Restart%Cloud_data(nc)%ice_size,    mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'nsum',             Restart%Cloud_data(nc)%nsum_out,    mandatory = .false.)
    endif

    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'uw_conv') then
      id_restart = register_restart_field(Til_restart, fname, 'shallow_cloud_area',     Restart%Cloud_data(nc)%cloud_area,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'shallow_liquid',         Restart%Cloud_data(nc)%liquid_amt,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'shallow_ice',            Restart%Cloud_data(nc)%ice_amt,        mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'shallow_droplet_number', Restart%Cloud_data(nc)%droplet_number, mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'shallow_ice_number',     Restart%Cloud_data(nc)%ice_number,     mandatory = .false.)
    endif
  enddo

    ! save large-scale clouds last to reproduce ulm code
    if (index_strat > 0) then
      nc = index_strat
      id_restart = register_restart_field(Til_restart, fname, 'lsc_cloud_area',     Restart%Cloud_data(nc)%cloud_area,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_liquid',         Restart%Cloud_data(nc)%liquid_amt,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_ice',            Restart%Cloud_data(nc)%ice_amt,        mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_droplet_number', Restart%Cloud_data(nc)%droplet_number, mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_ice_number',     Restart%Cloud_data(nc)%ice_number,     mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_snow',           Restart%Cloud_data(nc)%snow,           mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_rain',           Restart%Cloud_data(nc)%rain,           mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_snow_size',      Restart%Cloud_data(nc)%snow_size,      mandatory = .false.)
      id_restart = register_restart_field(Til_restart, fname, 'lsc_rain_size',      Restart%Cloud_data(nc)%rain_size,      mandatory = .false.)
    endif

end subroutine physics_driver_register_restart
! </SUBROUTINE>    
!#####################################################################
! <SUBROUTINE NAME="check_args">
!  <OVERVIEW>
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call check_args (lat, lon, area, p_half, p_full, z_half, z_full,&
!                        u, v, t, q, r, um, vm, tm, qm, rm,             &
!                        udt, vdt, tdt, qdt, rdt)
!  </TEMPLATE>
!  <IN NAME="lat" TYPE="real">
!   array of model latitudes at model points [radians]
!  </IN>
!  <IN NAME="lon" TYPE="real">
!   array of model longitudes at model points [radians]
!  </IN>
!  <IN NAME="area" TYPE="real">
!   grid box area - current not used
!  </IN>
!  <IN NAME="p_half" TYPE="real">
!   pressure at model interface levels (offset from t,q,u,v,r)
!  </IN>
!  <IN NAME="p_full" TPYE="real">
!   pressure at full levels
!  </IN>
!  <IN NAME="z_half" TYPE="real">
!   height at model interface levels
!  </IN>
!  <IN NAME="z_full" TPYE="real">
!   height at full levels
!  </IN>
!  <IN NAME="u" TYPE="real">
!   zonal wind at current time step
!  </IN>
!  <IN NAME="v" TYPE="real">
!   meridional wind at current time step
!  </IN>
!  <IN NAME="t" TYPE="real">
!   temperature at current time step
!  </IN>
!  <IN NAME="q" TYPE="real">
!   specific humidity at current time step
!  </IN>
!  <IN NAME="r" TPYE="real">
!   multiple 3d tracer fields at current time step
!  </IN>
!  <IN NAME="um" TYPE="real">
!   zonal wind at previous time step
!  </IN>
!  <IN NAME="vm" TYPE="real">
!   meridional wind at previous time step
!  </IN>
!  <IN NAME="tm" TYPE="real">
!   temperature at previous time step
!  </IN>
!  <IN NAME="qm" TYPE="real">
!   specific humidity at previous time step
!  </IN>
!  <IN NAME="rm" TPYE="real">
!   multiple 3d tracer fields at previous time step
!  </IN>
!  <IN NAME="udt" TYPE="real">
!   zonal wind tendency
!  </IN>
!  <IN NAME="vdt" TYPE="real">
!   meridional wind tendency
!  </IN>
!  <IN NAME="tdt" TYPE="real">
!   temperature tendency
!  </IN>
!  <IN NAME="qdt" TYPE="real">
!   moisture tracer tendencies
!  </IN>
!  <IN NAME="rdt" TYPE="real">
!   multiple tracer tendencies
!  </IN>
! </SUBROUTINE>
!
subroutine check_args (lat, lon, area, p_half, p_full, z_half, z_full,&
                        u, v, t, q, r, um, vm, tm, qm, rm,             &
                        udt, vdt, tdt, qdt, rdt, rdiag)

!----------------------------------------------------------------------
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!-----------------------------------------------------------------------

real,    dimension(:,:),    intent(in)           :: lat, lon, area
real,    dimension(:,:,:),  intent(in)           :: p_half, p_full,   &
                                                    z_half, z_full,   &
                                                    u, v, t, q, um, vm, &
                                                    tm, qm
real,    dimension(:,:,:,:),intent(in)           :: r, rm
real,    dimension(:,:,:),  intent(in)           :: udt, vdt, tdt, qdt
real,    dimension(:,:,:,:),intent(in)           :: rdt
real,    dimension(:,:,:,ntp+1:),intent(in)      :: rdiag

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous 
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency 
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!
!   intent(in), optional:
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      integer ::  id, jd, kd  ! model dimensions on the processor  
      integer ::  ierr        ! error flag

!--------------------------------------------------------------------
!    define the sizes that the arrays should be.
!--------------------------------------------------------------------
      id = size(u,1) 
      jd = size(u,2) 
      kd = size(u,3) 

!--------------------------------------------------------------------
!    check the dimensions of each input array. if they are incompat-
!    ible in size with the standard, the error flag is set to so
!    indicate.
!--------------------------------------------------------------------
      ierr = 0
      ierr = ierr + check_dim (lat, 'lat',  id,jd)
      ierr = ierr + check_dim (lon, 'lon',  id,jd)
      ierr = ierr + check_dim (area,'area', id,jd)

      ierr = ierr + check_dim (p_half,'p_half', id,jd,kd+1)
      ierr = ierr + check_dim (p_full,'p_full', id,jd,kd)
      ierr = ierr + check_dim (z_half,'z_half', id,jd,kd+1)
      ierr = ierr + check_dim (z_full,'z_full', id,jd,kd)

      ierr = ierr + check_dim (u, 'u',  id,jd,kd)
      ierr = ierr + check_dim (v, 'v',  id,jd,kd)
      ierr = ierr + check_dim (t, 't',  id,jd,kd)
      ierr = ierr + check_dim (q, 'q',  id,jd,kd)
      ierr = ierr + check_dim (um,'um', id,jd,kd)
      ierr = ierr + check_dim (vm,'vm', id,jd,kd)
      ierr = ierr + check_dim (tm,'tm', id,jd,kd)
      ierr = ierr + check_dim (qm,'qm', id,jd,kd)

      ierr = ierr + check_dim (udt,'udt', id,jd,kd)
      ierr = ierr + check_dim (vdt,'vdt', id,jd,kd)
      ierr = ierr + check_dim (tdt,'tdt', id,jd,kd)
      ierr = ierr + check_dim (qdt,'qdt', id,jd,kd)

      if (ntp > 0) then
        ierr = ierr + check_dim (r,  'r',   id,jd,kd,ntp)
        ierr = ierr + check_dim (rm, 'rm',  id,jd,kd,ntp)
      endif
      if (ntp > 0) then
        ierr = ierr + check_dim (rdt,'rdt', id,jd,kd,ntp)
      endif
      if (nt > ntp) then
        ierr = ierr + check_dim (rdiag,'rdiag', id,jd,kd,nt-ntp)
      endif

!--------------------------------------------------------------------
!    if any problems were detected, exit with an error message.
!--------------------------------------------------------------------
      if (ierr > 0) then
        call error_mesg ('physics_driver_mod', 'bad dimensions', FATAL)
      endif

!-----------------------------------------------------------------------


      end subroutine check_args


!#######################################################################
! <FUNCTION NAME="check_dim_2d">
!  <OVERVIEW>
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_2d (data,name,id,jd) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd" TYPE="integer">
!   expected i and j dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_2d (data,name,id,jd) result (ierr)

!--------------------------------------------------------------------
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

real,    intent(in), dimension(:,:) :: data
character(len=*), intent(in)        :: name
integer, intent(in)                 :: id, jd
integer                             :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd      expected i and j dimensions
!     
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
             'dimension 1 of argument ' //  &
              name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
           call error_mesg ('physics_driver_mod',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif

!----------------------------------------------------------------------

      end function check_dim_2d

!#######################################################################
! <FUNCTION NAME="check_dim_3d">
!  <OVERVIEW>
!    check_dim_3d compares the size of three-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_3d compares the size of three-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_3d (data,name,id,jd, kd) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd, kd" TYPE="integer">
!   expected i, j and k dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_3d (data,name,id,jd,kd) result (ierr)

!--------------------------------------------------------------------
!    check_dim_3d compares the size of thr1eedimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

real,    intent(in), dimension(:,:,:) :: data
character(len=*), intent(in)          :: name
integer, intent(in)                   :: id, jd, kd
integer  ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd,kd   expected i, j and k dimensions
!     
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
        call error_mesg ('physics_driver_mod',  &
              'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif

!---------------------------------------------------------------------


      end function check_dim_3d


!#######################################################################
! <FUNCTION NAME="check_dim_4d">
!  <OVERVIEW>
!    check_dim_4d compares the size of four-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </OVERVIEW>
!  <DESCRIPTION>
!    check_dim_4d compares the size of four-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!  </DESCRIPTION>
!  <TEMPLATE>
!    check_dim_4d (data,name,id,jd, kd, nt) result (ierr)
!  </TEMPLATE>
!  <IN NAME="data" TYPE="real">
!   array of data to be checked
!  </IN>
!  <IN NAME="name" TYPE="character">
!   name associated with array to be checked
!  </IN>
!  <IN NAME="id, jd, kd, nt" TYPE="integer">
!   expected i, j, k and 4th dimensions
!  </IN>
! </FUNCTION>
!
function check_dim_4d (data,name,id,jd,kd,nt) result (ierr)

!--------------------------------------------------------------------
!    check_dim_4d compares the size of four dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------
real,    intent(in), dimension(:,:,:,:) :: data
character(len=*), intent(in)            :: name
integer, intent(in)                     :: id, jd, kd, nt
integer                                 :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data          array to be checked
!     name          name associated with array to be checked
!     id,jd,kd,nt   expected i, j and k dimensions
!     
!  result variable:
!
!     ierr          set to 0 if ok, otherwise is a count of the number
!                   of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,4) /= nt) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 4 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif

!---------------------------------------------------------------------


      end function check_dim_4d


 
                end module physics_driver_mod
