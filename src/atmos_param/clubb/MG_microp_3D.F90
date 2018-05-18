module MG_microp_3D_mod

! Driver for the Morrison-Gettelman microphysics designed to be used in conjunction
! with the CLUBB cloud boundary layer turbulence.
!
! Authors:
!
! Huan Guo (Huan.Guo@noaa.gov)
! Chris Golaz (Chris.Golaz@noaa.gov)

  ! --- external modules ---

  use              fms_mod, only :  file_exist, open_namelist_file, close_file,    &
                                    error_mesg, FATAL, check_nml_error
  use        constants_mod, only :  cp_air, hlv, hlf, hls, tfreeze,  &
                                          rdgas, grav, rvgas
  use        cloud_rad_mod, only :  cloud_rad_init
  use     time_manager_mod, only :  time_type, get_time, set_date
  use   cldwat2m_micro_mod, only :  mmicro_pcond, ini_micro
  use micro_mg_mod,         only:   micro_mg_init, micro_mg_get_cols,&
                                          micro_mg_tend


  use strat_cloud_utilities_mod, only: strat_cloud_utilities_init, &
                                       diag_id_type, diag_pt_type, &
                                       strat_nml_type

  use strat_netcdf_mod,          only: strat_netcdf, strat_netcdf_init, &
                                       strat_netcdf_end
  implicit none


  ! --- available public interfaces ---
  public  MG_microp_3D_init, MG_microp_3D, MG_microp_3D_end
        
  ! --- namelist parameters ---

  !   parameter              definition                  unit
  !  ------------     -----------------------------   ---------------
  !
  !   N_land           fixed number of cloud drops     1/(m*m*m)
  !                    per unit volume in liquid 
  !                    clouds on land
  !
  !   N_ocean          fixed number of cloud drops     1/(m*m*m)
  !                    per unit volume in liquid 
  !                    clouds over ocean
  !
  !   qmin             minimum permissible value of    kg condensate/
  !                    cloud liquid, cloud ice,        kg air
  !                    saturated volume fraction,
  !                    or rain and snow areas
  !
  !                    NOTE: qmin should be chosen
  !                    such that the range of
  !                    {qmin, max(qa,ql,qi)} is
  !                    resolved by the precision 
  !                    of the numbers used.
  !
  !   do_liq_num       Should the prognostic droplet number
  !                    concentration be used?
  !
  !   overlap          value of the overlap parameter
  !                    from cloud rad
  !                    overlap = 1 is maximum-random
  !                    overlap = 2 is random
  !
  !   init_date        Option to delay start of microphysics by micro_begin_sec
  !   micro_begin_sec  seconds after init_data


  INTEGER               :: n_diag_4d, n_diag_4d_kp1
  TYPE(diag_id_type)    :: diag_id
  TYPE(diag_pt_type)    :: diag_pt
  

  real                  :: N_land           = 250.e+06
  real                  :: N_ocean          = 100.e+06
  real                  :: qmin             = 1.e-14
  integer               :: overlap          = 2
  logical               :: do_liq_num       = .true.
  logical               :: do_ice_num       = .false.
  logical               :: use_Meyers       = .false.
  logical               :: use_Cooper       = .false.
  integer, dimension(6) :: init_date        = (/ 1, 1, 1, 0, 0, 0 /)
  real                  :: micro_begin_sec  = 0.0
  integer               :: override_liq_num = 0
  integer               :: override_ice_num = 0

  logical               :: use_ndust = .true.
  real                  :: qcvar = 2.0
  real                  :: rhmini=0.80  ! minimum rh for ice cld fraction > 0
  logical               :: do_mg_ncar_microphys = .true.
  logical               :: do_ncar_microphys = .false.
  logical               :: microp_uniform = .false. 
                               ! .true. = configure uniform for sub-columns 
                               ! .false. = use w/o sub-columns (standard)
  logical               :: do_cldice = .true.
                               ! .true. = do all processes (standard)
                               ! .false. = skip all processes affecting
                               !           cloud ice
   
  namelist /MG_microp_3D_nml/ &
       N_land, N_ocean, qmin, overlap, do_liq_num, init_date, micro_begin_sec, &
       override_liq_num, override_ice_num, use_Meyers, use_Cooper, use_ndust, qcvar, do_ice_num, rhmini, &
       do_mg_ncar_microphys,  do_ncar_microphys,  microp_uniform, do_cldice

  integer, parameter    :: r8 = selected_real_kind(12)   ! 8 byte real
  logical               :: module_is_initialized = .false.
  integer, save         :: current_days0, current_sec0

contains

#ifdef CLUBB
!##############################################################################
subroutine MG_microp_3D_init(axes,Time,idim,jdim,kdim)

  ! --- calling arguments ---
  integer, intent (in)        :: axes(4)          ! x,y,z,z_half axes types
  integer, intent (in)        :: idim,jdim,kdim   ! dimensions
  type(time_type), intent(in) :: Time             ! time

  ! --- internal variables ---
  integer                     :: unit, io, ierr
  type(time_type)             :: Time_init
  character(len=128)          :: errstring ! Output status: non-blank for 
                                           ! error return


  ! --- is module already initialized ? ---
  if (module_is_initialized) then
    return
  else
    module_is_initialized = .true.
  endif

  ! --- read namelist ---
  if ( file_exist('input.nml')) then
    unit = open_namelist_file()
    ierr = 1
    do while (ierr /= 0)
      read(unit, nml=MG_microp_3D_nml, iostat=io, end=10)
      ierr = check_nml_error(io,'MG_microp_3D_nml')
    enddo
10  call close_file(unit)
  endif

  ! --- initialize cloud_rad_mod ---
  if (do_liq_num) then
    if( do_ice_num) then
      call cloud_rad_init( axes, Time, qmin_in=qmin,                         &
                           N_land_in=N_land, N_ocean_in=N_ocean,             &
                           prog_droplet_in=do_liq_num, overlap_out=overlap,  &
                           qcvar_in = qcvar,   prog_ice_num_in=do_ice_num )
    else
      call cloud_rad_init( axes, Time, qmin_in=qmin,                         &
                           N_land_in=N_land, N_ocean_in=N_ocean,             &
                           prog_droplet_in=do_liq_num, overlap_out=overlap )
    endif
  else
    call cloud_rad_init( axes, Time, qmin_in=qmin,                         &
                         N_land_in=N_land, N_ocean_in=N_ocean,             &
                         overlap_out=overlap )
  endif

  ! --- get namelist initial time to determine whether microphysics is active ---
  Time_init = set_date( init_date(1), init_date(2), init_date(3),  &
                        init_date(4), init_date(5), init_date(6) )
  call get_time( Time_init, current_sec0, current_days0)

  ! --- initialize diagnostics ---
  call strat_netcdf_init (axes, Time, diag_id, diag_pt, n_diag_4d, &
                          n_diag_4d_kp1)

  if( do_mg_ncar_microphys) then
    call ini_micro (grav, rdgas, rvgas, cp_air,   &
                                  tfreeze, hlv, hlf, rhmini)
  elseif ( do_ncar_microphys) then
    call micro_mg_init (r8, grav, rdgas, rvgas, cp_air, tfreeze, &
                            hlv, hlf, rhmini, microp_uniform, do_cldice, &
                            errstring )
    if (trim(errstring) /= '') then
       call error_mesg ('strat_cloud/microphysics/micro_mg_init', &
                                                errstring, FATAL)
    endif
  endif

end subroutine MG_microp_3D_init
!##############################################################################


!##############################################################################
subroutine MG_microp_3D( Time, is, ie, js, je, lon, lat, dtcloud,              &
                         pfull3d, phalf3d, zhalf3d, LAND,                      &
                         T3d, qv3d, ql3d, qi3d, qa3d, qn3d, qni3d, ahuco3d,    &
                         dcond_ls_liquid, dcond_ls_ice,                        &
                         Ndrop_act_CLUBB, Icedrop_act_CLUBB,                   &
                         ndust, rbar_dust,                                     &
                         ST3d, SQ3d, SL3d, SI3d, SA3d, SN3d, SNi3d,            &
                         f_snow_berg3d,                                        &
                         rain3d, snow3d, surfrain, surfsnow,                   &
                         do_clubb,  qcvar_clubb, MASK3d,                       &
                         lsc_snow, lsc_rain, lsc_snow_size, lsc_rain_size )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !       ------
  !       INPUT:
  !       ------
  !         variable              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !       Time           time type variable 
  !
  !       is,ie          starting and ending i indices 
  !                      for data window
  !
  !       js,je          starting and ending j indices 
  !                      for data window
  !
  !       lon,lat        cell longitude and latitude     
  !
  !       dtcloud        time between this call and      s
  !                      the next call to MG_microp_3D
  !
  !       pfull3d        pressure at full model levels   Pa
  !                      IMPORTANT NOTE: p(j)<p(j+1)
  !
  !       phalf3d        pressure at half model levels   Pa
  !                      phalf(j)<pfull(j)<phalf(j+1)
  !
  !       zhalf3d        height at half model levels     m
  !
  !       LAND           the fraction of the grid box    fraction
  !                      covered by land
  !                               
  !       T3d            temperature                     K
  !
  !       qv3d           specific humidity of water      kg vapor/kg air
  !                      vapor
  !
  !       ql3d           specific humidity of cloud      kg condensate/
  !                      liquid                          kg air
  !
  !       qi3d           specific humidity of cloud      kg condensate/
  !                      ice                             kg air
  !
  !       qa3d           saturated volume fraction       fraction
  !
  !       qn3d           cloud droplet number            #/kg air
  !
  !       qni3d          cloud droplet number            #/kg air
  !
  !       ahuco3d        fraction, cell+meso, from       fraction
  !                      donner_deep
  !                      index 1 nearest ground
  !
  !       dcond_ls_liquid   change in liquid due to      kg condensate/
  !                         non-convective condensation.    kg air
  !
  !       dcond_ls_ice   change in ice due to            kg condensate/
  !                      non-convective condensation.    kg air
  !
  !       Ndrop_act_CLUBB  activated droplet number      #/
  !                        concentration.                kg air
  !
  !       Icedrop_act_CLUBB   activated ice crystal      #/
  !                           number concentration.      kg air  
  !
  !       -------
  !       OUTPUT:
  !       -------
  !
  !       ST3d           temperature change due to       K
  !                      all stratiform processes
  !
  !       SQ3d           water vapor change due to       kg vapor/kg air
  !                      all stratiform processes
  !
  !       SL3d           cloud liquid change due to      kg condensate/
  !                      all stratiform processes        kg air
  !
  !       SI3d           cloud ice change due to         kg condensate/
  !                      all stratiform processes        kg air
  !
  !       SA3d           saturated volume fraction       fraction
  !                      change due to all stratiform 
  !                      processes
  !
  !       SN3d           cloud droplet number            #/kg air
  !                      change due to all stratiform 
  !                      processes
  !
  !       SNi3d          cloud droplet number            #/kg air
  !                      change due to all stratiform 
  !                      processes
  !
  !       f_snow_berg3d  fraction of Bergeron process in total snow production
  !                    
  !       rain3d         rain that falls through the     kg condensate/
  !                      edge of the model layer         (kg liquid/m2/sec)
  !
  !       snow3d         snow that falls through the     kg condensate/
  !                      edge of the model layer         (kg snow/m2/sec)
  !
  !       surfrain       rain that falls through the     kg condensate/
  !                      bottom of the column over       (m/sec)
  !                      the time dtcloud
  !
  !       surfsnow       snow that falls through the     kg condensate/
  !                      bottom of the column over       (m/sec)
  !                      the time dtcloud
  !
  !       ---------------
  !       optional INPUT:
  !       ---------------
  !
  !       MASK3d         real array indicating the 
  !                      point is above the surface
  !                      if equal to 1.0 and 
  !                      indicating the point is below
  !                      the surface if equal to 0.
  !
  !       ----------------
  !       optional OUTPUT:
  !       ----------------
  !
  !       lsc_snow       precipitation mass and size
  !       lsc_rain
  !       lsc_snow_size
  !       lsc_rain_size

  ! --- calling arguments ---
  type(time_type), intent (in)                         :: Time
  integer, intent (in)                                 :: is,ie,js,je
  real, intent (in),    dimension(:,:)                 :: lon,lat
  real, intent (in)                                    :: dtcloud
  real, intent (in),    dimension(:,:,:)               :: pfull3d,phalf3d
  real, intent (in),    dimension(:,:,:)               :: zhalf3d
  real, intent (in),    dimension(:,:)                 :: LAND
  real, intent (in),    dimension(:,:,:)               :: T3d,qv3d,ql3d,qi3d,qa3d,qn3d,qni3d
  real, intent (in),    dimension(:,:,:)               :: ahuco3d
  real, intent (in),    dimension(:,:,:)               :: dcond_ls_liquid,dcond_ls_ice
  real, intent (in),    dimension(:,:,:)               :: Ndrop_act_CLUBB,Icedrop_act_CLUBB
  real, intent (in),    dimension(:,:,:)               :: ndust, rbar_dust
  real, intent (out),   dimension(:,:,:)               :: ST3d,SQ3d,SL3d,SI3d,SA3d,SN3d,SNi3d
  real, intent (out),   dimension(:,:,:)               :: f_snow_berg3d  
  real, intent (out),   dimension(:,:,:)               :: rain3d,snow3d
  real, intent (out),   dimension(:,:)                 :: surfrain,surfsnow

  integer, intent (in),  optional                      :: do_clubb
  real, intent (in),  optional, dimension(:,:,:)       :: qcvar_clubb
  real, intent (in),  optional, dimension(:,:,:)       :: MASK3d
  real, intent (out), optional, dimension(:,:,:)       :: lsc_snow,      &
                                                          lsc_rain,      &
                                                          lsc_snow_size, &
                                                          lsc_rain_size

  ! --- internal variables ---
  integer                                              :: idim,jdim,kdim
  integer                                              :: i, j, k, nn, unit

  real, dimension(size(pfull3d,1),size(pfull3d,3))     :: pfull
  real, dimension(size(phalf3d,1),size(phalf3d,3))     :: phalf

  real, dimension(size(T3d,1),size(T3d,3))             :: T,qv,ql,qi,qa,qn,qni,tn,qvn
  real, dimension(size(T3d,1),size(T3d,3))             :: ahuco
  real, dimension(size(T3d,1),size(T3d,3))             :: ST,SQ,SL,SI,SA,SN,SNi

  real, dimension(size(T3d,1),size(T3d,3))             :: ST_micro,SQ_micro,SL_micro, &
                                                          SI_micro,SN_micro,SNi_micro

  logical                                              :: mass_cons = .true.
  real                                                 :: m1, m2, scalef

  real, dimension(size(T3d,1),size(T3d,3))             :: ql_upd,qi_upd,qa_upd,qn_upd,qni_upd
  real, dimension(size(T3d,1),size(T3d,3))             :: liqcldf, icecldf

  real, dimension(size(T3d,1),size(T3d,3))             :: drop2, crystal1
  real, dimension(size(T3d,1),size(T3d,3),4)           :: rbar_dust_4bin, ndust_4bin

  ! rain and snow mixing ratios from the Morrison scheme (kg/kg)
  !  2D  (i,k)
  real, dimension(size(T3d,1),size(T3d,3))             :: qsout2d_mg,     qrout2d_mg
  real, dimension(size(T3d,1),size(T3d,3))             :: reff_rain2d_mg, reff_snow2d_mg
  !  3D  (i,j,k)
  real, dimension(size(T3d,1),size(T3d,2),size(T3d,3)) :: reff_rain3d_mg, reff_snow3d_mg

  ! rain and snow fluxes from the Morrison scheme ( kg/m2/sec )
  !  2D  (i,k)
  real, dimension(size(T3d,1),size(T3d,3)+1)           :: rain2d_mg, snow2d_mg

  ! pressure difference across level (pa) as inputs for MG microphysics
  real, dimension(size(T3d,1),size(T3d,3))             :: pdel

  ! diagnostic variables
  real, allocatable, dimension(:,:,:,:)                :: diag_4d, diag_4d_kp1 
  real, allocatable, dimension(:,:,:)                  :: diag_3d

  ! set up time to determine whether microphysics is active or not
  real                                                 :: current_total_sec
  integer                                              :: current_sec, current_days
  real                                                 :: rho_air

  real, dimension(size(T3d,1),size(T3d,2))             :: enth_micro_col,  wat_micro_col
!--->h1g, 06-12-2013
  real, dimension(size(T3d,1),size(T3d,3))             :: delta_cf
  real, dimension(size(T3d,1),size(T3d,3))             :: D_eros_l4, nerosc4, D_eros_i4, nerosi4, dqcdt, dqidt
  real, dimension(size(T3d,1),size(T3d,3))             :: f_snow_berg, qa0, gamma_mg, SA_0, ssat_disposal
  type(strat_nml_type)                                 :: Nml
!<---h1g, 06-12-2013

!--->h1g, 04-12-2014
  real, dimension(size(T3d,1),size(T3d,3))             :: relvar
  integer                                              :: top_lev, nlev
  character(len=128)                                   :: errstring
  real, dimension(size(T3d,1),size(T3d,3))             :: accre_enhann, tnd_qsnown, &
                                                          tnd_nsnown, re_icen 
  integer,dimension(:),allocatable                     :: mgcols         
  integer                                              :: mgncol
!<---h1g, 04-12-2014


  ! --- is module initialized ? ---
  if (.not.module_is_initialized) then
    call error_mesg('MG_microp_3D','MG_microp_3D is not initialized',FATAL)
  endif

  ! --- initialization ---
  ! 3D rain/snow fluxes, (kg/m2/s).
  rain3d(:,:,:)    = 0.
  snow3d(:,:,:)    = 0.

  ! 2D rain/snow fluxes, (kg/m2/s).
  ! (i,k) from MG microphysics
  rain2d_mg(:,:)   = 0.
  snow2d_mg(:,:)   = 0.
  ! (i,j) surface rain/snow
  surfrain(:,:)    = 0.
  surfsnow(:,:)    = 0.

  ! rain/snow mixing ratio  (kg/kg)
  !  2D
  qrout2d_mg(:,:)  = 0.
  qsout2d_mg(:,:)  = 0.

  ! rain/snow effective radius (micron)
  !  3D
  reff_rain3d_mg(:,:,:)= 0.
  reff_snow3d_mg(:,:,:)= 0.
  !  2D
  reff_rain2d_mg(:,:)= 0.
  reff_snow2d_mg(:,:)= 0.

  ! 3D tendencies
  !  (T/sec)
  ST3d(:,:,:)      = 0.
  !  (kg/kg/sec)
  SQ3d(:,:,:)      = 0.
  SL3d(:,:,:)      = 0.
  SI3d(:,:,:)      = 0.
  SA3d(:,:,:)      = 0.
  SN3d(:,:,:)      = 0.
  SNi3d(:,:,:)     = 0.

  ! 2D (i,k) activated droplet and crystal number concentrations
  drop2(:,:)       = 0.
  crystal1(:,:)    = 0.

  ! --- determine dimensions of slab ---
  idim             = SIZE(T3d,1)
  jdim             = SIZE(T3d,2)
  kdim             = SIZE(T3d,3)

!------------------------------------------------------------------------
!    allocate and initialize the diagnostic variables.
!------------------------------------------------------------------------
  if (allocated(diag_3d)) deallocate (diag_3d)
      allocate (diag_3d(idim,jdim,0:n_diag_4d))
  if (allocated(diag_4d)) deallocate (diag_4d)
      allocate (diag_4d(idim,jdim,kdim,0:n_diag_4d))
  if (allocated(diag_4d_kp1)) deallocate (diag_4d_kp1)
      allocate (diag_4d_kp1(idim,jdim,kdim+1,0:n_diag_4d_kp1))
      diag_3d(:,:,0:) = 0.
      diag_4d(:,:,:,0:) = 0.
      diag_4d_kp1(:,:,:,0:) = 0.

  !-------------------
  !--- Main j loop ---
  !-------------------
  enth_micro_col(:,:) = 0.0
  wat_micro_col(:,:)  = 0.0
  relvar(:,:) = qcvar

  j_loop: do j = 1,jdim

    ! initialization of 2D (i,k) arrays for MG microphysics
    ST(:,:)  = 0.
    SQ(:,:)  = 0.
    SL(:,:)  = 0.
    SI(:,:)  = 0.
    SA(:,:)  = 0.
    SN(:,:)  = 0.
    SNi(:,:) = 0.

    if (present(qcvar_clubb)) &
    relvar(:,:) = qcvar_clubb(:,j,:)
    T(:,:)   = T3d(:,j,:)
    qv(:,:)  = qv3d(:,j,:)
    ql(:,:)  = ql3d(:,j,:)
    qi(:,:)  = qi3d(:,j,:)
    qa(:,:)  = qa3d(:,j,:)
    qn(:,:)  = qn3d(:,j,:)
    qni(:,:) = qni3d(:,j,:)

    ahuco(:,:) = ahuco3d(:,j,:)

    phalf(:,:) = phalf3d(:,j,:)
    pfull(:,:) = pfull3d(:,j,:)
    do k = 1, kdim
      pdel(:,k) = phalf3d(:,j,k+1)- phalf3d(:,j,k)
    enddo

    !-----------------------------------------------------------------------
    !       Account for the fact that other processes may have created
    !       negative tracer or extremely small values of tracer fields.
    !
    !       In this step any values of the prognostic variables which are 
    !       less than qmin are reset to zero, while conserving total 
    !       moisture.
    !----------------------------------------------------------------------

    where (qa.le. qmin)
      SA = SA - qa
      qa_upd = 0.
    elsewhere
      qa_upd = qa
    end where



    ! Total cloud fraction should be no greater than 1.0, i.e, 
    ! sum of large-scale cloud fraction (qa_upd) and convective cloud fraction (ahuco)
    ! qa_upd+ahuco <= 1.0
    where (qa_upd .gt. (1.0 - ahuco) )
      SA = SA + (1.0 - ahuco) - qa_upd
      qa_upd = (1.0 - ahuco)
    end where

    do k = 1,kdim
      do i = 1,idim
        if ( ql(i,k).le.qmin .or. qa(i,k).le.qmin .or. qn(i,k)*1.e-6.le.qmin ) then
          SL(i,k)     = SL(i,k) - ql(i,k)
          SQ(i,k)     = SQ(i,k) + ql(i,k)
          ST(i,k)     = ST(i,k) - hlv*ql(i,k)/cp_air
          SN(i,k)     = SN(i,k) - qn(i,k)
          ql_upd(i,k) = 0.
          qn_upd(i,k) = 0.

          if (diag_id%qdt_liquid_init > 0)  diag_4d(i,j,k,diag_pt%qdt_liquid_init)  =  ql(i,k)/dtcloud
          if (diag_id%qndt_fill  + diag_id%qn_fill_col + &
                    diag_id%qldt_fill + diag_id%ql_fill_col > 0 ) diag_4d(i,j,k,diag_pt%qndt_fill) = -qn(i,k)/dtcloud 

        else
          ql_upd(i,k) = ql(i,k)
          qn_upd(i,k) = qn(i,k)
        endif
      end do
    end do

    do k = 1,kdim
      do i = 1,idim
        if ( qi(i,k).le.qmin .or. qa(i,k).le.qmin .or. qni(i,k)*1.e-3.le.qmin ) then
           SI(i,k)     = SI(i,k)  - qi(i,k)
           SQ(i,k)     = SQ(i,k)  + qi(i,k)
           ST(i,k)     = ST(i,k)  - hls*qi(i,k)/cp_air
           SNi(i,k)    = SNi(i,k) - qni(i,k)
           qi_upd(i,k) = 0.
           qni_upd(i,k)= 0.

           if (diag_id%qdt_ice_init > 0)    diag_4d(i,j,k,diag_pt%qdt_ice_init )   =  qi(i,k)/dtcloud
           if (diag_id%qnidt_fill  + diag_id%qni_fill_col > 0 ) diag_4d(i,j,k,diag_pt%qnidt_fill) = -qni(i,k)/dtcloud
        else
           qi_upd(i,k) = qi (i,k)
           qni_upd(i,k)= qni (i,k)
        endif
      end do
    end do

    ! --- updated values of T, qv, ice crystal number for MG ---
    do k = 1,kdim
      do i = 1,idim
        crystal1(i,k) = Icedrop_act_CLUBB(i,j,k)   ! crystal1  in unit of (#/kg)
        tn(i,k)       = T(i,k)  + ST(i,k)
        qvn(i,k)      = qv(i,k) + SQ(i,k)
      end do
    end do

    ! Since droplet activation has been considered and droplet number has been 
    ! updated due to activation in CLUBB, we do not include activation in 
    ! MG microphysics, i.e drop2 = 0.0
    drop2(:,:) = 0.0

    ! Option to override ice crystal number concentrations
    ! override_ice_num == 1 : Meyers formula as in original strat_cloud
    if (override_ice_num == 1) then
      if( use_Meyers .and. use_Cooper ) &
        call error_mesg ('MG_microp_3D','use_Meyers and use_Cooper cannot both be true',FATAL)
      do k = 1,kdim
        do i = 1,idim
          if( use_Meyers ) then
             crystal1(i,k) = 1000.0 * exp((12.96*0.0125*(273.15-tn(i,k)))-0.639) 
             rho_air = pfull(i,k)/tn(i,k)/rdgas
             crystal1(i,k) = crystal1(i,k)/rho_air
          elseif( use_Cooper ) then
! cooper curve (factor of 1000 is to convert from L-1 to m-3)
             crystal1(i,k) = 0.005*exp(0.304*(273.15-tn(i,k)))*1000.0
! put limit on number of nucleated crystals, set to number at T=-35 C
! then convert from m-3 to kg-1.
             rho_air = pfull(i,k)/tn(i,k)/rdgas
             crystal1(i,k) = min( crystal1(i,k), 208.9e3)/rho_air
          else ! use a constant 0.5 /kg
             crystal1(i,k) = 1.0e6 * 0.5
          end if
        end do
      end do
    end if

    ! Determine whether MG is active
    call get_time( time, current_sec, current_days)
    current_total_sec = real(current_sec - current_sec0) + 86400.0*(current_days - current_days0)

    if (current_total_sec >= micro_begin_sec) then

      ! set liquid and ice cloud fraction to be the same as total large-scale cloud fraction
      liqcldf(:,:)  = qa_upd(:,:)
      icecldf(:,:)  = qa_upd(:,:)

      ! radius of 4 dust bins for contact freezing (in the unit of m)
    do k = 1,kdim
      do i = 1,idim
        rbar_dust_4bin(i,k,1) = 5.e-6
        rbar_dust_4bin(i,k,2) = 10.e-6
        rbar_dust_4bin(i,k,3) = 15.e-6
        rbar_dust_4bin(i,k,4) = 20.e-6
      ! number in 4 dust bins for contact freezing
        if ( use_ndust ) then
          ndust_4bin(i,k,1)     = ndust(i, j, k)
          ndust_4bin(i,k,2)     = ndust(i, j, k)
          ndust_4bin(i,k,3)     = ndust(i, j, k)
          ndust_4bin(i,k,4)     = ndust(i, j, k)
        else
          ndust_4bin(i,k,1)     = 0.0
          ndust_4bin(i,k,2)     = 0.0
          ndust_4bin(i,k,3)     = 0.0
          ndust_4bin(i,k,4)     = 0.0
        endif
      enddo
    enddo

      ! initialize 2D tendency due to MG microphysics
      ST_micro(:,:)         = 0.0
      SQ_micro(:,:)         = 0.0 
      SL_micro(:,:)         = 0.0
      SI_micro(:,:)         = 0.0
      SN_micro(:,:)         = 0.0
      SNI_micro(:,:)        = 0.0


! --->h1g, 06-12-2013
!     converge to Rick Hemler's NCAR microphysics
      delta_cf(:,:)        = 0.0
      D_eros_l4(:,:)       = 0.0
      nerosc4(:,:)         = 0.0
      D_eros_i4(:,:)       = 0.0
      nerosi4(:,:)         = 0.0
      dqcdt(:,:)           = 0.0
      dqidt(:,:)           = 0.0

      qa0(:,:)             = 0.0 
      gamma_mg(:,:)        = 0.0 
      SA_0 (:,:)           = 0.0 
      ssat_disposal(:,:)   = 0.0 
                                 ! disposition of supersaturation at end 
                                 ! of step; 0.= no ssat, 1.= liq, 2.=ice)
! <---h1g, 06-12-2013
  if( do_mg_ncar_microphys) then
     call mmicro_pcond (  .false.,   .false.,   & 
                          .false.,   .false.,   &
                          j, jdim, kdim,  idim,   idim,  dtcloud,  &
                          relvar,  tn,   &
                          qvn,      ql_upd,   qi_upd,  qn_upd,   qni_upd,  pfull,  pdel, phalf,  qa_upd, &
                          liqcldf,  icecldf,  delta_cf,                                     &
                          D_eros_l4, nerosc4, D_eros_i4, nerosi4, dqcdt,                    &
                          dqidt,  crystal1, drop2,    rbar_dust_4bin,  ndust_4bin,          &
                          ST_micro,      SQ_micro,  SL_micro,       SI_micro,       &
                          SN_micro,      SNI_micro,                                 &
                          surfrain(:,j), surfsnow(:,j),rain2d_mg,   snow2d_mg,      &
                          qrout2d_mg,    qsout2d_mg,   reff_rain2d_mg, reff_snow2d_mg,  &
                          f_snow_berg, Nml, qa0, gamma_mg, SA_0, SA, &
                          ssat_disposal,  n_diag_4d, diag_4d, diag_id, &
                          diag_pt, do_clubb=do_clubb, qcvar_clubb =qcvar_clubb(:,j,:) )
  elseif ( do_ncar_microphys) then
             nlev = kdim
             top_lev = 1
             accre_enhann(:,:) = 1.0  ! accretion enhancement factor
             tnd_qsnown(:,:) = 0.     
             tnd_nsnown(:,:) = 0.
             re_icen(:,:) = 0.
 
            call micro_mg_get_cols (idim, nlev, top_lev, &
               qvn(:,:), &
               ql_upd(:,:) + dqcdt(:,:)*dtcloud, &
               qi_upd(:,:) + dqidt(:,:)*dtcloud, &
                                                    mgncol, mgcols, .true.)


             call  micro_mg_tend (  .false.,   .false.,    .false.,   &
                                   j, jdim, mgncol,   mgcols,   nlev,    top_lev,   &
                                   dtcloud, tn,   &
                                   qvn,      ql_upd,   qi_upd,  qn_upd,   qni_upd, &
                                   relvar,   accre_enhann, &
                                   pfull,  pdel, phalf,  qa_upd, &
                                   liqcldf,  icecldf,  delta_cf, &
                                   D_eros_l4, nerosc4, D_eros_i4, nerosi4, dqcdt,  &
                                   dqidt,  crystal1, drop2,    rbar_dust_4bin,  ndust_4bin,  &
                                   ST_micro,      SQ_micro,  SL_micro,       SI_micro,       &
                                   SN_micro,      SNI_micro,                                 &
                                   surfrain(:,j), surfsnow(:,j), qsout2d_mg,  rain2d_mg,   snow2d_mg,      &
                                   qrout2d_mg,  reff_rain2d_mg, reff_snow2d_mg,  &
                                   tnd_qsnown(:,:),  tnd_nsnown(:,:), re_icen(:,:), &
                                   errstring, f_snow_berg, &
                                   Nml,  ssat_disposal, &
                                   n_diag_4d, diag_4d, diag_id, diag_pt, do_clubb=do_clubb, qcvar_clubbin=qcvar_clubb(:,j,:))

             if (trim(errstring) /= '') then
               call error_mesg ('strat_cloud/microphysics/micro_mg_tend', &
                                                         errstring, FATAL)
             endif
  endif



! ---> h1g, 2012-05-16, calculate column enthalpy and total water changes
! Note: in MG-microphysics, temperature tendency is multiplied by Cp_air.
      do i=1,idim
        do k=1,kdim
           enth_micro_col(i,j) = enth_micro_col(i,j)                    &
                + ( ST_micro(i,k) - HLV*SL_micro(i,k) - HLS*SI_micro(i,k) ) * pdel(i,k)/grav

           wat_micro_col(i,j) = wat_micro_col(i,j)                      &
               + ( SQ_micro(i,k) + SL_micro(i,k) + SI_micro(i,k) ) * pdel(i,k)/grav
        enddo

        enth_micro_col(i,j) = enth_micro_col(i,j)                                     &
                             + (-HLV*1000.0 * (surfrain(i,j)-surfsnow(i,j)) - HLS*1000.0 * surfsnow(i,j) )

        wat_micro_col(i,j) =   wat_micro_col(i,j) + surfrain(i,j) *1000.0
      enddo
! <--- h1g, 2012-05-16

      ! --- impose mass conservation ---
      if (mass_cons) then
        do i=1,idim
          m1 = 0.
          do k=1,kdim
            m1 = m1 + ( SQ_micro(i,k) + SL_micro(i,k) + SI_micro(i,k) ) * dtcloud * pdel(i,k)/grav
          end do
          m2 = 1.e3 * surfrain(i,j) * dtcloud
          if ( m2 > 1.e-12 ) then
            scalef = -m1/m2
            k=1
            if ( diag_id%rain_mass_conv > 0 ) &
            diag_4d(i,j,k,diag_pt%rain_mass_conv) = ( scalef*surfrain(i,j) - surfrain(i,j) ) * 1.e3
            if ( diag_id%snow_mass_conv > 0 ) &
            diag_4d(i,j,k,diag_pt%snow_mass_conv) = ( scalef*surfsnow(i,j) - surfsnow(i,j) ) * 1.e3
  
            surfrain(i,j) =  scalef * surfrain(i,j)
            surfsnow(i,j) =  scalef * surfsnow(i,j)
          end if
        end do
      end if

      ! surfrain from 'mmicro_pcond' is the total precipitation (rain+snow)
      surfrain(:,j) = surfrain(:,j) - surfsnow(:,j)
      ! convert units below from m/s to kg/m2
      surfrain(:,j) = surfrain(:,j)*1000.0 * dtcloud
      surfsnow(:,j) = surfsnow(:,j)*1000.0 * dtcloud

      ! update tendency due to microphysics
      SQ(:,:) = SQ(:,:) + SQ_micro(:,:)  * dtcloud
      ST(:,:) = ST(:,:) + ST_micro(:,:)/cp_air * dtcloud
      SL(:,:) = SL(:,:) + SL_micro(:,:)  * dtcloud
      SI(:,:) = SI(:,:) + SI_micro(:,:)  * dtcloud
      SN(:,:) = SN(:,:) + SN_micro(:,:)  * dtcloud
      SNI(:,:)= SNI(:,:)+ SNI_micro(:,:) * dtcloud

      ! --- cloud destruction ---
      do i=1,idim
        do k=1,kdim
          if ( (ql(i,k)+SL(i,k) .le. qmin) .and. (qi(i,k)+SI(i,k) .le. qmin) ) then
            if ( diag_id%qldt_destr > 0  .or. diag_id%ql_destr_col > 0 )  &
              diag_4d(i,j,k,diag_pt%qldt_destr) = - (ql(i,k) + SL(i,k))/dtcloud
            if ( diag_id%qidt_destr > 0  .or. diag_id%qi_destr_col > 0   ) &
              diag_4d(i,j,k,diag_pt%qidt_destr) = - (qi(i,k) + SI(i,k))/dtcloud
            if ( diag_id%qadt_destr > 0  .or. diag_id%qa_destr_col > 0   ) &
              diag_4d(i,j,k,diag_pt%qadt_destr) = - (qa(i,k) + SA(i,k))/dtcloud
            if ( diag_id%qndt_destr > 0  .or. diag_id%qn_destr_col > 0  ) &
              diag_4d(i,j,k,diag_pt%qndt_destr) = - (qn(i,k) + SN(i,k))/dtcloud
            if ( diag_id%qnidt_destr > 0   )                          &
              diag_4d(i,j,k,diag_pt%qnidt_destr)= - (qni(i,k) + SNI(i,k))/dtcloud
            if ( diag_id%qdt_destr > 0 )                              &
              diag_4d(i,j,k,diag_pt%qdt_destr)  = (ql(i,k)+SL(i,k) + qi(i,k)+SI(i,k))/dtcloud

            SQ(i,k) = SQ(i,k) + (ql(i,k)+SL(i,k)) + (qi(i,k)+SI(i,k))
            ST(i,k) = ST(i,k) - (hlv*(ql(i,k) + SL(i,k)) + hls*(qi(i,k) + SI(i,k)) )/cp_air
            SL(i,k) = SL(i,k) - (ql(i,k) + SL(i,k))
            SI(i,k) = SI(i,k) - (qi(i,k) + SI(i,k))
            SA(i,k) = SA(i,k) - (qa(i,k) + SA(i,k))
            SN(i,k) = SN(i,k) - (qn(i,k) + SN(i,k))
            SNI(i,k)= SNI(i,k)- (qni(i,k) + SNI(i,k))
          endif
        end do
      end do

      ! --- final clean up ---
      do i=1,idim
        do k=1,kdim
          if ( (ql(i,k)+SL(i,k) .le. qmin) ) then
            if ( diag_id%qdt_cleanup_liquid > 0 ) &
              diag_4d(i,j,k,diag_pt%qdt_cleanup_liquid) = (ql(i,k)+SL(i,k))/dtcloud

            SQ(i,k) = SQ(i,k) + (ql(i,k)+SL(i,k))
            ST(i,k) = ST(i,k) - (hlv*(ql(i,k)+SL(i,k)))/cp_air
            SL(i,k) = SL(i,k) - (ql(i,k)+SL(i,k))
            SN(i,k) = SN(i,k) - (qn(i,k)+SN(i,k))
            IF ( diag_id%qndt_cleanup > 0 ) &
              diag_4d(i,j,k,diag_pt%qndt_cleanup) = -(qn(i,k)+SN(i,k))/dtcloud
          endif
        end do
      end do

      do i=1,idim
        do k=1,kdim
          if ( (qi(i,k)+SI(i,k) .le. qmin) ) then
            if ( diag_id%qdt_cleanup_ice > 0 ) &
              diag_4d(i,j,k,diag_pt%qdt_cleanup_ice) = (qi(i,k)+SI(i,k))/dtcloud

            SQ(i,k) = SQ(i,k) + (qi(i,k)+SI(i,k))
            ST(i,k) = ST(i,k) - (hls*(qi(i,k)+SI(i,k)))/cp_air
            SI(i,k) = SI(i,k) - (qi(i,k)+SI(i,k))
            SNI(i,k) = SNI(i,k) - (qni(i,k)+SNI(i,k))
            if ( diag_id%qnidt_cleanup > 0 ) &
              diag_4d(i,j,k,diag_pt%qnidt_cleanup) = - (qni(i,k)+SNI(i,k))/dtcloud
          endif
        end do
      end do

      ! --- back to 3-D ---
      ST3d(:,j,:) = ST(:,:)
      SQ3d(:,j,:) = SQ(:,:)
      SL3d(:,j,:) = SL(:,:)
      SI3d(:,j,:) = SI(:,:)
      SA3d(:,j,:) = SA(:,:)
      SN3d(:,j,:) = SN(:,:)
      SNi3d(:,j,:) = SNI(:,:)

! ---> h1g, 2014-07-18, pass "f_now_berg" to moist_process 
!      for in-cloud scavenging in large-scale wet deposition
      f_snow_berg3d (:,j,:) = f_snow_berg(:,:)
! <--- h1g, 2014-07-18

      rain3d(:,j,:) = rain2d_mg(:,:)
      snow3d(:,j,:) = snow2d_mg(:,:)

      if (present(lsc_snow))      lsc_snow(:,j,:)      = qsout2d_mg(:,:)
      if (present(lsc_rain))      lsc_rain(:,j,:)      = qrout2d_mg(:,:)

      if ( do_mg_ncar_microphys) then
      ! in cldwat2m_micro.F90, the output is diameter
        if (present(lsc_snow_size)) lsc_snow_size(:,j,:) = reff_snow2d_mg(:,:)
        if (present(lsc_rain_size)) lsc_rain_size(:,j,:) = reff_rain2d_mg(:,:)
      elseif  ( do_ncar_microphys) then
      ! in MG1.5(micro_mg.F90), the output is radius
        if (present(lsc_snow_size)) lsc_snow_size(:,j,:) = 2*reff_snow2d_mg(:,:)
        if (present(lsc_rain_size)) lsc_rain_size(:,j,:) = 2*reff_rain2d_mg(:,:)
      endif

      if (diag_id%qrout > 0)     diag_4d(:,j,:,diag_pt%qrout) = qrout2d_mg(:,:)
      if (diag_id%qsout > 0)     diag_4d(:,j,:,diag_pt%qsout) = qsout2d_mg(:,:)

!------------------------------------------------------------------------
!    generate column integrated diagnostics.
!------------------------------------------------------------------------
      do nn=1, n_diag_4d
        do k =kdim,1, -1
          diag_3d(:,j,nn) = diag_3d(:,j,nn) &
                          + diag_4d(:,j,k,nn)*pdel(:,k)/grav
        enddo
      enddo

    endif !   current_total_sec >= micro_begin_sec
  end do j_loop

! ---> h1g, 06-14-2013, in order to reproduce bit-wise identical results as AM3-CLUBB
      ST3d = ST3d/dtcloud
      ST3d = ST3d*dtcloud

      SQ3d = SQ3d/dtcloud
      SQ3d = SQ3d*dtcloud

      SL3d = SL3d/dtcloud
      SL3d = SL3d*dtcloud

      SN3d = SN3d/dtcloud
      SN3d = SN3d*dtcloud

      SNI3d = SNI3d/dtcloud
      SNI3d = SNI3d*dtcloud
    
      SI3d = SI3d/dtcloud
      SI3d = SI3d*dtcloud
! <--- h1g, 06-14-2013

  call strat_netcdf (diag_id, diag_pt, diag_4d, diag_4d_kp1,  &
                     diag_3d, Time, is, js, kdim, MASK3d)

  deallocate ( diag_4d )
  deallocate ( diag_4d_kp1 )
  deallocate ( diag_3d )

end subroutine MG_microp_3D
!##############################################################################


!##############################################################################
subroutine MG_microp_3D_end()

  if (.not. module_is_initialized) return
  module_is_initialized = .false.

  call strat_netcdf_end

end subroutine MG_microp_3D_end
!##############################################################################

#else
! NULL routines return error if called but not compiled for clubb
subroutine MG_microp_3D_init(axes,Time,idim,jdim,kdim)

  ! --- calling arguments ---
  integer, intent (in)        :: axes(4)          ! x,y,z,z_half axes types
  integer, intent (in)        :: idim,jdim,kdim   ! dimensions
  type(time_type), intent(in) :: Time             ! time
    call error_mesg('MG_microp_3D_init','Not compiled with -DCLUBB',FATAL)
end subroutine MG_microp_3D_init
!##############################################################################


!##############################################################################
subroutine MG_microp_3D( Time, is, ie, js, je, lon, lat, dtcloud,              &
                         pfull3d, phalf3d, zhalf3d, LAND,                      &
                         T3d, qv3d, ql3d, qi3d, qa3d, qn3d, qni3d, ahuco3d,    &
                         dcond_ls_liquid, dcond_ls_ice,                        &
                         Ndrop_act_CLUBB, Icedrop_act_CLUBB,                   &
                         ndust, rbar_dust,                                     &
                         ST3d, SQ3d, SL3d, SI3d, SA3d, SN3d, SNi3d,            &
                         f_snow_berg3d,                                        &
                         rain3d, snow3d, surfrain, surfsnow,                   &
                         do_clubb,  qcvar_clubb, MASK3d,                       &
                         lsc_snow, lsc_rain, lsc_snow_size, lsc_rain_size )

  ! --- calling arguments ---
  type(time_type), intent (in)                         :: Time
  integer, intent (in)                                 :: is,ie,js,je
  real, intent (in),    dimension(:,:)                 :: lon,lat
  real, intent (in)                                    :: dtcloud
  real, intent (in),    dimension(:,:,:)               :: pfull3d,phalf3d
  real, intent (in),    dimension(:,:,:)               :: zhalf3d
  real, intent (in),    dimension(:,:)                 :: LAND
  real, intent (in),    dimension(:,:,:)               :: T3d,qv3d,ql3d,qi3d,qa3d,qn3d,qni3d
  real, intent (in),    dimension(:,:,:)               :: ahuco3d
  real, intent (in),    dimension(:,:,:)               :: dcond_ls_liquid,dcond_ls_ice
  real, intent (in),    dimension(:,:,:)               :: Ndrop_act_CLUBB,Icedrop_act_CLUBB
  real, intent (in),    dimension(:,:,:)               :: ndust, rbar_dust
  real, intent (out),   dimension(:,:,:)               :: ST3d,SQ3d,SL3d,SI3d,SA3d,SN3d,SNi3d
  real, intent (out),   dimension(:,:,:)               :: f_snow_berg3d  
  real, intent (out),   dimension(:,:,:)               :: rain3d,snow3d
  real, intent (out),   dimension(:,:)                 :: surfrain,surfsnow

  integer, intent (in),  optional                      :: do_clubb
  real, intent (in),  optional, dimension(:,:,:)       :: qcvar_clubb
  real, intent (in),  optional, dimension(:,:,:)       :: MASK3d
  real, intent (out), optional, dimension(:,:,:)       :: lsc_snow,      &
                                                          lsc_rain,      &
                                                          lsc_snow_size, &
                                                          lsc_rain_size
    call error_mesg('MG_microp_3D','Not compiled with -DCLUBB',FATAL)
end subroutine MG_microp_3D
!##############################################################################


!##############################################################################
subroutine MG_microp_3D_end()
    call error_mesg('MG_microp_3D_end','Not compiled with -DCLUBB',FATAL)
end subroutine MG_microp_3D_end
!##############################################################################

#endif
end module MG_microp_3D_mod
