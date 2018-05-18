module full_radiation_driver_mod

use                    mpp_mod, only: input_nml_file, mpp_pe, mpp_root_pe
use                 mpp_io_mod, only: MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE, mpp_get_field_index, mpp_get_info, mpp_read, &
                                      fieldtype, axistype, mpp_open, mpp_get_fields, mpp_get_axes, mpp_get_atts, mpp_get_axis_data
use                    fms_mod, only: open_namelist_file, close_file, stdlog, stdout, file_exist
use              constants_mod, only: SECONDS_PER_DAY, STEFAN, PI, TFREEZE, GRAV, RDGAS
use                    fms_mod, only: error_mesg, FATAL, NOTE
use           time_manager_mod, only: time_type, get_time, set_time, print_date, operator(+), operator(-)
use       press_and_geopot_mod, only: pressure_variables
use  radiation_driver_diag_mod, only: radiation_driver_diag_init, radiation_driver_diag_end
use  radiation_driver_diag_mod, only: update_rad_fields, produce_radiation_diagnostics
use radiation_driver_types_mod, only: radiation_control_type, astronomy_type, rad_output_type
use                sealw99_mod, only: sealw99_init , sealw99_time_vary , sealw99 , sealw99_endts , sealw99_end
use      solar_data_driver_mod, only: solar_data_driver_init , solar_data_driver_time_vary
use         longwave_types_mod, only: lw_output_type , lw_diagnostics_type
use        shortwave_types_mod, only: sw_output_type
use           esfsw_driver_mod, only: esfsw_driver_init , swresf , esfsw_driver_end , esfsw_number_of_bands
use              astronomy_mod, only: astronomy_init, diurnal_solar
use         sat_vapor_pres_mod, only: lookup_es
use           diag_manager_mod, only: register_diag_field, send_data
use           horiz_interp_mod, only: horiz_interp
use       aerosolrad_types_mod, only: aerosolrad_control_type
use  radiative_gases_types_mod, only: radiative_gases_type

implicit none
private

public :: full_radiation_driver_init, full_radiation_driver, full_radiation_driver_end

integer, parameter :: nlwaerb = 1, nswaerb = 1, nstoch = 1
integer :: id_qo3, qo3_file_nlev, ozone_option_i
integer :: nx, ny, nswbands, nlev, lwrad_alarm, swrad_alarm, flag_stoch = 0
type(radiation_control_type) :: Rad_control
type(rad_output_type) :: Rad_output
real, allocatable, dimension(:)         :: solflxband, qo3_file_pfull
real, allocatable, dimension(:,:)       :: pref, rad_lon, rad_lat, rad_lonb, rad_latb
real, allocatable, dimension(:,:,:)     :: deltaz, qo3, t_half, qo3_out
real, allocatable, dimension(:,:)       :: albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif
real, allocatable, dimension(:,:,:)     :: cmxolw
real, allocatable, dimension(:,:,:,:)   :: camtsw, crndlw
real, allocatable, dimension(:,:,:,:,:) :: cldsctsw, cldextsw, cldasymmsw, emrndlw, emmxolw
real, allocatable, dimension(:,:,:,:)   :: aerooptdep, aerooptdep_volc
real, allocatable, dimension(:,:,:,:)   :: aeroasymfac, aerosctopdep, aeroextopdep
real :: rrvco2, rrvch4, rrvn2o, rrvf11, rrvf12, rrvf113, rrvf22
logical :: use_ch4_for_tf_calc, use_n2o_for_tf_calc, use_co2_for_tf_calc
type(aerosolrad_control_type) :: Aerosolrad_control
type(radiative_gases_type)    :: Rad_gases

! ----------------- namelist -----------------
logical :: do_conserve_energy = .false.
real :: rch4  = 1.82120E-06
real :: rn2o  = 3.16000E-07
real :: rco2  = 3.69400E-04
real :: rf11  = 0.0
real :: rf12  = 0.0
real :: rf113 = 0.0
real :: rf22  = 0.0
logical :: do_ch4      = .true.
logical :: do_n2o      = .true.
logical :: do_co2      = .true.
logical :: do_cfc      = .false.
logical :: do_h2o_lw   = .true.
logical :: do_o3_lw    = .false.
logical :: do_co2_10um = .false.
integer :: lw_rad_time_step = 10800
integer :: sw_rad_time_step = 10800
logical :: renormalize_sw_fluxes = .true.
logical :: do_clear_sky_pass = .false.
real :: solar_constant = 1365.0
real :: albedo_sw = 0.265
character(len=16) :: ozone_option = 'from_ext_file' ! valid options are 'from_ext_file' and 'computed' 

namelist /full_radiation_driver_nml/ do_conserve_energy, rch4, rn2o, rco2, rf11, rf12, rf113, rf22, &
              do_ch4, do_n2o, do_co2, do_cfc, do_h2o_lw, do_o3_lw, do_co2_10um, lw_rad_time_step, sw_rad_time_step, &
              renormalize_sw_fluxes, do_clear_sky_pass, solar_constant, ozone_option, albedo_sw

contains
!=================================================================================================================================
subroutine full_radiation_driver_init(Time, is, ie, js, je, num_levels, axes, rad_lon_in, rad_lat_in, rad_lonb_in, rad_latb_in)
type(time_type), intent(in) :: Time
integer, intent(in) :: is, ie, js, je, num_levels, axes(:)
real,    intent(in) :: rad_lon_in(:,:), rad_lat_in(:,:), rad_lonb_in(:,:), rad_latb_in(:,:)
character(len=256) :: message
integer :: ierr, findex, num_dims, dimlens(4)
integer :: j,k, qo3_file_unit, ndims, nvar, natt, nrec, qo3_file_nlat
real, allocatable :: qo3_file_qo3(:,:), qo3_file_qo3_3d(:,:,:)
real, allocatable, dimension(:) :: lon_in, lat_in, lon_out, lat_out, dummy1, dummy2, dummy3, qo3_file_lat
type(fieldtype), allocatable, dimension(:) :: fields
type(axistype),  allocatable :: qo3_file_axes(:)

 nx = ie-is+1
 ny = je-js+1
 nlev = num_levels
 if(any((/size(rad_lon_in,1) /= nx , size(rad_lon_in,2) /= ny , size(rad_lat_in,1) /= nx , size(rad_lat_in,2) /= ny/))) then ! Remove this check later
   message = ''
   write(message(  1: 17),'(a)')    ' Extent mismatch:'
   write(message( 18: 43),'(a,i6)') ' size(rad_lon_in,1)=',size(rad_lon_in,1)
   write(message( 44: 69),'(a,i6)') ' size(rad_lon_in,2)=',size(rad_lon_in,2)
   write(message( 70: 95),'(a,i6)') ' size(rad_lat_in,1)=',size(rad_lat_in,1)
   write(message( 96:121),'(a,i6)') ' size(rad_lat_in,2)=',size(rad_lat_in,2)
   write(message(122:152),'(2(a,i6))') '  ie-is+1=',nx,' je-js+1=',ny
   call error_mesg('full_radiation_driver_init', trim(message), FATAL)
 endif

 read (input_nml_file, nml=full_radiation_driver_nml)
 write(stdlog(),nml=full_radiation_driver_nml)
 if(trim(ozone_option) == 'computed') then
   ozone_option_i = 1
 else if(trim(ozone_option) == 'from_ext_file') then
   ozone_option_i = 2
 else
   call error_mesg('full_radiation_driver_init', trim(ozone_option)//' is an invalid choice for ozone_option', FATAL)
 endif
 allocate(rad_lon(nx,ny), rad_lat(nx,ny))
 allocate(rad_lonb(nx+1,ny+1), rad_latb(nx+1,ny+1))
 rad_lon = rad_lon_in
 rad_lat = rad_lat_in
 rad_lonb = rad_lonb_in
 rad_latb = rad_latb_in
 allocate(dummy1(nlev+1), dummy2(nlev+1), pref(nlev+1,2), dummy3(nlev))
 call pressure_variables(dummy1, dummy2, pref(1:nlev,1), dummy3, 101325.)
 call pressure_variables(dummy1, dummy2, pref(1:nlev,2), dummy3,  81060.)
 deallocate(dummy1,dummy2,dummy3)
 pref(nlev+1,1) = 101325.
 pref(nlev+1,2) =  81060.

 call astronomy_init(rad_latb, rad_lonb)

 allocate(deltaz(nx,ny,nlev), qo3(nx,ny,nlev), t_half(nx,ny,nlev+1))

 lwrad_alarm = 1
 swrad_alarm = 1
 if (rch4 .eq. 0.0) do_ch4 = .false.
 if (rn2o .eq. 0.0) do_n2o = .false.
 if (rco2 .eq. 0.0) do_co2 = .false.
 if (rf11 .eq. 0.0 .and. rf12 .eq. 0.0 .and. rf22 .eq. 0.0 .and. rf113 .eq. 0.0) do_cfc = .false.
 use_ch4_for_tf_calc = do_ch4
 use_n2o_for_tf_calc = do_n2o
 use_co2_for_tf_calc = do_co2
 if (use_co2_for_tf_calc) then
   Rad_gases%co2_for_last_tf_calc = rco2
 else
   Rad_gases%co2_for_last_tf_calc = 0.0
 endif
 if (use_ch4_for_tf_calc) then
   Rad_gases%ch4_for_last_tf_calc = rch4
 else
   Rad_gases%ch4_for_last_tf_calc = 0.0
 endif
 if (use_n2o_for_tf_calc) then
   Rad_gases%n2o_for_last_tf_calc = rn2o
 else
   Rad_gases%n2o_for_last_tf_calc = 0.0
 endif
 Rad_gases%rrvco2 = rco2
 Rad_gases%rrvch4 = rch4
 Rad_gases%rrvn2o = rn2o
 Rad_gases%rrvf11 = rf11
 Rad_gases%rrvf12 = rf12
 Rad_gases%rrvf22 = rf22
 Rad_gases%rrvf113 = rf113
 call sealw99_init ( pref, do_h2o_lw, do_o3_lw, do_ch4, do_n2o, do_co2, do_co2_10um, do_cfc )

 call esfsw_driver_init

 call esfsw_number_of_bands (nswbands)
 allocate( solflxband(nswbands) )
 call solar_data_driver_init (nswbands, ierr)

 Rad_control%renormalize_sw_fluxes = renormalize_sw_fluxes
 Rad_control%do_totcld_forcing = do_clear_sky_pass
 Rad_control%using_restart_file = .false.
 Aerosolrad_control%do_aerosol = .false.  ! no aerosols in this run
 Aerosolrad_control%do_swaerosol = .false.
 Aerosolrad_control%do_lwaerosol = .false.
 Aerosolrad_control%do_swaerosol_forcing = .false.
 Aerosolrad_control%do_lwaerosol_forcing = .false.
 Aerosolrad_control%indx_swaf = 1
 Aerosolrad_control%indx_lwaf = 1
 Aerosolrad_control%volcanic_sw_aerosols = .false.
 Aerosolrad_control%volcanic_lw_aerosols = .false.
 call radiation_driver_diag_init(Time, nx, ny, nlev, axes, Rad_control, Aerosolrad_control)
 call Rad_output%alloc (nx, ny, nlev, Rad_control%do_totcld_forcing)
 call Rad_output%initvalues
 allocate(cmxolw         (nx,ny,nlev           )) ; cmxolw   = 0.0
 allocate(camtsw         (nx,ny,nlev,nstoch    )) ; camtsw   = 0.0
 allocate(crndlw         (nx,ny,nlev,nstoch    )) ; crndlw   = 0.0
 allocate(cldsctsw       (nx,ny,nlev,nswbands,1)) ; cldsctsw = 0.0
 allocate(cldextsw       (nx,ny,nlev,nswbands,1)) ; cldextsw = 0.0
 allocate(cldasymmsw     (nx,ny,nlev,nswbands,1)) ; cldasymmsw = 1.0
 allocate(emrndlw        (nx,ny,nlev,nswbands,1)) ; emrndlw    = 1.0
 allocate(emmxolw        (nx,ny,nlev,nswbands,1)) ; emmxolw    = 1.0
 allocate(aerooptdep     (nx,ny,nlev,nlwaerb   )) ; aerooptdep      = 0.0
 allocate(aerooptdep_volc(nx,ny,nlev,nlwaerb   )) ; aerooptdep_volc = 0.0
 allocate(aeroasymfac    (nx,ny,nlev,nswbands  )) ; aeroasymfac     = 0.0
 allocate(aerosctopdep   (nx,ny,nlev,nswbands  )) ; aerosctopdep    = 0.0
 allocate(aeroextopdep   (nx,ny,nlev,nswbands  )) ; aeroextopdep    = 0.0
 allocate(albedo_vis_dir (nx,ny)) ; albedo_vis_dir = albedo_sw
 allocate(albedo_nir_dir (nx,ny)) ; albedo_nir_dir = albedo_sw
 allocate(albedo_vis_dif (nx,ny)) ; albedo_vis_dif = albedo_sw
 allocate(albedo_nir_dif (nx,ny)) ; albedo_nir_dif = albedo_sw

 id_qo3 = register_diag_field('full_radiation_driver', 'qo3', axes(1:3), Time, 'Ozone mixing ratio', 'Kg ozone/Kg air')
 if(.not.file_exist('INPUT/ape_O3.nc')) then
   call error_mesg('full_radiation_driver_init', 'Ozone data file, INPUT/ape_O3.nc, does not exist', FATAL)
 endif
 call mpp_open(qo3_file_unit, 'INPUT/ape_O3.nc', action=MPP_RDONLY, form=MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)
 call mpp_get_info(qo3_file_unit,ndims,nvar,natt,nrec)
 allocate(fields(nvar))
 allocate(qo3_file_axes(ndims))
 call mpp_get_axes(qo3_file_unit, qo3_file_axes)
 call mpp_get_fields(qo3_file_unit,fields)
 findex = mpp_get_field_index(fields,'O3')
 if(findex <= 0) call error_mesg('full_radiation_driver_init', 'Field name "O3" not found in file INPUT/ape_O3.nc', FATAL)
 call mpp_get_atts(fields(findex), ndim=num_dims, siz=dimlens)
 qo3_file_nlat = dimlens(1)
 qo3_file_nlev = dimlens(2)
 allocate(qo3_file_lat(qo3_file_nlat), qo3_file_pfull(qo3_file_nlev))
 call mpp_get_axis_data(qo3_file_axes(1), qo3_file_lat)
 call mpp_get_axis_data(qo3_file_axes(2), qo3_file_pfull)
 qo3_file_lat = PI*qo3_file_lat/180.
 qo3_file_pfull = 100.*qo3_file_pfull
 allocate(qo3_file_qo3(qo3_file_nlat, qo3_file_nlev))
 call mpp_read(qo3_file_unit, fields(findex), qo3_file_qo3)
 qo3_file_qo3 = 1.0e-6*qo3_file_qo3
 allocate(qo3_file_qo3_3d(1,qo3_file_nlat, qo3_file_nlev), qo3_out(1,ny,qo3_file_nlev))
 qo3_file_qo3_3d(1,:,:) = qo3_file_qo3
 allocate(lon_in(2), lon_out(2))
 lon_in(1) = 0.0 ; lon_out(1) = 0.0 ; lon_in(2) = 2*PI ; lon_out(2) = 2*PI
 allocate(lat_in(qo3_file_nlat+1), lat_out(ny+1))
 lat_in(1) = -0.5*PI
 do j=2,qo3_file_nlat
   lat_in(j) = .5*(qo3_file_lat(j-1) + qo3_file_lat(j))
 enddo
 lat_in(qo3_file_nlat+1) = 0.5*PI
 do j=1,ny+1
   lat_out(j) = rad_latb(1,j)
 enddo
 do k=1,qo3_file_nlev
   call horiz_interp(qo3_file_qo3_3d(:,:,k), lon_in, lat_in, lon_out, lat_out, qo3_out(:,:,k))
 enddo
 deallocate(fields, qo3_file_axes, qo3_file_lat, qo3_file_qo3, qo3_file_qo3_3d, lon_in, lat_in, lon_out, lat_out)
 
end subroutine full_radiation_driver_init
!=================================================================================================================================
subroutine full_radiation_driver(Time, Time_step, p_half, p_full, tg, sphum, z_half, z_full, tdt, tdtlw, tdtsw, net_surf_sw_down, surf_lw_down)
type(time_type), intent(in) :: Time, Time_step
real, intent(in) :: p_half(:,:,:), p_full(:,:,:), tg(:,:,:), sphum(:,:,:), z_half(:,:,:), z_full(:,:,:)
real, intent(inout) :: tdt(:,:,:)
real, intent(out) :: tdtlw(:,:,:), tdtsw(:,:,:)
real, intent(out) :: net_surf_sw_down(:,:), surf_lw_down(:,:)

type(sw_output_type), dimension(1) :: Sw_output
type(lw_output_type), dimension(1) :: Lw_output
type(lw_diagnostics_type) :: Lw_diagnostics
type(astronomy_type) :: Astro_rad, Astro_phys
type(time_type) :: Time_diag
real, dimension(nx,ny) :: flux_ratio

! Sometimes the dynamic sphum can be infinitesimally negative; this causes problems in the radiation
! code.  Here we correct for this by setting any negative sphum values to be 1.0e-15.  This bug was discovered
! and fixed by Yi in the early going.  This radiation sphum is called `sphum_positive_definite` to indicate
! that it is based on the dynamic sphum, but positive definite.
real, dimension(size(sphum,1),size(sphum,2),size(sphum,3)) :: sphum_positive_definite
real :: SPHUM_MIN = 1.0e-13                                !adjust the minimum value by gc

integer :: days, seconds, time_step_int, i, j, k
character(len=24) :: label
logical :: used
 
 Time_diag = Time + Time_step
 call get_time(Time_step, seconds, days)
 time_step_int = seconds + seconds_per_day*days

 lwrad_alarm = lwrad_alarm - time_step_int
 if (lwrad_alarm <= 0) then
    Rad_control%do_lw_rad = .true.
 else
    Rad_control%do_lw_rad = .false.
 endif
 swrad_alarm = swrad_alarm - time_step_int
 if (swrad_alarm <= 0) then
    Rad_control%do_sw_rad = .true.
 else
    Rad_control%do_sw_rad = .false.
 endif
 if (Rad_control%do_lw_rad) then
    call sealw99_time_vary ( use_ch4_for_tf_calc, use_n2o_for_tf_calc, use_co2_for_tf_calc, rch4, rn2o, rco2 )
 endif
 if (Rad_control%do_sw_rad) then
    call solar_data_driver_time_vary ( Time, solar_constant, solflxband )
 endif
 if (Rad_control%do_sw_rad) then
    call diurnal_init (Astro_rad, Time-Time_step, 2*sw_rad_time_step)
 endif
 if (renormalize_sw_fluxes) then
    call diurnal_init (Astro_phys, Time, time_step_int)
 endif
 if (Rad_control%do_sw_rad .or. Rad_control%do_lw_rad) then
    label =''
    if (Rad_control%do_sw_rad) then
       label = trim(label)//' SW'
    else
       label = trim(label)//' '
    endif
    if (Rad_control%do_sw_rad .and. Rad_control%do_lw_rad) then
       label = trim(label)//'+LW'
    else if (Rad_control%do_lw_rad) then
       label = trim(label)//' LW'
    else
       label = trim(label)//' '
    endif
!  comment out by gc
!    call print_date(Time, label(1:15)//' Date = ')
 endif

 qo3 = get_qo3(z_full, p_full)
 if(id_qo3 > 0) used = send_data(id_qo3, qo3, Time)
 t_half = compute_t_half(tg,z_half,z_full)
 deltaz = compute_deltaz(z_half,z_full)

 ! [SKC] Populate sphum_positive_definite, and replace sphum with sphum_positive_definite in all sealw99 and swresf function calls.
 do k=1,size(sphum,3)
   do j=1,size(sphum,2)
     do i=1,size(sphum,1)
       if(sphum(i,j,k) .le. SPHUM_MIN) then
         sphum_positive_definite(i,j,k) = SPHUM_MIN
       else
         sphum_positive_definite(i,j,k) = sphum(i,j,k)
       endif
     enddo
   enddo
 enddo

 if (Rad_control%do_lw_rad) then
    call Lw_output(1)%alloc (nx, ny, nlev, Rad_control%do_totcld_forcing)
    call sealw99(p_full, p_half, tg, t_half, sphum_positive_definite, deltaz, qo3, rco2, rf11, rf12, rf113, rf22, &
         emrndlw, emmxolw, crndlw, cmxolw, aerooptdep, aerooptdep_volc, Lw_output(1), Lw_diagnostics, &
         flag_stoch, Rad_control%do_totcld_forcing, Aerosolrad_control%do_lwaerosol, &
         Aerosolrad_control%volcanic_lw_aerosols)
    Rad_output%tdtlw = Lw_output(1)%heatra/SECONDS_PER_DAY
    Rad_output%flxnet = Lw_output(1)%flxnet
    if (do_clear_sky_pass) then
       Rad_output%tdtlw_clr = Lw_output(1)%heatracf/SECONDS_PER_DAY
       Rad_output%flxnetcf = Lw_output(1)%flxnet
    endif
 endif
 if (Rad_control%do_sw_rad) then
    call Sw_output(1)%alloc (nx, ny, nlev, Rad_control%do_totcld_forcing)
    call swresf(p_full, p_half, tg, sphum_positive_definite, deltaz, albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif, qo3, rco2, rch4,&
                rn2o, solflxband, solar_constant, Astro_rad%rrsun, Astro_rad%cosz, Astro_rad%fracday, camtsw, cldsctsw, cldextsw,&
                cldasymmsw, aeroasymfac, aerosctopdep, aeroextopdep, Rad_control%do_totcld_forcing, flag_stoch, Sw_output(1))
    Rad_output%tdtsw = Sw_output(1)%hsw/SECONDS_PER_DAY
    Rad_output%ufsw = Sw_output(1)%ufsw
    Rad_output%dfsw = Sw_output(1)%dfsw
    if (do_clear_sky_pass) then
       Rad_output%tdtsw_clr = Sw_output(1)%hswcf/SECONDS_PER_DAY
       Rad_output%ufsw_clr = Sw_output(1)%ufswcf(:,:,:)
       Rad_output%dfsw_clr = Sw_output(1)%dfswcf
       Rad_output%flux_sw_down_total_dir_clr = Sw_output(1)%dfsw_dir_sfc_clr
       Rad_output%flux_sw_down_total_dif_clr = Sw_output(1)%dfsw_dif_sfc_clr
       Rad_output%flux_sw_down_vis_clr = Sw_output(1)%dfsw_vis_sfc_clr
    endif
 endif
 Rad_output%tdt_rad = (Rad_output%tdtsw+Rad_output%tdtlw)
 if (do_clear_sky_pass) then
    Rad_output%tdt_rad_clr = (Rad_output%tdtsw_clr + Rad_output%tdtlw_clr)
 endif
 if (Rad_control%do_lw_rad) then
    Rad_output%flux_lw_surf = STEFAN*t_half(:,:,nlev+1)**4 - Lw_output(1)%flxnet(:,:,nlev+1)
 endif
 if (Rad_control%do_sw_rad) then
    Rad_output%flux_sw_surf = Sw_output(1)%dfsw(:,:,nlev+1) - Sw_output(1)%ufsw(:,:,nlev+1)
    Rad_output%flux_sw_surf_dir = Sw_output(1)%dfsw_dir_sfc
    Rad_output%flux_sw_surf_refl_dir = Sw_output(1)%ufsw_dir_sfc
    Rad_output%flux_sw_surf_dif = Sw_output(1)%dfsw_dif_sfc - Sw_output(1)%ufsw_dif_sfc
    Rad_output%flux_sw_down_vis_dir = Sw_output(1)%dfsw_vis_sfc_dir
    Rad_output%flux_sw_down_vis_dif = Sw_output(1)%dfsw_vis_sfc_dif
    Rad_output%flux_sw_down_total_dir = Sw_output(1)%dfsw_dir_sfc
    Rad_output%flux_sw_down_total_dif = Sw_output(1)%dfsw_dif_sfc
    Rad_output%flux_sw_vis = Sw_output(1)%dfsw_vis_sfc - Sw_output(1)%ufsw_vis_sfc
    Rad_output%flux_sw_vis_dir = Sw_output(1)%dfsw_vis_sfc_dir
    Rad_output%flux_sw_refl_vis_dir = Sw_output(1)%ufsw_vis_sfc_dir
    Rad_output%flux_sw_vis_dif = Sw_output(1)%dfsw_vis_sfc_dif - Sw_output(1)%ufsw_vis_sfc_dif
 endif
 call update_rad_fields(1,nx,1,ny, Time_diag, Astro_rad, Astro_phys, Rad_control, Sw_output, Rad_output, flux_ratio)
 call produce_radiation_diagnostics(1,nx,1,ny, Time_diag, Time, rad_lat, pref(nlev+1,1), t_half(:,:,nlev+1), &
                                  p_half, p_half, albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif, &
                                  flux_ratio, Astro_rad, Astro_phys, Rad_output, Rad_gases, Rad_control, Lw_output, Sw_output)

 tdt = tdt + Rad_output%tdt_rad
 tdtlw = Rad_output%tdtlw
 tdtsw = Rad_output%tdtsw
 if (Rad_control%do_sw_rad) then
   net_surf_sw_down = Rad_output%flux_sw_surf
 endif
 if (Rad_control%do_lw_rad) then
   surf_lw_down = Rad_output%flux_lw_surf
 endif

 if (Rad_control%do_lw_rad) then
    call Lw_output(1)%dealloc
    call Lw_diagnostics%dealloc
 endif
 if (Rad_control%do_sw_rad) then
    call Sw_output(1)%dealloc
    call Astro_rad%dealloc
 endif
 call Astro_phys%dealloc
 use_ch4_for_tf_calc = .false.
 use_n2o_for_tf_calc = .false.
 use_co2_for_tf_calc = .false.
 if (Rad_control%do_lw_rad) then
    call sealw99_endts
    lwrad_alarm = lwrad_alarm + lw_rad_time_step
    Rad_control%do_lw_rad = .false.
 endif
 if (Rad_control%do_sw_rad) then
    swrad_alarm = swrad_alarm + sw_rad_time_step
    Rad_control%do_sw_rad = .false.
 endif

end subroutine full_radiation_driver
!=================================================================================================================================
subroutine full_radiation_driver_end

 call Rad_output%dealloc
 call radiation_driver_diag_end (Rad_control)
 call sealw99_end
 call esfsw_driver_end
 deallocate(solflxband)
 deallocate(pref, rad_lon, rad_lat, rad_lonb, rad_latb)
 deallocate(deltaz, qo3, t_half)
 deallocate(albedo_vis_dir, albedo_nir_dir, albedo_vis_dif, albedo_nir_dif)
 deallocate(cmxolw)
 deallocate(camtsw, crndlw)
 deallocate(cldsctsw, cldextsw, cldasymmsw, emrndlw, emmxolw)
 deallocate(aerooptdep, aerooptdep_volc)
 deallocate(aeroasymfac, aerosctopdep, aeroextopdep)

end subroutine full_radiation_driver_end
!=================================================================================================================================
 function get_qo3(z_full, p_full) result(ozone)
 real, intent(in) :: z_full(:,:,:), p_full(:,:,:)
 real :: ozone(size(z_full,1), size(z_full,2), nlev)
 real:: ozone_max = 10.2E-06, z0 = 36.0E3, width = 13.0E3
 integer :: i,j,k,L

 if(ozone_option_i == 1) then
    do k=1,nlev
      ozone(:,:,k) = ozone_max*exp(-((z_full(:,:,k)-z0)/width)**2)
    enddo
 else if(ozone_option_i == 2) then
    L = 1
    do k=1,nlev
    do j=1,size(z_full,2)
    do i=1,size(z_full,1)
      if(p_full(i,j,k) >= qo3_file_pfull(qo3_file_nlev)) then
         ozone(i,j,k) = qo3_out(1,j,qo3_file_nlev)            
      else if(p_full(i,j,k) <= qo3_file_pfull(1)) then
         ozone(i,j,k) = qo3_out(1,j,1)            
      else
         if(p_full(i,j,k) < qo3_file_pfull(L) .or. p_full(i,j,k) > qo3_file_pfull(L+1)) then
            do L=1,qo3_file_nlev
               if(p_full(i,j,k) >= qo3_file_pfull(L) .and. p_full(i,j,k) <= qo3_file_pfull(L+1)) exit
            enddo  
         endif
         ozone(i,j,k) = ((p_full(i,j,k)-qo3_file_pfull(L))*qo3_out(1,j,L+1) + (qo3_file_pfull(L+1)-p_full(i,j,k))*qo3_out(1,j,L))/(qo3_file_pfull(L+1)-qo3_file_pfull(L))
      endif
    enddo
    enddo
    enddo
 endif

 end function get_qo3
!=================================================================================================================================
 function compute_t_half(tg,z_half,z_full) result(t_half)
 real, intent(in) :: tg(:,:,:), z_half(:,:,:), z_full(:,:,:)
 real :: t_half(size(z_half,1), size(z_half,2), nlev+1)
 integer :: k

 t_half(:,:,1) = tg(:,:,1)
 do k=2,nlev
   t_half(:,:,k) = ((z_half(:,:,k)-z_full(:,:,k-1))*tg(:,:,k) + (z_full(:,:,k)-z_half(:,:,k))*tg(:,:,k-1))/(z_full(:,:,k)-z_full(:,:,k-1))
 enddo
 t_half(:,:,nlev+1) = tg(:,:,nlev)

 end function compute_t_half
!=================================================================================================================================
 function compute_deltaz(z_half,z_full) result(deltaz)
 real, intent(in) :: z_half(:,:,:), z_full(:,:,:)
 real :: deltaz(size(z_full,1), size(z_full,2), nlev)
 integer :: k

 deltaz(:,:,1) = z_full(:,:,1) - z_half(:,:,2) 
 do k=2,nlev
   deltaz(:,:,k) = z_half(:,:,k) - z_half(:,:,k+1)
 enddo

 end function compute_deltaz
!=================================================================================================================================
!subroutine state_init(p_half, pres, pflux, temp, tflux, rh2o, qo3, deltaz)

!real, intent(in), dimension(:,:,:) :: p_half
!real, intent(out), dimension(:,:,:) :: pres, pflux, temp, tflux, rh2o, qo3, deltaz
!real, dimension(8) :: h_ref, p_ref, t_ref
!real, dimension(7) :: lapse_ref
!real :: dhght, tmean, rtog, rgamog
!real, dimension(size(p_half,1),size(p_half,2),nlev+1) :: hght
!real, dimension(size(p_half,1),size(p_half,2)) :: lnptop, lnpbot, tempin, esat
!integer :: i, j, k, kr
!integer :: iunit
!real :: qo3k
!real :: log_p_at_top=2.0 
!character(len=64) :: err
!logical :: verbose = .false. 

!h_ref(1:8) = (/ 0., 11.e3, 20.e3, 32.e3, 47.e3, 51.e3, 71.e3, 110.e3 /)
!lapse_ref(1:7) = (/ -6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3 /)
!p_ref(1) = 101325.
!t_ref(1) = 15. + TFREEZE
!do k = 2, 8
!  dhght = h_ref(k) - h_ref(k-1)
!  t_ref(k) = t_ref(k-1) + dhght * lapse_ref(k-1)
!  tmean = 0.5*(t_ref(k-1)+t_ref(k))
!  p_ref(k) = p_ref(k-1) * exp(-dhght*GRAV/(RDGAS*tmean))
!enddo
!lnptop = log(p_half(:,:,1))
!do k = 1, nlev
!  lnpbot = log(p_half(:,:,k+1))
!  p_ref(k) = p_ref(k-1) * exp(-dhght*GRAV/(RDGAS*tmean))
!  pres(:,:,k) = (p_half(:,:,k+1)-p_half(:,:,k)) / (lnpbot - lnptop)
!  lnptop = lnpbot
!enddo

!do j = 1, ny
!do i = 1, nx
!do k = 1, nlev
!  if (pres(i,j,k) .ge. p_ref(1)) then
!     kr = 2
!  else if (pres(i,j,k) .le. p_ref(8)) then
!     kr = 7
!  else
!     do kr = 2, 8
!     if (pres(i,j,k) .lt. p_ref(kr-1) .and. pres(i,j,k) .ge. p_ref(kr)) then
!        exit
!     endif
!     enddo
!  endif
!  rgamog = -RDGAS*lapse_ref(kr-1)/GRAV
!  temp(i,j,k) = t_ref(kr-1) * (pres(i,j,k)/p_ref(kr-1)) ** rgamog
!enddo
!enddo
!enddo
!if (do_conserve_energy) then 
!   pflux = p_half
!else
!   pflux(:,:,1) = 0.0
!   do k = 2, nlev
!     pflux(:,:,k) = 0.5*(pres(:,:,k-1)+pres(:,:,k))
!   enddo
!   pflux(:,:,nlev+1) = p_half(:,:,nlev+1)
!endif
!tflux(:,:,nlev+1) = t_ref(1) * (p_half(:,:,nlev+1)/p_ref(1)) ** (-RDGAS*lapse_ref(1)/GRAV)
!do k = 2, nlev
!  tflux(:,:,k) = 0.5*(temp(:,:,k-1)+temp(:,:,k))
!enddo
!tflux(:,:,1) = temp(:,:,1) 
!do k = nlev, 1, -1
!  tempin = temp(:,:,k)
!  where (pres(:,:,k) .lt. p_ref(2)) tempin = t_ref(2)
!  call lookup_es(tempin, esat, err)
!  if (err .ne. '') then
!     call error_mesg('solo_rad_util_mod',err,FATAL)
!  end if
!  rh2o(:,:,k) = 0.50 * (0.622/pres(:,:,k)) * esat(:,:)
!  if (k .lt. nlev) then
!     where(rh2o(:,:,k) .gt. rh2o(:,:,k+1)) rh2o(:,:,k) = rh2o(:,:,k+1)
!  end if
!end do
!deltaz(:,:,1) = log_p_at_top*RDGAS*temp(:,:,1)*(1.0+0.608*rh2o(:,:,1))/GRAV
!do k = 2, nlev
!  deltaz(:,:,k) = log(pflux(:,:,k+1)/pflux(:,:,k))*RDGAS* temp(:,:,k)*(1.0+0.608*rh2o(:,:,k))/GRAV
!enddo
!qo3(:,:,1:nlev) = 0.0
!if (file_exist('INPUT/ozone.txt')) then
!   iunit = open_namelist_file ('INPUT/ozone.txt')
!   do k = 1, nlev
!     read (iunit,*,end=10) qo3k
!     qo3(:,:,k) = qo3k
!   enddo
!   10 call close_file (iunit)
!   call error_mesg('solo_rad_util_mod', 'Reading ozone profile from file INPUT/ozone.txt', NOTE)
!endif
!if (verbose) then
!   write(stdout(),9010)
!   do k = 1, nlev
!     write(stdout(),9011) k, pres(1,1,k), temp(1,1,k), rh2o(1,1,k)*1000., deltaz(1,1,k)
!   enddo
!   write(stdout(),9012) pflux(1,1,nlev+1), tflux(1,1,nlev+1)
!endif
!9010 format(2x,'k',5x,'pres',8x,'temp',8x,'rh2o',9x,'dz')
!9011 format(i3,4f12.2)
!9012 format(3x,2f12.2)
!end subroutine state_init
!=================================================================================================================================
subroutine diurnal_init ( Astro, Time, dt )
type(astronomy_type), intent(inout) :: Astro
type(time_type), intent(in) :: Time
integer, intent(in) :: dt

 call Astro%alloc(size(rad_lon,1), size(rad_lon,2))
 call diurnal_solar (rad_lat, rad_lon, Time, Astro%cosz, Astro%fracday, Astro%rrsun, dt_time=set_time(dt,0))
 Astro%fracday = min(Astro%fracday, 1.0)
 Astro%solar = Astro%cosz*Astro%fracday*Astro%rrsun

end subroutine diurnal_init
!=================================================================================================================================
end module full_radiation_driver_mod
