module tracer_source_sink_mod

use      constants_mod, only: RDGAS, RVGAS, GRAV, hlv, hls, cp_air
use            fms_mod, only: error_mesg, FATAL
use sat_vapor_pres_mod, only: sat_vapor_pres_init, lookup_es, compute_qs
use    tracer_type_mod, only: tracer_type
use       time_manager_mod, only: time_type
use strat_cloud_mod,            only: strat_cloud
use beta_dist_mod,        only: beta_dist_init, beta_dist_end, incomplete_beta 

implicit none

public :: radon_source_sink, water_vapor_source_sink
public :: beta_dist_init, incomplete_beta, beta_dist_end

real    :: rratio
logical :: sat_vapor_is_initialized = .false.

contains

!------------------------------------------------------------------------------------------------------
subroutine water_vapor_source_sink(Time, is, ie, js, je, attributes, tg, dp, p_full, p_half, delta_t, mix_rat, ql, qi, qa, condense, evap_rate, rh, rh2, rain, evap)
type(time_type), intent(in) :: Time
integer, intent(in) :: is, ie,js, je
type(tracer_type), intent(inout) :: attributes
real, intent(in), dimension(:,:,:) :: tg, dp, p_full, p_half
real, intent(in) :: delta_t
real, intent(inout), dimension(:,:,:) :: mix_rat, ql, qi, qa
real, intent(out),   dimension(:,:,:) :: condense, evap_rate, rh, rh2
real, intent(out),   dimension(:,:)   :: rain, evap
real, dimension(size(tg,1),size(tg,2),size(tg,3)) :: esat, sat_mix_rat, mix_rat_init, k_evap, tg2
real :: rh_crit, evap_time, x1
integer :: i, j, k, evap_layers, num_levels
character(len=256) :: err_msg
character(len=4) :: cstring

!yim: variables used by strat_cloud
real, dimension(size(tg,1),size(tg,2),size(tg,3)) :: omega, mc_full, &
                                                  diff_t,  convective_humidity_ratio,  &
                                                  convective_humidity_area, &
                                                  radturbten
logical                            :: limit_conv_cloud_frac = .false.
real, dimension(size(tg,1),size(tg,2)) :: land
real, dimension(size(tg,1),size(tg,2),size(tg,3)) :: ST, SQ, SL, SI, SA,  &
                                                  rain3d, snow3d,  &
                                                  snowclr3d, f_snow_berg
real, dimension(size(tg,1),size(tg,2)) :: rain2, snow
!  logical         :: module_is_initialized = .false.
  if (.not.sat_vapor_is_initialized) then
    call sat_vapor_pres_init
    rratio = RDGAS/RVGAS
    sat_vapor_is_initialized = .true.
  endif


!if (.not.module_is_initialized) then
!    call beta_dist_init
!endif

!yim: the original condensation scheme 
!isaac's way to compute es
!   call lookup_es( tg, esat, err_msg )
!   if(len_trim(err_msg) /= 0) call error_mesg('water_vapor_source_sink', trim(err_msg), FATAL)
!   sat_mix_rat = rratio*esat/(p_full-esat)
!yim's way (to be consistent with the pdf cloud scheme in strat_cloud_legacy.F90
    call compute_qs( tg-((hlv*ql+hls*qi)/cp_air), p_full, sat_mix_rat)
!yim's perturbation case 1
!    sat_mix_rat = 1.14*sat_mix_rat
!yim's perturbation case 2
!    sat_mix_rat = 1.14*sat_mix_rat*(293.15/tg)**3.

!   mix_rat_init = mix_rat
!   mix_rat = min(mix_rat_init, sat_mix_rat)
!   condense = mix_rat_init - mix_rat
   num_levels = size(tg,3)
!   rain = 0.0
!   do k =1,num_levels        
!     rain = rain + condense(:,:,k)*dp(:,:,k)/GRAV
!   enddo
!   rain = rain/delta_t

!yim: call cloud scheme
!radturbten omega mc_full diff_t land ttnd qtnd, q_tnd(3), f_snow_berg rain3d snow3d snowclr3d
!rain snow, convective_humidity_ratio, convective_humidity_area limit_conv_cloud_frac, mask
!initialize
  radturbten = 0.
  omega = 0.
  mc_full = 0.
  diff_t = 0.
  land = 0.
  ST = 0.
  SQ = 0.
  SL = 0.
  SI = 0.
  SA = 0.
  f_snow_berg = 0.
  rain3d = 0.
  snow3d = 0.
  snowclr3d = 0.
  rain2 = 0.
  snow = 0.
  convective_humidity_ratio = 0.
  convective_humidity_area = 0.

!yim perturbation increase temperature by 2K

!  tg2 = tg+2.

            call strat_cloud (Time, is, ie, js, je, delta_t, p_full, p_half,   & 
                          radturbten, tg, mix_rat, ql,        &
                          qi, qa,      &
                          omega, mc_full, diff_t, land, ST, SQ,    &
                          SL, SI,   &
                          SA, f_snow_berg,  &
                          rain3d, snow3d, snowclr3d,  &
                          rain2, snow, convective_humidity_ratio,  &
                          convective_humidity_area, &
                          limit_conv_cloud_frac)
!yim: update the cloud fields, intentionally leaving out ST
   qa=qa+SA
   ql=ql+SL
   qi=qi+SI
   mix_rat=mix_rat+SQ
   condense = SQ
   rain = (rain2+snow)/delta_t

   rh_crit = 1.0
   evap_time = 0.5*86400

   ! compute evaporation source at each level and then sum to compute total evaporation
   ! number of levels feeling evaporation is evap_layers
   ! evaporation is computed by relaxing exponentially to the critical relative humidity
    
   evap_layers = 3

   ! This evaporation scheme is crude: the user has to adjust it when the vertical resolution is changed. But I think it is
   ! simpler and safer this way for the time being -- it could get confusing when using generalized sigma-p vertical coordinate.
   if(num_levels /= 20) then
     write(cstring,'(i4)') num_levels
     call error_mesg('water_vapor_source_sink','Evaporation scheme is intended for use with 20 equally spaced pressure levels.'//&
                     ' You are using '//trim(cstring)//' levels', FATAL)
   endif
    
   x1 = 1-exp(-delta_t/evap_time)
   do j = 1, size(tg,2)
     do i = 1, size(tg,1)
       do k = num_levels+1-evap_layers, num_levels 
!mef: noEVAP_3030
!if ((js+j-1).gt.21.and.(js+j-1).lt.44) then
!  k_evap(i,j,k) = 0.0
!else
!mef: noEVAP_3090
!if ((js+j-1).lt.22) then
!  k_evap(i,j,k) = 0.0
!elseif ((js+j-1).gt.43) then
!  k_evap(i,j,k) = 0.0
!else
!mef: noEVAP_1010
!if ((js+j-1).gt.28.and.(js+j-1).lt.37) then
!  k_evap(i,j,k) = 0.0
!else
!mef: noEVAP_1020
!if ((js+j-1).gt.24.and.(js+j-1).lt.29) then
!  k_evap(i,j,k) = 0.0
!elseif ((js+j-1).gt.36.and.(js+j-1).lt.41) then
!  k_evap(i,j,k) = 0.0
!else
!mef: noEVAP_NH
!if ((js+j-1).gt.32) then
!  k_evap(i,j,k) = 0.0
!else
          k_evap(i,j,k) = (rh_crit*sat_mix_rat(i,j,k) - mix_rat(i,j,k))*x1
          !k_evap(i,j,k) = 0.0 !mef: noEVAP
          mix_rat(i,j,k) = mix_rat(i,j,k) + k_evap(i,j,k)
!endif
       enddo
     enddo
   enddo
     
   evap = 0.0
   evap_rate = 0.0
   do k = num_levels+1-evap_layers, num_levels
     evap = evap + (dp(:,:,k)/GRAV)*k_evap(:,:,k)
     evap_rate(:,:,k) = k_evap(:,:,k)/delta_t
   enddo
   evap = evap/delta_t
   rh = mix_rat/sat_mix_rat
   rh2 = (mix_rat+ql+qi)/sat_mix_rat

end subroutine water_vapor_source_sink
!------------------------------------------------------------------------------------------------------
subroutine radon_source_sink(attributes, radon, p_half, delta_t, dt_radon)
type(tracer_type), intent(inout) :: attributes
real, intent(in) :: radon(:,:,:), p_half(:,:,:)
real, intent(in) :: delta_t
real, intent(inout) :: dt_radon(:,:,:)
real :: rdamp, mass_air
integer :: i,j,k,kbot

if(attributes%flux == 0.0 .and. attributes%sink == 0.0) return

if(attributes%sink < 0.0) then ! negative value of sink is a flag that it's units are days
  rdamp = -1.0/(86400.*attributes%sink)
else if(attributes%sink == 0.0) then ! zero is a flag that there is no sink
  rdamp = 0.0
else
  rdamp = 1.0/attributes%sink ! sink is in seconds
endif

kbot = size(radon,3)

do j=1,size(radon,2)
do i=1,size(radon,1)
  mass_air = (p_half(i,j,kbot+1) - p_half(i,j,kbot))/GRAV
  dt_radon(i,j,kbot) = dt_radon(i,j,kbot) + (attributes%flux/mass_air - rdamp*radon(i,j,kbot))/(1.0 + delta_t*rdamp)
enddo
enddo

do k=1,kbot-1
do j=1,size(radon,2)
do i=1,size(radon,1)
  dt_radon(i,j,k) = dt_radon(i,j,k) - rdamp*radon(i,j,k)/(1.0 + delta_t*rdamp)
enddo
enddo
enddo

end subroutine radon_source_sink
!------------------------------------------------------------------------------------------------------

end module tracer_source_sink_mod
