module tracer_source_sink_mod

use      constants_mod, only: RDGAS, RVGAS, GRAV
use            fms_mod, only: error_mesg, FATAL
use sat_vapor_pres_mod, only: sat_vapor_pres_init, lookup_es
use    tracer_type_mod, only: tracer_type

implicit none

public :: radon_source_sink, water_vapor_source_sink

real    :: rratio
logical :: sat_vapor_is_initialized = .false.

contains

!------------------------------------------------------------------------------------------------------
subroutine water_vapor_source_sink(attributes, tg, dp, p_full, delta_t, mix_rat, condense, rh, rain, evap)
type(tracer_type), intent(inout) :: attributes
real, intent(in), dimension(:,:,:) :: tg, dp, p_full
real, intent(in) :: delta_t
real, intent(inout), dimension(:,:,:) :: mix_rat
real, intent(out),   dimension(:,:,:) :: condense, rh
real, intent(out),   dimension(:,:)   :: rain, evap
real, dimension(size(tg,1),size(tg,2),size(tg,3)) :: esat, sat_mix_rat, mix_rat_init, k_evap
real :: rh_crit, evap_time, x1
integer :: i, j, k, evap_layers, num_levels
character(len=256) :: err_msg
character(len=4) :: cstring

  if (.not.sat_vapor_is_initialized) then
    call sat_vapor_pres_init
    rratio = RDGAS/RVGAS
    sat_vapor_is_initialized = .true.
  endif

   call lookup_es( tg, esat, err_msg )
   if(len_trim(err_msg) /= 0) call error_mesg('water_vapor_source_sink', trim(err_msg), FATAL)
   sat_mix_rat = rratio*esat/(p_full-esat)

   mix_rat_init = mix_rat
   mix_rat = min(mix_rat_init, sat_mix_rat)
   condense = mix_rat_init - mix_rat
   num_levels = size(tg,3)
   rain = 0.0
   do k =1,num_levels        
     rain = rain + condense(:,:,k)*dp(:,:,k)/GRAV
   enddo
   rain = rain/delta_t

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
          k_evap(i,j,k) = (rh_crit*sat_mix_rat(i,j,k) - mix_rat(i,j,k))*x1
          mix_rat(i,j,k) = mix_rat(i,j,k) + k_evap(i,j,k)
       enddo
     enddo
   enddo
     
   evap = 0.0
   do k = num_levels+1-evap_layers, num_levels
     evap = evap + (dp(:,:,k)/GRAV)*k_evap(:,:,k)
   enddo
   evap = evap/delta_t
   rh = mix_rat/sat_mix_rat  

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
