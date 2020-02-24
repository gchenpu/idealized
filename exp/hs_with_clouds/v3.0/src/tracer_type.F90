module tracer_type_mod

#ifdef INTERNAL_FILE_NML
   use mpp_mod, only: input_nml_file
#else
   use fms_mod, only: open_namelist_file
#endif

use fms_mod, only: write_version_number, check_nml_error, close_file, error_mesg, NOTE, lowercase, stdlog

implicit none
private

public :: tracer_type, NO_TRACER, get_tracer_index, tracer_type_init, tracer_type_end

character(len=128) :: version = '$Id: tracer_type.F90,v 11.0.6.1 2014/12/12 17:38:54 pjp Exp $'
character(len=128) :: tagname = '$Name: post_ulm_pjp $'

type tracer_type
  character(len=32) :: name, longname, units, numerical_representation, advect_horiz, advect_vert, hole_filling
  real :: robert_coeff, flux, sink
end type

logical :: module_is_initialized = .false.
integer, parameter :: NO_TRACER = 0
integer, parameter :: max_num_tracers = 100
integer :: num_tracers
character(len=128), dimension(max_num_tracers) :: name, longname, units, numerical_representation, advect_vert, hole_filling
real,               dimension(max_num_tracers) :: robert_coeff, radon_flux, radon_sink

namelist /tracer_type_nml/ num_tracers, name, longname, units, numerical_representation, &
                           advect_vert, hole_filling, robert_coeff, radon_flux, radon_sink

contains
!---------------------------------------------------------------------------------------------------------------------------------
subroutine tracer_type_init(tracer_attributes)
type(tracer_type), allocatable, intent(out) :: tracer_attributes(:)
integer :: ierr, io, ntr

if(module_is_initialized) return
call write_version_number(version, tagname)

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=tracer_type_nml, iostat=io)
    ierr = check_nml_error(io, 'tracer_type_nml')
#else
    unit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read(unit, nml=tracer_type_nml, iostat=io, end=20)
      ierr = check_nml_error (io, 'tracer_type_nml')
    enddo
20  call close_file (unit)
#endif
write (stdlog(), nml=tracer_type_nml)

allocate(tracer_attributes(num_tracers))

do ntr=1,num_tracers
  tracer_attributes(ntr)%name     = lowercase(name(ntr))
  tracer_attributes(ntr)%longname = lowercase(longname(ntr))
  tracer_attributes(ntr)%units    = lowercase(units(ntr))
  tracer_attributes(ntr)%numerical_representation = numerical_representation(ntr)
  tracer_attributes(ntr)%advect_vert  = advect_vert(ntr)
  tracer_attributes(ntr)%hole_filling = hole_filling(ntr)
  tracer_attributes(ntr)%robert_coeff = robert_coeff(ntr)
  if(trim(name(ntr)) == 'radon') then
    tracer_attributes(ntr)%flux = radon_flux(ntr)
    tracer_attributes(ntr)%sink = radon_sink(ntr)
  else
    tracer_attributes(ntr)%flux = 0.0
    tracer_attributes(ntr)%sink = 0.0
  endif
    
  if(trim(tracer_attributes(ntr)%numerical_representation) == 'grid') then
    if(trim(tracer_attributes(ntr)%hole_filling) == 'on') then
      call error_mesg('tracer_type_init','Warning: hole_filling scheme = on will be ignored for grid tracer '// &
      tracer_attributes(ntr)%name,NOTE)
    endif
  endif
enddo

module_is_initialized = .true.

end subroutine tracer_type_init
!---------------------------------------------------------------------------------------------------------------------------------
integer function get_tracer_index(tracers,name)
type(tracer_type), intent(in) :: tracers(:)
character(len=*), intent(in)  :: name
integer :: ntr

get_tracer_index = NO_TRACER
do ntr=1,size(tracers)
  if(trim(tracers(ntr)%name) == trim(name)) then
    get_tracer_index = ntr
    exit
  endif
enddo

end function get_tracer_index
!---------------------------------------------------------------------------------------------------------------------------------
subroutine tracer_type_end(tracer_attributes)
type(tracer_type), allocatable, intent(out) :: tracer_attributes(:)

deallocate(tracer_attributes)
num_tracers = 0

end subroutine tracer_type_end
!---------------------------------------------------------------------------------------------------------------------------------
end module tracer_type_mod
