Program main
! The current code is designed to work on 1 processor, but can be extended to multiple processors.

! netcdf i/o 
  use  netcdf_file_mod, only:  define_axis,  define_time_axis,  &
                         define_variable, define_variable_tave,  &
                         define_tave,  define_file,             &
                         read_file_header,  write_file_header,  &
                         copy_file_header, get_axis_units,   &
                         get_axis_long_name, get_axis_cart_name,  &
                         get_axis_edges_name,               &
                         get_var_units, get_var_long_name,  &
                         get_missing_value, get_axis_direction,  &
                         get_calendar_type, get_valid_range
  use  ncread_write_mod, only:  write_time_coord, write_variable_tave, &
                                write_variable, read_time_coord,  &
                                read_variable, read_variable_tave, &
                                axis_values, time_axis_values,  &
                                variable_tave_info
  use  ncfile_access_mod, only:  file_name, num_axis, num_var,       &
                                 axis_name, axis_index, axis_length, &
                                 time_axis_index, time_axis_length,  &
                                 var_name, var_index, num_var_axis,  &
                                 var_axis_len, var_axis
  use       ncfile_mod, only:  ncfile_type, close_ncfile, empty_ncfile
  use   ncd_define_mod, only:  float_type, double_type

! spectral transform
 use         constants_mod, only: constants_init, radius, omega, pi, kappa, grav, RdGAS
 use               fms_mod, only: fms_init, fms_end, mpp_pe, mpp_root_pe, mpp_npes, &
                                  error_mesg, NOTE, FATAL, write_version_number, stdlog, &
				  lowercase, uppercase, check_nml_error, file_exist, &
				  set_domain, read_data, write_data, nullify_domain, &
                                  open_namelist_file, open_restart_file, open_file, close_file
 use            fms_io_mod, only: fms_io_init, fms_io_exit
 use        transforms_mod, only: transforms_init, transforms_end, &
                                  vor_div_from_uv_grid, uv_grid_from_vor_div,&
                                  trans_grid_to_spherical, trans_spherical_to_grid, &
                                  get_grid_domain, get_spec_domain, &
                                  grid_domain, spectral_domain,     &
				  compute_lat_deriv_cos, divide_by_cos, compute_lon_deriv_cos, compute_laplacian, &
				  get_eigen_laplacian, horizontal_advection, &
				  area_weighted_global_mean, get_deg_lat, get_deg_lon, get_grid_boundaries, get_cos_lat
 use       mpp_domains_mod, only: mpp_domains_init, domain2D,        &
                                  mpp_define_domains, mpp_get_compute_domain, mpp_global_field
 use  spectral_damping_mod, only: spectral_damping_init,   compute_spectral_damping
 use      horiz_interp_mod, only: horiz_interp
 use          leapfrog_mod, only: leapfrog
 use        rand_number_mod, only: gasdev
 
 implicit none
 
 character(len=128) :: exp_name    = "temp_adv"
 integer            :: num_time    = 6000                ! days
  
 real               :: theta_0=300., delta_theta=60., tau=10.
 real               :: phi_0=45., delta_phi=10, wavenumber=6., angular_speed=15., zeta_0=1.e-9, tau_s=2.
 logical            :: velocity_const = .false., total_gradient = .true.
 
 character(*), parameter :: SOURCE_ROOT  = "./"               ! work directory set by the run script
 character(*), parameter :: ARCHIVE_ROOT =  SOURCE_ROOT
 character(len=128) :: path, path2, path3

 !------------------ standard settings for T42 -----------------------------
 integer :: num_lon            = 128
 integer :: num_lat            = 64
 integer :: num_fourier        = 42
 integer :: num_spherical      = 43
 integer :: fourier_inc        = 1

 !------------------ standard settings for T85 -----------------------------
! integer :: num_lon          = 256
! integer :: num_lat          = 128
! integer :: num_fourier      = 85
! integer :: num_spherical    = 86
! integer :: fourier_inc      = 1
 
 integer            :: damping_order      = 4
 real               :: damping_coeff      = 1.15740741e-4     ! (one tenth day)**-1
 character(len=64)  :: damping_option     = 'resolution_dependent'

 real               :: robert_coeff = .04, dt_real = 1200.
 integer            :: num_step
 integer            :: previous, current, future, t
 
  !--------------------------------------------------------------------------
 logical :: root_pe
 integer :: npes, pe
 integer :: i, j, k, l, is, ie, js, je, n, ms, me, ns, ne, start, total_time, list, idum, unit, ierr, io
 integer :: iduminit  = -855287
 real    :: r1, r2, ks, delta_t
 
 type(NCFILE_TYPE) :: File

!--------------------------------------------------------------------------
 real,    allocatable, dimension(:,:):: ug, vg, theta_eq, forcing
 real,    allocatable, dimension(:)  :: lon, lat, lonb, latb, beta
 real,    allocatable, dimension(:,:,:) :: q_tracer, vorg
 complex, allocatable, dimension(:,:,:) :: qs_tracer, vors

 namelist /Tpdf_nml/ exp_name, num_time, theta_0, delta_theta, tau, phi_0, delta_phi,          &
                     wavenumber, angular_speed, zeta_0, tau_s, velocity_const, total_gradient, &
                     num_lon, num_lat, num_fourier, num_spherical, fourier_inc,                &
		     damping_order, damping_coeff, damping_option, dt_real

!--------------------------------------------------------------------------
    unit = open_namelist_file()
    ierr=1
    do while (ierr /= 0)
      read(unit, nml=Tpdf_nml, iostat=io, end=20)
      ierr = check_nml_error (io, 'Tpdf_nml')
    enddo
20  call close_file (unit)
write (stdlog(), nml=Tpdf_nml)

!--------------------------------------------------------------------------
 call fms_init ();  call fms_io_init;  call constants_init
 pe = mpp_pe();     npes = mpp_npes(); root_pe = (pe == mpp_root_pe())
 call transforms_init(radius, num_lat, num_lon, num_fourier, fourier_inc, num_spherical, &
         south_to_north=.true., triang_trunc=.true., longitude_origin=0.0 )	
 call get_grid_domain(is, ie, js, je)
 call get_spec_domain(ms, me, ns, ne)

 call allocate_fields

 call get_deg_lon(lon)
 call get_deg_lat(lat)
 call get_grid_boundaries(lonb, latb, global=.true.)
 
 beta = 2*omega*cos(lat*pi/180.)/radius

 num_step     = int(86400./dt_real)  ! 1 day/dt_real

!initialization of q tracers and damping
 call spectral_damping_init(damping_coeff, damping_order, damping_option, num_fourier, num_spherical, 1, 0., 0., 0.)

 previous = 1; current  = 1
 
 do i=1, num_lon
 do j=1, num_lat
   theta_eq(i,j)= theta_0 - delta_theta*sin(lat(j)*pi/180.)*sin(lat(j)*pi/180.)
 end do
 end do
 
 do j=1, size(lat)
   q_tracer(:,j,previous) = theta_0
 end do
 call trans_grid_to_spherical(q_tracer(:,:,previous), qs_tracer(:,:,previous))
 q_tracer(:,:,current) =q_tracer(:,:,previous)
 qs_tracer(:,:,current)=qs_tracer(:,:,previous)
 
 vorg = 0.; vors = 0.
 
 idum = iduminit

 if(tau_s /= 0.0) ks = 1.0/abs(tau_s)/86400.

!--------------------------------------------------------------------------
if(root_pe) print*, "--------------------------------------------------------------------------"
    path = trim(ARCHIVE_ROOT)//trim(exp_name)//".nc"
    call open_ncfile(path, File)
    call write_variable (File, 'theta_eq', (/1,1/), theta_eq )

    !ccccccccccccccccccccccccccccccc
    do l=1, num_time
      do t=1, num_step
        idum = idum - mpp_pe()*37
	if(velocity_const) then
	   r1 = 1.0; r2 = 0.0
	   else
  	   r1 = gasdev(idum)
	   r2 = gasdev(idum)
!	   print*, r1, r2
	endif

        if(previous == current) then
           delta_t = dt_real 
           future = 2
        else
           delta_t = 2*dt_real
           future = previous
        endif

        ! compute the linear forced-dissipiated vorticity equation
        call generate_stirring(vorg, vors, ug, vg, forcing)
	
        ! advection-diffusion scheme 
        call compute_tracer_advection(ug, vg, q_tracer, qs_tracer)
	
        previous = current
        current  = future
      end do
      
      call write_variable (File, 'ug',     (/1,1,l/), ug )
      call write_variable (File, 'vg',     (/1,1,l/), vg )
      call write_variable (File, 'forcing',(/1,1,l/), forcing )
      call write_variable (File, 'vorg',   (/1,1,l/), vorg(:,:,current) )
      call write_variable (File, 'theta',  (/1,1,l/), q_tracer(:,:,current) )
      call write_time_coord(File, l, real(l))
    end do   ! end of num_time
    call close_ncfile(File)
    !ccccccccccccccccccccccccccccccc
if(root_pe) print*, "--------------------------------------------------------------------------"

!--------------------------------------------------------------------------
 call transforms_end
 call fms_io_exit
 call fms_end

!--------------------------------------------------------------------------
contains

subroutine allocate_fields
allocate (lon(num_lon))
allocate (lat(num_lat))
allocate (lonb(num_lon+1))
allocate (latb(num_lat+1))
allocate (beta(num_lat))

allocate (ug(is:ie, js:je))
allocate (vg(is:ie, js:je))
allocate (theta_eq(is:ie, js:je))
allocate (forcing(is:ie, js:je))

allocate (q_tracer(is:ie, js:je,2))
allocate (qs_tracer(ms:me, ns:ne, 2))
allocate (vorg(is:ie, js:je,2))
allocate (vors(ms:me, ns:ne, 2))

end subroutine allocate_fields
  
!-------------------------------------------------------------------------
subroutine open_ncfile(filename, File)

    character(*), intent(in) :: filename
    type(NCFILE_TYPE), intent(inout) :: File
    
      call define_file(File, filename)
! axes
      call define_axis(File, 'lon',  lon,  'degrees_E', 'x', 'longitude');
      call define_axis(File, 'lonb', lonb, 'degrees_E', 'x', 'longitude edges');
      call define_axis(File, 'lat',  lat,  'degrees_N', 'y', 'latitude');
      call define_axis(File, 'latb', latb, 'degrees_N', 'y', 'latitude edges');
      call define_time_axis(File, 'time');
! variables
      call define_variable(File, 'ug',(/'lon','lat','time'/),'','ug')
      call define_variable(File, 'vg',(/'lon','lat','time'/),'','vg')
      call define_variable(File, 'forcing',(/'lon','lat','time'/),'','forcing')
      call define_variable(File, 'vorg',(/'lon','lat','time'/),'','vorg')
      call define_variable(File, 'theta',(/'lon','lat','time'/), '','theta')
      call define_variable(File, 'theta_eq',(/'lon','lat'/),'','theta_eq')

      call write_file_header(File)

end subroutine open_ncfile

!-------------------------------------------------------------------------
subroutine compute_tracer_advection(u, v, q, qs)
 real, intent(inout), dimension(is:ie, js:je, 2) :: q
 complex, intent(inout), dimension(ms:me, ns:ne, 2) :: qs
 real, intent(in), dimension(is:ie, js:je) :: u, v
 
 real, dimension(is:ie, js:je) :: dt_q
 complex, dimension(ms:me, ns:ne) :: dt_qs, qs_zm
 
 do j=js, je
   dt_q(:,j) = -(q(:,j,current)-theta_eq(:,j))/(86400.*tau)
 end do
 
 if(total_gradient) then
   call horizontal_advection (qs(:,:,current), u, v, dt_q)
 else
   ! use zonal mean tracer gradient
   if (ms==0) then
     qs_zm(ms,ns:ne)      = qs(ms,ns:ne,current)
     qs_zm(ms+1:me,ns:ne) = 0
   else
     qs_zm(ms:me,ns:ne) = 0
   endif
   call horizontal_advection (qs_zm, u, v, dt_q)
 endif
  
 call trans_grid_to_spherical (dt_q, dt_qs)

 call compute_spectral_damping(qs(:,:,previous), dt_qs, delta_t)
 
 call leapfrog(qs, dt_qs, previous, current, future, delta_t, robert_coeff)
 
 call trans_spherical_to_grid(qs(:,:,future), q(:,:,future))
  
end subroutine compute_tracer_advection

!===================================================================================
 
subroutine generate_stirring(vorg, vors, u, v, forcing)
 real, intent(inout), dimension(is:ie, js:je, 2) :: vorg
 complex, intent(inout), dimension(ms:me, ns:ne, 2) :: vors
 real, intent(out), dimension(is:ie, js:je) :: u, v, forcing
 
 real, dimension(is:ie, js:je) :: dvorg_dx
 real, dimension(is:ie, js:je) :: dt_vorg
 complex, dimension(ms:me, ns:ne) :: dt_vors
  
 dvorg_dx = compute_lon_deriv(vorg(:,:,current))
 call uv_grid_from_vor_div(vors(:,:,current), vors(:,:,current)*0.0, u, v)
 do i=is, ie
 do j=js, je
    forcing(i,j) = zeta_0*exp(-((abs(lat(j))-phi_0)/delta_phi)**2)  &
                         *(r1*cos(wavenumber*lon(i)*pi/180.) -r2*sin(wavenumber*lon(i)*pi/180.))
    dt_vorg(i,j) = -angular_speed*cos(lat(j)*pi/180.)*dvorg_dx(i,j) -v(i,j)*beta(j) -ks*vorg(i,j,previous) +forcing(i,j)
 end do
 end do
 
 call trans_grid_to_spherical(dt_vorg, dt_vors)

 call leapfrog(vors, dt_vors, previous, current, future, delta_t, robert_coeff)

 call trans_spherical_to_grid(vors(:,:,future), vorg(:,:,future))

return

end subroutine generate_stirring

!===================================================================================

function compute_lon_deriv(grid_1) result(grid_2)

 real,    dimension(is:ie, js:je) :: grid_1, grid_2
 complex, dimension(ms:me, ns:ne) :: spec_1, spec_2

 call trans_grid_to_spherical(grid_1, spec_1, do_truncation = .true.)
 spec_2 = compute_lon_deriv_cos(spec_1)
 call trans_spherical_to_grid(spec_2, grid_2)
 call divide_by_cos(grid_2)
 
end function compute_lon_deriv

end program main
