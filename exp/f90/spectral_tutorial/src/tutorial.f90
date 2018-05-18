program tutorial

! A TUTORIAL FOR TRANSFORMS_MOD, A MODULE FOR THE TRANSFORMATION OF 
! GRID FIELDS DEFINED ON THE SPHERE TO SPHERICAL HARMONICS, AND 
! MANIPULATIONS WITHIN THE SPHERICAL HARMONIC REPRESENTATION

!  RUN THE PROGRAM AND STUDY THE OUTPUT AND THE PROGRAM SIMULTANEOUSLY.
!  ALSO, COMPARE YOUR OUTPUT TO THE ONE THAT WE GENERATE WITH THIS
!  PROGRAM (TUTORIAL_OUTPUT WHICH WAS CREATED RUNNING ON 8 PROCESSORS)



! MODULES USED BY THIS TUTORIAL PROGRAM

use transforms_mod , only : transforms_init,          &
                            transforms_end,           &
                            get_lat_max,              &
                            get_sin_lat,              &
                            get_cos_lat,              &
                            get_wts_lat,              &
                            get_deg_lat,              &
                            get_deg_lon,              &
                            get_grid_boundaries,      &
                            get_fourier_wave,         &
                            get_spherical_wave,       &
                            get_eigen_laplacian,      &
                            rhomboidal_truncation,    &
                            triangular_truncation,    &
                            compute_legendre,         &
                            trans_grid_to_spherical,  &
                            trans_spherical_to_grid,  &
                            compute_lat_deriv_cos,    &
                            divide_by_cos,            &
                            divide_by_cos2,           &
                            compute_gradient_cos,     &
                            compute_laplacian,        &
                            compute_vor,              &
                            compute_div,              &
                            compute_vor_div,          &
                            get_grid_domain,          &
                            get_spec_domain,          &
                            uv_grid_from_vor_div,     &
                            vor_div_from_uv_grid,     &
                            compute_ucos_vcos,        &
                            horizontal_advection,     &
                            grid_domain,              &
                            spectral_domain

use mpp_mod, only : mpp_sync

use mpp_domains_mod, only : mpp_global_sum, BITWISE_EXACT_SUM

use fms_mod, only : write_version_number, mpp_pe, mpp_root_pe, mpp_npes, &
                    fms_init, fms_end

use constants_mod, only : constants_init

!======================================================================
implicit none

!======================================================================

! INPUT NEEDED FOR INITIALIZATION OF TRANSFORMS_MOD

character(len=128), parameter :: version = '$Id tutorial.f90 $'
character(len=128), parameter :: tagname = '$Name:  $'

real  :: radius            = 1.0 
real  :: longitude_origin  = 0.0 

integer, parameter :: num_lat       = 64
integer, parameter :: num_lon       = 128 
integer, parameter :: fourier_inc   =  1 
integer, parameter :: num_fourier   = 42 
integer, parameter :: num_spherical = 43

logical :: south_to_north  = .true.
logical :: triang_trunc    = .true.

! radius:  THE RADIUS OF THE SPHERE.  UNITS ARE UP TO THE USER. HERE WE
!   USE A UNIT SPHERE FOR SIMPLICITY.  IF CONSTRUCTING AN ATMOSPHERIC MODEL,
!   ONE CAN GET THE RADIUS OF THE EARTH FROM constants_mod.
  
! longitude origin:  GRID FIELDS ARE DEFINED ON A LON-LAT GRID; 
!   THE LONGITUDES ARE EQUALLY SPACED, THIS IS THE LOCATION 
!   OF THE FIRST LONGITUDINAL GRID POINT, IN RADIANS.
  
! fourier_inc: THE ONLY FOURIER WAVES RETAINED ARE MULTIPLES OF fourier_inc.  
!   FOR A STANDARD MODEL, SET fourier_inc = 1; 
!   fourier_inc = 3, FOR EXAMPLE, WOULD BE USED IF ONE WAS CONSTRUCTING A 
!   MODEL IN WHICH ALL FIELDS WERE ASSUMED TO BE THREE-FOLD
!   SYMMETRIC AROUND A LATITUDE CIRCLE
!   (THE OPTION OF SETTING fourier_inc TO SOMETHING OTHER THAN 1 
!    HAS NOT BEEN THOROUGHLY CHECKED RECENTLY, AND SHOULD BE USED WITH
!    CAUTION)

! num_fourier:  THE MAXIMUM FOURIER WAVENUMBER IN THE SPHERICAL HARMONIC 
!   BASIS.  THE RETAINED FOURIER MODES ARE (0, 1N, 2N, 3N, ..., num_fourier)
!   WHERE N = fourier_inc.  (num_fourier MUST BE DIVISIBLE BY fourier_inc)
!   PERFORMANCE IS BEST IF NUM_FOURIER IS A MULTIPLE OF FACTORS OF 
!   2,3 AND 5.  

! num_spherical: THE MAXIMUM MERIDIONAL WAVENUMBER IN THE SPHERICAL HARMONIC
!   BASIS.  FOR STANDARD "TRIANGULAR" AND "RHOMBOIDAL" TRUNCATIONS (SEE BELOW)
!   ALWAYS SET num_spherical = num_fourier + 1.  

!   (THE VALUES FOR num_lat, num_lon, num_fourier, num_spherical ABOVE ARE
!    THOSE FOR A STANDARD T21 TRUNCATION -- THE NUMBER OF GRID POINTS IS 
!    LARGE ENOUGH TO AVOID ALIASING OF QUADRATIC PRODUCTS, AND THE NUMBER 
!    OF LONGITUDINAL POINTS IS A NICE POWER OF TWO.  

! south_to_north:  THE LATITUDES RUN FROM SOUTH TO NORTH IF THIS IS
!   TRUE AND FROM NORTH TO SOUTH IF FALSE

! triang_trunc:  SOME ROUTINES AUTOMATICALLY TRUNCATE SPECTRAL FIELDS
!   BY ZERO'ING COEFFICIENTS OF HARMONICS LYING OUTSIDE OF THE TRUNCATION.
!   IF triang_trunc = .true. THIS TRUNCATION IS PERFORMED USING A TRIANGULAR
!   TRUNCATION (SEE BELOW)

!======================================================================

logical :: check_checking = .false.

!   TRANSFORMS_MOD HAS A NAMELIST THAT ALLOWS ONE TO TURN ON
!   AUTOMATIC CHECKING, TO SEE IF SPECTRAL FIELDS, WHICH ARE PRESUMED 
!   TO BE TRANSFORMS OF REAL DATA, HAVE ZERO IMAGINARY PARTS IN THE 
!   ZONALLY SYMMETRIC COMPONENTS OF THE FIELD
!   IF THIS IS NOT THE CASE, IT GENERALLY MEANS THAT THERE IS A SERIOUS
!   BUG IN ONE'S CODE.  AS THIS CHECKING CAN TAKE A NON-TRIVIAL AMOUNT OF TIME
!   IN SOME CONTEXTS, ONE GENERALLY TURNS IT OFF ONCE A PROGRAM IS 
!   DEBUGGED.  TO TURN IT ON, SET 
!   check = .true. IN transforms_nml .  THIS NAMELIST SHOULD BE COPIED OR APPENDED TO 
!   A FILE NAMED input.nml THAT IS PLACED IN THE SAME DIRECTORY AS THE 
!   EXECUTABLE FOR THIS PROGRAM.  IF YOU SET check_checking = .true.
!   THIS PROGRAM WILL TRY TO TRANSFORM AN UNACCEPTABLE FIELD AND WILL
!   BOMB IF YOU HAVE SET THE NAMELIST OPTION APPROPRIATELY, (AND IF 
!   YOU ARE READING IN THE NAMELIST PROPERLY!)


! SPACE NEEDED BY THIS PROGRAM -- THESE FIELDS WILL BE DESCRIBED
!   AS THEY ARE USED

real,    allocatable, dimension(:)     :: sin_lat, cos_lat, wts_lat
real,    allocatable, dimension(:)     :: sin_lat_global, wts_lat_global
real,    allocatable, dimension(:)     :: wts_lat_check
real,    allocatable, dimension(:)     :: lat_boundaries, lon_boundaries
real,    allocatable, dimension(:)     :: lat_boundaries_global
real,    allocatable, dimension(:)     :: deg_lat, deg_lon
real,    allocatable, dimension(:)     :: deg_lat_global

integer, allocatable, dimension(:,:)   :: mm, nn, ll
integer, allocatable, dimension(:,:)   :: mm_global, nn_global, ll_global
real   , allocatable, dimension(:,:)   :: eigen, eigen_global

complex, allocatable, dimension(:,:)   :: spec_1, spec_2, spec_3, spec_global
complex, allocatable, dimension(:,:)   :: u_cos, v_cos, u_cosm, v_cosm
complex, allocatable, dimension(:,:)   :: vor, div, vor1, div1, stream

real,    allocatable, dimension(:,:)   :: grid_1, grid_2, grid_3, grid_4, streamg, tendg
real,    allocatable, dimension(:,:)   :: u_cosg, v_cosg, vorg, divg, ug, vg, vorg1, divg1

complex, allocatable, dimension(:,:,:) :: nspec
real,    allocatable, dimension(:,:,:) :: ngrid

real,    dimension(0:num_fourier,   0:num_spherical, num_lat) :: legendre
real,    dimension(0:num_spherical, 0:num_spherical)          :: norm

complex, allocatable, dimension(:,:) :: work_s
real,    allocatable, dimension(:,:) :: work_g

!======================================================================

! THESE ARE IMPORTANT (SEE BELOW)

integer :: is, ie, js, je, ms, me, ns, ne
!======================================================================

! THESE ARE NOT IMPORTANT

integer :: pr, i, j, k, m, n, n1, p, pe, npes, global, power
real    :: sum, x, sum_pars_g, sum_pars_s, error, e1, e2, e3, e4, e5, e6
logical :: root_pe
!======================================================================

call fms_init()
call constants_init

call write_version_number(version, tagname)


! THE FIRST STEP IS TO INITIALIZE THE MODULE


call transforms_init(                         &
     radius, num_lat, num_lon,                &
     num_fourier, fourier_inc, num_spherical, &
     south_to_north    = south_to_north,      &
     triang_trunc      = triang_trunc,        &
     longitude_origin  = longitude_origin     )

!  transforms_mod INIITIALIZES mpp_mod AUTOMATICALLY, SO YOU 
!   DO NOT HAVE TO DO IT YOURSELF

!  THE DEFAULTS ON THE OPTIONAL ARGUMENTS ARE

!   south_to_north   = .true.
!   triang_trunc     = .true.
!   longitude_origin =  0.0

! IF THESE ARE OK, YOU COULD USE THE SIMPLER CALL
!  call initialize_transforms(radius, num_lat, num_lon, &
!	     num_fourier, fourier_inc, num_spherical)




! THE INITIALIZATION CALL WILL PRODUCE SOME OUTPUT DESCRIBING THE DOMAIN
!   DECOMPOSITION.  LET'S IGNORE IT.


! FROM mpp_mod, WE CAN OBTAIN THE NUMBER OF OUR PROCESSOR 	 
! AND THE TOTAL NUMBER OF PROCESSORS

pe   = mpp_pe()
npes = mpp_npes()
root_pe = (pe == mpp_root_pe())
call mpp_sync() ! THIS IS ALSO AVAILABLE FROM MPP_MOD, AND IS USED IN THIS 
                ! TUTORIAL PROGRAM ONLY TO INSURE THAT THE PRINTOUT IS LEGIBLE


!Output #1
!############## printing #############################################  
10000 format(1x,/,1x,'================', i3, ' ====================',/)
pr = 1    
call mpp_sync(); if(root_pe) write(*,10000) pr 
if(root_pe) then
  print*, 'total # of processors =', npes
endif
!############## printing #############################################


! WE ALSO HAVE AVAILABLE THE LIMITS OF THE LOCAL
! GRID AND SPECTRAL DOMAINS ON EACH PROCESSOR: 
!   is, ie, js, je, ms, me, ns,ne,
!   WHERE GRID FIELDS ARE DIMENSIONED (is:ie,js:je) 
! AND SPECTRAL FIELDS ARE DIMENSIONED (ms:me,ns:ne)

call get_grid_domain(is,ie,js,je)
call get_spec_domain(ms,me,ns,ne)

! THESE LIMITS ARE OFTEN USED TO DIMENSION ARRAYS WHICH HOLD DATA 
! WHICH IS DISTRIBUTED ACROSS PROCESSORS. FOR EXAMPLE, WE USE THEM
! HERE TO DIMENSION TWO WORK ARRAYS WHICH WILL BE USED LATER.

allocate ( work_g (is:ie, js:je))
allocate ( work_s (ms:me, ns:ne))

! DON'T WORRY ABOUT THE FORM OF THE PRINT STATEMENTS IN THIS PROGRAM -- 
!  IT IS DIFFICULT TO GET THE PROCESSORS TO PRINT OUT IN SUCCESSION TO 
!  GET NICE LOOKING OUTPUT -- 

!Output #2
!############## printing #############################################
pr = 2
call mpp_sync(); if(root_pe) write(*,10000) pr 
do p = 0, npes-1
  call mpp_sync()
  if(pe == p) write(*,1111) pe, is, ie, js, je, ms, me, ns, ne
enddo
1111 format(1x,'pe =', i4, 2x, 'grid(',i3,':',i3,',',i3,':',i3,')', &
                           2x, 'spectral(',i3,':',i3,',',i3,':',i3,')')

call mpp_sync()
!############## printing #############################################


! CHECKING YOUR OUTPUT, NOTE THAT THE DOMAIN DECOMPOSITION IN THE GRID 
! DOMAIN IS IN LATITUDE, WITH ALL POINTS ON A GIVEN LATITUDE CIRCLE ON 
! THE SAME PROCESSOR. THIS 1D DOMAIN DECOMPOSITION WILL BE GENERALIZED IN 
! LATER RELEASES

! IN THE SPECTRAL DOMAIN, THE DECOMPOSITION IS IN ZONAL WAVENUMBER,
! WITH ALL MERIDIONAL WAVES WITH THE SAME VALUE OF ZONAL WAVENUMBER 
! RESIDING ON THE SAME PROCESSOR



! THERE ARE A VARIETY OF SIMPLE CALLS THAT RETRIEVE VALUES OF
!  USEFUL PARAMETERS 

! FOR EXAMPLE, TO RETRIEVE THE NUMBER OF LATITUDES OVER THE ENTIRE
!  GLOBE, (FOR HISTORICAL REASONS, THE INTERFACE REFERS TO THIS NUMBER AS
! lat_max RATHER THAN num_lat)

call get_lat_max(n)


!Output #3
!############## printing #############################################
pr = 3
call mpp_sync(); if(root_pe) write(6,10000) pr
call mpp_sync(); if(root_pe) print*, 'num_lat =', n
!############## printing #############################################


! OTHER CALLS OF THIS TYPE ARE 

! get_num_fourier
! get_fourier_inc
! get_num_spherical
! get_triang_trunc
! get_lon_max
! get_longitude_origin

! MORE INTERESTINGLY, ONE CAN RETRIEVE 

! deg_lat = THE GAUSSIAN LATITUDES AT WHICH FIELDS ARE DEFINED (IN DEGREES)
! deg_lon 
! sin_lat
! cos_lat
! cosm_lat  = 1/cos
! cosm2_lat = 1/(cos*cos)
! wts_lat   = THE GAUSSIAN WEIGHTS (MORE ON THESE BELOW)
! lat_boundaries = LATITUDES OF BOUNDARIES OF GRID CELLS (IN RADIANS)
! lon_boundaries = LONGITUDES "" ""

! WE CAN RETRIEVE THE GLOBAL FIELDS OF THESE QUANTITIES OR JUST THE VALUES 
! FOR THE LATITUDES ON PROCESSOR

!FIRST ALLOCATE SOME FIELDS

allocate(sin_lat        (js:je     ))
allocate(sin_lat_global ( 1:num_lat))

allocate(wts_lat        (js:je     ))
allocate(wts_lat_global ( 1:num_lat))
allocate(wts_lat_check  ( 1:num_lat))

allocate(lat_boundaries        (js:je+1     ))
allocate(lat_boundaries_global ( 1:num_lat+1))
allocate(lon_boundaries        ( 1:num_lon+1))

allocate(deg_lat        (js:je    ))
allocate(deg_lat_global (1:num_lat))
allocate(deg_lon        (1:num_lon))


! USING THE MAGIC OF FORTRAN 90 WE CAN USE THE SAME CALL TO RETRIEVE 
!  THE GLOBAL OR THE LOCAL FIELDS 

call get_sin_lat(sin_lat)
call get_sin_lat(sin_lat_global)

call get_wts_lat(wts_lat)
call get_wts_lat(wts_lat_global)

call get_deg_lat(deg_lat)
call get_deg_lat(deg_lat_global)
call get_deg_lon(deg_lon)

call get_grid_boundaries(lon_boundaries, lat_boundaries)
call get_grid_boundaries(lon_boundaries, lat_boundaries_global, global= .true.)

! YOU MAY ASK WHY, IN THE LAST CASE, WE NEED TO SPECIFY THE OPTIONAL
! ARGUMENT global, RATHER THAN USING THE SAME "OVERLOADING" TRICK USED
! TO PROGRAM THE OTHERS -- BEATS ME!
! NOTE THAT THE DEFAULT VALUE OF GLOBAL IS .FALSE. 


!Output #4
!############## printing #############################################
pr = 4
call mpp_sync()
if(root_pe) then
  write(*,10000) pr
  write(*,1000)
  do j =1, num_lat
    write(*,1001) j, deg_lat_global(j), sin_lat_global(j), wts_lat_global(j)
  enddo
endif
1000 format(1x,/,1x, 'global latitudes',/,8x,'deg_lat',6x,'sin_lat', &
            6x, 'wts_lat',/)
1001 format(1x,i4,1x,f9.4,2(1x,e13.6))

call mpp_sync(); if(root_pe) write(*,1002)
do p = 0, 0
  call mpp_sync()
  if(pe == p) then
    do j = js, je
      write(*,1003) j, deg_lat(j), sin_lat(j), wts_lat(j)
    enddo
  endif
enddo
1002 format(1x,//,1x, 'local latitudes on processor 0',/,8x,'deg_lat',&
       6x,'sin_lat', 6x, 'wts_lat',/)
1003 format(1x,i4,1x,f9.4,2(1x,e13.6))
!############## printing #############################################




! THE GAUSSIAN WEIGHTS ARE USED IN ALL INTEGRATIONS OVER LATITUDE.
! (INTEGRAL OF F COS(LAT) DLAT = SUM OVER J OF F(J)*WTS_LAT(J) )
! THINK OF WTS_LAT(J) AS COS(LAT)*DEL_LAT WHERE DEL_LAT IS THE
! INCREMENT IN LATITUDE REPRESENTED BY THAT POINT, OR, MORE PRECISELY,
! WTS_LAT(J) = SIN(LAT_BOUNDARIES(J+1)) - SIN(LAT_BOUNDARIES(J)),
! WHERE LAT_BOUNDARIES(J+1) AND LAT_BOUNDARIES(J) ARE THE LATITUDES
! OF THE BOUNDARIES OF THE CELL J

! THUS, THE SUM OVER J OF THE GAUSSIAN WEIGHTS SHOULD EQUAL 2.0.
! LET'S CHECK THIS

sum = 0.0
do j = 1, num_lat
  sum = sum + wts_lat_global(j)
enddo


!Output #5
!############## printing #############################################
pr = 5
call mpp_sync(); if(root_pe) write(*,10000) pr
call mpp_sync(); if(root_pe) write(*,1004) sum, sum - 2.0
1004 format(1x,'sum of Gaussian wts = ', f10.6, 3x, 'sum - 2 = ',e13.4)
!############## printing #############################################


! WHILE WE'RE AT IT WE MIGHT AS WELL CHECK THE RELATION BETWEEN
! LAT_BOUNDARIES AND WTS_LAT


do j = 1, num_lat
  wts_lat_check(j) =  sin(lat_boundaries_global(j+1))  &
                    - sin(lat_boundaries_global(j))
enddo


!Output #6
!############## printing #############################################
pr = 6
call mpp_sync(); if(root_pe) write(*,10000) pr
if(root_pe) then
  write(*,1005)
  do j = 1,num_lat
    write(*,1006) j, wts_lat_global(j), wts_lat_check(j)
  enddo
endif
1005 format(1x, /,1x, 'checking consistency of weights and grid boundaries',/)
1006 format(1x,i4,2(1x,e13.6))
!############## printing #############################################



! WE MIGHT AS WELL PRINT OUT THE LONGITUDES AND LONGITUDE BOUNDARIES
! NOTE THAT LON_BOUNDARIES IS IN RADIANS AND DEG_LON IN DEGREES.
! THE FIRST GRID POINT IS LOCATED AT  longitude_origin , WHICH IS 
! PROVIDED TO THE INITIALIZATION ROUTINE IN RADIANS.
! YOU MIGHT WANT TO TRY TO REDEFINE  longitude_origin  SO THAT
! THE FIRST CELL HAS A WESTERN SIDE LOCATED AT 0.0 

! ALSO, IF YOU PLAY AROUND WITH  longitude_origin,  YOU WILL DISCOVER THAT
! GET_GRID_BOUNDARIES ALWAYS RETURNS A MONOTONIC SERIES OF LONGITUDES,
! SOME OF WHICH MAY BE NEGATIVE, 
! WHILE GET_DEG_LON RETURNS LONGITUDES BETWEEN 0.0 AND 360.0


!Output #7
!############## printing #############################################
pr = 7
call mpp_sync(); if(root_pe) write(*,10000)
if(root_pe) then
  write(*,1007)
  do i = 1, num_lon
    write(*,1008) i, lon_boundaries(i), deg_lon(i)
  end do
  write(*,1008) num_lon+1, lon_boundaries(num_lon+1)
endif
1007 format(1x, /,4x, 'lon bounds of cells(radians)   lon of grid points(deg)', /)
1008 format(1x,i4,2(8x,e13.6))
!############## printing #############################################



! NEXT WE EXAMINE HOW SPECTRAL FIELDS ARE LAID OUT
! WE FIRST GET THE ZONAL WAVENUMBER, M, THE TOTAL WAVENUMBER L = M + N, 
! AND (-1)*EIGENVALUE OF THE LAPLACIAN ON THE SPHERE = L*(L+1)/(A*A)
!  (WHERE A = radius)

! N IS THE NUMBER OF NODES OF THE MODE IN THE MERIDIONAL DIRECTION
! IF N IS EVEN, THE EIGENFUNCTION IS EVEN ABOUT THE EQUATOR
! L = N + M IS THE 'TOTAL HORIZONTAL WAVENUMBER' IN THE SENSE THAT
! IT DETERMINES THE EIGENVALUE OF THE LAPLACIAN
!  (CAUTION: ONE OFTEN SEES L HERE DENOTED BY THE SYMBOL N)

allocate(mm    (ms:me,ns:ne))
allocate(ll    (ms:me,ns:ne))
allocate(eigen (ms:me,ns:ne))
allocate(mm_global    (0:num_fourier,0:num_spherical))
allocate(ll_global    (0:num_fourier,0:num_spherical))
allocate(eigen_global (0:num_fourier,0:num_spherical))

call get_fourier_wave     (mm)
call get_spherical_wave   (ll)
call get_eigen_laplacian  (eigen) 

call get_fourier_wave     (mm_global)
call get_spherical_wave   (ll_global)
call get_eigen_laplacian  (eigen_global) 



!Output #7
!############## printing #############################################
pr = 7
call mpp_sync(); if( root_pe) write(*,10000) pr
if (root_pe) then
  write(*,1009) (m, m = 0, num_fourier)
  do n = 0,num_spherical
    write(*,1010) n, mm_global(:,n)
  enddo  
endif  
1009 format(1x,/,'zonal wavenumbers',/,3x,'n',1x, 'm =',22(1x,i3),/)
1010 format(1x, i3, 4x, 22(1x,i3))

if (root_pe) then
  write(*,1011) (m, m = 0, num_fourier)
  do n = 0,num_spherical
    write(*,1010) n, ll_global(:,n)
  enddo  
endif  
1011 format(1x,/,'total wavenumbers',/,3x,'n',1x, 'm =',22(1x,i3),/)

if (root_pe) then
  write(*,1012) (m, m = 0, num_fourier)
  do n = 0,num_spherical
    write(*,1013) n, eigen_global(:,n)
  enddo  
endif  
1012 format(1x,/,'negative of eigenvalues of Laplacian',/,3x,'n',1x, &
            'm =',22(3x,i3),/)
1013 format(1x, i3, 5x, 22(1x,f5.0))

if(root_pe) then
  write(*,1014)(m, m = ms, me)
 do n = ns,ne
    write(*,1013) n, eigen(:,n)
  enddo  
endif  
1014 format(1x,/,'eigenvalues on pe = 0',/,3x,'n',1x, &
            'm =',22(3x,i3),/)
!############## printing #############################################



! NOW LET'S DEFINE A GLOBAL SPECTRAL FIELD 
!   ( EVERTHING BELOW SHOULD ALSO WORK FOR LOCAL FIELDS, BUT IT
!     IS THEN MORE TROUBLE PRINTING OUT PROPERLY)

allocate(spec_global(0:num_fourier,0:num_spherical))

spec_global = ll_global


! LET'S SEE WHAT HAPPENS IF WE TRY THIS

call rhomboidal_truncation(spec_global)


!Output #8
!############## printing #############################################
pr = 8
call mpp_sync(); if(root_pe) write(*,10000) pr
if(root_pe) then
  write(*,1015) (m, m = 0, num_fourier)
  do n = 0, num_spherical
    write(*,1016)  n, real(spec_global(:,n))
  enddo
endif
1015 format(1x,/,'default rhomboidal truncation',/,/,&
            3x,'n',1x, 'm =', 22(1x,i3),/)
1016 format(1x, i3, 5x, 22(1x,f3.0))
!############## printing #############################################


! THE LAST MERIDIONAL MODE (N = NUM_SPHERICAL) HAS BEEN ZERO'ED, WHICH IS
! NOT VERY EXCITING.  A RHOMBOIDAL TRUNCATION (R21 IN THIS CASE OR R_N MORE
! GENERALLY) CONTAINS WAVES (m,n) WITH m = 0,N AND n = 0,N.  BUT
! THE SECOND DIMENSION OF OUR SPECTRAL FIELDS RUNS FROM 0 TO N+1.
! THIS EXTRA MERIDIONAL HARMONIC IS NEEDED BECAUSE THE EXACT MERIDIONAL
! DERIVATIVE OF A FIELD TRUNCATED AT n = N CONTAINS AN N+1 COMPONENT.
! SO THIS OPERATION COMES IN HANDY

! MORE GENERALLY, FOR EXAMPLE, WE HAVE


call rhomboidal_truncation(spec_global, 8, 4)


!Output #9
!############## printing #############################################
pr = 9
call mpp_sync(); if(root_pe) write(*,10000) pr
if(root_pe) then
  write(*,1017) (m, m = 0, num_fourier)
  do n = 0, num_spherical
    write(*,1016)  n, real(spec_global(:,n))
  enddo
endif
call mpp_sync()
1017 format(1x,/,'rhomboidal truncation at (8,4)',/,3x,'n',1x, 'm =', &
                  22(1x,i3),/)
!############## printing #############################################


! rhomboidal_truncation(x) IS EQUIVALENT TO 
! rhomboidal_truncation(x, num_fourier, num_spherical)
!  (THIS MORE GENERAL FORM IS RARELY USED)



! NOW LET'S TRY TRIANGULAR TRUNCATION

spec_global = ll_global
call triangular_truncation(spec_global)


!Output #10
!############## printing #############################################
pr = 10
call mpp_sync(); if(root_pe) write(*,10000) pr
if(root_pe) then
  write(*,1018) (m, m = 0, num_fourier)
  do n = 0, num_spherical
    write(*,1016)  n, real(spec_global(:,n))
  enddo
endif
1018 format(1x,/,'default triangular truncation',/,3x,'n',1x, 'm =', 22(1x,i3),/)
!############## printing #############################################



call triangular_truncation(spec_global,7)


!Output #11
!############## printing #############################################
pr = 11
call mpp_sync(); if(root_pe) write(*,10000) pr
if(root_pe) then
  write(*,1019) (m, m = 0, num_fourier)
  do n = 0, num_spherical
    write(*,1016)  n, real(spec_global(:,n))
  enddo
endif
1019 format(1x,/,'triangular truncation at 7',/,3x,'n',1x, 'm =', 22(1x,i3),/)
!############## printing #############################################


! triangular truncation(x) IS EQUIVALENT TO 
! triangular_truncation(x, num_fourier)


! TRIANGULAR TRUNCATION IS A SPHERICALLY ISOTROPIC TRUNCATION BECAUSE THE 
! THE ROTATION OF A FIELD ON THE SPHERE MIXES HARMONICS WITH THE SAME
! VALUES OF TOTAL WAVENUMBER L, BUT DOES NOT GENERATE ANY NEW VALUES OF L
! WITH A RHOMBOIDAL TRUNCATION, ONE INCLUDES SOME WAVES WITH PARTICULAR 
! VALUES OF L BUT NOT OTHERS WITH THE SAME VALUE OF L, SO A ROTATION
! OF A RHOMBOIDALLY TRUNCATED FIELD GEERNATED HARMONICS OUTSIDE OF THE 
! TRUNCATION



!  LET'S LOOK AT THE LEGENDRE POLYNOMIALS.
!  THE FOLLOWING CALL PROVIDES THE GLOBAL POLYNOMIALS ON EACH PROCESSOR
!  ONE DOES NOT NEED TO SEE THESE EXPLICITLY TO PERFORM TRANSFORMS,
!  BUT LET'S LOOK AT THEM ANYWAY


call compute_legendre(legendre, num_fourier, fourier_inc, num_spherical, &
        sin_lat_global, num_lat)


! CHECK THE ORTHOGONALITY AND NORMALIZATION BY COMPUTING 

!   SUM {P(M,N,J)*P(M'N',J)*WTS(J) 
!    J

! FIRST, LET'S CHECK M = 0

norm = 0.0
do j = 1, num_lat
  do n = 0, num_spherical
    do n1 = 0, num_spherical
      norm(n,n1) = norm(n,n1) &
           + legendre(0,n,j)*legendre(0,n1,j)*wts_lat_global(j)
    enddo
  enddo
end do

!Output #12
!############## printing #############################################
pr = 12
call mpp_sync(); if(root_pe) write(*,10000) pr
if(root_pe) then
  write(*,2300) (n1, n1=0, 7)
  do n = 0, num_spherical
    write(*,2303) n, (norm(n,n1), n1 = 0,7)
  enddo
  write(*,2301) (n1, n1=8, 15)
  do n = 0, num_spherical
    write(*,2303) n, (norm(n,n1), n1 = 8,15)
  enddo  
  write(*,2302) (n1, n1=16, 22)
  do n = 0, num_spherical
    write(*,2304) n, (norm(n,n1), n1 = 16,22)
  enddo
endif
2300 format(1x,/,3x, 'orthonormality for m = 0',/,1x, 11(8x, i3))
2301 format(1x,/,1x, 11(8x, i3))
2302 format(1x,/,1x, 10(8x, i3))
2303 format(1x, i3, 11(1x,f10.6))
2304 format(1x, i3, 10(1x,f10.6))

do n = 0, num_spherical
  norm(n,n) = norm(n,n) - 1.0
enddo
error = abs(maxval(norm))

call mpp_sync(); if(root_pe) write(*,2305) error
2305 format(1x,/,' maximum orthonormality error for m = 0:', e11.4)
!############## printing #############################################

! MIGHT AS WELL DO A NON-ZERO M WHILE WE'RE AT IT

m = 5
norm = 0.0
do j = 1, num_lat
  do n = 0, num_spherical
    do n1 = 0, num_spherical
      norm(n,n1) = norm(n,n1) &
           + legendre(m,n,j)*legendre(m,n1,j)*wts_lat_global(j)
    enddo
  enddo
end do
do n = 0, num_spherical
  norm(n,n) = norm(n,n) - 1.0
enddo
error = abs(maxval(norm))


!Output #13
!############## printing #############################################
pr = 13
call mpp_sync(); if(root_pe) write(*,10000) pr
call mpp_sync(); if(root_pe) write(*,2306) error
2306 format(1x,/,' maximum orthonormality error for m = 5:', e11.4)
!############## printing #############################################



!  BECAUSE SUM WTS(J) = 2
!            J

!  OUR NORMALIZATION IMPLIES THAT THE (0.,0.) MODE, WHICH IS A CONSTANT,
!  HAS THE VALUE 1/SQRT(2)

!  THIS IS A POTENTIAL SOURCE OF CONFUSION.  FOR EXAMPLE,
!  IF YOU WANT TO SET THE GLOBAL MEAN TEMPERATURE EQUAL TO T, SET THE
!  REAL PART OF THE (0.,0.) MODE OF THE COMPLEX TEMPERATURE FIELD EQUAL TO 
!  T X SQRT(2).

!  THE HARMONIC (0,1) IS PROPORTIONAL TO SIN(LAT). 
!    SINCE THE INTEGRAL OF SIN**2 X COS D(THETA)) = 2/3,
!    THE VALUES OF THESE HARMONICS SHOULD BE SQRT(3/2)*SIN(LAT)

!  LET'S CHECK THIS --


!Output #14
!############## printing #############################################
pr = 14
call mpp_sync(); if(root_pe) write(*,10000) pr
if(root_pe) then
  write(6,1200)
  do j = 1, num_lat
    write(*,1201) j, legendre(0,1,j), sqrt(3./2.)*sin_lat_global(j)
  end do
endif
1200 format(1x,/,1x,'checking normalization of (0,1)',/A)
1201 format(1x, i3, 2(3x,f13.9))
!############## printing #############################################



!  SO, IF YOU WANT TO ADD THE PLANETARY VORTICITY, 2*(OMEGA)*SIN(LAT),
!    TO THE RELATIVE VORTICITY, YOU CAN ADD 2*(OMEGA)*(SQRT(2./3.)) TO THE
!    (0,1) COMPONENT.

!  SINCE THIS NORMALIZATION IS QUITE CONFUSING, ONE SHOULD DEFINE
!  FIELDS IN PHYSICAL SPACE AND THEN TRANSFORM TO THE SPECTRAL DOMAIN
!  WHENEVER POSSIBLE

  

!  NOW DEFINE GRID AND SPECTRAL FIELDS, DISTRIBUTED ACROSS PROCESSORS

allocate(grid_1(is:ie, js:je))
allocate(grid_2(is:ie, js:je))
allocate(grid_3(is:ie, js:je))
allocate(grid_4(is:ie, js:je))
allocate(spec_1(ms:me, ns:ne))
allocate(spec_2(ms:me, ns:ne))
allocate(spec_3(ms:me, ns:ne))

! DEFINE THE GRID FIELD TO EQUAL 1 AT (4,10) AND 0 EVERYWHERE ELSE

grid_1 = 0.0
if(is <= 4 .and. ie >= 4 .and. js <= 4 .and. je >=4) then
  grid_1(4,4) = 1.0
endif


! TRANSFORM TO SPECTRAL DOMAIN AND THEN BACK


call trans_grid_to_spherical(grid_1,spec_1, do_truncation = .true.)
call trans_spherical_to_grid(spec_1,grid_1)

! THE OPTIONAL ARGUMENT IN trans_grid_to_spherical CONTROLS 
! WHETHER OR NOT THE OUTPUT IS TRUNCATED.  THE TRUNCATION IS TRIANGULAR IN
! OUR CASE BECAUSE triang_trunc = .true. IN THE INITALIZATION CALL
! THE DEFAULT IS do_truncation = .false.

! WE PRINT OUT THE FIRST 10 LONGITUDES TO LOOK AT WHAT WE HAVE


!Output #15
!############## printing #############################################
pr = 15
call mpp_sync(); if(root_pe) write(*,10000) pr
call mpp_sync(); if(root_pe) write(*,1205)
do p = 0, npes-1
  call mpp_sync()
  if(pe == p) then
    do j = js,je
      write(*,1206) pe, j, (grid_1(i,j), i = 1,10)
    enddo
  endif
enddo
1205 format(1x, 'set G = 0 except G(4,4) = 1; transform G => S => G1',/)
1206 format(1x, i3, 1x, i3, 10(1x,f10.4))
!############## printing #############################################



! NOW WE HAVE A GRID FIELD THAT IS TRUNCATED.  TRANSFORM BACK AND
! FORTH ONCE AGAIN AND WE SHOULD GET THE SAME FIELD BACK


call trans_grid_to_spherical(grid_1,spec_1, do_truncation = .true.)
call trans_spherical_to_grid(spec_1,grid_2)


!Output #16
!############## printing #############################################
pr = 16
call mpp_sync(); if(root_pe) write(*,10000) pr
call mpp_sync(); if(root_pe) write(*,1207)
do p = 0, npes-1
  call mpp_sync()
  if(pe == p) then
    do j = js,je
      write(*,1206) pe, j, ((grid_2(i,j) - grid_1(i,j))*1.e13, i = 1,10)
    enddo
  endif
enddo
1207 format(1x, 'G1 => S => G2;  printing (G2 - G1)*(10^13)',/)
!############## printing #############################################


! WE SHOULD ALSO CHECK PARSEVAL'S THEOREM.  

! NOTE HOW THE SPECTRAL HARMONICS ARE SUMMED (THE FACTOR OF 2 FOR THE 
! NON-ZONAL HARMONICS ARISES JUST AS IN FOURIER SERIES
! BECAUSE EACH WAVENUMBER M > 0 REPRESENTS BOTH POSITIVE AND 
! NEGATIVE M

! THESE SUMS INVOLVE GLOBAL FIELDS THAT ARE DISTRIBUTED ACROSS ALL 
! PROCESSORS. TO ENSURE THAT GLOBAL SUMS EXACTLY REPRODUCE, REGARDLESS
! OF THE NUMBER OF PROCESSORS USED, THE GLOBAL ARRAY ELEMENTS MUST ALWAYS
! BE SUMMED IN THE SAME SEQUENCE.
! THE FUNCTION mpp_global_sum DOES THIS WHEN WE PASS THE FLAG "BITWISE_EXACT_SUM"
! (BOTH mpp_global_sum AND BITWISE_EXACT_SUM ARE USED FROM mpp_mod)

do j = js,je
  do i = is,ie
    work_g(i,j) = grid_2(i,j)*grid_2(i,j)*wts_lat(j)
  enddo
enddo
sum_pars_g = mpp_global_sum(grid_domain, work_g, BITWISE_EXACT_SUM)/float(num_lon)

do n = ns,ne
  do m = ms,me
    x = 2.0
    if(m == 0) x = 1.0
    work_s(m,n) = x*real(spec_1(m,n)*conjg(spec_1(m,n)))
  enddo
enddo
!!sum_pars_s = mpp_global_sum(spectral_domain, work_s, BITWISE_EXACT_SUM)

!!error = (sum_pars_g - sum_pars_s)/sum_pars_g

!Output #17
!############## printing #############################################
pr = 17
call mpp_sync(); if(root_pe) write(*,10000) pr
call mpp_sync(); if(root_pe) write(*,1210)
print*, 'not supported in this version'
!!call mpp_sync(); if(root_pe) write(*,1211) sum_pars_g, sum_pars_s, error
1210 format(1x,'checking Parseval',/)
1211 format(1x,'grid_variance     = ', e17.10,/,1x, &
               'spectral variance = ', e17.10,/,1x, &
               'relative error    = ', e11.3)
!############## printing #############################################



! NOW COMPUTE THE LATITUDINAL DERIVATIVE OF COS^10
! TRANSFORM TO SPECTRAL DOMAIN, THEN COMPUTE THE DERIVATIVE
!    TIMES COS(LAT), THEN DIVIDE BY COS(LAT)

allocate(cos_lat(js:je))
call get_cos_lat(cos_lat)

power = 10
do j = js,je
  grid_1(:,j) = cos_lat(j)**power
  grid_2(:,j) = -float(power)*sin_lat(j)*(cos_lat(j)**(power-1))
       ! the correct answer
enddo

call trans_grid_to_spherical(grid_1, spec_1)
spec_2 = compute_lat_deriv_cos(spec_1)
call trans_spherical_to_grid(spec_2,grid_3)
call divide_by_cos(grid_3)



!Output #18
!############## printing #############################################
pr = 18
call mpp_sync(); if(root_pe) write(*,10000) pr
if(root_pe) write(*,1300)
do p = 0, npes-1
  call mpp_sync()
  if(pe == p) then
    do j = js, je
      write(*,1301) pe, j, grid_2(is,j), grid_3(is,j)
    enddo
  endif
enddo
1300 format(1x,/,1x,'checking latitudinal derivative of cos^10',/)
1301 format(1x, i2, 1x, i3, 2(1x, f12.8))
!############## printing #############################################



! compute_lon_deriv_cos WORKS IN THE SAME WAY;  ONE CAN ALSO USE 
! compute_gradient_cos(spec, ucos, vcos) TO DO BOTH SIMULTANEOUSLY



!  COMPUTING THE VORTICITY AND DIVEGENCE IS A LITTLE TRICKY
! LET'S COMPUTE THE LAPLACIAN IN TWO DIFFERENT WAYS

! THINK OF GRID_1 BELOW AS THE VELOCITY POTENTIAL FOR A 2D FLOW
! U_COS, V_COS ARE THE SPECTRAL FIELDS FOR COS(LAT) X VELOCITY
! U_COSG, V_COSG ARE THE CORRESPONDING GRID FIELDS
! THE INPUT NEEDED FOR THE COMPUTATION OF DIVERGENCE (AND VORTICITY)
! ARE THE SPECTRAL COMPONENTS OF U/COS AND V/COS !!!
! THE CALL TO DIVIDE_BY_COS2 COMES IN HANDY HERE


allocate(u_cos(ms:me,ns:ne))
allocate(v_cos(ms:me,ns:ne))
allocate(u_cosm(ms:me,ns:ne))
allocate(v_cosm(ms:me,ns:ne))
allocate(u_cosg(is:ie,js:je))
allocate(v_cosg(is:ie,js:je))


!  COMPUTING(DIVERGENCE(GRADIENT)): 

call trans_grid_to_spherical(grid_1, spec_1, do_truncation = .true.)
call compute_gradient_cos(spec_1, u_cos, v_cos)

call trans_spherical_to_grid(u_cos,u_cosg)
call trans_spherical_to_grid(v_cos,v_cosg)

call divide_by_cos2(u_cosg)
call divide_by_cos2(v_cosg)

call trans_grid_to_spherical(u_cosg,u_cosm)
call trans_grid_to_spherical(v_cosg,v_cosm)

spec_2 = compute_div(u_cosm, v_cosm)

! COMPUTING LAPLACIAN DIRECTLY IS SPHERICAL DOMAIN
!  (THIS JUST MULTIPLIES BY THE EIGENVALUES)

spec_3 = compute_laplacian(spec_1, power = 1)

!  POWER = 1 IS ACTUALLY THE DEFAULT
!  POWER = -1 IS THE INVERSE LAPLACIAN, POWER = 2 IS DEL**4 ETC

call trans_spherical_to_grid(spec_2, grid_2)
call trans_spherical_to_grid(spec_3, grid_3)


!Output #19
!############## printing #############################################
pr = 19
call mpp_sync(); if(root_pe) write(*,10000) pr
call mpp_sync(); if(root_pe) write(*,1302) 
do p = 0, npes-1
  call mpp_sync()
  if(pe == p) then
    do j = js, je
      write(*,1301) pe, j, grid_2(is,j), grid_3(is,j)
    enddo
  endif
enddo

1302 format(1x,'checking divergence ',/, &
1x, 'by computing laplacian in two different ways',/)
!############## printing #############################################


!  ONE CAN DO THE SAME THING FOR THE VORTICITY (compute_vor(u_cosm))
!  OR ONE CAN COMPUTE VORTICITY AND DIVERGENCE SIMULTANEOUSLY
!  call compute_vort_div(u_cosm, v_cosm, vor, div)



! WE ALSO NEED TO BE ABLE TO REGENERATE THE COMPONENTS OF A VELOCITY
!  FIELD FROM THE VORTICITY AND DIVERGENCE

allocate (vorg  (is:ie,js:je))
allocate (divg  (is:ie,js:je))
allocate (vorg1 (is:ie,js:je))
allocate (divg1 (is:ie,js:je))
allocate (ug    (is:ie,js:je))
allocate (vg    (is:ie,js:je))
allocate (vor   (ms:me,ns:ne))
allocate (div   (ms:me,ns:ne))
allocate (vor1  (ms:me,ns:ne))
allocate (div1  (ms:me,ns:ne))

! START BY GENERATING SOME ARBITRARY GRID FIELDS THAT WE CALL VORTICITY 
!   AND DIVERGENCE

vorg = 0.0
divg = 0.0
if(is .le. 10 .and. ie .ge. 10 .and. js .le. 4 .and. je .ge. 4) &
      vorg(10,4) = 1.0
if(is .le. 13 .and. ie .ge. 13 .and. js .le. 9 .and. je .ge. 9) &
      divg(13,9) = 1.0

call trans_grid_to_spherical(vorg, vor, do_truncation = .true.)
call trans_grid_to_spherical(divg, div, do_truncation = .true.)

if(ms == 0 .and. ns == 0) vor(0,0) = 0.0
if(ms == 0 .and. ns == 0) div(0,0) = 0.0

! NOTE THAT WE ARE CAREFUL TO TRUNCATE THEM, AND TO SET THE GLOBAL 
!  MEAN TO ZERO, AS MUST BE THE CASE FOR A REALIZEABLE VOR AND DIV

! WE ALSO GENERATE THE GRID FIELDS CONSISTENT WITH THESE TRUNCATED
!   SPECTRAL FIELDS

call trans_spherical_to_grid(vor, vorg)
call trans_spherical_to_grid(div, divg)


! THE CALL compute_ucos_vos GENERATES THE SPECTRAL COMPONENTS OF THE
!  VELOCITY FIELD MULTIPLIED BY COS(LAT)

call compute_ucos_vcos          (vor, div, u_cos, v_cos)
call trans_spherical_to_grid    (u_cos, ug)
call trans_spherical_to_grid    (v_cos, vg)
call divide_by_cos              (ug)
call divide_by_cos              (vg)

!  At this point we have the velocity field in the grid domain

! TO RETURN TO (VOR, DIV) WE FIRST HAVE TO DIVIDE BY COS(LAT) AGAIN

call divide_by_cos(ug)
call divide_by_cos(vg)
call trans_grid_to_spherical(ug, u_cosm, do_truncation=.false.)
call trans_grid_to_spherical(vg, v_cosm, do_truncation=.false.)

! NOTE THAT WHEN WE DO NOT TRUNCATE THE SPECTRAL U/COS AND V/COS !!
! YOU NEED TO RETAIN ONE HIGHER MERIDIONAL HARMONIC IN U/COS AND V/COS TO 
! REGENERATE THE CORRECT VOR AND DIV.

call compute_vor_div(u_cosm, v_cosm, vor1, div1)

! NOW WE CAN TRUNCATE

if(triang_trunc) then
  call triangular_truncation(vor1)
  call triangular_truncation(div1)
else
  call rhomboidal_truncation(vor1)
  call rhomboidal_truncation(div1)
endif

call trans_spherical_to_grid(vor1, vorg1)
call trans_spherical_to_grid(div1, divg1)

do j = js,je
  do i = is,ie
    work_g(i,j) = wts_lat(j)*(vorg1(i,j)-vorg(i,j))**2
  enddo
enddo
e1 = mpp_global_sum(grid_domain, work_g, BITWISE_EXACT_SUM)/float(num_lon)
do j = js,je
  do i = is,ie
    work_g(i,j) = wts_lat(j)*(divg1(i,j)-divg(i,j))**2
  enddo
enddo
e2 = mpp_global_sum(grid_domain, work_g, BITWISE_EXACT_SUM)/float(num_lon)

! LETS DO THE WHOLE THING OVER BUT MAKE THE MISTAKE OF TRUNCATING 
!  U/COS AND V/COS 

call trans_grid_to_spherical(ug, u_cosm, do_truncation=.true.)
call trans_grid_to_spherical(vg, v_cosm, do_truncation=.true.)
call compute_vor_div(u_cosm, v_cosm, vor1, div1)
if(triang_trunc) then
  call triangular_truncation(vor1)
  call triangular_truncation(div1)
else
  call rhomboidal_truncation(vor1)
  call rhomboidal_truncation(div1)
endif
call trans_spherical_to_grid(vor1, vorg1)
call trans_spherical_to_grid(div1, divg1)

do j = js,je
  do i = is,ie
    work_g(i,j) = wts_lat(j)*(vorg1(i,j)-vorg(i,j))**2
  enddo
enddo
e3 = mpp_global_sum(grid_domain, work_g, BITWISE_EXACT_SUM)/float(num_lon)
do j = js,je
  do i = is,ie
    work_g(i,j) = wts_lat(j)*(divg1(i,j)-divg(i,j))**2
  enddo
enddo
e4 = mpp_global_sum(grid_domain, work_g, BITWISE_EXACT_SUM)/float(num_lon)

! SINCE THIS IS A BIT INVOLVED, AN EASIER WAY IS PROVIDED 

! TO GET GRID VELOCITY FIELDS FROM SPECTRAL VORTICITY AND DIVERGENCE

call uv_grid_from_vor_div(vor, div, ug, vg)

! TO GET SPECTRAL VORTICITY AND DIVERGENCE FROM GRID U AND V

call vor_div_from_uv_grid(ug, vg, vor1, div1)

call trans_spherical_to_grid(vor1, vorg1)
call trans_spherical_to_grid(div1, divg1)

do j = js,je
  do i = is,ie
    work_g(i,j) = wts_lat(j)*(vorg1(i,j)-vorg(i,j))**2
  enddo
enddo
e5 = mpp_global_sum(grid_domain, work_g, BITWISE_EXACT_SUM)/float(num_lon)
do j = js,je
  do i = is,ie
    work_g(i,j) = wts_lat(j)*(divg1(i,j)-divg(i,j))**2
  enddo
enddo
e6 = mpp_global_sum(grid_domain, work_g, BITWISE_EXACT_SUM)/float(num_lon)

!Output #20
!############## printing #############################################
pr = 20
call mpp_sync(); if(root_pe) write(*,10000) pr
call mpp_sync(); if(root_pe) write(*,1350) e1, e2, e3, e4, e5, e6
1350 format(1x,'transforming (vor,div) to (u,v) and back',//,                       &
 1x,'the correct way:',/,1x,                                                        &
 'variance of error in vor = ', e13.6, 3x, 'variance of error in div = ', e13.6,//, &
 1x,'the incorrect way:',/,1x,                                                      &
 'variance of error in vor = ', e13.6, 3x, 'variance of error in div = ', e13.6,//, &
 1x,'using the shortcut:',/,1x,                                                     &
 'variance of error in vor = ', e13.6, 3x, 'variance of error in div = ', e13.6)
!############## printing #############################################




! IN FLUID MECHANICS WE ARE OFTEN COMPUTING THE HORIZONTAL ADVECTION OF A
! TRACER :  DT/DT (THIS IS A PARTIAL TIME DERIVATIVE) = -V.GRAD(T)
! SO WE PROVIDE A SIMPLE CALL TO COMPUTE -V.GRAD(T) .  TO CHECK IT, LET'S
! COMPUTE A NON-DIVERGENT VELOCITY FIELD FROM A STREAMFUNCTION, AND THEN
! TRY TO ADVECT THE STREAMFUNCTION ITSELF -- THIS SHOULD GIVE ZERO SINCE
! THE VECLOTY FIELD WILL BE PERPENDICULAR TO THE LINES OF CONSTANT STREAMFUNCTION

! SETTTING UP OUR FIELDS

allocate(streamg (is:ie, js:je))
allocate(tendg   (is:ie, js:je))
allocate(stream  (ms:me, ns:ne))
do j = js, je
  do i = is, ie
    streamg(i,j) = float(i)*float(j) ! an arbitrary field
  enddo
enddo
call trans_grid_to_spherical(streamg, stream, do_truncation = .true.)
vor = compute_laplacian(stream)
div = (0.,0.)
call uv_grid_from_vor_div(vor, div, ug, vg)

! NOW COMPUTE THE ADVECTION -- THIS CALL ACTUALLY ADDS THIS TERM TO AN EXISTING
!   TENDENCY FIELD -- ITS INPUT IS THE SPECTRAL FIELD OF THE TRACER AND THE 
!   GRID VELOCITY FIELDS  (A BIT SPECIALIZED ADMITTEDLY)

tendg = 0.0

call horizontal_advection(stream, ug, vg, tendg)

! THIS IS EQUIVALENT TO
!   call compute gradient(stream, spec_1, spec_2)
!   call trans_spherical_to_grid(spec_1, grid_1)
!   call trans_spherical_tp_grid(spec_2, grid_2)
!   call divide_by_cos(grid_1)
!   call divide_by_cos(grid_2)
!   tendg = tendg -  ug*grid_1 - vg*grid_2

do j = js,je
  do i = is,ie
    work_g(i,j) = wts_lat(j)*tendg(i,j)**2
  enddo
enddo
error = mpp_global_sum(grid_domain, work_g, BITWISE_EXACT_SUM)/float(num_lon)

!Output #21
!############## printing #############################################
pr = 21
call mpp_sync(); if(root_pe) write(*,10000) pr
call mpp_sync(); if(root_pe) write(*,1360) error
1360 format(1x,'advecting the streamfunction',//, 1x, &
     'variance of error = ', e13.6)
!############## printing #############################################





! FINALLY, ALL CALLS SUCH AS TRIANGULAR_TRUNCATION, COMPUTE_LAT_DERIV_COS,
!  COMPUTE_LAPLACIAN, TRANS_GRID_TO_SPHERICAL, DIVIDE_BY_COS, 
!  UV_GRID_FROM_VOR_DIV, HORIZONTAL_ADVECTION, ETC
!  WILL ALSO WORK IF THEY ARE CALLED WITH THREE DIMENSIONAL ARRAYS
!  IE, TO TRANSFORM 7 FIELDS SIMULTANEOUSLY, USE

allocate(ngrid(is:ie, js:je, 7))
allocate(nspec(ms:me, ns:ne, 7))

ngrid = 1.0
call trans_grid_to_spherical(ngrid, nspec)
call trans_spherical_to_grid(nspec, ngrid)

!Output #22
!############## printing #############################################
pr = 22
call mpp_sync(); if(root_pe) write(*,10000) pr
call mpp_sync(); if(root_pe) write(*,1370) (ngrid(is,js,k), k = 1,7)
1370 format(1x, &
 'transforming 7 fields simultaneously ',     &
 ' -- just checking if this flows properly ', /,&
  /,7(1x,f10.6))
call mpp_sync(); if(root_pe) write(*,10000)
!############## printing #############################################



if(check_checking) then
  spec_1 = (0.,0.)
  if(ms == 0) spec_1(ms,:) = (0.,1.)  ! this should crash the program
  call trans_spherical_to_grid(spec_1, grid_1)
endif
!======================================================================

call transforms_end

call fms_end()

end program tutorial

!======================================================================
