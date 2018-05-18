module microphysics_mod

use fms_mod,                      only :  error_mesg, FATAL, mpp_pe,   &
                                          mpp_root_pe, open_namelist_file, &
                                          check_nml_error, close_file,  &
                                          write_version_number, file_exist,&
                                          stdlog
use mpp_mod,                      only :  input_nml_file
use constants_mod,                only :  cp_air, hlv, hlf, hls, tfreeze,  &
                                          rdgas, grav, rvgas
use rotstayn_klein_mp_mod,        only :  rotstayn_klein_microp, &
                                          rotstayn_klein_microp_init,  &
                                          rotstayn_klein_microp_end
use morrison_gettelman_microp_mod, only : morrison_gettelman_microp,  &
                                          morrison_gettelman_microp_init, &
                                          morrison_gettelman_microp_end
use strat_cloud_utilities_mod,    only :  strat_cloud_utilities_init, &
                                          diag_id_type, diag_pt_type, &
                                          strat_nml_type, atmos_state_type,&
                                          cloud_state_type, particles_type,&
                                          strat_constants_type, &
                                          cloud_processes_type, &
                                          precip_state_type
use cldwat2m_micro_mod,           only :  ini_micro, mmicro_pcond,  &
                                          mmicro_end
use micro_mg_mod,              only:   micro_mg_init, micro_mg_get_cols,&
                                          micro_mg_tend

implicit none
private

!------------------------------------------------------------------------
!---interfaces-----------------------------------------------------------

public  microphysics, microphysics_init, microphysics_end    



!------------------------------------------------------------------------
!---version number-------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

!--------------------------------------------------------------------------
!---namelist---------------------------------------------------------------
 
! these variables used with ncar microphysics:
real     ::         rhmini=0.80  ! minimum rh for ice cld fraction > 0
logical  ::         microp_uniform = .false. 
                               ! .true. = configure uniform for sub-columns 
                               ! .false. = use w/o sub-columns (standard)

logical  ::         do_cldice = .true.
                               ! .true. = do all processes (standard)
                               ! .false. = skip all processes affecting
                               !           cloud ice

namelist / microphysics_nml /  rhmini, microp_uniform, do_cldice

!-------------------------------------------------------------------------

integer, parameter :: r8 = selected_real_kind(12)   ! 8 byte real
logical            :: module_is_initialized = .false.


CONTAINS




!##########################################################################

subroutine microphysics_init (Nml, Constants)

type(strat_nml_type), intent(in) :: Nml
type(strat_constants_type), intent(in) :: Constants

      integer :: unit, io, ierr, logunit
      character(len=128) :: errstring ! Output status: non-blank for 
                                      ! error return

      if (module_is_initialized) return

!-------------------------------------------------------------------------
!    process namelist.
!-------------------------------------------------------------------------
#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=microphysics_nml, iostat=io)
     ierr = check_nml_error(io,'microphysics_nml')
#else
     if ( file_exist('input.nml')) then
       unit = open_namelist_file ( )
       ierr=1; do while (ierr /= 0)
       read  (unit, nml=microphysics_nml, iostat=io, end=10)
       ierr = check_nml_error(io,'microphysics_nml')
       enddo
10      call close_file (unit)
     endif
#endif
!-------------------------------------------------------------------------
!    write version number and namelist to standard log.
!-------------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe()) &
                                  write (logunit, nml=microphysics_nml)

!-------------------------------------------------------------------------
!    make sure needed modules have been initialized.
!-------------------------------------------------------------------------
      call strat_cloud_utilities_init
      if (Constants%do_rk_microphys) then
        call rotstayn_klein_microp_init
      else if (Constants%do_mg_microphys) then
        call morrison_gettelman_microp_init (Nml%do_pdf_clouds, Nml%qcvar)
      else if (constants%do_mg_ncar_microphys) then
        call ini_micro (grav, rdgas, rvgas, cp_air,   &
                                         tfreeze, hlv, hlf, rhmini)
      else if (Constants%do_ncar_microphys) then
        call micro_mg_init (r8, grav, rdgas, rvgas, cp_air, tfreeze, &
                            hlv, hlf, rhmini, microp_uniform, do_cldice, &
                            errstring )
        if (trim(errstring) /= '') then
          call error_mesg ('strat_cloud/microphysics/micro_mg_init', &
                                                         errstring, FATAL)
        endif
      else
        call error_mesg ('strat_cloud/microphysics_init', &
              'invalid strat_cloud_nml microphys_scheme option', FATAL)
      endif

      module_is_initialized = .true.

end subroutine microphysics_init


!##########################################################################

subroutine microphysics &         
                    (idim, jdim, kdim, relvarn, Nml, Constants, N3D, Atmos_state, &
                     Cloud_state, Cloud_processes, Particles, n_diag_4d, &
                     diag_4d, diag_id, diag_pt, n_diag_4d_kp1, diag_4d_kp1,&
                     ST_out, SQ_out, Precip_state, otun, ncall, &
                     qa_upd_0, SA_0, nrefuse, isamp, jsamp, ksamp, &
                     debugo, debugo0, debugo1)    

!------------------------------------------------------------------------
type(strat_constants_type),        intent(inout) :: Constants     
type(atmos_state_type),            intent(inout) :: Atmos_state
type(cloud_state_type),            intent(inout) :: Cloud_state
type(precip_state_type),           intent(inout) :: Precip_state
type(cloud_processes_type),        intent(inout) :: Cloud_processes
type(particles_type),              intent(inout) :: Particles
integer,                           intent(inout) :: nrefuse         
integer,                           intent(in)    :: idim, jdim, kdim, &
                                                    n_diag_4d, otun, ncall,&
                                                    n_diag_4d_kp1, isamp, &
                                                    jsamp, ksamp
type(strat_nml_type),              intent(in)    :: Nml
logical,                           intent(in)    :: debugo, debugo0, debugo1
real, dimension (idim,jdim,kdim),  intent(in)    :: N3D, relvarn 
 real, dimension(idim,jdim,kdim),  intent(inout) :: ST_out,SQ_out,    &
                                                    qa_upd_0, SA_0   
real, dimension(idim,jdim,kdim,0:n_diag_4d),                &
                                   intent(inout) :: diag_4d
real, dimension(idim,jdim,kdim+1,0:n_diag_4d_kp1),               &
                                   intent(inout) :: diag_4d_kp1
type(diag_id_type),                intent(in)    :: diag_id
type(diag_pt_type),                intent(inout) :: diag_pt

!------------------------------------------------------------------------
!---local variables------------------------------------------------------

      real, dimension(idim,jdim,kdim,4) :: rbar_dust_4bin, ndust_4bin
      real, dimension(idim,jdim,kdim)   :: ST_micro, SQ_micro, SL_micro, &
                                           SI_micro, SN_micro, SNI_micro
      real, dimension(idim,jdim,kdim)   :: D_eros_l, D_eros_i,  &
                                           nerosc, nerosi,    &
                                           dqcdt, dqidt, &
                                           qa_new, &
                                           ssat_disposal, &
                                           ql_new,  qi_new,              &
                                           nctend, nitend, qn_new, qni_new
      real, dimension(idim,jdim,kdim)   :: rho, liqcldf, icecldf, tmp2s
      real, dimension(idim,jdim,kdim)   :: accre_enhann, tnd_qsnown, &
                                           tnd_nsnown, re_icen           
      real, dimension(idim,jdim)        :: m1, m2, scalef
      integer,dimension(:),allocatable  :: mgcols         
      integer                           :: mgncol
      integer                           :: i,j,k,n
      integer                           :: top_lev, nlev
      character(len=128)                :: errstring


!------------------------------------------------------------------------
!    call selected microphysics scheme.
!------------------------------------------------------------------------
      if (Constants%do_rk_microphys ) then

!------------------------------------------------------------------------
!    Rotstayn-Klein microphysics
!------------------------------------------------------------------------
        call rotstayn_klein_microp ( &
                         idim, jdim, kdim,  Nml, N3D,     &
! cjg: total activation for RK
                          Constants%total_activation,  &  
                         Constants%overlap, Constants%dtcloud,  &
                         Constants%inv_dtcloud, Atmos_state%pfull,&
                         Atmos_state%deltpg, Atmos_state%airdens,     &
                         Constants%mask_present, Constants%mask, &
                         Atmos_state%esat0, Cloud_state%ql_in,  &
                         Cloud_state%qi_in, Cloud_state%qa_in,   &
                         Cloud_state%ql_mean, Cloud_state%qa_mean, &
                         Cloud_state%qn_mean, Atmos_state%omega,  &
                         Atmos_state%T_in, Atmos_state%U_ca, &
                         Atmos_state%qv_in, Atmos_state%qs,  &
                         Cloud_processes%D_eros, Cloud_processes%dcond_ls, &
                         Cloud_processes% dcond_ls_ice,         &
                         Cloud_processes%qvg, Atmos_state%gamma,   &
                         Cloud_processes%delta_cf, Particles%drop1,    &
                         Particles%concen_dust_sub, Cloud_state%ql_upd,   &
                         Cloud_state%qi_upd, Cloud_state%qn_upd,       & 
                         Cloud_state%qi_mean, Cloud_state%qa_upd,   &
                         Atmos_state%ahuco, n_diag_4d, diag_4d, diag_id, &
                         diag_pt, n_diag_4d_kp1, diag_4d_kp1,  &
                         Constants%limit_conv_cloud_frac, &
                         Cloud_state%SA_out, Cloud_state%SN_out,        & 
                         ST_out, SQ_out, Cloud_state%SL_out,  &
                         Cloud_state%SI_out, Precip_state%rain3d,  &
                         Precip_state%snow3d, Precip_state%snowclr3d,    &
                         Precip_state%surfrain, Precip_state%surfsnow,  &
                         Cloud_processes%f_snow_berg, otun)                 
      else   !if (Constants%do_rk_microphys ) 

!--------------------------------------------------------------------------
!     for morrison-gettelman and NCAR, some additional fields are needed. 
!     for droplet activation Yi's drop1 is used. 
!--------------------------------------------------------------------------
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              Particles%drop2(i,j,k) = Particles%drop1(i,j,k)*  &
                                          1.e6/Atmos_state%airdens(i,j,k)
              Cloud_processes%dcond_ls_tot(i,j,k) =   &
                               Cloud_processes%dcond_ls(i,j,k) +   &
                                      Cloud_processes%dcond_ls_ice(i,j,k) 
              Atmos_state%tn(i,j,k) = Atmos_state%T_in(i,j,k) +   &
                                                             ST_out(i,j,k)
              Atmos_state%qvn(i,j,k) = Atmos_state%qv_in(i,j,k) + &
                                                              SQ_out(i,j,k)
              if (Constants%tiedtke_macrophysics) then
                D_eros_i(i,j,k) = -Cloud_state%qi_upd(i,j,k)* &
                                        Cloud_processes%D_eros(i,j,k)/ &
                                                     Constants%dtcloud
                D_eros_l(i,j,k) = -Cloud_state%ql_upd(i,j,k)* &
                                        Cloud_processes%D_eros(i,j,k)/ &
                                                     Constants%dtcloud
                if (Cloud_state%ql_upd(i,j,k) >= Nml%qmin) then
                  nerosc(i,j,k) = D_eros_l(i,j,k)/  &
                                      Cloud_state%ql_upd(i,j,k)* &
                               Cloud_state%qn_upd(i,j,k)/MAX(0.0001, &
                                               Cloud_state%qa_upd(i,j,k))
                else
                  nerosc(i,j,k) = 0.
                endif
                if (Cloud_state%qi_upd(i,j,k) >= Nml%qmin) then
                  nerosi(i,j,k) = D_eros_i(i,j,k)/   &
                                          Cloud_state%qi_upd(i,j,k)* &
                               Cloud_state%qni_upd(i,j,k)/MAX(0.0001, &
                                              Cloud_state%qa_upd(i,j,k))
                else
                  nerosi(i,j,k) = 0.
                endif
                if (Cloud_processes%dcond_ls_tot(i,j,k) > 0.) then
                  if (Atmos_state%tn(i,j,k) <= (tfreeze - 40.) ) then
                    dqcdt (i,j,k) = 0.
                    dqidt(i,j,k) = Cloud_processes%dcond_ls_tot(i,j,k)* &
                                                 Constants%inv_dtcloud
                  else
                    dqidt (i,j,k) = 0.
                    dqcdt(i,j,k) = Cloud_processes%dcond_ls_tot(i,j,k)*  &
                                                 Constants%inv_dtcloud
                  endif
                else
                  if (Atmos_state%tn(i,j,k) <= tfreeze) then
                    dqcdt(i,j,k) = MAX(Cloud_processes%dcond_ls_tot(i,j,k),&
                              -Cloud_state%ql_upd(i,j,k))
                    dqidt(i,j,k) = MAX(Cloud_processes%dcond_ls_tot(i,j,k) &
                              - dqcdt(i,j,k), -Cloud_state%qi_upd(i,j,k))
                    dqcdt(i,j,k) = dqcdt(i,j,k)* Constants%inv_dtcloud
                    dqidt(i,j,k) = dqidt(i,j,k)* Constants%inv_dtcloud
                  else
                    dqidt(i,j,k) = 0.
                    dqcdt(i,j,k) = MAX(Cloud_processes%dcond_ls_tot(i,j,k),&
                                           -Cloud_state%ql_upd(i,j,k))* &
                                                 Constants%inv_dtcloud
                  endif
                endif
              else
                dqidt(i,j,k) = 0.
                dqcdt(i,j,k) = 0.
                nerosi(i,j,k) = 0.
                nerosc(i,j,k) = 0.
                D_eros_l(i,j,k) = 0.                             
                D_eros_i(i,j,k) = 0.                             
              endif
            end do
          end do   
        end do   

        if (Constants%do_mg_microphys) then

          do j=1,jdim

!------------------------------------------------------------------------
!    if debugging is activated, output the temp tendency prior to
!    microphysics.
!------------------------------------------------------------------------
            if (debugo) then
              if ( j .eq. jsamp) then
                write(otun, *) " ST samp bef mg ", ST_out(isamp,jsamp,ksamp)
              end if
            endif
!-------------------------------------------------------------------------
!    call morrison-gettelman microphysics package.
!-------------------------------------------------------------------------
            call morrison_gettelman_microp( &
                  Constants%tiedtke_macrophysics, &
                  Constants%total_activation, Constants%dqa_activation, &
                  ncall, j ,idim, jdim, kdim, Nml, &
                  Constants%dtcloud, Atmos_state%pfull(:,j,:),  &
                  Atmos_state%delp(:,j,:),                          &
                  Atmos_state%tn(:,j,:),  Atmos_state%T_in(:,j,:),    &
                  Atmos_state%qvn(:,j,:), Atmos_state%qv_in(:,j,:),  &
                  Cloud_state%ql_upd(:,j,:), Cloud_state%qi_upd(:,j,:), &
                  Cloud_state%qn_upd(:,j,:), Cloud_state%qni_upd(:,j,:), &
                  Cloud_state%qa_upd(:,j,:),  &
                  dqcdt(:,j,:), dqidt(:,j,:), &
                  Particles%drop2(:,j,:), Particles%crystal1(:,j,:), &
                  Particles%rbar_dust(:,j,:), Particles%ndust(:,j,:),    &
                  Cloud_processes%delta_cf(:,j,:),    &
                  Cloud_state%qa_upd(:,j,:), &
                  qa_upd_0(:,j,:), SA_0(:,j,:), D_eros_l(:,j,:),  &
                  nerosc(:,j,:),  D_eros_i(:,j,:), nerosi(:,j,:), &
                  Atmos_state%gamma(:,j,:), &
                  Constants%inv_dtcloud, Cloud_state%qa_in(:,j,:),     &
                  ST_out(:,j,:), SQ_out(:,j,:), ssat_disposal(:,j,:), &
                  ST_micro(:,j,:), SQ_micro(:,j,:),&
                  SL_micro(:,j,:), SI_micro(:,j,:),  &
                  SN_micro(:,j,:), SNI_micro(:,j,:),&
                  Cloud_state%SA_out(:,j,:),  &
                  Precip_state%rain3d, Precip_state%snow3d,  &
                  Precip_state%surfrain(:,j), Precip_state%surfsnow(:,j), &
                  Precip_state%lsc_rain(:,j,:),   &
                  Precip_state%lsc_snow(:,j,:), &
                  Precip_state%lsc_rain_size(:,j,:),   &
                  Precip_state%lsc_snow_size(:,j,:), &
                  Cloud_processes%f_snow_berg(:,j,:), &
                  n_diag_4d, diag_4d, diag_id, &
                  diag_pt, nrefuse, debugo0, debugo1, otun)    

!------------------------------------------------------------------------
!    if debugging is activated, output the temp tendency after microphysics.
!------------------------------------------------------------------------
            IF (debugo) THEN
              if ( j .eq. jsamp) then
                write(otun, *) " ST samp aft mg ", ST_out(isamp,jsamp,ksamp)
              end if
            END IF
          end do   ! j loop

        else if (Constants%do_mg_ncar_microphys) then  !(do_mg_microphys) 
!-------------------------------------------------------------------------
!    use ncar maintained version of microphysics
!-------------------------------------------------------------------------
!  are these the actual bin centers, or arbitrary values ??
          rbar_dust_4bin(:,:,:,1) = 5.0e-6
          rbar_dust_4bin(:,:,:,2) = 10.0e-6
          rbar_dust_4bin(:,:,:,3) = 15.0e-6
          rbar_dust_4bin(:,:,:,4) = 20.0e-6
          ndust_4bin = 0.0
!         ndust_4bin(:,j,:,1) = 0.25*Particles%ndust(:,j,:) 
!         ndust_4bin(:,j,:,2) = 0.25*Particles%ndust(:,j,:) 
!         ndust_4bin(:,j,:,3) = 0.25*Particles%ndust(:,j,:) 
!         ndust_4bin(:,j,:,4) = 0.25*Particles%ndust(:,j,:) 

          

!-------------------------------------------------------------------------
!    call microphysics package as maintained at NCAR.
!-------------------------------------------------------------------------
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                rho(i,j,k) = Atmos_state%pfull(i,j,k)/   &
                                            (rdgas*Atmos_state%tn(i,j,k))
                liqcldf(i,j,k) = Cloud_state%qa_upd(i,j,k)
                icecldf(i,j,k) = Cloud_state%qa_upd(i,j,k)
              end do
            end do
          end do

          do j=1,jdim
            call mmicro_pcond( &
                   Constants%dqa_activation, &
                   Constants%total_activation, &
                   Constants%tiedtke_macrophysics, &
                   .false., j ,jdim, kdim, idim, idim,  &
                   Constants%dtcloud,    &
                   relvarn(:,j,:), &
                   Atmos_state%tn(:,j,:),     &
                   Atmos_state%qvn(:,j,:),   &
                   Cloud_state%ql_upd(:,j,:), Cloud_state%qi_upd(:,j,:), &
                   Cloud_state%qn_upd(:,j,:), Cloud_state%qni_upd(:,j,:), &
                   Atmos_state%pfull(:,j,:),  &
                   Atmos_state%delp(:,j,:),   &
                   Atmos_state%phalf(:,j,:), &
                   Cloud_state%qa_upd(:,j,:),  &
                   liqcldf(:,j,:)       , icecldf(:,j,:),   & 
                   Cloud_processes%delta_cf(:,j,:), &
                   D_eros_l(:,j,:), nerosc(:,j,:), &
                   D_eros_i(:,j,:), nerosi(:,j,:), &
                   dqcdt(:,j,:), dqidt(:,j,:), &
                   Particles%crystal1(:,j,:)/rho(:,j,:), &
                   Particles%drop2(:,j,:),   &
                   rbar_dust_4bin(:,j,:,:), ndust_4bin(:,j,:,:), &
                   ST_micro(:,j,:), SQ_micro(:,j,:), SL_micro(:,j,:), &
                   SI_micro(:,j,:), SN_micro(:,j,:), SNI_micro(:,j,:), &
                   Precip_state%surfrain(:,j),   &
                   Precip_state%surfsnow(:,j),   &
                   Precip_state%rain3d(:,j,:),   &
                   Precip_state%snow3d(:,j,:),   &
                   Precip_state%lsc_rain(:,j,:),   &
                   Precip_state%lsc_snow(:,j,:), &
                   Precip_state%lsc_rain_size(:,j,:),   &
                   Precip_state%lsc_snow_size(:,j,:), &
                   Cloud_processes%f_snow_berg(:,j,:), &

                   Nml, Cloud_state%qa_in(:,j,:), Atmos_state%gamma(:,j,:),&
                   SA_0(:,j,:), Cloud_state%SA_out(:,j,:),  &

                   ssat_disposal (:,j,:), &
                   n_diag_4d, diag_4d, diag_id, &
                   diag_pt)    
          end do

        else if (Constants%do_ncar_microphys) then  !(do_ncar_microphys) 
!-------------------------------------------------------------------------
!    use ncar maintained version of microphysics
!-------------------------------------------------------------------------
!  are these the actual bin centers, or arbitrary values ??
          rbar_dust_4bin(:,:,:,1) = 5.0e-6
          rbar_dust_4bin(:,:,:,2) = 10.0e-6
          rbar_dust_4bin(:,:,:,3) = 15.0e-6
          rbar_dust_4bin(:,:,:,4) = 20.0e-6
          ndust_4bin = 0.0
!         ndust_4bin(:,j,:,1) = 0.25*Particles%ndust(:,j,:) 
!         ndust_4bin(:,j,:,2) = 0.25*Particles%ndust(:,j,:) 
!         ndust_4bin(:,j,:,3) = 0.25*Particles%ndust(:,j,:) 
!         ndust_4bin(:,j,:,4) = 0.25*Particles%ndust(:,j,:) 

          

!-------------------------------------------------------------------------
!    call microphysics package as maintained at NCAR.
!-------------------------------------------------------------------------
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                rho(i,j,k) = Atmos_state%pfull(i,j,k)/   &
                                            (rdgas*Atmos_state%tn(i,j,k))
                liqcldf(i,j,k) = Cloud_state%qa_upd(i,j,k)
                icecldf(i,j,k) = Cloud_state%qa_upd(i,j,k)
              end do
            end do
          end do


          nlev = kdim
          top_lev = 1

          do j=1,jdim
            call micro_mg_get_cols (idim, nlev, top_lev, &
               Atmos_state%qvn(:,j,:), &
               Cloud_state%ql_upd(:,j,:) + dqcdt(:,j,:)*Constants%dtcloud, &
               Cloud_state%qi_upd(:,j,:) + dqidt(:,j,:)*Constants%dtcloud, &
                                                    mgncol, mgcols, .false.)

             accre_enhann(:,j,:) = 1.0  ! accretion enhancement factor
             tnd_qsnown(:,j,:) = 0.     
             tnd_nsnown(:,j,:) = 0.
             re_icen(:,j,:) = 0.

             call  micro_mg_tend ( &
                   Constants%dqa_activation, &
                   Constants%total_activation, &
                   Constants%tiedtke_macrophysics, &
                   j, jdim,   mgncol,   mgcols,   nlev,   top_lev,   &
                   Constants%dtcloud,  &
                   Atmos_state%tn(:,j,:),     &
                   Atmos_state%qvn(:,j,:),   &
                   Cloud_state%ql_upd(:,j,:), Cloud_state%qi_upd(:,j,:), &
                   Cloud_state%qn_upd(:,j,:), Cloud_state%qni_upd(:,j,:), &
                   relvarn(:,j,:), accre_enhann(:,j,:),                   &
                   Atmos_state%pfull(:,j,:),  &
                   Atmos_state%delp(:,j,:),   &
                   Atmos_state%phalf(:,j,:), &
                   Cloud_state%qa_upd(:,j,:),  &
                   liqcldf(:,j,:)       , icecldf(:,j,:),   & 
                   Cloud_processes%delta_cf(:,j,:), &
                   D_eros_l(:,j,:), nerosc(:,j,:), &
                   D_eros_i(:,j,:), nerosi(:,j,:), &
                   dqcdt(:,j,:), dqidt(:,j,:), &
                   Particles%crystal1(:,j,:)/rho(:,j,:), &
                   Particles%drop2(:,j,:),   &
                   rbar_dust_4bin(:,j,:,:), ndust_4bin(:,j,:,:), &
                   ST_micro(:,j,:), SQ_micro(:,j,:), SL_micro(:,j,:), &
                   SI_micro(:,j,:), SN_micro(:,j,:), SNI_micro(:,j,:), &
                   Precip_state%surfrain(:,j),   &
                   Precip_state%surfsnow(:,j),   &
                   Precip_state%lsc_snow(:,j,:), &
                   Precip_state%rain3d(:,j,:),   &
                   Precip_state%snow3d(:,j,:),   &
                   Precip_state%lsc_rain(:,j,:),   &
                   Precip_state%lsc_rain_size(:,j,:),  &
                   Precip_state%lsc_snow_size(:,j,:),   &
                   tnd_qsnown(:,j,:),  tnd_nsnown(:,j,:), re_icen(:,j,:), &
                   errstring, Cloud_processes%f_snow_berg(:,j,:), &
                   Nml,  ssat_disposal (:,j,:), &
                   n_diag_4d, diag_4d, diag_id, diag_pt)    

!Convert from effective radius to diameter for use in radiation.
! in old NCAR, diameter was returned from mmicro_pcond routine.
             Precip_state%lsc_rain_size(:,j,:) =     &
                                   2.0*Precip_state%lsc_rain_size(:,j,:)
             Precip_state%lsc_snow_size(:,j,:) =      &
                                  2.0*Precip_state%lsc_snow_size(:,j,:)

             if (trim(errstring) /= '') then
               call error_mesg ('strat_cloud/microphysics/micro_mg_tend', &
                                                         errstring, FATAL)
             endif
           end do

!-------------------------------------------------------------------------
!    exit with error if no valid microphysics scheme was specified.
!-------------------------------------------------------------------------
        else  
          call error_mesg ('strat_cloud/microphysics', &
              'invalid strat_cloud_nml microphys_scheme option', FATAL)
        endif    

!------ POST MICROPHYSICS ROUTINE CALL


        if (Nml%mass_cons) then
          do j=1,jdim
            do i=1,idim
              m1(i,j) = 0.
              do k=1,kdim
 
                m1(i,j) = m1(i,j) +   &
                   (SQ_micro(i,j,k) + SL_micro(i,j,k) + SI_micro(i,j,k))*  &
                             Constants%dtcloud*Atmos_state%delp(i,j,k)/grav
              end do
              m2(i,j) = 1.e3*Precip_state%surfrain(i,j)*Constants%dtcloud
              if (m2(i,j) .NE. 0.0) THEN
                scalef(i,j) = -m1(i,j)/m2(i,j)
 
                IF ( diag_id%rain_mass_conv > 0   ) &
                     diag_4d(i,j,1,diag_pt%rain_mass_conv) =   &
                       (scalef(i,j)*Precip_state%surfrain(i,j) -    &
                              Precip_state%surfrain(i,j))/Constants%dtcloud
                IF ( diag_id%snow_mass_conv > 0   ) &
                     diag_4d(i,j,1,diag_pt%snow_mass_conv) =   &
                        (scalef(i,j)*Precip_state%surfsnow(i,j) -  &
                              Precip_state%surfsnow(i,j))/Constants%dtcloud
 
                Precip_state%surfrain(i,j) =    &
                                     scalef(i,j)*Precip_state%surfrain(i,j)
                Precip_state%surfsnow(i,j) =    &
                                     scalef(i,j)*Precip_state%surfsnow(i,j)
              end if
            end do
          end do
        end if 

        if (diag_id%neg_rain > 0) &
              diag_4d(:,:,1,diag_pt%neg_rain) = 1.0e3*    &
               (Precip_state%surfrain(:,:) - Precip_state%surfsnow(:,:))* &
                                                        Constants%dtcloud
        if (diag_id%neg_snow > 0) &
             diag_4d(:,:,1,diag_pt%neg_snow) = 1.0e3*    &
                          (Precip_state%surfsnow(:,:))*Constants%dtcloud 

        Precip_state%surfrain = max(     &
           1.e3*(Precip_state%surfrain - Precip_state%surfsnow)*   &
                                                Constants%dtcloud , 0.0) 
        Precip_state%surfsnow = max(    &
              1.e3*Precip_state%surfsnow*Constants%dtcloud, 0.0)

        if (diag_id%neg_rain > 0) &
              diag_4d(:,:,1,diag_pt%neg_rain) =   &
                     -1.0*( (Precip_state%surfrain(:,:))  -   &
                                    diag_4d(:,:,1,diag_pt%neg_rain))/   &
                                                      (Constants%dtcloud) 
        if (diag_id%neg_snow > 0) &
             diag_4d(:,:,1,diag_pt%neg_snow) =   &
                      -1.0*( (Precip_state%surfsnow(:,:))  -   &
                                    diag_4d(:,:,1,diag_pt%neg_snow))/ &
                                                       (Constants%dtcloud) 

        ST_out = ST_out + ST_micro/cp_air*Constants%dtcloud
        SQ_out = SQ_out + SQ_micro*Constants%dtcloud
        Cloud_state%SL_out = Cloud_state%SL_out +   &
                                              SL_micro*Constants%dtcloud
        Cloud_state%SI_out = Cloud_state%SI_out +   &
                                              SI_micro*Constants%dtcloud
        Cloud_state%SN_out = Cloud_state%SN_out +   &
                                              SN_micro*Constants%dtcloud
        Cloud_state%SNI_out = Cloud_state%SNI_out +   &
                                              SNI_micro*Constants%dtcloud

!next time:  Get rid of nctend nitend.  simply update SN_out within the 
! loops below where nitend is updated.
        nctend = 0.
        nitend = 0.


!  tiedtke macrophysics --> cloud area adjustment required due to removing supersaturation 
        if (Constants%tiedtke_macrophysics .and.    &
                                         .not. Nml%do_pdf_clouds) then
          do k=1,kdim
            do j=1,jdim
              do i=1,idim
                rho(i,j,k) = Atmos_state%pfull(i,j,k)/   &
                                        (rdgas*Atmos_state%tn(i,j,k))
                if (Constants%limit_conv_cloud_frac) then
                  tmp2s(i,j,k) = Atmos_state%ahuco(i,j,k)
                else
                  tmp2s(i,j,k) = 0.
                endif
                if (ssat_disposal(i,j,k) == 2.) then
!  ming_activation --> additional aerosol activation proportional to change
!                      in cloud area due to supersaturation removal
                  nitend(i,j,k) = Particles%crystal1(i,j,k)/rho(i,j,k)*  &
                                  (1. - Cloud_state%qa_upd(i,j,k) -  &
                                           tmp2s(i,j,k))/Constants%dtcloud
                  if (diag_id%qnidt_super + diag_id%qni_super_col > 0 ) &
                      diag_4d(i,j,k,diag_pt%qnidt_super ) =    &
                         Particles%crystal1(i,j,k)/rho(i,j,k)*     &
                        (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))/  &
                                                          Constants%dtcloud
                else if (ssat_disposal(i,j,k) == 1.) then
                  nctend(i,j,k) = Particles%drop2(i,j,k)*   &
                                  (1. - Cloud_state%qa_upd(i,j,k) -   &
                                           tmp2s(i,j,k))/Constants%dtcloud
                  if (diag_id%qndt_super + diag_id%qn_super_col > 0 ) &
                      diag_4d(i,j,k,diag_pt%qndt_super ) =    &
                        Particles%drop2(i,j,k)*   &
                        (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))/  &
                                                           Constants%dtcloud
                endif
                if (ssat_disposal(i,j,k) > 0.0) then
                  if (max(diag_id%qadt_super,diag_id%qa_super_col) > 0) then
                    diag_4d(i,j,k,diag_pt%qadt_super ) =   &
                         (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))* &
                                                       Constants%inv_dtcloud
                  end if
                  Cloud_state%SA_out(i,j,k) = Cloud_state%SA_out(i,j,k) + &
                             (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))  
                  Cloud_state%qa_upd(i,j,k) = 1. - tmp2s(i,j,k)     
                endif
              end do
            end do
          end do
        end if 

        Cloud_state%SN_out = Cloud_state%SN_out + nctend*Constants%dtcloud
        Cloud_state%SNI_out = Cloud_state%SNI_out + nitend*Constants%dtcloud

!-----------------------------------------------------------------------
!       Cloud Destruction
!-----------------------------------------------------------------------
        ql_new  = Cloud_state%ql_in  + Cloud_state%SL_out 
        qi_new  = Cloud_state%qi_in  + Cloud_state%SI_out
        qn_new  = Cloud_state%qn_in  + Cloud_state%SN_out         
        qni_new = Cloud_state%qni_in + Cloud_state%SNi_out         
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              if ((ql_new(i,j,k) <= NMl%qmin  .and.    &
                   qi_new(i,j,k) <= Nml%qmin)   .or.   &
                  (Cloud_state%qa_upd(i,j,k) <= Nml%qmin)) then
                Cloud_state%SL_out(i,j,k) = Cloud_state%SL_out(i,j,k) -  &
                                                              ql_new(i,j,k)
                Cloud_state%SI_out(i,j,k) = Cloud_state%SI_out(i,j,k) -  &
                                                              qi_new(i,j,k)
                Cloud_state%SA_out(i,j,k) = Cloud_state%SA_out(i,j,k) -  &
                                                  Cloud_state%qa_upd(i,j,k)
                ST_out(i,j,k) = ST_out(i,j,k) -    &
                              (hlv*ql_new(i,j,k) + hls*qi_new(i,j,k))/cp_air
                SQ_out(i,j,k) = SQ_out(i,j,k) +    &
                                           (ql_new(i,j,k) + qi_new(i,j,k))
                Cloud_state%SN_out(i,j,k) = Cloud_state%SN_out(i,j,k) -  &
                                                             qn_new(i,j,k)
                Cloud_state%SNi_out(i,j,k) = Cloud_state%SNi_out(i,j,k) -  &
                                                            qni_new(i,j,k)
                if (diag_id%qldt_destr > 0 .or. diag_id%ql_destr_col > 0) &
                    diag_4d(i,j,k,diag_pt%qldt_destr) =    &
                                          - ql_new(i,j,k)/Constants%dtcloud
                if (diag_id%qidt_destr > 0 .or. diag_id%qi_destr_col > 0) &
                    diag_4d(i,j,k,diag_pt%qidt_destr) =    &
                                          - qi_new(i,j,k)/Constants%dtcloud
                if (diag_id%qadt_destr > 0 .or. diag_id%qa_destr_col > 0) &
                    diag_4d(i,j,k,diag_pt%qadt_destr) =     &
                              - Cloud_state%qa_upd(i,j,k)/Constants%dtcloud
                if (diag_id%qndt_destr > 0 .or. diag_id%qn_destr_col > 0) &
                    diag_4d(i,j,k,diag_pt%qndt_destr) =   &
                                        - qn_new(i,j,k)/Constants%dtcloud
                if (diag_id%qnidt_destr + diag_id%qni_destr_col > 0) &
                    diag_4d(i,j,k,diag_pt%qnidt_destr) =   &
                                      - qni_new(i,j,k)/Constants%dtcloud
                if (diag_id%qdt_destr + diag_id%q_destr_col > 0) &
                    diag_4d(i,j,k,diag_pt%qdt_destr) =   &
                         (ql_new(i,j,k) + qi_new(i,j,k))/Constants%dtcloud
              endif
            end do
          end do
        end do

        ql_new  =  Cloud_state%ql_in  + Cloud_state%SL_out 
        qi_new  =  Cloud_state%qi_in  + Cloud_state%SI_out
        qn_new  =  Cloud_state%qn_in  + Cloud_state%SN_out
        qni_new =  Cloud_state%qni_in + Cloud_state%SNI_out
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              if (abs(ql_new(i,j,k)) .le. Nml%qmin  .and.  &
                 Atmos_state%qv_in(i,j,k) + SQ_out(i,j,k) +   &
                                                 ql_new(i,j,k) > 0.0) then
                Cloud_state%SL_out(i,j,k) =  - Cloud_state%ql_in(i,j,k)
                SQ_out(i,j,k) = SQ_out(i,j,k) + ql_new(i,j,k)
                ST_out(i,j,k) = ST_out(i,j,k) - (hlv*ql_new(i,j,k))/cp_air
                if (diag_id%qdt_cleanup_liquid +    &
                    diag_id%q_cleanup_liquid_col > 0) &
                    diag_4d(i,j,k,diag_pt%qdt_cleanup_liquid) =   &
                                         ql_new(i,j,k)/Constants%dtcloud
                Cloud_state%SN_out(i,j,k) = Cloud_state%SN_out(i,j,k) -   &
                                                             qn_new(i,j,k) 
                if (diag_id%qndt_cleanup + diag_id%qn_cleanup_col > 0) &
                    diag_4d(i,j,k,diag_pt%qndt_cleanup) = - qn_new(i,j,k)/ &
                                                          Constants%dtcloud
              endif
            end do
          end do
        end do

        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              if (abs(qi_new(i,j,k)) .le. Nml%qmin  .and.  &
                  Atmos_state%qv_in(i,j,k) + SQ_out(i,j,k) +   &
                                                qi_new(i,j,k) > 0.0) then
                Cloud_state%SI_out(i,j,k) =  - Cloud_state%qi_in(i,j,k)
                SQ_out(i,j,k) = SQ_out(i,j,k) + qi_new(i,j,k)
                ST_out(i,j,k) = ST_out(i,j,k) - (hls*qi_new(i,j,k))/cp_air
                if (diag_id%qdt_cleanup_ice +    &
                    diag_id%q_cleanup_ice_col > 0) &
                    diag_4d(i,j,k,diag_pt%qdt_cleanup_ice) =   &
                                          qi_new(i,j,k)/Constants%dtcloud
                Cloud_state%SNI_out(i,j,k) = Cloud_state%SNI_out(i,j,k) -  &
                                                             qni_new(i,j,k) 
                if (diag_id%qnidt_cleanup + diag_id%qni_cleanup_col > 0) &
                    diag_4d(i,j,k,diag_pt%qnidt_cleanup) =     &
                                     - qni_new(i,j,k)/Constants%dtcloud
              endif
            end do
          end do
        end do
    
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              if (diag_id%qnidt_cleanup2 + diag_id%qni_cleanup2_col > 0) &
                  diag_4d(i,j,k,diag_pt%qnidt_cleanup2) =   &
                                               Cloud_state%SNi_out(i,j,k)
              Cloud_state%SNi_out(i,j,k) = MAX(Cloud_state%SNi_out(i,j,k), &
                                              - Cloud_state%qni_in(i,j,k))
              if (diag_id%qnidt_cleanup2 + diag_id%qni_cleanup2_col > 0) &
                  diag_4d(i,j,k,diag_pt%qnidt_cleanup2) =    &
                        (diag_4d(i,j,k,diag_pt%qnidt_cleanup2) -   &
                          Cloud_state%SNi_out(i,j,k))*Constants%inv_dtcloud
            end do
          end do
        end do

        do k=1,kdim
          do j=1,jdim
            do i=1,idim
              if (diag_id%qndt_cleanup2 + diag_id%qn_cleanup2_col > 0) &
                  diag_4d(i,j,k,diag_pt%qndt_cleanup2) =    &
                                                 Cloud_state%SN_out(i,j,k)
              Cloud_state%SN_out(i,j,k) = MAX(Cloud_state%SN_out(i,j,k),  &
                                              - Cloud_state%qn_in(i,j,k))
              if (diag_id%qndt_cleanup2 + diag_id%qn_cleanup2_col > 0) &
                  diag_4d(i,j,k,diag_pt%qndt_cleanup2) =    &
                         (diag_4d(i,j,k,diag_pt%qndt_cleanup2) -   &
                            Cloud_state%SN_out(i,j,k))*Constants%inv_dtcloud
            end do
          end do
        end do

        if (diag_id%qadt_destr + diag_id%qa_destr_col > 0)    &
            diag_4d(:,:,:,diag_pt%qadt_destr) =    &
                   diag_4d(:,:,:,diag_pt%qadt_destr) +    &
                                 Cloud_state%SA_out*Constants%inv_dtcloud

        qa_new = Cloud_state%qa_in + Cloud_state%SA_out
        where ( abs(qa_new) .le. Nml%qmin )
          Cloud_state%SA_out  = -Cloud_state%qa_in
        endwhere

        if (diag_id%qadt_destr + diag_id%qa_destr_col > 0)    &
            diag_4d(:,:,:,diag_pt%qadt_destr) =    &
                   diag_4d(:,:,:,diag_pt%qadt_destr) -     &
                                Cloud_state%SA_out*Constants%inv_dtcloud

        if (diag_id%qadt_limits + diag_id%qa_limits_col > 0)    &
            diag_4d(:,:,:,diag_pt%qadt_limits) = Cloud_state%SA_out(:,:,:)

        Cloud_state%SA_out = MAX(Cloud_state%SA_out,-Cloud_state%qa_in)
        Cloud_state%SA_out = MIN(Cloud_state%SA_out,   &
                                1. - Atmos_state%ahuco - Cloud_state%qa_in)
      endif ! rotstayn
!------------------------------------------------------------------------


end subroutine microphysics


!#######################################################################

subroutine microphysics_end (Nml)

type(strat_nml_type), intent(in) :: Nml

      if (.not. module_is_initialized) return

      call rotstayn_klein_microp_end 
      call morrison_gettelman_microp_end (Nml%do_pdf_clouds)

      module_is_initialized = .false.


end subroutine microphysics_end

!#######################################################################



end module microphysics_mod
