#include <fms_platform.h>
MODULE UW_CONV_MOD

  use           mpp_mod, only : mpp_pe, mpp_root_pe, stdlog
  use      Constants_Mod, ONLY: tfreeze,HLv,HLf,HLs,CP_AIR,GRAV,Kappa,rdgas,rvgas
  use   Diag_Manager_Mod, ONLY: register_diag_field, send_data
  use   Time_Manager_Mod, ONLY: time_type, get_time 
  use           mpp_mod, only : input_nml_file
  use           fms_mod, only : write_version_number, open_namelist_file, check_nml_error,&
                                FILE_EXIST, ERROR_MESG,  &
                                lowercase, &
                                CLOSE_FILE, FATAL, NOTE
  use  field_manager_mod, only: MODEL_ATMOS
  use  tracer_manager_mod, only: get_tracer_names, query_method, &
                                 get_tracer_index, NO_TRACER
  use  sat_vapor_pres_mod,only : sat_vapor_pres_init
  use atmos_tracer_utilities_mod, only : get_wetdep_param

  use  aerosol_types_mod, only : aerosol_type
  
  use  aer_ccn_act_mod, only :   aer_ccn_act_init
  use  conv_utilities_mod,only :   uw_params_init
  use  conv_utilities_k_mod,only : sd_init_k, sd_copy_k, sd_end_k,  &
                                   ac_init_k, ac_clear_k, ac_end_k, &
                                   pack_sd_k, adi_cloud_k, extend_sd_k,&
                                   exn_init_k, exn_end_k, findt_init_k,&
                                   findt_end_k, &
                                   check_tracer_realizability, &
                                   qt_parcel_k, qt_parcel_deep_k, &
                                   adicloud, sounding, uw_params

  use  conv_plumes_k_mod,only    : cp_init_k, cp_end_k, cp_clear_k, &
                                   ct_init_k, ct_end_k, ct_clear_k, &
                                   cumulus_tend_k, cumulus_plume_k, &
                                   cplume, ctend, cpnlist, cwetdep_type

  use  conv_closures_mod,only    : cclosure_bretherton,   &
                                   cclosure_relaxcbmf, &
                                   cclosure_relaxwfn,  &
                                   cclosure_implicit, cclosure

  use  deep_conv_mod,only        : deepc, cpn_copy, dpconv0, dpconv1, dpconv2, conv_forced
  use random_numbers_mod,    only: randomNumberStream, getRandomNumbers, &
                                   initializeRandomNumberStream, constructSeed
 
!---------------------------------------------------------------------
  implicit none
  private
!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

!---------------------------------------------------------------------
!-------  interfaces --------

  public  :: uw_conv, uw_conv_init, uw_conv_end

  real, parameter :: aday = 1.
  real, parameter :: mv = -999.
  logical         :: module_is_initialized = .false.

  character(len=7) :: mod_name = 'uw_conv'

  !namelist parameters for UW convection scheme
  integer :: iclosure = 0      ! 0: Bretherton UWShCu orginal / -CIN/TKE based
                               ! 1: Emanuel-Rayment: quasiequilibrium PBL
  real    :: rkm_sh1  = 10.0  
  real    :: rkm_sh   = 3.0    ! fractional lateral mixing rate for shallow
  real    :: cldhgt_max   = 50.e3
  real    :: landfact_m   = 0.5
  integer :: idpchoice = 0  
  logical :: do_deep = .false.
  logical :: do_relaxcape = .false.
  logical :: do_relaxwfn  = .false.
  logical :: do_coldT = .true.
  logical :: do_lands = .false.
  logical :: do_lclht = .false.
  logical :: do_uwcmt = .false.   
  logical :: do_fast  = .false.
  logical :: do_ice   = .true.
  logical :: do_ppen  = .true.
  logical :: do_forcedlifting = .false.
  real    :: atopevap = 0.
  logical :: apply_tendency = .true.
  logical :: prevent_unreasonable = .true.
  real    :: aerol = 1.e-12
  real    :: tkemin   = 1.e-6
  real    :: wmin_ratio = 0.05
  logical :: use_online_aerosol = .false.
  logical :: use_sub_seasalt = .true.
  logical :: do_auto_aero = .false.
  logical :: do_rescale   = .false.
  logical :: do_rescale_t = .false.
  logical :: do_debug     = .false.
!miz
  logical :: do_imposing_forcing = .false.
  real    :: tdt_rate = 0.0             
  real    :: qdt_rate = 0.0
  real    :: pres_min = 0.0
  real    :: pres_max = 0.0
  integer :: klevel = 10
  logical :: use_klevel   = .true.
  logical :: do_imposing_rad_cooling = .false.
  real    :: cooling_rate = -1.5 !K/day
  real    :: t_thresh = 207.5    !K
  real    :: t_strato = 200.0    !K
  real    :: tau_rad  = 5.0      !day
!miz
  integer :: cush_choice  = 0
  real    :: pcp_min      = 3e-5
  real    :: pcp_max      = 1.5e-3
  real    :: rh0          = 0.8
  real    :: cush_ref     = 0.
  real    :: plev_cin     = 60000.
  real    :: pblht0 = 500.
  real    :: lofactor0 = 1.
  integer :: lochoice  = 0
  real    :: wrel_min = 1.
  real    :: om_to_oc = 1.67
  real    :: sea_salt_scale = 0.1
  logical :: do_stime  = .false.
  logical :: do_dtime  = .false.
  logical :: do_qctflx_zero = .false.
  logical :: do_hlflx_zero  = .true.
  logical :: do_subcloud_flx = .false.
  logical :: do_detran_zero = .false.
  logical :: do_prog_tke  = .false.
  logical :: do_prog_gust = .false.
  logical :: do_gust_qt = .false.
  logical :: use_new_let = .false.
  logical :: use_lcl_only =.false.
  logical :: do_new_pevap =.false.
  logical :: stop_at_let  =.false.
  logical :: zero_out_conv_area = .false.
  integer :: src_choice = 0
  integer :: gqt_choice = 0
  real    :: plev_for   = 50000.
  real    :: tke0 = 0.1
  real    :: gama = 0.
  real    :: hgt0 = 500.
  real    :: duration = 10800
  real    :: pblrat0  = 2.0
  real    :: bfact  = 0.05
  real    :: tau_tke  = 7200
  real    :: tau_gust = 7200
  real    :: gfact  = 1.0
  real    :: gfact3 = 0.01
  real    :: gfact4 = 1
  real    :: cgust0 = 1.
  real    :: cgust_max = 10.
  real    :: sigma0 = 0.5
  real    :: stime0 = 0.5
  real    :: dtime0 = 0.5

  integer :: tracer_check_type = -999 !legacy
  !< select realizability checks to be applied to tracers
  !! -999 (default): apply min/max checks (with scaling of tendencies), no filling
  !!              1: omit min/max checks, apply filling (using sjl_fillz), no scaling
  !!              2: omit min/max checks, no filling (apply scaling to avoid negatives)

  logical :: use_turb_tke = .false.  !h1g, 2015-08-11

  NAMELIST / uw_conv_nml / iclosure, rkm_sh1, rkm_sh, cldhgt_max, plev_cin, &
       do_deep, idpchoice, do_relaxcape, do_relaxwfn, do_coldT, do_lands, do_uwcmt,       &
       do_fast, do_ice, do_ppen, do_forcedlifting, do_lclht, do_gust_qt, use_new_let, do_hlflx_zero, &
       atopevap, apply_tendency, prevent_unreasonable, aerol, tkemin, do_prog_tke, tau_tke, pblrat0, &
       wmin_ratio, use_online_aerosol, use_sub_seasalt, landfact_m, pblht0, tke0, lofactor0, lochoice, &
       do_auto_aero, do_rescale, do_rescale_t, wrel_min, om_to_oc, sea_salt_scale, bfact, gfact, gfact3, gfact4, &
       do_debug, cush_choice, pcp_min, pcp_max, cush_ref, do_prog_gust, tau_gust, cgust0, cgust_max, sigma0,&
       rh0, do_qctflx_zero, do_detran_zero, gama, hgt0, duration, do_stime, do_dtime, stime0, dtime0, &
       do_imposing_forcing, tdt_rate, qdt_rate, pres_min, pres_max, klevel, use_klevel, do_subcloud_flx,&
       do_imposing_rad_cooling, cooling_rate, t_thresh, t_strato, tau_rad, src_choice, gqt_choice,&
       zero_out_conv_area, tracer_check_type, use_turb_tke, use_lcl_only, do_new_pevap, plev_for, stop_at_let

  !namelist parameters for UW convective plume
  real    :: rle      = 0.10   ! for critical stopping distance for entrainment
  real    :: rpen     = 5.0    ! for entrainment efficiency
  real    :: rmaxfrac = 0.15   ! maximum allowable updraft fraction
  real    :: wmin     = 0.5    ! minimum vertical velocity for computing updraft fraction
  real    :: wmax     = 50     ! maximum allowable vertical velocity
  real    :: rbuoy    = 1.0    ! for nonhydrostatic pressure effects on updraft
  real    :: rdrag    = 1.0 
  real    :: frac_drs = 0.0    ! 
  real    :: bigc     = 0.7    ! for momentum transfer
  real    :: auto_th0 = 0.5e-3 ! threshold for precipitation
  real    :: auto_rate= 1.e-3
  real    :: tcrit    = -60.0  ! critical temperature 
  real    :: deltaqc0 = 0.5e-3 
  logical :: do_pdfpcp= .false.
  logical :: do_pmadjt= .false.
  logical :: do_emmax = .false.
  logical :: do_pnqv  = .false.
  logical :: do_tten_max = .false.
  real    :: rad_crit = 14.0   ! critical droplet radius
  real    :: emfrac_max = 1.0
  integer :: mixing_assumption = 0
  integer :: mp_choice = 1
  real    :: Nl_land   = 300.e6
  real    :: Nl_ocean  = 100.e6
  real    :: qi_thresh = 1.e-4
  real    :: r_thresh  = 12.e-6
  logical :: do_pevap = .false.
  real    :: cfrac     = 0.1
  real    :: hcevap    = 0.8
  real    :: pblfac    = 0.0
  real    :: ffldep    = 0.0
  logical :: do_weffect = .false.
  logical :: do_limit_wmax =.false.
  real    :: weffect    = 0.5
  real    :: peff_l     = 1.0
  real    :: peff_i     = 1.0
  real    :: t00        = 295
  real    :: tten_max   = 1000.

  NAMELIST / uw_plume_nml / rle, rpen, rmaxfrac, wmin, wmax, rbuoy, rdrag, frac_drs, bigc, ffldep, do_limit_wmax,&
       auto_th0, auto_rate, tcrit, deltaqc0, do_pdfpcp, do_pmadjt, do_emmax, do_pnqv, do_tten_max, rad_crit, emfrac_max, &
       mixing_assumption, mp_choice, Nl_land, Nl_ocean, qi_thresh, r_thresh, do_pevap, cfrac, hcevap, pblfac,&
       do_weffect, weffect, peff_l, peff_i, t00, tten_max
  !namelist parameters for UW convective closure
  integer :: igauss   = 1      ! options for cloudbase massflux closure
                               ! 1: cin/gaussian closure, using TKE to compute CIN.
                               ! 2: cin/gaussian closure, using W* to compute CIN.
                               ! 0: cin/tke mapse-style closure; 
  real    :: rkfre    = 0.05   ! vertical velocity variance as fraction of tke
  real    :: tau_sh   = 7200.  ! 
  real    :: wcrit_min= 0.
  real    :: mass_fact= 0.25
  logical :: do_old_cbmfmax = .true.

  NAMELIST / uw_closure_nml / igauss, rkfre, tau_sh, wcrit_min, mass_fact, do_old_cbmfmax


!========Option for deep convection=======================================
  real    :: cbmf0         = 0.0001
  real    :: rkm_dp1       = 10.
  real    :: rkm_dp2       = 1.
  real    :: cbmf_dp_frac1 = 0.
  real    :: cbmf_dp_frac2 = 1.
  real    :: crh_th_ocean  = 0.5
  real    :: crh_th_land   = 0.5
  real    :: cape_th       = 10.
  real    :: cin_th        = 5.
  real    :: cwfn_th       = 0.
  real    :: tau_dp        = 7200.
  real    :: rpen_d        = 5.0
  integer :: mixing_assumption_d = 0
  integer :: norder      = 1
  logical :: do_ppen_d   = .true.
  logical :: do_pevap_d  = .true.
  real    :: cfrac_d     = 0.05
  real    :: hcevap_d    = 0.8
  real    :: pblfac_d    = 0.0
  real    :: ffldep_d    = 0.0
  real    :: dcapedm_th  = 0
  real    :: dcwfndm_th  = 0
  real    :: frac_limit_d = 0.25
  real    :: lofactor_d   = 1.0
  real    :: auto_th0_d   = 1.0e-3
  real    :: tcrit_d      = -120
  real    :: peff_l_d     = 3.8e-5
  real    :: peff_i_d     = 3.8e-5
  integer :: src_choice_d = 0
  logical :: do_forcedlifting_d = .false.
  logical :: do_lod_rkm   = .false.
  logical :: do_lod_cfrac = .false.
  logical :: do_lod_tcrit = .false.
  logical :: do_lod_cape  = .false.
  logical :: do_lod_tau   = .false.
  logical :: do_lod_cush  = .false.
  logical :: do_stochastic_rkm = .false.
  logical :: do_cgust_dp  = .false.
  real    :: gustmax      = 3.  ! maximum gustiness wind (m/s)
  real    :: cpool_gust   = 10. ! constant for cool pool effect W/m2, default=10 W/m2
  logical :: do_forced_conv = .false.
  real    :: frac_rkm_pert = 1.
  real    :: tau_dp_fact = 10
  real    :: cin_fact    = 1
  real    :: wcrit_min_gust = 0.2
  integer :: cgust_choice = 0
  NAMELIST / deep_conv_nml / cbmf0, rkm_dp1, rkm_dp2, cbmf_dp_frac1, cbmf_dp_frac2, do_forced_conv, &
                 crh_th_ocean, crh_th_land, do_forcedlifting_d, frac_limit_d, wcrit_min_gust, cin_fact,&
                 cape_th, cin_th, cwfn_th, tau_dp, rpen_d, mixing_assumption_d, norder, dcwfndm_th, &
                 do_ppen_d, do_pevap_d, cfrac_d, hcevap_d, pblfac_d, ffldep_d, lofactor_d, dcapedm_th, &
                 auto_th0_d, tcrit_d, do_lod_rkm, do_lod_cfrac, do_lod_tcrit, do_lod_cape, &
		 peff_l_d, peff_i_d, do_lod_tau, do_lod_cush, cgust_choice, tau_dp_fact, &
                 do_stochastic_rkm, frac_rkm_pert, do_cgust_dp, gustmax, cpool_gust, src_choice_d
!========Option for deep convection=======================================

!------------------------------------------------------------------------

  integer :: nqv, nql, nqi, nqa ,nqn
  logical :: do_qn = .false.    ! use droplet number tracer field ?

  integer :: id_tdt_uwc, id_qdt_uwc, id_udt_uwc, id_vdt_uwc, id_prec_uwc, id_snow_uwc, &
       id_cin_uwc, id_cbmf_uwc, id_tke_uwc, id_tkep_uwc, id_plcl_uwc, id_zlcl_uwc, id_zinv_uwc,  &
       id_cush_uwc, id_pct_uwc, id_pcb_uwc, id_plfc_uwc, id_enth_uwc,  &
       id_qldt_uwc, id_qidt_uwc, id_qadt_uwc, id_qndt_uwc, id_cmf_uwc, id_cmf_uws, id_wu_uwc,   &
       id_fer_uwc,  id_fdr_uwc, id_fdrs_uwc, id_cqa_uwc, id_cql_uwc,   &
       id_cqi_uwc,  id_cqn_uwc, id_hlflx_uwc, id_qtflx_uwc, id_nqtflx_uwc, &
       id_cape_uwc, id_dcin_uwc, id_dcape_uwc, id_crh_uwc, id_pblht_uwc, &
       id_ocode_uwc, id_plnb_uwc, id_wrel_uwc, id_ufrc_uwc, id_qtmp_uwc,id_gust_uwc, &
       id_tdt_pevap_uwc, id_qdt_pevap_uwc, id_xpsrc_uwc, id_xhlsrc_uwc, id_xqtsrc_uwc,&
       id_qldet_uwc, id_qidet_uwc, id_qadet_uwc, id_qtdt_uwc, id_dting_uwc, &
       id_cfq_uwc, id_fdp_uwc, id_hmo_uwc, id_hms_uwc, id_abu_uwc, id_peo_uwc, &
       id_tten_rad_uwc, id_tdt_forc_uwc, id_qdt_forc_uwc, id_tdt_diss_uwc, &
       id_hm_vadv_uwc, id_pflx_uwc, id_lhflx_uwc, id_shflx_uwc, &
       id_tdt_rad_uwc, id_tdt_dyn_uwc, id_tdt_dif_uwc, id_qdt_dyn_uwc, id_qdt_dif_uwc, &
       id_tdt_rad_int, id_tdt_dyn_int, id_tdt_dif_int, id_qdt_dyn_int, id_qdt_dif_int, &
       id_tdt_rad_pbl, id_tdt_dyn_pbl, id_tdt_dif_pbl, id_qdt_dyn_pbl, id_qdt_dif_pbl, &
       id_tdt_rad_fre, id_tdt_dyn_fre, id_tdt_dif_fre, id_qdt_dyn_fre, id_qdt_dif_fre, &
       id_tdt_tot_pbl, id_tdt_tot_fre, id_cpool_uwc, id_bflux_uwc, &
       id_dgz_dyn_uwc, id_ddp_dyn_uwc, id_dgz_dyn_int, id_ddp_dyn_int, &
       id_hmint_uwc, id_hm_vadv0_uwc, id_hm_hadv0_uwc, id_hm_tot0_uwc, id_hm_total_uwc,&
       id_qtflx_up_uwc, id_qtflx_dn_uwc, id_omega_up_uwc, id_omega_dn_uwc, &
       id_omgmc_up_uwc, id_rkm_uwc, id_stime_uwc, id_scale_uwc, id_scaletr_uwc


  integer, allocatable :: id_tracerdt_uwc(:), id_tracerdt_uwc_col(:), &
                          id_tracerdtwet_uwc(:), id_tracerdtwet_uwc_col(:), &
                          id_tracerdt_uwc_nc(:), id_tracerdt_uwc_col_nc(:), id_rn(:)
  integer, allocatable :: id_trevp_uwc(:), id_trevp_uwd(:)

!========Option for deep convection=======================================
  integer :: id_tdt_uwd, id_qdt_uwd, id_qtdt_uwd, id_prec_uwd, id_snow_uwd,   &
       id_cbmf_uwd, id_enth_uwd, id_qldt_uwd, id_qidt_uwd,&
       id_qndt_uwd, id_qadt_uwd, id_cmf_uwd, id_wu_uwd, id_fer_uwd,    &
       id_fdr_uwd, id_fdrs_uwd, id_cqa_uwd, id_cql_uwd, id_cqi_uwd,    &
       id_cqn_uwd, id_hlflx_uwd, id_qtflx_uwd, id_nqtflx_uwd, id_dcin_uwd, &
       id_dcapedm_uwd, id_dcwfndm_uwd, id_ocode_uwd, id_cush_uwd,      &
       id_tdt_pevap_uwd, id_qdt_pevap_uwd, id_rkm_uwd, id_cbu_uwd,     &
       id_rand_uwd, id_taudp_uwd, id_pwfn_uwd, id_cwfn_uwd, id_dcwfndt_dpc, &
       id_dcwfndt_fre, id_dcwfndt_pbl, id_cwfn3d_uwd, id_cape3d_uwd, id_dtime_uwd
!========Option for deep convection=======================================

  type(cwetdep_type), dimension(:), allocatable :: wetdep
  type(uw_params),  save  :: Uw_p
  character(len=32), dimension(:), allocatable   :: tracername 
  character(len=32), dimension(:), allocatable   :: tracer_units 

contains

!#####################################################################
!#####################################################################

  SUBROUTINE UW_CONV_INIT(do_strat, axes, Time, kd, tracers_in_uw)
    logical,         intent(in) :: do_strat
    integer,         intent(in) :: axes(4), kd
    type(time_type), intent(in) :: Time
    logical,         intent(in) :: tracers_in_uw(:)
    
!---------------------------------------------------------------------
!  intent(in) variables:
!
!      tracers_in_uw 
!                   logical array indicating which of the activated 
!                   tracers are to be transported by UW convection
!
!-------------------------------------------------------------------

    integer   :: unit, io
    
    integer   :: ntracers, n, nn, ierr, logunit
    logical   :: flag
    character(len=200) :: text_in_scheme, control
     real :: frac_junk, frac_junk2
 
    ntracers = count(tracers_in_uw)

    call uw_params_init   (Uw_p)

!   Initialize lookup tables needed for findt and exn
!   sat_vapor_pres needs to be initialized if not already done
    call sat_vapor_pres_init     
    call exn_init_k (Uw_p)
    call findt_init_k (Uw_p)

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, nml=uw_closure_nml, iostat=io)
      ierr = check_nml_error(io,'uw_closure_nml')
      read (input_nml_file, nml=uw_conv_nml, iostat=io)
      ierr = check_nml_error(io,'uw_conv_nml')
      read (input_nml_file, nml=uw_plume_nml, iostat=io)
      ierr = check_nml_error(io,'uw_plume_nml')
      read (input_nml_file, nml=deep_conv_nml, iostat=io)
      ierr = check_nml_error(io,'deep_conv_nml')
#else   
    if( FILE_EXIST( 'input.nml' ) ) then
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_closure_nml, iostat = io, end = 10 )
          ierr = check_nml_error(io,'uw_closure_nml')
       end do
10     call close_file ( unit )

       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_conv_nml, iostat = io, end = 20 )
          ierr = check_nml_error(io,'uw_conv_nml')
       end do
20     call close_file ( unit )
       
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = uw_plume_nml, iostat = io, end = 30 )
          ierr = check_nml_error(io,'uw_plume_nml')
       end do
30     call close_file ( unit )

!========Option for deep convection=======================================
       unit = OPEN_NAMELIST_FILE ()
       io = 1
       do while ( io .ne. 0 )
          READ( unit, nml = deep_conv_nml, iostat = io, end = 40 )
          ierr = check_nml_error(io,'deep_conv_nml')
       end do
40     call close_file ( unit )
!========Option for deep convection=======================================
    end if
#endif
    call write_version_number (version, tagname)
    logunit = stdlog()
    WRITE( logunit, nml = uw_closure_nml )
    WRITE( logunit, nml = uw_conv_nml )
    WRITE( logunit, nml = uw_plume_nml )
    WRITE( logunit, nml = deep_conv_nml )

    if ( use_online_aerosol ) call aer_ccn_act_init

    nqv = get_tracer_index ( MODEL_ATMOS, 'sphum' )
    nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
    nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
    nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
    nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
    if (nqn /= NO_TRACER) do_qn = .true.
    if (ntracers > 0) then
      allocate ( tracername   (ntracers) )
      allocate ( tracer_units (ntracers) )
      allocate ( wetdep       (ntracers) )
      nn = 1
      do n=1,size(tracers_in_uw(:))
         if (tracers_in_uw(n)) then
             call get_tracer_names (MODEL_ATMOS, n,  &
                                    name = tracername(nn), &
                                    units = tracer_units(nn))
             flag = query_method( 'wet_deposition', MODEL_ATMOS, n, &
                                  text_in_scheme, control )
             call get_wetdep_param( text_in_scheme, control, &
                                    wetdep(nn)%scheme, &
                                    wetdep(nn)%Henry_constant, &
                                    wetdep(nn)%Henry_variable, &
                                    frac_junk, frac_junk2, &
                                    wetdep(nn)%alpha_r, &
                                    wetdep(nn)%alpha_s, &
                                    wetdep(nn)%Lwetdep, &
                                    wetdep(nn)%Lgas, &
                                    wetdep(nn)%Laerosol, &
                                    wetdep(nn)%Lice, &
                                    frac_in_cloud_uw=wetdep(nn)%frac_in_cloud )
             wetdep(nn)%scheme = lowercase( wetdep(nn)%scheme )
             nn = nn + 1
          endif
       end do
    endif

!---> h1g, 2015-08-11
    if ( do_prog_tke .and. use_turb_tke )  then
        call error_mesg ('uw_conv_mod',  &
                '  do_prog_tke and use_turb_tke cannot be true at the same time', FATAL)
    endif
!<--- h1g, 2015-08-11

    id_xpsrc_uwc  = register_diag_field (mod_name,'xpsrc_uwc', axes(1:2), Time, &
         'xpsrc', 'hPa' )
    id_xhlsrc_uwc = register_diag_field (mod_name,'xhlsrc_uwc', axes(1:2), Time, &
         'xhlsrc', 'J/kg' )
    id_xqtsrc_uwc = register_diag_field (mod_name,'xqtsrc_uwc', axes(1:2), Time, &
         'xqtsrc', 'kg/kg' )

    id_tdt_pevap_uwc = register_diag_field ( mod_name, 'tdt_pevap_uwc', axes(1:3), Time, &
         'Temperature tendency due to pevap from uw_conv', 'K/s', missing_value=mv )
    id_qdt_pevap_uwc = register_diag_field ( mod_name, 'qdt_pevap_uwc', axes(1:3), Time, &
         'Spec. humidity tendency due to pevap from uw_conv', 'kg/kg/s', missing_value=mv)

    id_tdt_uwc = register_diag_field ( mod_name, 'tdt_uwc', axes(1:3), Time, &
         'Temperature tendency from uw_conv', 'K/s', missing_value=mv )
    id_qdt_uwc = register_diag_field ( mod_name, 'qdt_uwc', axes(1:3), Time, &
         'Spec. humidity tendency from uw_conv', 'kg/kg/s', missing_value=mv)
    id_udt_uwc = register_diag_field ( mod_name, 'udt_uwc', axes(1:3), Time, &
         'U tendency from uw_conv', 'm/s2', missing_value=mv )
    id_vdt_uwc = register_diag_field ( mod_name, 'vdt_uwc', axes(1:3), Time, &
         'V tendency from uw_conv', 'm/s2', missing_value=mv)
    id_cmf_uwc = register_diag_field ( mod_name, 'cmf_uwc', axes(1:3), Time, &
         'Cloud vert. mass flux from uw_conv', 'kg/m2/s', missing_value=mv)
    id_cmf_uws = register_diag_field ( mod_name, 'cmf_uws', axes(1:3), Time, &
         'Cloud vert. mass flux from shallow plume uw_conv', 'kg/m2/s', missing_value=mv)
    id_cfq_uwc = register_diag_field ( mod_name, 'cfq_uwc', axes(1:3), Time,   &
         'Convective frequency', 'none', missing_value=mv)
    id_peo_uwc = register_diag_field ( mod_name, 'peo_uwc', axes(1:3), Time,   &
         'Convective precipitation efficiency', 'none', missing_value=mv)
    id_hmo_uwc = register_diag_field ( mod_name, 'hmo_uwc', axes(1:3), Time,   &
         'moist static energy', 'J/kg', missing_value=mv)
    id_hms_uwc = register_diag_field ( mod_name, 'hms_uwc', axes(1:3), Time,   &
         'moist static energy', 'J/kg', missing_value=mv)
    id_abu_uwc = register_diag_field ( mod_name, 'abu_uwc', axes(1:3), Time,   &
         'adiabatic buoyancy', 'K', missing_value=mv)
    id_wu_uwc = register_diag_field ( mod_name, 'wu_uwc', axes(1:3), Time,   &
         'Updraft vert. velocity from uw_conv', 'm/s', missing_value=mv)
    id_fer_uwc = register_diag_field ( mod_name, 'fer_uwc', axes(1:3), Time, &
         'Fractional entrainment rate from uw_conv', '1/Pa', missing_value=mv)
    id_fdr_uwc = register_diag_field ( mod_name, 'fdr_uwc', axes(1:3), Time, &
         'Fractional detrainment rate from uw_conv', '1/Pa', missing_value=mv)
    id_fdrs_uwc = register_diag_field (mod_name,'fdrs_uwc', axes(1:3), Time, &
         'Detrainment rate for sat. air from uw_conv', '1/Pa', missing_value=mv)
    id_cqa_uwc = register_diag_field ( mod_name, 'cqa_uwc', axes(1:3), Time, &
         'Updraft fraction from uw_conv', 'none', missing_value=mv)
    id_cql_uwc = register_diag_field ( mod_name, 'cql_uwc', axes(1:3), Time, &
         'Updraft liquid from uw_conv', 'kg/kg', missing_value=mv)
    id_cqi_uwc = register_diag_field ( mod_name, 'cqi_uwc', axes(1:3), Time, &
         'Updraft ice from uw_conv', 'kg/kg', missing_value=mv)
    id_cqn_uwc = register_diag_field ( mod_name, 'cqn_uwc', axes(1:3), Time, &
         'Updraft liquid drop from uw_conv', '/kg', missing_value=mv)
    id_hlflx_uwc=register_diag_field (mod_name,'hlflx_uwc',axes(1:3),Time, &
         'Liq.wat.pot.temp. flux from uw_conv', 'W/m2', missing_value=mv)
    id_qtflx_uwc = register_diag_field (mod_name,'qtflx_uwc',axes(1:3),Time, &
         'Total water flux from uw_conv', 'kg/m2/s', missing_value=mv)
    id_nqtflx_uwc = register_diag_field (mod_name,'nqtflx_uwc',axes(1:3),Time, &
         'net total water flux from uw_conv', 'kg/m2/s', missing_value=mv)
    id_qtflx_up_uwc = register_diag_field (mod_name,'qtflx_up_uwc',axes(1:3),Time, &
         'Total water flux from resolved upward flow', 'kg/m2/s', missing_value=mv)
    id_qtflx_dn_uwc = register_diag_field (mod_name,'qtflx_dn_uwc',axes(1:3),Time, &
         'Total water flux from resolved downward flow', 'kg/m2/s', missing_value=mv)
    id_omgmc_up_uwc = register_diag_field (mod_name,'omgmc_up_uwc',axes(1:3),Time, &
         'Total upward mass flux', 'kg/m2/s', missing_value=mv)
    id_omega_up_uwc = register_diag_field (mod_name,'omega_up_uwc',axes(1:3),Time, &
         'Total mass flux from resolved upward flow', 'kg/m2/s', missing_value=mv)
    id_omega_dn_uwc = register_diag_field (mod_name,'omega_dn_uwc',axes(1:3),Time, &
         'Total mass flux from resolved downward flow', 'kg/m2/s', missing_value=mv)
    id_hm_vadv_uwc = register_diag_field (mod_name,'hm_vadv_uwc',axes(1:3),Time, &
         'Vertical advection of MSE', 'W/m2', missing_value=mv)
    id_pflx_uwc = register_diag_field (mod_name,'pflx_uwc',axes(1:3),Time, &
         '3D precipitation flux', 'kg/m2/s', missing_value=mv)

    id_lhflx_uwc = register_diag_field (mod_name,'lhflx_uwc', axes(1:2), Time, &
         'surface latent heat flux from uw_conv', 'W/m2',                      &
         interp_method = "conserve_order1" )
    id_shflx_uwc = register_diag_field (mod_name,'shflx_uwc', axes(1:2), Time, &
         'surface buoyancy flux from uw_conv', 'W/m2',                       &
         interp_method = "conserve_order1" )
    id_hmint_uwc = register_diag_field (mod_name,'hmint_uwc', axes(1:2), Time, &
         'vertically averaged MSE', 'J/m2',                       &
         interp_method = "conserve_order1" )
    id_hm_vadv0_uwc = register_diag_field (mod_name,'hm_vadv0_uwc', axes(1:2), Time, &
         'vertical advection of MSE from uw_conv',   'W/m2',   &
         interp_method = "conserve_order1" )
    id_hm_hadv0_uwc = register_diag_field (mod_name,'hm_hadv0_uwc', axes(1:2), Time, &
         'horizontal advection of MSE from uw_conv', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_hm_tot0_uwc = register_diag_field (mod_name,'hm_tot0_uwc', axes(1:2), Time, &
         'total MSE tendency from uw_conv',   'W/m2',   &
         interp_method = "conserve_order1" )
    id_hm_total_uwc = register_diag_field (mod_name,'hm_total_uwc', axes(1:2), Time, &
         'total MSE tendency computed based on previous time step value', 'W/m2',   &
         interp_method = "conserve_order1" )

    id_tdt_rad_uwc = register_diag_field (mod_name,'tdt_rad_uwc',axes(1:3),Time, &
         'radiative cooling rate', 'K/s', missing_value=mv)
    id_tdt_dyn_uwc = register_diag_field (mod_name,'tdt_dyn_uwc',axes(1:3),Time, &
         'dynamics T tendency', 'K/s', missing_value=mv)
    id_tdt_dif_uwc = register_diag_field (mod_name,'tdt_dif_uwc',axes(1:3),Time, &
         'diffusion T tendency', 'K/s', missing_value=mv)
    id_qdt_dyn_uwc = register_diag_field (mod_name,'qdt_dyn_uwc',axes(1:3),Time, &
         'dynamics qv tendency', 'kg/kg/s', missing_value=mv)
    id_qdt_dif_uwc = register_diag_field (mod_name,'qdt_dif_uwc',axes(1:3),Time, &
         'diffusion qv tendency', 'kg/kg/s', missing_value=mv)
    id_dgz_dyn_uwc = register_diag_field (mod_name,'dgz_dyn_uwc',axes(1:3),Time, &
         'geopotential height tendency due to dynamics', 'm2/s2/s', missing_value=mv)
    id_ddp_dyn_uwc = register_diag_field (mod_name,'ddp_dyn_uwc',axes(1:3),Time, &
         'hm change due to mass change from the dynamics', 'm2/s2/s', missing_value=mv)

    id_tdt_rad_int = register_diag_field (mod_name,'tdt_rad_int', axes(1:2), Time, &
         'vertically integrated radiative T tedency', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_dyn_int = register_diag_field (mod_name,'tdt_dyn_int', axes(1:2), Time, &
         'vertically integrated dynamic T tendency', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_dif_int = register_diag_field (mod_name,'tdt_dif_int', axes(1:2), Time, &
         'vertically integrated diffusion T tendency', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_qdt_dyn_int = register_diag_field (mod_name,'qdt_dyn_int', axes(1:2), Time, &
         'vertically integrated dynamic q tendency', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_qdt_dif_int = register_diag_field (mod_name,'qdt_dif_int', axes(1:2), Time, &
         'vertically integrated diffusion q tendency', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_dgz_dyn_int = register_diag_field (mod_name,'dgz_dyn_int', axes(1:2), Time, &
         'geopotential height tendency due to dynamics', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_ddp_dyn_int = register_diag_field (mod_name,'ddp_dyn_int', axes(1:2), Time, &
         'hm change due to mass change from the dynamics', 'W/m2',   &
         interp_method = "conserve_order1" )

    id_tdt_rad_pbl = register_diag_field (mod_name,'tdt_rad_pbl', axes(1:2), Time, &
         'vertically integrated radiative T tedency with PBL', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_dyn_pbl = register_diag_field (mod_name,'tdt_dyn_pbl', axes(1:2), Time, &
         'vertically integrated dynamic T tendency with PBL', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_dif_pbl = register_diag_field (mod_name,'tdt_dif_pbl', axes(1:2), Time, &
         'vertically integrated diffusion T tendency with PBL', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_qdt_dyn_pbl = register_diag_field (mod_name,'qdt_dyn_pbl', axes(1:2), Time, &
         'vertically integrated dynamic q tendency with PBL', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_qdt_dif_pbl = register_diag_field (mod_name,'qdt_dif_pbl', axes(1:2), Time, &
         'vertically integrated diffusion q tendency with PBL', 'W/m2',   &
         interp_method = "conserve_order1" )

    id_tdt_rad_fre = register_diag_field (mod_name,'tdt_rad_fre', axes(1:2), Time, &
         'vertically integrated radiative T tedency within FRE', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_dyn_fre = register_diag_field (mod_name,'tdt_dyn_fre', axes(1:2), Time, &
         'vertically integrated dynamic T tendency within FRE', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_dif_fre = register_diag_field (mod_name,'tdt_dif_fre', axes(1:2), Time, &
         'vertically integrated diffusion T tendency within FRE', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_qdt_dyn_fre = register_diag_field (mod_name,'qdt_dyn_fre', axes(1:2), Time, &
         'vertically integrated dynamic q tendency within FRE', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_qdt_dif_fre = register_diag_field (mod_name,'qdt_dif_fre', axes(1:2), Time, &
         'vertically integrated diffusion q tendency within FRE', 'W/m2',   &
         interp_method = "conserve_order1" )

    id_tdt_tot_fre = register_diag_field (mod_name,'tdt_tot_fre', axes(1:2), Time, &
         'vertically integrated total T tedency within FRE', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_tdt_tot_pbl = register_diag_field (mod_name,'tdt_tot_pbl', axes(1:2), Time, &
         'vertically integrated total T tendency within FRE', 'W/m2',   &
         interp_method = "conserve_order1" )
    id_cpool_uwc = register_diag_field (mod_name,'cpool_uwc', axes(1:2), Time, &
         'cold pool', 'm2/s3', &
         interp_method = "conserve_order1" )
    id_bflux_uwc = register_diag_field (mod_name,'bflux_uwc', axes(1:2), Time, &
         'surface buoyancy flux', 'm2/s3', &
         interp_method = "conserve_order1" )

    id_prec_uwc = register_diag_field (mod_name,'prec_uwc', axes(1:2), Time, &
         'Precipitation rate from uw_conv', 'kg/m2/sec',                     &
         interp_method = "conserve_order1" )
    id_snow_uwc = register_diag_field (mod_name,'snow_uwc', axes(1:2), Time, &
         'Frozen precip. rate from uw_conv', 'kg/m2/sec',                       &
         interp_method = "conserve_order1" )
    id_cin_uwc = register_diag_field ( mod_name, 'cin_uwc', axes(1:2), Time, &
         'CIN from uw_conv', 'm2/s2' )
    id_cape_uwc= register_diag_field ( mod_name,'cape_uwc', axes(1:2), Time, &
         'CAPE from uw_conv', 'm2/s2' )
    id_gust_uwc= register_diag_field ( mod_name,'gust_uwc', axes(1:2), Time, &
         'gustiness from uw_conv', 'm2/s2' )
    id_crh_uwc= register_diag_field ( mod_name,'crh_uwc', axes(1:2), Time, &
         'Column RH from uw_conv', '%' )
    id_pblht_uwc= register_diag_field ( mod_name,'pblht_uwc', axes(1:2), Time, &
         'PBL height from uw_conv', 'm' )
    id_cbmf_uwc = register_diag_field (mod_name,'cbmf_uwc', axes(1:2), Time, &
         'Cloud-base mass flux from uw_conv', 'kg/m2/s' )
    id_wrel_uwc = register_diag_field (mod_name,'wrel_uwc', axes(1:2), Time, &
         'Release level vertical velocity from uw_conv', 'm/s' )
    id_ufrc_uwc = register_diag_field (mod_name,'ufrc_uwc', axes(1:2), Time, &
         'Release level updraft fraction from uw_conv', 'none' )
    id_tke_uwc = register_diag_field ( mod_name, 'tke_uwc', axes(1:2), Time, &
         'PBL mean TKE from uw_conv', 'm2/s2' )
    id_tkep_uwc = register_diag_field ( mod_name, 'tkep_uwc', axes(1:2), Time, &
         'prognostic estimate of PBL mean TKE from uw_conv', 'm2/s2' )
    id_plcl_uwc = register_diag_field (mod_name,'plcl_uwc', axes(1:2), Time, &
         'LCL pressure from uw_conv', 'hPa' )
    id_zlcl_uwc = register_diag_field (mod_name,'zlcl_uwc', axes(1:2), Time, &
         'LCL depth from uw_conv', 'm' )
    id_plfc_uwc = register_diag_field (mod_name,'plfc_uwc', axes(1:2), Time, &
         'LFC pressure from uw_conv', 'hPa' )
    id_plnb_uwc = register_diag_field (mod_name,'plnb_uwc', axes(1:2), Time, &
         'LNB pressure from uw_conv', 'hPa' )
    id_zinv_uwc = register_diag_field (mod_name,'zinv_uwc', axes(1:2), Time, &
         'Inversion pressure from uw_conv', 'm' )
    id_pct_uwc = register_diag_field ( mod_name, 'pct_uwc', axes(1:2), Time, &
         'Cloud-top pressure from uw_conv', 'hPa' )
    id_pcb_uwc = register_diag_field ( mod_name, 'pcb_uwc', axes(1:2), Time, &
         'Cloud-base pressure from uw_conv', 'hPa' )
    id_cush_uwc = register_diag_field (mod_name,'cush_uwc', axes(1:2), Time, &
         'Convective scale height from uw_conv', 'm' )
    id_dcin_uwc = register_diag_field (mod_name, 'dcin_uwc', axes(1:2), Time, &
         'dCIN/cbmf from uw_conv', 'm2/s2/(kg/m2/s)' )
    id_dcape_uwc= register_diag_field (mod_name, 'dcape_uwc', axes(1:2), Time, &
         'dCAPE/cbmf from uw_conv', 'm2/s2/(kg/m2/s)' )
    id_enth_uwc = register_diag_field (mod_name,'enth_uwc', axes(1:2), Time, &
         'Column-integrated enthalpy tendency from uw_conv', 'W/m2' )
    id_qtmp_uwc = register_diag_field (mod_name,'qtmp_uwc', axes(1:2), Time, &
         'Column-integrated water tendency from uw_conv', 'kg/m2/s' )
    id_dting_uwc = register_diag_field (mod_name,'dting_uwc', axes(1:2), Time, &
         'Column-integrated heating rate from uw_conv', 'W/m2' )
    id_ocode_uwc = register_diag_field (mod_name,'ocode_uwc', axes(1:2), Time, &
         'Out code from uw_conv', 'none' )
    id_fdp_uwc = register_diag_field (mod_name, 'fdp_uwc',   axes(1:2), Time,   &
         'Deep convective frequency', 'none', missing_value=mv)
    id_rkm_uwc = register_diag_field (mod_name, 'rkm_uwc', axes(1:2), Time, &
            'rkm for shallow_conv', 'none' )
    id_stime_uwc = register_diag_field (mod_name, 'stime_uwc', axes(1:2), Time, &
            'stime for shallow_conv', 's' )
    id_scale_uwc = register_diag_field (mod_name, 'scale_uwc', axes(1:2), Time, &
            'scale_uwc in shallow_conv', '' )
    id_scaletr_uwc = register_diag_field (mod_name, 'scaletr_uwc', axes(1:2), Time, &
            'scaletr_uwc in shallow_conv', '' )
    if ( do_strat ) then
       id_qldt_uwc= register_diag_field (mod_name,'qldt_uwc',axes(1:3),Time, &
            'Liquid water tendency from uw_conv', 'kg/kg/s', missing_value=mv)
       id_qidt_uwc= register_diag_field (mod_name,'qidt_uwc',axes(1:3),Time, &
            'Ice water tendency from uw_conv', 'kg/kg/s', missing_value=mv)
       id_qadt_uwc= register_diag_field (mod_name,'qadt_uwc',axes(1:3),Time, &
            'CLD fraction tendency from uw_conv', '1/s', missing_value=mv )
       id_qndt_uwc= register_diag_field (mod_name,'qndt_uwc',axes(1:3),Time, &
            'Cloud droplet number fraction tendency from uw_conv', '#/kg/s', missing_value=mv )
       id_qldet_uwc = register_diag_field (mod_name,'qldet_uwc',axes(1:3),Time, &
            'ql detrainment', 'kg/kg/s', missing_value=mv)
       id_qidet_uwc = register_diag_field (mod_name,'qidet_uwc',axes(1:3),Time, &
            'qi detrainment', 'kg/kg/s', missing_value=mv)
       id_qadet_uwc = register_diag_field (mod_name,'qadet_uwc',axes(1:3),Time, &
            'qa detrainment', '1/s', missing_value=mv)
       id_qtdt_uwc= register_diag_field (mod_name,'qtdt_uwc',axes(1:3),Time, &
            'Total water tendency from uw_conv', 'kg/kg/s', missing_value=mv)
    end if

    if (do_imposing_rad_cooling) then
       id_tten_rad_uwc = register_diag_field ( mod_name, 'tten_rad_uwc', axes(1:3), Time, &
         'Idealized radiative temperature tendency from uw_conv', 'K/s', missing_value=mv )
    end if
    if (do_imposing_forcing) then
       id_tdt_forc_uwc = register_diag_field ( mod_name, 'tdt_forc_uwc', axes(1:3), Time, &
         'Idealized temperature forcing from uw_conv', 'K/s', missing_value=mv )
       id_qdt_forc_uwc = register_diag_field ( mod_name, 'qdt_forc_uwc', axes(1:3), Time, &
         'Idealized humidity forcing from uw_conv', 'K/s', missing_value=mv )
    end if
    if (do_uwcmt) then
        id_tdt_diss_uwc = register_diag_field ( mod_name, 'tdt_diss_uwc', axes(1:3), Time, &
         'Temperature tendency due to dissipation from uw_conv', 'K/s', missing_value=mv )
    end if
!========Option for deep convection=======================================
    if (do_deep) then
       id_tdt_pevap_uwd = register_diag_field ( mod_name, 'tdt_pevap_uwd', axes(1:3), Time, &
            'Temperature tendency due to pevap from deep_conv', 'K/s', missing_value=mv )
       id_qdt_pevap_uwd = register_diag_field ( mod_name, 'qdt_pevap_uwd', axes(1:3), Time, &
            'Spec. humidity tendency due to pevap from deep_conv', 'kg/kg/s', missing_value=mv)

       id_tdt_uwd = register_diag_field ( mod_name, 'tdt_uwd', axes(1:3), Time, &
            'Temperature tendency from deep_conv', 'K/s', missing_value=mv )
       id_qdt_uwd = register_diag_field ( mod_name, 'qdt_uwd', axes(1:3), Time, &
            'Spec. humidity tendency from deep_conv', 'kg/kg/s', missing_value=mv)
       id_qtdt_uwd= register_diag_field ( mod_name, 'qtdt_uwd', axes(1:3), Time, &
            'Total water spec. humidity tendency from deep_conv', 'kg/kg/s', missing_value=mv)
       id_cmf_uwd = register_diag_field ( mod_name, 'cmf_uwd', axes(1:3), Time, &
            'Cloud vert. mass flux from deep_conv', 'kg/m2/s', missing_value=mv)
       id_wu_uwd = register_diag_field ( mod_name, 'wu_uwd', axes(1:3), Time,   &
            'Updraft vert. velocity from deep_conv', 'm/s', missing_value=mv)
       id_cbu_uwd= register_diag_field ( mod_name, 'cbu_uwd', axes(1:3), Time,   &
            'deep plume buoyancy', 'K', missing_value=mv)
       id_fer_uwd = register_diag_field ( mod_name, 'fer_uwd', axes(1:3), Time, &
         'Fractional entrainment rate from deep_conv', '1/Pa', missing_value=mv)
       id_fdr_uwd = register_diag_field ( mod_name, 'fdr_uwd', axes(1:3), Time, &
            'Fractional detrainment rate from deep_conv', '1/Pa', missing_value=mv)
       id_fdrs_uwd = register_diag_field (mod_name,'fdrs_uwd', axes(1:3), Time, &
            'Detrainment rate for sat. air from deep_conv', '1/Pa', missing_value=mv)
       id_cqa_uwd = register_diag_field ( mod_name, 'cqa_uwd', axes(1:3), Time, &
            'Updraft fraction from deep_conv', 'none', missing_value=mv)
       id_cql_uwd = register_diag_field ( mod_name, 'cql_uwd', axes(1:3), Time, &
         'Updraft liquid from deep_conv', 'kg/kg', missing_value=mv)
       id_cqi_uwd = register_diag_field ( mod_name, 'cqi_uwd', axes(1:3), Time, &
            'Updraft ice from deep_conv', 'kg/kg', missing_value=mv)
       id_cqn_uwd = register_diag_field ( mod_name, 'cqn_uwd', axes(1:3), Time, &
            'Updraft liquid drop from deep_conv', '/kg', missing_value=mv)
       id_hlflx_uwd=register_diag_field (mod_name,'hlflx_uwd',axes(1:3),Time, &
            'Liq.wat.pot.temp. flux from deep_conv', 'W/m2', missing_value=mv)
       id_qtflx_uwd = register_diag_field (mod_name,'qtflx_uwd',axes(1:3),Time, &
            'Total water flux from deep_conv', 'W/m2', missing_value=mv)
       id_nqtflx_uwd = register_diag_field (mod_name,'nqtflx_uwd',axes(1:3),Time, &
            'net total water flux from deep_conv', 'W/m2', missing_value=mv)
       id_prec_uwd = register_diag_field (mod_name,'prec_uwd', axes(1:2), Time, &
            'Precipitation rate from deep_conv', 'kg/m2/sec' )
       id_snow_uwd = register_diag_field (mod_name,'snow_uwd', axes(1:2), Time, &
            'Frozen precip. rate from deep_conv', 'kg/m2/sec' )
       id_cbmf_uwd = register_diag_field (mod_name,'cbmf_uwd', axes(1:2), Time, &
            'Cloud-base mass flux from deep_conv', 'kg/m2/s' )
       id_cwfn_uwd = register_diag_field (mod_name,'cwfn_uwd', axes(1:2), Time, &
            'Cloud work function from deep_conv', 'kg/m2/s' )
       id_dcapedm_uwd= register_diag_field (mod_name, 'dcapedm_uwd', axes(1:2), Time, &
            'dCAPE/cbmf from deep_conv', 'm2/s2/(kg/m2/s)' )
       id_dcwfndm_uwd= register_diag_field (mod_name, 'dcwfndm_uwd', axes(1:2), Time, &
            'dCWFN/cbmf from deep_conv', '(m2/s2)/(kg/m2/s)' )
       id_taudp_uwd= register_diag_field (mod_name, 'taudp_uwd', axes(1:2), Time, &
            'taudp from deep_conv', 's' )
       id_cush_uwd = register_diag_field (mod_name, 'cush_uwd',  axes(1:2), Time, &
            'convective depth from deep_conv', 'm' )
       id_enth_uwd = register_diag_field (mod_name,'enth_uwd', axes(1:2), Time, &
            'Column-integrated enthalpy tendency from deep_conv', 'K/s' )
       id_ocode_uwd = register_diag_field (mod_name,'ocode_uwd', axes(1:2), Time, &
            'Out code from deep_conv', 'none' )
       id_rkm_uwd = register_diag_field (mod_name,'rkm_uwd', axes(1:2), Time, &
            'rkm for deep_conv', 'none' )
       id_rand_uwd = register_diag_field (mod_name,'rand_uwd', axes(1:2), Time, &
         'rand_uwd', 'none' )
       id_dtime_uwd= register_diag_field (mod_name,'dtime_uwd',axes(1:2), Time, &
            'dtime for deep_conv', 's' )
       if ( do_strat ) then
          id_qldt_uwd= register_diag_field (mod_name,'qldt_uwd',axes(1:3),Time, &
               'Liquid water tendency from deep_conv', 'kg/kg/s', missing_value=mv)
          id_qidt_uwd= register_diag_field (mod_name,'qidt_uwd',axes(1:3),Time, &
               'Ice water tendency from deep_conv', 'kg/kg/s', missing_value=mv)
          id_qadt_uwd= register_diag_field (mod_name,'qadt_uwd',axes(1:3),Time, &
               'CLD fraction tendency from deep_conv', '1/s', missing_value=mv )
       end if
    end if
!========Option for deep convection=======================================


    if ( ntracers>0 ) then
      allocate(id_tracerdt_uwc(ntracers), id_tracerdt_uwc_col(ntracers) )
      allocate(id_tracerdt_uwc_nc(ntracers), id_tracerdt_uwc_col_nc(ntracers)) 
      allocate(id_rn(ntracers)) 
      allocate(id_tracerdtwet_uwc(ntracers), id_tracerdtwet_uwc_col(ntracers))
       allocate(id_trevp_uwc(ntracers))
       allocate(id_trevp_uwd(ntracers))
      do nn = 1,ntracers
         id_rn(nn) = &
            register_diag_field (mod_name, trim(tracername(nn))//'_rscale', &
                                    axes(1:3), Time, &
                                  trim(tracername(nn)) //' correction', &
                                  'none', missing_value=mv)

         id_tracerdt_uwc(nn) = &
            register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc', &
                                    axes(1:3), Time, &
                                  trim(tracername(nn)) //' tendency from uw_conv', &
                                  trim(tracer_units(nn))//'/s', missing_value=mv)
            id_tracerdt_uwc_col(nn) = &
              register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_col', &
                                     axes(1:2), Time, &
                                   trim(tracername(nn)) //' column tendency from uw_conv', &
                                   trim(tracer_units(nn))//'*(kg/m2)/s', missing_value=mv)
         id_tracerdt_uwc_nc(nn) = &
            register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_nc', &
                                    axes(1:3), Time, &
                                  trim(tracername(nn)) //' tendency from uw_conv before correction', &
                                  trim(tracer_units(nn))//'/s', missing_value=mv)
            id_tracerdt_uwc_col_nc(nn) = &
              register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_col_nc', &
                                     axes(1:2), Time, &
                                   trim(tracername(nn)) //' column tendency from uw_conv before correction', &
                                   trim(tracer_units(nn))//'*(kg/m2)/s', missing_value=mv)
           id_tracerdtwet_uwc(nn) = &
              register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_wet', &
                                    axes(1:3), Time, &
                                   trim(tracername(nn)) //' tendency from uw_conv wetdep', &
                                   trim(tracer_units(nn))//'/s', missing_value=mv)
            id_tracerdtwet_uwc_col(nn) = &
              register_diag_field (mod_name, trim(tracername(nn))//'dt_uwc_wet_col', &
                                   axes(1:2), Time, &
                                   trim(tracername(nn)) //' column tendency from uw_conv wetdep', &
                                   trim(tracer_units(nn))//'*(kg/m2)/s', missing_value=mv)
           id_trevp_uwc(nn) = &
            register_diag_field (mod_name, trim(tracername(nn))//'dt_evp_uwc', &
                                    axes(1:3), Time, &
                                  trim(tracername(nn)) //' tendency from uw_conv precip_revap', &
                                  trim(tracer_units(nn))//'/s', missing_value=mv)
           id_trevp_uwd(nn) = &
            register_diag_field (mod_name, trim(tracername(nn))//'dt_evp_uwd', &
                                    axes(1:3), Time, &
                                  trim(tracername(nn)) //' tendency from uw_conv precip_revap', &
                                  trim(tracer_units(nn))//'/s', missing_value=mv)
        end do
     end if

    select case (tracer_check_type)
       case(1)
          call error_mesg ('uw_conv', &
	     'tracer checks: no min/max check, yes filling (WARNING! non-conservative)', NOTE)
       case(2)
          call error_mesg ('uw_conv', 'tracer checks: no min/max check, no filling', NOTE)
       case DEFAULT
          call error_mesg ('uw_conv', 'tracer checks: DEFAULT = min/max check, no filling', NOTE)
    end select

    module_is_initialized = .true.

    
  end SUBROUTINE UW_CONV_INIT

!#####################################################################
!#####################################################################

  subroutine uw_conv_end
    call exn_end_k
    call findt_end_k
    module_is_initialized = .FALSE.
  end subroutine uw_conv_end

!#####################################################################
!#####################################################################

  SUBROUTINE uw_conv(is, js, Time, tb, qv, ub, vb, pmid, pint,zmid,  & !input
       zint, q, omega, delt, pblht, ustar, bstar, qstar, sflx, lflx, land, coldT,& !input
       asol, tdt_rad, tdt_dyn, qdt_dyn, dgz_dyn, ddp_dyn, tdt_dif, qdt_dif, hmint, lat, lon, & !input
       cush, do_strat,  skip_calculation, max_available_cf,          & !input
       tten, qvten, qlten, qiten, qaten, qnten,                      & !output
       uten, vten, rain, snow,                                       & !output
       cmf, hlflx, qtflx, pflx, liq_pflx, ice_pflx, cldql, cldqi, cldqa,cldqn, &
       cbmfo, gusto, tkep, pblhto, rkmo, taudpo, exist_shconv, exist_dpconv, tracers, trtend, uw_wetdep)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     SHALLOW CONVECTION SCHEME
!     Described in Bretherton et. al (MWR, April 2004)
!     For info contact Ming Zhao: ming.zhao@noaa.gov
!
!     Inputs: see below
!
!     Outputs: see below
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    implicit none

    type(time_type), intent(in)  :: Time
    integer,         intent(in)  :: is, js
    real,            intent(in)  :: delt 

    real, intent(in), dimension(:,:,:)   :: ub,vb !wind profile (m/s)
    real, intent(in), dimension(:,:,:)   :: zint  !height@model interfaces(m)
    real, intent(in), dimension(:,:,:)   :: pint  !pressure@model interfaces(pa)
    real, intent(in), dimension(:,:,:)   :: tb    !temperature profile (K)
    real, intent(in), dimension(:,:,:)   :: qv    !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:,:) :: q     !specific humidity profile (kg/kg)
    real, intent(in), dimension(:,:,:)   :: pmid  !pressure@model mid-levels (pa)
    real, intent(in), dimension(:,:,:)   :: zmid  !height@model mid-levels (m)
    real, intent(in), dimension(:,:,:)   :: omega !omega (Pa/s)
    real, intent(in), dimension(:,:)     :: land  !land fraction
    real, intent(in), dimension(:,:,:)   :: max_available_cf !  largest
                                     ! realizable value for uw cld frac
                                   ! after accounting for deep cld frac
    logical,intent(in), dimension(:,:)   :: skip_calculation ! do not
                                                 ! calculate where .true.
    logical,intent(in)                   :: do_strat !logical flag
    logical,intent(in), dimension(:,:)   :: coldT    !logical flag

    real, intent(in),    dimension(:,:)  :: pblht, ustar, bstar, qstar, sflx, lflx, lat, lon !pbl height...
    real, intent(inout), dimension(:,:)  :: cush  ! convective scale height (m) 
    real, intent(inout), dimension(:,:)  :: hmint ! vertically-averaged MSE (J/kg)
    real, intent(in),  dimension(:,:,:)  :: tdt_rad, tdt_dyn, qdt_dyn, dgz_dyn, ddp_dyn, tdt_dif, qdt_dif!miz

    integer, intent(inout), dimension(:,:,:)  :: exist_shconv, exist_dpconv

    type(aerosol_type),  intent (in)     :: asol
   
    real, intent(out), dimension(:,:,:)  :: tten,qvten              ! T,qv tendencies
    real, intent(out), dimension(:,:,:)  :: qlten,qiten,qaten,qnten ! q tendencies
    real, intent(out), dimension(:,:,:)  :: uten,vten               ! u,v tendencies
   
    real, intent(out), dimension(:,:,:)  :: cldql,cldqi,cldqa, cldqn!in-updraft q
    real, intent(out), dimension(:,:,:)  :: cmf    ! mass flux at level above layer (kg/m2/s)
    real, intent(out), dimension(:,:,:)  :: pflx   ! precipitation flux removed from a layer
    real, intent(out), dimension(:,:,:)  :: liq_pflx   ! liq precipitation flux removed from a layer
    real, intent(out), dimension(:,:,:)  :: ice_pflx   ! solid precipitation flux removed from a layer
    real, intent(out), dimension(:,:,:)  :: hlflx ! theta_l flux
    real, intent(out), dimension(:,:,:)  :: qtflx ! qt  flux
    real, intent(out), dimension(:,:)    :: rain, snow
    real, intent(inout), dimension(:,:)  :: cbmfo, gusto, tkep, pblhto, rkmo, taudpo  ! cloud-base mass flux
    real, intent(in),  dimension(:,:,:,:)  :: tracers         ! env. tracers
    real, intent(out), dimension(:,:,:,:)  :: trtend          ! calculated tracer tendencies
    real, intent(out), dimension(:,:,:)  :: uw_wetdep       ! calculated wet depostion for tracers

    integer i, j, k, kl, klm, nk, naer, na, n, ksrc, kinv

    real rhos0j, pblrat, pblht_old, pblht_cur
    real hlsrc, thcsrc, qctsrc, tmp, tmp1, lofactor, crh_th, tvs, qvs, gust_new, gust_dis
    real zsrc, psrc, cbmf_shallow, cbmf_old, cbmf_deep, rkm_shallow, rkm_dp, cbmf_dp_frac
    real del_crh, dcrh, dcrh0, dpsum
    real pblfact, numx
    real, dimension(size(tb,1),size(tb,2)) :: &
         plcl,       &     ! pressure of lifting condensation level (Pa)
         zlcl,       &     ! depth of lifting condensation level (m)
         plfc,       &     ! pressure of level of free convection (Pa)
         plnb,       &     ! pressure of level of neutral buoyancy (Pa)
         cino,       &     ! cin (m2/s2)
         capeo,      &     ! cape(m2/s2)
         tkeo,       &     ! tke (m2/s2)
         wrelo,      &     ! release level vertical velocity (m/s)
         ufrco,      &     ! cloud-base updraft fraction
         zinvo,      &     ! surface driven mixed-layer height
         einso,      &     ! estimated inversion strength (K)
         denth,      &     
         dqtmp,      &
         dting,      &
         dcino,      &     ! dcin (m2/s2)
         dcapeo,     &     ! dcape(m2/s2)
         stime,      &
         dtime,      &
         cwfno,      &     
         ocode,      &
         xpsrc,      &
         xhlsrc,     &
         xqtsrc,     &
         crho,       &
         rkm_s,      &
         rkm_d,      &
         cush_d,     &
         fdp,        &
	 rhos,       &
	 lhflx,      &
	 shflx,      &
         hmint_old,  &
         hm_vadv0,   &
         hm_hadv0,   &
         hm_tot0,    &
         hm_total,   &
	 tdt_rad_int, tdt_dyn_int, tdt_dif_int, qdt_dyn_int, qdt_dif_int, &  
	 tdt_rad_pbl, tdt_dyn_pbl, tdt_dif_pbl, qdt_dyn_pbl, qdt_dif_pbl, &
	 tdt_rad_fre, tdt_dyn_fre, tdt_dif_fre, qdt_dyn_fre, qdt_dif_fre, &
         tdt_tot_pbl, tdt_tot_fre, cpool, bflux,dgz_dyn_int, ddp_dyn_int

    real, dimension(size(qtflx,1),size(qtflx,2),size(qtflx,3)) :: qtflx_up, qtflx_dn, omega_up, omega_dn, hm_vadv
    real, dimension(size(qtflx,1),size(qtflx,2),size(qtflx,3)) :: omgmc_up, ddp_dyn_hm, nqtflx
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: wuo,fero,fdro,fdrso, tten_pevap, qvten_pevap
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: qldet, qidet, qadet, cfq, peo, hmo, hms, abu, cmf_s

    real, dimension(size(tb,1),size(tb,2))            :: scale_uw, scale_tr
    real :: tnew, qtin, dqt, temp_1, temp_max, temp_min
    
    !f1p
    real, dimension(size(tracers,1), size(tracers,2), size(tracers,3), size(tracers,4)) :: trtend_nc, rn_diag

!========Option for deep convection=======================================
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: uten_d, vten_d, tten_d, &
         qvten_d, qlten_d, qiten_d, qaten_d, qnten_d, cmf_d, cbu_d, pflx_d, hlflx_d, qtflx_d, nqtflx_d, qtten_d, &
         wuo_d, fero_d, fdro_d, fdrso_d, cldql_d, cldqi_d, cldqa_d, cldqn_d, tten_pevap_d, qvten_pevap_d
    real, dimension(size(tb,1),size(tb,2)) :: rain_d, snow_d, cwfn_d
    real, dimension(size(tb,1),size(tb,2)) :: dcapedm_d, dcwfndm_d, denth_d, dting_d, dqtmp_d, cbmf_d
    real, dimension(size(tracers,1),size(tracers,2),size(tracers,3),size(tracers,4)) :: trevp_d, trevp_s
    real, dimension(size(tracers,3),size(tracers,4)) :: trtend_t, trwet_t
!f1p
    real, dimension(size(tracers,3),size(tracers,4)) :: trtend_t_nc, trwet_t_nc, rn
!
    type(randomNumberStream), dimension(size(tb,1), size(tb,2)) :: streams
    real, dimension(size(tb,1),size(tb,2)) :: rand
    integer :: iseed
    real    :: seedwts(8) = (/3000.,1000.,300.,100.,30.,10.,3.,1./)
    real    :: tempseed(8)
    integer :: thisseed(8)
    integer :: seedperm = 0
 !========Option for deep convection=======================================

    real, dimension(size(tb,3)) :: am1, am2, am3, am4, am5, qntmp
    
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: pmass    ! layer mass (kg/m2)
    real, dimension(size(tb,1),size(tb,2))            :: tempdiag ! temporary diagnostic variable
    real, dimension(size(tracers,1),size(tracers,2),size(tracers,3),size(tracers,4))  :: trwet 
    ! calculated tracer wet deposition tendencies

    integer imax, jmax, kmax
    integer kd, ntracers
    integer ktop_tmp, kbot_tmp
    real :: tten_intg, qvten_intg, cp_inv, half_delt
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: tten_forc, qten_forc
    real, dimension(size(tb,3)) :: tten_tmp, qvten_tmp !f1p
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: tten_rad
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: dissipative_heat
    real, dimension(size(tb,1),size(tb,2),size(tb,3)) :: tdt_dif_l
    real  :: latx, lonx, lat1b, lat1e, lon1b, lon1e, lat2b, lat2e, lon2b, lon2e
    
    logical used
    type(sounding)          :: sd, sd1
    type(adicloud)          :: ac, ac1
    type(cclosure)          :: cc, cc1
    type(cplume)            :: cp, cp1
    type(ctend)             :: ct, ct1
    type(cpnlist)           :: cpn,dpn
    type(deepc)             :: dpc
    integer ::  ier
    character(len=256) :: ermesg

    kd = size(tracers,3)
    ntracers = size(tracers,4)
    call sd_init_k(kd,ntracers,sd);
    call sd_init_k(kd,ntracers,sd1);
    call ac_init_k(kd,ac);
    call ac_init_k(kd,ac1);
    call cp_init_k(kd,ntracers,cp)
    call cp_init_k(kd,ntracers,cp1)
    call ct_init_k(kd,ntracers,ct)
    call ct_init_k(kd,ntracers,ct1)
    !pack namelist parameters into plume and closure structure
    cpn % do_qctflx_zero = do_qctflx_zero
    cpn % do_hlflx_zero  = do_hlflx_zero
    cpn % do_subcloud_flx= do_subcloud_flx
    cpn % do_detran_zero = do_detran_zero
    cpn % rle       = rle
    cpn % rpen      = rpen
    cpn % rmaxfrac  = rmaxfrac
    cpn % wmin      = wmin
    cpn % wmax      = wmax
    cpn % rbuoy     = rbuoy
    cpn % rdrag     = rdrag  
    cpn % frac_drs  = frac_drs
    cpn % bigc      = bigc    
    cpn % auto_th0  = auto_th0
    cpn % deltaqc0  = deltaqc0
    cpn % do_pdfpcp = do_pdfpcp
    cpn % do_pmadjt = do_pmadjt
    cpn % do_emmax  = do_emmax
    cpn % do_pnqv   = do_pnqv
    cpn % do_tten_max   = do_tten_max
    cpn % tten_max  = tten_max
    cpn % emfrac_max= emfrac_max
    cpn % auto_rate = auto_rate
    cpn % tcrit     = tcrit  
    cpn % cldhgt_max= cldhgt_max
    cpn % do_ice    = do_ice
    cpn % do_ppen   = do_ppen
    cpn % do_pevap  = do_pevap
    cpn % hcevap    = hcevap
    cpn % cfrac     = cfrac
    cpn % pblfac    = pblfac
    cpn % ffldep    = ffldep
    cpn % mixing_assumption= mixing_assumption
    cpn % mp_choice = mp_choice
    cpn % Nl_land   = Nl_land
    cpn % Nl_ocean  = Nl_ocean
    cpn % qi_thresh = qi_thresh
    cpn % r_thresh  = r_thresh
    cpn % peff_l    = peff_l
    cpn % peff_i    = peff_i
    cpn % t00       = t00
    cpn % rh0       = rh0
    cpn % do_forcedlifting= do_forcedlifting
    cpn % atopevap  = atopevap
    cpn % wtwmin_ratio = wmin_ratio*wmin_ratio
    cpn % do_auto_aero = do_auto_aero
    cpn % rad_crit = rad_crit
    cpn % wrel_min = wrel_min
    cpn % do_weffect = do_weffect
    cpn % weffect    = weffect
    cpn % use_online_aerosol = use_online_aerosol
    cpn % use_new_let = use_new_let
    cpn % use_lcl_only= use_lcl_only
    cpn % do_new_pevap= do_new_pevap
    cpn % stop_at_let = stop_at_let
    cpn % do_limit_wmax= do_limit_wmax
    cpn % plev_for = plev_for
    if (ntracers > 0) then
      allocate ( cpn%tracername   (ntracers) )
      allocate ( cpn%tracer_units (ntracers) )
      allocate ( cpn%wetdep       (ntracers) )
      cpn%tracername(:) = tracername(:)
      cpn%tracer_units(:) = tracer_units(:)
      cpn%wetdep(:)%scheme = wetdep(:)%scheme
      cpn%wetdep(:)%Henry_constant = wetdep(:)%Henry_constant
      cpn%wetdep(:)%Henry_variable = wetdep(:)%Henry_variable
      cpn%wetdep(:)%frac_in_cloud = wetdep(:)%frac_in_cloud
      cpn%wetdep(:)%alpha_r = wetdep(:)%alpha_r
      cpn%wetdep(:)%alpha_s = wetdep(:)%alpha_s
      cpn%wetdep(:)%Lwetdep = wetdep(:)%Lwetdep
      cpn%wetdep(:)%Lgas = wetdep(:)%Lgas
      cpn%wetdep(:)%Laerosol = wetdep(:)%Laerosol
      cpn%wetdep(:)%Lice = wetdep(:)%Lice
      allocate ( dpn%tracername   (ntracers) )
      allocate ( dpn%tracer_units (ntracers) )
      allocate ( dpn%wetdep       (ntracers) )
      dpn%tracername(:) = tracername(:)
      dpn%tracer_units(:) = tracer_units(:)
      dpn%wetdep(:)%scheme = wetdep(:)%scheme
      dpn%wetdep(:)%Henry_constant = wetdep(:)%Henry_constant
      dpn%wetdep(:)%Henry_variable = wetdep(:)%Henry_variable
      dpn%wetdep(:)%frac_in_cloud = wetdep(:)%frac_in_cloud
      dpn%wetdep(:)%alpha_r = wetdep(:)%alpha_r
      dpn%wetdep(:)%alpha_s = wetdep(:)%alpha_s
      dpn%wetdep(:)%Lwetdep = wetdep(:)%Lwetdep
      dpn%wetdep(:)%Lgas = wetdep(:)%Lgas
      dpn%wetdep(:)%Laerosol = wetdep(:)%Laerosol
      dpn%wetdep(:)%Lice = wetdep(:)%Lice
    endif
    call cpn_copy(cpn, dpn)

    cc  % igauss    = igauss
    cc  % rkfre     = rkfre
    cc  % rmaxfrac  = rmaxfrac
    cc  % wcrit_min = wcrit_min
    cc  % mass_fact = mass_fact
    cc  % rbuoy     = rbuoy
    cc  % tau_sh    = tau_sh
    cc  % do_old_cbmfmax = do_old_cbmfmax
!========Option for deep convection=======================================
    dpc % cbmf0               = cbmf0
    dpc % rkm_dp1             = rkm_dp1
    dpc % rkm_dp2             = rkm_dp2
    dpc % cbmf_dp_frac1       = cbmf_dp_frac1
    dpc % cbmf_dp_frac2       = cbmf_dp_frac2
    dpc % crh_th_ocean        = crh_th_ocean
    dpc % crh_th_land         = crh_th_land
    dpc % cape_th             = cape_th
    dpc % cin_th              = cin_th
    dpc % cwfn_th             = cwfn_th
    dpc % tau_dp              = tau_dp
    dpc % mixing_assumption_d = mixing_assumption_d
    dpc % do_ppen_d           = do_ppen_d
    dpc % rpen_d              = rpen_d
    dpc % do_pevap_d          = do_pevap_d
    dpc % cfrac_d             = cfrac_d
    dpc % hcevap_d            = hcevap_d
    dpc % pblfac_d            = pblfac_d
    dpc % ffldep_d            = ffldep_d
    dpc % frac_limit_d        = frac_limit_d
    dpc % dcapedm_th          = dcapedm_th
    dpc % dcwfndm_th          = dcwfndm_th
    dpc % do_forcedlifting_d  = do_forcedlifting_d
    dpc % lofactor_d          = lofactor_d
    dpc % auto_th0_d          = auto_th0_d
    dpc % tcrit_d             = tcrit_d
    dpc % peff_l_d            = peff_l_d
    dpc % peff_i_d            = peff_i_d
    dpc % gama                = gama
    dpc % tke0                = tke0
    dpc % hgt0                = hgt0
    dpc % do_cgust_dp         = do_cgust_dp
    dpc % cgust_choice        = cgust_choice
    dpc % tau_dp_fact         = tau_dp_fact
    dpc % cin_fact            = cin_fact
    dpc % wcrit_min_gust      = wcrit_min_gust
    dpc % src_choice_d        = src_choice_d
    numx                      = duration/delt
!========Option for deep convection=======================================
    imax  = size( tb, 1 )
    jmax  = size( tb, 2 )
    kmax  = size( tb, 3 )
    sd % kmax=kmax
    sd % plev_cin=plev_cin


    kl=kmax-1
    klm=kl-1

   !initialize 3D variables outside the loop

    tten=0.; qvten=0.; qlten=0.; qiten=0.; qaten=0.; qnten=0.;
    uten=0.; vten =0.; rain =0.; snow =0.; plcl =0.; plfc=0.; plnb=0.;  
    cldqa=0.; cldql=0.; cldqi=0.; cldqn=0.; zlcl=0.;
    hlflx=0.; qtflx=0.; nqtflx=0.; pflx=0.; am1=0.; am2=0.; am3=0.; am4=0.;
    tten_pevap=0.; qvten_pevap=0.;
    ice_pflx = 0. ; liq_pflx = 0.; qtflx_up=0.; qtflx_dn=0.; 
    omega_up=0.; omega_dn=0.; omgmc_up=0.;  ddp_dyn_hm=0.;
    hm_vadv0=0.; hm_hadv0=0.; hm_tot0=0.; hm_total=0.; tdt_dif_l=0.; 

    tdt_rad_int=0.; tdt_dyn_int=0.; tdt_dif_int=0.; qdt_dyn_int=0.; qdt_dif_int=0.; 
    tdt_rad_pbl=0.; tdt_dyn_pbl=0.; tdt_dif_pbl=0.; qdt_dyn_pbl=0.; qdt_dif_pbl=0.; 
    tdt_rad_fre=0.; tdt_dyn_fre=0.; tdt_dif_fre=0.; qdt_dyn_fre=0.; qdt_dif_fre=0.; 
    tdt_tot_pbl=0.; tdt_tot_fre=0.; cpool=0.; bflux=0.; scale_uw=1.; scale_tr=1.;
    dgz_dyn_int=0.; ddp_dyn_int=0.;

    cino=0.; capeo=0.; tkeo=0.; wrelo=0.; ufrco=0.; zinvo=0.; einso=0.; wuo=0.; peo=0.; 
    fero=0.; fdro=0.; fdrso=0.; cmf=0.; denth=0.;  dqtmp=0.; ocode=0; cmf_s=0.;
    dcapeo=0.; dcino=0.; xpsrc=0.; xhlsrc=0.; xqtsrc=0.; fdp=0.; rkm_s=0.;
    trtend=0.; qldet=0.; qidet=0.; qadet=0.; crho=0.; hmo=0.; hms=0.; abu=0.;
    trwet = 0.
    dting = 0.
    dissipative_heat = 0.; rhos=0; lhflx=0; shflx=0; 
    hmint_old=hmint; hmint=0;

    cbmf_shallow = 0.

    naer = size(asol%aerosol,4)

!========Option for deep convection=======================================
    tten_d=0.; qvten_d=0.; qlten_d=0.; qiten_d=0.; qaten_d=0.; qnten_d=0.;
    uten_d=0.; vten_d =0.; rain_d =0.; snow_d =0.; qtten_d=0.;
    trevp_d=0.; trevp_s=0.; cush_d=0.;
    cldqa_d=0.; cldql_d=0.; cldqi_d=0.; cldqn_d=0.;
    hlflx_d=0.; qtflx_d=0.; nqtflx_d=0.; pflx_d=0.;
    wuo_d=0.; fero_d=0.; fdro_d=0.; fdrso_d=0.; 
    cmf_d=0.; cbu_d=0.;
    denth_d=0.; dting_d=0.; dqtmp_d=0.; cbmf_d=0.; cwfn_d=0.;
    dcapedm_d=0.; dcwfndm_d=0.;
    dcino=0.; 
    tten_pevap_d=0.; qvten_pevap_d=0.; rkm_d=0.; 

!========Option for deep convection=======================================
    if (do_stochastic_rkm) then
      do j = 1, jmax
        do i=1, imax
          tempseed = tb(i,j,1)*seedwts
	  thisseed = nint(tempseed)
          streams(i,j) = initializeRandomNumberStream(ishftc(thisseed, seedperm))
          call getRandomNumbers( streams(i,j), rand(i,j) )
        enddo
      enddo
      rand(:,:)=(rand(:,:)-0.5)*2.
    end if

    do j = 1, jmax
       do i=1, imax

	 stime(i,j)=0;
	 dtime(i,j)=0;
	 do k=1,numx
	    stime(i,j)=stime(i,j)+exist_shconv(i,j,k)
	    dtime(i,j)=dtime(i,j)+exist_dpconv(i,j,k)
	 end do
	 stime(i,j)=stime(i,j)/numx
	 dtime(i,j)=dtime(i,j)/numx

      	 do k=numx,2,-1
	    exist_shconv(i,j,k) = exist_shconv(i,j,k-1)
	    exist_dpconv(i,j,k) = exist_dpconv(i,j,k-1)
	 end do
         exist_shconv(i,j,1) = 0
         exist_dpconv(i,j,1) = 0

       	 trtend_t=0.; trwet_t=0.;
         do k=1,kmax
           pmass(i,j,k) = (pint(i,j,k+1) - pint(i,j,k))/GRAV
         enddo
    !relaxation TKE back to 0 with time-scale of disscale
    !tkeavg = ustar(i,j)*bstar(i,j)*disscale 
    !dissipate tke with length-scale of disscale
    !tkeavg=(ustar(i,j)*bstar(i,j)*disscale)**(2./3.)
    !below following Holtslag and Boville 1993

	 if (pblht(i,j).lt.0.) then
            temp_1=0.0
         elseif (pblht(i,j).gt.5000.) then
            temp_1=5000.
         else
            temp_1=pblht(i,j)
         endif
         bflux(i,j) = 0.5*(0.6*ustar(i,j)*bstar(i,j)*temp_1)**(2./3.)
         temp_1=ustar(i,j)**3.+0.6*ustar(i,j)*bstar(i,j)*temp_1
         if (temp_1 .gt. 0.) temp_1 = 0.5*temp_1**(2./3.)
         tkeo(i,j) = MAX (tkemin, temp_1)

	 if (do_prog_tke) then
    	    pblht_old=pblhto(i,j)
    	    pblht_cur=pblht (i,j)
	    pblht_cur=min(max(pblht_cur,10.),5000.)
     	    pblrat=pblht_old/pblht_cur
    	    pblhto(i,j)=pblht_cur
            !Buoyancy flux B=rhos(i,j)*ustar(i,j)*bstar(i,j); unit:kg/m/s3 or (kg/m3 * m2/s3)
	    bflux(i,j)=ustar(i,j)*bstar(i,j) !buoyancy flux in kinematic unit: m2/s3 or (m/s2) *(m/2)
	    !solve tke implicitly
	    if (pblrat.gt.pblrat0 .or. bflux(i,j).le.0) then
	       tkep(i,j)=tkep(i,j)/(1.+delt/tau_tke)
            else
	       tkep(i,j)=((bfact*bflux(i,j)*delt+tkep(i,j))*pblrat)/(1.+delt/tau_tke)
	    endif
            tkep(i,j) = MAX (tkemin, tkep(i,j))
       	 endif

         cbmf_shallow=0. ! Set cbmf_shallow to avoid usage before assignment.
         if (skip_calculation(i,j)) then
           ocode(i,j) = 6
           go to 100
         endif
         call clearit(ac, cc, cp, ct, cp1, ct1);

! restrict grid-box area available to shallow convection to that which 
! is not involved with deep convection
          cp%maxcldfrac = minval(max_available_cf(i,j,:))
          cc%maxcldfrac = cp%maxcldfrac

          cc%scaleh = cush(i,j); 
          cush(i,j) = -1.;
          if(cc%scaleh.le.0.0) cc%scaleh=1000.

          am1(:) = 0.; am2(:) = 0.; am3(:) = 0.; am4(:) = 0.; am5(:) = 0.;

          do k=1,kmax
            tmp=1. / (zint(i,j,k)-zint(i,j,k+1)) * 1.0e9 * 1.0e-12
            if(use_online_aerosol) then
              do na = 1,naer
                if(asol%aerosol_names(na) == 'so4' .or. &
                   asol%aerosol_names(na) == 'so4_anthro' .or. &
                   asol%aerosol_names(na) == 'so4_natural') then
                           am1(k)=am1(k)+asol%aerosol(i,j,k,na)*tmp
                else if(asol%aerosol_names(na) == 'omphilic' .or. &
                        asol%aerosol_names(na) == 'omphobic') then
                           am4(k)=am4(k)+asol%aerosol(i,j,k,na)*tmp
                else if(asol%aerosol_names(na) == 'bcphilic' .or. &
                        asol%aerosol_names(na) == 'bcphobic' .or. &
                        asol%aerosol_names(na) == 'dust1' .or. &
                        asol%aerosol_names(na) == 'dust2' .or. &
                        asol%aerosol_names(na) == 'dust3' .or. &
                        asol%aerosol_names(na) == 'dust_mode1_of_2') then   !h1g, 2015-09-19
                           am2(k)=am2(k)+asol%aerosol(i,j,k,na)*tmp
                else if(asol%aerosol_names(na) == 'seasalt1' .or. &
                        asol%aerosol_names(na) == 'seasalt2' .or. &
                        asol%aerosol_names(na) == 'seasalt_aitken' .or. &   !h1g, 2015-09-19
                        asol%aerosol_names(na) == 'seasalt_fine'  ) then    !h1g, 2015-09-19
                           am3(k)=am3(k)+asol%aerosol(i,j,k,na)*tmp
                else if(asol%aerosol_names(na) == 'seasalt3' .or. &
                        asol%aerosol_names(na) == 'seasalt4' .or. &
                        asol%aerosol_names(na) == 'seasalt5' .or. &
                        asol%aerosol_names(na) == 'seasalt_coarse') then    !h1g, 2015-09-19
                           am5(k)=am5(k)+asol%aerosol(i,j,k,na)*tmp
                end if
              end do
              am2(k)=am2(k)+am3(k)+am4(k)
              if(.not. use_sub_seasalt) am3(k)=am3(k)+am5(k)
            else
              am1(k)= asol%aerosol(i,j,k,2)*tmp
              am2(k)= asol%aerosol(i,j,k,1)*tmp
              am3(k)= sea_salt_scale*asol%aerosol(i,j,k,5)*tmp
              am4(k)= om_to_oc*asol%aerosol(i,j,k,3)*tmp
            endif
          end do

!========Pack column properties into a sounding structure====================

          if (do_qn) then
             qntmp(:)=q(i,j,:,nqn)
          else
             qntmp(:)=0.
          end if
	  ksrc=1
	  tdt_dif_l(i,j,:)=tdt_dif(i,j,:)-tdt_rad(i,j,:);
          call pack_sd_k(land(i,j), coldT(i,j), delt, pmid(i,j,:), pint(i,j,:),     &
               zmid(i,j,:), zint(i,j,:), ub(i,j,:), vb(i,j,:), omega(i,j,:), tb(i,j,:), &
               qv(i,j,:), q(i,j,:,nql), q(i,j,:,nqi), q(i,j,:,nqa), qntmp,       &
               am1(:), am2(:), am3(:), am4(:), &
	       tdt_rad(i,j,:), tdt_dyn(i,j,:), qdt_dyn(i,j,:), dgz_dyn(i,j,:), ddp_dyn(i,j,:), &
	       tdt_dif_l(i,j,:), qdt_dif(i,j,:), src_choice, tracers(i,j,:,:), sd, Uw_p)

!========Finite volume intepolation==========================================

          gusto(i,j)=max(gusto(i,j),tkemin)
          sd%do_gust_qt = do_gust_qt
          sd%gqt_choice = gqt_choice
          sd%cgust     = gusto(i,j)
          sd%cgust0    = cgust0
          sd%cgust_max = cgust_max
          sd%sigma0    = sigma0
          sd%tke       = tkeo(i,j)
          sd%lat       = lat(i,j)*180/3.1415926
          sd%lon       = lon(i,j)*180/3.1415926

	  if (do_prog_tke .or. use_turb_tke ) sd%tke = tkep(i,j)   !h1g, 2015-08-11

          call extend_sd_k(sd, pblht(i,j), do_ice, Uw_p)

          zinvo(i,j) = sd%zinv
          kinv       = sd%kinv
	  einso(i,j) = (sd%hl(kinv+1)-sd%hl(kinv))/(sd%z(kinv+1)-sd%z(kinv))/cp_air

	  tmp=sd%thvbot(1)*sd%exners(1)
          rhos(i,j)= sd%ps(0)/(rdgas*tmp)
!         lhflx(i,j)=rhos(i,j)*ustar(i,j)*qstar(i,j)
!	  qvs=sd%qv(1)+sd%ssqct(1)*(sd%ps(0)-sd%p(1))
!         tvs=tmp*(1+0.608*qvs)
!         shflx(i,j)=Uw_p%cp_air*(rhos(i,j)*ustar(i,j)*bstar(i,j)*tvs/Uw_p%grav-0.608*tmp*lhflx(i,j))&
!                       / (1+0.608*qvs);
!	  lhflx(i,j)=lhflx(i,j)*Uw_p%hlv

          shflx(i,j)=sflx(i,j)
	  lhflx(i,j)=lflx(i,j)*Uw_p%hlv

	  hmint(i,j) = sd%hmint;
	  tdt_rad_int(i,j) = sd%tdt_rad_int;
	  tdt_dyn_int(i,j) = sd%tdt_dyn_int;
	  tdt_dif_int(i,j) = sd%tdt_dif_int;
	  qdt_dyn_int(i,j) = sd%qdt_dyn_int;
	  qdt_dif_int(i,j) = sd%qdt_dif_int;
	  dgz_dyn_int(i,j) = sd%dgz_dyn_int;
	  ddp_dyn_int(i,j) = sd%ddp_dyn_int;
	  tdt_rad_pbl(i,j) = sd%tdt_rad_pbl;
	  tdt_dyn_pbl(i,j) = sd%tdt_dyn_pbl;
	  tdt_dif_pbl(i,j) = sd%tdt_dif_pbl;
	  qdt_dyn_pbl(i,j) = sd%qdt_dyn_pbl;
	  qdt_dif_pbl(i,j) = sd%qdt_dif_pbl;
 	  tdt_rad_fre(i,j) = sd%tdt_rad_fre;
	  tdt_dyn_fre(i,j) = sd%tdt_dyn_fre;
	  tdt_dif_fre(i,j) = sd%tdt_dif_fre;
	  qdt_dyn_fre(i,j) = sd%qdt_dyn_fre;
	  qdt_dif_fre(i,j) = sd%qdt_dif_fre;
	  tdt_tot_pbl(i,j) = sd%tdt_tot_pbl;
	  tdt_tot_fre(i,j) = sd%tdt_tot_fre;
          hm_vadv0(i,j) = sd%hm_vadv0
	  hm_hadv0(i,j) = sd%tdt_dyn_int + sd%qdt_dyn_int + sd%dgz_dyn_int + ddp_dyn_int(i,j) - sd%hm_vadv0
!note: qdt_dif_int = lhflx; tdt_dif_int approximately equal to shflx
!         hm_tot0 (i,j) = hm_hadv0(i,j)+hm_vadv0(i,j)+tdt_rad_int(i,j)+shflx(i,j)+lhflx(i,j)
          hm_tot0 (i,j) = hm_hadv0(i,j)+hm_vadv0(i,j)+tdt_rad_int(i,j)+tdt_dif_int(i,j)+qdt_dif_int(i,j)

          hm_total(i,j) = (hmint(i,j)-hmint_old(i,j))/delt

          do k = 1,kmax
             nk = kmax+1-k
             qtflx_up(i,j,nk) = sd%qtflx_up(k)
             qtflx_dn(i,j,nk) = sd%qtflx_dn(k)
             omega_up(i,j,nk) = sd%omega_up(k)
             omega_dn(i,j,nk) = sd%omega_dn(k)
             hm_vadv(i,j,nk)  = sd%hm_vadv(k)
             ddp_dyn_hm(i,j,nk) = sd%ddp_dyn(k)
          enddo

!========Find source air, and do adiabatic cloud lifting======================
	  ksrc  =sd%ksrc    
          zsrc  =sd%zsrc
          psrc  =sd%psrc
          thcsrc=sd%thcsrc
          qctsrc=sd%qctsrc
          hlsrc =sd%hlsrc

          rkm_shallow=rkm_sh

	  lat1b= 30.; lat1e=40.; lon1b=260.; lon1e=270.
      	  lat2b=-20.; lat2e=10.; lon2b=285.; lon2e=305.
      	  latx=lat(i,j)*180/3.1415926; lonx=lon(i,j)*180/3.1415926
      	  if (latx.gt.lat1b .and. latx.lt.lat1e .and. lonx.gt.lon1b .and. lonx.lt.lon1e) then
             tmp=1
      	  elseif (latx.gt.lat2b .and. latx.lt.lat2e .and. lonx.gt.lon2b .and. lonx.lt.lon2e) then
             tmp=2
      	  endif


          if (do_lands) then
            !wstar   = (ustar(i,j)*bstar(i,j)*pblht(i,j))**(1./3.)
             cpn % auto_th0 = auto_th0 * (1. + landfact_m * sd%land)
             call qt_parcel_k (sd%qs(1), qstar(i,j), pblht(i,j), sd%tke, sd%land, 0.0, &
                  pblht0, 1.0, lofactor0, lochoice, qctsrc, lofactor)
             rkm_shallow = rkm_sh   * lofactor
          end if

          call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, do_fast, do_ice, ac)
          ac % usrc = sd%u(sd%ktoppbl)
          ac % vsrc = sd%v(sd%ktoppbl)
          !if (ac%plfc.eq.0) ac%plfc=psrc
          !if (ac%plnb.eq.0) ac%plnb=psrc
          plcl (i,j) = ac%plcl
          zlcl (i,j) = ac%zlcl-sd%zs(0)
          plfc (i,j) = ac%plfc
          plnb (i,j) = ac%plnb
          cino (i,j) = ac%cin
          capeo(i,j) = ac%cape
          xpsrc(i,j) = sd%p(ksrc);
          xhlsrc(i,j)= ac%hlsrc; !xhlsrc(i,j)= sd%qct(ksrc); 
          xqtsrc(i,j)= ac%qctsrc; 
          crho(i,j)  = sd%crh;
          do k = 1,kmax
             nk = kmax+1-k
             hmo  (i,j,nk) = sd%hm(k);
             hms  (i,j,nk) = sd%hms(k);
             abu  (i,j,nk) = ac%buo(k);
          end do

!          if (do_lclht .and. sd%land.gt.0.5) then
!	      tmp = hgt0/zlcl(i,j) !tmp = hgt0/sd%pblht
!             tmp = min(max(tmp,0.1),10.)
!             rkm_shallow = rkm_sh * tmp
!             cpn % auto_th0 = auto_th0 * (1. + landfact_m * sd%land)
!          end if

	  if (do_stime) then
             tmp=rkm_sh1-stime(i,j)*(rkm_sh1-rkm_sh)
 	     if (stime(i,j).eq.0 .or. rkmo(i,j).gt.tmp) then
	     	rkmo(i,j) = tmp
             endif
	     rkm_shallow=rkmo(i,j)
	  endif

!	  if (do_prog_gust) then
!	     rkm_shallow = rkm_shallow/max(gusto(i,j),1.0)
!          endif

!	  if (do_tke_cgust) then
!	     cpn % do_forcedlifting = do_forcedlifting
!	     cpn % do_ppen  = do_ppen
!             cc  % rkfre    = rkfre
!             sd  % tke      = sd%cgust
!             sd  % tke      = sd%tke   * (1. - gfact2     * sd%land)
!             cpn % auto_th0 = auto_th0 * (1. + landfact_m * sd%land)
!	     if (sd%land.gt.0.5) then
! 	     	if (sd%cgust.gt.sd%cgust0 .and. ac%cape.gt.ac%cin) then
!             	   ac  % cin = 0
!             	   cpn % do_forcedlifting = .true.
!	     	   cpn % do_ppen = .false.
!	           rkm_shallow = rkm_sh
!              	   tmp = sd%cgust/(sd%cgust0+sd%cgust)
!	     	   tmp = (1.-sqrt(tmp))
!             	   rkm_shallow = rkm_shallow * tmp * gfact3
!                   cc  % rkfre = min(cc % rkfre / tmp * gfact4, 10.)
!                endif
! 	     endif
!          endif

	  rkm_s(i,j) = rkm_shallow

          if (do_fast) then
             if (ac%klcl.eq.0 .or. ac%plcl.eq.sd%ps(1) .or. ac%plcl.lt.20000.) then
                ocode(i,j)=1; cbmf_shallow=0.; goto 100
             end if
             if (ac%plfc.lt.50000.) then
                ocode(i,j)=2; cbmf_shallow=0.; goto 100
             end if
          end if

!========Cumulus closure to determine cloud base mass flux===================

          cbmf_old=cbmfo(i,j); cc%cbmf=cbmf_old;

          if (iclosure.eq.0) then
             call cclosure_bretherton(sd%tke, cpn, sd, Uw_p, ac, cc)
          else if (iclosure.eq.1) then
             call cclosure_implicit(sd%tke, cpn, sd, Uw_p, ac, cc, delt, rkm_shallow, &
                  do_coldT, sd1, ac1, cc1, cp, ct, ier, ermesg) 
             if (ier /= 0) then
               call error_mesg ('subroutine uw_conv iclosure=1 ', ermesg, FATAL)
             endif
          else if (iclosure.eq.2) then
             call cclosure_relaxwfn(sd%tke, cpn, sd, Uw_p, ac, cc, cp, ct, delt,  &
                  rkm_shallow, do_coldT, sd1, ac1, cc1, cp1, ct1, ier, ermesg)
             if (ier /= 0) then
               call error_mesg ('subroutine uw_conv iclosure=2 ', ermesg, FATAL)
             endif
          end if

          cbmfo(i,j) = cc%cbmf
          wrelo(i,j) = cc%wrel
          ufrco(i,j) = cc%ufrc

          if (ac%klcl.eq.0 .or. ac%plcl.eq.sd%ps(1) .or. ac%plcl.lt.20000.) then
             ocode(i,j)=1; cbmf_shallow=0.; goto 100
          end if
          if (ac%plfc.lt.50000.) then
             ocode(i,j)=2; cbmf_shallow=0.; goto 100
          end if
          if(cc%cbmf.lt.1.e-6 .or. cc%wrel.eq.0.) then
             ocode(i,j)=3; cbmf_shallow=0.; goto 100
          end if

!========Do shallow cumulus plume calculation================================

          cbmf_shallow = cc%cbmf
          cpn%do_ppen=do_ppen
          cpn%rpen   =rpen
          call cumulus_plume_k(cpn, sd, ac, cp, rkm_shallow, cbmf_shallow, cc%wrel, cc%scaleh, Uw_p, ier, ermesg)
          if (ier /= 0) then
            call error_mesg ('subroutine uw_conv', ermesg, FATAL)
          endif
          if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
             ocode(i,j)=4; cbmf_shallow=0.; goto 100
          end if
          if(cp%cldhgt.ge.cldhgt_max) then
             ocode(i,j)=5; cbmf_shallow=0.; goto 100
          end if

          !cpn%isdeep=.false.
          !if (cpn%isdeep .EQV. .true.) then 
          !   fdp(i,j) = 1
          !else
          !   fdp(i,j) = 0
          !end if

!========Calculate cumulus produced tendencies===============================

          call cumulus_tend_k(cpn, sd, Uw_p, cp, ct, do_coldT)

!========Unpack convective tendencies========================================
          do k = 1,cp%ltop
             nk = kmax+1-k
             uten  (i,j,nk) = ct%uten (k)
             vten  (i,j,nk) = ct%vten (k)
             qlten (i,j,nk) = ct%qlten(k)
             qiten (i,j,nk) = ct%qiten(k)
             qaten (i,j,nk) = ct%qaten(k)
             qnten (i,j,nk) = ct%qnten(k)
             qldet (i,j,nk) = ct%qldet(k)
             qidet (i,j,nk) = ct%qidet(k)
             qadet (i,j,nk) = ct%qadet(k)
             qvten (i,j,nk) = ct%qvten(k)
             pflx  (i,j,nk) = ct%pflx (k)
             ice_pflx(i,j,nk) = cp%ppti(k)
             liq_pflx(i,j,nk) = cp%pptr(k)
             tten  (i,j,nk) = ct%tten (k)
             rhos0j = sd%ps(k)/(rdgas*0.5*(cp%thvbot(k+1)+cp%thvtop(k))*sd%exners(k))
             hlflx(i,j,nk) = ct%hlflx(k)
             qtflx (i,j,nk)  = ct%qctflx(k-1)
             nqtflx(i,j,nk)  = ct%nqtflx(k-1)
             tten_pevap (i,j,nk) = ct%tevap (k)
             qvten_pevap(i,j,nk) = ct%qevap (k)
             
             cldqa (i,j,nk) = cp%ufrc(k)
             cldql (i,j,nk) = cp%qlu(k)
             cldqi (i,j,nk) = cp%qiu(k)
             cldqn (i,j,nk) = cp%qnu(k)
             cmf_s (i,j,nk) = cp%umf(k)
             cmf   (i,j,nk) = cp%umf(k)
             wuo   (i,j,nk) = cp%wu (k)
             peo   (i,j,nk) = cp%peff(k)
             fero  (i,j,nk) = cp%fer(k)
             fdro  (i,j,nk) = cp%fdr(k)
             fdrso (i,j,nk) = cp%fdrsat(k)*cp%umf(k)

             do n = 1, size(trtend,4)
              trevp_s(i,j,nk,n) = ct%trevp(k,n)
             enddo
          enddo
          cush  (i,j)  = cp%cush
          snow  (i,j)  = ct%snow
          rain  (i,j)  = ct%rain
          denth (i,j)  = ct%denth
          dqtmp (i,j)  = ct%dqtmp
          dting (i,j)  = ct%dting
          cpool (i,j)  = ct%cpool

! make sure the predicted tracer tendencies do not produce negative
! tracers due to convective tendencies. if necessary, adjust the 
! tendencies.
          if (do_deep) then
            trtend_t = ct%trten
            trwet_t  = ct%trwet
          else
!f1p
!            call check_tracer_realizability (kmax, size(trtend,4), delt, &
!                                           cp%tr, ct%trten, ct%trwet) 
             call check_tracer_realizability (kmax, size(trtend,4), delt, &
                                           cp%tr, ct%trten, ct%trwet, pmass(i,j,:), tracer_check_type, rn = rn )
            do k = 1,cp%ltop
              nk = kmax+1-k
              do n = 1, size(trtend,4)
                trtend(i,j,nk,n) = ct%trten(k,n) + ct%trwet(k,n)
                trwet(i,j,nk,n)  = ct%trwet(k,n)
                rn_diag(i,j,nk,n) = rn(k,n)
              enddo
            enddo
	  end if

	  !if (cp%cush .gt. cush_ref) then
100	  if (cbmf_shallow .gt. 0.0) then
	     exist_shconv(i,j,1) = 1
	  endif

!========Option for deep convection=======================================
          if (do_deep) then
	     cbmf_deep = 0.
	     rkm_dp = dpc%rkm_dp1 
	     crh_th = sd%land*dpc%crh_th_land+(1.-sd%land)*dpc%crh_th_ocean
             tmp = max(min (sd%crh, 1.0), 0.0)
	     del_crh = tmp - crh_th
             dcrh0  = 1.0001-crh_th
	     if (del_crh .gt. 0.) then
	        cbmf_deep = 0.0001 !first assuming existence of deep convective cloud base mass flux
	        dcrh = del_crh/dcrh0
		dcrh = dcrh**(1./norder)
	        rkm_dp  = dpc%rkm_dp1 + dcrh * (dpc%rkm_dp2 - dpc%rkm_dp1)
                lofactor= 1.- sd%land*(1.- dpc%lofactor_d) !option for introducing land difference
	        if (do_lod_rkm) then
               	   rkm_dp       = rkm_dp  * lofactor
	        elseif (do_lod_cfrac) then
		   dpc % cfrac_d= cfrac_d * lofactor
	        elseif (do_lod_tcrit) then
             	   dpc % tcrit_d= tcrit_d * lofactor
	        elseif (do_lod_cape) then
             	   dpc % cape_th= cape_th * lofactor
	        elseif (do_lod_tau) then
	      	   dpc % tau_dp = tau_dp  * lofactor
 	        end if

		if (do_stochastic_rkm) then
                 !rkm_dp = dpc%rkm_dp1 + rand(i,j)*(dpc%rkm_dp2 -dpc%rkm_dp1)
                 rkm_dp = rkm_dp * frac_rkm_pert * rand(i,j)
                end if

!		if (do_lclht .and. sd%land.gt.0.5) then
!		   tmp = hgt0/zlcl(i,j) !tmp = hgt0/sd%pblht
!             	   tmp = min(max(tmp,0.1),10.)
!                  rkm_dp = dpc%rkm_dp2*tmp
!                  dpc % tau_dp = tau_dp * tmp
!                  dpc % auto_th0_d = auto_th0_d * (1. + landfact_m * sd%land)
!		end if
	     end if

	     if (do_dtime) then
                tmp=tau_dp*tau_dp_fact
                tmp=tmp-dtime(i,j)*(tmp-tau_dp)
 	     	if (dtime(i,j).eq.0 .or. taudpo(i,j).gt.tmp) then
	     	   taudpo(i,j) = tmp
                endif
	        dpc%tau_dp = taudpo(i,j)
                if (dtime(i,j).lt.dtime0 .and. stime(i,j).lt.stime0) then
                   cbmf_deep = 0
                endif
	     endif

             dpn % do_ppen  = dpc % do_ppen_d
	     dpn % rpen     = dpc % rpen_d
             dpn % do_pevap = dpc % do_pevap_d
             dpn % cfrac    = dpc % cfrac_d
             dpn % hcevap   = dpc % hcevap_d
             dpn % pblfac   = dpc % pblfac_d
             dpn % ffldep   = dpc % ffldep_d
             dpn % tcrit    = dpc % tcrit_d
             dpn % auto_th0 = dpc % auto_th0_d
             dpn % peff_l   = dpc % peff_l_d
             dpn % peff_i   = dpc % peff_i_d
             dpn % mixing_assumption = dpc % mixing_assumption_d
             dpn % do_forcedlifting  = dpc % do_forcedlifting_d

             if (idpchoice.eq.0) then
                call  dpconv0(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                      rkm_dp, cbmf_deep, cp1, ct1, ocode(i,j), ier, ermesg)
             else if (idpchoice.eq.1) then
                call  dpconv1(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                      rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode(i,j), dcapedm_d(i,j), ier, ermesg)
             else if (idpchoice.eq.2) then
                call  dpconv2(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                      rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode(i,j), taudpo(i,j), &
		      dcwfndm_d(i,j), dcapedm_d(i,j), cwfn_d(i,j), lat(i,j), lon(i,j), ier, ermesg)
             end if
             if (ier /= 0) then
                call error_mesg ('uw_conv calling dpconv', ermesg, FATAL)
             endif

    	     if (cbmf_deep.eq.0) then
       	     	call cp_clear_k(cp1); cp1%maxcldfrac=1.; cp1%cush=-1;
       		call ct_clear_k(ct1);
		if (do_forced_conv) & !try if there is forced convection due to cold pool
	     	   call conv_forced(dpc, dpn, Uw_p, sd, ac, do_coldT, do_ice, rkm_dp, cbmf_deep,&
                                    cp1, ct1, lat(i,j), lon(i,j), ier, ermesg)
             end if
             if (cbmf_deep .eq.0) then 
             	fdp(i,j) = 0
             else
          	fdp(i,j) = 1
             end if
	     !if (cp1%cush .gt. cush_ref) then
	     if (cbmf_deep .gt. 0.0) then
	     	exist_dpconv(i,j,1) = 1
	     endif

             do k = 1,kmax !cp1%ltop
                nk = kmax+1-k
                uten_d  (i,j,nk) = ct1%uten (k)
                vten_d  (i,j,nk) = ct1%vten (k)
                qlten_d (i,j,nk) = ct1%qlten(k)
                qiten_d (i,j,nk) = ct1%qiten(k)
                qaten_d (i,j,nk) = ct1%qaten(k) 
                qnten_d (i,j,nk) = ct1%qnten(k) 
                qvten_d (i,j,nk) = ct1%qvten(k)
                qtten_d (i,j,nk) = ct1%qctten(k)
                pflx_d  (i,j,nk) = ct1%pflx (k)
                tten_d  (i,j,nk) = ct1%tten (k)
                hlflx_d (i,j,nk) = ct1%hlflx(k) 
                qtflx_d (i,j,nk) = ct1%qctflx(k-1)
                nqtflx_d(i,j,nk) = ct1%nqtflx(k-1)
                cldqa_d (i,j,nk) = cp1%ufrc(k)
                cldql_d (i,j,nk) = cp1%qlu(k)
                cldqi_d (i,j,nk) = cp1%qiu(k)
                cldqn_d (i,j,nk) = cp1%qnu(k)
                cmf_d   (i,j,nk) = cp1%umf(k) + cp1%emf(k)
                cbu_d   (i,j,nk) = cp1%buo(k)
                tten_pevap_d (i,j,nk) = ct1%tevap (k)
                qvten_pevap_d(i,j,nk) = ct1%qevap (k)
                wuo_d   (i,j,nk) = cp1%wu (k)
                fero_d  (i,j,nk) = cp1%fer(k)
                fdro_d  (i,j,nk) = cp1%fdr(k) 
                fdrso_d (i,j,nk) = cp1%fdrsat(k)*cp1%fdr(k)*cp1%umf(k)
                do n = 1, size(trtend,4)
                   trevp_d(i,j,nk,n) = ct1%trevp(k,n)
             	enddo
             enddo
             snow_d  (i,j)  = ct1%snow
             rain_d  (i,j)  = ct1%rain
             cbmf_d  (i,j)  = cbmf_deep
             denth_d (i,j)  = ct1%denth
             dting_d (i,j)  = ct1%dting
             dqtmp_d (i,j)  = ct1%dqtmp
             cush_d  (i,j)  = cp1%cush
             rkm_d   (i,j)  = rkm_dp;

             trtend_t = trtend_t+ct1%trten
             trwet_t  = trwet_t +ct1%trwet
!<f1p
             trtend_t_nc = trtend_t 
             trwet_t_nc  = trwet_t
!>
!
!             call check_tracer_realizability (kmax, size(trtend,4), delt, &
!                                              cp1%tr, trtend_t, trwet_t)
             call check_tracer_realizability (kmax, size(trtend,4), delt, &
                                              cp1%tr, trtend_t, trwet_t, pmass(i,j,:), tracer_check_type, rn = rn    )
             do k = 1,kmax!cp1%ltop
               nk = kmax+1-k
               do n = 1, size(trtend,4)
                 trtend(i,j,nk,n) = trtend_t(k,n) + trwet_t(k,n)
                 trwet(i,j,nk,n)  = trwet_t(k,n)
!f1p
                 trtend_nc(i,j,nk,n) = trtend_t_nc(k,n) + trwet_t_nc(k,n)
                 rn_diag(i,j,nk,n) = rn(k,n)
!>
	       enddo
             enddo
            
             uten  (i,j,:) = uten  (i,j,:) + uten_d  (i,j,:)
             vten  (i,j,:) = vten  (i,j,:) + vten_d  (i,j,:)
             qlten (i,j,:) = qlten (i,j,:) + qlten_d (i,j,:)
             qiten (i,j,:) = qiten (i,j,:) + qiten_d (i,j,:)
             qaten (i,j,:) = qaten (i,j,:) + qaten_d (i,j,:) 
             qnten (i,j,:) = qnten (i,j,:) + qnten_d (i,j,:) 
             qvten (i,j,:) = qvten (i,j,:) + qvten_d (i,j,:)
             pflx  (i,j,:) = pflx  (i,j,:) + pflx_d  (i,j,:)
             tten  (i,j,:) = tten  (i,j,:) + tten_d  (i,j,:)
             hlflx (i,j,:) = hlflx (i,j,:) + hlflx_d (i,j,:) 
             qtflx (i,j,:) = qtflx (i,j,:) + qtflx_d (i,j,:)
             nqtflx(i,j,:) = nqtflx(i,j,:) + nqtflx_d(i,j,:)
             cmf   (i,j,:) = cmf_s (i,j,:) + cmf_d   (i,j,:)
             tten_pevap (i,j,:)=tten_pevap (i,j,:) + tten_pevap_d (i,j,:) 
             qvten_pevap(i,j,:)=qvten_pevap(i,j,:) + qvten_pevap_d(i,j,:) 

             cldql (i,j,:) = cldql (i,j,:) + cldql_d(i,j,:)
             cldqi (i,j,:) = cldqi (i,j,:) + cldqi_d(i,j,:)
	     do k = 1,kmax
	     	cldqa (i,j,k) = max(cldqa (i,j,k),cldqa_d(i,j,k))
   	     end do

             snow  (i,j)  = snow  (i,j) + snow_d  (i,j)
             rain  (i,j)  = rain  (i,j) + rain_d  (i,j)
             denth (i,j)  = denth (i,j) + denth_d (i,j)
             dting (i,j)  = dting (i,j) + dting_d (i,j)
             dqtmp (i,j)  = dqtmp (i,j) + dqtmp_d (i,j)
             cpool (i,j)  = cpool (i,j) + ct1%cpool
             !cbmfo (i,j)  = cc%cbmf
             !cwfno (i,j)  = cc%cwfn
          end if
!========End of do_deep, Option for deep convection=======================================

	  if (do_prog_gust) then
	     gusto(i,j)=(gusto(i,j)+gfact*cpool(i,j)*delt)/(1+delt/tau_gust)

!	     tmp  =sd%thvbot(1)*sd%exners(1)
!	     tvs  =tmp*(1+0.608*sd%qv(1))
!             tmp=(sd%ps(0)-ac%plcl)/Uw_p%grav !unit:kg/m2
!             !Buoyancy flux B=rhos(i,j)*ustar(i,j)*bstar(i,j); unit:kg/m/s3 or (kg/m3 * m2/s3)
!	     bflux(i,j)=ustar(i,j)*bstar(i,j)  !buoyancy flux in kinematic unit: m2/s3 or (m/s2) *(m/2)
!	     bflux(i,j)=Uw_p%cp_air*rhos(i,j)*bflux(i,j)*tvs/Uw_p%grav !unit:W/m2; 
!	     bflux(i,j)=bflux(i,j)/tmp         !unit: W/kg or m2/s3
!	     gust_new=(cpool(i,j)+gfact1*bflux(i,j))*sd%delt !unit m2/s2
!	     gust_dis=gusto(i,j)/tau_gust*sd%delt
!             gusto(i,j)=gusto(i,j)+gust_new-gust_dis

!	     pblht_cur=min(max(pblht_cur,10.),5000.)
!     	     pblrat=pblht_old/pblht_cur
!	     !solve gusto implicitly
!	     !if (pblrat.gt.2. .or. pblht_cur.eq.10) then
!	     if (pblrat.gt.5. .or. bflux(i,j).le.0) then
!		 gusto(i,j)=gusto(i,j)/(1.+sd%delt/tau_gust)
!		 !gusto(i,j)=(gfact1*bflux(i,j)*sd%delt+gusto(i,j))*pblrat
!             else
!		 gusto(i,j)=((gfact1*bflux(i,j)*sd%delt+gusto(i,j))*pblrat)/(1.+sd%delt/tau_gust)
!	     endif
!
!             cbmfo(i,j)=pblht_cur
	  endif

!subtract parameterized convective mass flux
          do k = 1,kmax
	     tmp=cmf(i,j,k)*Uw_p%grav;
	     if ((-omega_up(i,j,k).gt.tmp) .and. (tmp.gt.0)) then
             	omgmc_up(i,j,k) = omega_up(i,j,k)+tmp;
 	     else
             	omgmc_up(i,j,k) = omega_up(i,j,k);
	     endif		 
          enddo

       enddo
    enddo

    call sd_end_k(sd)
    call sd_end_k(sd1)
    call ac_end_k(ac)
    call ac_end_k(ac1)
    call cp_end_k(cp)
    call cp_end_k(cp1)
    call ct_end_k(ct)
    call ct_end_k(ct1)
    if (_ALLOCATED ( cpn%tracername    ))  deallocate ( cpn%tracername    )
    if (_ALLOCATED ( cpn%tracer_units  ))  deallocate ( cpn%tracer_units  )
    if (_ALLOCATED ( cpn%wetdep        ))  deallocate ( cpn%wetdep        )
    if (_ALLOCATED ( dpn%tracername    ))  deallocate ( dpn%tracername    )
    if (_ALLOCATED ( dpn%tracer_units  ))  deallocate ( dpn%tracer_units  )
    if (_ALLOCATED ( dpn%wetdep        ))  deallocate ( dpn%wetdep        )

    if (do_uwcmt) then
        half_delt = delt*0.5
        cp_inv    = 1. / Uw_p%cp_air
	dissipative_heat(:,:,:) = -((ub(:,:,:) + half_delt*uten(:,:,:))*uten(:,:,:) + &
                                    (vb(:,:,:) + half_delt*vten(:,:,:))*vten(:,:,:))*cp_inv
	tten(:,:,:) = tten(:,:,:) + dissipative_heat(:,:,:)
    else 
    	uten=0.;
    	vten=0.;
    end if

    if ( prevent_unreasonable ) then
      temp_min=300.
      temp_max=200.
      scale_uw=HUGE(1.0)
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            if ((q(i,j,k,nqa) + qaten(i,j,k)*delt) .lt. 0. .and. (qaten(i,j,k).ne.0.)) then
              qaten(i,j,k) = -1.*q(i,j,k,nqa)/delt
            end if
            if ((q(i,j,k,nqa) + qaten(i,j,k)*delt) .gt. 1. .and. (qaten(i,j,k).ne.0.)) then
              qaten(i,j,k)= (1. - q(i,j,k,nqa))/delt
            end if
 
            if ((q(i,j,k,nql) + qlten(i,j,k)*delt) .lt. 0. .and. (qlten(i,j,k).ne.0.)) then
              tten (i,j,k) = tten(i,j,k) -(q(i,j,k,nql)/delt+qlten(i,j,k))*HLv/Cp_Air
              qvten(i,j,k) = qvten(i,j,k)+(q(i,j,k,nql)/delt+qlten(i,j,k))
              qlten(i,j,k) = qlten(i,j,k)-(q(i,j,k,nql)/delt+qlten(i,j,k))
            end if
 
            if ((q(i,j,k,nqi) + qiten(i,j,k)*delt) .lt. 0. .and. (qiten(i,j,k).ne.0.)) then
              tten (i,j,k) = tten(i,j,k) -(q(i,j,k,nqi)/delt+qiten(i,j,k))*HLs/Cp_Air
              qvten(i,j,k) = qvten(i,j,k)+(q(i,j,k,nqi)/delt+qiten(i,j,k))
    
              qiten(i,j,k) = qiten(i,j,k)-(q(i,j,k,nqi)/delt+qiten(i,j,k))
            end if

            if (do_qn) then
              if ((q(i,j,k,nqn) + qnten(i,j,k)*delt) .lt. 0. .and. (qnten(i,j,k).ne.0.)) then
                qnten(i,j,k) = qnten(i,j,k)-(q(i,j,k,nqn)/delt+qnten(i,j,k))
              end if
            endif
    !rescaling to prevent negative specific humidity for each grid point
            if (do_rescale) then
              qtin =  qv(i,j,k)
              dqt  =  qvten(i,j,k) * delt
              if ( dqt.lt.0 .and. qtin+dqt.lt.1.e-10 ) then
                temp_1 = max( 0.0, -(qtin-1.e-10)/dqt )
              else
                temp_1 = 1.0
              endif
    !scaling factor for each column is the minimum value within that column
              scale_uw(i,j) = min( temp_1, scale_uw(i,j))
            endif
    !rescaling to prevent excessive temperature tendencies
            if (do_rescale_t) then
              temp_max = max(temp_max,tb(i,j,k))
              temp_min = min(temp_min,tb(i,j,k))
              tnew  =  tb(i,j,k) + tten(i,j,k) * delt
              if ( tnew > 363.15 ) then
                temp_1 = 0.0
                print *, 'WARNING: setting scale_uw to zero to prevent large T tendencies in UW'
                print *, i,j,'lev=',k,'pressure=',pmid(i,j,k),'tb=',tb(i,j,k),'tten=',tten(i,j,k)*delt
                print *, 'lat=', sd%lat, 'lon=', sd%lon, 'land=',sd%land
              else
                temp_1 = 1.0
              endif
    !scaling factor for each column is the minimum value within that column
              scale_uw(i,j) = min( temp_1, scale_uw(i,j))
            endif
          enddo
        enddo
      enddo
!      if (temp_max > 350. .or. temp_min < 170.) then
!      	 print *, 'temp_min=',temp_min, 'temp_max=',temp_max
!      endif

!     where ((tracers(:,:,:,:) + trtend(:,:,:,:)*delt) .lt. 0.)
!        trtend(:,:,:,:) = -tracers(:,:,:,:)/delt
!     end where

      if (do_rescale .or. do_rescale_t) then
      !scale tendencies
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              uten (i,j,k)  = scale_uw(i,j) * uten (i,j,k)
              vten (i,j,k)  = scale_uw(i,j) * vten (i,j,k)
              tten (i,j,k)  = scale_uw(i,j) * tten (i,j,k)
              qvten(i,j,k)  = scale_uw(i,j) * qvten(i,j,k)
              qlten(i,j,k)  = scale_uw(i,j) * qlten(i,j,k)
              qiten(i,j,k)  = scale_uw(i,j) * qiten(i,j,k)
              qaten(i,j,k)  = scale_uw(i,j) * qaten(i,j,k)
              if (do_qn) qnten(i,j,k) = scale_uw(i,j) * qnten(i,j,k)
              if (k.eq.kmax) then
                rain(i,j) = scale_uw(i,j) * rain(i,j)
                snow(i,j) = scale_uw(i,j) * snow(i,j)
                rain_d(i,j)=scale_uw(i,j) * rain_d(i,j)
                snow_d(i,j)=scale_uw(i,j) * snow_d(i,j)
                cpool(i,j)= scale_uw(i,j) * cpool(i,j)
              endif
            end do
          end do
        end do
      end if
    endif


    do k=1,kmax
      do j=1,jmax
        do i=1,imax
          cfq(i,j,k) = 0
          if (wuo(i,j,k) .gt. 0.) then
            cfq(i,j,k) = 1
          endif
        enddo
      enddo
    enddo

    if (do_imposing_rad_cooling) then
      tten_rad(:,:,:) = 0.0
      do j = 1,jmax
         do i=1,imax
	   do k = 1,kmax
              if (tb(i,j,k) .gt. t_thresh) then
              	tten_rad (i,j,k) = cooling_rate/86400.
 	      else
		tten_rad (i,j,k) = (t_strato-tb(i,j,k))/(tau_rad*86400.)
              end if
           enddo
	 enddo
      enddo
!      used = send_data( id_tten_rad_uwc,tten_rad*aday,Time, is, js, 1)
    end if

    if (do_imposing_forcing) then
      tten_forc(:,:,:) = 0.0
      qten_forc(:,:,:) = 0.0
      do j = 1,jmax
      	do i=1,imax
	   tten_forc(i,j,:)=0; 
	   qten_forc(i,j,:)=0;
           if (use_klevel) then
             k = klevel
	     tten_forc(i,j,k)=tdt_rate/86400.
	     qten_forc(i,j,k)=qdt_rate/86400.
           else
	     kbot_tmp=1; 
 	     ktop_tmp=kmax;
 	     do k=1,kmax
	       if (pmid(i,j,k)>=7500) then
	         ktop_tmp=k
	       end if
	       if (pmid(i,j,k)>=85000) then
		 kbot_tmp=k
	       end if
	     enddo
	     do k = kbot_tmp,ktop_tmp
	     	if (pmid(i,j,k)>pres_min .and. pmid(i,j,k)<=pres_max) then
		   tten_forc(i,j,k)=tdt_rate/86400.
	    	   qten_forc(i,j,k)=qdt_rate/86400.
                end if
    	     enddo
	   end if
  	enddo
      enddo
      used = send_data( id_tdt_forc_uwc,tten_forc*aday,Time, is, js, 1)
      used = send_data( id_qdt_forc_uwc,qten_forc*aday,Time, is, js, 1)
    end if

    !diagnostic output
    used = send_data( id_xpsrc_uwc,        xpsrc,              Time, is, js)
    used = send_data( id_xhlsrc_uwc,       xhlsrc,             Time, is, js)
    used = send_data( id_xqtsrc_uwc,       xqtsrc,             Time, is, js)
    used = send_data( id_tdt_pevap_uwc,    tten_pevap*aday , Time, is, js, 1)
    used = send_data( id_qdt_pevap_uwc,    qvten_pevap*aday, Time, is, js, 1)

    used = send_data( id_tdt_uwc,    tten*aday , Time, is, js, 1)
    used = send_data( id_qdt_uwc,    qvten*aday, Time, is, js, 1)
    used = send_data( id_udt_uwc,    uten*aday , Time, is, js, 1)
    used = send_data( id_vdt_uwc,    vten*aday,  Time, is, js, 1)
    used = send_data( id_cmf_uwc,    cmf,          Time, is, js, 1)
    used = send_data( id_cmf_uws,    cmf_s,        Time, is, js, 1)
    used = send_data( id_cfq_uwc,    cfq,          Time, is, js, 1)
    used = send_data( id_wu_uwc,     wuo,          Time, is, js, 1)
    used = send_data( id_peo_uwc,    peo,          Time, is, js, 1)
    used = send_data( id_fer_uwc,    fero,         Time, is, js, 1)
    used = send_data( id_fdr_uwc,    fdro,         Time, is, js, 1)
    used = send_data( id_fdrs_uwc,   fdrso,        Time, is, js, 1)
    used = send_data( id_cqa_uwc,    cldqa,        Time, is, js, 1)
    used = send_data( id_cql_uwc,    cldql,        Time, is, js, 1)
    used = send_data( id_cqi_uwc,    cldqi,        Time, is, js, 1)
    used = send_data( id_cqn_uwc,    cldqn,        Time, is, js, 1)
    used = send_data( id_hlflx_uwc,  hlflx,        Time, is, js, 1)
    used = send_data( id_qtflx_uwc,  qtflx,        Time, is, js, 1)
    used = send_data( id_nqtflx_uwc, nqtflx,       Time, is, js, 1)
    used = send_data( id_qtflx_up_uwc, qtflx_up,   Time, is, js, 1)
    used = send_data( id_qtflx_dn_uwc, qtflx_dn,   Time, is, js, 1)
    used = send_data( id_omgmc_up_uwc, omgmc_up,   Time, is, js, 1)
    used = send_data( id_omega_up_uwc, omega_up,   Time, is, js, 1)
    used = send_data( id_omega_dn_uwc, omega_dn,   Time, is, js, 1)
    used = send_data( id_hm_vadv_uwc, hm_vadv,     Time, is, js, 1)
    used = send_data( id_pflx_uwc,   pflx,         Time, is, js, 1)
    used = send_data( id_hmo_uwc,    hmo,          Time, is, js, 1)
    used = send_data( id_hms_uwc,    hms,          Time, is, js, 1)
    used = send_data( id_abu_uwc,    abu,          Time, is, js, 1)

    used = send_data( id_tdt_rad_uwc,  tdt_rad,  Time, is, js, 1)
    used = send_data( id_tdt_dyn_uwc,  tdt_dyn,  Time, is, js, 1)
    used = send_data( id_tdt_dif_uwc,  tdt_dif,  Time, is, js, 1)
    used = send_data( id_qdt_dyn_uwc,  qdt_dyn,  Time, is, js, 1)
    used = send_data( id_qdt_dif_uwc,  qdt_dif,  Time, is, js, 1)
    used = send_data( id_dgz_dyn_uwc,  dgz_dyn,  Time, is, js, 1)
    used = send_data( id_ddp_dyn_uwc,  ddp_dyn_hm,  Time, is, js, 1)

    used = send_data( id_lhflx_uwc,    lhflx,    Time, is, js )
    used = send_data( id_shflx_uwc,    shflx,    Time, is, js )
    used = send_data( id_hmint_uwc,    hmint,    Time, is, js )
    used = send_data( id_hm_vadv0_uwc, hm_vadv0, Time, is, js )
    used = send_data( id_hm_hadv0_uwc, hm_hadv0, Time, is, js )
    used = send_data( id_hm_tot0_uwc,  hm_tot0,  Time, is, js )
    used = send_data( id_hm_total_uwc, hm_total, Time, is, js )

    used = send_data( id_tdt_rad_int,  tdt_rad_int, Time, is, js )
    used = send_data( id_tdt_dyn_int,  tdt_dyn_int, Time, is, js )
    used = send_data( id_tdt_dif_int,  tdt_dif_int, Time, is, js )
    used = send_data( id_qdt_dyn_int,  qdt_dyn_int, Time, is, js )
    used = send_data( id_qdt_dif_int,  qdt_dif_int, Time, is, js )
    used = send_data( id_dgz_dyn_int,  dgz_dyn_int, Time, is, js )
    used = send_data( id_ddp_dyn_int,  ddp_dyn_int, Time, is, js )

!    used = send_data( id_tdt_rad0_uwc, tdt_rad0, Time, is, js )
 
    used = send_data( id_prec_uwc, (rain+snow-rain_d-snow_d)*aday, Time, is, js )
    used = send_data( id_snow_uwc, (snow-snow_d)*aday,      Time, is, js )
    used = send_data( id_cin_uwc,  (cino),             Time, is, js )
    used = send_data( id_cape_uwc, (capeo),            Time, is, js )
    used = send_data( id_gust_uwc, (gusto),            Time, is, js )
    used = send_data( id_cpool_uwc,(cpool),            Time, is, js )
    used = send_data( id_bflux_uwc,(bflux),            Time, is, js )
    used = send_data( id_crh_uwc,  (crho),             Time, is, js )
    used = send_data( id_pblht_uwc,(pblht),            Time, is, js )
    used = send_data( id_tke_uwc,  (tkeo),             Time, is, js )
    used = send_data( id_tkep_uwc, (tkep),             Time, is, js )
    used = send_data( id_cbmf_uwc, (cbmfo),            Time, is, js )
    used = send_data( id_wrel_uwc, (wrelo),            Time, is, js )
    used = send_data( id_ufrc_uwc, (ufrco),            Time, is, js )
    used = send_data( id_plcl_uwc, (plcl*0.01),        Time, is, js )
    used = send_data( id_zlcl_uwc, (zlcl),             Time, is, js )
    used = send_data( id_plfc_uwc, (plfc*0.01),        Time, is, js )
    used = send_data( id_plnb_uwc, (plnb*0.01),        Time, is, js )
    used = send_data( id_zinv_uwc, (zinvo),            Time, is, js )
    used = send_data( id_cush_uwc, (cush),             Time, is, js )
    used = send_data( id_dcin_uwc, (dcino),            Time, is, js )
    used = send_data( id_dcape_uwc,(dcapeo),           Time, is, js )
!    used = send_data( id_dwfn_uwc, (dwfno),            Time, is, js )
    used = send_data( id_enth_uwc, (denth),            Time, is, js )
    used = send_data( id_qtmp_uwc, (dqtmp),            Time, is, js )
    used = send_data( id_dting_uwc,(dting),            Time, is, js )
    used = send_data( id_ocode_uwc,(ocode),            Time, is, js )
    used = send_data( id_fdp_uwc,  (fdp),              Time, is, js )
    used = send_data( id_rkm_uwc,  (rkm_s),            Time, is, js )
    used = send_data( id_stime_uwc,(stime),            Time, is, js )
    used = send_data( id_scale_uwc,(scale_uw),         Time, is, js )
    used = send_data( id_scaletr_uwc,(scale_tr),       Time, is, js )

    if ( do_uwcmt ) then
      used = send_data( id_tdt_diss_uwc,  dissipative_heat*aday , Time, is, js, 1)
    end if

    if ( do_strat ) then
       used = send_data( id_qldt_uwc, qlten*aday,    Time, is, js, 1)
       used = send_data( id_qidt_uwc, qiten*aday,    Time, is, js, 1)
       used = send_data( id_qadt_uwc, qaten*aday,    Time, is, js, 1)
       used = send_data( id_qndt_uwc, qnten*aday,    Time, is, js, 1)
       used = send_data( id_qldet_uwc,  qldet*aday,  Time, is, js, 1)
       used = send_data( id_qidet_uwc,  qidet*aday,  Time, is, js, 1)
       used = send_data( id_qadet_uwc,  qadet*aday,  Time, is, js, 1)
       used = send_data( id_qtdt_uwc,(qvten+qlten+qiten)*aday,Time, is, js, 1)
    end if
!f1p
    if ( allocated(id_rn) ) then
       do n = 1,size(id_rn)
          used = send_data( id_rn(n), rn_diag(:,:,:,n), Time, is, js, 1)
       end do
    end if

    if ( allocated(id_tracerdt_uwc) ) then
       do n = 1,size(id_tracerdt_uwc)
          used = send_data( id_tracerdt_uwc(n), trtend(:,:,:,n), Time, is, js, 1)
       end do
    end if
!f1p
    if ( allocated(id_tracerdt_uwc_nc) ) then
       do n = 1,size(id_tracerdt_uwc_nc)
          used = send_data( id_tracerdt_uwc_nc(n), trtend_nc(:,:,:,n), Time, is, js, 1)
       end do
    end if

    if ( allocated(id_tracerdt_uwc_col) ) then
       do n = 1,size(id_tracerdt_uwc_col)
          if ( id_tracerdt_uwc_col(n) > 0 ) then
            tempdiag = 0.
            do k = 1,kmax
               tempdiag(:,:) = tempdiag(:,:) + trtend(:,:,k,n) * pmass(:,:,k)
            end do
            used = send_data( id_tracerdt_uwc_col(n), tempdiag(:,:), Time, is, js)
          end if
       end do
    end if
!f1p
    if ( allocated(id_tracerdt_uwc_col_nc) ) then
       do n = 1,size(id_tracerdt_uwc_col_nc)
          if ( id_tracerdt_uwc_col_nc(n) > 0 ) then
            tempdiag = 0.
            do k = 1,kmax
               tempdiag(:,:) = tempdiag(:,:) + trtend_nc(:,:,k,n) * pmass(:,:,k)
            end do
            used = send_data( id_tracerdt_uwc_col_nc(n), tempdiag(:,:), Time, is, js)
          end if
       end do
    end if

    if ( allocated(id_tracerdtwet_uwc) ) then
       do n = 1,size(id_tracerdtwet_uwc)
          used = send_data( id_tracerdtwet_uwc(n), trwet(:,:,:,n), Time, is, js, 1)
       end do
    end if

!<<<fp
!this means that uw_wetdep = 0 if wet_uwc_col is not requested for tracer n

!     if ( allocated(id_tracerdtwet_uwc_col) ) then
!        uw_wetdep = 0.
!        do n = 1,size(id_tracerdtwet_uwc_col)
!           if ( id_tracerdtwet_uwc_col(n) > 0 ) then
!              tempdiag = 0.
!              do k = 1,kmax
!                tempdiag(:,:) = tempdiag(:,:) + trwet(:,:,k,n) * pmass(:,:,k)
!             end do
!             used = send_data( id_tracerdtwet_uwc_col(n), tempdiag(:,:), Time, is, js)
!             uw_wetdep(:,:,n) = tempdiag(:,:)
!           end if
!        end do
!     end if

!now calculated uw for all tracers

    uw_wetdep = 0.
    do n=1,ntracers
       tempdiag = 0.
       do k = 1,kmax
          tempdiag(:,:) = tempdiag(:,:) + trwet(:,:,k,n) * pmass(:,:,k)
       end do
       uw_wetdep(:,:,n) = tempdiag(:,:)       
    end do

     if ( allocated(id_tracerdtwet_uwc_col) ) then
        do n = 1,size(id_tracerdtwet_uwc_col)
           if ( id_tracerdtwet_uwc_col(n) > 0 ) then
             used = send_data( id_tracerdtwet_uwc_col(n), uw_wetdep(:,:,n), Time, is, js)
           end if
        end do
     end if

!>>>

    if ( allocated(id_trevp_uwc) ) then
       do n = 1,size(id_trevp_uwc)
          used = send_data( id_trevp_uwc(n), trevp_s(:,:,:,n)+trevp_d(:,:,:,n), Time, is, js, 1)
       end do
    end if

!========Option for deep convection=======================================
    if (do_deep) then
       used=send_data( id_tdt_pevap_uwd,    tten_pevap_d*aday , Time, is, js, 1)
       used=send_data( id_qdt_pevap_uwd,    qvten_pevap_d*aday, Time, is, js, 1)
       used=send_data( id_tdt_uwd,   tten_d*aday , Time, is, js, 1)
       used=send_data( id_qdt_uwd,   qvten_d*aday, Time, is, js, 1)
       used=send_data( id_qtdt_uwd,  qtten_d*aday, Time, is, js, 1)
       used=send_data( id_cmf_uwd,   cmf_d,          Time, is, js, 1)
       used=send_data( id_cbu_uwd,   cbu_d,          Time, is, js, 1)
       used=send_data( id_wu_uwd,    wuo_d,          Time, is, js, 1)
       used=send_data( id_fer_uwd,   fero_d,         Time, is, js, 1)
       used=send_data( id_fdr_uwd,   fdro_d,         Time, is, js, 1)
       used=send_data( id_fdrs_uwd,  fdrso_d,        Time, is, js, 1)
       used=send_data( id_cqa_uwd,   cldqa_d,        Time, is, js, 1)
       used=send_data( id_cql_uwd,   cldql_d,        Time, is, js, 1)
       used=send_data( id_cqi_uwd,   cldqi_d,        Time, is, js, 1)
       used=send_data( id_cqn_uwd,   cldqn_d,        Time, is, js, 1)
       used=send_data( id_hlflx_uwd, hlflx_d,        Time, is, js, 1)
       used=send_data( id_qtflx_uwd, qtflx_d,        Time, is, js, 1)
       used=send_data( id_nqtflx_uwd,nqtflx_d,       Time, is, js, 1)
!       used=send_data( id_trtend_uwd, trtend,        Time, is, js, 1)
!       used=send_data( id_trwet_uwd,  trwet,         Time, is, js, 1)
      
       used=send_data( id_prec_uwd, (rain_d+snow_d)*aday,Time, is, js )
       used=send_data( id_snow_uwd, (snow_d)*aday,       Time, is, js )
       used=send_data( id_cbmf_uwd, (cbmf_d),              Time, is, js )
       used=send_data( id_dcapedm_uwd,(dcapedm_d),         Time, is, js )
       used=send_data( id_dcwfndm_uwd,(dcwfndm_d),         Time, is, js )
       used=send_data( id_cwfn_uwd, (cwfn_d),              Time, is, js )
       used=send_data( id_taudp_uwd,(taudpo),              Time, is, js )
       used=send_data( id_cush_uwd, (cush_d),              Time, is, js )
       used=send_data( id_enth_uwd, (denth_d),             Time, is, js )
       used=send_data( id_rkm_uwd,  (rkm_d),               Time, is, js )
       used=send_data( id_rand_uwd, (rand),                Time, is, js )
       used=send_data( id_dtime_uwd,(dtime),               Time, is, js )
             
       if ( do_strat ) then
          used=send_data( id_qldt_uwd, qlten_d*aday,     Time, is, js, 1)
          used=send_data( id_qidt_uwd, qiten_d*aday,     Time, is, js, 1)
          used=send_data( id_qadt_uwd, qaten_d*aday,     Time, is, js, 1)
       end if
       if ( allocated(id_trevp_uwd) ) then
         do n = 1,size(id_trevp_uwd)
           used = send_data( id_trevp_uwd(n), trevp_d(:,:,:,n), Time, is, js, 1)
         end do
       end if
    end if
!========Option for deep convection=======================================


    if (.not.apply_tendency) then
       uten=0.; vten=0.; tten=0.; qvten=0.; cmf=0.; rain=0.; snow=0.;
       qlten=0.; qiten=0.; qaten=0.; qnten=0.;
    end if

    if (zero_out_conv_area) then
       cldql=0.; cldqi=0.; cldqa=0.;
    end if

    if (do_imposing_rad_cooling) then
      do j = 1,jmax
         do i=1,imax
    	   tten(i,j,:) = tten (i,j,:) + tten_rad (i,j,:)
	 enddo
       enddo
    end if

    if (do_imposing_forcing) then
       do j = 1,jmax
       	  do i=1,imax
             tten (i,j,:) = tten (i,j,:) + tten_forc(i,j,:)
             qvten(i,j,:) = qvten(i,j,:) + qten_forc(i,j,:)
  	  enddo
	enddo
    end if


  END SUBROUTINE UW_CONV

!#####################################################################
!#####################################################################

  subroutine clearit(ac, cc, cp, ct, cp1, ct1)

    type(adicloud), intent(inout) :: ac
    type(cclosure), intent(inout) :: cc
    type(cplume),   intent(inout) :: cp,cp1
    type(ctend),    intent(inout) :: ct,ct1

    call ac_clear_k(ac); 
    ac%klcl =0;  ac%klfc =0;  ac%klnb =0; ac%zlcl=0;

    cc%wrel=0.; cc%ufrc=0.; cc%scaleh=0.;

    call cp_clear_k(cp);  cp%maxcldfrac =1.;
    call ct_clear_k(ct);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);

  end subroutine clearit


!#####################################################################

end MODULE UW_CONV_MOD
