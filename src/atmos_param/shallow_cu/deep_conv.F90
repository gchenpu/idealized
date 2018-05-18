
MODULE DEEP_CONV_MOD

  use      fms_mod,         only : write_version_number
  use      Constants_Mod,   ONLY : tfreeze,HLv,HLf,HLs,CP_AIR,GRAV,Kappa,rdgas,rvgas
  use  conv_utilities_mod,  only : uw_params_init
  use  conv_utilities_k_mod,only : sd_init_k, sd_copy_k, sd_end_k,  &
                                   ac_init_k, ac_clear_k, ac_end_k, &
                                   pack_sd_k, adi_cloud_k, extend_sd_k,&
                                   exn_init_k, exn_end_k, findt_init_k,&
                                   findt_end_k, qt_parcel_deep_k, &
                                   adicloud, sounding, uw_params, findt_k, &
				   exn_k, qsat_k, erfccc

  use  conv_plumes_k_mod,   only : cp_init_k, cp_end_k, cp_clear_k, &
                                   ct_init_k, ct_end_k, ct_clear_k, &
                                   cumulus_tend_k, cumulus_plume_k, &
                                   cplume, ctend, cpnlist

  use  conv_closures_mod,   only : cclosure_bretherton,   &
                                   cclosure_relaxcbmf, &
                                   cclosure_relaxwfn,  &
                                   cclosure_implicit, cclosure

!---------------------------------------------------------------------
  implicit none
  private
!----------- ****** VERSION NUMBER ******* ---------------------------

  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

!-------  interfaces --------

  public  :: cpn_copy, dpconv0, dpconv1, dpconv2, dpconv3, dpconv4, conv_forced

  logical         :: module_is_initialized = .false.

  character(len=7) :: mod_name = 'deep_conv'

  public deepc
  type deepc
     real    :: cbmf0
     real    :: rkm_dp1
     real    :: rkm_dp2
     real    :: cbmf_dp_frac1
     real    :: cbmf_dp_frac2
     real    :: crh_th_land
     real    :: crh_th_ocean
     real    :: cape_th 
     real    :: cin_th 
     real    :: tau_dp 
     real    :: rpen_d
     integer :: mixing_assumption_d
     logical :: do_ppen_d
     logical :: do_pevap_d
     real    :: cfrac_d
     real    :: hcevap_d
     real    :: pblfac_d
     real    :: ffldep_d
     real    :: dcapedm_th
     real    :: dcwfndm_th
     real    :: frac_limit_d
     real    :: lofactor_d
     real    :: tcrit_d
     real    :: auto_th0_d
     real    :: peff_l_d
     real    :: peff_i_d
     integer :: src_choice_d
     real    :: gama
     real    :: tke0
     real    :: hgt0
     real    :: cwfn_th
!     real    :: cgust0
     real    :: cin_fact
     real    :: wcrit_min_gust
     real    :: tau_dp_fact
     integer :: cgust_choice
     logical :: do_cgust_dp
     logical :: do_forcedlifting_d
  end type deepc

contains

!#####################################################################
!#####################################################################

  subroutine cpn_copy(cpn, dpn)
    type(cpnlist), intent(in)    :: cpn
    type(cpnlist), intent(inout) :: dpn

    dpn % do_qctflx_zero     = cpn % do_qctflx_zero
    dpn % do_hlflx_zero      = cpn % do_hlflx_zero
    dpn % do_subcloud_flx    = cpn % do_subcloud_flx
    dpn % use_lcl_only       = cpn % use_lcl_only
    dpn % do_new_pevap       = cpn % do_new_pevap
    dpn % stop_at_let        = cpn % stop_at_let
    dpn % do_limit_wmax      = cpn % do_limit_wmax
    dpn % plev_for           = cpn % plev_for
    dpn % do_detran_zero     = cpn % do_detran_zero
    dpn % rle                = cpn % rle
    dpn % rpen               = cpn % rpen
    dpn % rmaxfrac           = cpn % rmaxfrac
    dpn % wmin               = cpn % wmin
    dpn % wmax               = cpn % wmax
    dpn % rbuoy              = cpn % rbuoy
    dpn % rdrag              = cpn % rdrag  
    dpn % frac_drs           = cpn % frac_drs
    dpn % bigc               = cpn % bigc    
    dpn % auto_th0           = cpn % auto_th0
    dpn % deltaqc0           = cpn % deltaqc0
    dpn % do_pdfpcp          = cpn % do_pdfpcp
    dpn % do_pmadjt          = cpn % do_pmadjt
    dpn % do_emmax           = cpn % do_emmax
    dpn % do_pnqv            = cpn % do_pnqv
    dpn % emfrac_max         = cpn % emfrac_max
    dpn % auto_rate          = cpn % auto_rate
    dpn % tcrit              = cpn % tcrit  
    dpn % cldhgt_max         = cpn % cldhgt_max
    dpn % do_ice             = cpn % do_ice
    dpn % do_ppen            = cpn % do_ppen
    dpn % do_pevap           = cpn % do_pevap
    dpn % hcevap             = cpn % hcevap
    dpn % cfrac              = cpn % cfrac
    dpn % pblfac             = cpn % pblfac
    dpn % ffldep             = cpn % ffldep
    dpn % mixing_assumption  = cpn % mixing_assumption
    dpn % mp_choice          = cpn % mp_choice
    dpn % Nl_land            = cpn % Nl_land
    dpn % Nl_ocean           = cpn % Nl_ocean
    dpn % qi_thresh          = cpn % qi_thresh
    dpn % r_thresh           = cpn % r_thresh
    dpn % peff_l             = cpn % peff_l
    dpn % peff_i             = cpn % peff_i
    dpn % peff               = cpn % peff
    dpn % t00                = cpn % t00
    dpn % rh0                = cpn % rh0
    dpn % do_forcedlifting   = cpn % do_forcedlifting
    dpn % atopevap           = cpn % atopevap
    dpn % wtwmin_ratio       = cpn % wtwmin_ratio
    dpn % do_auto_aero       = cpn % do_auto_aero
    dpn % rad_crit           = cpn % rad_crit
    dpn % wrel_min           = cpn % wrel_min
    dpn % do_weffect         = cpn % do_weffect
    dpn % weffect            = cpn % weffect
    dpn % use_online_aerosol = cpn % use_online_aerosol
    dpn % isdeep             = cpn % isdeep
    dpn % use_new_let        = cpn % use_new_let
    dpn % do_tten_max        = cpn % do_tten_max
    dpn % tten_max           = cpn % tten_max

  end subroutine cpn_copy

!#####################################################################
!#####################################################################

  subroutine dpconv0(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, cp1, ct1, ocode, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(in)     :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: zcldtop, dcrh, cbmf_dp_frac

    ier = 0
    ermesg = ' '
    if ( (ocode.ne.0 .and. ocode.ne.4) .or. (cbmf_deep.eq.0) ) then
       ocode=6; return
    end if

    zcldtop = sd%z(cp%ltop)
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, cc%wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       ocode=6; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if
 
  end subroutine dpconv0


!#####################################################################
!#####################################################################

  subroutine dpconv1(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode, dcapedm, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(inout)  :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode, dcapedm
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: zcldtop, wrel, cbmf0, cbmf_max, tmp
    real          :: zs0, ps0, hl0, thc0
    integer       :: ksrc
    real          :: zsrc, psrc, thcsrc, qctsrc, hlsrc

    ier = 0
    ermesg = ' '
    zcldtop = 2000 !sd%z(cp%ltop)
    wrel = max(cc%wrel, 0.1)

   if (dpc%do_cgust_dp) then
       !if ((ac%cape.gt.ac%cin) .and. (sd%land.gt.0.5)) then
       if ((ac%cape.gt.0) .and. (sd%land.gt.0.5)) then
       	  dpn%do_forcedlifting = .true.
          dpn%do_ppen = .false.
       endif
   endif

!    if ( (ocode.ne.0 .and. ocode.ne.4) .or. (cbmf_deep.eq.0) ) then
    if ( cbmf_deep.eq.0 ) then
       ocode=6; return
    end if
    if (ac%cape .lt. dpc%cape_th) then
       cbmf_deep = 0.; ocode=7; return
    end if

    cbmf0 = 0.0001
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       cbmf_deep = 0.; ocode=8; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if

    call sd_copy_k(sd, sd1)
    tmp      = 1.
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt * tmp
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt * tmp
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt * tmp
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt * tmp
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt * tmp
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt * tmp
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt * tmp
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt * tmp

    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)

    call ac_clear_k(ac1);

    ksrc  =sd1%ksrc
    zsrc  =sd1%zs (ksrc)
    psrc  =sd1%ps (ksrc)
    thcsrc=sd1%thc(ksrc)
    qctsrc=sd1%qct(ksrc)
    hlsrc =sd1%hl (ksrc)
    call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd1, Uw_p, .false., do_ice, ac1)

    dcapedm=(ac%cape-ac1%cape)/cbmf0

    if (dcapedm .lt. dpc%dcapedm_th) then
       cbmf_deep=0.; ocode=9; 
       call ct_clear_k(ct1);
       return
    else
       cbmf_deep= (ac%cape - dpc%cape_th) / dcapedm / (dpc%tau_dp/sd%delt)
    end if

    cbmf_max  = (sd%ps(0) - sd%ps(cp1%krel))*(dpc%frac_limit_d/sd%delt)/Grav
    cbmf_deep = max(min(cbmf_deep, cbmf_max), 0.)
 
    if(cbmf_deep.lt.1.e-10) then 
       cbmf_deep=0.; ocode=10; 
       call ct_clear_k(ct1);
       return
    end if

    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       cbmf_deep = 0.; ocode=11;
       call ct_clear_k(ct1);
       return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if
    
  end subroutine dpconv1

!#####################################################################
!#####################################################################

  subroutine dpconv2(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode, &
                     taudp, dcwfndm, dcapedm, cwfn, lat, lon, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(inout)  :: sd
    type(adicloud),  intent(inout)  :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode, taudp, dcwfndm, dcapedm
    real,            intent(out)    :: cwfn
    real,            intent(in)     :: lat, lon
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: zcldtop, wrel, wrel0, cbmf0, cwfn_new, cwfn_th, cbmf_max, tmp, ufrc
    integer       :: ksrc, k, let, krel
    real          :: latx, lonx, lat1b, lat1e, lon1b, lon1e, lat2b, lat2e, lon2b, lon2e
    real          :: zsrc, psrc, thcsrc, qctsrc, hlsrc, lofactor

    ier = 0
    ermesg = ' '
    cwfn = 0.; dcwfndm=0.; dcapedm=0.; wrel0=0.1;
    taudp = dpc%tau_dp;

    if ( cbmf_deep.eq.0 ) then
       ocode=6; return
    end if

    cwfn_th = dpc%cwfn_th;
    if (dpc%do_cgust_dp) then
       if (dpc%cgust_choice==0) then
       	  if (ac%cape.gt.dpc%cape_th .and. cc%wrel.le.0. .and.     &
	      ac%cape.gt.dpc%cin_fact*ac%cin .and. sd%cgust.gt.sd%cgust0) then
       	      dpn%do_forcedlifting = .true.
              dpn%do_ppen  = .false.;
              lofactor     = dpc%lofactor_d
              rkm_dp       = rkm_dp        * lofactor
              dpn % peff_l = dpn % peff_l  / lofactor
              dpn % peff_i = dpn % peff_i  / lofactor
 	  endif
       else if (dpc%cgust_choice==1 .and. sd%land.gt.0.5) then
       	  if (ac%cape.gt.dpc%cape_th .and. cc%wrel.le.0. .and.     &
	      ac%cape.gt.dpc%cin_fact*ac%cin .and. sd%cgust.gt.sd%cgust0) then
       	      dpn%do_forcedlifting = .true.
              dpn%do_ppen  = .false.;
              lofactor     = 1.- sd%land*(1.- dpc%lofactor_d)
              rkm_dp       = rkm_dp        * lofactor
              !dpn % peff_l = dpn % peff_l  / lofactor
              !dpn % peff_i = dpn % peff_i  / lofactor
 	  endif
       else if (dpc%cgust_choice==2 .and. sd%land.gt.0.5) then
       	  if (ac%cape.gt.dpc%cape_th .and. cc%wrel.le.0. .and.     &
	      ac%cape.gt.dpc%cin_fact*ac%cin .and. sd%cgust.gt.sd%cgust0) then
       	      dpn%do_forcedlifting = .true.
              dpn%do_ppen  = .false.;
              lofactor     = 1.- sd%land*(1.- dpc%lofactor_d)
              rkm_dp       = rkm_dp        * lofactor
              !dpn % peff_l = dpn % peff_l  / lofactor
              !dpn % peff_i = dpn % peff_i  / lofactor
 	  endif
       else if (dpc%cgust_choice==3 .and. sd%land.gt.0.5) then
       	  if (ac%cape.gt.dpc%cape_th .and. cc%wrel.le.0. .and.     &
	      ac%cape.gt.dpc%cin_fact*ac%cin .and. sd%cgust.gt.sd%cgust0) then
       	      dpn%do_forcedlifting = .true.
              dpn%do_ppen  = .false.;
              lofactor     = 1.- sd%land*(1.- dpc%lofactor_d)
              rkm_dp       = rkm_dp        * lofactor
              dpn % peff_l = dpn % peff_l  / lofactor
              dpn % peff_i = dpn % peff_i  / lofactor
 	  endif
       else if (dpc%cgust_choice==5) then !saved only for testing purpose
	  lat1b= 30.; lat1e=40.; lon1b=260.; lon1e=270.
	  lat2b=-20.; lat2e=10.; lon2b=285.; lon2e=305.
      	  latx=lat*180/3.1415926; lonx=lon*180/3.1415926
      	  if (latx.gt.lat1b .and. latx.lt.lat1e .and. lonx.gt.lon1b .and. lonx.lt.lon1e) then
             tmp=1
      	  elseif (latx.gt.lat2b .and. latx.lt.lat2e .and. lonx.gt.lon2b .and. lonx.lt.lon2e) then
             tmp=2
      	  endif
          sd%do_gust_qt = .true.
          sd%gqt_choice = 0
    	  call extend_sd_k(sd,sd%pblht, do_ice, Uw_p)
   	  ksrc  =sd%ksrc
	  zsrc  =sd%zsrc
	  psrc  =sd%psrc
	  thcsrc=sd%thcsrc
	  qctsrc=sd%qctsrc
	  hlsrc =sd%hlsrc
    	  call ac_clear_k(ac);
    	  call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, .false., do_ice, ac)
       endif
    end if

    zcldtop = 2000 !sd%z(cp%ltop)
    wrel  = max(cc%wrel, wrel0)
    cbmf0 = dpc%cbmf0
    call cp_clear_k(cp);  cp %maxcldfrac=1.;
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
       cbmf_deep = 0.; ocode=8; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp, ct1, do_coldT)
    end if

    call sd_copy_k(sd, sd1)
    tmp      = 1.
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt * tmp
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt * tmp
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt * tmp
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt * tmp
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt * tmp
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt * tmp
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt * tmp
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt * tmp

    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)

    ksrc  =sd1%ksrc
    zsrc  =sd1%zsrc
    psrc  =sd1%psrc
    thcsrc=sd1%thcsrc
    qctsrc=sd1%qctsrc
    hlsrc =sd1%hlsrc

    call ac_clear_k(ac1);
    call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd1, Uw_p, .false., do_ice, ac1)

    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)

    krel=max(cp%krel,cp1%krel)
    let =min(cp%let, cp1%let)

    cwfn=0.; !unit: m2/s2 or J/kg
    do k = krel,let-1
       cwfn = cwfn + cp%buo(k)/cp%thvtop(k)*cp%dp(k)/sd%rho(k)
    end do

    cwfn_new=0.;
    do k = krel,let-1
       cwfn_new = cwfn_new + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)/sd1%rho(k)
    end do

    dcapedm=(ac%cape-ac1%cape)/cbmf0
    dcwfndm =(cwfn - cwfn_new)/cbmf0

    if (dcwfndm .lt. dpc%dcwfndm_th .or. cwfn.lt.0) then
       cbmf_deep=0.; ocode=9; return
    else
       cbmf_deep= (cwfn-cwfn_th) / dcwfndm / (taudp/sd%delt)
    end if

    tmp = sd%ps(0)-sd%ps(krel)
!    if (krel.gt.sd%kinv) then !elevated convection
!       tmp = sd%dp(krel)*0.1
!    end if
    cbmf_max = tmp*dpc%frac_limit_d/sd%delt/Grav

    cbmf_deep = max(min(cbmf_deep, cbmf_max), 0.)
 
    if(cbmf_deep.lt.1.e-10) then 
       cbmf_deep=0.; ocode=10; return
    end if

    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
       cbmf_deep = 0.; ocode=11; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    end if
    
  end subroutine dpconv2
!#####################################################################
!#####################################################################


!#####################################################################
!#####################################################################
  subroutine dpconv3(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode, taudp,  &
                     dcwfndt_fre, dcwfndt_pbl, dcwfndt_dpc, cwfn, cwfn3d, cape3d, &
                     lat, lon, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(inout)  :: sd
    type(adicloud),  intent(inout)  :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode, taudp
    real,            intent(inout)  :: dcwfndt_fre, dcwfndt_pbl, dcwfndt_dpc, cwfn
    real,            intent(inout),  dimension(:)  :: cwfn3d, cape3d
    real,            intent(in)     :: lat, lon
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    integer       :: ksrc, k, k1, k2
    real          :: zcldtop, wrel, cbmf0, tmp, tmp1, avg, fact
    real          :: cape0, cape, cape_new, cwfn0, cwfn_new
    real          :: dcapedt_pbl, dcapedt_fre, dcapedt_tot, dcapedt_dpc
    real          :: dcwfndt_tot, dcwfndt_fre1, dcapedt_dpc_dm, dcwfndt_dpc_dm
    real          :: cbmf_deep1, cbmf_max, scale_fact
    real          :: hltmp, qttmp, thj, qvj, qlj, qij, qse, thvj, nu, exnj
    real          :: latx, lonx, lat1b, lat1e, lon1b, lon1e, lat2b, lat2e, lon2b, lon2e
    real, dimension(size(cwfn3d)) :: temp, qcon

    ier = 0
    ermesg = ' '
    zcldtop = 2000.
    wrel = max(cc%wrel, 0.1)
    taudp=dpc%tau_dp;

   tmp = 1.
   if (tmp.lt.0.) then
    lat1b= 10.; lat1e=18.; lon1b=0.;   lon1e=22.
    lat2b=-20.; lat2e=10.; lon2b=285.; lon2e=305.
    latx=lat*180/3.1415926; lonx=lon*180/3.1415926
    if (sd%land.ge.1) then
      if (latx.gt.lat1b .and. latx.lt.lat1e .and. lonx.gt.lon1b .and. lonx.lt.lon1e) then
        call compute_cwfn3d(dpc,dpn,Uw_p,sd,ac1,cp1,do_coldT,do_ice,rkm_dp,cwfn3d,cape3d,ier,ermesg)
      elseif (latx.gt.lat2b .and. latx.lt.lat2e .and. lonx.gt.lon2b .and. lonx.lt.lon2e) then
        call compute_cwfn3d(dpc,dpn,Uw_p,sd,ac1,cp1,do_coldT,do_ice,rkm_dp,cwfn3d,cape3d,ier,ermesg)
	if (sd%tke.lt.0.5) then
	   tmp=1
	elseif(sd%tke.lt.1.) then	
	   tmp=2
	elseif(sd%tke.lt.1.5) then	
	   tmp=3
        elseif(sd%tke.lt.2.0) then	
	   tmp=4
	elseif(sd%tke.lt.2.5) then	
	   tmp=5
	endif
      endif
    endif
   endif

!create the non-updated sounding (sd1) by removing tendencies from other processes
    call sd_copy_k(sd, sd1)
    do k=1,sd%kmax
       sd1%t (k)=sd1%t (k)-(sd%tdt_rad(k)+sd%tdt_dyn(k)+sd%tdt_dif(k)) * sd%delt
       sd1%qv(k)=sd1%qv(k)-(sd%qdt_dyn(k)+sd%qdt_dif(k)) * sd%delt
       sd1%qv(k)=max(sd1%qv(k),0.)
    end do
    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)

!lauch both adiabatic and entraining plumes (ac1, cp1) based on sd1 to compute cape0 and cwfn0
!will need addtional adiabatic and entraining plumes to evaluate dcwfn and dcape due to PBL tendencies
    call ac_clear_k(ac1);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call adi_cloud_k(sd1%zsrc, sd1%psrc, sd1%hlsrc, sd1%thcsrc, sd1%qctsrc, &
         sd1, Uw_p, .false., do_ice, ac1)
    cape0=ac1%cape

!make base mass flux cbmf0 dependent on cape so that first guess is closer to the actual cbmf_deep
!restrict cbmf0 to be between 0.001 to 0.01, this cbmf0 is used before true cbmf_deep is diagnosed
    cbmf0 = 0.001+(0.01-0.001)*cape0*0.001 
    cbmf0 = min(cbmf0,0.01)
    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn0=0.
    do k = cp1%krel,cp1%let
       cwfn0 = cwfn0 + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)
    end do

!compute free troposphere generation of dcwfndt due to other processes; 
!sd1 contain previous sounding, sd contain updated sounding;
!ct contain shallow_cu tendencies; deep convective depth has been diagnosed in cp1
    dcwfndt_fre1= 0.;
    dcwfndt_fre = 0.; !unit:W/m2,dcwfndt from processes other than deep convection;
    do k=sd%kinv, cp1%let
       tmp = sd%tdt_rad(k)+sd%tdt_dyn(k)+sd%tdt_dif(k)!+ct%tten(k)
       dcwfndt_fre1=dcwfndt_fre1-tmp/sd%t(k)*sd%dp(k)
       hltmp=sd%hl (k)+ct%hlten (k)*sd%delt
       qttmp=sd%qct(k)+ct%qctten(k)*sd%delt
       call findt_k(sd%z(k),sd%p(k),hltmp,qttmp,thj,qvj,qlj,qij,qse,thvj,dpn%do_ice,Uw_p)
       temp(k)=thj*exn_k(sd%p(k),Uw_p); qcon(k)=qlj+qij;
       tmp = (temp(k)-sd1%t(k))/sd%delt
       dcwfndt_fre=dcwfndt_fre-tmp/sd%t(k)*sd%dp(k)
    end do
    dcwfndt_fre=dcwfndt_fre1

! additional plumes using the updated sounding (sd) to compute wfn and cape 
! differences wfn-wfn0 and cape-cape0 give total cape tendencies over the timestep
! adiabatic plume ac which contains source air properties were already computed outside the routine
! convective tendencies are not needed at this time
    call ac_clear_k(ac);
    call cp_clear_k(cp);  cp %maxcldfrac=1.;
    call adi_cloud_k(sd%zsrc,sd%psrc,sd%hlsrc,sd%thcsrc,sd%qctsrc, sd, Uw_p, .false., do_ice, ac)
    call cumulus_plume_k(dpn, sd, ac, cp, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn=0.; !unit: (kg/m3)(m2/s2) or N/m2/, mass weighted, more readily related to mass flux
    do k = cp%krel,cp%let
       cwfn = cwfn + cp%buo(k)/cp%thvtop(k)*cp%dp(k)
    end do
    cape=ac%cape;

!compute dcwfn_dt between updated (sd) and previous sounding (sd1)
    dcapedt_tot=(cape-cape0)/sd%delt !unit:W/m2, dcape_dt
    dcwfndt_tot=(cwfn-cwfn0)/sd%delt !unit:W/m2, dcwfn_dt due to PBL tendencies

!compute dcwfn_dt due to PBL as the difference between total and free
    dcwfndt_pbl=dcwfndt_tot-dcwfndt_fre

!make taudp depend on pbl and free troposphere processes
    if (sd%land.ge.1) then
      if (dcwfndt_fre.gt.0 .and. dcwfndt_pbl.gt.0) then
    	taudp=(dcwfndt_fre*dpc%tau_dp+dcwfndt_pbl*dpc%tau_dp*dpc%tau_dp_fact) &
                 /(dcwfndt_fre+dcwfndt_pbl)
      elseif (dcwfndt_fre.gt.0 .and. dcwfndt_pbl.lt.0) then
	taudp=dpc%tau_dp
      elseif (dcwfndt_fre.lt.0 .and. dcwfndt_pbl.gt.0) then
        taudp=dpc%tau_dp*dpc%tau_dp_fact
      endif
    endif

    if (dpc%do_cgust_dp .and. sd%land.gt.0.5) then
       tmp=sd%cgust/(sd%cgust0+sd%cgust)
       tmp=(tmp-0.5)/sqrt(2.*sd%sigma0*sd%sigma0)
       fact=(1.-0.5*erfccc(tmp))*dpc%tau_dp_fact
       taudp=dpc%tau_dp/fact
    endif

!now integration may exit if deep convection criteria are not satisfied
    call ct_clear_k(ct1);
    if ( cbmf_deep.eq.0 ) then
       ocode=6;
       return
    end if
    if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
       cbmf_deep = 0.; ocode=7; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp, ct1, do_coldT)
    end if

!if deep convection exist, compute dcwfndt_dpc and dcapedt_dpc due to deep convection by updating the
!convective tendencies (in ct1) based on cbmf0 and evaluate new cwfn and cape
    call sd_copy_k(sd, sd1)
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt
    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)
    call ac_clear_k(ac1);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call adi_cloud_k(sd1%zsrc,sd1%psrc,sd1%hlsrc,sd1%thcsrc,sd1%qctsrc,sd1,Uw_p,.false.,do_ice,ac1)
    cape_new = ac1%cape
    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn_new=0.;
    do k = cp1%krel,cp1%let !new convective depth used here
       cwfn_new = cwfn_new + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)
    end do

!compute changes in cape and cwfn due to convective mass flux cbmf0)
!dcapedt_dpc and dcwfndt_dpc have unit N/m2/s or W/m2
    dcapedt_dpc = (ac%cape - ac1%cape)/sd%delt !CAPE destructive rate by convection
    dcwfndt_dpc = (cwfn    - cwfn_new)/sd%delt !CWFN destructive rate by convection

!normalized by cbmf0 so that they are changes due to unit cloud base mass flux
    dcwfndt_dpc_dm = dcwfndt_dpc / cbmf0            !per unit mass flux
    dcapedt_dpc_dm = dcapedt_dpc / cbmf0            !per unit mass flux

!if convection reduce cwfn and cwfn is greater than a threshold, generate deep convection
    if ((cwfn .gt. dpc%cwfn_th) .and. dcwfndt_dpc_dm .gt. 0.) then
       cbmf_deep= (cwfn-dpc%cwfn_th) / taudp / dcwfndt_dpc_dm
    else
       call ct_clear_k(ct1);
       cbmf_deep=0.; ocode=8; return
    end if

!restrict cbmf_deep with a maximum allowable value so that CFL condition satisfied
    cbmf_max  = (sd%ps(0) - sd%ps(cp%krel))*(dpc%frac_limit_d/sd%delt)/Grav
    cbmf_deep = max(min(cbmf_deep, cbmf_max), 0.)

!rescale convective tendencies generated in ct1 due to cbmf0 based on cbmf_deep/cbmf0
!if no precip reevaporation, the system should be very linear and rescaling nearly exactly
!reproduce the actual tendencies if one lauch the plume with cbmf_deep
    scale_fact = cbmf_deep/cbmf0 
    if (scale_fact.lt.1.e-10) scale_fact = 0.
    call scale_ct (ct1, sd, Uw_p, scale_fact)

!compute convective destruction of cape and cwfn over this time-step, rescaling may not be 
!linear, so an additional computation of the actual change would be necessary
    dcapedt_dpc = dcapedt_dpc_dm * cbmf_deep
    dcwfndt_dpc = dcwfndt_dpc_dm * cbmf_deep

!check the actual change in cwfn due to the scaled tendencies in ct1
    call sd_copy_k(sd, sd1)
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt
    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)
    call ac_clear_k(ac1);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call adi_cloud_k(sd1%zsrc,sd1%psrc,sd1%hlsrc,sd1%thcsrc,sd1%qctsrc,sd1,Uw_p,.false.,do_ice,ac1)
    cape_new=ac1%cape
    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn_new=0.;
    do k = cp1%krel,cp1%let
       cwfn_new = cwfn_new + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)
    end do
    dcapedt_dpc = (cape - cape_new)/sd%delt !actual dcape_dt due to convection tendencies
    dcwfndt_dpc = (cwfn - cwfn_new)/sd%delt !actual dcwfn_dt due to convective tendencies

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!check how much differences if we actually lauch the plume with cbmf_deep (to be deleted)
   tmp = 1.
   if (tmp.lt.0.) then
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    call sd_copy_k(sd, sd1)
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt
    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)
    call ac_clear_k(ac1);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    ksrc=sd1%ksrc
    call adi_cloud_k(sd1%zs(ksrc), sd1%ps(ksrc), sd1%hl(ksrc), sd1%thc(ksrc), sd1%qct(ksrc), &
         sd1, Uw_p, .false., do_ice, ac1)
    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn_new=0.;
    do k = cp1%krel,cp1%let
       cwfn_new = cwfn_new + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)
    end do
    dcwfndt_dpc = (cwfn - cwfn_new)/sd%delt
   end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine dpconv3
!#####################################################################
!#####################################################################

!#####################################################################
!#####################################################################
  subroutine dpconv4(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode, taudp,  &
                     dcwfndt_fre, dcwfndt_pbl, dcwfndt_dpc, cwfn, cwfn3d, cape3d, &
                     lat, lon, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(inout)  :: sd
    type(adicloud),  intent(inout)  :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode, taudp
    real,            intent(inout)  :: dcwfndt_fre, dcwfndt_pbl, dcwfndt_dpc, cwfn
    real,            intent(inout),  dimension(:)  :: cwfn3d, cape3d
    real,            intent(in)     :: lat, lon
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    integer       :: ksrc, k, k1, k2
    real          :: zcldtop, wrel, cbmf0, tmp, tmp1, avg, fact
    real          :: cape0, cape, cape_new, cwfn0, cwfn_new
    real          :: dcapedt_pbl, dcapedt_fre, dcapedt_tot, dcapedt_dpc
    real          :: dcwfndt_tot, dcwfndt_fre1, dcapedt_dpc_dm, dcwfndt_dpc_dm
    real          :: cbmf_deep1, cbmf_max, scale_fact
    real          :: hltmp, qttmp, thj, qvj, qlj, qij, qse, thvj, nu, exnj
    real          :: latx, lonx, lat1b, lat1e, lon1b, lon1e, lat2b, lat2e, lon2b, lon2e
    real, dimension(size(cwfn3d)) :: temp, qcon

    ier = 0
    ermesg = ' '
    zcldtop = 2000.
    wrel = max(cc%wrel, 0.1)
    taudp=dpc%tau_dp;

    if (dpc%do_cgust_dp .and. (sd%land.gt.0.5)) then
       if (sd%cgust.gt.sd%cgust0 .and. ac%cape.gt.ac%cin) then
       	  dpn%do_forcedlifting = .true.
          dpn%do_ppen = .false.
	  rkm_dp = dpc%rkm_dp2
          tmp = sd%cgust/(sd%cgust0+sd%cgust)
	  tmp = (1.-sqrt(tmp))
          rkm_dp = rkm_dp * tmp
       endif
    end if

!create the non-updated sounding (sd1) by removing tendencies from other processes
    call sd_copy_k(sd, sd1)
    do k=1,sd%kmax
       sd1%t (k)=sd1%t (k)-(sd%tdt_rad(k)+sd%tdt_dyn(k)+sd%tdt_dif(k)) * sd%delt
       sd1%qv(k)=sd1%qv(k)-(sd%qdt_dyn(k)+sd%qdt_dif(k)) * sd%delt
       sd1%qv(k)=max(sd1%qv(k),0.)
    end do
    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)

!lauch both adiabatic and entraining plumes (ac1, cp1) based on sd1 to compute cape0 and cwfn0
!will need addtional adiabatic and entraining plumes to evaluate dcwfn and dcape due to PBL tendencies
    call ac_clear_k(ac1);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call adi_cloud_k(sd1%zsrc, sd1%psrc, sd1%hlsrc, sd1%thcsrc, sd1%qctsrc, &
         sd1, Uw_p, .false., do_ice, ac1)
    cape0=ac1%cape

!make base mass flux cbmf0 dependent on cape so that first guess is closer to the actual cbmf_deep
!restrict cbmf0 to be between 0.001 to 0.01, this cbmf0 is used before true cbmf_deep is diagnosed
    cbmf0 = 0.001+(0.01-0.001)*cape0*0.001 
    cbmf0 = min(cbmf0,0.01)
    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn0=0.
    do k = cp1%krel,cp1%let
       cwfn0 = cwfn0 + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)
    end do

!compute free troposphere generation of dcwfndt due to other processes; 
!sd1 contain previous sounding, sd contain updated sounding;
!ct contain shallow_cu tendencies; deep convective depth has been diagnosed in cp1
    dcwfndt_fre = 0.; !unit:W/m2,dcwfndt from processes other than deep convection;
    do k=sd%kinv, cp1%let
       tmp = sd%tdt_rad(k)+sd%tdt_dyn(k)+sd%tdt_dif(k) !+ct%tten(k)
       dcwfndt_fre=dcwfndt_fre-tmp/sd%t(k)*sd%dp(k)
    end do

! additional plumes using the updated sounding (sd) to compute wfn and cape 
! differences wfn-wfn0 and cape-cape0 give total cape tendencies over the timestep
! adiabatic plume ac which contains source air properties were already computed outside the routine
! convective tendencies are not needed at this time
    call ac_clear_k(ac);
    call cp_clear_k(cp);  cp %maxcldfrac=1.;
    call adi_cloud_k(sd%zsrc,sd%psrc,sd%hlsrc,sd%thcsrc,sd%qctsrc, sd, Uw_p, .false., do_ice, ac)
    call cumulus_plume_k(dpn, sd, ac, cp, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn=0.; !unit: (kg/m3)(m2/s2) or N/m2/, mass weighted, more readily related to mass flux
    do k = cp%krel,cp%let
       cwfn = cwfn + cp%buo(k)/cp%thvtop(k)*cp%dp(k)
    end do
    cape=ac%cape;

!compute dcwfn_dt between updated (sd) and previous sounding (sd1)
    dcapedt_tot=(cape-cape0)/sd%delt !unit:W/m2, dcape_dt
    dcwfndt_tot=(cwfn-cwfn0)/sd%delt !unit:W/m2, dcwfn_dt due to PBL tendencies

!compute dcwfn_dt due to PBL as the difference between total and free
    dcwfndt_pbl=dcwfndt_tot-dcwfndt_fre

!make taudp depend on pbl and free troposphere processes
    if (sd%land.ge.0.5) then
      if (dcwfndt_fre.gt.0 .and. dcwfndt_pbl.gt.0) then
    	taudp=(dcwfndt_fre*dpc%tau_dp+dcwfndt_pbl*dpc%tau_dp*dpc%tau_dp_fact) &
                 /(dcwfndt_fre+dcwfndt_pbl)
      elseif (dcwfndt_fre.gt.0 .and. dcwfndt_pbl.le.0) then
	taudp=dpc%tau_dp
      elseif (dcwfndt_fre.le.0 .and. dcwfndt_pbl.gt.0) then
        taudp=dpc%tau_dp*dpc%tau_dp_fact
      endif
    endif

!now integration may exit if deep convection criteria are not satisfied
    call ct_clear_k(ct1);
    if ( cbmf_deep.eq.0 ) then
       ocode=6; cp1%cush=-1;
       return
    end if
    if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
       cbmf_deep = 0.; ocode=7; cp1%cush=-1; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp, ct1, do_coldT)
    end if

!if deep convection exist, compute dcwfndt_dpc and dcapedt_dpc due to deep convection by updating the
!convective tendencies (in ct1) based on cbmf0 and evaluate new cwfn and cape
    call sd_copy_k(sd, sd1)
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt
    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)
    call ac_clear_k(ac1);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call adi_cloud_k(sd1%zsrc,sd1%psrc,sd1%hlsrc,sd1%thcsrc,sd1%qctsrc,sd1,Uw_p,.false.,do_ice,ac1)
    cape_new = ac1%cape
    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn_new=0.;
    do k = cp1%krel,cp1%let !new convective depth used here
       cwfn_new = cwfn_new + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)
    end do

!compute changes in cape and cwfn due to convective mass flux cbmf0)
!dcapedt_dpc and dcwfndt_dpc have unit N/m2/s or W/m2
    dcapedt_dpc = (ac%cape - ac1%cape)/sd%delt !CAPE destructive rate by convection
    dcwfndt_dpc = (cwfn    - cwfn_new)/sd%delt !CWFN destructive rate by convection

!normalized by cbmf0 so that they are changes due to unit cloud base mass flux
    dcwfndt_dpc_dm = dcwfndt_dpc / cbmf0            !per unit mass flux
    dcapedt_dpc_dm = dcapedt_dpc / cbmf0            !per unit mass flux

!if convection reduce cwfn and cwfn is greater than a threshold, generate deep convection
    if ((cwfn .gt. dpc%cwfn_th) .and. dcwfndt_dpc_dm .gt. 0.) then
       cbmf_deep= (cwfn-dpc%cwfn_th) / taudp / dcwfndt_dpc_dm
    else
       call ct_clear_k(ct1);
       cbmf_deep=0.; ocode=8; cp1%cush=-1; return
    end if

!restrict cbmf_deep with a maximum allowable value so that CFL condition satisfied
    cbmf_max  = (sd%ps(0) - sd%ps(cp%krel))*(dpc%frac_limit_d/sd%delt)/Grav
    cbmf_deep = max(min(cbmf_deep, cbmf_max), 0.)

!rescale convective tendencies generated in ct1 due to cbmf0 based on cbmf_deep/cbmf0
!if no precip reevaporation, the system should be very linear and rescaling nearly exactly
!reproduce the actual tendencies if one lauch the plume with cbmf_deep
    scale_fact = cbmf_deep/cbmf0 
    if (scale_fact.lt.1.e-10) scale_fact = 0.
    call scale_ct (ct1, sd, Uw_p, scale_fact)

  end subroutine dpconv4

!#####################################################################
!#####################################################################

  subroutine scale_ct (ct, sd, Uw_p, scale)
    implicit none
    type(ctend),     intent(inout) :: ct
    type(sounding),  intent(in)    :: sd
    type(uw_params), intent(inout) :: Uw_p
    real,            intent(in)    :: scale
    integer  :: i, k

    ct%uten   =ct%uten  *scale;
    ct%vten   =ct%vten  *scale;
    ct%tten   =ct%tten  *scale;
    ct%qvten  =ct%qvten *scale;
    ct%qlten  =ct%qlten *scale;
    ct%qiten  =ct%qiten *scale;
    ct%qaten  =ct%qaten *scale;
    ct%qnten  =ct%qnten *scale;
    ct%qldet  =ct%qldet *scale;
    ct%qidet  =ct%qidet *scale;
    ct%qadet  =ct%qadet *scale;
    ct%qndet  =ct%qndet *scale;
    ct%hlten  =ct%hlten *scale;
    ct%thcten =ct%thcten*scale;
    ct%qctten =ct%qctten*scale;
    ct%qvdiv  =ct%qvdiv *scale;
    ct%qldiv  =ct%qldiv *scale;
    ct%qidiv  =ct%qidiv *scale;
    ct%hlflx  =ct%hlflx *scale;
    ct%thcflx =ct%thcflx*scale;
    ct%qctflx =ct%qctflx*scale;
    ct%qvflx  =ct%qvflx *scale;
    ct%qlflx  =ct%qlflx *scale;
    ct%qiflx  =ct%qiflx *scale;
    ct%qaflx  =ct%qaflx *scale;
    ct%qnflx  =ct%qnflx *scale;
    ct%umflx  =ct%umflx *scale;
    ct%vmflx  =ct%vmflx *scale;
    ct%pflx   =ct%pflx  *scale;
    ct%tevap  =ct%tevap *scale;
    ct%qevap  =ct%qevap *scale;
    ct%trflx  =ct%trflx *scale; 
    ct%trten  =ct%trten *scale;
    ct%trwet  =ct%trwet *scale;
    ct%trevp  =ct%trevp *scale;
    ct%dtring =ct%dtring*scale;
    ct%snow   =ct%snow  *scale;
    ct%rain   =ct%rain  *scale;
    ct%denth  =ct%denth *scale;
    ct%dting  =ct%dting *scale;
    ct%dqtmp  =ct%dqtmp *scale;

    ct%dting=0.;
    ct%denth=0.; ct%dqtmp=0.; ct%uav=0.;ct%vav=0.; 
    do k = 1,sd%kmax
       ct%denth = ct%denth + (Uw_p%cp_air*ct%tten(k)-Uw_p%HLv*ct%qlten(k)-Uw_p%HLs*ct%qiten(k))*sd%dp(k)/Uw_p%grav
       ct%dqtmp = ct%dqtmp + (ct%qvten(k) + ct%qlten(k) + ct%qiten(k))*sd%dp(k)/Uw_p%grav
       ct%uav   = ct%uav   + ct%uten(k)*sd%dp(k)/Uw_p%grav
       ct%vav   = ct%vav   + ct%vten(k)*sd%dp(k)/Uw_p%grav
       ct%dting = ct%dting + Uw_p%cp_air*ct%tten(k)*sd%dp(k)/Uw_p%grav
    end do
    ct%denth = ct%denth - Uw_p%HLv*ct%rain - Uw_p%HLs*ct%snow
    ct%dqtmp = ct%dqtmp + ct%rain + ct%snow
 
    do i = 1,size(sd%tr,2)
       ct%dtring(i)=0.;
       do k = 1,sd%kmax
         ct%dtring(i) = ct%dtring(i) + Uw_p%cp_air*ct%trten(k,i)*sd%dp(k)/Uw_p%grav
       enddo
     end do

  end subroutine scale_ct

!#####################################################################
!#####################################################################
  subroutine dpconv5(dpc, dpn, Uw_p, sd, ac, cc, cp, ct, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, sd1, ac1, cp1, ct1, ocode, &
                     taudp, dcwfn_con, dcwfn_ndp, dcwfn_pbl, pwfn, dwfn, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(in)     :: ac
    type(cclosure),  intent(in)     :: cc
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    type(cplume),    intent(inout)  :: cp
    type(ctend),     intent(inout)  :: ct
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(sounding),  intent(inout)  :: sd1
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(inout)  :: ocode, taudp, dcwfn_con, dcwfn_ndp, dcwfn_pbl, pwfn, dwfn
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    integer       :: ksrc, k, k1, k2
    real          :: zcldtop, wrel, cbmf0, tmp
    real          :: cape, cwfn, cape_pbl, cwfn_pbl, cape_new, cwfn_new, pwfn_new
    real          :: dcape_pbl, dcwfn_con_dm, dcape_con, dcape_con_dm, dcwfn_ndp1
    real          :: cbmf_deep1, cbmf_max, scale_fact
    real          :: hltmp, qttmp, thj, qvj, qlj, qij, qse, thvj, nu, exnj, temp

    ier = 0
    ermesg = ' '
    zcldtop = 2000 !sd%z(cp%ltop)
    wrel = max(cc%wrel, 0.1)
    taudp=dpc%tau_dp;

!creat a new sounding with update of only PBL tendencies from other processes
    call sd_copy_k(sd, sd1)
    k=sd%kinv-1;
    sd1%t (1:k)=sd1%t (1:k)+(sd%tdt_rad(1:k)+sd%tdt_dyn(1:k)+sd%tdt_dif(1:k)+ct%tten(1:k)) * sd%delt
    sd1%qv(1:k)=sd1%qv(1:k)+(sd%qdt_dyn(1:k)+sd%qdt_dif(1:k)+ct%qvten(1:k)) * sd%delt
    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)

!lauch both adiabatic and entraining plumes to compute cape and cwfn 
!will need addtional adiabatic and entraining plumes to evaluate dcwfn and dcape due to PBL tendencies
    call ac_clear_k(ac1);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    ksrc=sd%ksrc
    call adi_cloud_k(sd1%zs(ksrc), sd1%ps(ksrc), sd1%hl(ksrc), sd1%thc(ksrc), sd1%qct(ksrc), &
         sd1, Uw_p, .false., do_ice, ac1)
    cape_pbl=ac1%cape

!make base mass flux cbmf0 dependent on cape so that first guess is closer to the actual cbmf_deep
!restrict cbmf0 to be between 0.001 to 0.01, this cbmf0 is used before true cbmf_deep is diagnosed
    cbmf0 = 0.001+(0.01-0.001)*ac%cape*0.001 
    cbmf0 = min(cbmf0,0.01)
    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn_pbl=0.
    do k = cp1%krel,cp1%let
       cwfn_pbl = cwfn_pbl + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)
    end do

! additional plumes using non-updated sd to compute wfn and cape 
! differences wfn_pbl-wfn and cape_pbl-cape are changes due to PBL tendencies over the timestep
! note adiabatic plume ac which contain source air properties were already computed outside the routine
! convective tendencies are not needed for this evaluation
    call cp_clear_k(cp);  cp %maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn=0.; !unit: (kg/m3)(m2/s2) or N/m2/, mass weighted, more readily related to mass flux
    do k = cp%krel,cp%let
       cwfn = cwfn + cp%buo(k)/cp%thvtop(k)*cp%dp(k)
    end do
    cape=ac%cape;

!now compute dcwfn and dcape due to PBL tendencies
    dcape_pbl=(cape_pbl-cape)/sd%delt !unit:W/m2, dcape_dt due to PBL tendencies
    dcwfn_pbl=(cwfn_pbl-cwfn)/sd%delt !unit:W/m2, dcwfn_dt due to PBL tendencies

!compute dcwfn due to other processes; ct contain shallow_cu tendencies at this time
!at this time, deep convective depth at current state has been diagnosed in cp
    dcwfn_ndp1= 0.;
    dcwfn_ndp = 0.; !unit:W/m2,dcwfn_dt from processes other than deep convection;
    do k=sd%kinv, cp%let
       tmp = sd%tdt_rad(k)+sd%tdt_dyn(k)+sd%tdt_dif(k)+ct%tten(k)
       dcwfn_ndp1=dcwfn_ndp1-tmp/sd%t(k)*sd%dp(k)

       hltmp=sd%hl (k)+(sd%tdt_dyn(k)+sd%tdt_dif(k)+ct%tten(k)+sd%tdt_rad(k))*sd%delt*Uw_p%cp_air
       qttmp=sd%qct(k)+(sd%qdt_dyn(k)+sd%qdt_dif(k)+ct%qvten(k))*sd%delt
       call findt_k (sd%z(k),sd%p(k),hltmp,qttmp,thj,qvj,qlj,qij,qse,thvj,dpn%do_ice,Uw_p)
       temp=thj*exn_k(sd%p(k),Uw_p)
       temp=(temp-sd%t(k))/sd%delt
       dcwfn_ndp=dcwfn_ndp-temp/sd%t(k)*sd%dp(k)
    end do

!prognose pwfn due to processes other than deep convection
!integration should not exit before this point in any case (whether or not deep convection exist)
!if no deep convection supported and no positive tendency for cwfn reset pwfn 
    pwfn_new = pwfn + (dcwfn_ndp + dcwfn_pbl)*sd%delt
    if ((cwfn .le. 0.) .and. (dcwfn_ndp + dcwfn_pbl).le.0) then
       pwfn_new=0.0
       dcwfn_ndp=0.0
       dcwfn_pbl=0.0
    end if
    pwfn = pwfn_new

    if (pwfn_new .gt. 10000.) then
       tmp=pwfn_new-pwfn
    end if
    if (cape .gt. 1000.) then
       tmp=pwfn_new-pwfn
    end if


!make taudp depend on pbl and free troposphere processes
    if (sd%land.eq.1) then
       taudp = (dcwfn_ndp*dpc%tau_dp + dcwfn_pbl*dpc%tau_dp*dpc%tau_dp_fact)/(dcwfn_ndp+dcwfn_pbl)
    end if

!integration now may exit if deep convection criteria are not satisfied
    if ( cbmf_deep.eq.0 ) then
       call ct_clear_k(ct1);
       ocode=6;
       return
    end if

!Shoud we use updated sounding to compute convective destruction of wfn? currently not!
    if(cp%ltop.lt.cp%krel+2 .or. cp%let.le.cp%krel+1) then
       call ct_clear_k(ct1);
       cbmf_deep = 0.; ocode=7; return
    else
       call cumulus_tend_k(dpn, sd, Uw_p, cp, ct1, do_coldT)
    end if

!if deep convection exist, compute dcwfn and dcape due to deep convection by updating the
!convective tendencies based on cbmf0 (contained in ct1) and evaluate new cwfn and cape
    call sd_copy_k(sd, sd1)
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt
    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)
    call ac_clear_k(ac1);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    ksrc=sd1%ksrc
    call adi_cloud_k(sd1%zs(ksrc), sd1%ps(ksrc), sd1%hl(ksrc), sd1%thc(ksrc), sd1%qct(ksrc), &
         sd1, Uw_p, .false., do_ice, ac1)
    cape_new = ac1%cape
    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn_new=0.;
    do k = cp1%krel,cp1%let !new convective depth used here
       cwfn_new = cwfn_new + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)
    end do

!compute changes in cape and cwfn due to convective mass flux cbmf0)
!dcape_con and dcwfn_con have unit N/m2/s or W/m2
    dcape_con    = (ac%cape - ac1%cape)/sd%delt !CAPE destructive rate by convection
    dcwfn_con    = (cwfn    - cwfn_new)/sd%delt !CWFN destructive rate by convection

!normalized by cbmf0 so that they are changes due to unit cloud base mass flux
    dcwfn_con_dm = dcwfn_con / cbmf0            !per unit mass flux
    dcape_con_dm = dcape_con / cbmf0            !per unit mass flux

!if convection reduce cwfn and cwfn is greater than a threshold, generate deep convection
    if ((cwfn .gt. dpc%cwfn_th) .and. dcwfn_con_dm .gt. 0.) then
       cbmf_deep= (cwfn-dpc%cwfn_th) / taudp / dcwfn_con_dm
    else
       call ct_clear_k(ct1);
       cbmf_deep=0.; ocode=8; return
    end if

    if ((pwfn_new .gt. dpc%cwfn_th) .and. dcwfn_con_dm .gt. 0.) then
       cbmf_deep1= (pwfn_new-dpc%cwfn_th) / taudp / dcwfn_con_dm
       !cbmf_deep1= dcwfn_ndp / dcwfn_con_dm
    else
       call ct_clear_k(ct1);
       cbmf_deep1=0.; ocode=8; return
    end if
    cbmf_deep=cbmf_deep1;

!restrict cbmf_deep with a maximum allowable value so that CFL-like condition satisfied
    cbmf_max  = (sd%ps(0) - sd%ps(cp%krel))*(dpc%frac_limit_d/sd%delt)/Grav
    cbmf_deep = max(min(cbmf_deep, cbmf_max), 0.)

!rescale convective tendencies generated in ct1 due to cbmf0 based on cbmf_deep/cbmf0
!if no precip reevaporation, the system should be very linear and rescaling nearly exactly
!reproduce the actual tendencies if one lauch the plume with cbmf_deep
    scale_fact = cbmf_deep/cbmf0 
    if (scale_fact.lt.1.e-10) scale_fact = 0.
    call scale_ct (ct1, sd, Uw_p, scale_fact)

!compute convective destruction of cape and cwfn over this time-step, rescaling may not be 
!linear, so an additional computation of the actual change would be necessary
    dcape_con = dcape_con_dm * cbmf_deep
    dcwfn_con = dcwfn_con_dm * cbmf_deep

!check the actual change in cwfn due to the scaled tendencies in ct1
    call sd_copy_k(sd, sd1)
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt
    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)
    call ac_clear_k(ac1);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    ksrc=sd1%ksrc
    call adi_cloud_k(sd1%zs(ksrc), sd1%ps(ksrc), sd1%hl(ksrc), sd1%thc(ksrc), sd1%qct(ksrc), &
         sd1, Uw_p, .false., do_ice, ac1)
    cape_new=ac1%cape
    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn_new=0.;
    do k = cp1%krel,cp1%let
       cwfn_new = cwfn_new + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)
    end do
    dcape_con = (cape - cape_new)/sd%delt !actual dcape_dt due to convection tendencies
    dcwfn_con = (cwfn - cwfn_new)/sd%delt !actual dcwfn_dt due to convective tendencies
    dcwfn_con_dm = dcwfn_con / cbmf_deep  !actual value per unit mass flux
    dcape_con_dm = dcape_con / cbmf_deep  !actual value per unit mass flux

!update prognostic cwfn with dcwfn
    dwfn = cwfn
    pwfn = pwfn_new - dcwfn_con*sd%delt
!    pwfn = max(pwfn,0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!check how much differences if we actually lauch the plume with cbmf_deep (to be deleted)
   tmp = 1.
   if (tmp.lt.0.) then
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    call ct_clear_k(ct1);
    call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
    call sd_copy_k(sd, sd1)
    sd1 % t  = sd1 % t  + ct1%tten  * sd%delt
    sd1 % qv = sd1 % qv + ct1%qvten * sd%delt
    sd1 % ql = sd1 % ql + ct1%qlten * sd%delt
    sd1 % qi = sd1 % qi + ct1%qiten * sd%delt
    sd1 % qa = sd1 % qa + ct1%qaten * sd%delt
    sd1 % qn = sd1 % qn + ct1%qnten * sd%delt
    sd1 % u  = sd1 % u  + ct1%uten  * sd%delt
    sd1 % v  = sd1 % v  + ct1%vten  * sd%delt
    call extend_sd_k(sd1,sd%pblht, do_ice, Uw_p)
    call ac_clear_k(ac1);
    call cp_clear_k(cp1); cp1%maxcldfrac=1.;
    ksrc=sd1%ksrc
    call adi_cloud_k(sd1%zs(ksrc), sd1%ps(ksrc), sd1%hl(ksrc), sd1%thc(ksrc), sd1%qct(ksrc), &
         sd1, Uw_p, .false., do_ice, ac1)
    call cumulus_plume_k(dpn, sd1, ac1, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    cwfn_new=0.;
    do k = cp1%krel,cp1%let
       cwfn_new = cwfn_new + cp1%buo(k)/cp1%thvtop(k)*cp1%dp(k)
    end do
    dcwfn_con = (cwfn - cwfn_new)/sd%delt
   end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine dpconv5
!#####################################################################
!#####################################################################


!#####################################################################
!#####################################################################
  subroutine compute_cwfn3d(dpc, dpn, Uw_p, sd, ac1, cp1, do_coldT, do_ice, &
                     rkm_dp, cwfn3d, cape3d, ier, ermesg)
    implicit none

    type(deepc),     intent(in)     :: dpc
    type(cpnlist),   intent(in)     :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(inout)  :: sd
    type(adicloud),  intent(inout)  :: ac1
    type(cplume),    intent(inout)  :: cp1
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    real,            intent(in)     :: rkm_dp
    real,            intent(inout),  dimension(:)  :: cwfn3d, cape3d
    integer,            intent(inout) :: ier
    character(len=256), intent(inout) :: ermesg

    integer       :: ksrc, i, k, km
    real          :: zcldtop, wrel, cbmf0, tmp, zsrc, psrc, hlsrc, thcsrc, qctsrc,cwfnmax

    ier = 0
    ermesg = ' '
    zcldtop = 2000 
    wrel = 0.1
    cbmf0=0.001
    km=15
    cape3d(:)=0.; cwfn3d(:)=0.;
    do k=1,km
       call ac_clear_k(ac1);
       call cp_clear_k(cp1); cp1%maxcldfrac=1.;
       ksrc=k
       zsrc  =sd%zs (ksrc)
       psrc  =sd%ps (ksrc)
       hlsrc =sd%hl (ksrc)
       thcsrc=sd%thc(ksrc)
       qctsrc=sd%qct(ksrc)
       call adi_cloud_k(zsrc, psrc, hlsrc, thcsrc, qctsrc, sd, Uw_p, .false., do_ice, ac1)
       cape3d(k)=ac1%cape
       call cumulus_plume_k(dpn, sd, ac1, cp1, rkm_dp, cbmf0, wrel, zcldtop, Uw_p, ier, ermesg)
       cwfn3d(k)=0.
       do i = cp1%krel,cp1%let
       	  cwfn3d(k) = cwfn3d(k) + cp1%buo(i)/cp1%thvtop(i)*cp1%dp(i)
       end do
    end do

!    ksrc=1
!    cwfnmax=cwfn3d(1)
!    do k=1,km
!     if (cwfn3d(k) .gt. cwfnmax) then
!       cwfnmax=cwfn3d(k)
!       ksrc=k
!     end if
!    end do
!    sd%ksrc  =ksrc;

  end subroutine compute_cwfn3d
!#####################################################################
!#####################################################################

  subroutine conv_forced(dpc, dpn, Uw_p, sd, ac, do_coldT, do_ice, &
                     rkm_dp, cbmf_deep, cp1, ct1, lat, lon, ier, ermesg)
    implicit none

    type(deepc),     intent(inout)  :: dpc
    type(cpnlist),   intent(inout)  :: dpn
    type(uw_params), intent(inout)  :: Uw_p
    type(sounding),  intent(in)     :: sd
    type(adicloud),  intent(in)     :: ac
    logical,         intent(in)     :: do_coldT
    logical,         intent(in)     :: do_ice
    real,            intent(inout)  :: rkm_dp, cbmf_deep
    type(cplume),    intent(inout)  :: cp1
    type(ctend),     intent(inout)  :: ct1
    real,            intent(in)     :: lat, lon
    integer,            intent(out) :: ier
    character(len=256), intent(out) :: ermesg

    real          :: latx, lonx, lat1b, lat1e, lon1b, lon1e, lat2b, lat2e, lon2b, lon2e
    real          :: zcldtop, wrel, ufrc, tmp
    integer       :: k

    ier = 0
    ermesg = ' '
    zcldtop = 2000

    if (sd%land.gt.0.5) then
      if (ac%cape.gt.dpc%cape_th .and. sd%cgust.gt.sd%cgust0) then
	 lat1b= 30.; lat1e=40.; lon1b=260.; lon1e=270.
      	 lat2b=-20.; lat2e=10.; lon2b=285.; lon2e=305.
      	 latx=lat*180/3.1415926; lonx=lon*180/3.1415926
      	 if (latx.gt.lat1b .and. latx.lt.lat1e .and. lonx.gt.lon1b .and. lonx.lt.lon1e) then
            tmp=1
      	 elseif (latx.gt.lat2b .and. latx.lt.lat2e .and. lonx.gt.lon2b .and. lonx.lt.lon2e) then
            tmp=2
      	 endif
       	 dpn%do_forcedlifting = .true.
         dpn%do_ppen = .false. !dpn%do_pevap = .false.
         call cclosure_gust(sd%cgust, sd, Uw_p, ac, dpc%wcrit_min_gust, cbmf_deep, wrel, ufrc)
	 if(cbmf_deep.lt.1.e-10) then 
	    cbmf_deep=0.;
       	    call cp_clear_k(cp1); cp1%maxcldfrac=1.; cp1%cush=-1;
       	    call ct_clear_k(ct1);
       	    return
    	 end if
    	 call cumulus_plume_k(dpn, sd, ac, cp1, rkm_dp, cbmf_deep, wrel, zcldtop, Uw_p, ier, ermesg)
    	 if(cp1%ltop.lt.cp1%krel+2 .or. cp1%let.le.cp1%krel+1) then
     	    cbmf_deep = 0.; 
       	    call cp_clear_k(cp1); cp1%maxcldfrac=1.; cp1%cush=-1;
       	    call ct_clear_k(ct1);
       	    return
    	 else
	    call cumulus_tend_k(dpn, sd, Uw_p, cp1, ct1, do_coldT)
            return
    	 end if
      endif
    endif
  end subroutine conv_forced
!#####################################################################
!#####################################################################


!#####################################################################
!#####################################################################
  subroutine cclosure_gust(cgust, sd, Uw_p, ac, wcrit_min, cbmf, wrel, ufrc)
  
    implicit none
    real,            intent(in)    :: cgust
    type(sounding),  intent(in)    :: sd
    type(uw_params), intent(in)    :: Uw_p
    type(adicloud),  intent(in)    :: ac
    real,            intent(in)    :: wcrit_min !wcrit_min=0.2
    real,            intent(out)   :: cbmf, wrel, ufrc
    
    integer :: kkk
    real    :: sigmaw, wcrit, erfarg, wexp, wtw, cbmf0, tmp

    cbmf=0.; wrel=0.; ufrc=0.;
    wcrit  = sqrt(2. * ac % cin)
    sigmaw = sqrt(cgust)
    wcrit = max(wcrit, wcrit_min*sigmaw)
    cbmf   = ac % rho0lcl * sigmaw / 2.5066 * exp(-0.5*((wcrit/sigmaw)**2.))

    !Diagnose updraft fraction sqrt(2.) = 1.4142
    erfarg=wcrit / (1.4142 * sigmaw)
    if(erfarg.lt.20.)then
       ufrc = min(0.15, 0.5*erfccc(erfarg))
    else
       ufrc = 0.
    endif

    if(ufrc.gt.0.0) then !Diagnose expected value of cloud base vertical velocity
       wexp = cbmf / ac % rho0lcl / ufrc
    else
       wexp = 0.
       cbmf = 0.
    endif

    wtw = wexp * wexp - 2 * ac % cin 
    if(wtw.le.0.) then
       wrel=0.; 
    else
       wrel=sqrt(wtw)
    end if
    wrel=min(wrel, 25.)

    cbmf0 = sd%dp(sd%ksrc) / sd%delt / Uw_p%GRAV
    kkk = max(max(ac%klcl, sd%kinv),1)
    tmp = (sd%ps(0) - sd%ps(kkk)) * 0.25 
    cbmf0 = tmp / sd%delt / Uw_p%GRAV

    cbmf = min (cbmf, cbmf0)
    if (wrel .gt. 0.) then
      ufrc=cbmf / wrel /ac % rho0lcl
    else
      ufrc=0.
    end if

    return

  end subroutine cclosure_gust
!#####################################################################
!#####################################################################


end MODULE DEEP_CONV_MOD
