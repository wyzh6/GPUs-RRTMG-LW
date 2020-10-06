
       module rrtmg_lw_rad


! -------- Modules --------
      use spmd_utils, only:iam
      use shr_kind_mod, only: r8 => shr_kind_r8
      use ppgrid,       only: pcols, begchunk, endchunk
      use perf_mod,          only: t_startf, t_stopf

      use rrlw_vsn
      use mcica_subcol_gen_lw, only: mcica_subcol_lw
      !use rrtmg_lw_cldprmc, only: cldprmc
      use rrtmg_lw_setcoef
      use rrlw_vsn, only: hvrclc, hnamclc,hvrset, hnamset,hvrtau, hnamtau,hvrrtc, hnamrtc
      !----------------------------------

      use parrrtm, only : nbndlw, mg, maxxsec, mxmol,ngptlw
     
      use rrlw_ref
      !----------------------------------

!      use parkind, only : im => kind_im, rb => kind_r8

      use rrlw_con, only: fluxfac, heatfac

  

      use rrlw_tbl, only: ntbl,tblint, bpade, tau_tbl, exp_tbl, tfn_tbl


      use parrrtm, only : ng1
      use rrlw_kg01, only : fracrefa01 => fracrefa,fracrefb01 => fracrefb,absa01 => absa,absb01 => absb, &
                        ka_mn201 => ka_mn2,kb_mn201 => kb_mn2,selfref01 => selfref,forref01 => forref

      use parrrtm, only : ng2, ngs1
      use rrlw_kg02, only : fracrefa02 => fracrefa,fracrefb02 => fracrefb,absa02 => absa,absb02 => absb, &
                          selfref02 => selfref,forref02 => forref

      use parrrtm, only : ng3, ngs2
      use rrlw_ref, only : chi_mls
      use rrlw_kg03, only : fracrefa03 => fracrefa,fracrefb03 => fracrefb,absa03 => absa,absb03 => absb,&
                          ka_mn2o03 => ka_mn2o,kb_mn2o03 => kb_mn2o,selfref03 => selfref,forref03 => forref

      use parrrtm, only : ng4, ngs3
      use rrlw_kg04, only : fracrefa04 => fracrefa,fracrefb04 => fracrefb,absa04 => absa,absb04 => absb, &
                           selfref04 => selfref,forref04 => forref

      use parrrtm, only : ng5, ngs4
      use rrlw_kg05, only : fracrefa05 => fracrefa,fracrefb05 => fracrefb,absa05 => absa,absb05 => absb, &
                          ka_mo305 => ka_mo3,selfref05 => selfref,forref05 => forref,ccl405 => ccl4

      use parrrtm, only : ng6, ngs5
      use rrlw_kg06, only : fracrefa06 => fracrefa,absa06 => absa,ka_mco206 => ka_mco2, &
                           selfref06 => selfref,forref06 => forref,cfc11adj06 => cfc11adj,cfc1206 => cfc12

      use parrrtm, only : ng7, ngs6
      use rrlw_kg07, only : fracrefa07 => fracrefa,fracrefb07 => fracrefb,absa07 => absa,absb07 => absb, &
                          ka_mco207 => ka_mco2,kb_mco207 => kb_mco2,selfref07 => selfref,forref07 => forref

      use parrrtm, only : ng8, ngs7
      use rrlw_kg08, only : fracrefa08 => fracrefa,fracrefb08 => fracrefb,absa08 => absa,absb08 => absb, &
                           ka_mco208 => ka_mco2,ka_mn2o08 => ka_mn2o,ka_mo308 => ka_mo3,kb_mco208 => kb_mco2,kb_mn2o08 => kb_mn2o, &
                          selfref08 => selfref,forref08 => forref,cfc1208 => cfc12,cfc22adj08 => cfc22adj

      use parrrtm, only : ng9, ngs8
      use rrlw_kg09, only : fracrefa09 => fracrefa,fracrefb09 => fracrefb,absa09 => absa,absb09 => absb, &
                          ka_mn2o09 => ka_mn2o,kb_mn2o09 => kb_mn2o,selfref09 => selfref,forref09 => forref

      use parrrtm, only : ng10, ngs9
      use rrlw_kg10, only : fracrefa10 => fracrefa,fracrefb10 => fracrefb,absa10 => absa,absb10 => absb, &
                           selfref10 => selfref,forref10 => forref

      use parrrtm, only : ng11, ngs10
      use rrlw_kg11, only : fracrefa11 => fracrefa,fracrefb11 => fracrefb,absa11 => absa,absb11 => absb, &
                          ka_mo211 => ka_mo2,kb_mo211 => kb_mo2,selfref11 => selfref,forref11 => forref

      use parrrtm, only : ng12, ngs11
      use rrlw_kg12, only : fracrefa12 => fracrefa,absa12 => absa, &
                           selfref12 => selfref,forref12 => forref

      use parrrtm, only : ng13, ngs12
      use rrlw_kg13, only : fracrefa13 => fracrefa,fracrefb13 => fracrefb,absa13 => absa, &
                           ka_mco213 => ka_mco2,ka_mco13 => ka_mco,kb_mo313 => kb_mo3,selfref13 => selfref,forref13 => forref

      use parrrtm, only : ng14, ngs13
      use rrlw_kg14, only : fracrefa14 => fracrefa,fracrefb14 => fracrefb,absa14 => absa,absb14 => absb, &
                           selfref14 => selfref,forref14 => forref

      use parrrtm, only : ng15, ngs14
      use rrlw_kg15, only : fracrefa15 => fracrefa,absa15 => absa, &
                        ka_mn215 => ka_mn2,selfref15 => selfref,forref15 => forref

      use parrrtm, only : ng16, ngs15
      use rrlw_kg16, only : fracrefa16 => fracrefa,fracrefb16 => fracrefb,absa16 => absa,absb16 => absb, &
                           selfref16 => selfref,forref16 => forref




      use cudafor
      use kernel
      

        public :: rrtmg_lw,rtrnmc
!---------------------------------------------------------------
      contains

               subroutine rrtmg_lw &
            (lchnk   ,ncol    ,nlay    ,icld    ,                   &
             play    ,plev    ,tlay    ,tlev    ,tsfc    ,h2ovmr  , &
             o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  ,&
             cfc11vmr,cfc12vmr, &
             cfc22vmr,ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
             cldfmcl ,taucmcl ,ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc)

! -------- Description --------

! --------- Modules ----------

      use parrrtm, only : nbndlw, ngptlw, maxxsec, mxmol,nmol,maxinpx
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi,grav,avogad
      use rrlw_wvn, only: ng, ngb, wavenum1, wavenum2, delwave
      use rrlw_cld, only: abscld1, absliq0, absliq1, &
                          absice0, absice1, absice2, absice3

       use rrlw_wvn, only: totplnk, totplk16,ngb,delwave, ngs,ixindx,nspa, nspb



        implicit none
! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: lchnk                      ! chunk identifier
      integer, intent(in) :: ncol                       ! Number of horizontal columns
      integer, intent(in) :: nlay                       ! Number of model layers
      integer, intent(inout) :: icld                    ! Cloud overlap method
                                                        !    0: Clear only
                                                        !    1: Random
 real(kind=r8), intent(in) :: play(:,:)            ! Layer pressures (hPa, mb)
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: plev(:,:)            ! Interface pressures (hPa, mb)
                                                        !    Dimensions: (pcols,nlay+1)
      real(kind=r8), intent(in) :: tlay(:,:)            ! Layer temperatures (K)
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: tlev(:,:)            ! Interface temperatures (K)
                                                        !    Dimensions: (pcols,nlay+1)
      real(kind=r8), intent(in) :: tsfc(:)              ! Surface temperature (K)
                                                        !    Dimensions: (pcols)
      real(kind=r8), intent(in) :: h2ovmr(:,:)          ! H2O volume mixing ratio
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: o3vmr(:,:)           ! O3 volume mixing ratio
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: co2vmr(:,:)          ! CO2 volume mixing ratio
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: ch4vmr(:,:)          ! Methane volume mixing ratio
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: cfc11vmr(:,:)        ! CFC11 volume mixing ratio
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: cfc12vmr(:,:)        ! CFC12 volume mixing ratio
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: cfc22vmr(:,:)        ! CFC22 volume mixing ratio
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: ccl4vmr(:,:)         ! CCL4 volume mixing ratio
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: emis(:,:)            ! Surface emissivity
                                                        !    Dimensions: (pcols,nbndlw)

      integer, intent(in) :: inflglw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflglw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflglw                   ! Flag for liquid droplet specification
      integer, parameter :: nsubclw = ngptlw
      real(kind=r8), intent(in) :: cldfmcl(:, :, :)       ! Cloud fraction
                                                        !    Dimensions: (ngptlw,pcols,nlay)
      real(kind=r8), intent(in) :: ciwpmcl(:, :, :)       ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ngptlw,pcols,nlay)
      real(kind=r8), intent(in) :: clwpmcl(:, :, :)       ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ngptlw,pcols,nlay)
      real(kind=r8), intent(in) :: reicmcl(:,:)         ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: relqmcl(:,:)         ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(in) :: taucmcl(:, :, :)       ! Cloud optical depth
                                                        !    Dimensions: (ngptlw,pcols,nlay)

      real(kind=r8), intent(in) :: tauaer(:,:,:)        ! aerosol optical depth
                                                        !   at mid-point of LW spectral bands
                                                        !    Dimensions: (pcols,nlay,nbndlw)


! ----- Output -----

      real(kind=r8), intent(out) :: uflx(pcols,nlay+2)           ! Total sky longwave upward flux (W/m2)
                                                        !    Dimensions: (pcols,nlay+1)
      real(kind=r8), intent(out) :: dflx(pcols,nlay+2)           ! Total sky longwave downward flux (W/m2)
                                                        !    Dimensions: (pcols,nlay+1)
      real(kind=r8), intent(out) :: hr(pcols,nlay+1)             ! Total sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (pcols,nlay)
      real(kind=r8), intent(out) :: uflxc(pcols,nlay+2)          ! Clear sky longwave upward flux (W/m2)
                                                        !    Dimensions: (pcols,nlay+1)
      real(kind=r8), intent(out) :: dflxc(pcols,nlay+2)          ! Clear sky longwave downward flux (W/m2)
                                                        !    Dimensions: (pcols,nlay+1)
      real(kind=r8), intent(out) :: hrc(pcols,nlay+1)            ! Clear sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (pcols,nlay)

! ----- Local -----

! Control
      !integer :: nlayers    !by zy                      ! total number of layers
      integer :: istart                         ! beginning band of calculation
      integer :: iend                           ! ending band of calculation
      integer :: iout                           ! output option flag (inactive)
      integer :: iaer                           ! aerosol option flag
     
      integer :: imca                           ! flag for mcica [0=off, 1=on]
      integer :: ims                            ! value for changing mcica permute seed
      integer :: nlayers
      integer :: i,j,k
!***********************************************
integer :: istat
!-------inatm out
real(kind=r8) :: cldfmc_c(140,2048,52)
real(kind=r8) :: ciwpmc_c(140,2048,52)
real(kind=r8) :: clwpmc_c(140,2048,52)
real(kind=r8) :: relqmc_c(2048,52)
real(kind=r8) :: reicmc_c(2048,52)
real(kind=r8) :: dgesmc_c(2048,52)
real(kind=r8) :: taucmc_c(140,2048,52)
real(kind=r8):: taua_c(2048,52,16)
real(kind=r8) :: pavel_c(2048,52)
real(kind=r8) :: tavel_c(2048,52)
real(kind=r8) :: pz_c(2048,0:52)
real(kind=r8) :: tz_c(2048,0:52)
real(kind=r8) :: coldry_c(2048,52)
real(kind=r8) :: wbrodl_c(2048,52)
real(kind=r8) :: wkl_c(2048,38,52)
real(kind=r8) :: wx_c(2048,4,52)
real(kind=r8) :: semiss_c(2048,16)
real(kind=r8) :: tbound_c(2048)            
real(kind=r8) :: pwvcm_c(2048)  

!-----------------------cldprmc out
 !real(kind=r8) :: taucmc_c(pcols,140,52)
 integer :: ncbands_c(2048)

!---------------------setcoef out
 integer :: laytrop_c(2048)
integer :: jp_c(2048,52)
integer :: jt_c(2048,52)
integer :: jt1_c(2048,52)
real(kind=r8) :: planklay_c(2048,52,16)
real(kind=r8) :: planklev_c(2048,0:52,16)
real(kind=r8) :: plankbnd_c(2048,16)
real(kind=r8) :: colh2o_c(2048,52)
real(kind=r8) :: colco2_c(2048,52)
real(kind=r8) :: colo3_c(2048,52)
real(kind=r8) :: coln2o_c(2048,52)
real(kind=r8) :: colco_c(2048,52)
real(kind=r8) :: colch4_c(2048,52)
real(kind=r8) :: colo2_c(2048,52)
real(kind=r8) :: colbrd_c(2048,52)
integer :: indself_c(2048,52)
integer :: indfor_c(2048,52)
real(kind=r8) :: selffac_c(2048,52)
real(kind=r8) :: selffrac_c(2048,52)
real(kind=r8) :: forfac_c(2048,52)
real(kind=r8) :: forfrac_c(2048,52)
integer :: indminor_c(2048,52)
real(kind=r8) :: minorfrac_c(2048,52)
real(kind=r8) :: scaleminor_c(2048,52)
real(kind=r8) :: scaleminorn2_c(2048,52)
real(kind=r8) :: fac00_c(2048,52)
real(kind=r8) :: fac01_c(2048,52)
real(kind=r8) :: fac10_c(2048,52)
real(kind=r8) :: fac11_c(2048,52)

real(kind=r8) :: rat_h2oco2_c(2048,52)
real(kind=r8) :: rat_h2oco2_1_c(2048,52)
real(kind=r8) :: rat_h2oo3_c(2048,52)
real(kind=r8) :: rat_h2oo3_1_c(2048,52)
real(kind=r8) :: rat_h2on2o_c(2048,52)
real(kind=r8) :: rat_h2on2o_1_c(2048,52)
real(kind=r8) :: rat_h2och4_c(2048,52)
real(kind=r8) :: rat_h2och4_1_c(2048,52)
real(kind=r8) :: rat_n2oco2_c(2048,52)
real(kind=r8) :: rat_n2oco2_1_c(2048,52)
real(kind=r8) :: rat_o3co2_c(2048,52)
real(kind=r8) :: rat_o3co2_1_c(2048,52)

!-------------------------------taumol out
real(kind=r8) :: fracs_c(2048,52,140)                                                   
real(kind=r8) :: taug_c(2048,52,140)         
real(kind=r8) :: taut_c(2048,52,140)   

!--------------------------------

 !real(kind=r8) :: a0_c(16)
 ! real(kind=r8) :: a1_c(16)
 !  real(kind=r8) :: a2_c(16)
!--------------------------------
integer,device :: inflag
integer,device :: iceflag
integer,device :: liqflag
!---------------------------
      type(dim3) :: grid1
      type(dim3) :: tBlock1 

      type(dim3) :: grid2
      type(dim3) :: tBlock2

      type(dim3) :: grid4
      type(dim3) :: tBlock4
!  This is the main longitude/column loop within RRTMG.
      oneminus = 1._r8 - 1.e-6_r8
      pi = 2._r8 * asin(1._r8)
      fluxfac = pi * 2.e4_r8                    ! orig:   fluxfac = pi * 2.d4
      istart = 1
      iend = 16
      iout = 0
      ims = 1
      nlayers = nlay + 1
 

  tBlock1=dim3(128,4,1)
  grid1=dim3(ceiling(real(ncol)/tBlock1%x),ceiling(real(nlayers)/tBlock1%y),1)
  tBlock2=dim3(2,128,2)
  grid2=dim3(ceiling(real(ngptlw)/tBlock2%x),ceiling(real(ncol)/tBlock2%y),ceiling(real(nlayers)/tBlock2%z))
 
  tBlock4=dim3(128,2,2)
  grid4=dim3(ceiling(real(ncol)/tBlock4%x),ceiling(real(nlayers)/tBlock4%y),ceiling(real(ngptlw)/tBlock4%z))

 if (icld.lt.0.or.icld.gt.3) icld = 2


      iaer = 10






    play_d=play
    plev_d=plev
    tlay_d=tlay
    tlev_d=tlev
    tsfc_d=tsfc
    h2ovmr_d=h2ovmr
    o3vmr_d=o3vmr
    co2vmr_d=co2vmr
    ch4vmr_d=ch4vmr
    o2vmr_d=o2vmr
    n2ovmr_d=n2ovmr
    cfc11vmr_d=cfc11vmr
    cfc12vmr_d=cfc12vmr
    cfc22vmr_d=cfc22vmr
    ccl4vmr_d=ccl4vmr
    emis_d=emis
    cldfmcl_d=cldfmcl
    ciwpmcl_d=ciwpmcl
    clwpmcl_d=clwpmcl
    reicmcl_d=reicmcl
    relqmcl_d=relqmcl
    taucmcl_d=taucmcl
    tauaer_d=tauaer

!--------------------------------------
    ixindx_d=ixindx
!-----------------------------------

    absliq1_d=absliq1
    absice0_d=absice0
    absice1_d=absice1
    absice2_d=absice2
    absice3_d=absice3
    ngb_d=ngb
!--------------------------------------------------
    totplnk_d=totplnk
    totplk16_d=totplk16
    preflog_d=preflog
    tref_d=tref
    chi_mls_d=chi_mls

!----------------------------------------------
    fracrefa01_d = fracrefa01
    fracrefb01_d = fracrefb01
    absa01_d = absa01
    absb01_d = absb01
    ka_mn201_d = ka_mn201
    kb_mn201_d = kb_mn201
    selfref01_d = selfref01
    forref01_d = forref01
    fracrefa02_d = fracrefa02
    fracrefb02_d = fracrefb02
    absa02_d = absa02
    absb02_d = absb02
    selfref02_d = selfref02
    forref02_d = forref02
    fracrefa03_d = fracrefa03
    fracrefb03_d = fracrefb03
    absa03_d = absa03
    absb03_d = absb03
    ka_mn2o03_d = ka_mn2o03
    kb_mn2o03_d = kb_mn2o03
    selfref03_d = selfref03
    forref03_d = forref03
    fracrefa04_d = fracrefa04
    fracrefb04_d = fracrefb04
    absa04_d = absa04
    absb04_d = absb04
    selfref04_d = selfref04
    forref04_d = forref04
    fracrefa05_d = fracrefa05
    fracrefb05_d = fracrefb05
    absa05_d = absa05
    absb05_d = absb05
    ka_mo305_d = ka_mo305
    selfref05_d = selfref05
    forref05_d = forref05
    ccl405_d = ccl405
    fracrefa06_d = fracrefa06
    absa06_d = absa06
    ka_mco206_d = ka_mco206
    selfref06_d = selfref06
    forref06_d = forref06
    cfc11adj06_d = cfc11adj06
    cfc1206_d = cfc1206
    fracrefa07_d = fracrefa07
    fracrefb07_d = fracrefb07
    absa07_d = absa07
    absb07_d = absb07
    ka_mco207_d = ka_mco207
    kb_mco207_d = kb_mco207
    selfref07_d = selfref07
    forref07_d = forref07
    fracrefa08_d = fracrefa08
    fracrefb08_d = fracrefb08
    absa08_d = absa08
    absb08_d = absb08
    ka_mco208_d = ka_mco208
    ka_mn2o08_d = ka_mn2o08
    ka_mo308_d = ka_mo308
    kb_mco208_d = kb_mco208
    kb_mn2o08_d = kb_mn2o08
    selfref08_d = selfref08
    forref08_d = forref08
    cfc1208_d = cfc1208
    cfc22adj08_d = cfc22adj08
    fracrefa09_d = fracrefa09
    fracrefb09_d = fracrefb09
    absa09_d = absa09
    absb09_d = absb09
    ka_mn2o09_d = ka_mn2o09
    kb_mn2o09_d = kb_mn2o09
    selfref09_d = selfref09
    forref09_d = forref09
    fracrefa10_d = fracrefa10
    fracrefb10_d = fracrefb10
    absa10_d = absa10
    absb10_d = absb10
    selfref10_d = selfref10
    forref10_d = forref10
    fracrefa11_d = fracrefa11
    fracrefb11_d = fracrefb11
    absa11_d = absa11
    absb11_d = absb11
    ka_mo211_d = ka_mo211
    kb_mo211_d = kb_mo211
    selfref11_d = selfref11
    forref11_d = forref11
    fracrefa12_d = fracrefa12
    absa12_d = absa12
    selfref12_d = selfref12
    forref12_d = forref12
    fracrefa13_d = fracrefa13
    fracrefb13_d = fracrefb13
    absa13_d = absa13
    ka_mco213_d = ka_mco213
    ka_mco13_d = ka_mco13
    kb_mo313_d = kb_mo313
    selfref13_d = selfref13
    forref13_d = forref13
    fracrefa14_d = fracrefa14
    fracrefb14_d = fracrefb14
    absa14_d = absa14
    absb14_d = absb14
    selfref14_d = selfref14
    forref14_d = forref14
    fracrefa15_d = fracrefa15
    absa15_d = absa15
    ka_mn215_d = ka_mn215
    selfref15_d = selfref15
    forref15_d = forref15
    fracrefa16_d = fracrefa16
    fracrefb16_d = fracrefb16
    absa16_d = absa16
    absb16_d = absb16
    selfref16_d = selfref16
    forref16_d = forref16
    nspa_d = nspa
    nspb_d = nspb
    !------------------------------
     delwave_d=delwave
    tau_tbl_d=tau_tbl
    exp_tbl_d=exp_tbl
    tfn_tbl_d=tfn_tbl
    ngs_d=ngs

 call t_startf('kernel')

call inatm_d1<<<grid1,tBlock1>>>(ncol,nlayers,nbndlw)
call inatm_d2<<<grid1,tBlock1>>>(nlayers,ncol,ngptlw,nbndlw)

call inatm_d3<<<grid2,tBlock2>>>(ncol,nlayers,ngptlw,icld)


call inatm_d4<<<ceiling(real(ncol)/512),512>>>(nlayers,ncol,avogad,grav,nmol,nbndlw)
call inatm_d5<<<grid1,tBlock1>>>(ncol,nlayers,nbndlw,icld,iaer,inflglw,iceflglw,liqflglw,inflag,iceflag,liqflag)


! cldfmc_c=cldfmc
! ciwpmc_c=ciwpmc
! clwpmc_c=clwpmc
! relqmc_c=relqmc
! reicmc_c=reicmc
! dgesmc_c=dgesmc
! taucmc_c=taucmc
! taua_c=taua
! pavel_c=pavel
! tavel_c=tavel
! pz_c=pz
! tz_c=tz
! coldry_c=coldry
! wbrodl_c=wbrodl
! wkl_c=wkl
! wx_c=wx
! semiss_c=semiss
! tbound_c=tbound
! pwvcm_c=pwvcm






! do j=1,140
!  do i=1,2048
!    do k=1,52
! write(*,*)'cldfmc_c=',cldfmc_c(40,500,20),cldfmc_c(140,2048,52)
!    enddo
!  enddo
! enddo


!do j=1,140
!   do i=1,2048
!  do k=1,52
!write(*,*)'ciwpmc_c=',ciwpmc_c(40,500,20),ciwpmc_c(140,2048,52)
!  enddo
!enddo
!enddo


!do j=1,140
!   do i=1,2048
!       do k=1,52
!write(*,*)'clwpmc_c=',clwpmc_c(40,500,20),clwpmc_c(140,2048,52)
!  enddo
!enddo
!enddo

!do i=1,2048
!do j=1,52
!write(*,*)'relqmc_c=',relqmc_c(500,20),relqmc_c(2048,52)
!enddo
!enddo




!do i=1,2048
!do j=1,52
!write(*,*)'reicmc_c=',reicmc_c(500,20),reicmc_c(2048,52)
!enddo
!enddo

!do i=1,2048
!do j=1,52
!write(*,*)'dgesmc_c=',dgesmc_c(500,20),dgesmc_c(2048,52)
!enddo
!enddo

!do j=1,140
!  do i=1,2048
!  do k=1,52
!write(*,*)'taucmc_c=',taucmc_c(40,500,20),taucmc_c(140,2048,52)
!  enddo
!enddo
!enddo


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!do i=1,2048
!do j=1,52
!  do k=1,16
!write(*,*)'taua_c=',taua_c(500,20,6),taua_c(2048,52,16)
!   enddo
!enddo
!enddo

! do i=1,2048
!do j=1,52
!write(*,*)'pavel_c=',pavel_c(500,20),pavel_c(2048,52)
!enddo
!enddo


! do i=1,2048
!do j=1,52
!write(*,*)'tavel_c=',tavel_c(500,20),tavel_c(2048,52)
!enddo
!enddo




!do i=1,2048
!do j=1,52
!write(*,*)'pz_c=',pz_c(500,20),pz_c(2048,52)
!enddo
!enddo
!enddo



!do i=1,2048
!do j=1,52
!write(*,*)'tz_c=',tz_c(500,20),tz_c(2048,52)
!enddo
!enddo


! do i=1,2048
!do j=1,52
!write(*,*)'coldry_c=',coldry_c(500,20),coldry_c(2048,52)
!enddo
!enddo

! do i=1,2048
!do j=1,52
!write(*,*)'wbrodl_c=',wbrodl_c(500,20),wbrodl_c(2048,52)
!enddo
!enddo

!do i=1,2048
!do j=1,38
!  do k=1,52
!write(*,*)'wkl_c=',wkl_c(500,20,20),wkl_c(2048,38,52)
!  enddo
!enddo
!enddo

!do i=1,2048
!do j=1,4
!  do k=1,52
!write(*,*)'wx_c=',wx_c(500,2,20),wx_c(2048,4,52)
!  enddo
!  enddo
!enddo
!enddo

 !do i=1,2048
!do j=1,16
!write(*,*)'semiss_c=',semiss_c(500,6),semiss_c(2048,16)
!enddo
!enddo

! do i=1,2048
!write(*,*)'tbound_c=',tbound_c(500),tbound_c(2048)
!enddo

! do i=1,2048
!write(*,*)'pwvcm_c=',pwvcm_c(500),pwvcm_c(2048)
!enddo

call cldprmc_d<<<ceiling(real(ncol)/512),512>>>(ncol,nlay,nlayers,ngptlw,absliq0, inflag, iceflag, liqflag)

!taucmc_c=taucmc
!ncbands_c=ncbands

!do i=1,2048

!do j=1,140
!  do k=1,52
!write(100+iam,*)i,j,k,'taucmc_c2=',taucmc_c(i,j,k)
!  enddo
!enddo


!write(100+iam,*)i,'ncbands_c=',ncbands_c(i)

!enddo
!iplon


call setcoef_d1<<<grid1,tBlock1>>>(ncol,nlay,nlayers, istart)
call setcoef_d2<<<ceiling(real(ncol)/512),512>>>(ncol,nlayers)


!laytrop_c=laytrop
!jp_c=jp
!jt_c=jt
!jt1_c=jt1
!planklay_c=planklay
!planklev_c=planklev
!plankbnd_c=plankbnd
!colh2o_c=colh2o
!colco2_c=colco2
!colo3_c=colo3
!coln2o_c=coln2o
!colco_c=colco
!colch4_c=colch4
!colo2_c=colo2
!colbrd_c=colbrd
!indself_c=indself
!indfor_c=indfor
!selffac_c=selffac
!selffrac_c=selffrac
!forfac_c=forfac
!forfrac_c=forfrac
!indminor_c=indminor
!minorfrac_c=minorfrac
!scaleminor_c=scaleminor
!scaleminorn2_c=scaleminorn2
!fac00_c=fac00
!fac01_c=fac01
!fac10_c=fac10
!fac11_c=fac11

!rat_h2oco2_c=rat_h2oco2
!rat_h2oco2_1_c=rat_h2oco2_1
!rat_h2oo3_c=rat_h2oo3
!rat_h2oo3_1_c=rat_h2oo3_1
!rat_h2on2o_c=rat_h2on2o
!rat_h2on2o_1_c=rat_h2on2o_1
!rat_h2och4_c=rat_h2och4
!rat_h2och4_1_c=rat_h2och4_1
!rat_n2oco2_c=rat_n2oco2
!rat_n2oco2_1_c=rat_n2oco2_1
!rat_o3co2_c=rat_o3co2
!rat_o3co2_1_c=rat_o3co2_1

!do i=1,2048

!write(100+iam,*)i,'laytrop_c',laytrop_c(i)

!do j=1,52
!write(100+iam,*)i,j,'jp_c',jp_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'jt_c',jt_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'jt1_c',jt1_c(i,j)
!enddo

!do j=1,52
!  do k=1,16
!write(100+iam,*)i,j,k,'planklay_c',planklay_c(i,j,k)
!  enddo
!enddo

!do j=1,52
!  do k=1,16
!write(100+iam,*)i,j,k,'planklev_c',planklev_c(i,j,k)
!  enddo
!enddo

!do j=1,16
!write(100+iam,*)i,j,'plankbnd_c',plankbnd_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'colh2o_c',colh2o_c(i,j)
!enddo


!do j=1,52
!write(100+iam,*)i,j,'colco2_c',colco2_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'colo3_c',colo3_c(i,j)
!enddo

!()()()((()(()())((()()()()()(()()()()()()()()())))))
!do j=1,52
!write(100+iam,*)i,j,'coln2o_c',coln2o_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'colco_c',colco_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'colch4_c',colch4_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'colo2_c',colo2_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'colbrd_c',colbrd_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'indself_c',indself_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'indfor_c',indfor_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'selffac_c',selffac_c(i,j)
!enddo


!do j=1,52

!write(100+iam,*)i,j,'selffrac_c',selffrac_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'forfac_c',forfac_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'forfrac_c',forfrac_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'indminor_c',indminor_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'minorfrac_c',minorfrac_c(i,j)
!enddo


!do j=1,52
!write(100+iam,*)i,j,'scaleminor_c',scaleminor_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'scaleminorn2_c',scaleminorn2_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'fac00_c',fac00_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'fac01_c',fac01_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'fac10_c',fac10_c(i,j)
!enddo


!do j=1,52
!write(100+iam,*)i,j,'fac11_c',fac11_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'rat_h2oco2_c',rat_h2oco2_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'rat_h2oco2_1_c',rat_h2oco2_1_c(i,j)
!enddo


!do j=1,52
!write(100+iam,*)i,j,'rat_h2oo3_c',rat_h2oo3_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'rat_h2oo3_1_c',rat_h2oo3_1_c(i,j)
!enddo



!do j=1,52
!write(100+iam,*)i,j,'rat_h2on2o_c',rat_h2on2o_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'rat_h2on2o_1_c',rat_h2on2o_1_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'rat_h2och4_c',rat_h2och4_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'rat_h2och4_1_c',rat_h2och4_1_c(i,j)
!enddo
!do j=1,52
!write(100+iam,*)i,j,'rat_n2oco2_c',rat_n2oco2_c(i,j)
!enddo
!do j=1,52
!write(100+iam,*)i,j,'rat_n2oco2_1_c',rat_n2oco2_1_c(i,j)
!enddo

!do j=1,52
!write(100+iam,*)i,j,'rat_o3co2_c',rat_o3co2_c(i,j)
!enddo


!do j=1,52
!write(100+iam,*)i,j,'rat_o3co2_1_c',rat_o3co2_1_c(i,j)
!enddo


!enddo
!iplon

call taumol_d<<<grid1,tBlock1>>>(ncol,nlayers,oneminus)
!--------------------------------------------------------------------
call taumol_d2<<<grid4,tBlock4>>>(ncol,nlayers,ngptlw)



!taut_c=taut
!pz_c=pz
!pwvcm_c=pwvcm
!semiss_c=semiss
!planklay_c=planklay
!planklev_c=planklev
!plankbnd_c=plankbnd
!fracs_c=fracs
!ncbands_c=ncbands
!cldfmc_c=cldfmc
!taucmc_c=taucmc



!do i=1,2048


  !do j=0,52
!write(100+iam,*)i,j,'pz=',pz_c(i,j)
  !enddo




!write(100+iam,*)i,'pwvcm=',pwvcm_c(i)



 ! do j=1,16
!write(100+iam,*)i,j,'semiss=',semiss_c(i,j)
 ! enddo



  !do j=1,52
   ! do k=1,16
!write(100+iam,*)i,j,k,'planklay=',planklay_c(i,j,k)
 !   enddo
  !enddo



!  do j=0,52
 !   do k=1,16
!write(100+iam,*)i,j,k,'planklev=',planklev_c(i,j,k)
 !   enddo
 ! enddo



 ! do j=1,16
!write(100+iam,*)i,j,'plankbnd=',plankbnd_c(i,j)
 ! enddo



  !do j=1,52
   ! do k=1,140
!write(100+iam,*)i,j,k,'fracs=',fracs_c(i,j,k)
 !   enddo
  !enddo




  !do j=1,52
  !  do k=1,140
!write(200+iam,*)i,j,k,'taut=',taut_c(i,j,k)
  !  enddo
  !enddo

!do j=1,nlayers
 ! do k=1,ngptlw
!write(100+iam,*)i,j,k,'taug=',taug_c(i,j,k)
 ! enddo
!enddo

!do j=1,nlayers
  !do k=1,nbndlw
!write(100+iam,*)i,j,k,'taua=',taua_c(i,j,k)
  !enddo
!enddo

!write(100+iam,*)i,'ncbands=',ncbands_c(i)


!write(100+iam,*)i,'laytrop=',laytrop_c(i)
  





 ! do j=1,140
  !  do k=1,52
!write(100+iam,*)i,j,k,'cldfmc=',cldfmc_c(i,j,k)
 !   enddo
 ! enddo



  !do j=1,140
   ! do k=1,52
!write(100+iam,*)i,j,k,'taucmc=',taucmc_c(i,j,k)
 !   enddo
  !enddo


!enddo
!iplon


!----------------------------------------------------------------------

!fracs_c=fracs
!taug_c=taug
!taut_c=taut

!write(*,*)'fracs_c',fracs_c(1,1,1),fracs_c(2048,52,140)
!write(*,*)'taug_c',taug_c(1,1,1),taug_c(2048,52,140)
!write(*,*)'taut_c',taut_c(1,1,1),taut_c(2048,52,140)

!a0=a0_c
!a1=a1_c
!a2=a2_c
call rtrnmc_d<<<ceiling(real(pcols)/64),64>>>(ncol,nlayers,istart,iend,iout,nbndlw,ngptlw,tblint,bpade,fluxfac,heatfac)


istat=cudaDeviceSynchronize()
 call t_stopf('kernel')

uflx=uflx_d
dflx=dflx_d
uflxc=uflxc_d
dflxc=dflxc_d
hr=hr_d
hrc=hrc_d




! do i=1,2048


! do j=1,53
! write(1000+iam,*)i,j,'uflx=',uflx(i,j)
! enddo

! do j=1,53
! write(1000+iam,*)i,j,'dflx=',dflx(i,j)
! enddo

! do j=1,53
! write(1000+iam,*)i,j,'uflxc=',uflxc(i,j)
! enddo

! do j=1,53
! write(1000+iam,*)i,j,'dflxc=',dflxc(i,j)
! enddo

! do j=1,52
! write(1000+iam,*)i,j,'hrc=',hrc(i,j)
! enddo

! do j=1,52
! write(1000+iam,*)i,j,'hr=',hr(i,j)
! enddo


! enddo




  


end subroutine rrtmg_lw

end module
