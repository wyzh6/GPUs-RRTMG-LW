
module kernel
	use shr_kind_mod, only: r8 => shr_kind_r8
	use cudafor
	!------inatm in
real(kind=r8),device :: play_d(2048,51)
real(kind=r8),device :: plev_d(2048,52)
real(kind=r8),device :: tlay_d(2048,51)
real(kind=r8),device :: tlev_d(2048,52)
real(kind=r8),device :: tsfc_d(2048)
real(kind=r8),device :: h2ovmr_d(2048,51)
real(kind=r8),device :: o3vmr_d(2048,51)
real(kind=r8),device :: co2vmr_d(2048,51)
real(kind=r8),device :: ch4vmr_d(2048,51)
real(kind=r8),device :: o2vmr_d(2048,51)
real(kind=r8),device :: n2ovmr_d(2048,51)
real(kind=r8),device :: cfc11vmr_d(2048,51)
real(kind=r8),device :: cfc12vmr_d(2048,51)
real(kind=r8),device :: cfc22vmr_d(2048,51)
real(kind=r8),device :: ccl4vmr_d(2048,51)
real(kind=r8),device :: emis_d(2048,16)
real(kind=r8),device :: cldfmcl_d(140, 2048, 51)
real(kind=r8),device :: ciwpmcl_d(140, 2048, 51)
real(kind=r8),device :: clwpmcl_d(140, 2048, 51)
real(kind=r8),device :: reicmcl_d(2048,51)
real(kind=r8),device :: relqmcl_d(2048,51)
real(kind=r8),device :: taucmcl_d(140, 2048, 51)
real(kind=r8),device :: tauaer_d(2048,51,16)


!-------inatm use
integer,device :: ixindx_d(38)
!integer,device :: ngb_d(140)

!-------inatm out
real(kind=r8),device :: cldfmc(140,2048,52)
real(kind=r8),device :: ciwpmc(140,2048,52)
real(kind=r8),device :: clwpmc(140,2048,52)
real(kind=r8),device :: relqmc(2048,52)
real(kind=r8),device :: reicmc(2048,52)
real(kind=r8),device :: dgesmc(2048,52)
real(kind=r8),device :: taucmc(140,2048,52)
real(kind=r8),device :: taua(2048,52,16)
real(kind=r8),device :: pavel(2048,52)
real(kind=r8),device :: tavel(2048,52)
real(kind=r8),device :: pz(2048,0:52)
real(kind=r8),device :: tz(2048,0:52)
real(kind=r8),device :: coldry(2048,52)
real(kind=r8),device :: wbrodl(2048,52)
real(kind=r8),device :: wkl(2048,38,52)
real(kind=r8),device :: wx(2048,4,52)
real(kind=r8),device :: semiss(2048,16)
real(kind=r8),device :: tbound(2048)            
real(kind=r8),device :: pwvcm(2048)   

!*******************************************************8
!----------------------cldprmc use
real(kind=r8),device :: absliq1_d(58,16)
real(kind=r8),device :: absice0_d(2)
real(kind=r8),device :: absice1_d(2,5)
real(kind=r8),device :: absice2_d(43,16)
real(kind=r8),device :: absice3_d(46,16)
integer,device :: ngb_d(140)

!-----------------------cldprmc in
 !real(kind=r8),device :: cldfmc(pcols,140,52)
!real(kind=r8),device :: ciwpmc(pcols,140,52)
!real(kind=r8),device :: clwpmc(pcols,140,52)
!real(kind=r8),device :: relqmc(pcols,52)
!real(kind=r8),device :: reicmc(pcols,52)
!real(kind=r8),device :: dgesmc(pcols,52)
!real(kind=r8),device :: taucmc(pcols,140,52)

!-----------------------cldprmc out
 !real(kind=r8),device :: taucmc(pcols,140,52)
 integer, device :: ncbands(2048)
!**************************************************************
!--------------------setcoef use
 real(kind=r8),device :: totplnk_d(181,16)
 real(kind=r8),device :: totplk16_d(181)
 real(kind=r8),device :: preflog_d(59)
 real(kind=r8),device :: tref_d(59)
 real(kind=r8),device :: chi_mls_d(7,59)
!--------------------setcoef in
 !real(kind=r8),device :: pavel(2048,52)
 !real(kind=r8),device :: tavel(2048,52)
 !real(kind=r8),device :: tavel(2048,52)          
   !real(kind=r8),device :: tz(2048,0:52)                                                
 !real(kind=r8),device :: tbound(2048)      
 ! real(kind=r8),device :: coldry(2048,52)                                                   
 !real(kind=r8),device :: wbrodl(2048,52)    
  !real(kind=r8),device :: wkl(2048,38,52)  
 !real(kind=r8),device :: semiss(2048,16)                                                      
    
 !---------------------setcoef out
 integer,device :: laytrop(2048)
integer,device :: jp(2048,52)
integer,device :: jt(2048,52)
integer,device :: jt1(2048,52)
real(kind=r8),device :: planklay(2048,52,16)
real(kind=r8),device :: planklev(2048,0:52,16)
real(kind=r8),device :: plankbnd(2048,16)
real(kind=r8),device :: colh2o(2048,52)
real(kind=r8),device :: colco2(2048,52)
real(kind=r8),device :: colo3(2048,52)
real(kind=r8),device :: coln2o(2048,52)
real(kind=r8),device :: colco(2048,52)
real(kind=r8),device :: colch4(2048,52)
real(kind=r8),device :: colo2(2048,52)
real(kind=r8),device :: colbrd(2048,52)
integer,device :: indself(2048,52)
integer,device :: indfor(2048,52)
real(kind=r8),device :: selffac(2048,52)
real(kind=r8),device :: selffrac(2048,52)
real(kind=r8),device :: forfac(2048,52)
real(kind=r8),device :: forfrac(2048,52)
integer,device :: indminor(2048,52)
real(kind=r8),device :: minorfrac(2048,52)
real(kind=r8),device :: scaleminor(2048,52)
real(kind=r8),device :: scaleminorn2(2048,52)
real(kind=r8),device :: fac00(2048,52)
real(kind=r8),device :: fac01(2048,52)
real(kind=r8),device :: fac10(2048,52)
real(kind=r8),device :: fac11(2048,52)

real(kind=r8),device :: rat_h2oco2(2048,52)
real(kind=r8),device :: rat_h2oco2_1(2048,52)
real(kind=r8),device :: rat_h2oo3(2048,52)
real(kind=r8),device :: rat_h2oo3_1(2048,52)
real(kind=r8),device :: rat_h2on2o(2048,52)
real(kind=r8),device :: rat_h2on2o_1(2048,52)
real(kind=r8),device :: rat_h2och4(2048,52)
real(kind=r8),device :: rat_h2och4_1(2048,52)
real(kind=r8),device :: rat_n2oco2(2048,52)
real(kind=r8),device :: rat_n2oco2_1(2048,52)
real(kind=r8),device :: rat_o3co2(2048,52)
real(kind=r8),device :: rat_o3co2_1(2048,52)

!************************************************************************
!-------------------------taumol_d use
real(kind=r8),device :: fracrefa01_d(10)
real(kind=r8),device :: fracrefb01_d(10)
real(kind=r8),device :: absa01_d(65,10)
real(kind=r8),device :: absb01_d(235,10)
real(kind=r8),device :: ka_mn201_d(19,10)
real(kind=r8),device :: kb_mn201_d(19,10)
real(kind=r8),device :: selfref01_d(10,10)
real(kind=r8),device :: forref01_d(4,10)

real(kind=r8),device :: fracrefa02_d(12)
real(kind=r8),device :: fracrefb02_d(12)
real(kind=r8),device :: absa02_d(65,12)
real(kind=r8),device :: absb02_d(235,12)
real(kind=r8),device :: selfref02_d(10,12)
real(kind=r8),device :: forref02_d(4,12)


real(kind=r8),device :: fracrefa03_d(16,10)
real(kind=r8),device :: fracrefb03_d(16,5)
real(kind=r8),device :: absa03_d(585,16)
real(kind=r8),device :: absb03_d(1175,16)
real(kind=r8),device :: ka_mn2o03_d(9,19,16)
real(kind=r8),device :: kb_mn2o03_d(5,19,16)
real(kind=r8),device :: selfref03_d(10,16)
real(kind=r8),device :: forref03_d(4,16)

real(kind=r8),device :: fracrefa04_d(14,9)
real(kind=r8),device :: fracrefb04_d(14,6)
real(kind=r8),device :: absa04_d(585,14)
real(kind=r8),device :: absb04_d(1175,14)
real(kind=r8),device :: selfref04_d(10,14)
real(kind=r8),device :: forref04_d(4,14)

real(kind=r8),device :: fracrefa05_d(16,9)
real(kind=r8),device :: fracrefb05_d(16,5)
real(kind=r8),device :: absa05_d(585,16)
real(kind=r8),device :: absb05_d(1175,16)
real(kind=r8),device :: ka_mo305_d(9,19,16)
real(kind=r8),device :: selfref05_d(10,16)
real(kind=r8),device :: forref05_d(4,16)
real(kind=r8),device :: ccl405_d(16)

real(kind=r8),device :: fracrefa06_d(8)
real(kind=r8),device :: absa06_d(65,8)
real(kind=r8),device :: ka_mco206_d(19,8)
real(kind=r8),device :: selfref06_d(10,8)
real(kind=r8),device :: forref06_d(4,8)
real(kind=r8),device :: cfc11adj06_d(8)
real(kind=r8),device :: cfc1206_d(8)

real(kind=r8),device :: fracrefa07_d(12,9)
real(kind=r8),device :: fracrefb07_d(12)
real(kind=r8),device :: absa07_d(585,12)
real(kind=r8),device :: absb07_d(235,12)
real(kind=r8),device :: ka_mco207_d(9,19,12)
real(kind=r8),device :: kb_mco207_d(19,12)
real(kind=r8),device :: selfref07_d(10,12)
real(kind=r8),device :: forref07_d(4,12)

real(kind=r8),device :: fracrefa08_d(8) 
real(kind=r8),device :: fracrefb08_d(8)
real(kind=r8),device :: absa08_d(65,8)
real(kind=r8),device :: absb08_d(235,8)
real(kind=r8),device :: ka_mco208_d(19,8)
real(kind=r8),device :: ka_mn2o08_d(19,8)
real(kind=r8),device :: ka_mo308_d(19,8)
real(kind=r8),device :: kb_mco208_d(19,8)
real(kind=r8),device :: kb_mn2o08_d(19,8)
real(kind=r8),device :: selfref08_d(10,8)
real(kind=r8),device :: forref08_d(4,8)
real(kind=r8),device :: cfc1208_d(8)
real(kind=r8),device :: cfc22adj08_d(8)

real(kind=r8),device :: fracrefa09_d(12,9)
real(kind=r8),device :: fracrefb09_d(12)
real(kind=r8),device :: absa09_d(585,12)
real(kind=r8),device :: absb09_d(235,12)
real(kind=r8),device :: ka_mn2o09_d(9,19,12)
real(kind=r8),device :: kb_mn2o09_d(19,12)
real(kind=r8),device :: selfref09_d(10,12)
real(kind=r8),device :: forref09_d(4,12)


real(kind=r8),device :: fracrefa10_d(6)
real(kind=r8),device :: fracrefb10_d(6)
real(kind=r8),device :: absa10_d(65,6)
real(kind=r8),device :: absb10_d(235,6)
real(kind=r8),device :: selfref10_d(10,6)
real(kind=r8),device :: forref10_d(4,6)

real(kind=r8),device :: fracrefa11_d(8)
real(kind=r8),device :: fracrefb11_d(8)
real(kind=r8),device :: absa11_d(65,8)
real(kind=r8),device :: absb11_d(235,8)
real(kind=r8),device :: ka_mo211_d(19,8)
real(kind=r8),device :: kb_mo211_d(19,8)
real(kind=r8),device :: selfref11_d(10,8)
real(kind=r8),device :: forref11_d(4,8)

real(kind=r8),device :: fracrefa12_d(8,9)
real(kind=r8),device :: absa12_d(585,8)
real(kind=r8),device :: selfref12_d(10,8)
real(kind=r8),device :: forref12_d(4,8)


real(kind=r8),device :: fracrefa13_d(4,9)
real(kind=r8),device :: fracrefb13_d(4)
real(kind=r8),device :: absa13_d(585,4)
real(kind=r8),device :: ka_mco213_d(9,19,4)
real(kind=r8),device :: ka_mco13_d(9,19,4)
real(kind=r8),device :: kb_mo313_d(19,4)
real(kind=r8),device :: selfref13_d(10,4)
real(kind=r8),device :: forref13_d(4,4)


real(kind=r8),device :: fracrefa14_d(2)
real(kind=r8),device :: fracrefb14_d(2)
real(kind=r8),device :: absa14_d(65,2)
real(kind=r8),device :: absb14_d(235,2)
real(kind=r8),device :: selfref14_d(10,2)
real(kind=r8),device :: forref14_d(4,2)


real(kind=r8),device :: fracrefa15_d(2,9)
real(kind=r8),device :: absa15_d(585,2)
real(kind=r8),device :: ka_mn215_d(9,19,2)
real(kind=r8),device :: selfref15_d(10,2)
real(kind=r8),device :: forref15_d(4,2)


real(kind=r8),device :: fracrefa16_d(2,9)
real(kind=r8),device :: fracrefb16_d(2) 
real(kind=r8),device :: absa16_d(585,2)
real(kind=r8),device :: absb16_d(235,2)
real(kind=r8),device :: selfref16_d(10,2)
real(kind=r8),device :: forref16_d(4,2)

integer,device :: nspa_d(16)
integer,device :: nspb_d(16)

!-------------------------------taumol in

!-------------------------------taumol out
real(kind=r8), device :: fracs(2048,52,140)                                                   
real(kind=r8), device :: taug(2048,52,140)         
real(kind=r8),device :: taut(2048,52,140)                                                    
!-----------------------------rtrnmc use
real(kind=r8),device :: delwave_d(16)
integer,device :: ngs_d(16)
real(kind=r8),device :: tau_tbl_d(0:10000)
real(kind=r8),device :: exp_tbl_d(0:10000)
real(kind=r8),device :: tfn_tbl_d(0:10000)
!!!!integer,device :: ngb_d(140)
!-----------------------------rtrnmc in
 !real(kind=r8),device :: pz(2048,0:52)
!real(kind=r8),device :: pwvcm(2048) 
!real(kind=r8),device :: semiss(2048,16)
!real(kind=r8),device :: planklay(2048,52,16)
!real(kind=r8),device :: planklev(2048,0:52,16)
!real(kind=r8),device :: plankbnd(2048,16)
!real(kind=r8), device :: fracs(2048,52,140) 
!real(kind=r8),device :: taut(2048,52,140)
!integer, device :: ncbands(2048)
!real(kind=r8),device :: cldfmc(2048,140,52)
!real(kind=r8),device :: taucmc(2048,140,52)


!------------------------------rtrnmc out 
real(kind=r8),device :: fnet(2048,0:52)
real(kind=r8),device :: fnetc(2048,0:52)
real(kind=r8),device :: totuflux(2048,0:52)
real(kind=r8),device :: totdflux(2048,0:52)
real(kind=r8),device :: htr(2048,0:52)
real(kind=r8),device :: totuclfl(2048,0:52)
real(kind=r8),device :: totdclfl(2048,0:52)
real(kind=r8),device :: htrc(2048,0:52)

!-----------------------------rtrnmc local
 real(kind=r8),device :: abscld(2048,52,140)
 real(kind=r8),device :: odcld(2048,52,140)
 real(kind=r8),device :: efclfrac(2048,52,140)


 !real(kind=r8),device :: a0(16)
  !real(kind=r8),device :: a1(16)
  ! real(kind=r8),device :: a2(16)
!------------------------------out!!!!
real(kind=r8),device :: uflx_d(2048,53)
real(kind=r8),device :: dflx_d(2048,53)
real(kind=r8),device :: uflxc_d(2048,53)
real(kind=r8),device :: dflxc_d(2048,53)
real(kind=r8),device :: hr_d(2048,52)
real(kind=r8),device :: hrc_d(2048,52)

!----------------------------------
real(kind=r8),device :: aaa1
real(kind=r8),device :: aaa2
real(kind=r8),device :: aaa3
real(kind=r8),device :: aaa4
real(kind=r8),device :: aaa5
contains
attributes(global)subroutine inatm_d1(ncol,nlayers,nbndlw)
implicit none
 
     integer, value :: nlayers 
     integer, value :: ncol
     integer, value :: nbndlw


 integer :: iplon,lay
 integer :: i,j
 iplon=(blockIdx%x-1)*blockDim%x+threadIdx%x
     lay=(blockIdx%y-1)*blockDim%y+threadIdx%y
if((iplon>=1 .and. iplon <=ncol).and. (lay>=1 .and. lay<=nlayers))then 


do i=1,38
wkl(iplon,i,lay)=0.0_r8
end do
 
       

do i=1,4
wx(iplon,i,lay) = 0.0_r8
end do

reicmc(iplon,lay) = 0.0_r8



dgesmc(iplon,lay) = 0.0_r8
     
      

relqmc(iplon,lay) = 0.0_r8
     
      
do j=1,nbndlw

taua(iplon,lay,j) = 0.0_r8

end do     


endif

if((iplon>=1 .and. iplon <=ncol).and. (lay>=1 .and. lay<=nlayers+1))then 


 pz(iplon,lay-1)= 0.0_r8


 tz(iplon,lay-1)= 0.0_r8

endif
aaa1=10
end subroutine inatm_d1

  attributes(global)subroutine inatm_d2(nlayers,ncol,ngptlw,nbndlw)
implicit none
! ------- Declarations -------
     
     integer, value :: nlayers 
     integer, value :: ncol
     integer,value :: ngptlw
     integer,value :: nbndlw

      integer :: iplon,lay
  
     

     
! Add one to nlayers here to include extra model layer at top of atmosphere
     

!  Initialize all molecular amounts and cloud properties to zero here, then pass input amounts
!  into RRTM arrays below.
 iplon=(blockIdx%x-1)*blockDim%x+threadIdx%x
     lay=(blockIdx%y-1)*blockDim%y+threadIdx%y    


if  ((iplon>=1 .and. iplon<=ncol) .and. (lay>=1 .and. lay<=nlayers-1))then
   
         pavel(iplon,lay) = play_d(iplon,nlayers-lay)
         tavel(iplon,lay) = tlay_d(iplon,nlayers-lay)
         pz(iplon,lay) = plev_d(iplon,nlayers-lay)
         tz(iplon,lay) = tlev_d(iplon,nlayers-lay)
         wkl(iplon,1,lay) = h2ovmr_d(iplon,nlayers-lay)
         wkl(iplon,2,lay) = co2vmr_d(iplon,nlayers-lay)
         wkl(iplon,3,lay) = o3vmr_d(iplon,nlayers-lay)
         wkl(iplon,4,lay) = n2ovmr_d(iplon,nlayers-lay)
         wkl(iplon,6,lay) = ch4vmr_d(iplon,nlayers-lay)
         wkl(iplon,7,lay) = o2vmr_d(iplon,nlayers-lay)
         wx(iplon,1,lay) = ccl4vmr_d(iplon,nlayers-lay)
         wx(iplon,2,lay) = cfc11vmr_d(iplon,nlayers-lay)
         wx(iplon,3,lay) = cfc12vmr_d(iplon,nlayers-lay)
         wx(iplon,4,lay) = cfc22vmr_d(iplon,nlayers-lay)
         
endif
aaa2=10
end subroutine inatm_d2


attributes(global)subroutine inatm_d3(ncol,nlayers,ngptlw,icld)
implicit none

integer,value :: ncol
integer,value :: nlayers
integer,value :: ngptlw
integer,value :: icld

integer :: iplon,lay,ig

  ig=(blockIdx%x-1)*blockDim%x+threadIdx%x
     iplon=(blockIdx%y-1)*blockDim%y+threadIdx%y 
    lay=(blockIdx%z-1)*blockDim%z+threadIdx%z 

if  ((iplon>=1 .and. iplon<=ncol).and.(ig>=1 .and. ig<=ngptlw).and.(lay>=1 .and. lay<=nlayers-1))then

 if (icld .ge. 1) then
cldfmc(ig,iplon,lay) = cldfmcl_d(ig,iplon,nlayers-lay)
taucmc(ig,iplon,lay) = taucmcl_d(ig,iplon,nlayers-lay)
ciwpmc(ig,iplon,lay) = ciwpmcl_d(ig,iplon,nlayers-lay)
clwpmc(ig,iplon,lay) = clwpmcl_d(ig,iplon,nlayers-lay) 
endif
end if
aaa3=10

end subroutine inatm_d3



attributes(global)subroutine inatm_d4(nlayers,ncol,avogad,grav,nmol,nbndlw)
implicit none
! ------- Declarations -------
     
     integer, value :: nlayers 
     integer, value :: ncol
                        
     real(kind=r8),value :: avogad 
     real(kind=r8),value :: grav
     integer,value :: nmol 
     
     integer,value :: nbndlw



      integer :: iplon,lay
      integer :: ix, n, imol,i,j           ! Loop indices
      real(kind=r8) :: amm, amttl, wvttl, wvsh, summol


      
! ----- Local -----
      real(kind=r8), parameter :: amd = 28.9660_r8      ! Effective molecular weight of dry air (g/mol)
      real(kind=r8), parameter :: amw = 18.0160_r8      ! Molecular weight of water vapor (g/mol)
      real(kind=r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
      real(kind=r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
      real(kind=r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
      real(kind=r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
      real(kind=r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
      real(kind=r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
      real(kind=r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12
      real(kind=r8), parameter :: sbc = 5.67e-08_r8     ! Stefan-Boltzmann constant (W/m2K4)
      
amttl = 0.0_r8
wvttl = 0.0_r8

iplon=(blockIdx%x-1)*blockDim%x+threadIdx%x



if  (iplon>=1 .and. iplon<=ncol)  then


!  Set surface temperature.
tbound(iplon) = tsfc_d(iplon)
pz(iplon,0) = plev_d(iplon,nlayers)
tz(iplon,0) = tlev_d(iplon,nlayers)

           



  do lay=1,nlayers-1
          
            amm = (1._r8 - wkl(iplon,1,lay)) * amd + wkl(iplon,1,lay) * amw
         coldry(iplon,lay) = (pz(iplon,lay-1)-pz(iplon,lay)) * 1.e3_r8 * avogad / &
                     (1.e2_r8 * grav * amm * (1._r8 + wkl(iplon,1,lay)))
 
enddo

      pavel(iplon,nlayers) = 0.5_r8 * pz(iplon,nlayers-1)
      tavel(iplon,nlayers) = tavel(iplon,nlayers-1)
      pz(iplon,nlayers) = 1.e-4_r8
      tz(iplon,nlayers-1) = 0.5_r8 * (tavel(iplon,nlayers)+tavel(iplon,nlayers-1))
      tz(iplon,nlayers) = tz(iplon,nlayers-1)
      wkl(iplon,1,nlayers) = wkl(iplon,1,nlayers-1)
      wkl(iplon,2,nlayers) = wkl(iplon,2,nlayers-1)
      wkl(iplon,3,nlayers) = wkl(iplon,3,nlayers-1)
      wkl(iplon,4,nlayers) = wkl(iplon,4,nlayers-1)
      wkl(iplon,6,nlayers) = wkl(iplon,6,nlayers-1)
      wkl(iplon,7,nlayers) = wkl(iplon,7,nlayers-1)
      amm = (1._r8 - wkl(iplon,1,nlayers-1)) * amd + wkl(iplon,1,nlayers-1) * amw
      coldry(iplon,nlayers) = (pz(iplon,nlayers-1)) * 1.e3_r8 * avogad / &
                        (1.e2_r8 * grav * amm * (1._r8 + wkl(iplon,1,nlayers-1)))
      wx(iplon,1,nlayers) = wx(iplon,1,nlayers-1)
      wx(iplon,2,nlayers) = wx(iplon,2,nlayers-1)
      wx(iplon,3,nlayers) = wx(iplon,3,nlayers-1)
      wx(iplon,4,nlayers) = wx(iplon,4,nlayers-1)


do lay=1,nlayers
       summol = 0.0_r8
         do imol = 2,nmol
            summol = summol + wkl(iplon,imol,lay)
         enddo
         wbrodl(iplon,lay) = coldry(iplon,lay) * (1._r8 - summol)
         do imol = 1, nmol
            wkl(iplon,imol,lay) = coldry(iplon,lay) * wkl(iplon,imol,lay)
         enddo
         amttl = amttl + coldry(iplon,lay)+wkl(iplon,1,lay)
         wvttl = wvttl + wkl(iplon,1,lay)
         do ix = 1,4
            if (ixindx_d(ix) .ne. 0) then
               wx(iplon,ixindx_d(ix),lay) = coldry(iplon,lay) * wx(iplon,ix,lay) * 1.e-20_r8
            endif
         enddo
enddo

      wvsh = (amw * wvttl) / (amd * amttl)
      pwvcm(iplon) = wvsh * (1.e3_r8 * pz(iplon,0)) / (1.e2_r8 * grav)

! Set spectral surface emissivity for each longwave band.

      do n=1,nbndlw
         semiss(iplon,n) = emis_d(iplon,n)

      enddo


endif

aaa4=10
end subroutine inatm_d4

attributes(global)subroutine inatm_d5(ncol,nlayers,nbndlw,icld,iaer,inflglw,iceflglw,liqflglw,&
                                       inflag,iceflag,liqflag)
implicit none

integer, value :: icld  
integer, value :: iaer
integer, value :: inflglw                    
integer, value :: iceflglw                  
integer, value :: liqflglw   

integer,value :: ncol
integer,value :: nlayers
integer,value :: nbndlw

integer, intent(out) :: inflag                    ! flag for cloud property method
integer, intent(out) :: iceflag                   ! flag for ice cloud properties
integer, intent(out) :: liqflag 



integer :: iplon,lay
integer :: ib

iplon=(blockIdx%x-1)*blockDim%x+threadIdx%x
     lay=(blockIdx%y-1)*blockDim%y+threadIdx%y  
! Transfer aerosol optical properties to RRTM variable;
! modify to reverse layer indexing here if necessary.

     

  if  ((iplon>=1 .and. iplon<=ncol).and.(lay>=1 .and. lay<=nlayers-1))then
     if (iaer .ge. 1) then
            do ib = 1, nbndlw
               taua(iplon,lay,ib) = tauaer_d(iplon,nlayers-lay,ib)
               
            enddo
     endif
      
 !Transfer cloud fraction and cloud optical properties to RRTM variables,
! modify to reverse layer indexing here if necessary.

      if (icld .ge. 1) then
         inflag = inflglw
         iceflag = iceflglw
         liqflag = liqflglw

   
            reicmc(iplon,lay) = reicmcl_d(iplon,nlayers-lay)
            if (iceflag .eq. 3) then
               dgesmc(iplon,lay) = 1.5396_r8 * reicmcl_d(iplon,nlayers-lay)
            endif
            relqmc(iplon,lay) = relqmcl_d(iplon,nlayers-lay)
      endif



endif 
end subroutine inatm_d5



attributes(global) subroutine cldprmc_d(ncol,nlay,nlayers,ngptlw,absliq0, inflag, iceflag, liqflag)
                                      
! ------------------------------------------------------------------------------

! Purpose:  Compute the cloud optical depth(s) for each cloudy layer.
implicit none
! ------- Input -------

                         ! total number of layers
      integer,intent(in) :: inflag                     ! see definitions
      integer,intent(in) :: iceflag                    ! see definitions
      integer,intent(in) :: liqflag
  
      integer,value :: ncol
      integer,value :: nlay
      integer,value :: nlayers 
      integer,value :: ngptlw                    ! see definitions
      real(kind=r8),value :: absliq0
      

   

      
     
! ------- Local -------
      integer :: iplon
      integer :: lay                            ! Layer index
      integer :: ib                             ! spectral band index
      integer :: ig                             ! g-point interval index
      integer :: index

      real(kind=r8) :: abscoice(140)         
      real(kind=r8) :: abscoliq(140) 
      real(kind=r8) :: cwp                      ! cloud water path
      real(kind=r8) :: radice                   ! cloud ice effective radius (microns)
      real(kind=r8) :: dgeice                   ! cloud ice generalized effective size
      real(kind=r8) :: factor                   !
      real(kind=r8) :: fint                     !
      real(kind=r8) :: radliq                   ! cloud liquid droplet radius (microns)
      
      real(kind=r8), parameter :: cldmin = 1.e-80_r8 ! minimum value for cloud quantities
     

iplon=(blockIdx%x-1)*blockDim%x+threadIdx%x
    
if(iplon>=1 .and. iplon <=ncol) then


      ncbands(iplon) = 1

     ! Main layer loop
      do lay = 1, nlayers

        do ig = 1, ngptlw
          cwp = ciwpmc(ig,iplon,lay) + clwpmc(ig,iplon,lay)
          if (cldfmc(ig,iplon,lay) .ge. cldmin .and. &
             (cwp .ge. cldmin .or. taucmc(ig,iplon,lay) .ge. cldmin)) then

! Ice clouds and water clouds combined.
            if (inflag .eq. 0) then
! Cloud optical depth already defined in taucmc, return to main program
               return

            !elseif(inflag .eq. 1) then
                !stop 'INFLAG = 1 OPTION NOT AVAILABLE WITH MCICA'
!
! Separate treatement of ice clouds and water clouds.
            elseif(inflag .eq. 2) then
               radice = reicmc(iplon,lay)

! Calculation of absorption coefficients due to ice clouds.
               if (ciwpmc(ig,iplon,lay) .eq. 0.0_r8) then
                  abscoice(ig) = 0.0_r8

               elseif (iceflag .eq. 0) then
                  !if (radice .lt. 10.0_r8) stop 'ICE RADIUS TOO SMALL'


                  abscoice(ig) = absice0_d(1) + absice0_d(2)/radice

               elseif (iceflag .eq. 1) then
! mji - turn off limits to mimic CAM3
!
                 ncbands(iplon) = 5
                  ib = ngb_d(ig)
                  abscoice(ig) = absice1_d(1,ib) + absice1_d(2,ib)/radice

! For iceflag=2 option, combine with iceflag=0 option to handle out of bounds
! particle sizes.
! Use iceflag=2 option for ice particle effective radii from 5.0 and 131.0 microns
! and use iceflag=0 option for ice particles greater than 131.0 microns.
! *** NOTE: Transition between two methods has not been smoothed.

               elseif (iceflag .eq. 2) then
                  !if (radice .lt. 5.0_r8) stop 'ICE RADIUS OUT OF BOUNDS'
                  if (radice .ge. 5.0_r8 .and. radice .le. 131._r8) then
                     ncbands(iplon) = 16
                     factor = (radice - 2._r8)/3._r8
                     index = int(factor)
                     if (index .eq. 43) index = 42
                     fint = factor - float(index)
                     ib = ngb_d(ig)
                     abscoice(ig) = &
                         absice2_d(index,ib) + fint * &
                         (absice2_d(index+1,ib) - (absice2_d(index,ib)))
                  elseif (radice .gt. 131._r8) then
                     abscoice(ig) = absice0_d(1) + absice0_d(2)/radice
                  endif

! For iceflag=3 option, combine with iceflag=0 option to handle large particle sizes.
! Use iceflag=3 option for ice particle effective radii from 3.2 and 91.0 microns
! (generalized effective size, dge, from 5 to 140 microns), and use iceflag=0 option
! for ice particle effective radii greater than 91.0 microns (dge = 140 microns).
! *** NOTE: Fu parameterization requires particle size in generalized effective size.
! *** NOTE: Transition between two methods has not been smoothed.

               elseif (iceflag .eq. 3) then
                  dgeice = dgesmc(iplon,lay)
                 ! if (dgeice .lt. 5.0_r8) stop 'ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'
                  if (dgeice .ge. 5.0_r8 .and. dgeice .le. 140._r8) then
                     ncbands(iplon) = 16
                     factor = (dgeice - 2._r8)/3._r8
                     index = int(factor)
                     if (index .eq. 46) index = 45
                     fint = factor - float(index)
                     ib = ngb_d(ig)
                     abscoice(ig) = &
                         absice3_d(index,ib) + fint * &
                         (absice3_d(index+1,ib) - (absice3_d(index,ib)))
                  elseif (dgeice .gt. 140._r8) then
                     abscoice(ig) = absice0_d(1) + absice0_d(2)/radice
                  endif

               endif

! Calculation of absorption coefficients due to water clouds.
               if (clwpmc(ig,iplon,lay) .eq. 0.0_r8) then
                  abscoliq(ig) = 0.0_r8

               elseif (liqflag .eq. 0) then
                   abscoliq(ig) = absliq0

               elseif (liqflag .eq. 1) then
                  radliq = relqmc(iplon,lay)
                  !if (radliq .lt. 1.5_r8 .or. radliq .gt. 60._r8) stop &
                       !'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
                  index = radliq - 1.5_r8
                  if (index .eq. 58) index = 57
                  if (index .eq. 0) index = 1
                  fint = radliq - 1.5_r8 - index
                  ib = ngb_d(ig)
                  abscoliq(ig) = &
                        absliq1_d(index,ib) + fint * &
                        (absliq1_d(index+1,ib) - (absliq1_d(index,ib)))
               endif

               taucmc(ig,iplon,lay) = ciwpmc(ig,iplon,lay) * abscoice(ig) + &
                                clwpmc(ig,iplon,lay) * abscoliq(ig)

            endif
         endif
         enddo
      enddo



endif

end subroutine cldprmc_d


attributes(global) subroutine setcoef_d1(ncol,nlay,nlayers, istart)
                            
implicit none

! ----- Input -----

      integer,value :: ncol
      integer,value :: nlay
      integer,value:: nlayers         ! total number of layers
      integer,value:: istart          ! beginning band of calculation
      !integer,value :: nbndlw 
    
      

! ----- Local -----
      integer :: iplon
      integer :: indbound, indlev0
      integer :: lay, indlay, indlev, iband
      integer :: jp1
      real(kind=r8) :: stpfac, tbndfrac, t0frac, tlayfrac, tlevfrac
      real(kind=r8) :: dbdtlev, dbdtlay
      real(kind=r8) :: plog, fp, ft, ft1, water, scalefac, factor, compfp

    

 iplon=(blockIdx%x-1)*blockDim%x+threadIdx%x
   lay=(blockIdx%y-1)*blockDim%y+threadIdx%y  


      stpfac = 296._r8/1013._r8
if  ((iplon>=1 .and. iplon<=ncol).and.(lay==1)) then
      indbound = tbound(iplon) - 159._r8
      if (indbound .lt. 1) then
         indbound = 1
      elseif (indbound .gt. 180) then
         indbound = 180
      endif
      tbndfrac = tbound(iplon) - 159._r8 - float(indbound)
      indlev0 = tz(iplon,0) - 159._r8
      if (indlev0 .lt. 1) then
         indlev0 = 1
      elseif (indlev0 .gt. 180) then
         indlev0 = 180
      endif
      t0frac = tz(iplon,0) - 159._r8 - float(indlev0)
     
end if


if  ((iplon>=1 .and. iplon<=ncol) .and. (lay>=1 .and. lay<=nlayers))then
         indlay = tavel(iplon,lay) - 159._r8
         if (indlay .lt. 1) then
            indlay = 1
         elseif (indlay .gt. 180) then
            indlay = 180
         endif
         tlayfrac = tavel(iplon,lay) - 159._r8 - float(indlay)
         indlev = tz(iplon,lay) - 159._r8
         if (indlev .lt. 1) then
            indlev = 1
         elseif (indlev .gt. 180) then
            indlev = 180
         endif
         tlevfrac = tz(iplon,lay) - 159._r8 - float(indlev)

! Begin spectral band loop

         do iband = 1, 15
  if  ((iplon>=1 .and. iplon<=ncol) .and. (lay==1))then

           
               dbdtlev = totplnk_d(indbound+1,iband) - totplnk_d(indbound,iband)
               plankbnd(iplon,iband) = semiss(iplon,iband) * &
                   (totplnk_d(indbound,iband) + tbndfrac * dbdtlev)
               dbdtlev = totplnk_d(indlev0+1,iband)-totplnk_d(indlev0,iband)
               planklev(iplon,0,iband) = totplnk_d(indlev0,iband) + t0frac * dbdtlev
   endif

            dbdtlev = totplnk_d(indlev+1,iband) - totplnk_d(indlev,iband)
            dbdtlay = totplnk_d(indlay+1,iband) - totplnk_d(indlay,iband)
            planklay(iplon,lay,iband) = totplnk_d(indlay,iband) + tlayfrac * dbdtlay
            planklev(iplon,lay,iband) = totplnk_d(indlev,iband) + tlevfrac * dbdtlev
         enddo

         iband = 16
         if (istart .eq. 16) then
 if  ((iplon>=1 .and. iplon<=ncol) .and. (lay==1))then
               dbdtlev = totplk16_d(indbound+1) - totplk16_d(indbound)
               plankbnd(iplon,iband) = semiss(iplon,iband) * &
                    (totplk16_d(indbound) + tbndfrac * dbdtlev)
               dbdtlev = totplnk_d(indlev0+1,iband)-totplnk_d(indlev0,iband)
               planklev(iplon,0,iband) = totplk16_d(indlev0) + &
                    t0frac * dbdtlev
            endif
            dbdtlev = totplk16_d(indlev+1) - totplk16_d(indlev)
            dbdtlay = totplk16_d(indlay+1) - totplk16_d(indlay)
            planklay(iplon,lay,iband) = totplk16_d(indlay) + tlayfrac * dbdtlay
            planklev(iplon,lay,iband) = totplk16_d(indlev) + tlevfrac * dbdtlev
         else
if  ((iplon>=1 .and. iplon<=ncol) .and. (lay==1))then
               dbdtlev = totplnk_d(indbound+1,iband) - totplnk_d(indbound,iband)
               plankbnd(iplon,iband) = semiss(iplon,iband) * &
                    (totplnk_d(indbound,iband) + tbndfrac * dbdtlev)
               dbdtlev = totplnk_d(indlev0+1,iband)-totplnk_d(indlev0,iband)
               planklev(iplon,0,iband) = totplnk_d(indlev0,iband) + t0frac * dbdtlev
            endif
            dbdtlev = totplnk_d(indlev+1,iband) - totplnk_d(indlev,iband)
            dbdtlay = totplnk_d(indlay+1,iband) - totplnk_d(indlay,iband)
            planklay(iplon,lay,iband) = totplnk_d(indlay,iband) + tlayfrac * dbdtlay
            planklev(iplon,lay,iband) = totplnk_d(indlev,iband) + tlevfrac * dbdtlev
         endif


         plog = dlog(pavel(iplon,lay))
         jp(iplon,lay) = int(36._r8 - 5*(plog+0.04_r8))
         if (jp(iplon,lay) .lt. 1) then
            jp(iplon,lay) = 1
         elseif (jp(iplon,lay) .gt. 58) then
            jp(iplon,lay) = 58
         endif
         jp1 = jp(iplon,lay) + 1
         fp = 5._r8 *(preflog_d(jp(iplon,lay)) - plog)


         jt(iplon,lay) = int(3._r8 + (tavel(iplon,lay)-tref_d(jp(iplon,lay)))/15._r8)
         if (jt(iplon,lay) .lt. 1) then
            jt(iplon,lay) = 1
         elseif (jt(iplon,lay) .gt. 4) then
            jt(iplon,lay) = 4
         endif
         ft = ((tavel(iplon,lay)-tref_d(jp(iplon,lay)))/15._r8) - float(jt(iplon,lay)-3)
         jt1(iplon,lay) = int(3._r8 + (tavel(iplon,lay)-tref_d(jp1))/15._r8)
         if (jt1(iplon,lay) .lt. 1) then
            jt1(iplon,lay) = 1
         elseif (jt1(iplon,lay) .gt. 4) then
            jt1(iplon,lay) = 4
         endif
         ft1 = ((tavel(iplon,lay)-tref_d(jp1))/15._r8) - float(jt1(iplon,lay)-3)
         water = wkl(iplon,1,lay)/coldry(iplon,lay)
         scalefac = pavel(iplon,lay) * stpfac / tavel(iplon,lay)

!  If the pressure is less than ~100mb, perform a different
!  set of species interpolations.
         if (plog .le. 4.56_r8) go to 5300
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7
        ! laytrop(iplon) =  laytrop(iplon) + 1

         forfac(iplon,lay) = scalefac / (1.+water)
         factor = (332.0_r8-tavel(iplon,lay))/36.0_r8
         indfor(iplon,lay) = min(2, max(1, int(factor)))
         forfrac(iplon,lay) = factor - float(indfor(iplon,lay))

!  Set up factors needed to separately include the water vapor
!  self-continuum in the calculation of absorption coefficient.
         selffac(iplon,lay) = water * forfac(iplon,lay)
         factor = (tavel(iplon,lay)-188.0_r8)/7.2_r8
         indself(iplon,lay) = min(9, max(1, int(factor)-7))
         selffrac(iplon,lay) = factor - float(indself(iplon,lay) + 7)

!  Set up factors needed to separately include the minor gases
!  in the calculation of absorption coefficient
         scaleminor(iplon,lay) = pavel(iplon,lay)/tavel(iplon,lay)
         scaleminorn2(iplon,lay) = (pavel(iplon,lay)/tavel(iplon,lay)) &
             *(wbrodl(iplon,lay)/(coldry(iplon,lay)+wkl(iplon,1,lay)))
         factor = (tavel(iplon,lay)-180.8_r8)/7.2_r8
         indminor(iplon,lay) = min(18, max(1, int(factor)))
         minorfrac(iplon,lay) = factor - float(indminor(iplon,lay))

!  Setup reference ratio to be used in calculation of binary
!  species parameter in lower atmosphere.
         rat_h2oco2(iplon,lay)=chi_mls_d(1,jp(iplon,lay))/chi_mls_d(2,jp(iplon,lay))
         rat_h2oco2_1(iplon,lay)=chi_mls_d(1,jp(iplon,lay)+1)/chi_mls_d(2,jp(iplon,lay)+1)

         rat_h2oo3(iplon,lay)=chi_mls_d(1,jp(iplon,lay))/chi_mls_d(3,jp(iplon,lay))
         rat_h2oo3_1(iplon,lay)=chi_mls_d(1,jp(iplon,lay)+1)/chi_mls_d(3,jp(iplon,lay)+1)

         rat_h2on2o(iplon,lay)=chi_mls_d(1,jp(iplon,lay))/chi_mls_d(4,jp(iplon,lay))
         rat_h2on2o_1(iplon,lay)=chi_mls_d(1,jp(iplon,lay)+1)/chi_mls_d(4,jp(iplon,lay)+1)

         rat_h2och4(iplon,lay)=chi_mls_d(1,jp(iplon,lay))/chi_mls_d(6,jp(iplon,lay))
         rat_h2och4_1(iplon,lay)=chi_mls_d(1,jp(iplon,lay)+1)/chi_mls_d(6,jp(iplon,lay)+1)

         rat_n2oco2(iplon,lay)=chi_mls_d(4,jp(iplon,lay))/chi_mls_d(2,jp(iplon,lay))
         rat_n2oco2_1(iplon,lay)=chi_mls_d(4,jp(iplon,lay)+1)/chi_mls_d(2,jp(iplon,lay)+1)

!  Calculate needed column amounts.
         colh2o(iplon,lay) = 1.e-20_r8 * wkl(iplon,1,lay)
         colco2(iplon,lay) = 1.e-20_r8 * wkl(iplon,2,lay)
         colo3(iplon,lay) =  1.e-20_r8 * wkl(iplon,3,lay)
         coln2o(iplon,lay) = 1.e-20_r8 * wkl(iplon,4,lay)
         colco(iplon,lay) = 1.e-20_r8 * wkl(iplon,5,lay)
         colch4(iplon,lay) = 1.e-20_r8 * wkl(iplon,6,lay)
         colo2(iplon,lay) = 1.e-20_r8 * wkl(iplon,7,lay)
         if (colco2(iplon,lay) .eq. 0._r8) colco2(iplon,lay) = 1.e-32_r8 * coldry(iplon,lay)
        ! if (colo3(iplon,lay) .eq. 0._r8) colo3(iplon,lay) = 1.e-32_r8 * coldry(iplon,lay)
         if (coln2o(iplon,lay) .eq. 0._r8) coln2o(iplon,lay) = 1.e-32_r8 * coldry(iplon,lay)
         if (colco(iplon,lay) .eq. 0._r8) colco(iplon,lay) = 1.e-32_r8 * coldry(iplon,lay)
         if (colch4(iplon,lay) .eq. 0._r8) colch4(iplon,lay) = 1.e-32_r8 * coldry(iplon,lay)
         colbrd(iplon,lay) = 1.e-20_r8 * wbrodl(iplon,lay)
         go to 5400

!  Above laytrop.
 5300    continue
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
         forfac(iplon,lay) = scalefac / (1.+water)
         factor = (tavel(iplon,lay)-188.0_r8)/36.0_r8
         indfor(iplon,lay) = 3
         forfrac(iplon,lay) = factor - 1.0_r8

!  Set up factors needed to separately include the water vapor
!  self-continuum in the calculation of absorption coefficient.
         selffac(iplon,lay) = water * forfac(iplon,lay)

!  Set up factors needed to separately include the minor gases
!  in the calculation of absorption coefficient
         scaleminor(iplon,lay) = pavel(iplon,lay)/tavel(iplon,lay)
         scaleminorn2(iplon,lay) = (pavel(iplon,lay)/tavel(iplon,lay)) &
             * (wbrodl(iplon,lay)/(coldry(iplon,lay)+wkl(iplon,1,lay)))
         factor = (tavel(iplon,lay)-180.8_r8)/7.2_r8
         indminor(iplon,lay) = min(18, max(1, int(factor)))
         minorfrac(iplon,lay) = factor - float(indminor(iplon,lay))

!  Setup reference ratio to be used in calculation of binary
!  species parameter in upper atmosphere.
         rat_h2oco2(iplon,lay)=chi_mls_d(1,jp(iplon,lay))/chi_mls_d(2,jp(iplon,lay))
         rat_h2oco2_1(iplon,lay)=chi_mls_d(1,jp(iplon,lay)+1)/chi_mls_d(2,jp(iplon,lay)+1)

         rat_o3co2(iplon,lay)=chi_mls_d(3,jp(iplon,lay))/chi_mls_d(2,jp(iplon,lay))
         rat_o3co2_1(iplon,lay)=chi_mls_d(3,jp(iplon,lay)+1)/chi_mls_d(2,jp(iplon,lay)+1)

!  Calculate needed column amounts.
         colh2o(iplon,lay) = 1.e-20_r8 * wkl(iplon,1,lay)
         colco2(iplon,lay) = 1.e-20_r8 * wkl(iplon,2,lay)
         colo3(iplon,lay) = 1.e-20_r8 * wkl(iplon,3,lay)
         coln2o(iplon,lay) = 1.e-20_r8 * wkl(iplon,4,lay)
         colco(iplon,lay) = 1.e-20_r8 * wkl(iplon,5,lay)
         colch4(iplon,lay) = 1.e-20_r8 * wkl(iplon,6,lay)
         colo2(iplon,lay) = 1.e-20_r8 * wkl(iplon,7,lay)
         if (colco2(iplon,lay) .eq. 0._r8) colco2(iplon,lay) = 1.e-32_r8 * coldry(iplon,lay)
         !if (colo3(iplon,lay) .eq. 0._r8) colo3(iplon,lay) =  coldry(iplon,lay)
         if (coln2o(iplon,lay) .eq. 0._r8) coln2o(iplon,lay) = 1.e-32_r8 * coldry(iplon,lay)
         if (colco(iplon,lay)  .eq. 0._r8) colco(iplon,lay) = 1.e-32_r8 * coldry(iplon,lay)
         if (colch4(iplon,lay) .eq. 0._r8) colch4(iplon,lay) = 1.e-32_r8 * coldry(iplon,lay)
         colbrd(iplon,lay) = 1.e-20_r8 * wbrodl(iplon,lay)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

 5400    continue

         compfp = 1. - fp
         fac10(iplon,lay) = compfp * ft
         fac00(iplon,lay) = compfp * (1._r8 - ft)
         fac11(iplon,lay) = fp * ft1
         fac01(iplon,lay) = fp * (1._r8 - ft1)

!  Rescale selffac and forfac for use in taumol
         selffac(iplon,lay) = colh2o(iplon,lay)*selffac(iplon,lay)
         forfac(iplon,lay) = colh2o(iplon,lay)*forfac(iplon,lay)

! End layer loop
endif

    end subroutine setcoef_d1


attributes(global) subroutine setcoef_d2(ncol,nlayers)
implicit none
integer,value :: ncol
   
integer,value:: nlayers 

integer :: iplon,lay
real(kind=r8) :: plog2

iplon=(blockIdx%x-1)*blockDim%x+threadIdx%x
laytrop(iplon) = 0
if  (iplon>=1 .and. iplon<=ncol) then
  laytrop(iplon) = 0
  do lay=1,nlayers
 plog2 = dlog(pavel(iplon,lay))
   if (plog2 .GT. 4.56_r8) then

   laytrop(iplon) =  laytrop(iplon) + 1
  endif
   enddo



endif

end subroutine setcoef_d2

 
attributes(global) subroutine taumol_d(ncol,nlayers,oneminus)

implicit none
! ----- Input -----
      integer,value :: ncol
      real(kind=r8),value :: oneminus
      integer,value :: nlayers 
      integer,parameter :: ngptlw = 140
     !  integer,value :: iaer


 !__________LOCAL____________
      integer :: iplon,lay
      integer :: ig

iplon=(blockIdx%x-1)*blockDim%x+threadIdx%x
  lay=(blockIdx%y-1)*blockDim%y+threadIdx%y   
 
if((iplon>=1 .and. iplon<=ncol).and.(lay>=1 .and. lay<=nlayers)) then
      call taugb1_d(ncol,iplon,lay,nlayers)

      call taugb2_d(ncol,iplon,lay,nlayers)
                   
      call taugb3_d(ncol,iplon,lay,nlayers,oneminus)

      call taugb4_d(ncol,iplon,lay,nlayers,oneminus)

      call taugb5_d(ncol,iplon,lay,nlayers,oneminus)

      call taugb6_d(ncol,iplon,lay,nlayers)
                    
      call taugb7_d(ncol,iplon,lay,nlayers,oneminus)
                   
      call taugb8_d(ncol,iplon,lay,nlayers)

      call taugb9_d(ncol,iplon,lay,nlayers, oneminus)
                   
      call taugb10_d(ncol,iplon,lay,nlayers)
                     
      call taugb11_d(ncol,iplon,lay,nlayers)

      call taugb12_d(ncol,iplon,lay,nlayers,oneminus)

      call taugb13_d(ncol,iplon,lay,nlayers,oneminus)
                   
      call taugb14_d(ncol,iplon,lay,nlayers)
                   
      call taugb15_d(ncol,iplon,lay,nlayers,oneminus)
                 
      call taugb16_d(ncol,iplon,lay,nlayers,oneminus)
       
endif



 end subroutine taumol_d


attributes(global) subroutine taumol_d2(ncol,nlayers,ngptlw)

implicit none
! ----- Input -----
      integer,value :: ncol
      integer,value :: ngptlw
       integer,value :: nlayers

       !__________LOCAL____________
      integer :: iplon,lay,ig 

 iplon=(blockIdx%x-1)*blockDim%x+threadIdx%x
  lay=(blockIdx%y-1)*blockDim%y+threadIdx%y  
  ig=(blockIdx%z-1)*blockDim%z+threadIdx%z 


 

if((iplon>=1 .and. iplon<=ncol).and.(lay>=1 .and. lay<=nlayers).and.(ig>=1 .and. ig<=ngptlw)) then


 taut(iplon,lay,ig) = taug(iplon,lay,ig) + taua(iplon,lay,ngb_d(ig))

endif

 end subroutine taumol_d2


attributes(device)subroutine taugb1_d(ncol,iplon,lay,nlayers)
 implicit none                                

!     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                          (high key - h2o; high minor - n2)
!
!     note: previous versions of rrtm band 1:
!           10-250 cm-1 (low - h2o; high - h2o)
!----------------------------------------------------------------------------
! ------- Declarations -------
! ----- Input -----
      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol    
!------- Local----------
      integer ::  ind0, ind1, inds, indf, indm, ig
      real(kind=r8) :: pp, corradj, scalen2, tauself, taufor, taun2


      integer, parameter :: ng1  = 10

! Lower atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then
     

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(1) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(1) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)
         pp = pavel(iplon,lay)
         corradj =  1.
         if (pp .lt. 250._r8) then
            corradj = 1._r8 - 0.15_r8 * (250._r8-pp) / 154.4_r8
         endif

         scalen2 = colbrd(iplon,lay) * scaleminorn2(iplon,lay)
         do ig = 1, ng1
            tauself = selffac(iplon,lay) * (selfref01_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref01_d(inds+1,ig) - selfref01_d(inds,ig)))
            taufor =  forfac(iplon,lay) * (forref01_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref01_d(indf+1,ig) -  forref01_d(indf,ig)))
            taun2 = scalen2*(ka_mn201_d(indm,ig) + &
                 minorfrac(iplon,lay) * (ka_mn201_d(indm+1,ig) - ka_mn201_d(indm,ig)))
            taug(iplon,lay,ig) = corradj * (colh2o(iplon,lay) * &
                (fac00(iplon,lay) * absa01_d(ind0,ig) + &
                 fac10(iplon,lay) * absa01_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absa01_d(ind1,ig) + &
                 fac11(iplon,lay) * absa01_d(ind1+1,ig)) &
                 + tauself + taufor + taun2)
             fracs(iplon,lay,ig) = fracrefa01_d(ig)
         enddo



endif

! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(1) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(1) + 1
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)
         pp = pavel(iplon,lay)
         corradj =  1._r8 - 0.15_r8 * (pp / 95.6_r8)

         scalen2 = colbrd(iplon,lay) * scaleminorn2(iplon,lay)
         do ig = 1, ng1
            taufor = forfac(iplon,lay) * (forref01_d(indf,ig) + &
                 forfrac(iplon,lay) * (forref01_d(indf+1,ig) - forref01_d(indf,ig)))
            taun2 = scalen2*(kb_mn201_d(indm,ig) + &
                 minorfrac(iplon,lay) * (kb_mn201_d(indm+1,ig) - kb_mn201_d(indm,ig)))
            taug(iplon,lay,ig) = corradj * (colh2o(iplon,lay) * &
                (fac00(iplon,lay) * absb01_d(ind0,ig) + &
                 fac10(iplon,lay) * absb01_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absb01_d(ind1,ig) + &
                 fac11(iplon,lay) * absb01_d(ind1+1,ig)) &
                 + taufor + taun2)
            fracs(iplon,lay,ig) = fracrefb01_d(ig)
         enddo





endif
end subroutine taugb1_d

attributes(device)subroutine taugb2_d(ncol,iplon,lay,nlayers)
 implicit none                                       
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!
!     note: previous version of rrtm band 2:
!           250 - 500 cm-1 (low - h2o; high - h2o)
!----------------------------------------------------------------------------

! ------- Declarations -------
! ----- Input -----
      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 


! Local
      integer :: ind0, ind1, inds, indf, ig
      real(kind=r8) :: pp, corradj, tauself, taufor

      integer, parameter :: ng2  = 12
      integer, parameter :: ngs1  = 10

! Lower atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then 

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(2) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(2) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         pp = pavel(iplon,lay)
         corradj = 1._r8 - .05_r8 * (pp - 100._r8) / 900._r8
         do ig = 1, ng2
            tauself = selffac(iplon,lay) * (selfref02_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref02_d(inds+1,ig) - selfref02_d(inds,ig)))
            taufor =  forfac(iplon,lay) * (forref02_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref02_d(indf+1,ig) - forref02_d(indf,ig)))
            taug(iplon,lay,ngs1+ig) = corradj * (colh2o(iplon,lay) * &
                (fac00(iplon,lay) * absa02_d(ind0,ig) + &
                 fac10(iplon,lay) * absa02_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absa02_d(ind1,ig) + &
                 fac11(iplon,lay) * absa02_d(ind1+1,ig)) &
                 + tauself + taufor)
            fracs(iplon,lay,ngs1+ig) = fracrefa02_d(ig)
         enddo



endif

! Upper atmosphere loop
 if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(2) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(2) + 1
         indf = indfor(iplon,lay)
         do ig = 1, ng2
            taufor =  forfac(iplon,lay) * (forref02_d(indf,ig) + &
                 forfrac(iplon,lay) * (forref02_d(indf+1,ig) - forref02_d(indf,ig)))
            taug(iplon,lay,ngs1+ig) = colh2o(iplon,lay) * &
                (fac00(iplon,lay) * absb02_d(ind0,ig) + &
                 fac10(iplon,lay) * absb02_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absb02_d(ind1,ig) + &
                 fac11(iplon,lay) * absb02_d(ind1+1,ig)) &
                 + taufor
            fracs(iplon,lay,ngs1+ig) = fracrefb02_d(ig)
         enddo



endif
end subroutine taugb2_d

attributes(device)subroutine taugb3_d(ncol,iplon,lay,nlayers,oneminus)
 implicit none                                    
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                           (high key - h2o,co2; high minor - n2o)
!----------------------------------------------------------------------------
! ------- Declarations -------
! ----- Input -----


  
      real(kind=r8),intent(in) :: oneminus

      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 

! Local
      integer :: ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmn2o, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, &
                       fmn2o, fmn2omf, chi_n2o, ratn2o, adjfac, adjcoln2o
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor, n2om1, n2om2, absn2o
      real(kind=r8) :: refrat_planck_a, refrat_planck_b, refrat_m_a, refrat_m_b
      real(kind=r8) :: tau_major, tau_major1


      integer, parameter :: ng3  = 16
      integer, parameter :: ngs2  = 22

!  P = 212.725 mb
      refrat_planck_a = chi_mls_d(1,9)/chi_mls_d(2,9)

!  P = 95.58 mb
      refrat_planck_b = chi_mls_d(1,13)/chi_mls_d(2,13)

!  P = 706.270mb
      refrat_m_a = chi_mls_d(1,3)/chi_mls_d(2,3)

!  P = 95.58 mb
      refrat_m_b = chi_mls_d(1,13)/chi_mls_d(2,13)


! Lower atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then

         speccomb = colh2o(iplon,lay) + rat_h2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(iplon,lay) + rat_h2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mn2o = colh2o(iplon,lay) + refrat_m_a*colco2(iplon,lay)
         specparm_mn2o = colh2o(iplon,lay)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminus) specparm_mn2o = oneminus
         specmult_mn2o = 8._r8*specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0_r8)
         fmn2omf = minorfrac(iplon,lay)*fmn2o

         chi_n2o = coln2o(iplon,lay)/coldry(iplon,lay)
         ratn2o = 1.e20_r8*chi_n2o/chi_mls_d(4,jp(iplon,lay)+1)
         if (ratn2o .gt. 1.5_r8) then
            adjfac = 0.5_r8+(ratn2o-0.5_r8)**0.65_r8
            adjcoln2o = adjfac*chi_mls_d(4,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20_r8
         else
            adjcoln2o = coln2o(iplon,lay)
         endif

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colco2(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(3) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(3) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1._r8 - fs) * fac00(iplon,lay)
            fac010 = (1._r8 - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif
         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1._r8 - fs1) * fac01(iplon,lay)
            fac011 = (1._r8 - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng3
            tauself = selffac(iplon,lay)* (selfref03_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref03_d(inds+1,ig) - selfref03_d(inds,ig)))
            taufor = forfac(iplon,lay) * (forref03_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref03_d(indf+1,ig) - forref03_d(indf,ig)))
            n2om1 = ka_mn2o03_d(jmn2o,indm,ig) + fmn2o * &
                 (ka_mn2o03_d(jmn2o+1,indm,ig) - ka_mn2o03_d(jmn2o,indm,ig))
            n2om2 = ka_mn2o03_d(jmn2o,indm+1,ig) + fmn2o * &
                 (ka_mn2o03_d(jmn2o+1,indm+1,ig) - ka_mn2o03_d(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(iplon,lay) * (n2om2 - n2om1)

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa03_d(ind0,ig) + &
                    fac100 * absa03_d(ind0+1,ig) + &
                    fac200 * absa03_d(ind0+2,ig) + &
                    fac010 * absa03_d(ind0+9,ig) + &
                    fac110 * absa03_d(ind0+10,ig) + &
                    fac210 * absa03_d(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa03_d(ind0-1,ig) + &
                    fac100 * absa03_d(ind0,ig) + &
                    fac000 * absa03_d(ind0+1,ig) + &
                    fac210 * absa03_d(ind0+8,ig) + &
                    fac110 * absa03_d(ind0+9,ig) + &
                    fac010 * absa03_d(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa03_d(ind0,ig) + &
                    fac100 * absa03_d(ind0+1,ig) + &
                    fac010 * absa03_d(ind0+9,ig) + &
                    fac110 * absa03_d(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa03_d(ind1,ig) + &
                    fac101 * absa03_d(ind1+1,ig) + &
                    fac201 * absa03_d(ind1+2,ig) + &
                    fac011 * absa03_d(ind1+9,ig) + &
                    fac111 * absa03_d(ind1+10,ig) + &
                    fac211 * absa03_d(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa03_d(ind1-1,ig) + &
                    fac101 * absa03_d(ind1,ig) + &
                    fac001 * absa03_d(ind1+1,ig) + &
                    fac211 * absa03_d(ind1+8,ig) + &
                    fac111 * absa03_d(ind1+9,ig) + &
                    fac011 * absa03_d(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa03_d(ind1,ig) +  &
                    fac101 * absa03_d(ind1+1,ig) + &
                    fac011 * absa03_d(ind1+9,ig) + &
                    fac111 * absa03_d(ind1+10,ig))
            endif

            taug(iplon,lay,ngs2+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracs(iplon,lay,ngs2+ig) = fracrefa03_d(ig,jpl) + fpl * &
                 (fracrefa03_d(ig,jpl+1)-fracrefa03_d(ig,jpl))
         enddo




endif

! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then

         speccomb = colh2o(iplon,lay) + rat_h2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(iplon,lay) + rat_h2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 4._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         fac000 = (1._r8 - fs) * fac00(iplon,lay)
         fac010 = (1._r8 - fs) * fac10(iplon,lay)
         fac100 = fs * fac00(iplon,lay)
         fac110 = fs * fac10(iplon,lay)
         fac001 = (1._r8 - fs1) * fac01(iplon,lay)
         fac011 = (1._r8 - fs1) * fac11(iplon,lay)
         fac101 = fs1 * fac01(iplon,lay)
         fac111 = fs1 * fac11(iplon,lay)

         speccomb_mn2o = colh2o(iplon,lay) + refrat_m_b*colco2(iplon,lay)
         specparm_mn2o = colh2o(iplon,lay)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminus) specparm_mn2o = oneminus
         specmult_mn2o = 4._r8*specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0_r8)
         fmn2omf = minorfrac(iplon,lay)*fmn2o

         chi_n2o = coln2o(iplon,lay)/coldry(iplon,lay)
         ratn2o = 1.e20*chi_n2o/chi_mls_d(4,jp(iplon,lay)+1)
         if (ratn2o .gt. 1.5_r8) then
            adjfac = 0.5_r8+(ratn2o-0.5_r8)**0.65_r8
            adjcoln2o = adjfac*chi_mls_d(4,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20_r8
         else
            adjcoln2o = coln2o(iplon,lay)
         endif

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_b*colco2(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 4._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(3) + js
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(3) + js1
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         do ig = 1, ng3
            taufor = forfac(iplon,lay) * (forref03_d(indf,ig) + &
                 forfrac(iplon,lay) * (forref03_d(indf+1,ig) - forref03_d(indf,ig)))
            n2om1 = kb_mn2o03_d(jmn2o,indm,ig) + fmn2o * &
                 (kb_mn2o03_d(jmn2o+1,indm,ig)-kb_mn2o03_d(jmn2o,indm,ig))
            n2om2 = kb_mn2o03_d(jmn2o,indm+1,ig) + fmn2o * &
                 (kb_mn2o03_d(jmn2o+1,indm+1,ig)-kb_mn2o03_d(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(iplon,lay) * (n2om2 - n2om1)
            taug(iplon,lay,ngs2+ig) = speccomb * &
                (fac000 * absb03_d(ind0,ig) + &
                fac100 * absb03_d(ind0+1,ig) + &
                fac010 * absb03_d(ind0+5,ig) + &
                fac110 * absb03_d(ind0+6,ig)) &
                + speccomb1 * &
                (fac001 * absb03_d(ind1,ig) +  &
                fac101 * absb03_d(ind1+1,ig) + &
                fac011 * absb03_d(ind1+5,ig) + &
                fac111 * absb03_d(ind1+6,ig))  &
                + taufor &
                + adjcoln2o*absn2o
            fracs(iplon,lay,ngs2+ig) = fracrefb03_d(ig,jpl) + fpl * &
                (fracrefb03_d(ig,jpl+1)-fracrefb03_d(ig,jpl))
         enddo







endif

end subroutine taugb3_d

attributes(device)subroutine taugb4_d(ncol,iplon ,lay,nlayers,oneminus)
  implicit none                                   
!     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
!----------------------------------------------------------------------------
! ------- Declarations -------

! ----- Input -----


    real(kind=r8),intent(in) :: oneminus 
     
  
      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 





! Local
      integer ::  ind0, ind1, inds, indf, ig
      integer :: js, js1, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor
      real(kind=r8) :: refrat_planck_a, refrat_planck_b
      real(kind=r8) :: tau_major, tau_major1


      integer, parameter :: ng4  = 14
      integer, parameter :: ngs3  = 38
     

! P =   142.5940 mb
      refrat_planck_a = chi_mls_d(1,11)/chi_mls_d(2,11)

! P = 95.58350 mb
      refrat_planck_b = chi_mls_d(3,13)/chi_mls_d(2,13)

! Compute the optical depth by interpolating in ln(pressure) and
! temperature, and appropriate species.  Below laytrop, the water
! vapor self-continuum and foreign continuum is interpolated (in temperature)
! separately.

! Lower atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then


         speccomb = colh2o(iplon,lay) + rat_h2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(iplon,lay) + rat_h2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colco2(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(4) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(4) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1._r8 - fs) * fac00(iplon,lay)
            fac010 = (1._r8 - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1._r8 - fs1) * fac01(iplon,lay)
            fac011 = (1._r8 - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng4
            tauself = selffac(iplon,lay)* (selfref04_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref04_d(inds+1,ig) - selfref04_d(inds,ig)))
            taufor =  forfac(iplon,lay) * (forref04_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref04_d(indf+1,ig) - forref04_d(indf,ig)))

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa04_d(ind0,ig) + &
                    fac100 * absa04_d(ind0+1,ig) + &
                    fac200 * absa04_d(ind0+2,ig) + &
                    fac010 * absa04_d(ind0+9,ig) + &
                    fac110 * absa04_d(ind0+10,ig) + &
                    fac210 * absa04_d(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa04_d(ind0-1,ig) + &
                    fac100 * absa04_d(ind0,ig) + &
                    fac000 * absa04_d(ind0+1,ig) + &
                    fac210 * absa04_d(ind0+8,ig) + &
                    fac110 * absa04_d(ind0+9,ig) + &
                    fac010 * absa04_d(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa04_d(ind0,ig) + &
                    fac100 * absa04_d(ind0+1,ig) + &
                    fac010 * absa04_d(ind0+9,ig) + &
                    fac110 * absa04_d(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa04_d(ind1,ig) +  &
                    fac101 * absa04_d(ind1+1,ig) + &
                    fac201 * absa04_d(ind1+2,ig) + &
                    fac011 * absa04_d(ind1+9,ig) + &
                    fac111 * absa04_d(ind1+10,ig) + &
                    fac211 * absa04_d(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa04_d(ind1-1,ig) + &
                    fac101 * absa04_d(ind1,ig) + &
                    fac001 * absa04_d(ind1+1,ig) + &
                    fac211 * absa04_d(ind1+8,ig) + &
                    fac111 * absa04_d(ind1+9,ig) + &
                    fac011 * absa04_d(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa04_d(ind1,ig) + &
                    fac101 * absa04_d(ind1+1,ig) + &
                    fac011 * absa04_d(ind1+9,ig) + &
                    fac111 * absa04_d(ind1+10,ig))
            endif

            taug(iplon,lay,ngs3+ig) = tau_major + tau_major1 &
                 + tauself + taufor
            fracs(iplon,lay,ngs3+ig) = fracrefa04_d(ig,jpl) + fpl * &
                 (fracrefa04_d(ig,jpl+1)-fracrefa04_d(ig,jpl))
         enddo



endif

! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then

         speccomb = colo3(iplon,lay) + rat_o3co2(iplon,lay)*colco2(iplon,lay)
         specparm = colo3(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colo3(iplon,lay) + rat_o3co2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colo3(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 4._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         fac000 = (1._r8 - fs) * fac00(iplon,lay)
         fac010 = (1._r8 - fs) * fac10(iplon,lay)
         fac100 = fs * fac00(iplon,lay)
         fac110 = fs * fac10(iplon,lay)
         fac001 = (1._r8 - fs1) * fac01(iplon,lay)
         fac011 = (1._r8 - fs1) * fac11(iplon,lay)
         fac101 = fs1 * fac01(iplon,lay)
         fac111 = fs1 * fac11(iplon,lay)

         speccomb_planck = colo3(iplon,lay)+refrat_planck_b*colco2(iplon,lay)
         specparm_planck = colo3(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 4._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(4) + js
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(4) + js1

         do ig = 1, ng4
            taug(iplon,lay,ngs3+ig) =  speccomb * &
                (fac000 * absb04_d(ind0,ig) + &
                fac100 * absb04_d(ind0+1,ig) + &
                fac010 * absb04_d(ind0+5,ig) + &
                fac110 * absb04_d(ind0+6,ig)) &
                + speccomb1 * &
                (fac001 * absb04_d(ind1,ig) +  &
                fac101 * absb04_d(ind1+1,ig) + &
                fac011 * absb04_d(ind1+5,ig) + &
                fac111 * absb04_d(ind1+6,ig))
            fracs(iplon,lay,ngs3+ig) = fracrefb04_d(ig,jpl) + fpl * &
                (fracrefb04_d(ig,jpl+1)-fracrefb04_d(ig,jpl))
         enddo

! Empirical modification to code to improve stratospheric cooling rates
! for co2.  Revised to apply weighting for g-point reduction in this band.

         taug(iplon,lay,ngs3+8)=taug(iplon,lay,ngs3+8)*0.92
         taug(iplon,lay,ngs3+9)=taug(iplon,lay,ngs3+9)*0.88
         taug(iplon,lay,ngs3+10)=taug(iplon,lay,ngs3+10)*1.07
         taug(iplon,lay,ngs3+11)=taug(iplon,lay,ngs3+11)*1.1
         taug(iplon,lay,ngs3+12)=taug(iplon,lay,ngs3+12)*0.99
         taug(iplon,lay,ngs3+13)=taug(iplon,lay,ngs3+13)*0.88
         taug(iplon,lay,ngs3+14)=taug(iplon,lay,ngs3+14)*0.943





endif
end subroutine taugb4_d

attributes(device)subroutine taugb5_d(ncol,iplon,lay,nlayers,oneminus)
 implicit none                       
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!                           (high key - o3,co2)
!----------------------------------------------------------------------------

! ------- Declarations -------
! ----- Input -----

     real(kind=r8),intent(in) :: oneminus
     
     
      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 

! Local
      integer ::  ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmo3, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mo3, specparm_mo3, specmult_mo3, fmo3
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor, o3m1, o3m2, abso3
      real(kind=r8) :: refrat_planck_a, refrat_planck_b, refrat_m_a
      real(kind=r8) :: tau_major, tau_major1


      integer, parameter :: ng5  = 16
      integer, parameter :: ngs4  = 52

! P = 473.420 mb
      refrat_planck_a = chi_mls_d(1,5)/chi_mls_d(2,5)

! P = 0.2369 mb
      refrat_planck_b = chi_mls_d(3,43)/chi_mls_d(2,43)

! P = 317.3480
      refrat_m_a = chi_mls_d(1,7)/chi_mls_d(2,7)

! Compute the optical depth by interpolating in ln(pressure) and
! temperature, and appropriate species.  Below laytrop, the
! water vapor self-continuum and foreign continuum is
! interpolated (in temperature) separately.

! Lower atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then


         speccomb = colh2o(iplon,lay) + rat_h2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(iplon,lay) + rat_h2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mo3 = colh2o(iplon,lay) + refrat_m_a*colco2(iplon,lay)
         specparm_mo3 = colh2o(iplon,lay)/speccomb_mo3
         if (specparm_mo3 .ge. oneminus) specparm_mo3 = oneminus
         specmult_mo3 = 8._r8*specparm_mo3
         jmo3 = 1 + int(specmult_mo3)
         fmo3 = mod(specmult_mo3,1.0_r8)

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colco2(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(5) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(5) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1._r8 - fs) * fac00(iplon,lay)
            fac010 = (1._r8 - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1._r8 - fs1) * fac01(iplon,lay)
            fac011 = (1._r8 - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng5
            tauself = selffac(iplon,lay) * (selfref05_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref05_d(inds+1,ig) - selfref05_d(inds,ig)))
            taufor =  forfac(iplon,lay) * (forref05_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref05_d(indf+1,ig) - forref05_d(indf,ig)))
            o3m1 = ka_mo305_d(jmo3,indm,ig) + fmo3 * &
                 (ka_mo305_d(jmo3+1,indm,ig)-ka_mo305_d(jmo3,indm,ig))
            o3m2 = ka_mo305_d(jmo3,indm+1,ig) + fmo3 * &
                 (ka_mo305_d(jmo3+1,indm+1,ig)-ka_mo305_d(jmo3,indm+1,ig))
            abso3 = o3m1 + minorfrac(iplon,lay)*(o3m2-o3m1)

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa05_d(ind0,ig) + &
                    fac100 * absa05_d(ind0+1,ig) + &
                    fac200 * absa05_d(ind0+2,ig) + &
                    fac010 * absa05_d(ind0+9,ig) + &
                    fac110 * absa05_d(ind0+10,ig) + &
                    fac210 * absa05_d(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa05_d(ind0-1,ig) + &
                    fac100 * absa05_d(ind0,ig) + &
                    fac000 * absa05_d(ind0+1,ig) + &
                    fac210 * absa05_d(ind0+8,ig) + &
                    fac110 * absa05_d(ind0+9,ig) + &
                    fac010 * absa05_d(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa05_d(ind0,ig) + &
                    fac100 * absa05_d(ind0+1,ig) + &
                    fac010 * absa05_d(ind0+9,ig) + &
                    fac110 * absa05_d(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa05_d(ind1,ig) + &
                    fac101 * absa05_d(ind1+1,ig) + &
                    fac201 * absa05_d(ind1+2,ig) + &
                    fac011 * absa05_d(ind1+9,ig) + &
                    fac111 * absa05_d(ind1+10,ig) + &
                    fac211 * absa05_d(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa05_d(ind1-1,ig) + &
                    fac101 * absa05_d(ind1,ig) + &
                    fac001 * absa05_d(ind1+1,ig) + &
                    fac211 * absa05_d(ind1+8,ig) + &
                    fac111 * absa05_d(ind1+9,ig) + &
                    fac011 * absa05_d(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa05_d(ind1,ig) + &
                    fac101 * absa05_d(ind1+1,ig) + &
                    fac011 * absa05_d(ind1+9,ig) + &
                    fac111 * absa05_d(ind1+10,ig))
            endif

            taug(iplon,lay,ngs4+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + abso3*colo3(iplon,lay) &
                 + wx(iplon,1,lay) * ccl405_d(ig)
            fracs(iplon,lay,ngs4+ig) = fracrefa05_d(ig,jpl) + fpl * &
                 (fracrefa05_d(ig,jpl+1)-fracrefa05_d(ig,jpl))
         enddo





endif

! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then

         speccomb = colo3(iplon,lay) + rat_o3co2(iplon,lay)*colco2(iplon,lay)
         specparm = colo3(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 4._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colo3(iplon,lay) + rat_o3co2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colo3(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 4._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         fac000 = (1._r8 - fs) * fac00(iplon,lay)
         fac010 = (1._r8 - fs) * fac10(iplon,lay)
         fac100 = fs * fac00(iplon,lay)
         fac110 = fs * fac10(iplon,lay)
         fac001 = (1._r8 - fs1) * fac01(iplon,lay)
         fac011 = (1._r8 - fs1) * fac11(iplon,lay)
         fac101 = fs1 * fac01(iplon,lay)
         fac111 = fs1 * fac11(iplon,lay)

         speccomb_planck = colo3(iplon,lay)+refrat_planck_b*colco2(iplon,lay)
         specparm_planck = colo3(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 4._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(5) + js
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(5) + js1

         do ig = 1, ng5
            taug(iplon,lay,ngs4+ig) = speccomb * &
                (fac000 * absb05_d(ind0,ig) + &
                fac100 * absb05_d(ind0+1,ig) + &
                fac010 * absb05_d(ind0+5,ig) + &
                fac110 * absb05_d(ind0+6,ig)) &
                + speccomb1 * &
                (fac001 * absb05_d(ind1,ig) + &
                fac101 * absb05_d(ind1+1,ig) + &
                fac011 * absb05_d(ind1+5,ig) + &
                fac111 * absb05_d(ind1+6,ig))  &
                + wx(iplon,1,lay) * ccl405_d(ig)
            fracs(iplon,lay,ngs4+ig) = fracrefb05_d(ig,jpl) + fpl * &
                (fracrefb05_d(ig,jpl+1)-fracrefb05_d(ig,jpl))
         enddo






endif

end subroutine taugb5_d

attributes(device)subroutine taugb6_d(ncol,iplon,lay,nlayers)
 implicit none                        
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                           (high key - nothing; high minor - cfc11, cfc12)
!----------------------------------------------------------------------------

! ------- Declarations -------

! ----- Input -----
     integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 

     



! Local
      integer ::  ind0, ind1, inds, indf, indm, ig
      real(kind=r8) :: chi_co2, ratco2, adjfac, adjcolco2
      real(kind=r8) :: tauself, taufor, absco2


      integer, parameter :: ng6  = 8
      integer, parameter :: ngs5  = 68


! Lower atmosphere loop
 if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then


         chi_co2 = colco2(iplon,lay)/(coldry(iplon,lay))
         ratco2 = 1.e20_r8*chi_co2/chi_mls_d(2,jp(iplon,lay)+1)
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 2.0_r8+(ratco2-2.0_r8)**0.77_r8
            adjcolco2 = adjfac*chi_mls_d(2,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20_r8
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(6) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(6) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         do ig = 1, ng6
            tauself = selffac(iplon,lay) * (selfref06_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref06_d(inds+1,ig) - selfref06_d(inds,ig)))
            taufor =  forfac(iplon,lay) * (forref06_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref06_d(indf+1,ig) - forref06_d(indf,ig)))
            absco2 =  (ka_mco206_d(indm,ig) + minorfrac(iplon,lay) * &
                 (ka_mco206_d(indm+1,ig) - ka_mco206_d(indm,ig)))
            taug(iplon,lay,ngs5+ig) = colh2o(iplon,lay) * &
                (fac00(iplon,lay) * absa06_d(ind0,ig) + &
                 fac10(iplon,lay) * absa06_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absa06_d(ind1,ig) +  &
                 fac11(iplon,lay) * absa06_d(ind1+1,ig))  &
                 + tauself + taufor &
                 + adjcolco2 * absco2 &
                 + wx(iplon,2,lay) * cfc11adj06_d(ig) &
                 + wx(iplon,3,lay) * cfc1206_d(ig)
            fracs(iplon,lay,ngs5+ig) = fracrefa06_d(ig)
         enddo




endif

! Upper atmosphere loop
! Nothing important goes on above laytrop in this band.
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then

         do ig = 1, ng6
            taug(iplon,lay,ngs5+ig) = 0.0_r8 &
                 + wx(iplon,2,lay) * cfc11adj06_d(ig) &
                 + wx(iplon,3,lay) * cfc1206_d(ig)
            fracs(iplon,lay,ngs5+ig) = fracrefa06_d(ig)
         enddo






endif



end subroutine taugb6_d

attributes(device)subroutine taugb7_d(ncol,iplon,lay,nlayers,oneminus)
 implicit none                           


!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!                            (high key - o3; high minor - co2)
!----------------------------------------------------------------------------
! ------- Declarations -------

! ----- Input -----

      
     integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 

  
      real(kind=r8),intent(in) :: oneminus 
      

      
      


! Local
      integer :: ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmco2, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor, co2m1, co2m2, absco2
      real(kind=r8) :: chi_co2, ratco2, adjfac, adjcolco2
      real(kind=r8) :: refrat_planck_a, refrat_m_a
      real(kind=r8) :: tau_major, tau_major1

      integer, parameter :: ng7  = 12
      integer, parameter :: ngs6  = 76


! P = 706.2620 mb
      refrat_planck_a = chi_mls_d(1,3)/chi_mls_d(3,3)

! P = 706.2720 mb
      refrat_m_a = chi_mls_d(1,3)/chi_mls_d(3,3)

! Lower atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then

         speccomb = colh2o(iplon,lay) + rat_h2oo3(iplon,lay)*colo3(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(iplon,lay) + rat_h2oo3_1(iplon,lay)*colo3(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mco2 = colh2o(iplon,lay) + refrat_m_a*colo3(iplon,lay)
         specparm_mco2 = colh2o(iplon,lay)/speccomb_mco2
         if (specparm_mco2 .ge. oneminus) specparm_mco2 = oneminus
         specmult_mco2 = 8._r8*specparm_mco2

         jmco2 = 1 + int(specmult_mco2)
         fmco2 = mod(specmult_mco2,1.0_r8)

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor
!  to obtain the proper contribution.
         chi_co2 = colco2(iplon,lay)/(coldry(iplon,lay))
         ratco2 = 1.e20*chi_co2/chi_mls_d(2,jp(iplon,lay)+1)
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 3.0_r8+(ratco2-3.0_r8)**0.79_r8
            adjcolco2 = adjfac*chi_mls_d(2,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20_r8
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colo3(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(7) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(7) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1._r8 - fs) * fac00(iplon,lay)
            fac010 = (1._r8 - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif
         if (specparm .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1._r8 - fs1) * fac01(iplon,lay)
            fac011 = (1._r8 - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng7
            tauself = selffac(iplon,lay)* (selfref07_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref07_d(inds+1,ig) - selfref07_d(inds,ig)))
            taufor = forfac(iplon,lay) * (forref07_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref07_d(indf+1,ig) - forref07_d(indf,ig)))
            co2m1 = ka_mco207_d(jmco2,indm,ig) + fmco2 * &
                 (ka_mco207_d(jmco2+1,indm,ig) - ka_mco207_d(jmco2,indm,ig))
            co2m2 = ka_mco207_d(jmco2,indm+1,ig) + fmco2 * &
                 (ka_mco207_d(jmco2+1,indm+1,ig) - ka_mco207_d(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(iplon,lay) * (co2m2 - co2m1)

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa07_d(ind0,ig) + &
                    fac100 * absa07_d(ind0+1,ig) + &
                    fac200 * absa07_d(ind0+2,ig) + &
                    fac010 * absa07_d(ind0+9,ig) + &
                    fac110 * absa07_d(ind0+10,ig) + &
                    fac210 * absa07_d(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa07_d(ind0-1,ig) + &
                    fac100 * absa07_d(ind0,ig) + &
                    fac000 * absa07_d(ind0+1,ig) + &
                    fac210 * absa07_d(ind0+8,ig) + &
                    fac110 * absa07_d(ind0+9,ig) + &
                    fac010 * absa07_d(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa07_d(ind0,ig) + &
                    fac100 * absa07_d(ind0+1,ig) + &
                    fac010 * absa07_d(ind0+9,ig) + &
                    fac110 * absa07_d(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa07_d(ind1,ig) + &
                    fac101 * absa07_d(ind1+1,ig) + &
                    fac201 * absa07_d(ind1+2,ig) + &
                    fac011 * absa07_d(ind1+9,ig) + &
                    fac111 * absa07_d(ind1+10,ig) + &
                    fac211 * absa07_d(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa07_d(ind1-1,ig) + &
                    fac101 * absa07_d(ind1,ig) + &
                    fac001 * absa07_d(ind1+1,ig) + &
                    fac211 * absa07_d(ind1+8,ig) + &
                    fac111 * absa07_d(ind1+9,ig) + &
                    fac011 * absa07_d(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa07_d(ind1,ig) +  &
                    fac101 * absa07_d(ind1+1,ig) + &
                    fac011 * absa07_d(ind1+9,ig) + &
                    fac111 * absa07_d(ind1+10,ig))
            endif

            taug(iplon,lay,ngs6+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcolco2*absco2
            fracs(iplon,lay,ngs6+ig) = fracrefa07_d(ig,jpl) + fpl * &
                 (fracrefa07_d(ig,jpl+1)-fracrefa07_d(ig,jpl))
         enddo





endif

! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then   

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor
!  to obtain the proper contribution.
         chi_co2 = colco2(iplon,lay)/(coldry(iplon,lay))
         ratco2 = 1.e20*chi_co2/chi_mls_d(2,jp(iplon,lay)+1)
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 2.0_r8+(ratco2-2.0_r8)**0.79_r8
            adjcolco2 = adjfac*chi_mls_d(2,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20_r8
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(7) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(7) + 1
         indm = indminor(iplon,lay)

         do ig = 1, ng7
            absco2 = kb_mco207_d(indm,ig) + minorfrac(iplon,lay) * &
                 (kb_mco207_d(indm+1,ig) - kb_mco207_d(indm,ig))
            taug(iplon,lay,ngs6+ig) = colo3(iplon,lay) * &
                 (fac00(iplon,lay) * absb07_d(ind0,ig) + &
                 fac10(iplon,lay) * absb07_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absb07_d(ind1,ig) + &
                 fac11(iplon,lay) * absb07_d(ind1+1,ig)) &
                 + adjcolco2 * absco2
            fracs(iplon,lay,ngs6+ig) = fracrefb07_d(ig)
         enddo

! Empirical modification to code to improve stratospheric cooling rates
! for o3.  Revised to apply weighting for g-point reduction in this band.

         taug(iplon,lay,ngs6+6)=taug(iplon,lay,ngs6+6)*0.92_r8
         taug(iplon,lay,ngs6+7)=taug(iplon,lay,ngs6+7)*0.88_r8
         taug(iplon,lay,ngs6+8)=taug(iplon,lay,ngs6+8)*1.07_r8
         taug(iplon,lay,ngs6+9)=taug(iplon,lay,ngs6+9)*1.1_r8
         taug(iplon,lay,ngs6+10)=taug(iplon,lay,ngs6+10)*0.99_r8
         taug(iplon,lay,ngs6+11)=taug(iplon,lay,ngs6+11)*0.855_r8





endif
end subroutine taugb7_d

attributes(device)subroutine taugb8_d(ncol,iplon,lay,nlayers)
 implicit none                          

!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!                             (high key - o3; high minor - co2, n2o)
!----------------------------------------------------------------------------
! ------- Declarations -------

! ----- Input -----
      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 

    
      
     

! Local
      integer ::  ind0, ind1, inds, indf, indm, ig
      real(kind=r8) :: tauself, taufor, absco2, abso3, absn2o
      real(kind=r8) :: chi_co2, ratco2, adjfac, adjcolco2

      integer, parameter :: ng8  = 8
      integer, parameter :: ngs7  = 88


! Lower atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then


         chi_co2 = colco2(iplon,lay)/(coldry(iplon,lay))
         ratco2 = 1.e20_r8*chi_co2/chi_mls_d(2,jp(iplon,lay)+1)
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 2.0_r8+(ratco2-2.0_r8)**0.65_r8
            adjcolco2 = adjfac*chi_mls_d(2,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20_r8
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(8) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(8) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         do ig = 1, ng8
            tauself = selffac(iplon,lay) * (selfref08_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref08_d(inds+1,ig) - selfref08_d(inds,ig)))
            taufor = forfac(iplon,lay) * (forref08_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref08_d(indf+1,ig) - forref08_d(indf,ig)))
            absco2 =  (ka_mco208_d(indm,ig) + minorfrac(iplon,lay) * &
                 (ka_mco208_d(indm+1,ig) - ka_mco208_d(indm,ig)))
            abso3 =  (ka_mo308_d(indm,ig) + minorfrac(iplon,lay) * &
                 (ka_mo308_d(indm+1,ig) - ka_mo308_d(indm,ig)))
            absn2o =  (ka_mn2o08_d(indm,ig) + minorfrac(iplon,lay) * &
                 (ka_mn2o08_d(indm+1,ig) - ka_mn2o08_d(indm,ig)))
            taug(iplon,lay,ngs7+ig) = colh2o(iplon,lay) * &
                 (fac00(iplon,lay) * absa08_d(ind0,ig) + &
                 fac10(iplon,lay) * absa08_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absa08_d(ind1,ig) +  &
                 fac11(iplon,lay) * absa08_d(ind1+1,ig)) &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colo3(iplon,lay) * abso3 &
                 + coln2o(iplon,lay) * absn2o &
                 + wx(iplon,3,lay) * cfc1208_d(ig) &
                 + wx(iplon,4,lay) * cfc22adj08_d(ig)
            fracs(iplon,lay,ngs7+ig) = fracrefa08_d(ig)
         enddo




endif

! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor
!  to obtain the proper contribution.
         chi_co2 = colco2(iplon,lay)/coldry(iplon,lay)
         ratco2 = 1.e20_r8*chi_co2/chi_mls_d(2,jp(iplon,lay)+1)
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 2.0_r8+(ratco2-2.0_r8)**0.65_r8
            adjcolco2 = adjfac*chi_mls_d(2,jp(iplon,lay)+1) * coldry(iplon,lay)*1.e-20_r8
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(8) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(8) + 1
         indm = indminor(iplon,lay)

         do ig = 1, ng8
            absco2 =  (kb_mco208_d(indm,ig) + minorfrac(iplon,lay) * &
                 (kb_mco208_d(indm+1,ig) - kb_mco208_d(indm,ig)))
            absn2o =  (kb_mn2o08_d(indm,ig) + minorfrac(iplon,lay) * &
                 (kb_mn2o08_d(indm+1,ig) - kb_mn2o08_d(indm,ig)))
            taug(iplon,lay,ngs7+ig) = colo3(iplon,lay) * &
                 (fac00(iplon,lay) * absb08_d(ind0,ig) + &
                 fac10(iplon,lay) * absb08_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absb08_d(ind1,ig) + &
                 fac11(iplon,lay) * absb08_d(ind1+1,ig)) &
                 + adjcolco2*absco2 &
                 + coln2o(iplon,lay)*absn2o &
                 + wx(iplon,3,lay) * cfc1208_d(ig) &
                 + wx(iplon,4,lay) * cfc22adj08_d(ig)
            fracs(iplon,lay,ngs7+ig) = fracrefb08_d(ig)
         enddo






endif
end subroutine taugb8_d

attributes(device)subroutine taugb9_d(ncol,iplon,lay,nlayers, oneminus)
implicit none                           
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!                             (high key - ch4; high minor - n2o)
!----------------------------------------------------------------------------
! ------- Declarations -------

! ----- Input -----
      real(kind=r8),intent(in) :: oneminus
   

      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 


! Local
      integer :: ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmn2o, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, fmn2o
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor, n2om1, n2om2, absn2o
      real(kind=r8) :: chi_n2o, ratn2o, adjfac, adjcoln2o
      real(kind=r8) :: refrat_planck_a, refrat_m_a
      real(kind=r8) :: tau_major, tau_major1


      integer, parameter :: ng9  = 12
      integer, parameter :: ngs8  = 96

! P = 212 mb
      refrat_planck_a = chi_mls_d(1,9)/chi_mls_d(6,9)

! P = 706.272 mb
      refrat_m_a = chi_mls_d(1,3)/chi_mls_d(6,3)

! Lower atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then

         speccomb = colh2o(iplon,lay) + rat_h2och4(iplon,lay)*colch4(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(iplon,lay) + rat_h2och4_1(iplon,lay)*colch4(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mn2o = colh2o(iplon,lay) + refrat_m_a*colch4(iplon,lay)
         specparm_mn2o = colh2o(iplon,lay)/speccomb_mn2o
         if (specparm_mn2o .ge. oneminus) specparm_mn2o = oneminus
         specmult_mn2o = 8._r8*specparm_mn2o
         jmn2o = 1 + int(specmult_mn2o)
         fmn2o = mod(specmult_mn2o,1.0_r8)

!  In atmospheres where the amount of N2O is too great to be considered
!  a minor species, adjust the column amount of N2O by an empirical factor
!  to obtain the proper contribution.
         chi_n2o = coln2o(iplon,lay)/(coldry(iplon,lay))
         ratn2o = 1.e20_r8*chi_n2o/chi_mls_d(4,jp(iplon,lay)+1)
         if (ratn2o .gt. 1.5_r8) then
            adjfac = 0.5_r8+(ratn2o-0.5_r8)**0.65_r8
            adjcoln2o = adjfac*chi_mls_d(4,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20_r8
         else
            adjcoln2o = coln2o(iplon,lay)
         endif

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colch4(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(9) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(9) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1._r8 - fs) * fac00(iplon,lay)
            fac010 = (1._r8 - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1._r8 - fs1) * fac01(iplon,lay)
            fac011 = (1._r8 - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng9
            tauself = selffac(iplon,lay)* (selfref09_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref09_d(inds+1,ig) - selfref09_d(inds,ig)))
            taufor = forfac(iplon,lay) * (forref09_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref09_d(indf+1,ig) - forref09_d(indf,ig)))
            n2om1 = ka_mn2o09_d(jmn2o,indm,ig) + fmn2o * &
                 (ka_mn2o09_d(jmn2o+1,indm,ig) - ka_mn2o09_d(jmn2o,indm,ig))
            n2om2 = ka_mn2o09_d(jmn2o,indm+1,ig) + fmn2o * &
                 (ka_mn2o09_d(jmn2o+1,indm+1,ig) - ka_mn2o09_d(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(iplon,lay) * (n2om2 - n2om1)

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa09_d(ind0,ig) + &
                    fac100 * absa09_d(ind0+1,ig) + &
                    fac200 * absa09_d(ind0+2,ig) + &
                    fac010 * absa09_d(ind0+9,ig) + &
                    fac110 * absa09_d(ind0+10,ig) + &
                    fac210 * absa09_d(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa09_d(ind0-1,ig) + &
                    fac100 * absa09_d(ind0,ig) + &
                    fac000 * absa09_d(ind0+1,ig) + &
                    fac210 * absa09_d(ind0+8,ig) + &
                    fac110 * absa09_d(ind0+9,ig) + &
                    fac010 * absa09_d(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa09_d(ind0,ig) + &
                    fac100 * absa09_d(ind0+1,ig) + &
                    fac010 * absa09_d(ind0+9,ig) + &
                    fac110 * absa09_d(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa09_d(ind1,ig) + &
                    fac101 * absa09_d(ind1+1,ig) + &
                    fac201 * absa09_d(ind1+2,ig) + &
                    fac011 * absa09_d(ind1+9,ig) + &
                    fac111 * absa09_d(ind1+10,ig) + &
                    fac211 * absa09_d(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa09_d(ind1-1,ig) + &
                    fac101 * absa09_d(ind1,ig) + &
                    fac001 * absa09_d(ind1+1,ig) + &
                    fac211 * absa09_d(ind1+8,ig) + &
                    fac111 * absa09_d(ind1+9,ig) + &
                    fac011 * absa09_d(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa09_d(ind1,ig) + &
                    fac101 * absa09_d(ind1+1,ig) + &
                    fac011 * absa09_d(ind1+9,ig) + &
                    fac111 * absa09_d(ind1+10,ig))
            endif

            taug(iplon,lay,ngs8+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracs(iplon,lay,ngs8+ig) = fracrefa09_d(ig,jpl) + fpl * &
                 (fracrefa09_d(ig,jpl+1)-fracrefa09_d(ig,jpl))
         enddo





endif

! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then

         chi_n2o = coln2o(iplon,lay)/(coldry(iplon,lay))
         ratn2o = 1.e20_r8*chi_n2o/chi_mls_d(4,jp(iplon,lay)+1)
         if (ratn2o .gt. 1.5_r8) then
            adjfac = 0.5_r8+(ratn2o-0.5_r8)**0.65_r8
            adjcoln2o = adjfac*chi_mls_d(4,jp(iplon,lay)+1)*coldry(iplon,lay)*1.e-20_r8
         else
            adjcoln2o = coln2o(iplon,lay)
         endif

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(9) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(9) + 1
         indm = indminor(iplon,lay)

         do ig = 1, ng9
            absn2o = kb_mn2o09_d(indm,ig) + minorfrac(iplon,lay) * &
                (kb_mn2o09_d(indm+1,ig) - kb_mn2o09_d(indm,ig))
            taug(iplon,lay,ngs8+ig) = colch4(iplon,lay) * &
                 (fac00(iplon,lay) * absb09_d(ind0,ig) + &
                 fac10(iplon,lay) * absb09_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absb09_d(ind1,ig) +  &
                 fac11(iplon,lay) * absb09_d(ind1+1,ig)) &
                 + adjcoln2o*absn2o
            fracs(iplon,lay,ngs8+ig) = fracrefb09_d(ig)
         enddo





endif

end subroutine taugb9_d

attributes(device)subroutine taugb10_d(ncol,iplon,lay,nlayers)
implicit none                                   
!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
! ------- Declarations -------

! ----- Input -----
      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 


  
 
! Local
      integer :: ind0, ind1, inds, indf, ig
      real(kind=r8) :: tauself, taufor

      integer, parameter :: ng10 = 6
      integer, parameter :: ngs9  = 108

! Lower atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(10) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(10) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)

         do ig = 1, ng10
            tauself = selffac(iplon,lay) * (selfref10_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref10_d(inds+1,ig) - selfref10_d(inds,ig)))
            taufor = forfac(iplon,lay) * (forref10_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref10_d(indf+1,ig) - forref10_d(indf,ig)))
            taug(iplon,lay,ngs9+ig) = colh2o(iplon,lay) * &
                 (fac00(iplon,lay) * absa10_d(ind0,ig) + &
                 fac10(iplon,lay) * absa10_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absa10_d(ind1,ig) + &
                 fac11(iplon,lay) * absa10_d(ind1+1,ig))  &
                 + tauself + taufor
            fracs(iplon,lay,ngs9+ig) = fracrefa10_d(ig)
         enddo




endif

! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(10) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(10) + 1
         indf = indfor(iplon,lay)

         do ig = 1, ng10
            taufor = forfac(iplon,lay) * (forref10_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref10_d(indf+1,ig) - forref10_d(indf,ig)))
            taug(iplon,lay,ngs9+ig) = colh2o(iplon,lay) * &
                 (fac00(iplon,lay) * absb10_d(ind0,ig) + &
                 fac10(iplon,lay) * absb10_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absb10_d(ind1,ig) +  &
                 fac11(iplon,lay) * absb10_d(ind1+1,ig)) &
                 + taufor
            fracs(iplon,lay,ngs9+ig) = fracrefb10_d(ig)
         enddo




endif
end subroutine taugb10_d

attributes(device)subroutine taugb11_d(ncol,iplon,lay,nlayers)
 implicit none                                   

!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!----------------------------------------------------------------------------
! ------- Declarations -------

! ----- Input -----
      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 

     

! Local
      integer :: ind0, ind1, inds, indf, indm, ig
      real(kind=r8) :: scaleo2, tauself, taufor, tauo2

      integer, parameter :: ng11 = 8
      integer, parameter :: ngs10 = 114

! Lower atmosphere loop
  if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then 

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(11) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(11) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)
         scaleo2 = colo2(iplon,lay)*scaleminor(iplon,lay)
         do ig = 1, ng11
            tauself = selffac(iplon,lay) * (selfref11_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref11_d(inds+1,ig) - selfref11_d(inds,ig)))
            taufor = forfac(iplon,lay) * (forref11_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref11_d(indf+1,ig) - forref11_d(indf,ig)))
            tauo2 =  scaleo2 * (ka_mo211_d(indm,ig) + minorfrac(iplon,lay) * &
                 (ka_mo211_d(indm+1,ig) - ka_mo211_d(indm,ig)))
            taug(iplon,lay,ngs10+ig) = colh2o(iplon,lay) * &
                 (fac00(iplon,lay) * absa11_d(ind0,ig) + &
                 fac10(iplon,lay) * absa11_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absa11_d(ind1,ig) + &
                 fac11(iplon,lay) * absa11_d(ind1+1,ig)) &
                 + tauself + taufor &
                 + tauo2
            fracs(iplon,lay,ngs10+ig) = fracrefa11_d(ig)
         enddo




endif

! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then

         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(11) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(11) + 1
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)
         scaleo2 = colo2(iplon,lay)*scaleminor(iplon,lay)
         do ig = 1, ng11
            taufor = forfac(iplon,lay) * (forref11_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref11_d(indf+1,ig) - forref11_d(indf,ig)))
            tauo2 =  scaleo2 * (kb_mo211_d(indm,ig) + minorfrac(iplon,lay) * &
                 (kb_mo211_d(indm+1,ig) - kb_mo211_d(indm,ig)))
            taug(iplon,lay,ngs10+ig) = colh2o(iplon,lay) * &
                 (fac00(iplon,lay) * absb11_d(ind0,ig) + &
                 fac10(iplon,lay) * absb11_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absb11_d(ind1,ig) + &
                 fac11(iplon,lay) * absb11_d(ind1+1,ig))  &
                 + taufor &
                 + tauo2
            fracs(iplon,lay,ngs10+ig) = fracrefb11_d(ig)
         enddo




endif
end subroutine taugb11_d

attributes(device)subroutine taugb12_d(ncol,iplon,lay,nlayers,oneminus)
 implicit none                                    
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!----------------------------------------------------------------------------
! ------- Declarations -------

! ----- Input -----
      
       real(kind=r8),intent(in) :: oneminus
       

       integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 

! Local
      integer :: ind0, ind1, inds, indf, ig
      integer :: js, js1, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor
      real(kind=r8) :: refrat_planck_a
      real(kind=r8) :: tau_major, tau_major1



      integer, parameter :: ng12 = 8
      integer, parameter :: ngs11 = 122

! P =   174.164 mb
      refrat_planck_a = chi_mls_d(1,10)/chi_mls_d(2,10)



! Lower atmosphere loop

if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then

         speccomb = colh2o(iplon,lay) + rat_h2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(iplon,lay) + rat_h2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colco2(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(12) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(12) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1._r8 - fs) * fac00(iplon,lay)
            fac010 = (1._r8 - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1._r8 - fs1) * fac01(iplon,lay)
            fac011 = (1._r8 - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng12
            tauself = selffac(iplon,lay)* (selfref12_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref12_d(inds+1,ig) - selfref12_d(inds,ig)))
            taufor = forfac(iplon,lay) * (forref12_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref12_d(indf+1,ig) - forref12_d(indf,ig)))

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa12_d(ind0,ig) + &
                    fac100 * absa12_d(ind0+1,ig) + &
                    fac200 * absa12_d(ind0+2,ig) + &
                    fac010 * absa12_d(ind0+9,ig) + &
                    fac110 * absa12_d(ind0+10,ig) + &
                    fac210 * absa12_d(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa12_d(ind0-1,ig) + &
                    fac100 * absa12_d(ind0,ig) + &
                    fac000 * absa12_d(ind0+1,ig) + &
                    fac210 * absa12_d(ind0+8,ig) + &
                    fac110 * absa12_d(ind0+9,ig) + &
                    fac010 * absa12_d(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa12_d(ind0,ig) + &
                    fac100 * absa12_d(ind0+1,ig) + &
                    fac010 * absa12_d(ind0+9,ig) + &
                    fac110 * absa12_d(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa12_d(ind1,ig) + &
                    fac101 * absa12_d(ind1+1,ig) + &
                    fac201 * absa12_d(ind1+2,ig) + &
                    fac011 * absa12_d(ind1+9,ig) + &
                    fac111 * absa12_d(ind1+10,ig) + &
                    fac211 * absa12_d(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa12_d(ind1-1,ig) + &
                    fac101 * absa12_d(ind1,ig) + &
                    fac001 * absa12_d(ind1+1,ig) + &
                    fac211 * absa12_d(ind1+8,ig) + &
                    fac111 * absa12_d(ind1+9,ig) + &
                    fac011 * absa12_d(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa12_d(ind1,ig) + &
                    fac101 * absa12_d(ind1+1,ig) + &
                    fac011 * absa12_d(ind1+9,ig) + &
                    fac111 * absa12_d(ind1+10,ig))
            endif

            taug(iplon,lay,ngs11+ig) = tau_major + tau_major1 &
                 + tauself + taufor
            fracs(iplon,lay,ngs11+ig) = fracrefa12_d(ig,jpl) + fpl * &
                 (fracrefa12_d(ig,jpl+1)-fracrefa12_d(ig,jpl))
         enddo




endif

! Upper atmosphere loop
    
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then


         do ig = 1, ng12
            taug(iplon,lay,ngs11+ig) = 0.0_r8
            fracs(iplon,lay,ngs11+ig) = 0.0_r8
         enddo



endif
end subroutine taugb12_d

attributes(device)subroutine taugb13_d(ncol,iplon,lay,nlayers,oneminus)
 implicit none                                   
!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
!----------------------------------------------------------------------------

! ----- Input -----
      
      real(kind=r8),intent(in) :: oneminus
     

        integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 

     

! ------- Declarations -------

! Local
      integer :: ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmco2, jmco, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real(kind=r8) :: speccomb_mco, specparm_mco, specmult_mco, fmco
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor, co2m1, co2m2, absco2
      real(kind=r8) :: com1, com2, absco, abso3
      real(kind=r8) :: chi_co2, ratco2, adjfac, adjcolco2
      real(kind=r8) :: refrat_planck_a, refrat_m_a, refrat_m_a3
      real(kind=r8) :: tau_major, tau_major1


      integer, parameter :: ng13 = 4
      integer, parameter :: ngs12 = 130

! P = 473.420 mb (Level 5)
      refrat_planck_a = chi_mls_d(1,5)/chi_mls_d(4,5)

! P = 1053. (Level 1)
      refrat_m_a = chi_mls_d(1,1)/chi_mls_d(4,1)

! P = 706. (Level 3)
      refrat_m_a3 = chi_mls_d(1,3)/chi_mls_d(4,3)

! Lower atmosphere loop

if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then

         speccomb = colh2o(iplon,lay) + rat_h2on2o(iplon,lay)*coln2o(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(iplon,lay) + rat_h2on2o_1(iplon,lay)*coln2o(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mco2 = colh2o(iplon,lay) + refrat_m_a*coln2o(iplon,lay)
         specparm_mco2 = colh2o(iplon,lay)/speccomb_mco2
         if (specparm_mco2 .ge. oneminus) specparm_mco2 = oneminus
         specmult_mco2 = 8._r8*specparm_mco2
         jmco2 = 1 + int(specmult_mco2)
         fmco2 = mod(specmult_mco2,1.0_r8)

!  In atmospheres where the amount of CO2 is too great to be considered
!  a minor species, adjust the column amount of CO2 by an empirical factor
!  to obtain the proper contribution.
         chi_co2 = colco2(iplon,lay)/(coldry(iplon,lay))
         ratco2 = 1.e20_r8*chi_co2/3.55e-4_r8
         if (ratco2 .gt. 3.0_r8) then
            adjfac = 2.0_r8+(ratco2-2.0_r8)**0.68_r8
            adjcolco2 = adjfac*3.55e-4*coldry(iplon,lay)*1.e-20_r8
         else
            adjcolco2 = colco2(iplon,lay)
         endif

         speccomb_mco = colh2o(iplon,lay) + refrat_m_a3*coln2o(iplon,lay)
         specparm_mco = colh2o(iplon,lay)/speccomb_mco
         if (specparm_mco .ge. oneminus) specparm_mco = oneminus
         specmult_mco = 8._r8*specparm_mco
         jmco = 1 + int(specmult_mco)
         fmco = mod(specmult_mco,1.0_r8)

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*coln2o(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(13) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(13) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1._r8 - fs) * fac00(iplon,lay)
            fac010 = (1._r8 - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1._r8 - fs1) * fac01(iplon,lay)
            fac011 = (1._r8 - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng13
            tauself = selffac(iplon,lay)* (selfref13_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref13_d(inds+1,ig) - selfref13_d(inds,ig)))
            taufor = forfac(iplon,lay) * (forref13_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref13_d(indf+1,ig) - forref13_d(indf,ig)))
            co2m1 = ka_mco213_d(jmco2,indm,ig) + fmco2 * &
                 (ka_mco213_d(jmco2+1,indm,ig) - ka_mco213_d(jmco2,indm,ig))
            co2m2 = ka_mco213_d(jmco2,indm+1,ig) + fmco2 * &
                 (ka_mco213_d(jmco2+1,indm+1,ig) - ka_mco213_d(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(iplon,lay) * (co2m2 - co2m1)
            com1 = ka_mco13_d(jmco,indm,ig) + fmco * &
                 (ka_mco13_d(jmco+1,indm,ig) - ka_mco13_d(jmco,indm,ig))
            com2 = ka_mco13_d(jmco,indm+1,ig) + fmco * &
                 (ka_mco13_d(jmco+1,indm+1,ig) - ka_mco13_d(jmco,indm+1,ig))
            absco = com1 + minorfrac(iplon,lay) * (com2 - com1)

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa13_d(ind0,ig) + &
                    fac100 * absa13_d(ind0+1,ig) + &
                    fac200 * absa13_d(ind0+2,ig) + &
                    fac010 * absa13_d(ind0+9,ig) + &
                    fac110 * absa13_d(ind0+10,ig) + &
                    fac210 * absa13_d(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa13_d(ind0-1,ig) + &
                    fac100 * absa13_d(ind0,ig) + &
                    fac000 * absa13_d(ind0+1,ig) + &
                    fac210 * absa13_d(ind0+8,ig) + &
                    fac110 * absa13_d(ind0+9,ig) + &
                    fac010 * absa13_d(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa13_d(ind0,ig) + &
                    fac100 * absa13_d(ind0+1,ig) + &
                    fac010 * absa13_d(ind0+9,ig) + &
                    fac110 * absa13_d(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa13_d(ind1,ig) + &
                    fac101 * absa13_d(ind1+1,ig) + &
                    fac201 * absa13_d(ind1+2,ig) + &
                    fac011 * absa13_d(ind1+9,ig) + &
                    fac111 * absa13_d(ind1+10,ig) + &
                    fac211 * absa13_d(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa13_d(ind1-1,ig) + &
                    fac101 * absa13_d(ind1,ig) + &
                    fac001 * absa13_d(ind1+1,ig) + &
                    fac211 * absa13_d(ind1+8,ig) + &
                    fac111 * absa13_d(ind1+9,ig) + &
                    fac011 * absa13_d(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa13_d(ind1,ig) + &
                    fac101 * absa13_d(ind1+1,ig) + &
                    fac011 * absa13_d(ind1+9,ig) + &
                    fac111 * absa13_d(ind1+10,ig))
            endif

            taug(iplon,lay,ngs12+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colco(iplon,lay)*absco
            fracs(iplon,lay,ngs12+ig) = fracrefa13_d(ig,jpl) + fpl * &
                 (fracrefa13_d(ig,jpl+1)-fracrefa13_d(ig,jpl))
         enddo
  


  endif

! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then

         indm = indminor(iplon,lay)
         do ig = 1, ng13
            abso3 = kb_mo313_d(indm,ig) + minorfrac(iplon,lay) * &
                 (kb_mo313_d(indm+1,ig) - kb_mo313_d(indm,ig))
            taug(iplon,lay,ngs12+ig) = colo3(iplon,lay)*abso3
            fracs(iplon,lay,ngs12+ig) =  fracrefb13_d(ig)
         enddo





endif
end subroutine taugb13_d

attributes(device)subroutine taugb14_d(ncol,iplon,lay,nlayers)
implicit none                                
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)
!----------------------------------------------------------------------------

! ----- Input -----
      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 


    

! ------- Declarations -------

! Local
      integer ::  ind0, ind1, inds, indf, ig
      real(kind=r8) :: tauself, taufor


      integer, parameter :: ng14 = 2
      integer, parameter :: ngs13 = 134

! Lower atmosphere loop
  if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then 

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(14) + 1
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(14) + 1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         do ig = 1, ng14
            tauself = selffac(iplon,lay) * (selfref14_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref14_d(inds+1,ig) - selfref14_d(inds,ig)))
            taufor =  forfac(iplon,lay) * (forref14_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref14_d(indf+1,ig) - forref14_d(indf,ig)))
            taug(iplon,lay,ngs13+ig) = colco2(iplon,lay) * &
                 (fac00(iplon,lay) * absa14_d(ind0,ig) + &
                 fac10(iplon,lay) * absa14_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absa14_d(ind1,ig) + &
                 fac11(iplon,lay) * absa14_d(ind1+1,ig)) &
                 + tauself + taufor
            fracs(iplon,lay,ngs13+ig) = fracrefa14_d(ig)
         enddo





endif
! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then


         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(14) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(14) + 1
         do ig = 1, ng14
            taug(iplon,lay,ngs13+ig) = colco2(iplon,lay) * &
                 (fac00(iplon,lay) * absb14_d(ind0,ig) + &
                 fac10(iplon,lay) * absb14_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absb14_d(ind1,ig) + &
                 fac11(iplon,lay) * absb14_d(ind1+1,ig))
            fracs(iplon,lay,ngs13+ig) = fracrefb14_d(ig)
         enddo


endif
end subroutine taugb14_d

attributes(device)subroutine taugb15_d(ncol,iplon,lay,nlayers,oneminus)
implicit none                                 
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
!----------------------------------------------------------------------------
! ------- Declarations -------

! ----- Input -----
      
     real(kind=r8),intent(in) :: oneminus
   

      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 


! Local
      integer ::  ind0, ind1, inds, indf, indm, ig
      integer :: js, js1, jmn2, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_mn2, specparm_mn2, specmult_mn2, fmn2
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: scalen2, tauself, taufor, n2m1, n2m2, taun2
      real(kind=r8) :: refrat_planck_a, refrat_m_a
      real(kind=r8) :: tau_major, tau_major1

      integer, parameter :: ng15 = 2
      integer, parameter :: ngs14 = 136

      

! P = 1053. mb (Level 1)
      refrat_planck_a = chi_mls_d(4,1)/chi_mls_d(2,1)

! P = 1053.
      refrat_m_a = chi_mls_d(4,1)/chi_mls_d(2,1)


! Lower atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then



         speccomb = coln2o(iplon,lay) + rat_n2oco2(iplon,lay)*colco2(iplon,lay)
         specparm = coln2o(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = coln2o(iplon,lay) + rat_n2oco2_1(iplon,lay)*colco2(iplon,lay)
         specparm1 = coln2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_mn2 = coln2o(iplon,lay) + refrat_m_a*colco2(iplon,lay)
         specparm_mn2 = coln2o(iplon,lay)/speccomb_mn2
         if (specparm_mn2 .ge. oneminus) specparm_mn2 = oneminus
         specmult_mn2 = 8._r8*specparm_mn2
         jmn2 = 1 + int(specmult_mn2)
         fmn2 = mod(specmult_mn2,1.0_r8)

         speccomb_planck = coln2o(iplon,lay)+refrat_planck_a*colco2(iplon,lay)
         specparm_planck = coln2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(15) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(15) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)
         indm = indminor(iplon,lay)

         scalen2 = colbrd(iplon,lay)*scaleminor(iplon,lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1._r8 - fs) * fac00(iplon,lay)
            fac010 = (1._r8 - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif
         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1._r8 - fs1) * fac01(iplon,lay)
            fac011 = (1._r8 - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng15
            tauself = selffac(iplon,lay)* (selfref15_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref15_d(inds+1,ig) - selfref15_d(inds,ig)))
            taufor =  forfac(iplon,lay) * (forref15_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref15_d(indf+1,ig) - forref15_d(indf,ig)))
            n2m1 = ka_mn215_d(jmn2,indm,ig) + fmn2 * &
                 (ka_mn215_d(jmn2+1,indm,ig) - ka_mn215_d(jmn2,indm,ig))
            n2m2 = ka_mn215_d(jmn2,indm+1,ig) + fmn2 * &
                 (ka_mn215_d(jmn2+1,indm+1,ig) - ka_mn215_d(jmn2,indm+1,ig))
            taun2 = scalen2 * (n2m1 + minorfrac(iplon,lay) * (n2m2 - n2m1))

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa15_d(ind0,ig) + &
                    fac100 * absa15_d(ind0+1,ig) + &
                    fac200 * absa15_d(ind0+2,ig) + &
                    fac010 * absa15_d(ind0+9,ig) + &
                    fac110 * absa15_d(ind0+10,ig) + &
                    fac210 * absa15_d(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa15_d(ind0-1,ig) + &
                    fac100 * absa15_d(ind0,ig) + &
                    fac000 * absa15_d(ind0+1,ig) + &
                    fac210 * absa15_d(ind0+8,ig) + &
                    fac110 * absa15_d(ind0+9,ig) + &
                    fac010 * absa15_d(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa15_d(ind0,ig) + &
                    fac100 * absa15_d(ind0+1,ig) + &
                    fac010 * absa15_d(ind0+9,ig) + &
                    fac110 * absa15_d(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa15_d(ind1,ig) + &
                    fac101 * absa15_d(ind1+1,ig) + &
                    fac201 * absa15_d(ind1+2,ig) + &
                    fac011 * absa15_d(ind1+9,ig) + &
                    fac111 * absa15_d(ind1+10,ig) + &
                    fac211 * absa15_d(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa15_d(ind1-1,ig) + &
                    fac101 * absa15_d(ind1,ig) + &
                    fac001 * absa15_d(ind1+1,ig) + &
                    fac211 * absa15_d(ind1+8,ig) + &
                    fac111 * absa15_d(ind1+9,ig) + &
                    fac011 * absa15_d(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa15_d(ind1,ig) + &
                    fac101 * absa15_d(ind1+1,ig) + &
                    fac011 * absa15_d(ind1+9,ig) + &
                    fac111 * absa15_d(ind1+10,ig))
            endif

            taug(iplon,lay,ngs14+ig) = tau_major + tau_major1 &
                 + tauself + taufor &
                 + taun2
            fracs(iplon,lay,ngs14+ig) = fracrefa15_d(ig,jpl) + fpl * &
                 (fracrefa15_d(ig,jpl+1)-fracrefa15_d(ig,jpl))
         enddo




endif

! Upper atmosphere loop
  
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then


         do ig = 1, ng15
            taug(iplon,lay,ngs14+ig) = 0.0_r8
            fracs(iplon,lay,ngs14+ig) = 0.0_r8
         enddo




endif
end subroutine taugb15_d

attributes(device)subroutine taugb16_d(ncol,iplon,lay,nlayers,oneminus)
 implicit none                                    

!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!----------------------------------------------------------------------------

! ----- Input -----
       real(kind=r8),intent(in) :: oneminus
      

      integer,intent(in) :: iplon
      integer,intent(in) :: lay
      integer,intent(in) :: nlayers         ! total number of layers
      integer,intent(in) :: ncol 

! ------- Declarations -------

! Local
      integer ::  ind0, ind1, inds, indf, ig
      integer :: js, js1, jpl
      real(kind=r8) :: speccomb, specparm, specmult, fs
      real(kind=r8) :: speccomb1, specparm1, specmult1, fs1
      real(kind=r8) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(kind=r8) :: p, p4, fk0, fk1, fk2
      real(kind=r8) :: fac000, fac100, fac200, fac010, fac110, fac210
      real(kind=r8) :: fac001, fac101, fac201, fac011, fac111, fac211
      real(kind=r8) :: tauself, taufor
      real(kind=r8) :: refrat_planck_a
      real(kind=r8) :: tau_major, tau_major1


      integer, parameter :: ng16 = 2
      integer, parameter :: ngs15 = 138

! P = 387. mb (Level 6)
      refrat_planck_a = chi_mls_d(1,6)/chi_mls_d(6,6)


! Lower atmosphere loop
  if((iplon>=1 .and. iplon <=ncol).and.(lay>=1 .and. lay<=laytrop(iplon))) then

         speccomb = colh2o(iplon,lay) + rat_h2och4(iplon,lay)*colch4(iplon,lay)
         specparm = colh2o(iplon,lay)/speccomb
         if (specparm .ge. oneminus) specparm = oneminus
         specmult = 8._r8*(specparm)
         js = 1 + int(specmult)
         fs = mod(specmult,1.0_r8)

         speccomb1 = colh2o(iplon,lay) + rat_h2och4_1(iplon,lay)*colch4(iplon,lay)
         specparm1 = colh2o(iplon,lay)/speccomb1
         if (specparm1 .ge. oneminus) specparm1 = oneminus
         specmult1 = 8._r8*(specparm1)
         js1 = 1 + int(specmult1)
         fs1 = mod(specmult1,1.0_r8)

         speccomb_planck = colh2o(iplon,lay)+refrat_planck_a*colch4(iplon,lay)
         specparm_planck = colh2o(iplon,lay)/speccomb_planck
         if (specparm_planck .ge. oneminus) specparm_planck=oneminus
         specmult_planck = 8._r8*specparm_planck
         jpl= 1 + int(specmult_planck)
         fpl = mod(specmult_planck,1.0_r8)

         ind0 = ((jp(iplon,lay)-1)*5+(jt(iplon,lay)-1))*nspa_d(16) + js
         ind1 = (jp(iplon,lay)*5+(jt1(iplon,lay)-1))*nspa_d(16) + js1
         inds = indself(iplon,lay)
         indf = indfor(iplon,lay)

         if (specparm .lt. 0.125_r8) then
            p = fs - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else if (specparm .gt. 0.875_r8) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac000 = fk0*fac00(iplon,lay)
            fac100 = fk1*fac00(iplon,lay)
            fac200 = fk2*fac00(iplon,lay)
            fac010 = fk0*fac10(iplon,lay)
            fac110 = fk1*fac10(iplon,lay)
            fac210 = fk2*fac10(iplon,lay)
         else
            fac000 = (1._r8 - fs) * fac00(iplon,lay)
            fac010 = (1._r8 - fs) * fac10(iplon,lay)
            fac100 = fs * fac00(iplon,lay)
            fac110 = fs * fac10(iplon,lay)
         endif

         if (specparm1 .lt. 0.125_r8) then
            p = fs1 - 1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else if (specparm1 .gt. 0.875_r8) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1 - p - 2.0_r8*p4
            fk2 = p + p4
            fac001 = fk0*fac01(iplon,lay)
            fac101 = fk1*fac01(iplon,lay)
            fac201 = fk2*fac01(iplon,lay)
            fac011 = fk0*fac11(iplon,lay)
            fac111 = fk1*fac11(iplon,lay)
            fac211 = fk2*fac11(iplon,lay)
         else
            fac001 = (1._r8 - fs1) * fac01(iplon,lay)
            fac011 = (1._r8 - fs1) * fac11(iplon,lay)
            fac101 = fs1 * fac01(iplon,lay)
            fac111 = fs1 * fac11(iplon,lay)
         endif

         do ig = 1, ng16
            tauself = selffac(iplon,lay)* (selfref16_d(inds,ig) + selffrac(iplon,lay) * &
                 (selfref16_d(inds+1,ig) - selfref16_d(inds,ig)))
            taufor =  forfac(iplon,lay) * (forref16_d(indf,ig) + forfrac(iplon,lay) * &
                 (forref16_d(indf+1,ig) - forref16_d(indf,ig)))

            if (specparm .lt. 0.125_r8) then
               tau_major = speccomb * &
                    (fac000 * absa16_d(ind0,ig) + &
                    fac100 * absa16_d(ind0+1,ig) + &
                    fac200 * absa16_d(ind0+2,ig) + &
                    fac010 * absa16_d(ind0+9,ig) + &
                    fac110 * absa16_d(ind0+10,ig) + &
                    fac210 * absa16_d(ind0+11,ig))
            else if (specparm .gt. 0.875_r8) then
               tau_major = speccomb * &
                    (fac200 * absa16_d(ind0-1,ig) + &
                    fac100 * absa16_d(ind0,ig) + &
                    fac000 * absa16_d(ind0+1,ig) + &
                    fac210 * absa16_d(ind0+8,ig) + &
                    fac110 * absa16_d(ind0+9,ig) + &
                    fac010 * absa16_d(ind0+10,ig))
            else
               tau_major = speccomb * &
                    (fac000 * absa16_d(ind0,ig) + &
                    fac100 * absa16_d(ind0+1,ig) + &
                    fac010 * absa16_d(ind0+9,ig) + &
                    fac110 * absa16_d(ind0+10,ig))
            endif

            if (specparm1 .lt. 0.125_r8) then
               tau_major1 = speccomb1 * &
                    (fac001 * absa16_d(ind1,ig) + &
                    fac101 * absa16_d(ind1+1,ig) + &
                    fac201 * absa16_d(ind1+2,ig) + &
                    fac011 * absa16_d(ind1+9,ig) + &
                    fac111 * absa16_d(ind1+10,ig) + &
                    fac211 * absa16_d(ind1+11,ig))
            else if (specparm1 .gt. 0.875_r8) then
               tau_major1 = speccomb1 * &
                    (fac201 * absa16_d(ind1-1,ig) + &
                    fac101 * absa16_d(ind1,ig) + &
                    fac001 * absa16_d(ind1+1,ig) + &
                    fac211 * absa16_d(ind1+8,ig) + &
                    fac111 * absa16_d(ind1+9,ig) + &
                    fac011 * absa16_d(ind1+10,ig))
            else
               tau_major1 = speccomb1 * &
                    (fac001 * absa16_d(ind1,ig) + &
                    fac101 * absa16_d(ind1+1,ig) + &
                    fac011 * absa16_d(ind1+9,ig) + &
                    fac111 * absa16_d(ind1+10,ig))
            endif

            taug(iplon,lay,ngs15+ig) = tau_major + tau_major1 &
                 + tauself + taufor
            fracs(iplon,lay,ngs15+ig) = fracrefa16_d(ig,jpl) + fpl * &
                 (fracrefa16_d(ig,jpl+1)-fracrefa16_d(ig,jpl))
         enddo


endif
! Upper atmosphere loop
if((iplon>=1 .and. iplon <=ncol).and.(lay>=laytrop(iplon)+1 .and. lay<=nlayers)) then


         ind0 = ((jp(iplon,lay)-13)*5+(jt(iplon,lay)-1))*nspb_d(16) + 1
         ind1 = ((jp(iplon,lay)-12)*5+(jt1(iplon,lay)-1))*nspb_d(16) + 1
         do ig = 1, ng16
            taug(iplon,lay,ngs15+ig) = colch4(iplon,lay) * &
                 (fac00(iplon,lay) * absb16_d(ind0,ig) + &
                 fac10(iplon,lay) * absb16_d(ind0+1,ig) + &
                 fac01(iplon,lay) * absb16_d(ind1,ig) + &
                 fac11(iplon,lay) * absb16_d(ind1+1,ig))
            fracs(iplon,lay,ngs15+ig) = fracrefb16_d(ig)
         enddo
 



 endif
end subroutine taugb16_d
attributes(global)subroutine rtrnmc_d(ncol,nlayers,istart,iend,iout,nbndlw,ngptlw,tblint,bpade,fluxfac,heatfac)
implicit none                          
! ------- Declarations -------

! ----- Input -----
      integer,value :: ncol
      integer,value :: nlayers                    ! total number of layers
      integer,value :: istart                     ! beginning band of calculation
      integer,value :: iend                       ! ending band of calculation
      integer,value :: iout 
      integer,value :: nbndlw 
      integer,value :: ngptlw
      !integer,value :: ntbl 
      real(kind=r8),value :: tblint                      ! output option flag
      real(kind=r8),value :: bpade
      real(kind=r8),value :: fluxfac
      real(kind=r8),value :: heatfac

     

! Clouds
     
! ----- Local -----
! Declarations for radiative transfer

  
      real(kind=r8) :: atot(52)
      real(kind=r8) :: atrans(52)
      real(kind=r8) :: bbugas(52)
      real(kind=r8) :: bbutot(52)
      real(kind=r8) :: clrurad(0:52)
      real(kind=r8) :: clrdrad(0:52)

      real(kind=r8) :: uflux(0:52)
      real(kind=r8) :: dflux(0:52)
      real(kind=r8) :: urad(0:52)
      real(kind=r8) :: drad(0:52)
      real(kind=r8) :: uclfl(0:52)
      real(kind=r8) :: dclfl(0:52)
      
      integer :: icldlyr(52)

      real(kind=r8) :: secdiff(16)                 

    
       real(kind=r8) :: a0(16),a1(16),a2(16)  ! diffusivity angle adjustment coefficients
      real(kind=r8) :: wtdiff, rec_6
      real(kind=r8) :: transcld, radld, radclrd, plfrac, blay, dplankup, dplankdn
      real(kind=r8) :: odepth, odtot, odepth_rec, odtot_rec, gassrc
      real(kind=r8) :: tblind, tfactot, bbd, bbdtot, tfacgas, transc, tausfac
      real(kind=r8) :: rad0, reflect, radlu, radclru

      integer :: iplon,k                     
      integer :: ibnd, ib, iband, lay, lev, l, ig        ! loop indices
      integer :: igc                                     ! g-point interval counter
      integer :: iclddn                                  ! flag for cloud in down path
      integer :: ittot, itgas, itr

iplon=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if (iplon>=1 .and. iplon<=ncol) then
       wtdiff =0.5_r8
       rec_6 =0.166667_r8


a0(1)=1.66_r8
a0(2)=1.55_r8
a0(3)=1.58_r8
a0(4)=1.66_r8
a0(5)=1.54_r8
a0(6)=1.454_r8
a0(7)=1.89_r8
a0(8)=1.33_r8
a0(9)=1.668_r8
a0(10)=1.66_r8
a0(11)=1.66_r8
a0(12)=1.66_r8
a0(13)=1.66_r8
a0(14)=1.66_r8
a0(15)=1.66_r8
a0(16)=1.66_r8

a1(1)=0.00_r8
a1(2)=0.25_r8
a1(3)=0.22_r8
a1(4)=0.00_r8
a1(5)=0.13_r8
a1(6)=0.446_r8
a1(7)=-0.10_r8
a1(8)=0.40_r8
a1(9)=-0.006_r8
a1(10)=0.00_r8
a1(11)=0.00_r8
a1(12)=0.00_r8
a1(13)=0.00_r8
a1(14)=0.00_r8
a1(15)=0.00_r8
a1(16)=0.00_r8

a2(1)=0.00_r8
a2(2)=-12.0_r8
a2(3)=-11.7_r8
a2(4)=0.00_r8
a2(5)=-0.72_r8
a2(6)=-0.243_r8
a2(7)=0.19_r8
a2(8)=-0.062_r8
a2(9)=0.414_r8
a2(10)=0.00_r8
a2(11)=0.00_r8
a2(12)=0.00_r8
a2(13)=0.00_r8
a2(14)=0.00_r8
a2(15)=0.00_r8
a2(16)=0.00_r8




! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
! and 1.80) as a function of total column water vapor.  The function
! has been defined to minimize flux and cooling rate errors in these bands
! over a wide range of precipitable water values.
   

      do ibnd = 1,nbndlw
         if (ibnd.eq.1 .or. ibnd.eq.4 .or. ibnd.ge.10) then
           secdiff(ibnd) = 1.66_r8
         else
           secdiff(ibnd) = a0(ibnd) + a1(ibnd)*exp(a2(ibnd)*pwvcm(iplon))
         endif
      enddo
      if (pwvcm(iplon).lt.1.0) secdiff(6) = 1.80_r8
      if (pwvcm(iplon).gt.7.1) secdiff(7) = 1.50_r8

      urad(0) = 0.0_r8
      drad(0) = 0.0_r8
      totuflux(iplon,0) = 0.0_r8
      totdflux(iplon,0) = 0.0_r8
      clrurad(0) = 0.0_r8
      clrdrad(0) = 0.0_r8
      totuclfl(iplon,0) = 0.0_r8
      totdclfl(iplon,0) = 0.0_r8

      do lay = 1, nlayers
         urad(lay) = 0.0_r8
         drad(lay) = 0.0_r8
         totuflux(iplon,lay) = 0.0_r8
         totdflux(iplon,lay) = 0.0_r8
         clrurad(lay) = 0.0_r8
         clrdrad(lay) = 0.0_r8
         totuclfl(iplon,lay) = 0.0_r8
         totdclfl(iplon,lay) = 0.0_r8
         icldlyr(lay) = 0

! Change to band loop?
         do ig = 1, ngptlw
            if (cldfmc(ig,iplon,lay) .eq. 1._r8) then
               ib = ngb_d(ig)
               odcld(iplon,lay,ig) = secdiff(ib) * taucmc(ig,iplon,lay)
               transcld = exp(-odcld(iplon,lay,ig))
               abscld(iplon,lay,ig) = 1._r8 - transcld
               efclfrac(iplon,lay,ig) = abscld(iplon,lay,ig) * cldfmc(ig,iplon,lay)
               icldlyr(lay) = 1
            else
               odcld(iplon,lay,ig) = 0.0_r8
               abscld(iplon,lay,ig) = 0.0_r8
               efclfrac(iplon,lay,ig) = 0.0_r8
            endif
         enddo

      enddo

      igc = 1
! Loop over frequency bands.
      do iband = istart,iend

! Reinitialize g-point counter for each band if output for each band is requested.
         if (iout.gt.0.and.iband.ge.2) igc = ngs_d(iband-1)+1

! Loop over g-channels.
 1000    continue

! Radiative transfer starts here.
         radld = 0._r8
         radclrd = 0._r8
         iclddn = 0

! Downward radiative transfer loop.

         do lev = nlayers, 1, -1
               plfrac = fracs(iplon,lev,igc)
               blay = planklay(iplon,lev,iband)
               dplankup = planklev(iplon,lev,iband) - blay
               dplankdn = planklev(iplon,lev-1,iband) - blay
               odepth = secdiff(iband) * taut(iplon,lev,igc)
               if (odepth .lt. 0.0_r8) odepth = 0.0_r8
!  Cloudy layer
               if (icldlyr(lev).eq.1) then
                  iclddn = 1
                  odtot = odepth + odcld(iplon,lev,igc)
                  if (odtot .lt. 0.06_r8) then
                     atrans(lev) = odepth - 0.5_r8*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     atot(lev) =  odtot - 0.5_r8*odtot*odtot
                     odtot_rec = rec_6*odtot
                     bbdtot =  plfrac * (blay+dplankdn*odtot_rec)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                     radld = radld - radld * (atrans(lev) + &
                         efclfrac(iplon,lev,igc) * (1. - atrans(lev))) + &
                         gassrc + cldfmc(igc,iplon,lev) * &
                         (bbdtot * atot(lev) - gassrc)
                     drad(lev-1) = drad(lev-1) + radld

                     bbugas(lev) =  plfrac * (blay+dplankup*odepth_rec)
                     bbutot(lev) =  plfrac * (blay+dplankup*odtot_rec)

                  elseif (odepth .le. 0.06_r8) then
                     atrans(lev) = odepth - 0.5_r8*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     odtot = odepth + odcld(iplon,lev,igc)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_r8
                     tfactot = tfn_tbl_d(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                     atot(lev) = 1. - exp_tbl_d(ittot)

                     radld = radld - radld * (atrans(lev) + &
                         efclfrac(iplon,lev,igc) * (1._r8 - atrans(lev))) + &
                         gassrc + cldfmc(igc,iplon,lev) * &
                         (bbdtot * atot(lev) - gassrc)
                     drad(lev-1) = drad(lev-1) + radld

                     bbugas(lev) = plfrac * (blay + dplankup*odepth_rec)
                     bbutot(lev) = plfrac * (blay + tfactot * dplankup)

                  else

                     tblind = odepth/(bpade+odepth)
                     itgas = tblint*tblind+0.5_r8
                     odepth = tau_tbl_d(itgas)
                     atrans(lev) = 1._r8 - exp_tbl_d(itgas)
                     tfacgas = tfn_tbl_d(itgas)
                     gassrc = atrans(lev) * plfrac * (blay + tfacgas*dplankdn)

                     odtot = odepth + odcld(iplon,lev,igc)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_r8
                     tfactot = tfn_tbl_d(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+tfacgas*dplankdn)
                     atot(lev) = 1._r8 - exp_tbl_d(ittot)

                  radld = radld - radld * (atrans(lev) + &
                    efclfrac(iplon,lev,igc) * (1._r8 - atrans(lev))) + &
                    gassrc + cldfmc(igc,iplon,lev) * &
                    (bbdtot * atot(lev) - gassrc)
                  drad(lev-1) = drad(lev-1) + radld
                  bbugas(lev) = plfrac * (blay + tfacgas * dplankup)
                  bbutot(lev) = plfrac * (blay + tfactot * dplankup)
                  endif
!  Clear layer
               else
                  if (odepth .le. 0.06_r8) then
                     atrans(lev) = odepth-0.5_r8*odepth*odepth
                     odepth = rec_6*odepth
                     bbd = plfrac*(blay+dplankdn*odepth)
                     bbugas(lev) = plfrac*(blay+dplankup*odepth)
                  else
                     tblind = odepth/(bpade+odepth)
                     itr = tblint*tblind+0.5_r8
                     transc = exp_tbl_d(itr)
                     atrans(lev) = 1._r8-transc
                     tausfac = tfn_tbl_d(itr)
                     bbd = plfrac*(blay+tausfac*dplankdn)
                     bbugas(lev) = plfrac * (blay + tausfac * dplankup)
                  endif
                  radld = radld + (bbd-radld)*atrans(lev)
                  drad(lev-1) = drad(lev-1) + radld
               endif
!  Set clear sky stream to total sky stream as long as layers
!  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
!  and clear sky stream must be computed separately from that point.
                  if (iclddn.eq.1) then
                     radclrd = radclrd + (bbd-radclrd) * atrans(lev)
                     clrdrad(lev-1) = clrdrad(lev-1) + radclrd
                  else
                     radclrd = radld
                     clrdrad(lev-1) = drad(lev-1)
                  endif
            enddo


         rad0 = fracs(iplon,1,igc) * plankbnd(iplon,iband)
!  Add in specular reflection of surface downward radiance.
         reflect = 1._r8 - semiss(iplon,iband)
         radlu = rad0 + reflect * radld
         radclru = rad0 + reflect * radclrd


! Upward radiative transfer loop.
         urad(0) = urad(0) + radlu
         clrurad(0) = clrurad(0) + radclru

         do lev = 1, nlayers
!  Cloudy layer
            if (icldlyr(lev) .eq. 1) then
               gassrc = bbugas(lev) * atrans(lev)
               radlu = radlu - radlu * (atrans(lev) + &
                   efclfrac(iplon,lev,igc) * (1._r8 - atrans(lev))) + &
                   gassrc + cldfmc(igc,iplon,lev) * &
                   (bbutot(lev) * atot(lev) - gassrc)
               urad(lev) = urad(lev) + radlu
!  Clear layer
            else
               radlu = radlu + (bbugas(lev)-radlu)*atrans(lev)
               urad(lev) = urad(lev) + radlu
            endif

               if (iclddn.eq.1) then
                  radclru = radclru + (bbugas(lev)-radclru)*atrans(lev)
                  clrurad(lev) = clrurad(lev) + radclru
               else
                  radclru = radlu
                  clrurad(lev) = urad(lev)
               endif
         enddo

! Increment g-point counter
         igc = igc + 1
! Return to continue radiative transfer for all g-channels in present band
         if (igc .le. ngs_d(iband)) go to 1000

         do lev = nlayers, 0, -1
            uflux(lev) = urad(lev)*wtdiff
            dflux(lev) = drad(lev)*wtdiff
            urad(lev) = 0.0_r8
            drad(lev) = 0.0_r8
            totuflux(iplon,lev) = totuflux(iplon,lev) + uflux(lev) * delwave_d(iband)
            totdflux(iplon,lev) = totdflux(iplon,lev) + dflux(lev) * delwave_d(iband)
            uclfl(lev) = clrurad(lev)*wtdiff
            dclfl(lev) = clrdrad(lev)*wtdiff
            clrurad(lev) = 0.0_r8
            clrdrad(lev) = 0.0_r8
            totuclfl(iplon,lev) = totuclfl(iplon,lev) + uclfl(lev) * delwave_d(iband)
            totdclfl(iplon,lev) = totdclfl(iplon,lev) + dclfl(lev) * delwave_d(iband)
         enddo

! End spectral band loop
      enddo

! Calculate fluxes at surface
      totuflux(iplon,0) = totuflux(iplon,0) * fluxfac
      totdflux(iplon,0) = totdflux(iplon,0) * fluxfac
      fnet(iplon,0) = totuflux(iplon,0) - totdflux(iplon,0)
      totuclfl(iplon,0) = totuclfl(iplon,0) * fluxfac
      totdclfl(iplon,0) = totdclfl(iplon,0) * fluxfac
      fnetc(iplon,0) = totuclfl(iplon,0) - totdclfl(iplon,0)

! Calculate fluxes at model levels
      do lev = 1, nlayers
         totuflux(iplon,lev) = totuflux(iplon,lev) * fluxfac
         totdflux(iplon,lev) = totdflux(iplon,lev) * fluxfac
         fnet(iplon,lev) = totuflux(iplon,lev) - totdflux(iplon,lev)
         totuclfl(iplon,lev) = totuclfl(iplon,lev) * fluxfac
         totdclfl(iplon,lev) = totdclfl(iplon,lev) * fluxfac
         fnetc(iplon,lev) = totuclfl(iplon,lev) - totdclfl(iplon,lev)
         l = lev - 1

! Calculate heating rates at model layers
         htr(iplon,l)=heatfac*(fnet(iplon,l)-fnet(iplon,lev))/(pz(iplon,l)-pz(iplon,lev))
         htrc(iplon,l)=heatfac*(fnetc(iplon,l)-fnetc(iplon,lev))/(pz(iplon,l)-pz(iplon,lev))
      enddo

! Set heating rate to zero in top layer
      htr(iplon,nlayers) = 0.0_r8
      htrc(iplon,nlayers) = 0.0_r8
      
        do k = 0, nlayers
            uflx_d(iplon,k+1) = totuflux(iplon,k)
            dflx_d(iplon,k+1) = totdflux(iplon,k)
            uflxc_d(iplon,k+1) = totuclfl(iplon,k)
            dflxc_d(iplon,k+1) = totdclfl(iplon,k)
         enddo
         do k = 0, nlayers-1
            hr_d(iplon,k+1) = htr(iplon,k)
            hrc_d(iplon,k+1) = htrc(iplon,k)
         enddo

endif
end subroutine rtrnmc_d




end module
