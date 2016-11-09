!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: avhrr_pixel_aerosol.f90 (src)
!       AEROSOL_PROPERTIES (program)
!
! PURPOSE: a module for pixel level aerosol properties from AVHRR
!
! DESCRIPTION: this is the single channel algorithm used in PATMOS for
!              aerosol optical in ch1 over the ocean
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Lookup tables provided by A. Ignatov (02/2005)
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! NOTES:
! Input and output to these routines is through public arrays
! through the PIXEL_COMMON module. 
!
! Public routines used in this module:
!  READ_AER_CH123A_REF_LUTS - read the lookup tables (LUTS)
!  PIXEL_AER_RET_OCEAN - perform the estimation of optical depth
!
! Private routines used in this module:
!  AER_RET - routine used in estimation aerosol optical depth.
!
! File I/O 
! Logical units are opened and closed during LUT read.
!
! public variables used
!  Ref_Ch1 - channel 1 relflectance (%)
!  Ref_Ch2 - channel 2 reflectance (%)
!  Ref_Ch6 - channel 3a reflectance (%)
!  satzen - sensor zenith angle (deg)
!  Solzen - solar zenith angle (deg)
!  Relaz - relative azimuth angle (deg)
!
! output passed through PIXEL_COMMON
! aot1 - channel 1 aerosol optical depth
! aot2 - channel 2 aerosol optical depth
! aot3a - channel 3a aerosol optical depth
!
! aerosol lookup table variables
! tau_Lut - optical depth used in tables - dimension (Ntau)
! Solzen_Lut - solar zenith angles used in tables - dimension(Nsolzen)
! zen_Lut - sensor zenith angles used in tables -dimension (nzen)
! az_Lut - relative azimuth angles used in tables - dimension(Naz)
! Ref_Lut_Ch1_aer - ch1 reflectances (0-1) of lookup table - 
!                    dimension(Nsolzen,nzen,Ntau,Naz)
! Ref_Lut_Ch2_aer - ch2 reflectances (0-1) of lookup table
! Ref_Lut_Ch3a_aer - ch3a reflectances (0-1) of lookup table
!
! Note, aerosol product quality flags are computed elsewhere.
!
!--------------------------------------------------------------------------------------
module AVHRR_PIXEL_AEROSOL
   
   use CONSTANTS,only: &
      sym &
      , real4 &
      , Missing_Value_Real4 &
      , DTOR
   
   use FILE_TOOLS,only: &
    get_lun
   
   
   use PIXEL_COMMON,only: &
      image &
      , geo &
      , bad_pixel_mask &
      , sfc &
      , ch &
      , aot1 &
      , aot2 &
      , aot3a &
      , sensor &
      , ch3a_on_avhrr

   implicit none
   private
   public:: READ_AER_CH123A_REF_LUTS,&
          PIXEL_AER_RET_OCEAN
  

!------------------------------------------------------------------------
 integer, parameter, private:: nzen = 15, Naz = 19, Nsolzen = 15, Ntau = 7

!--- aerosol lut variables
 real, dimension(Ntau), save, private::tau_Lut
 real, dimension(Nsolzen), save, private::Solzen_Lut
 real, dimension(nzen), save, private::zen_Lut
 real, dimension(Naz), save, private::az_Lut
 real, save, public:: lambda_ref
 real, save, private:: delta_Solzen,delta_zen,delta_az,delta_tau
 real, dimension(Nsolzen,nzen,Ntau,Naz), save, private::  Ref_Lut_Ch1_aer, &
                                                 Ref_Lut_Ch2_aer,Ref_Lut_Ch3a_aer

!--- other module-wide variables
 real, dimension(Ntau), save, private::temp_Ref_Lut

 contains

!----------------------------------------------------------------------------
! grid_cell aerosol retrieval
!
! This uses linear in tau and angle interpolation - should improve this base
! on scheme used for operational aerosol retrieval.
!
! temp_Ref_Lut is a reflectance vector of dimension(Ntau) that is already
!         interpolated in the angular dimensions.
!
!----------------------------------------------------------------------------
subroutine PIXEL_AER_RET_OCEAN(jmin,jmax)

   integer, intent(in):: jmin, jmax

   integer:: i,j,it,iu,iv,ifail
   real(kind=real4):: v,t,u,missing_value
 
   missing_value = -999.0
   ifail = 0

!--------------------------------------------------------------------------
! loop over pixels in scanlines
!--------------------------------------------------------------------------
scan_loop: do j = jmin, jmax + jmin - 1

pixel_loop:  do i = 1, Image%Number_Of_Elements

!--- check for a bad pixel
if (Bad_Pixel_Mask(i,j) == sym%YES) then 
    cycle
endif

!--- check for valid illumination and surface

if ((Geo%Solzen(i,j) > 87.0) .or. (Sfc%Sfc_Type(i,j) /= sym%WATER_SFC) .or. (Sfc%Snow(i,j) /= sym%NO_SNOW)) then

   Aot1(i,j) = Missing_Value_Real4
   Aot2(i,j) = Missing_Value_Real4
   Aot3a(i,j) = Missing_Value_Real4


 else

  it = max(1,min(Nsolzen-1,1+int( (Geo%Solzen(i,j)-Solzen_Lut(1)) / (Delta_Solzen))))
  iu = max(1,min(Nzen-1,1+int( (Geo%Satzen(i,j)-Zen_Lut(1)) / (Delta_Zen))))
  iv = max(1,min(Naz-1,1+int( min(Geo%Relaz(i,j),360.0-Geo%Relaz(i,j)) / (delta_az)) ))
  u = (Geo%Satzen(i,j) - Zen_Lut(iu))/(Zen_Lut(iu+1)-Zen_Lut(iu))
  v = (min(Geo%Relaz(i,j),360.0-Geo%Relaz(i,j)) - az_Lut(iv))/(az_Lut(iv+1)-az_Lut(iv))
  t =  (Geo%Solzen(i,j) - Solzen_Lut(it))/(Solzen_Lut(it+1)-Solzen_Lut(it))

  ifail = 0
 
 
  !--  channel 1
  if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
    temp_Ref_Lut = (1.0-t)*(1.0-u)*(1.0-v)*Ref_Lut_Ch1_aer(it,iu,:,iv) +  &
                 (1.0-t)*(1.0-u)*(v)*Ref_Lut_Ch1_aer(it,iu,:,iv+1) +  &
                 (t)*(1.0-u)*(1.0-v)*Ref_Lut_Ch1_aer(it+1,iu,:,iv) +  &
                 (t)*(1.0-u)*(v)*Ref_Lut_Ch1_aer(it+1,iu,:,iv+1) +  &
                 (1.0-t)*(u)*(1.0-v)*Ref_Lut_Ch1_aer(it,iu+1,:,iv) +  &
                 (1.0-t)*(u)*(v)*Ref_Lut_Ch1_aer(it,iu+1,:,iv+1) +  &
                 (t)*(u)*(1.0-v)*Ref_Lut_Ch1_aer(it+1,iu+1,:,iv) +  &
                 (t)*(u)*(v)*Ref_Lut_Ch1_aer(it+1,iu+1,:,iv+1)
    call AER_RET(ch(1)%Ref_Toa(i,j)/100.0, aot1(i,j))
  endif

  !--  channel 2
  if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then
     temp_Ref_Lut = (1.0-t)*(1.0-u)*(1.0-v)*Ref_Lut_Ch2_aer(it,iu,:,iv) +  &
                 (1.0-t)*(1.0-u)*(v)*Ref_Lut_Ch2_aer(it,iu,:,iv+1) +  &
                 (t)*(1.0-u)*(1.0-v)*Ref_Lut_Ch2_aer(it+1,iu,:,iv) +  &
                 (t)*(1.0-u)*(v)*Ref_Lut_Ch2_aer(it+1,iu,:,iv+1) +  &
                 (1.0-t)*(u)*(1.0-v)*Ref_Lut_Ch2_aer(it,iu+1,:,iv) +  &
                 (1.0-t)*(u)*(v)*Ref_Lut_Ch2_aer(it,iu+1,:,iv+1) +  &
                 (t)*(u)*(1.0-v)*Ref_Lut_Ch2_aer(it+1,iu+1,:,iv) +  &
                 (t)*(u)*(v)*Ref_Lut_Ch2_aer(it+1,iu+1,:,iv+1)
     call AER_RET(ch(2)%Ref_Toa(i,j)/100.0, aot2(i,j))
   endif


   if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
     if (Ch3a_On_AVHRR(j) == sym%YES) then
     !--  channel 3a
       temp_Ref_Lut = (1.0-t)*(1.0-u)*(1.0-v)*Ref_Lut_Ch3a_aer(it,iu,:,iv) +  &
                 (1.0-t)*(1.0-u)*(v)*Ref_Lut_Ch3a_aer(it,iu,:,iv+1) +  &
                 (t)*(1.0-u)*(1.0-v)*Ref_Lut_Ch3a_aer(it+1,iu,:,iv) +  &
                 (t)*(1.0-u)*(v)*Ref_Lut_Ch3a_aer(it+1,iu,:,iv+1) +  &
                 (1.0-t)*(u)*(1.0-v)*Ref_Lut_Ch3a_aer(it,iu+1,:,iv) +  &
                 (1.0-t)*(u)*(v)*Ref_Lut_Ch3a_aer(it,iu+1,:,iv+1) +  &
                 (t)*(u)*(1.0-v)*Ref_Lut_Ch3a_aer(it+1,iu+1,:,iv) +  &
                 (t)*(u)*(v)*Ref_Lut_Ch3a_aer(it+1,iu+1,:,iv+1)
      call AER_RET(ch(6)%Ref_Toa(i,j)/100.0, aot3a(i,j))
    endif
   endif

 endif            !end of check for appropriate retrieval conditions

 end do pixel_loop
end do scan_loop
 
end subroutine PIXEL_AER_RET_OCEAN

!------------------------------------------------------------
! Actual routine to estimate aerosol optical depth from LUT
!------------------------------------------------------------
subroutine AER_RET(y,aot)

real, intent(in):: y
real, intent(out):: aot
integer:: is
real:: s

!--- linear interp
  call LOCATE(temp_Ref_Lut,Ntau,y,is)
  is = max(1,min(Ntau-1,is))
  s = (y - temp_Ref_Lut(is))/(temp_Ref_Lut(is+1) - temp_Ref_Lut(is))
  aot = (1.0-s)*tau_Lut(is) + s*tau_Lut(is+1) 
end subroutine AER_RET

!------------------------------------------------------------------------------
! READ IN THE NESDIS OPERATIONAL LOOKUP TABLES FOR CH1, CH2 and CH3A OCEANIC AEROSOL
!------------------------------------------------------------------------------
 subroutine READ_AER_CH123A_REF_LUTS(ancil_data_dir,wmonum)
  character (len=*), intent(in):: ancil_data_dir
  integer, intent(in):: wmonum
  character (len=1020):: ch1_aer_Lut_file,ch2_aer_Lut_file,ch3a_aer_Lut_file
  character (len=1020):: header,basename
  integer:: iSolzen,izen,itau
  integer:: nzen1,Naz1,Nsolzen1,Ntau1
  integer:: nzen2,Naz2,Nsolzen2,Ntau2
  integer:: nzen3,Naz3,Nsolzen3,Ntau3
  real(kind=real4):: cos_Solzen
  integer:: ch1_lun
  integer:: ch2_lun
  integer:: ch3a_lun


  if (wmonum == 708) then
      basename = "Ch1F.lut.dat-NOAATN"
  elseif (wmonum == 706) then
      basename = "Ch1F.lut.dat-NOAA06"
  elseif(wmonum == 707) then
      basename = "Ch1F.lut.dat-NOAA07"
  elseif(wmonum == 200) then
      basename = "Ch1F.lut.dat-NOAA08"
  elseif(wmonum == 201) then
      basename = "Ch1F.lut.dat-NOAA09"
  elseif(wmonum == 202) then
      basename = "Ch1F.lut.dat-NOAA10"
  elseif(wmonum == 203) then
      basename = "Ch1F.lut.dat-NOAA11"
  elseif(wmonum == 204) then
      basename = "Ch1F.lut.dat-NOAA12"
  elseif(wmonum == 205) then
      basename = "Ch1F.lut.dat-NOAA14"
  elseif(wmonum == 206) then
      basename = "Ch1F.lut.dat-NOAA15"
  elseif(wmonum == 207) then
      basename = "Ch1F.lut.dat-NOAA16"
  elseif(wmonum == 208) then
      basename = "Ch1F.lut.dat-NOAA17"
  elseif(wmonum == 209) then
      basename = "Ch1F.lut.dat-NOAA18"
  elseif(wmonum == 223) then
      basename = "Ch1F.lut.dat-NOAA19"
  elseif(wmonum == 4) then
      basename = "Ch1F.lut.dat-MetOpA"
  elseif(wmonum == 3) then
      basename = "Ch1F.lut.dat-MetOpB"
  elseif(wmonum == 5) then
      basename = "Ch1F.lut.dat-MetOpC"
  else
      print *, "no aerosol lookup tables for satellite number = ", wmonum
      print *, "stopping"
      stop
  endif

  ch1_aer_Lut_file = trim(ancil_data_dir)//"static/luts/aerosol/"//trim(basename)//"_1"
  ch2_aer_Lut_file = trim(ancil_data_dir)//"static/luts/aerosol/"//trim(basename)//"_2"
  ch3a_aer_Lut_file = trim(ancil_data_dir)//"static/luts/aerosol/"//trim(basename)//"_3"

  if (wmonum >= 200 .and. wmonum <= 205) then
    ch3a_aer_Lut_file = trim(ancil_data_dir)//"static/luts/aerosol/"//"Ch1F.lut.dat-NOAA15_3"
  endif
  if (wmonum >= 706 .and. wmonum <= 708) then
    ch3a_aer_Lut_file = trim(ancil_data_dir)//"static/luts/aerosol/"//"Ch1F.lut.dat-NOAA15_3"
  endif

  ch1_lun = GET_LUN()
  open(unit=ch1_lun,file=ch1_aer_Lut_file,status="old",position="rewind",action="read")
  ch2_lun = GET_LUN()
  open(unit=ch2_lun,file=ch2_aer_Lut_file,status="old",position="rewind",action="read")
  ch3a_lun = GET_LUN()
  open(unit=ch3a_lun,file=ch3a_aer_Lut_file,status="old",position="rewind",action="read")

   read(unit=ch1_lun,fmt="(4i5)") nzen1, Naz1, Nsolzen1,Ntau1
   read(unit=ch2_lun,fmt="(4i5)") nzen2, Naz2, Nsolzen2,Ntau2
   read(unit=ch3a_lun,fmt="(4i5)") nzen3, Naz3, Nsolzen3,Ntau3

   Ref_Lut_Ch1_aer = 0.0
   Ref_Lut_Ch2_aer = 0.0
   Ref_Lut_Ch3a_aer = 0.0
   tau_Lut = 0.00
   Solzen_Lut = 0.00
   zen_Lut = 0.00
   az_Lut = 0.00
   lambda_ref = 0.63
   delta_tau = 0.15

   temp_Ref_Lut = 0.0

   delta_Solzen = 90.0 / Nsolzen
   do iSolzen = 1,Nsolzen
    Solzen_Lut(iSolzen) = (iSolzen-1)*delta_Solzen
   enddo


   read(unit=ch1_lun, fmt=*) tau_Lut
   read(unit=ch2_lun, fmt=*) tau_Lut
   read(unit=ch3a_lun, fmt=*) tau_Lut

   tau_Lut = tau_Lut * delta_tau

   do iSolzen = 1,Nsolzen
       cos_Solzen = cos(Solzen_Lut(iSolzen)*dtor)
     do itau = 1, Ntau
       read(unit=ch1_lun,fmt="(a)") header
       read(unit=ch2_lun,fmt="(a)") header
       read(unit=ch3a_lun,fmt="(a)") header

       read(unit=ch1_lun,fmt="(10f6.0,/,9f6.0)") az_Lut
       read(unit=ch2_lun,fmt="(10f6.0,/,9f6.0)") az_Lut
       read(unit=ch3a_lun,fmt="(10f6.0,/,9f6.0)") az_Lut
      
       do izen = 1, nzen
        read(unit=ch1_lun,fmt=*) zen_Lut(izen), Ref_Lut_Ch1_aer(iSolzen,izen,itau,1:5)
        read(unit=ch1_lun,fmt=*) Ref_Lut_Ch1_aer(iSolzen,izen,itau,6:10)
        read(unit=ch1_lun,fmt=*) Ref_Lut_Ch1_aer(iSolzen,izen,itau,11:15)
        read(unit=ch1_lun,fmt=*) Ref_Lut_Ch1_aer(iSolzen,izen,itau,16:19)

        read(unit=ch2_lun,fmt=*) zen_Lut(izen), Ref_Lut_Ch2_aer(iSolzen,izen,itau,1:5)
        read(unit=ch2_lun,fmt=*) Ref_Lut_Ch2_aer(iSolzen,izen,itau,6:10)
        read(unit=ch2_lun,fmt=*) Ref_Lut_Ch2_aer(iSolzen,izen,itau,11:15)
        read(unit=ch2_lun,fmt=*) Ref_Lut_Ch2_aer(iSolzen,izen,itau,16:19)

        read(unit=ch3a_lun,fmt=*) zen_Lut(izen), Ref_Lut_Ch3a_aer(iSolzen,izen,itau,1:5)
        read(unit=ch3a_lun,fmt=*) Ref_Lut_Ch3a_aer(iSolzen,izen,itau,6:10)
        read(unit=ch3a_lun,fmt=*) Ref_Lut_Ch3a_aer(iSolzen,izen,itau,11:15)
        read(unit=ch3a_lun,fmt=*) Ref_Lut_Ch3a_aer(iSolzen,izen,itau,16:19)

       enddo 
     enddo
     Ref_Lut_Ch1_aer(iSolzen,:,:,:) = Ref_Lut_Ch1_aer(iSolzen,:,:,:)/cos_Solzen
     Ref_Lut_Ch2_aer(iSolzen,:,:,:) = Ref_Lut_Ch2_aer(iSolzen,:,:,:)/cos_Solzen
     Ref_Lut_Ch3a_aer(iSolzen,:,:,:) = Ref_Lut_Ch3a_aer(iSolzen,:,:,:)/cos_Solzen
   enddo

  close(unit=ch1_lun)
  close(unit=ch2_lun)
  close(unit=ch3a_lun)

   Delta_Zen = Zen_Lut(2) - Zen_Lut(1)
   Delta_Az = Az_Lut(2) - Az_Lut(1)

 end subroutine READ_AER_CH123A_REF_LUTS

end module AVHRR_PIXEL_AEROSOL
