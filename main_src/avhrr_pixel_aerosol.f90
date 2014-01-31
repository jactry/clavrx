!$Id: avhrr_pixel_aerosol.f90,v 1.11.2.2 2014/01/26 04:48:31 heidinger Exp $
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
! tau_lut - optical depth used in tables - dimension (ntau)
! Solzen_lut - solar zenith angles used in tables - dimension(nSolzen)
! zen_lut - sensor zenith angles used in tables -dimension (nzen)
! az_lut - relative azimuth angles used in tables - dimension(naz)
! Ref_lut_Ch1_aer - ch1 reflectances (0-1) of lookup table - 
!                    dimension(nSolzen,nzen,ntau,naz)
! Ref_lut_Ch2_aer - ch2 reflectances (0-1) of lookup table
! Ref_lut_Ch3a_aer - ch3a reflectances (0-1) of lookup table
!
! Note, aerosol product quality flags are computed elsewhere.
!
!--------------------------------------------------------------------------------------
module AEROSOL_PROPERTIES
 use PIXEL_COMMON
 use CONSTANTS
 use FILE_UTILITY
 use NUMERICAL_ROUTINES
 implicit none
 private
 public:: READ_AER_CH123A_REF_LUTS,&
          PIXEL_AER_RET_OCEAN
 private:: AER_RET

!------------------------------------------------------------------------
 integer, parameter, private:: nzen = 15, naz = 19, nSolzen = 15, ntau = 7

!--- aerosol lut variables
 real, dimension(ntau), save, private::tau_lut
 real, dimension(nSolzen), save, private::Solzen_lut
 real, dimension(nzen), save, private::zen_lut
 real, dimension(naz), save, private::az_lut
 real, save, public:: lambda_ref
 real, save, private:: delta_Solzen,delta_zen,delta_az,delta_tau
 real, dimension(nSolzen,nzen,ntau,naz), save, private::  Ref_lut_Ch1_aer, &
                                                 Ref_lut_Ch2_aer,Ref_lut_Ch3a_aer

!--- other module-wide variables
 real, dimension(ntau), save, private::temp_Ref_lut

 contains

!----------------------------------------------------------------------------
! grid_cell aerosol retrieval
!
! This uses linear in tau and angle interpolation - should improve this base
! on scheme used for operational aerosol retrieval.
!
! temp_Ref_lut is a reflectance vector of dimension(ntau) that is already
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

pixel_loop:  do i = 1, num_pix

!--- check for a bad pixel
if (bad_pixel_mask(i,j) == sym%YES) then 
    cycle
endif

!--- check for valid illumination and surface

if ((Solzen(i,j) > 87.0) .or. (sfc_type(i,j) /= sym%WATER_SFC) .or. (snow(i,j) /= sym%NO_SNOW)) then

   aot1(i,j) = missing_value_real4
   aot2(i,j) = missing_value_real4
   aot3a(i,j) = missing_value_real4


 else

  it = max(1,min(nSolzen-1,1+int( (Solzen(i,j)-Solzen_lut(1)) / (delta_Solzen))))
  iu = max(1,min(nzen-1,1+int( (satzen(i,j)-zen_lut(1)) / (delta_zen))))
  iv = max(1,min(naz-1,1+int( min(Relaz(i,j),360.0-Relaz(i,j)) / (delta_az)) ))
  u = (satzen(i,j) - zen_lut(iu))/(zen_lut(iu+1)-zen_lut(iu))
  v = (min(Relaz(i,j),360.0-Relaz(i,j)) - az_lut(iv))/(az_lut(iv+1)-az_lut(iv))
  t =  (Solzen(i,j) - Solzen_lut(it))/(Solzen_lut(it+1)-Solzen_lut(it))

  ifail = 0
 
 
  !--  channel 1
  temp_Ref_lut = (1.0-t)*(1.0-u)*(1.0-v)*Ref_lut_Ch1_aer(it,iu,:,iv) +  &
                 (1.0-t)*(1.0-u)*(v)*Ref_lut_Ch1_aer(it,iu,:,iv+1) +  &
                 (t)*(1.0-u)*(1.0-v)*Ref_lut_Ch1_aer(it+1,iu,:,iv) +  &
                 (t)*(1.0-u)*(v)*Ref_lut_Ch1_aer(it+1,iu,:,iv+1) +  &
                 (1.0-t)*(u)*(1.0-v)*Ref_lut_Ch1_aer(it,iu+1,:,iv) +  &
                 (1.0-t)*(u)*(v)*Ref_lut_Ch1_aer(it,iu+1,:,iv+1) +  &
                 (t)*(u)*(1.0-v)*Ref_lut_Ch1_aer(it+1,iu+1,:,iv) +  &
                 (t)*(u)*(v)*Ref_lut_Ch1_aer(it+1,iu+1,:,iv+1)
  call AER_RET(ch(1)%Ref_Toa(i,j)/100.0, aot1(i,j))

  !--  channel 2
  temp_Ref_lut = (1.0-t)*(1.0-u)*(1.0-v)*Ref_lut_Ch2_aer(it,iu,:,iv) +  &
                 (1.0-t)*(1.0-u)*(v)*Ref_lut_Ch2_aer(it,iu,:,iv+1) +  &
                 (t)*(1.0-u)*(1.0-v)*Ref_lut_Ch2_aer(it+1,iu,:,iv) +  &
                 (t)*(1.0-u)*(v)*Ref_lut_Ch2_aer(it+1,iu,:,iv+1) +  &
                 (1.0-t)*(u)*(1.0-v)*Ref_lut_Ch2_aer(it,iu+1,:,iv) +  &
                 (1.0-t)*(u)*(v)*Ref_lut_Ch2_aer(it,iu+1,:,iv+1) +  &
                 (t)*(u)*(1.0-v)*Ref_lut_Ch2_aer(it+1,iu+1,:,iv) +  &
                 (t)*(u)*(v)*Ref_lut_Ch2_aer(it+1,iu+1,:,iv+1)
  call AER_RET(ch(2)%Ref_Toa(i,j)/100.0, aot2(i,j))


  if (Ch3a_On_AVHRR(j) == sym%YES) then
   !--  channel 3a
   temp_Ref_lut = (1.0-t)*(1.0-u)*(1.0-v)*Ref_lut_Ch3a_aer(it,iu,:,iv) +  &
                 (1.0-t)*(1.0-u)*(v)*Ref_lut_Ch3a_aer(it,iu,:,iv+1) +  &
                 (t)*(1.0-u)*(1.0-v)*Ref_lut_Ch3a_aer(it+1,iu,:,iv) +  &
                 (t)*(1.0-u)*(v)*Ref_lut_Ch3a_aer(it+1,iu,:,iv+1) +  &
                 (1.0-t)*(u)*(1.0-v)*Ref_lut_Ch3a_aer(it,iu+1,:,iv) +  &
                 (1.0-t)*(u)*(v)*Ref_lut_Ch3a_aer(it,iu+1,:,iv+1) +  &
                 (t)*(u)*(1.0-v)*Ref_lut_Ch3a_aer(it+1,iu+1,:,iv) +  &
                 (t)*(u)*(v)*Ref_lut_Ch3a_aer(it+1,iu+1,:,iv+1)
   call AER_RET(ch(6)%Ref_Toa(i,j)/100.0, aot3a(i,j))
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
  call LOCATE(temp_Ref_lut,ntau,y,is)
  is = max(1,min(ntau-1,is))
  s = (y - temp_Ref_lut(is))/(temp_Ref_lut(is+1) - temp_Ref_lut(is))
  aot = (1.0-s)*tau_lut(is) + s*tau_lut(is+1) 
end subroutine AER_RET

!------------------------------------------------------------------------------
! READ IN THE NESDIS OPERATIONAL LOOKUP TABLES FOR CH1, CH2 and CH3A OCEANIC AEROSOL
!------------------------------------------------------------------------------
 subroutine READ_AER_CH123A_REF_LUTS(ancil_data_dir,wmonum)
  character (len=*), intent(in):: ancil_data_dir
  integer, intent(in):: wmonum
  character (len=128):: ch1_aer_lut_file,ch2_aer_lut_file,ch3a_aer_lut_file
  character (len=128):: header,basename
  integer:: iSolzen,izen,itau
  integer:: nzen1,naz1,nSolzen1,ntau1
  integer:: nzen2,naz2,nSolzen2,ntau2
  integer:: nzen3,naz3,nSolzen3,ntau3
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

  ch1_aer_lut_file = trim(ancil_data_dir)//"luts/aerosol/"//trim(basename)//"_1"
  ch2_aer_lut_file = trim(ancil_data_dir)//"luts/aerosol/"//trim(basename)//"_2"
  ch3a_aer_lut_file = trim(ancil_data_dir)//"luts/aerosol/"//trim(basename)//"_3"

  if (wmonum >= 200 .and. wmonum <= 205) then
    ch3a_aer_lut_file = trim(ancil_data_dir)//"luts/aerosol/"//"Ch1F.lut.dat-NOAA15_3"
  endif
  if (wmonum >= 706 .and. wmonum <= 708) then
    ch3a_aer_lut_file = trim(ancil_data_dir)//"luts/aerosol/"//"Ch1F.lut.dat-NOAA15_3"
  endif

  ch1_lun = GET_LUN()
  open(unit=ch1_lun,file=ch1_aer_lut_file,status="old",position="rewind",action="read")
  ch2_lun = GET_LUN()
  open(unit=ch2_lun,file=ch2_aer_lut_file,status="old",position="rewind",action="read")
  ch3a_lun = GET_LUN()
  open(unit=ch3a_lun,file=ch3a_aer_lut_file,status="old",position="rewind",action="read")

   read(unit=ch1_lun,fmt="(4i5)") nzen1, naz1, nSolzen1,ntau1
   read(unit=ch2_lun,fmt="(4i5)") nzen2, naz2, nSolzen2,ntau2
   read(unit=ch3a_lun,fmt="(4i5)") nzen3, naz3, nSolzen3,ntau3

   Ref_lut_Ch1_aer = 0.0
   Ref_lut_Ch2_aer = 0.0
   Ref_lut_Ch3a_aer = 0.0
   tau_lut = 0.00
   Solzen_lut = 0.00
   zen_lut = 0.00
   az_lut = 0.00
   lambda_ref = 0.63
   delta_tau = 0.15

   temp_Ref_lut = 0.0

   delta_Solzen = 90.0 / nSolzen
   do iSolzen = 1,nSolzen
    Solzen_lut(iSolzen) = (iSolzen-1)*delta_Solzen
   enddo


   read(unit=ch1_lun, fmt=*) tau_lut
   read(unit=ch2_lun, fmt=*) tau_lut
   read(unit=ch3a_lun, fmt=*) tau_lut

   tau_lut = tau_lut * delta_tau

   do iSolzen = 1,nSolzen
       cos_Solzen = cos(Solzen_lut(iSolzen)*dtor)
     do itau = 1, ntau
       read(unit=ch1_lun,fmt="(a)") header
       read(unit=ch2_lun,fmt="(a)") header
       read(unit=ch3a_lun,fmt="(a)") header

       read(unit=ch1_lun,fmt="(10f6.0,/,9f6.0)") az_lut
       read(unit=ch2_lun,fmt="(10f6.0,/,9f6.0)") az_lut
       read(unit=ch3a_lun,fmt="(10f6.0,/,9f6.0)") az_lut
      
       do izen = 1, nzen
        read(unit=ch1_lun,fmt=*) zen_lut(izen), Ref_lut_Ch1_aer(iSolzen,izen,itau,1:5)
        read(unit=ch1_lun,fmt=*) Ref_lut_Ch1_aer(iSolzen,izen,itau,6:10)
        read(unit=ch1_lun,fmt=*) Ref_lut_Ch1_aer(iSolzen,izen,itau,11:15)
        read(unit=ch1_lun,fmt=*) Ref_lut_Ch1_aer(iSolzen,izen,itau,16:19)

        read(unit=ch2_lun,fmt=*) zen_lut(izen), Ref_lut_Ch2_aer(iSolzen,izen,itau,1:5)
        read(unit=ch2_lun,fmt=*) Ref_lut_Ch2_aer(iSolzen,izen,itau,6:10)
        read(unit=ch2_lun,fmt=*) Ref_lut_Ch2_aer(iSolzen,izen,itau,11:15)
        read(unit=ch2_lun,fmt=*) Ref_lut_Ch2_aer(iSolzen,izen,itau,16:19)

        read(unit=ch3a_lun,fmt=*) zen_lut(izen), Ref_lut_Ch3a_aer(iSolzen,izen,itau,1:5)
        read(unit=ch3a_lun,fmt=*) Ref_lut_Ch3a_aer(iSolzen,izen,itau,6:10)
        read(unit=ch3a_lun,fmt=*) Ref_lut_Ch3a_aer(iSolzen,izen,itau,11:15)
        read(unit=ch3a_lun,fmt=*) Ref_lut_Ch3a_aer(iSolzen,izen,itau,16:19)

       enddo 
     enddo
     Ref_lut_Ch1_aer(iSolzen,:,:,:) = Ref_lut_Ch1_aer(iSolzen,:,:,:)/cos_Solzen
     Ref_lut_Ch2_aer(iSolzen,:,:,:) = Ref_lut_Ch2_aer(iSolzen,:,:,:)/cos_Solzen
     Ref_lut_Ch3a_aer(iSolzen,:,:,:) = Ref_lut_Ch3a_aer(iSolzen,:,:,:)/cos_Solzen
   enddo

  close(unit=ch1_lun)
  close(unit=ch2_lun)
  close(unit=ch3a_lun)

   delta_zen = zen_lut(2) - zen_lut(1)
   delta_az = az_lut(2) - az_lut(1)
 end subroutine READ_AER_CH123A_REF_LUTS

end module AEROSOL_PROPERTIES
