! $Id: $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: clavrx_olr_module.f90 (src)
!       clavrx_olr_module.f90 (program)
!
! PURPOSE: 
!       compute Outgoing Longwave Radiation (OLR) 
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
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
! File I/O: None, note this a clavrx module makes use of the clavrx global
!           memory in PIXEL_COMMON and CONSTANTS
!
! Public Routines:
!  COMPUTE_OLR 
!
! Private Routines:
!  AVHRR_OLR
!  GOES_OLR
!
!--------------------------------------------------------------------------------------
module CLAVRX_OLR_MODULE
use CONSTANTS
use PIXEL_COMMON, only: ch, sensor, geo, olr
use ALGORITHM_CONSTANTS

implicit none

public:: COMPUTE_OLR, SETUP_OLR
private:: SPLIT_WINDOW_OLR

real, dimension(10), save, private:: olr_coef

contains

!-----------------------------------------------------------------
!  set olr_coef vector for the sensor
!
!-----------------------------------------------------------------
subroutine SETUP_OLR()


!--- initialize
olr_coef = MISSING_VALUE_REAL4

!--- 
select case (Sensor%WMO_Id)

      case(4) !METOP-A
         olr_coef(1:4) = (/ -60.14468, 0.92820,  0.00532, 0.59607/)
         olr_coef(5:9) = (/ 108.31544, 0.76409, -5.04253, 1.50100, -0.00285/)

      case(3) !METOP-B
         olr_coef(1:4) = (/ -60.04060, 0.92773,  0.00535, 0.56909/)
         olr_coef(5:9) = (/ 108.48708, 0.76272, -4.89100, 1.45674, -0.00272/)

      case(5) !METOP-C

      case(55) !MSG-8

      case(56) !MSG-9

      case(57) !MSG-10

      case(171) !MTSAT-1R

      case(172) !MTSAT-2

      case(173) !AHI

      case(200) !NOAA-8 AVHRR

      case(201) !NOAA-9 AVHRR
         olr_coef(1:4) = (/ -59.94189, 0.92706,  0.00549, 0.63098/)
         olr_coef(5:9) = (/ 109.02869, 0.75937, -5.79289, 1.77363, -0.00270/)

      case(202) !NOAA-10 AVHRR

      case(203) !NOAA-11 AVHRR
         olr_coef(1:4) = (/-60.02730, 0.92756, 0.00540, 0.61423/)
         olr_coef(5:9) = (/108.63747, 0.76197, -5.46845, 1.64940, -0.00279/)

      case(204) !NOAA-12 AVHRR
         olr_coef(1:4) = (/-60.24969, 0.92879,  0.00522, 0.59537/)
         olr_coef(5:9) = (/107.55368, 0.76867, -4.90049, 1.41152, -0.00285/) 

      case(205) !NOAA-14 AVHRR
         olr_coef(1:4) = (/-60.36504, 0.92936, 0.00517, 0.54885/)
         olr_coef(5:9) = (/107.45860, 0.76914,-4.37248, 1.24385, -0.00278/) 

      case(206) !NOAA-15 AVHRR
         olr_coef(1:4) = (/-60.07516, 0.92778,  0.00539, 0.61436/)
         olr_coef(5:9) = (/108.56152, 0.76234, -5.30072, 1.59413, -0.00277/) 

      case(207) !NOAA-16 AVHRR
         olr_coef(1:4) = (/-60.42642, 0.92973,  0.00510, 0.55280/)
         olr_coef(5:9) = (/106.94872, 0.77204, -4.27699, 1.18709, -0.00271/) 

      case(208) !NOAA-17 AVHRR
         olr_coef(1:4) = (/-60.07043, 0.92775,  0.00540, 0.60629/)
         olr_coef(5:9) = (/108.59656, 0.76195, -5.20385, 1.56456, -0.00271/) 

      case(209) !NOAA-18 AVHRR
         olr_coef(1:4) = (/-60.48009, 0.92993,  0.00510, 0.53741/)
         olr_coef(5:9) = (/107.19191, 0.77077, -4.15499, 1.16309, -0.00276/) 

      case(223) !NOAA-19 AVHRR
         olr_coef(1:4) = (/-60.75243, 0.93132,  0.00497, 0.50949/)
         olr_coef(5:9) = (/106.38530, 0.77575, -3.74003, 1.00011, -0.00272/) 

      case(224) !NPP-VIIRS 
         olr_coef(1:4) = (/-59.97073, 0.92727,  0.00545, 0.60347/)
         olr_coef(5:9) = (/108.71840, 0.76120, -5.44672, 1.64156, -0.00269/) 

      case(252) !GOES-8

      case(253) !GOES-9

      case(254) !GOES-10

      case(255) !GOES-11

      case(256) !GOES-12

      case(257) !GOES-13

      case(258) !GOES-14

      case(259) !GOES-15

      case(706) !NOAA-6

      case(707) !NOAA-7
         olr_coef(1:4) = (/-60.09697, 0.92794,  0.00535, 0.59990/)
         olr_coef(5:9) = (/108.26270, 0.76421, -5.18247, 1.53830, -0.00279/) 

      case(708) !NOAA-5

      case(783) !TERRA-MODIS 
         olr_coef(1:4) = (/-60.10922, 0.92743,  0.00570, 0.75807/)
         olr_coef(5:9) = (/110.84430, 0.74922, -7.00103, 2.31716, -0.00280/) 

      case(784) !AQUA-MODIS 
         olr_coef(1:4) = (/-60.11890, 0.92744,  0.00572, 0.76198/)
         olr_coef(5:9) = (/110.90846, 0.74883, -7.06467, 2.34522,-0.00280/) 

      case(810) !COMS

      case(514) !FY2D      

      case(515) !FY2E          

      case default
         print*,'sensor for WMO number not found in OLR  ', Sensor%WMO_id
         print*,'OLR will not be computed'
      end select

end subroutine SETUP_OLR
!-----------------------------------------------------------------
! 
!-----------------------------------------------------------------
subroutine COMPUTE_OLR()

      Olr = MISSING_VALUE_REAL4

      select case (Sensor%WMO_Id)

      !--- AVHRR/2 and AVHRR/3 sensors
      case(3:5,201,203:209,223)
 
      if (Sensor%Chan_On_Flag_Default(31)==sym%YES .and. &
          Sensor%Chan_On_Flag_Default(32)==sym%YES) then

          Olr = SPLIT_WINDOW_OLR(ch(31)%Bt_Toa,ch(32)%Bt_Toa,Geo%Seczen)

          where(ch(31)%Bt_Toa == MISSING_VALUE_REAL4 .or.  &
                ch(32)%Bt_Toa == MISSING_VALUE_REAL4) 
                Olr = MISSING_VALUE_REAL4
          endwhere
 
      endif

      case default

      end select

end subroutine COMPUTE_OLR

!----------------------------------------------------------------------
! function to compute olr from split window observations
!----------------------------------------------------------------------
real elemental function split_window_olr(bt11,bt12,seczen)

  real, intent(in):: bt11
  real, intent(in):: bt12
  real, intent(in):: seczen

  !this makes the window eff flux temp (K) (win_Olr = sigma*win_temp^4)
  split_window_olr = olr_coef(1) +  &
                     olr_coef(2)*bt12  + &
                     olr_coef(3)*bt12*seczen  + &
                     olr_coef(4)*(bt11-bt12)*seczen 

  !--- this make the Olr temp (K)
  split_window_olr = olr_coef(5) + olr_coef(6)*split_window_olr + &
                     olr_coef(7)*(bt11-bt12) + &
                     olr_coef(8)*(bt11-bt12)*seczen + &
                     olr_coef(9)*split_window_olr**seczen

 !--- this make the Olr  flux (W/m^2)
 split_window_olr = stefan_boltzmann_constant * (split_window_olr)**4

end function split_window_olr

end module CLAVRX_OLR_MODULE
