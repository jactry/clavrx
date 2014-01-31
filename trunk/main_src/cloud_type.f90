!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: cloud_type.f90 (src)
!       CLOUD_TYPING (program)
!
! PURPOSE: This module performs a cloud typing decision on pixel by pixel basis
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Michael Pavolonis (NOAA/NESDIS)
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
! REVISION HISTORY:
!   October 2006, Added retype routine - Heidinger
!   March 2012 - Added check for solar contamination flag
!
! Subroutines included in module: 
!	CLOUD_TYPE
!
! DEPENDENCIES:
!	CONSTANTS
!	PIXEL_COMMON
!
! SIDE EFFECTS: None
!--------------------------------------------------------------------------------------

module CLOUD_TYPING
 use CONSTANTS
 use PIXEL_COMMON
 use RTM_COMMON

 implicit none
 private
 public:: CLOUD_TYPE,  &
          CLOUD_RETYPE
 
 integer, parameter, private:: n_box=10

 contains

!------------------------------------------------------------------------
! SUBROUTINE NAME: CLOUD_TYPE
!
! This subroutine performs a cloud typing decision on pixel by pixel basis
! the resulting cloud mask codes are
!
! cloud type codes
!  0 - clear or partly cloudy
!  1 - fog
!  2 - liquid cloud
!  3 - mixed phase cloud
!  4 - glaciated
!  5 - cirrus cloud
!  6 - cirrus over lower
!
! INPUTS:
!  j1 - the first scan index to process
!  j2 - the last scan index to process
!
! OUTPUTS:
!  
! CALLING SEQUENCE:  call CLOUD_TYPE(j_min,num_scans_read)
!			(called from PROCESS_AVHRR_CLAVR)
!
! REFERENCES:
!
! DEPENDENCIES:
!	CONSTANTS
!	PIXEL_COMMON
!
! SIDE EFFECTS: None
!
! MODIFICATIONS: 
!   april 2005 - added spatial filtering to prevent ct =5,6 for isolated pixels
!                in regions where there are no cold pixels
!
!   August 30, 2005 - cleaned-up code and eliminated cirrus quality-based
!                     filter; only Bt_Ch31 filter is now used - mpav
!
!   April 14, 2006 - Modified to process on an arbitrary range of scans
!------------------------------------------------------------------------
!
! GLOBAL VARIABLES:
! lat				= latitude of pixel
! lon				= longitude of pixel
! sfc_type			= GFS land type
! cld_mask 			= cloud mask output
! cirrus_quality 		= quality of cirrus flag
! Solzen			= solar zenith angle
! satzen			= satellite zenith angle
! coszen			= cosine of the solar zenith angle
! bad_pixel_mask                = pixel quality flag
! ch3a_on_avhrr			= whether or not AVHRR channel 3a is used (0 = False, 1 = True)
! Bt_Ch31				= Channel 4 Brightness Tempertature
! Ref_Ch20				= channel 3b reflectance
! Ref_Ch6 			= channel 3a reflectance
! ems_Ch20 			= channel 3b emissivity
! Btd_Ch20_Ch31				= Brightness Tempurature channel 3b - Brightness Tempurature channel 4
! Btd_Ch31_Ch32 				= Brightness Tempurature channel 4 - Brightness Tempurature channel 5
!
! LOCAL VARIABLES:
!
! NIR_PHASE_THRES		= 1.6 micron phase threshold
! NIR_CIRRUS_THRES		= 1.6 micron cirrus threshold
! NIR_OVER_THRES		= Minimum 1.6 um reflectance allowed for cloud overlap over SNOW/ICE.
! BTD3811_PHASE_THRES		= 3.75um - 11um thresholds used for phase determination.
! EMS38_PHASE_THRES		= 3.75 um thresholds used for phase determination.
! BTD1112_DOVERLAP_THRES	= 11um - 12um cloud overlap threshold
! BTD1112_CIRRUS_THRES		= 11um - 12um cirrus threshold
! BTD1112_NOVERLAP_THRES_L	= Split window nighttime low cloud overlap threshold
! BTD1112_NOVERLAP_THRES_H	= Split window nighttime high cloud overlap threshold
! EMS38_NOVERLAP_THRES_L	= EMS38 nighttime low cloud overlap threshold
! EMS38_NOVERLAP_THRES_H	= EMS38 nighttime high cloud overlap threshold
! MIN_BTD1112_DOVERLAP		= The minimum 11um - 12um BTD allowed for overlap detection. 
! MIN_BTD1112_NOVERLAP		= The minimum allowed Bt_Ch31 - Bt_Ch32 allowed for nighttime overlap 
! A1				= Coefficient needed to determine the 11um - 12um BTD for cirrus detection
! B1				= Coefficient needed to determine the 11um - 12um BTD for cirrus detection
! C1				= Coefficient needed to determine the 11um - 12um BTD for cirrus detection
! D1				= Coefficient needed to determine the 11um - 12um BTD for cirrus detection
! E1				= Coefficient needed to determine the 11um - 12um BTD for cirrus detection
! A2				= Coefficient needed to determine the 3.75um - 11um BTD thresholds that
!		differentiate between ice and water as a function of 11um - 12um BTD.
! B2				= Coefficient needed to determine the 3.75um - 11um BTD thresholds that
!		differentiate between ice and water as a function of 11um - 12um BTD.
! C2				= Coefficient needed to determine the 3.75um - 11um BTD thresholds that
!		differentiate between ice and water as a function of 11um - 12um BTD.
! D2				= Coefficient needed to determine the 3.75um - 11um BTD thresholds that
!		differentiate between ice and water as a function of 11um - 12um BTD.
! E2				= Coefficient needed to determine the 3.75um - 11um BTD thresholds that
!		differentiate between ice and water as a function of 11um - 12um BTD.
! A3				= Coefficient needed to determine the 11um - 12um BTD used to find cloud
!		overlap as a function of 0.65 um reflectance.
! B3				= Coefficient needed to determine the 11um - 12um BTD used to find cloud
!		overlap as a function of 0.65 um reflectance.
! C3				= Coefficient needed to determine the 11um - 12um BTD used to find cloud
!		overlap as a function of 0.65 um reflectance.
! D3				= Coefficient needed to determine the 11um - 12um BTD used to find cloud
!		overlap as a function of 0.65 um reflectance.
! E3				= Coefficient needed to determine the 11um - 12um BTD used to find cloud
!		overlap as a function of 0.65 um reflectance.
!
! i				= pixel counter
! j				= scanline counter
! j1				= start scanline
! j2				= ending scanline
! index1			= viewing zenith angle bin
! index2			= solar zenith angle bin
! wflg				= IR window flag
! start_line			= line to start filtering 
! end_line			= line to end filtering
! start_pix			= pixel to start filtering
! end_pix			= pixel to end filtering
! npix				= number of pixels to filer
!
! day				= day/night flag
!
! t4_filter_thresh		= BT4 threshold, accounting for atmospheric effects, for filtering
! nir_ref			= channel 3a/b reflectance

! Not used: A4, B4, C4, D4, E4, n
!
!------------------------------------------------------------------------
 subroutine CLOUD_TYPE(j1,j2)

  integer, intent(in):: j1, j2

 !----------------------------------------------------------------------------
 ! Declare some variables to hold various thresholds.
 !----------------------------------------------------------------------------
 
 real (kind=real4):: NIR_PHASE_THRES, NIR_CIRRUS_THRES, NIR_OVER_THRES,&
                     BTD3811_PHASE_THRES, EMS38_PHASE_THRES,&
                     BTD1112_DOVERLAP_THRES, BTD1112_CIRRUS_THRES,&
                     BTD1112_NOVERLAP_THRES_L, BTD1112_NOVERLAP_THRES_H,&
                     EMS38_NOVERLAP_THRES_L, EMS38_NOVERLAP_THRES_H
     
 real(kind=real8), dimension(7):: A1
 real(kind=real8), dimension(7):: B1
 real(kind=real8), dimension(7):: C1
 real(kind=real8), dimension(7):: D1
 real(kind=real8), dimension(7):: E1
 real(kind=real8), dimension(7):: A2
 real(kind=real8), dimension(7):: B2
 real(kind=real8), dimension(7):: C2
 real(kind=real8), dimension(7):: D2
 real(kind=real8), dimension(7):: E2
 real(kind=real8), dimension(7,8) :: A3
 real(kind=real8), dimension(7,8) :: B3
 real(kind=real8), dimension(7,8) :: C3
 real(kind=real8), dimension(7,8) :: D3
 real(kind=real8), dimension(7,8) :: E3
 real(kind=real8), dimension(7,8) :: MIN_BTD1112_DOVERLAP
 real(kind=real8), dimension(7) :: MIN_BTD1112_NOVERLAP

 !----------------------------------------------------------------------------
 ! Declare some miscelaneous variables.
 !----------------------------------------------------------------------------
 
 integer:: i, j, index1, index2, wflg,&
           start_line, end_line, start_pix, end_pix, npix
 logical:: day
 real:: t4_filter_thresh, nir_ref

! -------------------------------------------------------------------------
! 					EXECUTABLE CODE
! -------------------------------------------------------------------------

! Fill coefficents

!----------------------------------------------------------------------------		
! Coefficients needed to determine the 11um - 12um BTD for cirrus detection 
! thresholds as a function of 11 um BT.  All coefficients are for 7
! different viewing zenith angles.
!----------------------------------------------------------------------------

 A1 = reshape((/-3.21578298543e+03, -2.94034788795e+03, -3.21256258031e+03,   &
                -3.47061397845e+03, -3.50486299804e+03, -5.08846627537e+03,   &
                -5.09507382017e+03/),(/7/))
           
 B1 = reshape((/ 4.88462861927e+01,  4.47332404792e+01,  4.86993920820e+01,   &
                 5.27678258867e+01,  5.32849278382e+01,  7.75358901756e+01,   &
                 7.80030571000e+01/),(/7/))
           
 C1 = reshape((/-2.76528129603e-01, -2.53525879101e-01, -2.75138895479e-01,   &
                -2.99072421299e-01, -3.01969990521e-01, -4.40956004022e-01,   &
                -4.45699597692e-01/),(/7/))
          
 D1 = reshape((/ 6.90692856276e-04,  6.33594236649e-04,  6.85787480971e-04,   &
                 7.48047893283e-04,  7.55160035017e-04,  1.10843453838e-03,   &
                 1.12561057802e-03/),(/7/))
           
 E1 = reshape((/-6.41179656790e-07, -5.88096462349e-07, -6.35205859105e-07,   &
                -6.95627948240e-07, -7.02035486892e-07, -1.03800396763e-06,   &
                -1.05900437811e-06/),(/7/))    
 
!----------------------------------------------------------------------------
! Coefficients needed to determine the 3.75um - 11um BTD thresholds that 
! differentiate between ice and water as a function of 11um - 12um BTD.  
! All coefficients are for 7 different viewing zenith angles.
!----------------------------------------------------------------------------

 A2 = reshape((/ 3.67107643481e-01,  3.36458606683e-01,  2.41676270282e-01,   &
                 1.74658763419e-01,  2.87704447664e-01,  2.17330552664e-02,   &
                -7.31011031866e-01/),(/7/))
           
 B2 = reshape((/ 3.78078114458e+00,  3.76731335795e+00,  3.76913819531e+00,   &
                 4.11018993353e+00,  4.47071677565e+00,  4.55202437527e+00,   &
                 4.97669171136e+00/),(/7/))
           
 C2 = reshape((/ 1.68542169247e+00,  1.69898191715e+00,  1.71234157395e+00,   &
                 1.50137445131e+00,  1.08453312818e+00,  1.02913912350e+00,   &
                 1.47244717791e+00/),(/7/))
          
 D2 = reshape((/-1.29320690316e+00, -1.27347975779e+00, -1.25712279459e+00,   &
                -1.24807113787e+00, -1.09555849488e+00, -1.00754783102e+00,   &
                -1.31340559947e+00/),(/7/))
           
 E2 = reshape((/ 1.82935664832e-01,  1.72289757442e-01,  1.66813773502e-01,   &
                 1.76165053220e-01,  1.53279656708e-01,  1.28047566835e-01,   &
                 1.74999237696e-01/),(/7/))

!----------------------------------------------------------------------------
! Coefficients needed to determine the 11um - 12um BTD used to find cloud
! overlap as a function of 0.65 um reflectance.  All coefficients are for 7
! different viewing zenith angles and 8 different solar zenith angles.
!----------------------------------------------------------------------------

 A3 = reshape((/&
       2.94e+00, 3.14e+00, 3.15e+00, 3.03e+00, 3.27e+00, 3.77e+00, 3.77e+00,  &
       2.94e+00, 3.14e+00, 3.15e+00, 3.03e+00, 3.27e+00, 3.77e+00, 3.77e+00,  &
       2.76e+00, 3.04e+00, 3.14e+00, 3.20e+00, 3.23e+00, 3.25e+00, 3.25e+00,  &
       2.95e+00, 2.75e+00, 3.03e+00, 3.15e+00, 3.34e+00, 3.48e+00, 3.48e+00,  &
       2.62e+00, 2.71e+00, 2.65e+00, 2.80e+00, 2.80e+00, 2.97e+00, 2.97e+00,  &
       2.26e+00, 2.59e+00, 2.33e+00, 2.43e+00, 2.62e+00, 3.01e+00, 3.01e+00,  &
       1.94e+00, 1.29e+00, 1.65e+00, 1.65e+00, 1.88e+00, 6.49e-01, 6.49e-01,  &
      -2.33e+00,-1.83e+00, 4.17e-01,-2.67e+00,-7.20e-01, 2.34e-01, 2.34e-01/),&
      (/7,8/))

 B3 = reshape((/&
       9.36e-01,-3.25e+00,-2.60e+00, 1.71e+00,-7.43e-01,-8.27e+00,-8.27e+00,  &
       9.36e-01,-3.25e+00,-2.60e+00, 1.71e+00,-7.43e-01,-8.27e+00,-8.27e+00,  &
       4.48e+00,-1.20e+00,-2.31e+00,-1.98e+00, 1.48e-01, 2.65e+00, 2.65e+00,  &
       3.65e-01, 4.94e+00,-2.40e-01,-1.20e+00,-2.60e+00,-2.09e+00,-2.09e+00,  &
       6.62e+00, 4.96e+00, 6.72e+00, 5.76e+00, 8.27e+00, 9.71e+00, 9.71e+00,  &
       1.21e+01, 6.67e+00, 1.12e+01, 1.05e+01, 9.62e+00, 7.24e+00, 7.24e+00,  &
       1.67e+01, 2.45e+01, 1.99e+01, 1.99e+01, 1.81e+01, 3.37e+01, 3.37e+01,  &
       6.92e+01, 6.26e+01, 3.53e+01, 6.57e+01, 4.58e+01, 3.56e+01, 3.56e+01/),&
       (/7,8/))

 C3 = reshape((/&
      -4.12e+01,-2.41e+01,-2.77e+01,-4.57e+01,-3.60e+01,-8.56e+00,-8.56e+00,  &
      -4.12e+01,-2.41e+01,-2.77e+01,-4.57e+01,-3.60e+01,-8.56e+00,-8.56e+00,  &
      -5.47e+01,-3.20e+01,-2.77e+01,-2.96e+01,-3.80e+01,-5.29e+01,-5.29e+01,  &
      -3.76e+01,-5.51e+01,-3.41e+01,-3.01e+01,-2.40e+01,-3.02e+01,-3.02e+01,  &
      -6.08e+01,-5.30e+01,-5.77e+01,-5.38e+01,-6.41e+01,-7.67e+01,-7.67e+01,  &
      -8.14e+01,-5.96e+01,-7.29e+01,-6.59e+01,-6.21e+01,-4.88e+01,-4.88e+01,  &
      -1.02e+02,-1.27e+02,-1.06e+02,-1.00e+02,-9.37e+01,-1.38e+02,-1.38e+02,  &
      -3.09e+02,-2.80e+02,-1.69e+02,-2.56e+02,-1.86e+02,-1.23e+02,-1.23e+02/),&
      (/7,8/))

 D3 = reshape((/&
       8.58e+01, 6.05e+01, 6.69e+01, 9.35e+01, 7.91e+01, 4.24e+01, 4.24e+01,  &
       8.58e+01, 6.05e+01, 6.69e+01, 9.35e+01, 7.91e+01, 4.24e+01, 4.24e+01,  &
       1.05e+02, 7.15e+01, 6.50e+01, 6.78e+01, 7.92e+01, 1.07e+02, 1.07e+02,  &
       7.88e+01, 1.03e+02, 7.16e+01, 6.47e+01, 5.46e+01, 6.82e+01, 6.82e+01,  &
       1.11e+02, 9.79e+01, 1.01e+02, 9.37e+01, 1.08e+02, 1.33e+02, 1.33e+02,  &
       1.41e+02, 1.08e+02, 1.20e+02, 1.03e+02, 9.63e+01, 7.01e+01, 7.01e+01,  &
       1.78e+02, 2.08e+02, 1.70e+02, 1.53e+02, 1.43e+02, 1.87e+02, 1.87e+02,  &
       5.08e+02, 4.55e+02, 2.75e+02, 3.69e+02, 2.66e+02, 1.40e+02, 1.40e+02/),&
       (/7,8/))


 E3 = reshape((/&
      -5.09e+01,-3.83e+01,-4.20e+01,-5.52e+01,-4.80e+01,-3.16e+01,-3.16e+01,  &
      -5.09e+01,-3.83e+01,-4.20e+01,-5.52e+01,-4.80e+01,-3.16e+01,-3.16e+01,  &
      -6.02e+01,-4.36e+01,-4.01e+01,-4.14e+01,-4.67e+01,-6.36e+01,-6.36e+01,  &
      -4.67e+01,-5.81e+01,-4.20e+01,-3.79e+01,-3.25e+01,-4.11e+01,-4.11e+01,  &
      -6.25e+01,-5.49e+01,-5.42e+01,-4.99e+01,-5.63e+01,-7.13e+01,-7.13e+01,  &
      -7.72e+01,-6.02e+01,-6.28e+01,-5.13e+01,-4.70e+01,-3.22e+01,-3.22e+01,  &
      -1.00e+02,-1.12e+02,-8.96e+01,-7.65e+01,-7.04e+01,-8.40e+01,-8.40e+01,  &
      -2.85e+02,-2.52e+02,-1.49e+02,-1.82e+02,-1.28e+02,-5.21e+01,-5.21e+01/),&
      (/7,8/))

       
!----------------------------------------------------------------------------
! The minimum 11um - 12um BTD allowed for overlap detection.  All values are 
! for 7 different viewing zenith angles and 8 different solar zenith angles.
!----------------------------------------------------------------------------
	     
  MIN_BTD1112_DOVERLAP = reshape((/&
             0.70,0.70,0.70,0.70,0.75,0.80,0.80,  &
             0.70,0.70,0.70,0.70,0.75,0.80,0.80,  &
             0.70,0.70,0.70,0.70,0.75,0.80,0.80,  &
             0.70,0.70,0.70,0.70,0.75,0.80,0.80,  &
             0.70,0.70,0.70,0.70,0.75,0.80,0.80,  &
             0.70,0.70,0.70,0.70,0.75,0.90,0.90,  &
             0.75,0.75,0.75,0.80,0.80,0.90,0.90,  &
             0.75,0.75,0.75,0.80,0.80,0.90,0.90/),&
	     (/7,8/))

  
!----------------------------------------------------------------------------
! The minimum allowed Bt_Ch31 - Bt_Ch32 allowed for nighttime overlap for 7
! viewing angles.
!----------------------------------------------------------------------------

  MIN_BTD1112_NOVERLAP = reshape((/0.57948485,0.57948485,0.57948485, &
                             0.57948485,0.57948485,0.57948485, &
                             0.57948485/),(/7/))
 
! Begin cloud typing.

! loop over all pixels
i_loop: do  i = 1,num_pix

!loop over scanlines
 j_loop: do j = j1,j2 + j1 - 1
   
!---- check for bad pixels quality
    if (bad_pixel_mask(i,j) == sym%YES) then
     cycle
    endif

!print *, "Calling cloud_type ", Ch3a_On(j), ch(6)%Ref_Toa(i,j), ch(20)%Ref_Toa(i,j), ch(1)%Ref_Toa(i,j)
!---- check for lines without ch3a or ch3b
 if (Ch3a_On_AVHRR(j) == -1) then
     cycle
 endif

!-----
nir_ref = Missing_Value_Real4

! Determine the viewing zenith angle bin.
   index1 = min(7,max(1,int(satzen(i,j)/10.0) + 1))
   
! Determine the solar zenith angle bin.
   index2 = min(8,max(1,int(Solzen(i,j)/10.0) + 1))
     
! Set 11um - 12um cirrus thresholds.
   BTD1112_CIRRUS_THRES = A1(index1) + B1(index1)*ch(31)%Bt_Toa(i,j) + C1(index1)*ch(31)%Bt_Toa(i,j)**2 + &
                      D1(index1)*ch(31)%Bt_Toa(i,j)**3 + E1(index1)*ch(31)%Bt_Toa(i,j)**4

   BTD1112_CIRRUS_THRES = max(1.0,min(4.0,BTD1112_CIRRUS_THRES))

! Check if daytime or nighttime algorithm is to be used.
   day = .false.

   if (Solzen(i,j) < 88.0) then
     day = .true.
   endif 
   
! Check for clear or partly cloudy grid cells.
   if (cld_mask(i,j) == sym%CLEAR) then
     cld_type(i,j) = sym%CLEAR_TYPE
   endif
   if (cld_mask(i,j) == sym%PROB_CLEAR) then
     cld_type(i,j) = sym%PROB_CLEAR_TYPE
   endif
   
   cirrus_quality(i,j) = 0 
   
! If daytime, use daytime algorithm. So, check to see if it is day
   if ((cld_mask(i,j) > 1) .and. (day .eqv. .true.)) then
   
! Set 11um - 12um overlap thresholds.
! Note that the channel reflectance must be in % form.
     if ((ch(1)%Ref_Toa(i,j) >= 35.0) .and. (ch(1)%Ref_Toa(i,j) <= 60.0)) then
       BTD1112_DOVERLAP_THRES = max(((A3(index1,index2) + B3(index1,index2)*ch(1)%Ref_Toa(i,j)*0.01 + &
                                 C3(index1,index2)*(ch(1)%Ref_Toa(i,j)*0.01)**2 + &
                                 D3(index1,index2)*(ch(1)%Ref_Toa(i,j)*0.01)**3 + &
                                 E3(index1,index2)*(ch(1)%Ref_Toa(i,j)*0.01)**4) - 0.1),MIN_BTD1112_DOVERLAP(index1,index2) - 0.1)
     elseif (ch(1)%Ref_Toa(i,j) > 60.0 .and. ch(1)%Ref_Toa(i,j) < 90.0) then
       BTD1112_DOVERLAP_THRES = MIN_BTD1112_DOVERLAP(index1,index2) - 0.1
     else
       BTD1112_DOVERLAP_THRES = 9999.0
     endif

! In the high latitudes, the 3.75 um reflectance must be less than 20% to
! prevent single layer water clouds from being typed as overlap.
   if ((lat(i,j) > 65.0 .or. lat(i,j) < -65.0) .and. ch(20)%Ref_Toa(i,j) > 20.0) then
     BTD1112_DOVERLAP_THRES = 9999.0
   endif
   
! Perform initial IR window brightness temperature-based typing.
     wflg = 1
     
     if (ch(31)%Bt_Toa(i,j) <= 233.16) then
       cld_type(i,j) = sym%OPAQUE_ICE_TYPE
       wflg = 0
     elseif ((ch(31)%Bt_Toa(i,j) > 233.16) .and. (ch(31)%Bt_Toa(i,j) <= 253.16)) then
       cld_type(i,j) = sym%OPAQUE_ICE_TYPE
     elseif ((ch(31)%Bt_Toa(i,j) > 253.16) .and. (ch(31)%Bt_Toa(i,j) <= 273.16)) then
       cld_type(i,j) =sym%SUPERCOOLED_TYPE
     else
       cld_type(i,j) = sym%WATER_TYPE
     endif
     
! Use 1.6 um algorithm, if available.
     if ((ch3a_on_avhrr(j) == 1) .and. (ch(1)%Ref_Toa(i,j) > 0.0)) then 
   
!----------------------------------------------------------------------------
! Set some 1.6 um thresholds used in phase identification and cirrus detection.
! 0:  WATER
! 14: SNOW/ICE
! 12: DESERT
!----------------------------------------------------------------------------

       if ((sfc_type(i,j) == 0) .or. (sfc_type(i,j) == 14) ) then
         NIR_CIRRUS_THRES = 20.0
         NIR_PHASE_THRES = 17.0
       elseif (sfc_type(i,j) == 12) then
         NIR_CIRRUS_THRES = 55.0
         NIR_PHASE_THRES = 32.0 
       else
         NIR_CIRRUS_THRES = 33.0
         NIR_PHASE_THRES = 32.0
       endif
       
!----------------------------------------------------------------------------
! Set the minimum 1.6 um reflectance allowed for cloud overlap over SNOW/ICE.
! 14: SNOW/ICE
!----------------------------------------------------------------------------

       if (sfc_type(i,j) == 14) then
         NIR_OVER_THRES = 17.0
       else
         NIR_OVER_THRES = 0.0
       endif

! The reflectance used in the typing tests to 1.65 um when ch3a is on.
       nir_ref = ch(6)%Ref_Toa(i,j)
     
! Use 3.75 um channel (used only if 1.6um is not available).
     elseif ((ch3a_on_avhrr(j) == 0) .and. (ch(1)%Ref_Toa(i,j) > 0.0)) then
       
!----------------------------------------------------------------------------
! Set some 3.75 um thresholds used in phase identification and cirrus detection.
! 0:  WATER
! 14: SNOW/ICE
! 12: DESERT
!----------------------------------------------------------------------------
       if ((sfc_type(i,j) == 0) .or. (sfc_type(i,j) == 14) ) then
         NIR_CIRRUS_THRES = 12.0
         NIR_PHASE_THRES = 6.0
       elseif (sfc_type(i,j) == 12) then
         NIR_CIRRUS_THRES = 40.0
         NIR_PHASE_THRES = 6.0   
       else
         NIR_CIRRUS_THRES = 12.0
         NIR_PHASE_THRES = 6.0
       endif
       
!----------------------------------------------------------------------------
! Set the minimum 3.75 um reflectance allowed for cloud overlap over SNOW/ICE.
! 14: SNOW/ICE
!----------------------------------------------------------------------------
       if (sfc_type(i,j) == 14) then
         NIR_OVER_THRES = 6.0
       else
         NIR_OVER_THRES = 0.0
       endif

! The reflectance used in the typing tests to 3.75 um when ch3b is on.
       nir_ref = ch(20)%Ref_Toa(i,j) 
       
     endif
     
! Perform the NIR reflectance bulk cloud phase test.
     if (cld_type(i,j) == sym%SUPERCOOLED_TYPE .and. nir_ref <= NIR_PHASE_THRES .and. ch(31)%Bt_Toa(i,j) < 263.16) then
       cld_type(i,j) = sym%OPAQUE_ICE_TYPE
     endif

     if (cld_type(i,j) == sym%OPAQUE_ICE_TYPE .and. wflg == 1 .and. nir_ref > NIR_PHASE_THRES) then
       cld_type(i,j) = sym%SUPERCOOLED_TYPE
     endif
       
!----------------------------------------------------------------------------
! Perform the cloud overlap test - not used over DESERT surfaces.
! 12: DESERT
!----------------------------------------------------------------------------
     if((Btd_Ch31_Ch32(i,j) > BTD1112_DOVERLAP_THRES) .and. (ch(31)%Bt_Toa(i,j) < 270.0) .and. &
         sfc_type(i,j) /= 12 .and. nir_ref > NIR_OVER_THRES .and. ch(31)%Bt_Toa(i,j) > 210.0) then
       cld_type(i,j) = sym%OVERLAP_TYPE
     endif

!----------------------------------------------------------------------------
! Look for cirrus clouds.
!----------------------------------------------------------------------------
!-- note, akh modified so that nir_ref test only applied when Solzen < 70
     if ((cld_type(i,j) /= sym%OVERLAP_TYPE) .and. (Btd_Ch31_Ch32(i,j) > (BTD1112_CIRRUS_THRES-0.2)) .and. &
         ch(31)%Bt_Toa(i,j) < 295.0) then

         if (Solzen(i,j) < 70.0) then
            if (nir_ref < NIR_CIRRUS_THRES) then
             cld_type(i,j) = sym%CIRRUS_TYPE
             cirrus_quality(i,j) = 1
            endif

         else
           cld_type(i,j) = sym%CIRRUS_TYPE
           cirrus_quality(i,j) = 0   !note, this is a low quality

         endif

     endif   

!----------------------------------------------------------------------------
! Look for fog - not used over DESERT surfaces.
! The ref3b/ref1 threshold is used to prevent near terminator and 
! sunglint pixels from being classified as fog.
! THERE IS CURRENTLY NO FOG DETECTION WHEN CH3A IS ON!
! 12: DESERT
!----------------------------------------------------------------------------
   if (ch3a_on_avhrr(j) == 0 .and. ch(20)%Ref_Toa(i,j) >= 25.0 .and. &
       sfc_type(i,j) /= 12 .and. ch(31)%Bt_Toa(i,j) > 240.0 .and. &
       ch(20)%Ref_Toa(i,j)/ch(1)%Ref_Toa(i,j) < 0.6) then

     cld_type(i,j) = sym%FOG_TYPE

   endif
         
! If nighttime, use tri-spectral nighttime algorithm.
   elseif((cld_mask(i,j) > 1) .and. (day .eqv. .false.)) then
   
! Set 3.75um - 11um thresholds used for phase determination.
     BTD3811_PHASE_THRES = A2(index1) + B2(index1)*Btd_Ch31_Ch32(i,j) + C2(index1)*Btd_Ch31_Ch32(i,j)**2 + &
                           D2(index1)*Btd_Ch31_Ch32(i,j)**3 + E2(index1)*Btd_Ch31_Ch32(i,j)**4

     BTD3811_PHASE_THRES = min(8.0,max(-2.0,BTD3811_PHASE_THRES))
     
!----------------------------------------------------------------------------
! Set the 3.75 um thresholds used for phase determination.
! 0: WATER
!----------------------------------------------------------------------------
     if (ch(31)%Bt_Toa(i,j) <= 245.0) then

       if (sfc_type(i,j) == 0) then
         EMS38_PHASE_THRES = 0.9
       else
         EMS38_PHASE_THRES = 0.9
       endif

     else

       if (sfc_type(i,j) == 0) then
         EMS38_PHASE_THRES = 1.12
       else
         EMS38_PHASE_THRES = 1.12
       endif

     endif

!----------------------------------------------------------------------------
! Set the split window and EMS3b thresholds used in nighttime cloud
! overlap detection.
! 0:  WATER
!----------------------------------------------------------------------------
     
     if (Btd_Ch20_Ch31(i,j) > 0.0) then

       if (sfc_type(i,j) == 0) then

         if (lat(i,j) > -30.0 .and. lat(i,j) < 30.0) then
           BTD1112_NOVERLAP_THRES_H = 2.5
	   BTD1112_NOVERLAP_THRES_L = MIN_BTD1112_NOVERLAP(index1)+0.2
           EMS38_NOVERLAP_THRES_H = 5.0
	   EMS38_NOVERLAP_THRES_L = 1.1

         else
           BTD1112_NOVERLAP_THRES_H = 2.0
	   BTD1112_NOVERLAP_THRES_L = MIN_BTD1112_NOVERLAP(index1)
           EMS38_NOVERLAP_THRES_H = 2.5
	   EMS38_NOVERLAP_THRES_L = 1.05

         endif

       else

         if (lat(i,j) > -30.0 .and. lat(i,j) < 30.0) then
	   BTD1112_NOVERLAP_THRES_H = 2.5
	   BTD1112_NOVERLAP_THRES_L = MIN_BTD1112_NOVERLAP(index1)+0.2
           EMS38_NOVERLAP_THRES_H = 5.0
	   EMS38_NOVERLAP_THRES_L = 1.1

         else
           BTD1112_NOVERLAP_THRES_H = 2.0
	   BTD1112_NOVERLAP_THRES_L = MIN_BTD1112_NOVERLAP(index1)
           EMS38_NOVERLAP_THRES_H = 2.0
	   EMS38_NOVERLAP_THRES_L = 1.0

	   endif

       endif

     else
       BTD1112_NOVERLAP_THRES_H = -99.0
       BTD1112_NOVERLAP_THRES_L = 999.0
       EMS38_NOVERLAP_THRES_H = -99.0
       EMS38_NOVERLAP_THRES_L = 999.0

     endif 

     if (sfc_type(i,j) == 12) then
       BTD1112_NOVERLAP_THRES_H = -99.0
       BTD1112_NOVERLAP_THRES_L = 999.0
       EMS38_NOVERLAP_THRES_H = -99.0
       EMS38_NOVERLAP_THRES_L = 999.0
     endif

! Perform initial IR window brightness temperature-based typing.
     wflg = 1
     
     if (ch(31)%Bt_Toa(i,j) <= 233.16) then
       cld_type(i,j) = sym%OPAQUE_ICE_TYPE
       wflg = 0
     elseif ((ch(31)%Bt_Toa(i,j) > 233.16) .and. (ch(31)%Bt_Toa(i,j) <= 253.16)) then
       cld_type(i,j) = sym%OPAQUE_ICE_TYPE
     elseif ((ch(31)%Bt_Toa(i,j) > 253.16) .and. (ch(31)%Bt_Toa(i,j) <= 273.16)) then
       cld_type(i,j) = sym%SUPERCOOLED_TYPE
     else
       cld_type(i,j) = sym%WATER_TYPE
     endif  
       
! Perform the EMS 3.75 um test for bulk cloud phase determination.
     if (Solar_Contamination_Mask(i,j) == sym%NO) then
      if ((cld_type(i,j) == sym%SUPERCOOLED_TYPE .and. ems_Ch20(i,j) >= EMS38_PHASE_THRES .and. ch(31)%Bt_Toa(i,j) < 263.16)) then
       cld_type(i,j) = sym%OPAQUE_ICE_TYPE
      endif

      if ((cld_type(i,j) == sym%OPAQUE_ICE_TYPE .and. wflg == 1 .and. ems_Ch20(i,j) < EMS38_PHASE_THRES)) then
       cld_type(i,j) = sym%SUPERCOOLED_TYPE
      endif
     endif

! Perform 3.75um/11um/12um trispectral test with the EMS 3.75 um test for
! bulk cloud phase determination.
! CURRENTLY NOT USED!
!     if ((cld_type(i,j) == 3 .and. Btd_Ch20_Ch31(i,j) >= BTD3811_PHASE_THRES .and. ch(31)%Bt_Toa(i,j) < 263.16) .or. &
!         (cld_type(i,j) == 3 .and. ems_Ch20(i,j) >= 0.9 .and. ch(31)%Bt_Toa(i,j) < 263.16)) then
!       cld_type(i,j) = 4
!     endif
!
!     if ((cld_type(i,j) == 4 .and. wflg == 1 .and. Btd_Ch20_Ch31(i,j) < BTD3811_PHASE_THRES) .or. &
!         (cld_type(i,j) == 4 .and. wflg == 1 .and. ems_Ch20(i,j) < 0.9)) then
!       cld_type(i,j) = 3
!     endif

! Nighttime cloud overlap test.
     if (Btd_Ch31_Ch32(i,j) > BTD1112_NOVERLAP_THRES_L .and. Btd_Ch31_Ch32(i,j) < BTD1112_NOVERLAP_THRES_H .and. &
         ems_Ch20(i,j) > EMS38_NOVERLAP_THRES_L .and. ems_Ch20(i,j) < EMS38_NOVERLAP_THRES_H .and. &
	 ch(31)%Bt_Toa(i,j) > 210.0 .and. ch(31)%Bt_Toa(i,j) < 283.0) then

       cld_type(i,j) = sym%OVERLAP_TYPE

     endif
     
! Look for cirrus clouds using the split window test.   
     if ((cld_type(i,j) /= sym%OVERLAP_TYPE) .and. (Btd_Ch31_Ch32(i,j) > BTD1112_CIRRUS_THRES) .and. &
          ems_Ch20(i,j) > 1.3) then

       cld_type(i,j) = sym%CIRRUS_TYPE

     endif
     
     if (ems_Ch20(i,j) > 1.6 .and. ch(31)%Bt_Toa(i,j) < 300.0) then
         cirrus_quality(i,j) = 1
     endif
     
! Look for cirrus clouds using the EMS 3.75 um test.   
     if (Solar_Contamination_Mask(i,j) == sym%NO) then

      if (cld_type(i,j) /= sym%OVERLAP_TYPE .and. cld_type(i,j) /= sym%OPAQUE_ICE_TYPE .and. &
         ems_Ch20(i,j) > 1.10 .and. ch(31)%Bt_Toa(i,j) < 300.0) then

       cld_type(i,j) = sym%CIRRUS_TYPE

      endif
     
      if ((ems_Ch20(i,j) > 1.6 .and. ch(31)%Bt_Toa(i,j) < 300.0) .or. &
         (ems_Ch20(i,j) > 1.4 .and. ch(31)%Bt_Toa(i,j) < 300.0 .and. Btd_Ch31_Ch32(i,j) > BTD1112_CIRRUS_THRES)) then

       cirrus_quality(i,j) = 1
 
      endif

     endif

!----------------------------------------------------------------------------
! Look for fog - not used over DESERT surfaces.
! 12: DESERT
!----------------------------------------------------------------------------
     
     if (Solzen(i,j) >= 90.0 .and. ems_Ch20(i,j) <= 0.90 .and. &
         ch(31)%Bt_Toa(i,j) > 240.0 .and. sfc_type(i,j) /= 12) then

       cld_type(i,j) = sym%FOG_TYPE

     endif
       
   endif

!end scanline and pixel loops        
 end do j_loop
   
end do i_loop

!----------------------------------------------------------------------------
! END OF TYPING SPECTRAL TESTS.
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
! Start spatial filter used on cirrus and overlap pixels.
! A 2n_box x 2n_box pixel filter is employed.
!----------------------------------------------------------------------------

!loop over pixels and scanlines
i_loop2: do  i = 1,num_pix
  j_loop2: do j = j1,j2

!---- check for bad pixels
    if (bad_pixel_mask(i,j) == sym%YES) then
     cycle
    endif

!--- determine box for filtering
    start_line = max(j1,j - n_box)
    end_line = min(j2,j + n_box)
    start_pix = max(1,i - n_box)
    end_pix = min(num_pix,i + n_box)
    npix = ((end_line - start_line)+1)*((end_pix - start_pix)+1)

! At least one pixel in the 2n_box x 2n_box array must have a Bt_Ch31 < 295 K and
! the average ems 3.75 um must be < 1.2 for low quality cirrus; otherwise 
! the pixel is reset to water or mixed.  
    if (cld_type(i,j) == sym%CIRRUS_TYPE .and. cirrus_quality(i,j) == 0) then

!account for atmospheric effects
      t4_filter_thresh = 295.0 - 12.0*(1.0-coszen(i,j))

!filter pixels in region
!-- added criterion for missing values in ems_Ch20 - akh 1/07
      if ( (minval(ch(31)%Bt_Toa(start_pix:end_pix,start_line:end_line)) > t4_filter_thresh) .or. &          
           ( (sum(ems_Ch20(start_pix:end_pix,start_line:end_line))/npix < 1.2) .and. &
             (minval(ems_Ch20(start_pix:end_pix,start_line:end_line)) > 0.0)) )then

        if (ch(31)%Bt_Toa(i,j) <= 273.16) then
          cld_type(i,j) = sym%SUPERCOOLED_TYPE

        else
          cld_type(i,j) = sym%WATER_TYPE

        endif

      endif

    endif
    
! At least one pixel in the 2n x 2n array must have a Bt_Ch31 < 275 K for overlap;
! otherwise the pixel is reset to water or mixed.

    if (cld_type(i,j) == sym%OVERLAP_TYPE .and. Solzen(i,j) > 90.0) then
!account for atmospheric effects in BT4 threshold for filtering
      t4_filter_thresh = 273.0 - 12.0*(1.0-coszen(i,j))

      if (minval(ch(31)%Bt_Toa(start_pix:end_pix,start_line:end_line)) > t4_filter_thresh) then

        if (ch(31)%Bt_Toa(i,j) <= 273.16) then
           cld_type(i,j) = sym%SUPERCOOLED_TYPE
        else
           cld_type(i,j) = sym%WATER_TYPE
        endif

      endif

    endif

!end scanline and pixel loop
  end do j_loop2
end do i_loop2

!----------------------------------------------------------------------------
! End of subroutine CLOUD_TYPE
!----------------------------------------------------------------------------
  end subroutine CLOUD_TYPE


!----------------------------------------------------------------------------
! Beginning of subroutine CLOUD_RETYPE
!  This routine modifies the cloud_type array based on spatial
!  analysis.  Its goal is to reduce the edges of stratus clouds
!  from being typed as cirrus
!----------------------------------------------------------------------------
  subroutine CLOUD_RETYPE(jmin,numj,cld_type_array)
   integer (kind=int1), intent(inout), dimension(:,:):: cld_type_array
   integer, intent(in):: jmin
   integer, intent(in):: numj
   integer:: ilrc
   integer:: jlrc
   integer:: i
   integer:: j

   i_loop: do i = 1, num_pix

      j_loop: do j = jmin, numj + jmin - 1

      ilrc = i_lrc(i,j)
      jlrc = j_lrc(i,j)

      !--- check for ice-phase pixels where the local radiative center
      !--- is a water cloud  - retype these as water

      if ((cld_type_array(i,j) == sym%CIRRUS_TYPE)  .or. &
          (cld_type_array(i,j) == sym%OVERLAP_TYPE)) then

!         ilrc = i_min_Bt_Ch31_3x3(i,j)
!         jlrc = j_min_Bt_Ch31_3x3(i,j)
          ilrc = i_lrc(i,j)
          jlrc = j_lrc(i,j)

          !-- skip this if no lrc is available
          if (ilrc < 1 .or. jlrc < 1) then
             cycle
          endif

          if ((cld_type_array(ilrc,jlrc) == sym%FOG_TYPE) .or. &
              (cld_type_array(ilrc,jlrc) == sym%WATER_TYPE) .or. &
              (cld_type_array(ilrc,jlrc) == sym%SUPERCOOLED_TYPE)) then

              cld_type_array(i,j) = cld_type_array(ilrc,jlrc)

          endif
      endif

      !--------------------------------------------
      !check for water clouds on the edge of cirrus
      !--------------------------------------------
      if ((cld_type_array(i,j) == sym%FOG_TYPE)  .or. &
          (cld_type_array(i,j) == sym%WATER_TYPE)) then

          ilrc = i_lrc(i,j)
          jlrc = j_lrc(i,j)

          if (ilrc < 1 .or. jlrc < 1) then
                 cycle
          endif

          if ((cld_type_array(ilrc,jlrc) == sym%CIRRUS_TYPE) .or. &
              (cld_type_array(ilrc,jlrc) == sym%OVERLAP_TYPE) .or. &
              (cld_type_array(ilrc,jlrc) == sym%OPAQUE_ICE_TYPE)) then

              cld_type_array(i,j) = sym%CIRRUS_TYPE

          endif
      endif

      end do j_loop

   end do i_loop

  end subroutine CLOUD_RETYPE

!----------------------------------------------------------------------------
! End of module CLOUD_TYPING
!----------------------------------------------------------------------------
end module CLOUD_TYPING
