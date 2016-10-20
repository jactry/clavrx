!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: algorithm_constants.f90 (src)
!       ALGORITHM_CONSTANTS (program)
!
! PURPOSE: This module serves as a common block for passing the 
!          non-cloud algorithm coefficients
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
! NOTES:
!
! File I/O - none
!
! Public routines: none
!
! Private routines: none
!
!-----------------------------------------------------------------------
module ALGORITHM_CONSTANTS

  use CONSTANTS

  implicit none

  real (kind=real4), parameter, public:: GLINT_ZEN_THRESH = 40.0
  real (kind=real4), parameter, public:: SCAT_ANGLE_THRESH = 90.0
  real (kind=real4), parameter, public:: STEFAN_BOLTZMANN_CONSTANT = 5.670e-08  !W/m^2/K^4
  real (kind=real4), parameter, public:: SOLAR_CONSTANT = 1360.0   !W/m^2

  !--- default modis white sky values for Snow and water
  real, parameter, public:: REF_SFC_WHITE_SKY_WATER = 5.0 

!--- End of module

end module ALGORITHM_CONSTANTS
