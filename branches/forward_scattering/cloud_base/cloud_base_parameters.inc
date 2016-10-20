!$Id: awg_cld_hght_include_1.inc 405 2014-06-10 20:28:10Z heidinger $
!----------------------------------------------------------------------
!
! Name:
!   cloud_base_parameters.inc
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
! Description:
!   Include file that declares a number of physical constants,
!   apriori estimates and uncertainty values
!
!----------------------------------------------------------------------

 !--- set starting index value
 integer(kind=INT4), parameter, private:: Element_Idx_Min = 1
 integer(kind=INT4), parameter, private:: Line_Idx_Min = 1
 integer(kind=INT4), parameter, private:: LRC_Meander_Flag = 1

 !---------------------------------------------------------------------
 ! retrieval specific parameters
 !---------------------------------------------------------------------
 integer(kind=INT4), private, parameter :: Num_Levels_Rtm_Prof = 101

 !--- parameters that allow for cloud geometrical thickness estimation
 real(kind=real4), parameter, private:: WATER_EXTINCTION = 25.00   ! 1/km

 real(kind=real4), parameter, private:: ICE_EXTINCTION1 = 1.71
 real(kind=real4), parameter, private:: ICE_EXTINCTION2 = 1.87
 real(kind=real4), parameter, private:: ICE_EXTINCTION3 = 2.24
 real(kind=real4), parameter, private:: ICE_EXTINCTION4 = 2.88
 real(kind=real4), parameter, private:: ICE_EXTINCTION5 = 4.74

 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION1 = 0.13
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION2 = 0.25
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION3 = 0.39
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION4 = 0.55
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION5 = 0.67
