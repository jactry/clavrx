! $Id: avhrr_module.f90,v 1.5.2.2 2014/01/26 04:48:31 heidinger Exp $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: avhrr_module.f90 (src)
!       AVHRR_MODULE (program)
!
! PURPOSE: AVHRR calibration and navigation routines
!
! DESCRIPTION: EXPERMINENTAL VERSION = 
!              READS LEVEL-1b one segment at a time - not 2 scans
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
! Routines in this module and their purpose:
!
! DEFINE_1B_DATA - based on level 1b data type, define some parameters
!                  needed to read and unpack the avhrr data
! DETERMINE_AVHRR_FILE_TYPE - interrogate header and determine the file type
! UNPACK_AVHRR_HEADER_RECORD - unpack level 1b header for pre-AVHRR_KLM_Flag data
! UNPACK_AVHRR_DATA_RECORD - unpack level 1b data for pre-AVHRR_KLM_Flag data
! LAGRANGIAN_ANCHOR_INTERP - perform Lagrangian interpolatation of Anchors 
! LINEAR_ANCHOR_INTERP - perform linear interpolation of Anchors
! COMPUTE_ANGLE_ANCHORS - compute Satzen and Relaz Anchors for pre-AVHRR_KLM_Flag data
! UNPACK_AVHRR_HEADER_RECORD_KLM - unpack level 1b header for AVHRR_KLM_Flag+ data
! UNPACK_AVHRR_DATA_RECORD_KLM - unpack level 1b data for AVHRR_KLM_Flag+ data
! i4word_to_string - convert a 4 byte integer to a string
! READ_AVHRR_INSTR_CONSTANTS - read avhrr instrument constant files
! REF_CAL - perform the reflectance calibration
! THERM_CAL - perform the thermal calibration
! COMPUTE_NEW_THERM_CAL_COEF - compute thermal calibration coefficients
! WRITE_HEADER_1BX - write a header to the level 1bx files
! READ_AVHRR_LEVEL1B_DATA - read level 1b data records
! READ_AVHRR_LEVEL1B_HEADER - read level 1b header records
! CALCULATE_ASC_DES - Calculates ascending/descending flag if error presennt
!                       in level1b file
!
!--------------------------------------------------------------------------------------
module AVHRR_MODULE

  use PIXEL_COMMON
  use CALIBRATION_CONSTANTS
  use CONSTANTS
  use FILE_UTILITY
  use VIEWING_GEOMETRY_MODULE
  use PLANCK

  implicit none

  public:: WRITE_HEADER_1BX,&
           ASSIGN_AVHRR_SAT_ID_NUM_INTERNAL, &
           READ_AVHRR_INSTR_CONSTANTS, &
           READ_AVHRR_LEVEL1B_DATA,  &
           READ_AVHRR_LEVEL1B_HEADER, &
           DETERMINE_AVHRR_FILE_TYPE, &
           DETERMINE_AVHRR_1, &
           DEFINE_1B_DATA, &
           CALCULATE_ASC_DES

  private::  &
           LAGRANGIAN_ANCHOR_INTERP, &
           LINEAR_ANCHOR_INTERP, &
           GNOMIC_ANCHOR_INTERP, &
           COMPUTE_ANGLE_ANCHORS, &
           i4word_to_string, &
           REF_CAL,  &
           THERM_CAL, &
           COMPUTE_NEW_THERM_CAL_COEF,  &
           MAKE_I4WORD, &
           MAKE_I2WORD, &
           UNPACK_AVHRR_HEADER_RECORD, &
           UNPACK_AVHRR_HEADER_RECORD_KLM, &
           UNPACK_AVHRR_DATA_RECORD,  &
           UNPACK_AVHRR_DATA_RECORD_KLM, &
           UNPACK_AVHRR_DATA_RECORD_AAPP, &
           CREATE_AVHRR_ARRAYS, &
           RESET_AVHRR_ARRAYS, &
           DESTROY_AVHRR_ARRAYS, &
           CONVERT_AVHRR_COUNTS_SINGLE_GAIN

  !--- variable declaration
  integer(kind=int2), dimension(5,10), private:: Space_Count_Temp
  integer(kind=int2), dimension(3:5,10),private:: Blackbody_Count
  integer(kind=int2), dimension(103),private:: Telemetry_Word
  integer(kind=int1), public, save:: Prt_Idx=0, Cal_Idx = 0
  integer(kind=int4), parameter, private:: Ncal = 10 
  integer(kind=int2), private, save:: sum_Prt
  real(kind=real4), private:: Mean_Prt
  real(kind=real4), private:: NBB_3b
  real(kind=real4), private:: NBB_4
  real(kind=real4), private:: NBB_5
  real(kind=real4), dimension(Ncal), private:: T_Prt
  logical, public, save:: Valid_Therm_Cal
  logical, public, save:: Spinup_New_Therm_Cal
  logical, public, save:: Valid_Ref_Cal
  real(kind=real4), dimension(Ncal),private,save:: Space_Count_Temp_1
  real(kind=real4), dimension(Ncal),private,save:: Space_Count_Temp_2
  real(kind=real4), dimension(Ncal),private,save:: Space_Count_Temp_3
  real(kind=real4), dimension(Ncal),private,save:: Space_Count_Temp_4
  real(kind=real4), dimension(Ncal),private,save:: Space_Count_Temp_5
  real(kind=real4), dimension(Ncal),private,save:: BB_Count_3
  real(kind=real4), dimension(Ncal),private,save:: BB_Count_4
  real(kind=real4), dimension(Ncal),private,save:: BB_Count_5
  real(kind=real4),public,save::Mean_Space_Count_1=0.0
  real(kind=real4),public,save::Mean_Space_Count_2=0.0
  real(kind=real4),public,save::Mean_Space_Count_3=0.0
  real(kind=real4),public,save::Mean_Space_Count_4=0.0
  real(kind=real4),public,save::Mean_Space_Count_5=0.0
  real(kind=real4),public,save::Mean_T_Prt=0.0
  real(kind=real4),public,save::Mean_BB_Count_3=0.0
  real(kind=real4),public,save::Mean_BB_Count_4=0.0
  real(kind=real4),public,save::Mean_BB_Count_5=0.0
  real(kind=real4), private, parameter:: Cal_Smoothing_Weight=0.2
  integer,public, save:: Calc_Asc_Des

  character(24), parameter, private :: MOD_PROMPT = " AVHRR_MODULE: "

  !--- values of iostatus indicative of an eof.  These values are not standard
  !--- on all fortran compilers.  These values come from xlf90, lahey, intel,
  !--- g95 and gfortran - there must be a more elegant solution
  integer, parameter, private:: Num_eof_ios = 6
  integer, dimension(Num_eof_ios), parameter, private::  &
           eof_iostatus = (/-1,-1,36,213,80,5002/)


  integer(kind=int4), parameter:: Biggest_I4 = HUGE(Biggest_I4)
  integer(kind=int2), parameter:: Biggest_I2 = HUGE(Biggest_I2)

  real(kind=real4), private:: AVHRR_ALTITUDE = 870.0   !954.0 !km
  real(kind=real4), private:: AVHRR_ALTITUDE_MORNING = 833.0 !km
  real(kind=real4), private:: AVHRR_ALTITUDE_AFTERNOON = 870.0 !km
  real(kind=real4), private:: AVHRR_ALTITUDE_METOP = 817.0 !km

  !---------------------------------------------------------------------
  ! arrays and variables formerly in pixel common)
  !---------------------------------------------------------------------
  integer (kind=int1), dimension(:), allocatable, private, save:: Header_Buffer_AVHRR
  integer (kind=int1), dimension(:), allocatable, private, save:: Buffer_AVHRR
  integer (kind=int1), dimension(:), allocatable, private, save:: Segment_Buffer_AVHRR
  integer (kind=int1), dimension(:), allocatable, private, save:: Fatal_AVHRR
  integer (kind=int1), dimension(:), allocatable, private, save:: Resync_AVHRR
  integer (kind=int1), dimension(:), allocatable, private, save:: Clavr_Status_1b_AVHRR  
  integer (kind=int2), dimension(:,:,:), allocatable, private, save:: Chan_Counts_Avhrr
  real(kind=int4), dimension(:,:), allocatable, private, save:: Scan_Space_Counts_Avhrr
  real(kind=int4), dimension(:,:), allocatable, private, save:: Scan_BB_Counts_Avhrr
  integer (kind=int2), dimension(:,:,:), allocatable, private, save:: Chan_Counts_Avhrr_Sg
  real (kind=real4), dimension(:,:), allocatable, private, save:: IR_Coef_1_1b,IR_Coef_2_1b,IR_Coef_3_1b
  real (kind=real4), dimension(:,:), allocatable, private, save:: IR_Linear_Slope_1b
  real (kind=real4), dimension(:,:), allocatable, private, save:: IR_Linear_Intercept_1b
  real (kind=real4), dimension(:,:), allocatable, private, save:: IR_Linear_Slope_New
  real (kind=real4), dimension(:,:), allocatable, private, save:: IR_Linear_Intercept_New
  real (kind=real4), dimension(:,:), allocatable, private, save:: Vis_Slope_1
  real (kind=real4), dimension(:,:), allocatable, private, save:: Vis_Slope_2
  real (kind=real4), dimension(:,:), allocatable, private, save:: Vis_Intercept_1
  real (kind=real4), dimension(:,:), allocatable, private, save:: Vis_Intercept_2
  real (kind=real4), dimension(:,:), allocatable, private, save:: Vis_Intersection

  contains
!-------------------------------------------------------------------
!-- define a Sc_Id that is unique value for each satellite
!-------------------------------------------------------------------
 subroutine ASSIGN_AVHRR_SAT_ID_NUM_INTERNAL(Sc_Id,Sc_Id_Internal,Sc_Id_Char)

  integer(kind=int2), intent(in):: Sc_Id
  integer(kind=int4), intent(out):: Sc_Id_Internal
  character(len=6), intent(out):: Sc_Id_Char

  Sc_Id_Internal = 0
  if (AVHRR_KLM_Flag == sym%NO) then
    if(Sc_Id == 1 .and. avhrr_1_Flag == sym%YES) then  !TIROS-N
       Sc_Id_Internal = 5
       Sc_Id_WMO = 708
       Sc_Id_Char = 'tirosn'
       Platform_Name_Attribute = 'TIROS-N'
       Sensor_Name_Attribute = 'AVHRR/1'
!      Sensor_Name_Attribute = 'TIROS-N : AVHRR/1'
    endif
    if(Sc_Id == 2 .and. avhrr_1_Flag == sym%YES) then  !NOAA-6
       Sc_Id_Internal = 6
       Sc_Id_WMO = 706
       Sc_Id_Char = 'noaa06'
       Platform_Name_Attribute = 'NOAA-6'
       Sensor_Name_Attribute = 'AVHRR/1'
!      Sensor_Name_Attribute = 'NOAA-6 : AVHRR/1'
    endif
    if(Sc_Id == 4 .and. AVHRR_KLM_Flag == sym%NO) then  !NOAA-7
       Sc_Id_Internal = 7
       Sc_Id_WMO = 707
       Sc_Id_Char = 'noaa07'
       Platform_Name_Attribute = 'NOAA-7'
       Sensor_Name_Attribute = 'AVHRR/2'
!      Sensor_Name_Attribute = 'NOAA-7 : AVHRR/2'
    endif
    if(Sc_Id == 6 .and. avhrr_1_Flag == sym%YES) then  !NOAA-8
       Sc_Id_Internal = 8
       Sc_Id_WMO = 200
       Sc_Id_Char = 'noaa08'
       Platform_Name_Attribute = 'NOAA-8'
       Sensor_Name_Attribute = 'AVHRR/1'
!      Sensor_Name_Attribute = 'NOAA-8 : AVHRR/1'
    endif
    if(Sc_Id == 7 .and. avhrr_1_Flag == sym%NO) then  !NOAA-9
       Sc_Id_Internal = 9
       Sc_Id_WMO = 201
       Sc_Id_Char = 'noaa09'
       Platform_Name_Attribute = 'NOAA-9'
       Sensor_Name_Attribute = 'AVHRR/2'
!      Sensor_Name_Attribute = 'NOAA-9 : AVHRR/2'
    endif
    if(Sc_Id == 8 .and. avhrr_1_Flag == sym%YES) then  !NOAA-10
       Sc_Id_Internal = 10
       Sc_Id_WMO = 202
       Sc_Id_Char = 'noaa10'
       Platform_Name_Attribute = 'NOAA-10'
       Sensor_Name_Attribute = 'AVHRR/1'
!      Sensor_Name_Attribute = 'NOAA-10 : AVHRR/1'
    endif
    if(Sc_Id == 1 .and. avhrr_1_Flag == sym%NO) then  !NOAA-11
       Sc_Id_Internal = 11
       Sc_Id_WMO = 203
       Sc_Id_Char = 'noaa11'
       Platform_Name_Attribute = 'NOAA-11'
       Sensor_Name_Attribute = 'AVHRR/2'
!      Sensor_Name_Attribute = 'NOAA-11 : AVHRR/2'
    endif
    if(Sc_Id == 5 .and. AVHRR_KLM_Flag==sym%NO) then             !NOAA-12
       Sc_Id_Internal = 12
       Sc_Id_WMO = 204
       Sc_Id_Char = 'noaa12'
       Platform_Name_Attribute = 'NOAA-12'
       Sensor_Name_Attribute = 'AVHRR/2'
!      Sensor_Name_Attribute = 'NOAA-12 : AVHRR/2'
    endif
    if(Sc_Id == 3 .and. AVHRR_KLM_Flag == sym%NO) then           !NOAA-14
       Sc_Id_Internal = 14
       Sc_Id_WMO = 205
       Sc_Id_Char = 'noaa14'
       Platform_Name_Attribute = 'NOAA-14'
       Sensor_Name_Attribute = 'AVHRR/2'
!      Sensor_Name_Attribute = 'NOAA-14 : AVHRR/2'
    endif
  else
    if(Sc_Id == 4 .and. AVHRR_KLM_Flag == sym%YES) then          !NOAA-15
       Sc_Id_Internal = 15
       Sc_Id_WMO = 206
       Sc_Id_Char = 'noaa15'
       Platform_Name_Attribute = 'NOAA-15'
       Sensor_Name_Attribute = 'AVHRR/3'
!      Sensor_Name_Attribute = 'NOAA-15 : AVHRR/3'
    endif
    if(Sc_Id == 2 .and. AVHRR_KLM_Flag == sym%YES) then          !NOAA-16
       Sc_Id_Internal = 16
       Sc_Id_WMO = 207
       Sc_Id_Char = 'noaa16'
       Platform_Name_Attribute = 'NOAA-16'
       Sensor_Name_Attribute = 'AVHRR/3'
!      Sensor_Name_Attribute = 'NOAA-16 : AVHRR/3'
    endif
    if(Sc_Id == 6 .and. AVHRR_KLM_Flag == sym%YES) then          !NOAA-17
       Sc_Id_Internal = 17
       Sc_Id_WMO = 208
       Sc_Id_Char = 'noaa17'
       Platform_Name_Attribute = 'NOAA-17'
       Sensor_Name_Attribute = 'AVHRR/3'
!      Sensor_Name_Attribute = 'NOAA-17 : AVHRR/3'
    endif
    if(Sc_Id == 7 .and. AVHRR_KLM_Flag == sym%YES) then          !NOAA-18
       Sc_Id_Internal = 18
       Sc_Id_WMO = 209
       Sc_Id_Char = 'noaa18'
       Platform_Name_Attribute = 'NOAA-18'
       Sensor_Name_Attribute = 'AVHRR/3'
!      Sensor_Name_Attribute = 'NOAA-18 : AVHRR/3'
    endif
    if(Sc_Id == 8 .and. AVHRR_KLM_Flag == sym%YES) then         !NOAA-19
       Sc_Id_Internal = 19
       Sc_Id_WMO = 223
       Sc_Id_Char = 'noaa19'
       Platform_Name_Attribute = 'NOAA-19'
       Sensor_Name_Attribute = 'AVHRR/3'
!      Sensor_Name_Attribute = 'NOAA-19 : AVHRR/3'
    endif
    if(Sc_Id == 12 .and. AVHRR_KLM_Flag == sym%YES) then         !Metop-A
       Sc_Id_Internal =  2
       Sc_Id_WMO = 4
       Sc_Id_Char = 'metopa'
       Platform_Name_Attribute = 'METOP-A'
       Sensor_Name_Attribute = 'AVHRR/3'
!      Sensor_Name_Attribute = 'METOP-A : AVHRR/3'
    endif

    if(Sc_Id == 11 .and. AVHRR_KLM_Flag == sym%YES) then         !Metop-B
       Sc_Id_Internal =  1
       Sc_Id_WMO = 3
       Sc_Id_Char = 'metopb'
       Platform_Name_Attribute = 'METOP-B'
       Sensor_Name_Attribute = 'AVHRR/3'
!      Sensor_Name_Attribute = 'METOP-B : AVHRR/3'
    endif
    !----
    !------- Metop-C Sc_Id numbers are not unknown at this time
    !----
    if(Sc_Id == 13 .and. AVHRR_KLM_Flag == sym%YES) then         !Metop-C
       Sc_Id_Internal =  4
       Sc_Id_WMO = 5
       Sc_Id_Char = 'metopc'
       Platform_Name_Attribute = 'METOP-C'
       Sensor_Name_Attribute = 'AVHRR/3'
!      Sensor_Name_Attribute = 'METOP-C : AVHRR/3'
    endif

  endif

  end subroutine ASSIGN_AVHRR_SAT_ID_NUM_INTERNAL

!---------------------------------------------------------------
! read the avhrr constants into memory
!-----------------------------------------------------------------
subroutine READ_AVHRR_INSTR_CONSTANTS(Instr_Const_file)
 character(len=*), intent(in):: Instr_Const_file
 integer:: ios0, erstat
 integer:: Instr_Const_lun

 Instr_Const_lun = GET_LUN()

 open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

 erstat = 0
 if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, "Error opening AVHRR constants file, ios0 = ", ios0
    print*, 'file : ',trim(Instr_Const_file)
    stop 19
 endif

  read(unit=Instr_Const_lun,fmt="(a3)") sat_name
  read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
  read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
  read(unit=Instr_Const_lun,fmt=*) Space_Rad_3b, b0_3b, b1_3b, b2_3b
  read(unit=Instr_Const_lun,fmt=*) Space_Rad_4, b0_4, b1_4, b2_4
  read(unit=Instr_Const_lun,fmt=*) Space_Rad_5, b0_5, b1_5, b2_5
  read(unit=Instr_Const_lun,fmt=*) nu_20
  read(unit=Instr_Const_lun,fmt=*) a1_20
  read(unit=Instr_Const_lun,fmt=*) a2_20
  read(unit=Instr_Const_lun,fmt=*) nu_31
  read(unit=Instr_Const_lun,fmt=*) a1_31
  read(unit=Instr_Const_lun,fmt=*) a2_31
  read(unit=Instr_Const_lun,fmt=*) nu_32
  read(unit=Instr_Const_lun,fmt=*) a1_32
  read(unit=Instr_Const_lun,fmt=*) a2_32
  read(unit=Instr_Const_lun,fmt=*) Ch1_Dark_Count
  read(unit=Instr_Const_lun,fmt=*) Ch2_Dark_Count
  read(unit=Instr_Const_lun,fmt=*) Ch3a_Dark_Count
  read(unit=Instr_Const_lun,fmt=*) Ch1_Gain_Low_0,Ch1_Degrad_Low_1, Ch1_Degrad_Low_2
  read(unit=Instr_Const_lun,fmt=*) Ch1_Gain_High_0,Ch1_Degrad_High_1, Ch1_Degrad_High_2
  read(unit=Instr_Const_lun,fmt=*) Ch2_Gain_Low_0,Ch2_Degrad_Low_1, Ch2_Degrad_Low_2
  read(unit=Instr_Const_lun,fmt=*) Ch2_Gain_High_0,Ch2_Degrad_High_1, Ch2_Degrad_High_2
  read(unit=Instr_Const_lun,fmt=*) Ch3a_Gain_Low_0,Ch3a_Degrad_Low_1, Ch3a_Degrad_Low_2
  read(unit=Instr_Const_lun,fmt=*) Ch3a_Gain_High_0,Ch3a_Degrad_High_1, Ch3a_Degrad_High_2
  read(unit=Instr_Const_lun,fmt=*) launch_date
  read(unit=Instr_Const_lun,fmt=*) Ch1_Switch_Count
  read(unit=Instr_Const_lun,fmt=*) Ch2_Switch_Count
  read(unit=Instr_Const_lun,fmt=*) Ch3a_Switch_Count
  read(unit=Instr_Const_lun,fmt=*) Prt_Coef(:,1)
  read(unit=Instr_Const_lun,fmt=*) Prt_Coef(:,2)
  read(unit=Instr_Const_lun,fmt=*) Prt_Coef(:,3)
  read(unit=Instr_Const_lun,fmt=*) Prt_Coef(:,4)
  read(unit=Instr_Const_lun,fmt=*) prt_weight(:)
  read(unit=Instr_Const_lun,fmt=*) b1_day_mask,b2_day_mask,b3_day_mask,b4_day_mask
  close(unit=Instr_Const_lun)

!-- convert solar flux in channel 3b to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

!----- define planck constants here

end subroutine READ_AVHRR_INSTR_CONSTANTS


!----------------------------------------------------------------------
! DEFINE_1B_DATA - based on the AVHRR_GAC_Flag and AVHRR_KLM_Flag flags, this routine defines
!                  the parameters necessary to read in the header and
!                  data records from the various 1B data types. This allows
!                  defines quantities to allow the clavr cloud mask to be
!                  written into the AVHRR_KLM_Flag 1b files.
!
! input
!      AVHRR_GAC_Flag - T if AVHRR_GAC_Flag, F if hrpt or lac
!      AVHRR_KLM_Flag - T if NOAA-AVHRR_KLM_Flag, F if NOAA-J or earlier
!      AVHRR_AAPP_Flag - T if AAPP, F if NESDIS
! output
!      Num_Pix - the number of pixels per scanline
!      l1b_rec_length - the number of bytes per record of the 1b file
!      Nrec_Header - the number of records used for the header
!      Nword_Clavr_Start - byte number of first clavr byte on 1b
!      Nword_Clavr - number of bytes used for clavr cloud mask on 1b
!-------------------------------------------------------------------------------
 subroutine DEFINE_1B_DATA(AVHRR_GAC_Flag,AVHRR_KLM_Flag,AVHRR_AAPP_Flag,Nrec_Header,Nword_Clavr_Start,Nword_Clavr)

 integer(kind=int4), intent(in):: AVHRR_GAC_Flag
 integer(kind=int4), intent(in):: AVHRR_KLM_Flag
 integer(kind=int4), intent(in):: AVHRR_AAPP_Flag
 integer(kind=int4), intent(out):: Nword_Clavr
 integer(kind=int4), intent(out):: Nword_Clavr_Start
 integer(kind=int4), intent(out):: Nrec_Header

  if (AVHRR_GAC_Flag == sym%YES) then
      Num_Pix = 409              !number of pixels per scan line
      if (AVHRR_KLM_Flag == sym%YES) then
       l1b_rec_length  = 4608
       Nrec_Header = 1
       Nword_Clavr_Start = 4049
       Nword_Clavr = 2*52
      else
       l1b_rec_length  = 3220
       Nrec_Header = 2
       Nword_Clavr_Start = 0
       Nword_Clavr = 0
      endif
   else
     if (AVHRR_KLM_Flag == sym%YES) then
       Num_Pix = 2048
       Nrec_Header = 1
       l1b_rec_length = 15872
       Nword_Clavr_Start = 14977
       Nword_Clavr = 2*256
     else
      Num_Pix = 2048
      Nrec_Header = 1
      l1b_rec_length = 14800
      Nword_Clavr_Start = 0
      Nword_Clavr = 0
      endif
   endif

!---- AVHRR_AAPP_Flag - assumed to be 1 km  and no clavr bits
   if (AVHRR_AAPP_Flag == sym%YES) then
       Num_Pix = 2048
       Nrec_Header = 1
       l1b_rec_length = 22016
       Nword_Clavr_Start = 0
       Nword_Clavr = 0
    endif

!--- number of geolocation Anchors
  Num_Anchors = 51


 end subroutine DEFINE_1B_DATA
!======================================================================
! Determine what type of file based on the header
!======================================================================
!-------------------------------------------------------
! read header to get data type and version of 1b data
! inorder to get AVHRR_KLM_Flag and AVHRR_GAC_Flag which determine
! the data format for subsequent read statements
! note, data_type and ver_1b are read in again later
!-------------------------------------------------------
   subroutine DETERMINE_AVHRR_FILE_TYPE(file_1b_local,AVHRR_GAC_Flag,AVHRR_KLM_Flag,AVHRR_AAPP_Flag,ver_1b,data_type,Byte_Swap_1b)

    character(len=*), intent(in):: file_1b_local
    integer(kind=int4), intent(out):: AVHRR_GAC_Flag
    integer(kind=int4), intent(out):: AVHRR_KLM_Flag
    integer(kind=int4), intent(out):: AVHRR_AAPP_Flag
    integer(kind=int2), intent(out):: ver_1b
    integer(kind=int2), intent(out):: data_type
    integer(kind=int4), intent(out):: Byte_Swap_1b

    integer(kind=int1), dimension(100):: Header_Buffer_Temp
    integer:: word_start
    integer:: bytes_per_word
    integer:: Number_Of_Words
    integer:: Number_Of_Words_Read

    !--- use c-routine
    bytes_per_word = 1
    word_start = 0
    Number_Of_Words = 100
    call mreadf_int(file_1b_local//CHAR(0),word_start,bytes_per_word, &
                    Number_Of_Words,Number_Of_Words_Read,Header_Buffer_Temp)

    !--- use instrinsic fortran read
!   word_start = 1
!   read(unit=lun_level1b,pos=word_start) Header_Buffer_Temp


    !--- look at first byte and determine AVHRR_KLM_Flag flag - A. Jelenak Suggestion  - check this for AAPP!
     if (Header_Buffer_Temp(1) <= 10) then
        AVHRR_KLM_Flag = sym%NO
     elseif (Header_Buffer_Temp(1) >= 65) then
        AVHRR_KLM_Flag = sym%YES
     else
        print *, EXE_PROMPT, MOD_PROMPT, "Error diagnosing AVHRR_KLM_Flag from first byte, assuming AVHRR_KLM_Flag = 0"
        AVHRR_KLM_Flag = sym%NO
     endif

!--- assuming AVHRR_KLM_Flag, unpack ver1b, Valid values are (1-5)
     ver_1b = MAKE_I2WORD(Header_Buffer_Temp(5:6),sym%UNSIGNED,sym%NOSWAP)

!--- if version number unrealistic, assume AVHRR_AAPP_Flag
     AVHRR_AAPP_Flag = sym%NO
     Byte_Swap_1b = sym%NO
     if (AVHRR_KLM_Flag == sym%YES) then
        if ((ver_1b > 10).or.(ver_1b < 1)) then
          AVHRR_AAPP_Flag = sym%YES
          Byte_Swap_1b = sym%YES
          ver_1b = MAKE_I2WORD(Header_Buffer_Temp(5:6),sym%UNSIGNED,Byte_Swap_1b)
        endif
     endif

!--- read spacedraft id
     Sc_Id = MAKE_I2WORD(Header_Buffer_Temp(73:74),sym%UNSIGNED,Byte_Swap_1b)

!--- fix noaa-15 bug
     if ((ver_1b == 1) .and. (Sc_Id == 4)) then
       ver_1b = 2
     endif

!--- assign pre-AVHRR_KLM_Flag to version 1
     if (AVHRR_KLM_Flag == sym%NO) then
       ver_1b = 1
     endif

!--- based on AVHRR_KLM_Flag flag, choose bytes to construct data_type
     if (AVHRR_KLM_Flag == sym%YES) then
         data_type = MAKE_I2WORD(Header_Buffer_Temp(77:78),sym%UNSIGNED,Byte_Swap_1b)
     else
         data_type = ishft(Header_Buffer_Temp(2),-4)
     endif

!--- based on data_type, set AVHRR_GAC_Flag
     AVHRR_GAC_Flag = sym%NO
     if (data_type == 2) then
       AVHRR_GAC_Flag = sym%YES
     endif

end subroutine DETERMINE_AVHRR_FILE_TYPE

!==================================================================================================
!--- determine if this data is from AVHRR/1 which has no channel 5
!==================================================================================================
subroutine DETERMINE_AVHRR_1(year, AVHRR_KLM_Flag_Flag, avhrr_1_Flag)

    integer(kind=int2), intent(in):: year       !year of this data-set
    integer(kind=int4), intent(in):: AVHRR_KLM_Flag_Flag   !AVHRR_KLM_Flag flag (yes/no)
    integer(kind=int4), intent(out):: avhrr_1_Flag  !AVHRR/1 flag (yes/no)

     !--- initialize to no
     avhrr_1_Flag = sym%NO

     !--- check Sc_Id for pre-AVHRR_KLM_Flag data
     if (AVHRR_KLM_Flag_Flag == sym%NO) then
      if (Sc_Id == 1 .and. year < 1985) avhrr_1_Flag = sym%YES   !TIROS-N (NOAA-11 repeat) 
      if (Sc_Id == 2) avhrr_1_Flag = sym%YES   !NOAA-6
      if (Sc_Id == 6) avhrr_1_Flag = sym%YES   !NOAA-8
      if (Sc_Id == 8) avhrr_1_Flag = sym%YES   !NOAA-10
     endif

end subroutine DETERMINE_AVHRR_1

!----------------------------------------------------------------------------
! This routine reads in a specifed number of scans from a level 1b file
! it calibrates, navigates and quality controls the data
!
! this processes the level1b data two scanline lines at a time
!
!----------------------------------------------------------------------------
  subroutine READ_AVHRR_LEVEL1B_DATA(file_1b_local, & 
                                     AVHRR_KLM_Flag, &
                                     AVHRR_AAPP_Flag, &
                                     Therm_Cal_1b, &
                                     Time_Since_Launch,  &
                                     Nrec_Header, &
                                     Seg_Idx)

    character(len=*), intent(in):: file_1b_local
    integer, intent(in):: AVHRR_KLM_Flag
    integer, intent(in):: AVHRR_AAPP_Flag
    integer, intent(in):: Therm_Cal_1b
    real(kind=real4), intent(in):: Time_Since_Launch
    integer, intent(in):: Nrec_Header
    integer, intent(in):: Seg_Idx
    integer:: Line_Idx
    integer:: Irec_start
    integer:: Irec_end
    integer:: Irec
    integer:: Error_Code

    integer(kind=int4):: word_start
    integer(kind=int4):: word_end
    integer(kind=int4):: Number_Of_Words
    integer(kind=int4):: Number_Of_Words_Read
    integer(kind=int4):: bytes_per_word

    !--- initialize Calc_Asc_Des flag
    Calc_Asc_Des = sym%NO

    !--- based on segment number, compute starting and ending record in level-1b
    !--- for this segment
    Irec_start = Nrec_Header + (Seg_Idx-1)*Num_Scans_per_Segment + 1
    Irec_end =   Irec_start + Num_Scans_per_Segment - 1

    !--- initialize number of scans-read index
    Num_Scans_Read = 0

    !--- create interal arrays needed to process this segment
    call CREATE_AVHRR_ARRAYS(Num_Pix, Num_Scans_Per_Segment, L1b_Rec_Length)

    !--- reset interal arrays needed to process this segment
    call RESET_AVHRR_ARRAYS()

    !-------------------------------------------------------------------
    ! read one segment of data from level1b file
    !-------------------------------------------------------------------
    Irec = Irec_start
    Line_Idx = Line_Idx_Min_Segment

    !------
    !-- try to read in segment as a whole
    !----

    !--- initialize
    Segment_Buffer_AVHRR = 0


    !--- use external c-routine for read
    bytes_per_word = 1
    word_start = l1b_rec_length* ( 1 + (Seg_Idx-1)*Num_Scans_per_Segment)
    Number_Of_Words = l1b_rec_length*Num_Scans_per_Segment
    call mreadf_int(file_1b_local//CHAR(0),word_start,bytes_per_word, &
                   Number_Of_Words,Number_Of_Words_Read,Segment_Buffer_AVHRR)
    !--- update number of scans read
    Num_Scans_Read = Number_Of_Words_Read / l1b_rec_length

!   !--- use instrinsic fortran read
!   word_start = l1b_rec_length* ( 1 + (Seg_Idx-1)*Num_Scans_per_Segment) + 1
!   read(unit=lun_level1b,pos=word_start,iostat=ios_l1b) Segment_Buffer_AVHRR
!   if (ios_l1b == 0) then
!      Num_Scans_Read = Num_Scans_per_Segment
!   else
!      print *, "error on read from l1b", ios_l1b
!   stop
!   endif
 
  

   !----------------------------------------------------------------------
   !  reset thermal calibration stats for this orbit
   !----------------------------------------------------------------------
   if (Seg_Idx == 1) then
          Valid_Ref_cal = .false.
          Valid_therm_cal = .false.
          Spinup_New_therm_cal = .true.
          Prt_Idx = 0
          Cal_Idx = 0
          Mean_T_Prt = 0.0
          Mean_Space_Count_1 = 0.0
          Mean_Space_Count_2 = 0.0
          Mean_Space_Count_3 = 0.0
          Mean_Space_Count_4 = 0.0
          Mean_Space_Count_5 = 0.0
          Mean_BB_Count_3 = 0.0
          Mean_BB_Count_4 = 0.0
          Mean_BB_Count_5 = 0.0
          Ch1_Gain_Low = 0.0
          Ch1_Gain_High = 0.0
          Ch2_Gain_Low = 0.0
          Ch2_Gain_High = 0.0
          Ch3a_Gain_Low = 0.0
          Ch3a_Gain_High = 0.0
   endif

   !------------------------------------------------------------------------------------------------------------
   !--- UNPACK Level-1b Data
   !------------------------------------------------------------------------------------------------------------
   segment_loop: do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1

       !extract level-1b information for this scanline
       word_start = (Line_Idx-1)*l1b_rec_length + 1
       word_end = word_start +  l1b_rec_length - 1

       Buffer_AVHRR = Segment_Buffer_AVHRR(word_start:word_end)

       !--- unpack the level-1b data
       if (AVHRR_AAPP_Flag == sym%YES) then
           call UNPACK_AVHRR_DATA_RECORD_AAPP(Line_Idx)
       else
         if (AVHRR_KLM_Flag == sym%YES) then
           call UNPACK_AVHRR_DATA_RECORD_KLM(Line_Idx)
         else
           call UNPACK_AVHRR_DATA_RECORD(Line_Idx)
           Ch3a_On_AVHRR(Line_Idx) = 0   !set for pre-AVHRR_KLM_Flag data
         endif
       endif


       !----------------------------------------------------------------------
       ! generate new calibration coefficients
       !----------------------------------------------------------------------
       if (Therm_Cal_1b == sym%NO) then

          call COMPUTE_NEW_THERM_CAL_COEF(Line_Idx)

          !--- if sufficient scans lines have not been processed, this is fatal
          if (Spinup_New_Therm_Cal .eqv. .true.) then
            Fatal_AVHRR(Line_Idx) = sym%YES
          endif
        else
           Valid_Therm_Cal = .true.
       endif

       !--- compute reflectance calibration
       call REF_CAL(Time_Since_Launch,Line_Idx)

       !---- determine if Valid coefficients are made
       if (Valid_Therm_Cal .eqv. .false.) then
         Fatal_AVHRR(Line_Idx) = sym%YES
       endif
       if (Valid_Ref_Cal .eqv. .false.) then
         Fatal_AVHRR(Line_Idx) = sym%YES
       endif

      !------------------------------------------------------------------------------
      ! Compute latitude and longitude
      !------------------------------------------------------------------------------
      if (maxval(abs(Lat_Anchor_1b(:,Line_Idx))) > 85.0) then
         call GNOMIC_ANCHOR_INTERP(Lon_Anchor_1b(:,Line_Idx),Lat_Anchor_1b(:,Line_Idx), &
                                   Lon_1b(:,Line_Idx),Lat_1b(:,Line_Idx))
      else
         call LAGRANGIAN_ANCHOR_INTERP(Lat_Anchor_1b(:,Line_Idx),Lat_1b(:,Line_Idx))
         call LAGRANGIAN_ANCHOR_INTERP(Lon_Anchor_1b(:,Line_Idx),Lon_1b(:,Line_Idx))
      endif

      !--- constrain lat and lon to physical values
      if ((minval(Lat_1b(:,Line_Idx)) < -90.0).or.(maxval(Lat_1b(:,Line_Idx)) > 90.0)) then
           Lat_1b(:,Line_Idx) = min(90.0,max(-90.0,Lat_1b(:,Line_Idx)))
      endif
      if ((minval(Lon_1b(:,Line_Idx)) < -180.0).or.(maxval(Lon_1b(:,Line_Idx)) > 180.0)) then
           Lon_1b(:,Line_Idx) = min(180.0,max(-180.0,Lon_1b(:,Line_Idx)))
      endif 

      !----------------------------------------------------------------
      !--- Compute viewing geometry
      !----------------------------------------------------------------

!     !--- compute Anchor points for pre-AVHRR_KLM_Flag series (not included in level-1b)
!     if (AVHRR_KLM_Flag == sym%NO) then
        call COMPUTE_ANGLE_ANCHORS(Scan_Day(Line_Idx),        &
                                   Scan_Time(Line_Idx),       &
                                   Lon_Anchor_1b(:,Line_Idx), &
                                   Lat_Anchor_1b(:,Line_Idx), &
                                   Satzen_Anchor(:,Line_Idx), &
                                   Solzen_Anchor(:,Line_Idx), &
                                   Relaz_Anchor(:,Line_Idx),  &
                                   Solaz_Anchor(:,Line_Idx),  &
                                   Sataz_Anchor(:,Line_Idx),  &
                                   Glintzen_Anchor(:,Line_Idx),  &
                                   Scatangle_Anchor(:,Line_Idx))
!     else
!       !--- modify AVHRR_KLM_Flag Relaz to expected convention
!       Relaz_Anchor(:,Line_Idx) = 180.0 - abs(Relaz_Anchor(:,Line_Idx))
!     endif

      !--- interpolate pixel values from Anchor values
      call LINEAR_ANCHOR_INTERP(Satzen_Anchor(:,Line_Idx),Satzen(:,Line_Idx))
      call LINEAR_ANCHOR_INTERP(Solzen_Anchor(:,Line_Idx),Solzen(:,Line_Idx))
      call LINEAR_ANCHOR_INTERP(Relaz_Anchor(:,Line_Idx),Relaz(:,Line_Idx))
      call LINEAR_ANCHOR_INTERP(Solaz_Anchor(:,Line_Idx),Solaz(:,Line_Idx))
      call LINEAR_ANCHOR_INTERP(Sataz_Anchor(:,Line_Idx),Sataz(:,Line_Idx))
      call LINEAR_ANCHOR_INTERP(Glintzen_Anchor(:,Line_Idx),Glintzen(:,Line_Idx))
      call LINEAR_ANCHOR_INTERP(Scatangle_Anchor(:,Line_Idx),Scatangle(:,Line_Idx))

      !--- constrain sensor zenith
      if ((minval(Satzen(:,Line_Idx)) < 0.0).or.(maxval(Satzen(:,Line_Idx)) > 89.9)) then
          Satzen(:,Line_Idx) = min(89.9,max(0.0,Satzen(:,Line_Idx)))
      endif

      !------------------------------------------------------------------------------
      ! check for bad scan lines here
      !------------------------------------------------------------------------------
       Bad_Scan_Flag(Line_Idx) = sym%NO
       Error_Code = 0

       !--- any scan with a fatal error is considered bad
       if (Fatal_AVHRR(Line_Idx) == sym%YES) then
           Bad_Scan_Flag(Line_Idx) = sym%YES
           Error_Code = 1
       endif

       !--- if in spinup mode, set all scans as bad 
       !--- note, these are also set to have a fatal code, but this overwrites
       !--- it to limit printing to screen of excessive reports
       if ((Therm_Cal_1b == sym%NO) .and. (Spinup_New_Therm_Cal .eqv. .true.)) then
           Bad_Scan_Flag(Line_Idx) = sym%YES
           Error_Code = 2
       endif

       !---- check for missing navigation
       if ( (sum(Lat_Anchor_1b(:,Line_Idx)) == 0.00).or. (sum(Lon_Anchor_1b(:,Line_Idx)) == 0.00)) then
           Bad_Scan_Flag(Line_Idx) = sym%YES
           Error_Code = 3
       endif

       !--- check erroneous asc/des flag 
       !--- erroneous values are common on AVHRR/1 data
       !--- so this is no longer considered a fatal error
       if ((Ascend(Line_Idx) < 0) .or. (Ascend(Line_Idx)  > 1)) then
           Calc_Asc_Des = sym%YES
       endif

       !---- check angle Anchors on AVHRR_KLM_Flag data
       if (AVHRR_KLM_Flag == sym%YES) then
         if ((minval(Satzen_Anchor(:,Line_Idx)) < 0.0) .or. (maxval(Satzen_Anchor(:,Line_Idx)) > 90.0)) then
           Bad_Scan_Flag(Line_Idx) = sym%YES
           Error_Code  = 5
         endif
         if ((minval(Solzen_Anchor(:,Line_Idx)) < 0.0) .or. (maxval(Solzen_Anchor(:,Line_Idx)) > 180.0)) then
           Bad_Scan_Flag(Line_Idx) = sym%YES
           Error_Code = 5
         endif
         if ((minval(Relaz_Anchor(:,Line_Idx)) < -180.0) .or. (maxval(Relaz_Anchor(:,Line_Idx)) > 180.0)) then
           Bad_Scan_Flag(Line_Idx) = sym%YES
           Error_Code = 5
         endif
       endif

       !----- report bad scan code to standard output (but not bad scans during
       !----- spin-up of thermal cal
       if (Bad_Scan_Flag(Line_Idx) == sym%YES .and. Error_Code /= 2) then
         write(unit=6,fmt="(a,a,a,i6,a,i2)") EXE_PROMPT, MOD_PROMPT,  &
               "BAD SCAN: Numbers = ", Scan_Number(Line_Idx), " error code = ", Error_Code
       endif 

       !---- set geolocation to missing for bad-scans to avoid ancil-data interp
       if (Bad_Scan_Flag(Line_Idx) == sym%YES) then
               Lat_1b(:,Line_Idx) = Missing_Value_Real4
               Lon_1b(:,Line_Idx) = Missing_Value_Real4
       endif

  end do segment_loop

  !------------------------------------------------------------------------------
  ! Calculate ascending/descending flag if needed for this segment
  !------------------------------------------------------------------------------
  if (Calc_Asc_Des == sym%YES) then
      call CALCULATE_ASC_DES(Line_Idx_Min_Segment,Num_Scans_Read)
  endif

  !------------------------------------------------------------------------------
  ! Calibrate segment
  !------------------------------------------------------------------------------

  !--- RADIANCE CALIBRATION (3b,4,5)
  call THERM_CAL(Num_Scans_Read)


  !--- convert AVHRR counts to single gain and store in global count arrays
  call CONVERT_AVHRR_COUNTS_SINGLE_GAIN(AVHRR_KLM_Flag,Line_Idx_Min_Segment,Num_Scans_Read)

  call DESTROY_AVHRR_ARRAYS()
end subroutine READ_AVHRR_LEVEL1B_DATA


!========================================================================
! Read in the AVHRR Level1b header
!
! input and output passed through avhrr_Pixel_common
!
! handles AVHRR_GAC_Flag,lac, hrpt for pre and post AVHRR_KLM_Flag data - and AVHRR_AAPP_Flag
!========================================================================
 subroutine READ_AVHRR_LEVEL1B_HEADER(file_1b_local)

   character(len=*), intent(in):: file_1b_local
   integer(kind=int4):: word_start
   integer(kind=int4):: Number_Of_Words
   integer(kind=int4):: Number_Of_Words_Read
   integer(kind=int4):: bytes_per_word

    !--- allocate header Buffer - taken from pixel_common
    allocate(Header_Buffer_AVHRR(l1b_rec_length))
 
    !--- use external c
    bytes_per_word = 1
    word_start = 0
    Number_Of_Words = l1b_rec_length
    call mreadf_int(file_1b_local//CHAR(0),word_start,bytes_per_word, &
                    Number_Of_Words,Number_Of_Words_Read,Header_Buffer_AVHRR)

!   !--- use internal fortran
!   word_start = 1
!   read(unit=lun_level1b,pos=word_start), Header_Buffer_AVHRR

    if (AVHRR_GAC_Flag == sym%YES) then

     if (AVHRR_KLM_Flag == sym%YES) then

      call UNPACK_AVHRR_HEADER_RECORD_KLM(Sc_Id,data_type,Start_Year, &
                   Start_day,Start_time,Num_Scans,End_Year,End_day,End_time, &
                   tip_parity,aux_sync,ramp_auto_Cal,proc_block_Id,ver_1b)

     else

      call UNPACK_AVHRR_HEADER_RECORD(Sc_Id,data_type,Start_Year, &
                   Start_day,Start_time,Num_Scans,End_Year,End_day,End_time, &
                   tip_parity,aux_sync,ramp_auto_Cal,proc_block_Id,ver_1b)

      !--- pre AVHRR_KLM_Flag used a 2 digit year
      if (End_Year > 50) then
       End_Year = End_Year + 1900
      else
        End_Year = End_Year + 2000
      endif
      if (Start_Year > 50) then
        Start_Year = Start_Year + 1900
      else
        Start_Year = Start_Year + 2000
      endif
                                                                                                                                          
     endif

   else     !if not AVHRR_GAC_Flag

     if (AVHRR_KLM_Flag == sym%YES) then

      call UNPACK_AVHRR_HEADER_RECORD_KLM(Sc_Id,data_type,Start_Year, &
                  Start_day,Start_time,Num_Scans,End_Year,End_day,End_time, &
                  tip_parity,aux_sync,ramp_auto_Cal,proc_block_Id,ver_1b)

     else

      call UNPACK_AVHRR_HEADER_RECORD(Sc_Id,data_type,Start_Year, &
                Start_day,Start_time,Num_Scans,End_Year,End_day,End_time, &
                tip_parity,aux_sync,ramp_auto_Cal,proc_block_Id,ver_1b)

!--- pre AVHRR_KLM_Flag used a 2 digit year
      if (End_Year > 50) then
        End_Year = End_Year + 1900
      else
        End_Year = End_Year + 2000
      endif
      if (Start_Year > 50) then
        Start_Year = Start_Year + 1900
      else
        Start_Year = Start_Year + 2000
      endif
     endif
   endif


   deallocate(Header_Buffer_AVHRR)

end subroutine READ_AVHRR_LEVEL1B_HEADER

!---------------------------------------------------------------------
! when creating a 1bx file, this routine write the original level 1b
! header to it
!---------------------------------------------------------------------
 subroutine WRITE_HEADER_1BX(file_1b)

    character(len=*), intent(in):: file_1b
    character(len=120):: file_1bx
    integer:: ios0,erstat

    file_1bx = trim(file_1b)
    if (file_1bx(1:3) == 'NSS') then
       file_1bx(1:3) = 'TCL'
    endif
    print *,EXE_PROMPT, "PROCESS_AVHRR unit 11 out 1bx==>",dir_1bx
    print *,EXE_PROMPT, "creating file ",trim(file_1bx)

    open(unit=lun_level1b,file=trim(dir_1bx)//trim(file_1bx),  &
         status="replace",form="unformatted", &
         access="dIrect",recl=l1b_rec_length,action="write",iostat=ios0)

         if (ios0 /= 0) then
           erstat = 12
           print *, EXE_PROMPT, "Error opening 1Bx file, ios0 = ",ios0
           stop 12
         endif

      write(unit=lun_level1b,rec=1,iostat=ios0) Header_Buffer_AVHRR

      if (ios0 /= 0) then
       print *, EXE_PROMPT, "Error writing 1Bx file header, ios0 = ",ios0
       stop 13
      endif

end subroutine WRITE_HEADER_1BX
!*******************************************************************


!-----------------------------------------------------------
! REFLECTANCE CALIBRATION ROUTINE
!
! comments given for channel 1 only but apply to channel 2,3a
!
! input
!   Time_Temp_Since_Launch - the Time_Temp in year_Temps since the launch
!   Ref_Cal_1b - flag controlling use of 1b coefficients
!
!  output
!   none - only through common memory
!
! Parameters defined in 'avhrr_Pixel_common.f' module
!  Ch1_Gain_low - calibration slope or gain in low count region
!                  with degradation accounted for
!  Ch1_Gain_High - calibration slope or gain in high count region
!                  with degradation accounted for
!
! internal parameters
!  Ref_Ch1_Switch - the albedo at count = count1_Switch
!
! arrays used from common memory
!  Chan_Counts_Avhrr_Temp(1,:,:) - two scanlines of channel 1 counts
!
! arrays modified in common memory
!  Ref_Ch1_Temp - channel 1 albedo for scanline pair, not this is
!         not corrected for sun-earth distance or for solar zenith angle
!
!  variables used from common memory
!   Ch1_Gain_Low_0 - calibration slope or gain in low count region at launch
!   Ch1_Gain_High_0 - calibration slope or gain in high count region at launch
!   Ch1_Switch_Count - count separating low from high count region
!   Ch1_Dark_Count_Cal - the dark count used for calibration - may differ from the
!                         the true dark count for any given slope,intercept pair
!
! Discussion
!  This code computes the coefficients to compute the reflectance calibration 
!  do to the following (Method 1)
!   R_low = slope_low * (C - C_dark)
!   R_High = R_Switch + slope_High*(C - C_Switch)
! 
! The method in the POD guide is as follows (Method 2)
!   R_low = Slope_low * C + intercept_low
!   R_High = Slope_High * C + intercept_High
!
! If Ref_Cal_1b is true, then the coefficients read in from the 1b 
! file are converted into forms for Method 1 (done in REF_CAL()).
! This transformation is exact
!
! Revision History
! Coded January 2003 - A. Heidinger
!
! 2004-09-03 - A. Jelenak
!     Original REF_CAL sub split into two to allow computation of channel gain
!     coeffs and their storage in the OBS file. REF_CAL_COEFFS sub does this
!     task now, and REF_CAL just applies them to the Chan_Counts_Avhrr_Temp() array.
!
! 2004-12-04 - A. Heidinger
!     Added capability to use 1b reflectance coeffcients
!
! 2008-12-24 - A. Heidinger
!     Now operates on single scan-line
!
! 2014-01-04 - A. Heidinger
!     Reintegrated REF_CAL_COEFFS into REF_CAL for simplicity
!
!
! Note, if Therm_Cal_1b == sym%NO, this routine must be called
!       after Compute_New_Therm_Cal_Coeffs
!-----------------------------------------------------------
 subroutine REF_CAL(Time_Temp_Since_Launch,Line_Idx)

   real(kind=real4), intent(in):: Time_Temp_Since_Launch
   integer, intent(in):: Line_Idx

   Valid_Ref_Cal = .true.

   !--- compute coefficient if using level-1b values
   if (Ref_Cal_1b == sym%YES) then

    Ch1_Gain_low = vis_slope_1(1,Line_Idx)
    Ch1_Dark_Count_Cal = -1.0*vis_intercept_1(1,Line_Idx)/vis_slope_1(1,Line_Idx)
    Ch1_Gain_High = vis_slope_2(1,Line_Idx)
    Ch1_Switch_Count_Cal = vis_intersection(1,Line_Idx)

    Ch2_Gain_low = vis_slope_1(2,Line_Idx)
    Ch2_Dark_Count_Cal = -1.0*vis_intercept_1(2,Line_Idx)/vis_slope_1(2,Line_Idx)
    Ch2_Gain_High = vis_slope_2(2,Line_Idx)
    Ch2_Switch_Count_Cal = vis_intersection(2,Line_Idx)

    Ch3a_Gain_low = vis_slope_1(3,Line_Idx)
    Ch3a_Dark_Count_Cal = -1.0*vis_intercept_1(3,Line_Idx)/vis_slope_1(3,Line_Idx)
    Ch3a_Gain_High = vis_slope_2(3,Line_Idx)
    Ch3a_Switch_Count_Cal = vis_intersection(3,Line_Idx)

    !---- compute reflectances at switch count (define consistent with)
    Ref_Ch1_Switch = vis_slope_2(1,Line_Idx) * vis_intersection(1,Line_Idx) + vis_intercept_2(1,Line_Idx)
    Ref_Ch2_Switch = vis_slope_2(2,Line_Idx) * vis_intersection(2,Line_Idx) + vis_intercept_2(2,Line_Idx)
    Ref_Ch6_Switch = vis_slope_2(3,Line_Idx) * vis_intersection(3,Line_Idx) + vis_intercept_2(3,Line_Idx)

   else

     !--- adjust gains to account for degradation
     if (Therm_Cal_1b == sym%YES) then
       Ch1_Dark_Count_Cal = Ch1_Dark_Count
       Ch2_Dark_Count_Cal = Ch2_Dark_Count
       Ch3a_Dark_Count_Cal = Ch3a_Dark_Count
     else
        Ch1_Dark_Count_Cal =  Scan_Space_Counts_Avhrr(1,Line_Idx) 
        Ch2_Dark_Count_Cal = Scan_Space_Counts_Avhrr(2,Line_Idx)
        Ch3a_Dark_Count_Cal = Scan_Space_Counts_Avhrr(3,Line_Idx)
     endif

     !-- set switch counts to values in instrument files
     Ch1_Switch_Count_Cal = Ch1_Switch_Count
     Ch2_Switch_Count_Cal = Ch2_Switch_Count
     Ch3a_Switch_Count_Cal = Ch3a_Switch_Count

     Ch1_Gain_low = Ch1_Gain_Low_0*(100.0+Ch1_Degrad_Low_1*Time_Temp_Since_Launch+ &
                                  Ch1_Degrad_Low_2*Time_Temp_Since_Launch**2)/100.0
     Ch1_Gain_High = Ch1_Gain_High_0*(100.0+Ch1_Degrad_High_1*Time_Temp_Since_Launch+ &
                                  Ch1_Degrad_High_2*Time_Temp_Since_Launch**2)/100.0
     Ch2_Gain_low = Ch2_Gain_Low_0*(100.0+Ch2_Degrad_Low_1*Time_Temp_Since_Launch+ &
                                  Ch2_Degrad_Low_2*Time_Temp_Since_Launch**2)/100.0
     Ch2_Gain_High = Ch2_Gain_High_0*(100.0+Ch2_Degrad_High_1*Time_Temp_Since_Launch+ &
                                  Ch2_Degrad_High_2*Time_Temp_Since_Launch**2)/100.0
     Ch3a_Gain_low = Ch3a_Gain_Low_0*(100.0 + Ch3a_Degrad_Low_1 * Time_Temp_Since_Launch + &
                                  Ch3a_Degrad_Low_2 * Time_Temp_Since_Launch**2)/100.0
     Ch3a_Gain_High = Ch3a_Gain_High_0*(100.0+Ch3a_Degrad_High_1*Time_Temp_Since_Launch + &
                                  Ch3a_Degrad_High_2 * Time_Temp_Since_Launch**2)/100.0

     !---- compute reflectances at switch count
     Ref_Ch1_Switch = Ch1_Gain_low * (Ch1_Switch_Count_Cal - Ch1_Dark_Count_Cal)
     Ref_Ch2_Switch = Ch2_Gain_low * (Ch2_Switch_Count_Cal - Ch2_Dark_Count_Cal)
     Ref_Ch6_Switch = Ch3a_Gain_low * (Ch3a_Switch_Count_Cal - Ch3a_Dark_Count_Cal)

     !--- test for Valid dark count
     if (abs(Ch1_Dark_Count_Cal - Ch1_Dark_Count) > 5) then
        Valid_Ref_Cal = .false.
     endif

   endif

   !---- determine if coefficients are Valid
   if ((Ch1_Gain_low < 0.0) .or. (Ch1_Gain_High < 0.0)) then
        Valid_Ref_Cal = .false.
   endif


   !--- Apply Calibration
   if (Valid_Ref_Cal .eqv. .true.) then

     !--- channel 1
     where (Chan_Counts_Avhrr(1,:,Line_Idx) <= Ch1_Switch_Count_Cal)
       ch(1)%Ref_Toa(:,Line_Idx) = Ch1_Gain_low*(Chan_Counts_Avhrr(1,:,Line_Idx) - Ch1_Dark_Count_Cal) 
     elsewhere
       ch(1)%Ref_Toa(:,Line_Idx) = Ref_Ch1_Switch + Ch1_Gain_High*(Chan_Counts_Avhrr(1,:,Line_Idx) - Ch1_Switch_Count_Cal) 
     end where

     !--- channel 2
     where (Chan_Counts_Avhrr(2,:,Line_Idx) <= Ch2_Switch_Count_Cal)
       ch(2)%Ref_Toa(:,Line_Idx) = Ch2_Gain_low*(Chan_Counts_Avhrr(2,:,Line_Idx) - Ch2_Dark_Count_Cal) 
     elsewhere
       ch(2)%Ref_Toa(:,Line_Idx) = Ref_Ch2_Switch + Ch2_Gain_High*(Chan_Counts_Avhrr(2,:,Line_Idx) - Ch2_Switch_Count_Cal) 
     end where

     !--- channel 3a
     if (Ch3a_On_AVHRR(Line_Idx) == sym%NO) then
       ch(6)%Ref_Toa(:,Line_Idx) = Missing_Value_Real4 
     else 
       where (Chan_Counts_Avhrr(3,:,Line_Idx) <= Ch3a_Switch_Count_Cal)
         ch(6)%Ref_Toa(:,Line_Idx) = Ch3a_Gain_low*(Chan_Counts_Avhrr(3,:,Line_Idx) - Ch3a_Dark_Count_Cal) 
       elsewhere
         ch(6)%Ref_Toa(:,Line_Idx) = Ref_Ch6_Switch + Ch3a_Gain_High*(Chan_Counts_Avhrr(3,:,Line_Idx) - Ch3a_Switch_Count_Cal) 
       end where
     endif

  else

    ch(1)%Ref_Toa(:,Line_Idx) = Missing_Value_Real4 
    ch(2)%Ref_Toa(:,Line_Idx) = Missing_Value_Real4 
    ch(6)%Ref_Toa(:,Line_Idx) = Missing_Value_Real4 

  endif

end subroutine REF_CAL

!-----------------------------------------------------------
! THERMAL CALIBRATION
!
! variables
!    Therm_Cal_1b - use 1b cal. coefficients (true)
!    AVHRR_KLM_Flag - logical flag set to true if this AVHRR_KLM_Flag data
!    Ch3a_On_AVHRR_Temp - integer flag set to 0 if Ch3a is off
!    count3 - vector of channel 3 counts
!    count4 - vector of channel 4 counts
!    count5 - vector of channel 5 counts
!    d_coef - the "d" coefficient as illustrated in AVHRR_KLM_Flag-POD
!    a_coef - the "a" coefficient as illustrated in AVHRR_KLM_Flag-POD
!    b_coef - the "b" coefficient as illustrated in AVHRR_KLM_Flag-POD
!    linear_slope - slope in linear radiance computation
!    linear_intercept - intercept in linear radiance computation
!    nu_20, a1_20, a2_20 - radiance to temperature conversion parameters
!                        (see AVHRR_KLM_Flag-PODUG for example)
!    nu_31, a1_31, a2_31 - radiance to temperature conversion parameters
!    nu_32, a1_32, a2_32 - radiance to temperature conversion parameters
!
!
! output 
!   Rad_Ch20_Temp, Rad_Ch31_Temp, Rad_Ch32_Temp - vector of channel radiance in mW/m^2/str/cm^-1
!   Bt_Ch20_Temp, Bt_Ch31_Temp, Bt_Ch32_Temp - vector of brightness temperature in K
!
! Some Notes
! 1- pre-AVHRR_KLM_Flag data is first converted to a radiance assuming a linear calibration
!    then the non-linear adjustment is made
! 2- Does not incorporate the correct procedure for data
!    before NOAA-14.  Need to review original CLAVR-1 for this.
! 3 - AVHRR_KLM_Flag data, the non-linear correction is added to the linear radiance 
!     but for non-AVHRR_KLM_Flag data is the non-linear correction is additive (see pods)
!     this is reflected in the difference is the pre-AVHRR_KLM_Flag and AVHRR_KLM_Flag coefs.
!
! Revision History
!  Coded January 2003 - A. Heidinger
!        February 2003 - added Therm_Cal_1b argument and new calibration
!
!
! some clarification
!  the AVHRR_KLM_Flag-POD nomenclature has changed.  
!
!---------------------------------------------------------------------------
subroutine THERM_CAL(Num_Scans_Read)

  integer, intent(in):: Num_Scans_Read
  integer(kind=int4):: i
  integer(kind=int4):: j
  integer, parameter:: dx = 2
  integer, parameter:: dy = 0
  integer:: Num_Valid
  integer:: Elem_Idx_Min
  integer:: Elem_Idx_Max
  integer:: Line_Idx_Max
  integer:: Line_Idx_Min
  integer:: Elem_Idx
  integer:: Line_Idx
  integer, parameter:: Ch20_Counts_Filter_Thresh = 900

  !------------------------------------------------------------------------
  ! Perform Filter on Ch20 counts for appropriate sensors
  ! Do this only for TIROS-N, NOAA-5,6,7,8,9,10
  !------------------------------------------------------------------------
  if ((Sc_Id_WMO >= 706 .and. Sc_Id_WMO <= 708) .or. (Sc_Id_WMO >= 200 .and. Sc_Id_WMO <= 202)) then

       Temp_Pix_Array = Chan_Counts_Avhrr(3,:,:)
       Ch20_Counts_Filtered = Chan_Counts_Avhrr(3,:,:)
       One_Byte_Temp = 1
       where(Temp_Pix_Array <= 100 .or. Temp_Pix_Array > 1024)
              One_Byte_Temp = 0
       endwhere

       do Elem_Idx = 1, Num_Pix

          Elem_Idx_Min = min(Num_Pix,max(1,Elem_Idx - dx))
          Elem_Idx_Max = min(Num_Pix,max(1,Elem_Idx + dx))

         do Line_Idx = 1, Num_Scans_Read

           if (Chan_Counts_Avhrr(3,Elem_Idx,Line_Idx) > Ch20_Counts_Filter_Thresh) then
             Line_Idx_Min = min(Num_Scans_Read,max(1,Line_Idx - dy))
             Line_Idx_Max = min(Num_Scans_Read,max(1,Line_Idx + dy))
             Num_Valid = sum(One_Byte_Temp(Elem_Idx_Min:Elem_Idx_Max,Line_Idx_Min:Line_Idx_Max))
             if (Num_Valid > 1) then
              Ch20_Counts_Filtered(Elem_Idx,Line_Idx) =  &
                    sum(One_Byte_Temp(Elem_Idx_Min:Elem_Idx_Max,Line_Idx_Min:Line_Idx_Max)* &
                        Temp_Pix_Array(Elem_Idx_Min:Elem_Idx_Max,Line_Idx_Min:Line_Idx_Max)) / Num_Valid
             endif
           endif
          enddo
       enddo
       Chan_Counts_Avhrr(3,:,:) = int(Ch20_Counts_Filtered,kind=int2)

  endif 

  !------------------------------------------------------------------------
  ! Loop through each line and calibrate
  !------------------------------------------------------------------------
  do j = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1

    if (Bad_Scan_Flag(j) == sym%YES) then
       cycle
    endif 

    !----  convert counts to radiance
    if (Therm_Cal_1b == sym%YES) then

       if (AVHRR_KLM_Flag == sym%YES) then       !apply level-1b non-linear correction

         if (Ch3a_On_AVHRR(j) == 0 .and. Chan_On_Flag_Default(20) == sym%YES) then
           ch(20)%Rad_Toa(:,j) = ir_coef_1_1b(3,j) + ir_coef_2_1b(3,j)*Chan_Counts_Avhrr(3,:,j) +  &
                        ir_coef_3_1b(3,j)*Chan_Counts_Avhrr(3,:,j)**2
         endif
         if (Chan_On_Flag_Default(31) == sym%YES) then
             ch(31)%Rad_Toa(:,j) = ir_coef_1_1b(4,j) + ir_coef_2_1b(4,j)*Chan_Counts_Avhrr(4,:,j) +  &
                       ir_coef_3_1b(4,j)*Chan_Counts_Avhrr(4,:,j)**2
         endif
         if (Chan_On_Flag_Default(32) == sym%YES) then
            ch(32)%Rad_Toa(:,j) = ir_coef_1_1b(5,j) + ir_coef_2_1b(5,j)*Chan_Counts_Avhrr(5,:,j) +  &
                            ir_coef_3_1b(5,j)*Chan_Counts_Avhrr(5,:,j)**2
         endif

       else                         !compute level-1b linear correction for pre-AVHRR_KLM_Flag

         !--- Eqn 7.1.2.4-7
         if (Chan_On_Flag_Default(20) == sym%YES) then
           ch(20)%Rad_Toa(:,j) = ir_linear_slope_1b(3,j)*Chan_Counts_Avhrr(3,:,j) + ir_linear_intercept_1b(3,j)
         endif
         if (Chan_On_Flag_Default(31) == sym%YES) then
           ch(31)%Rad_Toa(:,j) = ir_linear_slope_1b(4,j)*Chan_Counts_Avhrr(4,:,j) + ir_linear_intercept_1b(4,j)
         endif
         if (Chan_On_Flag_Default(32) == sym%YES) then
           ch(32)%Rad_Toa(:,j) = ir_linear_slope_1b(5,j)*Chan_Counts_Avhrr(5,:,j) + ir_linear_intercept_1b(5,j)
         endif

       endif

    else    !generate linear radiance using internally computed coefficients

      if (Ch3a_On_AVHRR(j) == 0 .and. Chan_On_Flag_Default(20) == sym%YES) then
         ch(20)%Rad_Toa(:,j) = IR_Linear_Slope_New(3,j)*Chan_Counts_Avhrr(3,:,j) + IR_Linear_Intercept_New(3,j)
      endif
      if (Chan_On_Flag_Default(31) == sym%YES) then
        ch(31)%Rad_Toa(:,j) = IR_Linear_Slope_New(4,j)*Chan_Counts_Avhrr(4,:,j) + IR_Linear_Intercept_New(4,j)
      endif
      if (Chan_On_Flag_Default(32) == sym%YES) then
        ch(32)%Rad_Toa(:,j) = IR_Linear_Slope_New(5,j)*Chan_Counts_Avhrr(5,:,j) + IR_Linear_Intercept_New(5,j)
      endif
    endif

   !--- apply non-linear correction where appropriate (non-AVHRR_KLM_Flag or any if thermal cal redone)
   if ((AVHRR_KLM_Flag == sym%NO) .or. (Therm_Cal_1b == sym%NO)) then 

     if (Ch3a_On_AVHRR(j) == 0 .and. Chan_On_Flag_Default(20) == sym%YES) then
       ch(20)%Rad_Toa(:,j) = ch(20)%Rad_Toa(:,j) + b0_3b + b1_3b*ch(20)%Rad_Toa(:,j) + b2_3b*ch(20)%Rad_Toa(:,j)**2
     endif
     if (Chan_On_Flag_Default(31) == sym%YES) then
       ch(31)%Rad_Toa(:,j) = ch(31)%Rad_Toa(:,j) + b0_4 + b1_4*ch(31)%Rad_Toa(:,j) + b2_4*ch(31)%Rad_Toa(:,j)**2
     endif
     if (Chan_On_Flag_Default(32) == sym%YES) then
       ch(32)%Rad_Toa(:,j) = ch(32)%Rad_Toa(:,j) + b0_5 + b1_5*ch(32)%Rad_Toa(:,j) + b2_5*ch(32)%Rad_Toa(:,j)**2
     endif
   endif

   !---- convert radiance to brightness temperature

    do i = 1, Num_Pix

        if (Chan_On_Flag_Default(20) == sym%YES) then
         if (ch(20)%Rad_Toa(i,j) > 0.0) then
          ch(20)%Bt_Toa(i,j) = PLANCK_TEMP_FAST(20,ch(20)%Rad_Toa(i,j)) 
         else
          ch(20)%Rad_Toa(i,j) = Missing_Value_Real4
          ch(20)%Bt_Toa(i,j) = Missing_Value_Real4
         endif
        endif

        if (Chan_On_Flag_Default(31) == sym%YES) then
         if (ch(31)%Rad_Toa(i,j) > 0.0) then
           ch(31)%Bt_Toa(i,j) = PLANCK_TEMP_FAST(31,ch(31)%Rad_Toa(i,j)) 
         else
           ch(31)%Rad_Toa(i,j) = Missing_Value_Real4
           ch(31)%Bt_Toa(i,j) = Missing_Value_Real4
         endif
        endif

        if (Chan_On_Flag_Default(32) == sym%YES) then
         if (ch(32)%Rad_Toa(i,j) > 0.0) then
           ch(32)%Bt_Toa(i,j) = PLANCK_TEMP_FAST(32,ch(32)%Rad_Toa(i,j)) 
         else
           ch(32)%Rad_Toa(i,j) = Missing_Value_Real4
           ch(32)%Bt_Toa(i,j) = Missing_Value_Real4
         endif
        endif

    end do

  enddo

  !----- additional quality check on Bt_Ch31
  do j = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1
    if (minval(ch(31)%Bt_Toa(:,j)) < 0.0) then
      Bad_Scan_Flag(j) = sym%YES
    endif
  enddo

end subroutine THERM_CAL

!--------------------------------------------------------------------------
! compute parameters needed to derive new thermal calibration
! coefficients using a method similar to that described in the
! PODUG for HPRT users.  
!
! The methodology is 
!  average all space, blackbody counts and PRT's rom one scan together
!  average all quantities together for Ncal scanlines
!  every Ncal scanlines, update a running mean of all these values
!  the running mean is defined as
!   new_running_mean = (1-Cal_Smoothing_Weight)*old_running_mean + new_value 
!
!  experimentation has shown Ncal = 10 and Cal_Smoothing_Weight = 1/5 result in smooth curves
!  the variation 1b slope and intercept is approximated with Ncal=5, Cal_Smoothing_Weight = 1.0
!  unsure if smoother is better (TBD).
!
! Author: Andrew Heidinger
! Date: February 2003
!
! Revision - added check for missing space and blackbody counts - noticed
!           this occured for D95001 - D95005 data for NOAA-14
!
!  - modified for one scan at a time processing - akh
!
! Reference:
! http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/c7/sec7-1.htm
!--------------------------------------------------------------------------
subroutine COMPUTE_NEW_THERM_CAL_COEF(Line_Idx)

 integer, intent(in):: Line_Idx
 integer:: iword
 integer:: Byte_Idx
 integer:: j
 integer:: Pix_Idx
 integer:: Chan_Idx
 integer(kind=int4):: i4word
 real(kind=real4):: Mean_Space_Count_1_Temp
 real(kind=real4):: Mean_Space_Count_2_Temp
 real(kind=real4):: Mean_Space_Count_3_Temp
 real(kind=real4):: Mean_Space_Count_4_Temp
 real(kind=real4):: Mean_Space_Count_5_Temp
 real(kind=real4):: Mean_BB_Count_3_Temp
 real(kind=real4):: Mean_BB_Count_4_Temp
 real(kind=real4):: Mean_BB_Count_5_Temp
 real(kind=real4):: Space_Count_Temp_4_test

 Valid_Therm_Cal = .true.

!-----------------------------------------------------------------------
! unpack PRT readings
!-----------------------------------------------------------------------

!-- reset prt readings

  Sum_Prt = 0

  !-----------------------------------------------------------------------
  ! unpack telemetry words for AVHRR_KLM_Flag
  !-----------------------------------------------------------------------
  if (AVHRR_KLM_Flag == sym%YES) then
    Byte_Idx = 1091 
    Sum_Prt = MAKE_I2WORD(Buffer_AVHRR(Byte_Idx:Byte_Idx+1),sym%SIGNED,Byte_Swap_1b) +  &
              MAKE_I2WORD(Buffer_AVHRR(Byte_Idx+2:Byte_Idx+3),sym%SIGNED,Byte_Swap_1b) + &
              MAKE_I2WORD(Buffer_AVHRR(Byte_Idx+4:Byte_Idx+5),sym%SIGNED,Byte_Swap_1b)
  else
  !-----------------------------------------------------------------------
  ! unpack telemetry words for pre-AVHRR_KLM_Flag
  !-----------------------------------------------------------------------
    iword = 1
    do Byte_Idx = 309,448,4
      i4word = MAKE_I4WORD(Buffer_AVHRR(Byte_Idx:Byte_Idx+3),sym%UNSIGNED,Byte_Swap_1b)
      do j = 1,3
       Telemetry_Word(iword) = i4word / (1024**(3-j))
       i4word = i4word - Telemetry_Word(iword)*(1024**(3-j))
       iword = iword + 1
       if (iword > 103) then
        exit
       endif
      enddo
       if (iword > 103) then
        exit
       endif
    enddo

    !--- prt readings
    Sum_Prt = Telemetry_Word(18)+Telemetry_Word(19)+Telemetry_Word(20)

  endif

  !--- increment prt counter
  Prt_Idx = Prt_Idx + 1

  !---- check for reset of PRTs (Reset by Prt_Idx = 0 )
  if (Sum_Prt < 50) then
     Prt_Idx = 0
  endif

  !--------------------------------------------------
  ! test to see if any coefficients can be made
  !--------------------------------------------------
  if ((Prt_Idx == 0) .and. (Mean_T_Prt <= 0.0)) then
        print *, EXE_PROMPT, "Error: PRT Temperatures bad "
        Valid_Therm_Cal = .false.
        return
  endif

  !---- check to see if out of cycle
  if (Prt_Idx > 4) then
      print *, EXE_PROMPT, "Error: PRTs out of order "
      Valid_Therm_Cal = .false.
      Prt_Idx = 0    ! added akh
      return
  endif

  !----------------------------------------------------
  !---- unpack target view
  !----------------------------------------------------
  if (AVHRR_KLM_Flag == sym%YES) then

    !--- the internal target data from AVHRR_KLM_Flag
    do j = 1,10
      Blackbody_Count(3,j) = MAKE_I2WORD(Buffer_AVHRR(1101+6*(j-1):1102+6*(j-1)),sym%SIGNED,Byte_Swap_1b)
      Blackbody_Count(4,j) = MAKE_I2WORD(Buffer_AVHRR(1103+6*(j-1):1104+6*(j-1)),sym%SIGNED,Byte_Swap_1b)
      Blackbody_Count(5,j)  = MAKE_I2WORD(Buffer_AVHRR(1105+6*(j-1):1106+6*(j-1)),sym%SIGNED,Byte_Swap_1b)
    enddo
  else

    !--- the internal target data from pre-AVHRR_KLM_Flag
    Chan_Idx = 3
    Pix_Idx = 1
    do iword = 23,52,1
       Blackbody_Count(Chan_Idx, Pix_Idx) = Telemetry_Word(iword)
       Chan_Idx = Chan_Idx + 1
       if (Chan_Idx > 5) then
          Chan_Idx = 3
          Pix_Idx = Pix_Idx + 1
       endif
     enddo

  endif

  !------------------------------------------------------
  !---- unpack space view
  !------------------------------------------------------

  if (AVHRR_KLM_Flag == sym%YES) then

    !---space data from AVHRR_KLM_Flag data
    do j = 1,10
      Space_Count_Temp(1,j)  = MAKE_I2WORD(Buffer_AVHRR(1161+10*(j-1):1162+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
      Space_Count_Temp(2,j) = MAKE_I2WORD(Buffer_AVHRR(1163+10*(j-1):1164+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
      Space_Count_Temp(3,j) = MAKE_I2WORD(Buffer_AVHRR(1165+10*(j-1):1166+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
      Space_Count_Temp(4,j) = MAKE_I2WORD(Buffer_AVHRR(1167+10*(j-1):1168+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
      Space_Count_Temp(5,j)  = MAKE_I2WORD(Buffer_AVHRR(1169+10*(j-1):1170+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
    enddo

  else

    !---space data from pre-AVHRR_KLM_Flag data
    Chan_Idx = 1
    Pix_Idx = 1
    do iword = 53,102,1
       Space_Count_Temp(Chan_Idx, Pix_Idx) = Telemetry_Word(iword)
       Chan_Idx = Chan_Idx + 1
       if (Chan_Idx > 5) then
          Chan_Idx = 1
          Pix_Idx = Pix_Idx + 1
       endif
    enddo

  endif



   !---- check for rapid change in ch4 space count (ie during a moon view)
   Space_Count_Temp_4_test = sum(Space_Count_Temp(4,:))/10.0
   if ((Mean_Space_Count_4 > 0) .and.  &
       (abs(Space_Count_Temp_4_test - Mean_Space_Count_4) > 4.0))  then
       print *, EXE_PROMPT, "Error: rapid change in ch4 space count "
       Valid_Therm_Cal = .false.
       return
   endif

   !---- check for bad space or blackbody counts
   if ((minval(Space_Count_Temp(4:5,:)) < 50) .or. (minval(Blackbody_Count(4:5,:)) < 50)) then
       print *, EXE_PROMPT, "Error: ch4 -ch5 space count threshold exceeded "
       Valid_Therm_Cal = .false.
       return
   endif

    !----- assume a Valid calibration can be made
    if (Valid_Therm_Cal .eqv. .true.) then
     Cal_Idx = Cal_Idx + 1
    endif

    !--- compute blackbody temperature if Valid PRT reading (Eq 7.1.2.4-1 in NOAA KLM User's Guide)
    if ((Prt_Idx > 0) .and. (Valid_Therm_Cal .eqv. .true.)) then
       Mean_Prt = Sum_Prt/3.0
       T_Prt(Cal_Idx) = Prt_Coef(0,Prt_Idx) +  &
                     Prt_Coef(1,Prt_Idx)*Mean_Prt + &
                     Prt_Coef(2,Prt_Idx)*Mean_Prt**2 + &
                     Prt_Coef(3,Prt_Idx)*Mean_Prt**3 + &
                     Prt_Coef(4,Prt_Idx)*Mean_Prt**4

    endif

    !--- if no PRT reading, fill in with mean value
    if ((Prt_Idx == 0) .and. (Valid_Therm_Cal .eqv. .true.)) then
       T_Prt(Cal_Idx) = Mean_T_Prt
    endif

    !--- average space and bb counts for this scan
    if (Cal_Idx > 0) then
       Space_Count_Temp_1(Cal_Idx) = sum(Space_Count_Temp(1,:))/10.0
       Space_Count_Temp_2(Cal_Idx) = sum(Space_Count_Temp(2,:))/10.0
       Space_Count_Temp_3(Cal_Idx) = sum(Space_Count_Temp(3,:))/10.0
       Space_Count_Temp_4(Cal_Idx) = sum(Space_Count_Temp(4,:))/10.0
       Space_Count_Temp_5(Cal_Idx) = sum(Space_Count_Temp(5,:))/10.0

       BB_Count_3(Cal_Idx) = sum(Blackbody_Count(3,:))/10.0
       BB_Count_4(Cal_Idx) = sum(Blackbody_Count(4,:))/10.0
       BB_Count_5(Cal_Idx) = sum(Blackbody_Count(5,:))/10.0
    endif

!-------------------------------------------------------------------
! compute running means of space, blackbody and prt temps
!-------------------------------------------------------------------

    !--- if initial spinup period, sum all to get mean
    if (Spinup_New_Therm_Cal .eqv. .true.) then

      if (Cal_Idx > 0) then
       Mean_T_Prt = sum(T_Prt(1:Cal_Idx))/Cal_Idx                      ! Ref KLM-POD 7.1.2.4-2
       Mean_Space_Count_1 = sum(Space_Count_Temp_1(1:Cal_Idx))/Cal_Idx
       Mean_Space_Count_2 = sum(Space_Count_Temp_2(1:Cal_Idx))/Cal_Idx
       Mean_Space_Count_3 = sum(Space_Count_Temp_3(1:Cal_Idx))/Cal_Idx
       Mean_Space_Count_4 = sum(Space_Count_Temp_4(1:Cal_Idx))/Cal_Idx
       Mean_Space_Count_5 = sum(Space_Count_Temp_5(1:Cal_Idx))/Cal_Idx
       Mean_BB_Count_3 = sum(BB_Count_3(1:Cal_Idx))/Cal_Idx
       Mean_BB_Count_4 = sum(BB_Count_4(1:Cal_Idx))/Cal_Idx
       Mean_BB_Count_5 = sum(BB_Count_5(1:Cal_Idx))/Cal_Idx
      else
       Mean_T_Prt = Missing_Value_Real4
       Mean_Space_Count_1 = Missing_Value_Real4
       Mean_Space_Count_2 = Missing_Value_Real4
       Mean_Space_Count_3 = Missing_Value_Real4
       Mean_Space_Count_4 = Missing_Value_Real4
       Mean_Space_Count_5 = Missing_Value_Real4
       Mean_BB_Count_3 = Missing_Value_Real4
       Mean_BB_Count_4 = Missing_Value_Real4
       Mean_BB_Count_5 = Missing_Value_Real4
       Valid_Therm_Cal = .false.
      endif

      if (Cal_Idx == Ncal) then
        Spinup_New_Therm_Cal = .false.
        Cal_Idx = 0
        T_Prt = 0.0
        Space_Count_Temp_1 = 0.0
        Space_Count_Temp_2 = 0.0
        Space_Count_Temp_3 = 0.0
        Space_Count_Temp_4 = 0.0
        Space_Count_Temp_5 = 0.0
        BB_Count_3 = 0.0
        BB_Count_4 = 0.0
        BB_Count_5 = 0.0
      endif

    else

      !--- after Ncal lines, update space and bb count values 
      if (Cal_Idx == Ncal) then

       !--- compute running-mean of the black-body temperature
       Mean_T_Prt = (1.0-Cal_Smoothing_Weight)*Mean_T_Prt + Cal_Smoothing_Weight*(sum(T_Prt)/Ncal)

       !--- compute instantaneous mean of the space counts
       Mean_Space_Count_1_Temp = sum(Space_Count_Temp_1)/Ncal
       Mean_Space_Count_2_Temp = sum(Space_Count_Temp_2)/Ncal
       Mean_Space_Count_3_Temp = sum(Space_Count_Temp_3)/Ncal
       Mean_Space_Count_4_Temp = sum(Space_Count_Temp_4)/Ncal
       Mean_Space_Count_5_Temp = sum(Space_Count_Temp_5)/Ncal

       !--- compute instantaneous mean of the black-body counts
       Mean_BB_Count_3_Temp = sum(BB_Count_3)/Ncal
       Mean_BB_Count_4_Temp = sum(BB_Count_4)/Ncal
       Mean_BB_Count_5_Temp = sum(BB_Count_5)/Ncal

       !--- compute running mean of space counts
       Mean_Space_Count_1 = (1.0-Cal_Smoothing_Weight)*Mean_Space_Count_1 + Cal_Smoothing_Weight*Mean_Space_Count_1_Temp
       Mean_Space_Count_2 = (1.0-Cal_Smoothing_Weight)*Mean_Space_Count_2 + Cal_Smoothing_Weight*Mean_Space_Count_2_Temp
       Mean_Space_Count_3 = (1.0-Cal_Smoothing_Weight)*Mean_Space_Count_3 + Cal_Smoothing_Weight*Mean_Space_Count_3_Temp
       Mean_Space_Count_4 = (1.0-Cal_Smoothing_Weight)*Mean_Space_Count_4 + Cal_Smoothing_Weight*Mean_Space_Count_4_Temp
       Mean_Space_Count_5 = (1.0-Cal_Smoothing_Weight)*Mean_Space_Count_5 + Cal_Smoothing_Weight*Mean_Space_Count_5_Temp

       !--- compute running mean of black-body counts
       Mean_BB_Count_3 = (1.0-Cal_Smoothing_Weight)*Mean_BB_Count_3 + Cal_Smoothing_Weight*Mean_BB_Count_3_Temp
       Mean_BB_Count_4 = (1.0-Cal_Smoothing_Weight)*Mean_BB_Count_4 + Cal_Smoothing_Weight*Mean_BB_Count_4_Temp
       Mean_BB_Count_5 = (1.0-Cal_Smoothing_Weight)*Mean_BB_Count_5 + Cal_Smoothing_Weight*Mean_BB_Count_5_Temp


       !--- handle ch3 differently due to Ch3a/ch3b switching
       !--- if a switch detected, reset running means with instantataneous values

       if ((Ch3a_On_AVHRR(Line_Idx) == 1) .and. (Mean_Space_Count_3 > 50)) then
          Mean_Space_Count_3 = Mean_Space_Count_3_Temp
          Mean_BB_Count_3 =  Mean_BB_Count_3_Temp
       endif

       if ((Ch3a_On_AVHRR(Line_Idx) == 0) .and. (Mean_Space_Count_3 < 950)) then
          Mean_Space_Count_3 = Mean_Space_Count_3_Temp
          Mean_BB_Count_3 =  Mean_BB_Count_3_Temp
       endif

       Cal_Idx = 0
       T_Prt = 0.0
       Space_Count_Temp_1 = 0.0
       Space_Count_Temp_2 = 0.0
       Space_Count_Temp_3 = 0.0
       Space_Count_Temp_4 = 0.0
       Space_Count_Temp_5 = 0.0
       BB_Count_3 = 0.0
       BB_Count_4 = 0.0
       BB_Count_5 = 0.0

      endif

    endif


    !---  look for space counts that inValid for Ch3a/3b
    if ((Mean_Space_Count_3 > 50.0) .and. (Mean_Space_Count_3 < 950.0)) then
             Valid_Therm_Cal = .false.
    endif

    if (Valid_Therm_Cal .eqv. .true.) then

       !--- blackbody radiance
       NBB_3b = PLANCK_RAD_FAST(20,Mean_T_Prt)
       NBB_4 = PLANCK_RAD_FAST(31,Mean_T_Prt)
       NBB_5 = PLANCK_RAD_FAST(32,Mean_T_Prt)

       !--- compute new linear slopes     (Ref 7.1.2.4-5)
       if (Ch3a_On_AVHRR(Line_Idx) == 0) then
        IR_Linear_Slope_New(3,Line_Idx) = (Space_Rad_3b - NBB_3b) / &
                                       (Mean_Space_Count_3 - Mean_BB_Count_3)
        IR_Linear_Intercept_New(3,Line_Idx) = Space_Rad_3b - &
                  IR_Linear_Slope_New(3,Line_Idx)*Mean_Space_Count_3
       endif


       IR_Linear_Slope_New(4,Line_Idx) = (Space_Rad_4 - NBB_4) / &
                                (Mean_Space_Count_4 - Mean_BB_Count_4)

       IR_Linear_Intercept_New(4,Line_Idx) = Space_Rad_4 - &
                  IR_Linear_Slope_New(4,Line_Idx)*Mean_Space_Count_4

       IR_Linear_Slope_New(5,Line_Idx) = (Space_Rad_5 - NBB_5) / &
                                (Mean_Space_Count_5 - Mean_BB_Count_5)
       IR_Linear_Intercept_New(5,Line_Idx) = Space_Rad_5 - &
                  IR_Linear_Slope_New(5,Line_Idx)*Mean_Space_Count_5

     else
        Valid_Therm_Cal = .false.
     endif

     !--- store these into public arrays for access by other routines
     Scan_Space_Counts_Avhrr(1,Line_Idx) = Mean_Space_Count_1
     Scan_Space_Counts_Avhrr(2,Line_Idx) = Mean_Space_Count_2
     Scan_Space_Counts_Avhrr(3,Line_Idx) = Mean_Space_Count_3
     Scan_Space_Counts_Avhrr(4,Line_Idx) = Mean_Space_Count_4
     Scan_Space_Counts_Avhrr(5,Line_Idx) = Mean_Space_Count_5
     Scan_BB_Counts_Avhrr(3,Line_Idx) = Mean_BB_Count_3
     Scan_BB_Counts_Avhrr(4,Line_Idx) = Mean_BB_Count_4
     Scan_BB_Counts_Avhrr(5,Line_Idx) = Mean_BB_Count_5

end subroutine COMPUTE_NEW_THERM_CAL_COEF

!-----------------------------------------------------------------
! functions to convert 1 byte integers into 2 or 4 byte integers
! if the integeris signed (isign = 0) then the first byte
! remains signed
!---------------------------------------------------------------- 
 function MAKE_I2WORD(byte1,isign,iswap) result(i2word)
  integer(kind=int1), intent(in), dimension(:):: byte1
  integer(kind=int4), intent(in):: isign
  integer(kind=int4), intent(in):: iswap
  integer(kind=int2):: i2word
  integer(kind=int2), dimension(:), allocatable:: byte2
  integer:: i,n

  n = size(byte1)

  allocate(byte2(n))

  !--- swap bytes if necessary
  if (iswap == sym%YES) then
    do i = 1, n
      byte2(i) = byte1(n-i+1) 
    enddo
  else
    byte2 = byte1
  endif

  !--- account for unsigned integers
  if ((byte2(1) < 0) .and. (isign == 0)) then
     byte2(1) = byte2(1) + 256
  endif
  do i = 2,n
     if ((byte2(i) < 0)) then
        byte2(i) = byte2(i) + 256
     end if
  end do

!--- make word
   i2word = byte2(1)*256 + byte2(2)

   deallocate(byte2)

 end function MAKE_I2WORD

!-----
 function MAKE_I4WORD(byte1,isign,iswap) result(i4word)
  integer(kind=int1), intent(in), dimension(:):: byte1
  integer(kind=int4), intent(in):: isign
  integer(kind=int4), intent(in):: iswap
  integer(kind=int4):: i4word
  integer(kind=int2), dimension(:), allocatable:: byte2
  integer:: i,n

   n = size(byte1)

   allocate(byte2(n))

!--- swap bytes if necessary
  if (iswap == sym%YES) then
   do i = 1, n
    byte2(i) = byte1(n-i+1) 
   enddo
  else
   byte2 = byte1
  endif

!--- account for unsigned integers
   if ((byte2(1) < 0) .and. (isign == sym%NO)) then
     byte2(1) = byte2(1) + 256
   endif
   do i = 2,n
    if ((byte2(i) < 0)) then
     byte2(i) = byte2(i) + 256
    end if
   end do
   i4word = 0
   do i = n,1,-1
     i4word = i4word + byte2(i)*(256**(n-i))
   enddo

   deallocate(byte2)

 end function MAKE_I4WORD


 function i4word_to_string(i4word,n) result(char_word)
    integer(kind=int4), intent(in):: i4word
    integer, intent(in):: n
    character(len=n):: char_word
    integer(kind=int4):: i4temp
    integer:: i

    i4temp = i4word
    char_word = ""
    do i = 1,n
     char_word(i:i) = char(i4temp/(10**(n-i))+48)
     i4temp = i4temp - 10**(n-i)*(i4temp / (10**(n-i)))
    enddo
 end function i4word_to_string

!======================================================================
! unpack avhrr header record - routine should work for lac and AVHRR_GAC_Flag on
! all pre-AVHRR_KLM_Flag sensors
!
! Modified December 2008 for segment processing
!======================================================================
subroutine UNPACK_AVHRR_HEADER_RECORD(Sc_Id,data_type,Start_Year_Temp,&
              Start_day_Temp,Start_Time_Temp,Num_Scans,End_Year_Temp,End_day_Temp,End_Time_Temp, &
              tip_parity,aux_sync,ramp_auto_Cal,proc_block_Id,ver_1b)

   integer(kind=int2), intent(out):: Sc_Id,data_type,Start_Year_Temp,Start_day_Temp, &
                                     End_Year_Temp,End_day_Temp, &
                                     tip_parity,aux_sync,ramp_auto_Cal,ver_1b
   integer(kind=int4), intent(out):: Start_Time_Temp,End_Time_Temp,Num_Scans
   character(len=*), intent(out):: proc_block_Id
   integer(kind=int4):: i4word
   integer(kind=int2), dimension(2):: data_byte
   integer:: i
     
   !-----------------------------------------------------------------------
   ! begin executable code
   !-----------------------------------------------------------------------

   Sc_Id = Header_Buffer_AVHRR(1)
   data_type = ishft(Header_Buffer_AVHRR(2),-4)

   tip_parity = 0
   aux_sync = 0
   ramp_auto_Cal = 0
   ver_1b = 1

   !----------------------------------------------------------------------
   ! unpack starting Time_Temp
   !----------------------------------------------------------------------
   data_byte(1) = Header_Buffer_AVHRR(3)
   data_byte(2) = Header_Buffer_AVHRR(4)
   where (data_byte < 0)
         data_byte = data_byte + 256
   end where
   i4word =  data_byte(1)*256 + data_byte(2)
   Start_Year_Temp = i4word / (2**9)
   Start_day_Temp = i4word - Start_Year_Temp * (2**9)

   i4word = MAKE_I4WORD(Header_Buffer_AVHRR(5:8),sym%UNSIGNED,Byte_Swap_1b)
   Start_Time_Temp = i4word - ((i4word/2**27)*(2**27))

   !---------------------------------------------------------------------
   ! unpack processing block id
   !---------------------------------------------------------------------
   do i = 1,len(proc_block_Id)
    proc_block_Id(i:i) =  char(Header_Buffer_AVHRR(17+(i-1)))
   end do 
       
   !----------------------------------------------------------------------
   ! unpack  number of scan lines
   !----------------------------------------------------------------------
   Num_Scans = MAKE_I4WORD(Header_Buffer_AVHRR(9:10),sym%UNSIGNED,Byte_Swap_1b)

   !----------------------------------------------------------------------
   ! unpack end time
   !----------------------------------------------------------------------
   data_byte(1) = Header_Buffer_AVHRR(11)
   data_byte(2) = Header_Buffer_AVHRR(12)
   where (data_byte < 0)
       data_byte = data_byte + 256
   end where
   i4word =  data_byte(1)*256 + data_byte(2)

   End_Year_Temp = i4word / (2**9)
   End_day_Temp = i4word - Start_Year_Temp * (2**9)

   i4word = MAKE_I4WORD(Header_Buffer_AVHRR(13:16),sym%UNSIGNED,Byte_Swap_1b)
   End_Time_Temp = i4word - ((i4word/2**27)*(2**27))

   !----------------------------------------------------------------------
   ! unpack quality indicators
   !-----------------------------------------------------------------------
!  tip_parity = MAKE_I2WORD(Header_Buffer_AVHRR(139:140),sym%UNSIGNED,Byte_Swap_1b)
!  aux_sync = MAKE_I2WORD(Header_Buffer_AVHRR(141:142),sym%UNSIGNED,Byte_Swap_1b)
!  ramp_auto_Cal =  MAKE_I2WORD(Header_Buffer_AVHRR(187:188),sym%UNSIGNED,Byte_Swap_1b)

 end subroutine UNPACK_AVHRR_HEADER_RECORD
!======================================================================
! NOAA 15+ unpack avhrr header record - routine should work for lac and AVHRR_GAC_Flag
!
! Modified January 2005 for NOAA-N
!
! Modified December 2008 for segment processing
!======================================================================
subroutine UNPACK_AVHRR_HEADER_RECORD_KLM(Sc_Id,data_type,Start_Year_Temp,&
              Start_day_Temp,Start_Time_Temp,Num_Scans,End_Year_Temp,End_day_Temp,End_Time_Temp, &
              tip_parity,aux_sync,ramp_auto_Cal,proc_block_Id,ver_1b)

   integer(kind=int2), intent(out):: Sc_Id,data_type,Start_Year_Temp,Start_day_Temp, &
                                     End_Year_Temp,End_day_Temp, &
                                     tip_parity,aux_sync,ramp_auto_Cal,ver_1b
   integer(kind=int4), intent(out):: Start_Time_Temp,End_Time_Temp,Num_Scans
   character(len=*), intent(out):: proc_block_Id
   integer:: i

!-----------------------------------------------------------------------
! begin executable code
!-----------------------------------------------------------------------
      Sc_Id = MAKE_I2WORD(Header_Buffer_AVHRR(73:74),sym%UNSIGNED,Byte_Swap_1b)
      data_type = MAKE_I2WORD(Header_Buffer_AVHRR(77:78),sym%UNSIGNED,Byte_Swap_1b)

      do i = 1,len(proc_block_Id)
       proc_block_Id(i:i) =  char(Header_Buffer_AVHRR(66+(i-1)))
      end do 

!----------------------------------------------------------------------
! unpack 1b version number
!----------------------------------------------------------------------
      ver_1b = MAKE_I2WORD(Header_Buffer_AVHRR(5:6),sym%UNSIGNED,Byte_Swap_1b)

     !--- correct for noaa-15 bug in version number
     if ((ver_1b == 1) .and. (Sc_Id== 4)) then
       ver_1b = 2
     endif  

     !----------------------------------------------------------------------
     ! unpack starting Time_Temp
     !----------------------------------------------------------------------
     Start_Year_Temp = MAKE_I2WORD(Header_Buffer_AVHRR(85:86),sym%UNSIGNED,Byte_Swap_1b)
     Start_day_Temp = MAKE_I2WORD(Header_Buffer_AVHRR(87:88),sym%UNSIGNED,Byte_Swap_1b)
     Start_Time_Temp = MAKE_I4WORD(Header_Buffer_AVHRR(89:92),sym%UNSIGNED,Byte_Swap_1b)
   
    !----------------------------------------------------------------------
    ! unpack  number of scan lines
    ! note, the POD stores this number as a two byte unsigned integer
    !----------------------------------------------------------------------
    Num_Scans = MAKE_I4WORD(Header_Buffer_AVHRR(129:130),sym%UNSIGNED,Byte_Swap_1b)

    !----------------------------------------------------------------------
    ! unpack end Time_Temp
    !----------------------------------------------------------------------
    End_Year_Temp = MAKE_I2WORD(Header_Buffer_AVHRR(97:98),sym%UNSIGNED,Byte_Swap_1b)
    End_day_Temp = MAKE_I2WORD(Header_Buffer_AVHRR(99:100),sym%UNSIGNED,Byte_Swap_1b)
    End_Time_Temp = MAKE_I4WORD(Header_Buffer_AVHRR(101:104),sym%UNSIGNED,Byte_Swap_1b)

    !----------------------------------------------------------------------
    ! unpack quality indicators - check this
    !-----------------------------------------------------------------------
    tip_parity = MAKE_I2WORD(Header_Buffer_AVHRR(139:140),sym%UNSIGNED,Byte_Swap_1b)
    aux_sync = MAKE_I2WORD(Header_Buffer_AVHRR(141:142),sym%UNSIGNED,Byte_Swap_1b)
    ramp_auto_Cal =  MAKE_I2WORD(Header_Buffer_AVHRR(187:188),sym%UNSIGNED,Byte_Swap_1b)

 end subroutine UNPACK_AVHRR_HEADER_RECORD_KLM

!---------------------------------------------------------------------
! EXPLANATION OF LAC TO GAC GEOLOCATION CONVERSION:
!
! A LAC scan line is converted to GAC line in the following way:
!
!  1. GAC processing starts from the 3rd LAC pixel; 1st and 2nd LAC pixels are
!     not used.
!
!  2. The 1st GAC pixel is formed from the 3rd, 4th, 5th and 6th LAC pixel.
!
!  3. The 7th LAC pixel is skipped and the processing for the 2nd GAC pixel uses
!     the 8th, 9th, 10th and 11th LAC pixel.
!
!  4. LAC pixels contributing to a GAC pixel can be described with a formula:
!
!        l, l+1, l+2, l+3 for l = 5*g-2 where g = 1,..,409
!
!     `l' and `g' being LAC and GAC pixel index, respectively, in a scan line.
!
!  5. Geolocation of the third LAC pixel of the four used for a GAC pixel is
!     assigned to that GAC pixel. For the 1st and 2nd GAC pixels these are the
!     5th and 10th LAX pixels. Or in short:
!
!        l = 5*g (g=1,..,409)
!
!  6. The 2048th LAC pixel is not used for GAC pixel processing. The processing
!     ends with the 2047th LAX pixels, for a total of 2045 LAX pixels used
!     (2045/5=409).
!
! This subroutine performs a five-point Lagrangian interpolation for 2048 LAC
! pixels regardless of the actual length of the scan line (GAC or LAC). If it's
! a GAC scan line, the geolocations of the 2nd and 3rd LAC pixels for each GAC
! pixel are averaged to yield the final GAC pixel's geolocation.
!---------------------------------------------------------------------
 subroutine LAGRANGIAN_ANCHOR_INTERP(Anchor,y)
   real (kind=ipre), dimension(:), intent(in):: Anchor
   real (kind=ipre), dimension(:), intent(out):: y
   real (kind=ipre), dimension(2048):: temp_y
   real(kind=ipre), dimension(:), allocatable:: temp_Anchor
   real (kind=ipre):: x,L,xa_i,xa_j
   integer:: m,i,j,ia,ia1,ia2,Pix_Idx,Num_Anchor,Num_skip,Num_start
   integer:: discon

   m = 2    !2 means 5 point interp
   Num_Anchor = size(Anchor,1)  !should be 51
   allocate(Temp_Anchor(Num_Anchor))

! I'm leaving `Num_skip' and `Num_start' here as useful constants although they
! are not used in this subroutine.
   if (Num_Pix == 409) then
      ! It's a GAC scan line.
      Num_start = 5 !number of pixels between Anchors
      Num_skip = 8  !pixel number for first Anchor point
   else
      ! Must be a LAC scan line.
      Num_start = 25
      Num_skip = 40
   endif

!--- initialize
   y = 0.0
   temp_y = 0.0
   temp_Anchor = Anchor

!-------- check for a discontinuity (ie. crossing the 180ºW meridian)
    discon = 0
    if ((minval(temp_Anchor) < -140.0).and.(maxval(temp_Anchor) > 140.0)) then
     discon = 1
     where (temp_Anchor < 0.0)
        temp_Anchor =  temp_Anchor + 360.0
     endwhere
    endif

!-- perform 2m+1 point Lagrangian interpolation
      do Pix_Idx = 1, 2048
         ia = (Pix_Idx-25)/40+1
         ia1 =  min(Num_Anchor,max(1,ia-m))
         ia2 =  max(1,min(Num_Anchor,ia+m))
         x = Pix_Idx
         do i = ia1,ia2
           L = 1.0
           xa_i = (i-1)*(40)+25
           do j = ia1,ia2
            xa_j = (j-1)*(40)+25
            if (i /= j) L = L*(x-xa_j)/(xa_i-xa_j)
          enddo
          temp_y(Pix_Idx) = temp_y(Pix_Idx)+L*temp_Anchor(i)
         enddo
      enddo

! If it's a GAC scan line perform averaging of `temp_y' for the 2nd and 3rd LAC
! pixel going into every GAC pixel for the final value of `y'.
if (Num_Pix == 409) then
   do Pix_Idx = 1, 409
      y(Pix_Idx) = 0.5*(temp_y(5*Pix_Idx)+temp_y(5*Pix_Idx-1))
   enddo
else
   y = temp_y
endif

!----- account for any discontinuities
   where(y > 180.0)
    y = y - 360.0
   endwhere
   where(y < -180.0)
    y = 360.0 + y
   endwhere

   deallocate(Temp_Anchor)

 end subroutine LAGRANGIAN_ANCHOR_INTERP

!---------------------------------------------------------------------
! gnomic navigation Anchor interpolation routine
! note, this is only Valid for lat and lon, do not apply to angles
!---------------------------------------------------------------------
 subroutine GNOMIC_ANCHOR_INTERP(Lon_Anchor_deg,Lat_Anchor_deg,longitude,latitude)

   real (kind=ipre), dimension(:), intent(in):: Lon_Anchor_deg,Lat_Anchor_deg
   real (kind=ipre), dimension(:), intent(out):: longitude,latitude
   real (kind=ipre), dimension(:), allocatable:: k,x_Anchor,y_Anchor, &
                                                 Lon_Anchor,Lat_Anchor
   real (kind=ipre), dimension(Num_Pix):: x,y
   real (kind=ipre):: yy,xx,xa1,xa2,rho,v,v_prev,c,Lat_center,Lon_center,offset
   integer:: ia1,ia2,Pix_Idx,Num_Anchor,Num_skip,Num_start


   Num_Anchor = size(Lon_Anchor_deg)  !should be 51
   Num_skip = Num_Pix / Num_Anchor
   Num_start = (Num_Pix - Num_skip*(Num_Anchor-1))/2 + 1

   allocate(k(Num_Anchor))
   allocate(x_Anchor(Num_Anchor))
   allocate(y_Anchor(Num_Anchor))
   allocate(Lon_Anchor(Num_Anchor))
   allocate(Lat_Anchor(Num_Anchor))

!--- initialize
     latitude = 0.0
     longitude = 0.0

!--- convert Anchors to radians
     Lat_Anchor = Lat_Anchor_deg * dtor
     Lon_Anchor = Lon_Anchor_deg * dtor

!--- pick center for projection (in radians)
     Lat_center = Lat_Anchor(Num_Anchor/2)
     Lon_center = Lon_Anchor(Num_Anchor/2)

!--- convert to gnomic space
     k = sin(Lat_center) * sin(Lat_Anchor) +  &
         cos(Lat_center) * cos(Lat_Anchor) * cos(Lon_Anchor - Lon_center)
     k = 1.0 / k

     x_Anchor = k * cos(Lat_Anchor) * sin(Lon_Anchor - Lon_center)
     y_Anchor = k * (cos(Lat_center) * sin(Lat_Anchor) -  &
                     sin(Lat_center) * cos(Lat_Anchor) * cos(Lon_Anchor - Lon_center))

!--- linear interp in gnomic space
      do Pix_Idx = 1, Num_Pix

!--- determine two closest Anchor points
         ia1 = max(1,min(Num_Anchor-1,(Pix_Idx - Num_start)/Num_skip + 1))
         ia2 = ia1 + 1

!--- determine weights of two closest Anchor points
         xa1 = (ia1-1)*(Num_skip) + Num_start
         xa2 = (ia2-1)*(Num_skip) + Num_start

!--- determine point to interpolate to
         xx = Pix_Idx

!--- interp x at xx
         yy = (xx-xa1)*(x_Anchor(ia2)-x_Anchor(ia1))/(xa2-xa1) +  &
                   x_Anchor(ia1)
         x(Pix_Idx) = yy
         
!--- interp y at xx
         yy = (xx-xa1)*(y_Anchor(ia2)-y_Anchor(ia1))/(xa2-xa1) +  &
                   y_Anchor(ia1)
         y(Pix_Idx) = yy

      enddo

!--- convert other quantities
      v_prev = 0.0
      offset = 0
      do Pix_Idx = 1,Num_Pix

!--- some needed constants
        rho = sqrt(x(Pix_Idx)**2 + y(Pix_Idx)**2)

!--- handle center point where rho = 0.0
        if (rho == 0.0) then 
           latitude(Pix_Idx) = Lat_center
           longitude(Pix_Idx) = Lon_center
           cycle
        endif
        c = atan(rho)

!--- v is the longitude in radians relative to the center longitude
        v = atan( x(Pix_Idx) * sin(c) /  &
          (rho * cos(Lat_center) * cos(c) - y(Pix_Idx) * sin(Lat_center) * sin(c)) )

!--- look for a jump and determine offset
        if (Pix_Idx == 1) then
          v_prev = v
        endif

!       dv = v - v_prev
!       offset = 0.0
!       if (abs(dv) > 0.9*pi) then
!         if (dv < 0.0) then
!           offset = pi
!         else
!           offset = -1.0*pi
!         endif
!       endif

!       if (Pix_Idx == 205) then 
!          print *, dv, 0.9*pi, offset
!       endif
!       v = v + offset

        latitude(Pix_Idx) = asin(cos(c) * sin(Lat_center) + (y(Pix_Idx)*sin(c)*cos(Lat_center) / rho))
        longitude(Pix_Idx) = Lon_center + v

        v_prev = v

      enddo

!--- convert to degrees
      latitude = latitude / dtor
      longitude = longitude / dtor
     

!----- account for discontinuity if there was one
  where (longitude > 180.0)
    longitude = longitude - 360.0
  endwhere
  where (longitude < -180.0)
    longitude = longitude + 360.0
  endwhere


  deallocate(k)
  deallocate(x_Anchor)
  deallocate(y_Anchor)
  deallocate(Lon_Anchor)
  deallocate(Lat_Anchor)

  end subroutine GNOMIC_ANCHOR_INTERP

!---------------------------------------------------------------------
! linear interpolation routine
!---------------------------------------------------------------------
 subroutine LINEAR_ANCHOR_INTERP(Anchor,y)
   real (kind=ipre), dimension(:), intent(in):: Anchor
   real (kind=ipre), dimension(:), intent(out):: y
   real (kind=ipre), dimension(:), allocatable:: temp_Anchor
   real (kind=ipre):: x,xa1,xa2
   integer:: ia1,ia2,Pix_Idx,Num_Anchor,Num_skip,Num_start
   integer:: discon

   Num_Anchor = size(Anchor)  !should be 51
   Num_skip = Num_Pix / Num_Anchor
   Num_start = (Num_Pix - Num_skip*(Num_Anchor-1))/2 + 1

   allocate(Temp_Anchor(Num_Anchor))

!--- initialize
     y = 0.0
     temp_Anchor = Anchor

!-------- check for a discontinuity (ie. crossing dateline)
    discon = 0
    if ((minval(temp_Anchor) < -175.0).and.(maxval(temp_Anchor) > 175.0)) then
     discon = 1
     where( temp_Anchor < 0.0)
        temp_Anchor =  temp_Anchor + 360.0
     endwhere
    endif

!-----------------------------------------------------------
      do Pix_Idx = 1, Num_Pix
         ia1 = max(1,min(Num_Anchor-1,(Pix_Idx - Num_start)/Num_skip + 1))
         ia2 = ia1 + 1

         xa1 = (ia1-1)*(Num_skip) + Num_start
         xa2 = (ia2-1)*(Num_skip) + Num_start

         x = Pix_Idx

         y(Pix_Idx) = (x-xa1)*(temp_Anchor(ia2)-temp_Anchor(ia1))/(xa2-xa1) +  &
                   temp_Anchor(ia1)

      enddo

!----- account for discontinuity if there was one
  if (discon == 1) then
   where (y > 180.0)
    y = y - 360.0
   endwhere
  endif

  deallocate(Temp_Anchor)

  end subroutine LINEAR_ANCHOR_INTERP

!-------------------------------------------------------------------------------
! COMPUTE_ANGLE_ANCHORS
!
! input
!   jday_Temp -  Julian Day
!   Time_Temp - Time_Temp of day_Temp
!   Lon_Temp - vector of longitude (degrees)
!   Lat_Temp - vector of latitude (degrees)
!
! output
!   Satzen_Temp - vector of satellite view zenith angles (degrees)
!   Solzen_Temp - vector of solar zenith angles (degrees)
!   Relaz_Temp - vector of relative azimuth angles
!   Solaz_Temp - vector of solar azimuth angles
! 
! Revision History
!   Jan 2003 - Coded A. Heidinger (Based on CLAVR-1 code)   
!
! note - Relaz is the relative azimuth defined as this
!   Relaz = 180 means the sun is behind you (backscatter)
!   Relaz = 0 mean the sun is infront of you (forward scatter)
!   this is convention used in noaa-AVHRR_KLM_Flag Relaz but not in clavr-1
!-------------------------------------------------------------------------------
   subroutine COMPUTE_ANGLE_ANCHORS(jday_Temp,   &
                                    Time_Temp,   &
                                    Lon_Temp,    &
                                    Lat_Temp,    &
                                    Satzen_Temp, &
                                    Solzen_Temp, &
                                    Relaz_Temp,  &
                                    Solaz_Temp,  &
                                    Sataz_Temp,  &
                                    Glint_Temp,  &
                                    Scat_Temp)


     integer(kind=int2), intent(in):: jday_Temp
     integer(kind=int4), intent(in):: Time_Temp
     real(kind=ipre), dimension(:), intent(in):: Lon_Temp
     real(kind=ipre), dimension(:), intent(in):: Lat_Temp
     real(kind=ipre), dimension(:), intent(out):: Satzen_Temp
     real(kind=ipre), dimension(:), intent(out):: Solzen_Temp
     real(kind=ipre), dimension(:), intent(out):: Relaz_Temp
     real(kind=ipre), dimension(:), intent(out):: Solaz_Temp
     real(kind=ipre), dimension(:), intent(out):: Sataz_Temp
     real(kind=ipre), dimension(:), intent(out):: Glint_Temp
     real(kind=ipre), dimension(:), intent(out):: Scat_Temp
     integer:: i
     integer:: n
     integer:: jday_Temp4
     real(kind=real4):: Frac_Time_Temp
     real(kind=real4):: Lat_Temp_Center
     real(kind=real4):: Lon_Temp_Center
     real(kind=real4):: Solzen_Temp_Center

!    real(kind=real8):: Cosgeo
!    real(kind=real8):: Psix
!    real(kind=real8):: Geox
!    real(kind=real8):: Numor
!    real(kind=real8):: Denom
!    real(kind=real8):: Cos_Solzen_Temp

     real(kind=real4):: Geox


      !--- convert time to fraction hours
      Frac_Time_Temp = Time_Temp / (60.0*60.0*1000.0)

      !--- determine size of the arrays to process
      n = size(Lon_Temp)

      !--- make a local int4 julian day copy
      Jday_Temp4 = Jday_Temp

     !-------------------------------------------------
     ! solar zenith angle
     !-------------------------------------------------
      do i = 1, n
       call POSSOL(Jday_Temp4, &
                   Frac_Time_Temp, &
                   Lon_Temp(i), &
                   Lat_Temp(i), &
                   Solzen_Temp(i), &
                   Solaz_Temp(i))
      enddo

      !-------------------------------------------------
      !---- relative azimuth and satellite zenith angle
      !-------------------------------------------------
      !-- sub-satellite point
      Lat_Temp_Center = Lat_Temp(26)
      Lon_Temp_Center = Lon_Temp(26)
      Solzen_Temp_Center = Solzen_Temp(26)

      do i = 1,n

       !-------------------------------------
       ! original CLAVR-1 Method
       !-------------------------------------

       Geox = great_circle_angle(Lon_Temp_Center,Lat_Temp_Center,Lon_Temp(i),Lat_Temp(i))
       Satzen_Temp(i) = sensor_zenith_avhrr_anchor(Geox,i)
       Relaz_Temp(i) = relative_azimuth_avhrr_anchor(Geox,Solzen_Temp_Center,Solzen_Temp(i))

       !-------------------------------------
       ! New Generic Method is not generating sate
       !-------------------------------------
!-->   Satzen_Temp(i) = SENSOR_ZENITH(AVHRR_ALTITUDE,Lon_Temp_Center,Lat_Temp_Center,Lon_Temp(i),Lat_Temp(i)) 

       Sataz_Temp(i) = SENSOR_AZIMUTH(Lon_Temp_Center,Lat_Temp_Center,Lon_Temp(i),Lat_Temp(i)) 

!-->   Relaz_Temp(i) = Relative_AZIMUTH(Solaz_Temp(i), Sataz_Temp(i)) 
!      print *, 'New = ', Satzen_Temp(i), Relaz_Temp(i)

       Glint_Temp(i) = GLINT_ANGLE(Solzen_Temp(i),Satzen_Temp(i),Relaz_Temp(i)) 

       Scat_Temp(i) = SCATTERING_ANGLE(Solzen_Temp(i),Satzen_Temp(i),Relaz_Temp(i)) 

      enddo

 end subroutine COMPUTE_ANGLE_ANCHORS

!-------------------------------------------------------------------------------
! to use for AVHRR_GAC_Flag data:
!  1. set AVHRR_GAC_Flag to true
!  2. read the AVHRR 1B data records dIrectly into Buffer
!
! to use for lac data:
!  1. set AVHRR_GAC_Flag to false (0)
!  2. combine the two  data-records from the AVHRR 1B file  into Buffer
!-------------------------------------------------------------------------------
 subroutine UNPACK_AVHRR_DATA_RECORD(Line_Idx)

   integer, intent(in):: Line_Idx
   integer(kind=int4):: i4word
   integer(kind=int2):: i2word
   integer(kind=int4):: Num_data_bytes
   integer(kind=int4):: Num_Chan
   integer:: Chan_Idx
   integer:: Pix_Idx
   integer:: j
   integer:: Byte_Idx
   real, dimension(5):: slope
   real, dimension(5):: intercept
   integer(kind=int2), dimension(2):: data_byte

   !-------------------------------------------------------------------------------
   ! set parameters for AVHRR_GAC_Flag or lac processing
   !-------------------------------------------------------------------------------
   if (AVHRR_GAC_Flag == sym%YES) then
      Num_data_bytes = 2728
    else
      Num_data_bytes = 6952 + 6704
   endif
   Num_Chan = size(Chan_Counts_Avhrr,1)

   !----------------------------------------------------------------
   ! unpack scan line number
   !----------------------------------------------------------------
   Scan_Number(Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(1:2),sym%UNSIGNED,Byte_Swap_1b)

   !----------------------------------------------------------------------
   ! unpack Time_Temp code
   !----------------------------------------------------------------------
   data_byte(1) = Buffer_AVHRR(3)
   data_byte(2) = Buffer_AVHRR(4)
   where (data_byte < 0)
         data_byte = data_byte + 256
   end where
   i4word =  data_byte(1)*256 + data_byte(2)
   Scan_Year(Line_Idx) = i4word / (2**9)
   Scan_Day(Line_Idx) = i4word - Scan_Year(Line_Idx) * (2**9)

   i4word = MAKE_I4WORD(Buffer_AVHRR(5:8),sym%UNSIGNED,Byte_Swap_1b)
   Scan_Time(Line_Idx) = i4word - ((i4word/2**27)*(2**27))
   if (Scan_Time(Line_Idx) > Biggest_I4) Scan_Time(Line_Idx) = Missing_Value_Int4
   if (Scan_Time(Line_Idx) < -1*Biggest_I4) Scan_Time(Line_Idx) = Missing_Value_Int4

   !----------------------------------------------------------------------
   ! unpack quality indicators - not done fully
   !----------------------------------------------------------------------
   i2word = Buffer_AVHRR(9)
   if (i2word < 0) then
      i2word = i2word + 256
   end if
   Fatal_AVHRR(Line_Idx) = ishft(i2word, -7)

   i2word = Buffer_AVHRR(9)
   if (i2word < 0) then
      i2word = i2word + 256
   end if
   Ascend(Line_Idx) = ishft(ishft(i2word, 6),-7)

   !--- this needs to be verified
   i2word = Buffer_AVHRR(10)
   if (i2word < 0) then
      i2word = i2word + 256
   end if
   Resync_AVHRR(Line_Idx) = ishft(ishft(i2word, 0),-7)
   
!----------------------------------------------------------------------
! unpack calibration coefficients
!----------------------------------------------------------------------
    Chan_Idx = 1
    do Byte_Idx = 13,52,8
     slope(Chan_Idx) = MAKE_I4WORD(Buffer_AVHRR(Byte_Idx:Byte_Idx+3),sym%SIGNED,Byte_Swap_1b) / (2.0**30)
     intercept(Chan_Idx) = MAKE_I4WORD(Buffer_AVHRR(Byte_Idx+4:Byte_Idx+7),sym%SIGNED,Byte_Swap_1b) / (2.0**22)
     Chan_Idx = Chan_Idx + 1
    enddo

!----------------------------------------------
! mimic noaa AVHRR_KLM_Flag type calibration coefficients
!----------------------------------------------
    vis_slope_1(1,Line_Idx) = slope(1)   !high gain is same as low gain
    vis_slope_1(2,Line_Idx) = slope(2)   !high gain is same as low gain
    vis_slope_1(3,Line_Idx) = 0.0        !no channel 3a
    vis_slope_2(:,Line_Idx) = vis_slope_1(:,Line_Idx)   !high gain is same as low gain
    vis_intercept_2(:,Line_Idx) = vis_intercept_1(:,Line_Idx)   !high gain intercept same as low
    vis_intersection(:,Line_Idx) = 1024     !these values are ignored

    ir_linear_slope_1b(3,Line_Idx) = slope(3)
    ir_linear_slope_1b(4,Line_Idx) = slope(4)
    ir_linear_slope_1b(5,Line_Idx) = slope(5)

    ir_linear_intercept_1b(3,Line_Idx) = intercept(3)
    ir_linear_intercept_1b(4,Line_Idx) = intercept(4)
    ir_linear_intercept_1b(5,Line_Idx) = intercept(5)

    !----------------------------------------------------------------------
    ! unpack geolocation coefficients
    !----------------------------------------------------------------------
    Num_loc = Buffer_AVHRR(53)

    !----------------------------------------------------------------------
    ! unpack solar zenith angle data - treating as unsigned integer
    !----------------------------------------------------------------------
    do Byte_Idx = 54,104
        j = Byte_Idx - 53
        if (Buffer_AVHRR(Byte_Idx) < 0) then
             Solzen_Anchor(j,Line_Idx) = (256.0 + Buffer_AVHRR(Byte_Idx)) / 2.0
        else
             Solzen_Anchor(j,Line_Idx) = Buffer_AVHRR(Byte_Idx) / 2.0
        endif
    enddo

    !----------------------------------------------------------------------
    ! unpack earth location
    !----------------------------------------------------------------------
    j = 1
    do Byte_Idx = 105,308,4
        Lat_Anchor_1b(j,Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(Byte_Idx:Byte_Idx+1),sym%SIGNED,Byte_Swap_1b) / 128.0
        Lon_Anchor_1b(j,Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(Byte_Idx+2:Byte_Idx+3),sym%SIGNED,Byte_Swap_1b) / 128.0
        j = j + 1
    enddo

    !-- make sure longitudes are correct ( ?? check this)
    where (Lon_Anchor_1b(:,Line_Idx) > 180)
        Lon_Anchor_1b(:,Line_Idx) = Lon_Anchor_1b(:,Line_Idx) - 360.0
    end where

    !----------------------------------------------------------------------
    ! unpack channel  data
    !----------------------------------------------------------------------
    Chan_Idx = 1
    Pix_Idx = 1

    do Byte_Idx = 449,448+Num_data_bytes,4

      i4word = MAKE_I4WORD(Buffer_AVHRR(Byte_Idx:Byte_Idx+3),sym%UNSIGNED,Byte_Swap_1b)

      do j = 1,3

        if (Pix_Idx > Num_Pix) then
            exit
        endif

        Chan_Counts_Avhrr(Chan_Idx,Pix_Idx,Line_Idx) = ishft(ishft(i4word,2+10*(j-1)),-22)

        Chan_Idx = Chan_Idx + 1

        if (Chan_Idx > Num_Chan) then
          Chan_Idx = 1
          Pix_Idx = Pix_Idx + 1
        endif

      enddo

      if (Pix_Idx > Num_Pix) then
            exit
      endif

    enddo

  end subroutine UNPACK_AVHRR_DATA_RECORD

!-------------------------------------------------------------------------------
! to use for AVHRR_GAC_Flag data:
!  1. set AVHRR_GAC_Flag to true (1)
!  2. read the AVHRR 1B data records dIrectly into Buffer
!
! to use for lac data:
!  1. set AVHRR_GAC_Flag to false (0)
!  2. combine the two  data-records from the AVHRR 1B file  into Buffer
!
!
! Modifications
!  January 2005 - Modified for NOAA-N
!
!  Modified for one scan-line at a time
!-------------------------------------------------------------------------------
 subroutine UNPACK_AVHRR_DATA_RECORD_KLM(Line_Idx)

   integer, intent(in):: Line_Idx
   integer(kind=int4):: i4word
   integer(kind=int2):: i2word
   integer(kind=int4):: Num_data_bytes
   integer(kind=int4):: Num_Chan
   integer(kind=int4):: Num_clavr_bytes
   integer(kind=int4):: Num_digtel_bytes
   integer(kind=int4):: Num_house_bytes
   integer(kind=int4):: Start_Pixel
   integer(kind=int4):: End_Pixel
   integer:: Chan_Idx
   integer:: Pix_Idx
   integer:: i
   integer:: j
   integer:: Byte_Idx
   integer:: iword
   integer:: Start_word
   integer:: End_word
   integer:: word_clavr_start
   integer(kind=int1):: onebyte

!-------------------------------------------------------------------------------
! set parameters for AVHRR_GAC_Flag or lac processing
!-------------------------------------------------------------------------------
   if (AVHRR_GAC_Flag == sym%YES) then
      Num_data_bytes = 4000 - 1265 + 1
      Num_digtel_bytes =  4016 - 4001 + 1
      Num_house_bytes =  4048 - 4017 + 1
      Num_clavr_bytes = 4160 - 4049 + 1
    else
      Num_data_bytes = 14928 - 1265
      Num_digtel_bytes = 14944 - 14929
      Num_house_bytes = 14976 - 14971
      Num_clavr_bytes = 15496 - 14977
   endif
   Num_Chan = size(Chan_Counts_Avhrr,1)

   Valid_Therm_Cal = .true.

   !----------------------------------------------------------------
   ! unpack scan line number
   !----------------------------------------------------------------
   Scan_Number(Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(1:2),sym%UNSIGNED,Byte_Swap_1b)

   !----------------------------------------------------------------------
   ! unpack Time_Temp code
   !----------------------------------------------------------------------
   Scan_Year(Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(3:4),sym%UNSIGNED,Byte_Swap_1b)
   Scan_Day(Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(5:6),sym%UNSIGNED,Byte_Swap_1b)
   Scan_Time(Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(9:12),sym%UNSIGNED,Byte_Swap_1b)
   if (Scan_Time(Line_Idx) > Biggest_I4) Scan_Time(Line_Idx) = Missing_Value_Int4
   if (Scan_Time(Line_Idx) < -1*Biggest_I4) Scan_Time(Line_Idx) = Missing_Value_Int4

   !----------------------------------------------------------------------
   ! unpack quality indicators - not done fully
   !----------------------------------------------------------------------
   i2word  = Buffer_AVHRR(25)
   if (i2word < 0) then
      i2word = i2word + 256
   end if
   Fatal_AVHRR(Line_Idx) = ishft(i2word, -7)

   i2word  = Buffer_AVHRR(13)
   if (i2word < 0) then
      i2word = i2word + 256
   end if
   Ascend(Line_Idx) = ishft(i2word, -7)

   i2word  = Buffer_AVHRR(14)
   if (i2word < 0) then
      i2word = i2word + 256
   end if
   Ch3a_On_AVHRR(Line_Idx) = ishft(ishft(i2word, 6),-6)
   Ch6_On_Pixel_Mask(:,Line_Idx) = Ch3a_On_AVHRR(Line_Idx)

   i2word  = Buffer_AVHRR(28)
   if (i2word < 0) then
      i2word = i2word + 256
   end if
   Resync_AVHRR(Line_Idx) = ishft(ishft(i2word, 6),-7)

!---------------------------------------------------------------------
! unpack space view
!---------------------------------------------------------------------
  do j = 1,10
   Space_Count_Temp(1,j) = MAKE_I2WORD(Buffer_AVHRR(1161+10*(j-1):1162+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
   Space_Count_Temp(2,j) = MAKE_I2WORD(Buffer_AVHRR(1163+10*(j-1):1164+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
   Space_Count_Temp(3,j) = MAKE_I2WORD(Buffer_AVHRR(1165+10*(j-1):1166+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
   Space_Count_Temp(4,j) = MAKE_I2WORD(Buffer_AVHRR(1167+10*(j-1):1168+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
   Space_Count_Temp(5,j) = MAKE_I2WORD(Buffer_AVHRR(1169+10*(j-1):1170+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
  enddo
  
!----------------------------------------------------------------------
! unpack calibration coefficients
!----------------------------------------------------------------------
  do Chan_Idx = 1, 3

    vis_slope_1(Chan_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(49+(Chan_Idx-1)*60:52+(Chan_Idx-1)*60), &
                                     sym%SIGNED,Byte_Swap_1b) / (10.0**7)

    vis_intercept_1(Chan_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(53+(Chan_Idx-1)*60:56+(Chan_Idx-1)*60), &
                                         sym%SIGNED,Byte_Swap_1b) / (10.0**6)

    vis_slope_2(Chan_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(57+(Chan_Idx-1)*60:60+(Chan_Idx-1)*60), &
                                     sym%SIGNED,Byte_Swap_1b) / (10.0**7)

    vis_intercept_2(Chan_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(61+(Chan_Idx-1)*60:64+(Chan_Idx-1)*60), &
                                     sym%SIGNED,Byte_Swap_1b) / (10.0**6)

    vis_intersection(Chan_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(65+(Chan_Idx-1)*60:68+(Chan_Idx-1)*60), &
                                             sym%SIGNED,Byte_Swap_1b)

 enddo

 do Chan_Idx = 1,3
    ir_coef_1_1b(Chan_Idx+2,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(229+(Chan_Idx-1)*24:232+(Chan_Idx-1)*24), &
                                              sym%SIGNED,Byte_Swap_1b)/(10.0**6)
    ir_coef_2_1b(Chan_Idx+2,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(233+(Chan_Idx-1)*24:236+(Chan_Idx-1)*24), &
                                              sym%SIGNED,Byte_Swap_1b) / (10.0**6)
    ir_coef_3_1b(Chan_Idx+2,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(237+(Chan_Idx-1)*24:240+(Chan_Idx-1)*24), &
                                              sym%SIGNED,Byte_Swap_1b) / (10.0**6)
 enddo

!-----------------------------------------------------------------------------
! handling NOAA-N change to scaling of ir_coef_3_1b for ch4 and ch5 
! the scaling changed from 10**6 to 10**7 - account for extra
! factor of 10 in NOAA-N format here
!------------------------------------------------------------------------------
  if (ver_1b >= 3) then
    ir_coef_3_1b(4,Line_Idx) = ir_coef_3_1b(4,Line_Idx) / 10.0
    ir_coef_3_1b(5,Line_Idx) = ir_coef_3_1b(5,Line_Idx) / 10.0
  endif

!----------------------------------------------------------------------
! unpack number of Anchor points
!----------------------------------------------------------------------
!   Num_loc = Buffer_Temp(53)
    Num_loc = 51

!----------------------------------------------------------------------
! unpack angle data
!----------------------------------------------------------------------
      do j = 1, 51
        Byte_Idx = 329 + (j-1)*6
        Solzen_Anchor(j,Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(Byte_Idx:Byte_Idx+1),sym%SIGNED,Byte_Swap_1b)/(10.0**2)
        Satzen_Anchor(j,Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(Byte_Idx+2:Byte_Idx+3),sym%SIGNED,Byte_Swap_1b)/(10.0**2)
        Relaz_Anchor(j,Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(Byte_Idx+4:Byte_Idx+5),sym%SIGNED,Byte_Swap_1b)/(10.0**2)
      enddo

!----------------------------------------------------------------------
! unpack earth location
!----------------------------------------------------------------------
      j = 1
      do j = 1,51
        Byte_Idx = 641 + (j-1)*8 
        Lat_Anchor_1b(j,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(Byte_Idx:Byte_Idx+3),sym%SIGNED,Byte_Swap_1b)/(10.0**4)
        Lon_Anchor_1b(j,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(Byte_Idx+4:Byte_Idx+7),sym%SIGNED,Byte_Swap_1b)/(10.0**4)
      enddo

!-- make sure longitudes are correct ( ?? check this)
      where (Lon_Anchor_1b(:,Line_Idx) > 180)
         Lon_Anchor_1b(:,Line_Idx) = Lon_Anchor_1b(:,Line_Idx) - 360.0
      end where

!----------------------------------------------------------------------
! unpack channel  data
!----------------------------------------------------------------------
   Chan_Idx = 1
   Pix_Idx = 1

   do Byte_Idx = 1265,1264+Num_data_bytes,4

      i4word = MAKE_I4WORD(Buffer_AVHRR(Byte_Idx:Byte_Idx+3),sym%UNSIGNED,Byte_Swap_1b)

      do j = 1,3
       if (Pix_Idx > Num_Pix) then
           exit
       endif

!--- old
!      Chan_Counts_Avhrr(Chan_Idx,Pix_Idx,Line_Idx) = i4word/ (1024**(3-j))
!      i4word = i4word - Chan_Counts_Avhrr(Chan_Idx,Pix_Idx,Line_Idx)*(1024**(3-j))
!--- new
       Chan_Counts_Avhrr(Chan_Idx,Pix_Idx,Line_Idx) = ishft(ishft(i4word,2+10*(j-1)),-22)

       Chan_Idx = Chan_Idx + 1
        if (Chan_Idx > Num_Chan) then
          Chan_Idx = 1
          Pix_Idx = Pix_Idx + 1
        endif
      enddo

      if (Pix_Idx > Num_Pix) then
            exit
      endif

   enddo

    !-------------------------------------------------------------------------
    ! read in cloud mask 
    !-----------------------------------------------------------------------
    Pix_Idx = 0
    word_clavr_start = 1264 + Num_data_bytes + Num_digtel_bytes + Num_house_bytes + 1
 
    !------ clavr status flag
    onebyte = Buffer_AVHRR(word_clavr_start + 3)
    Clavr_Status_1b_AVHRR(Line_Idx) = ishft(ishft(onebyte,7),-7)

    !--- cloud mask
    Start_word = word_clavr_start + 8
    End_word = Start_word + Num_clavr_bytes

    do iword = Start_word, End_word
      onebyte = Buffer_AVHRR(iword)
      Pix_Idx = 1
      Start_Pixel = (iword - Start_word)*4 + 1
      End_Pixel = min(Num_Pix, Start_Pixel+3_int2)

      do Pix_Idx = Start_Pixel, End_Pixel
         if (Pix_Idx > Num_Pix) then
           exit
         endif
         i = (Pix_Idx - Start_Pixel) + 1
         Cld_Mask_Aux(Pix_Idx,Line_Idx) = ishft(ishft(onebyte,2*(i-1)),-6)
      enddo
    enddo

    !--- set flag to communicate that mask was read in
    Cloud_Mask_Aux_Read_Flag = sym%YES

  end subroutine UNPACK_AVHRR_DATA_RECORD_KLM

!----------------------------------------------------------------------------------------
!  Read data records for AAPP
!----------------------------------------------------------------------------------------
 subroutine UNPACK_AVHRR_DATA_RECORD_AAPP(Line_Idx)

   integer, intent(in):: Line_Idx
   integer(kind=int2):: i2word
   integer(kind=int4):: Num_Chan
   integer(kind=int4):: Num_clavr_bytes
   integer(kind=int4):: Num_digtel_bytes
   integer(kind=int4):: Num_house_bytes
   integer:: Chan_Idx
   integer:: Pix_Idx
   integer:: Anchor_Idx
   integer:: j
   integer:: Byte_Idx
   integer(kind=int4):: first_data_byte
   integer(kind=int1):: onebyte

   !-------------------------------------------------------------------------------
   ! set parameters for AVHRR_GAC_Flag or lac processing
   !-------------------------------------------------------------------------------
   if (AVHRR_GAC_Flag == sym%YES) then
      print *, EXE_PROMPT, "ERROR: Cannot handle GAC AAPP, stopping"
      stop
   endif

   first_data_byte = 1265
   Num_digtel_bytes = 14944 - 14929
   Num_house_bytes = 14976 - 14971
   Num_clavr_bytes = 0

   Num_Chan = size(Chan_Counts_Avhrr,1)

   Valid_Therm_Cal = .true.

   !----------------------------------------------------------------
   ! unpack scan line number
   !----------------------------------------------------------------
   Scan_Number(Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(1:2),sym%UNSIGNED,Byte_Swap_1b)
 
   !----------------------------------------------------------------------
   ! unpack Time_Temp code
   !----------------------------------------------------------------------
   Scan_Year(Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(3:4),sym%UNSIGNED,Byte_Swap_1b)
   Scan_Day(Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(5:6),sym%UNSIGNED,Byte_Swap_1b)
   Scan_Time(Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(9:12),sym%UNSIGNED,Byte_Swap_1b)
   if (Scan_Time(Line_Idx) > Biggest_I4) Scan_Time(Line_Idx) = Missing_Value_Int4
   if (Scan_Time(Line_Idx) < -1*Biggest_I4) Scan_Time(Line_Idx) = Missing_Value_Int4

   !----------------------------------------------------------------------
   ! unpack quality indicators - not done fully
   !----------------------------------------------------------------------
   i2word  = Buffer_AVHRR(25)
   if (i2word < 0) then
      i2word = i2word + 256
   end if
   Fatal_AVHRR(Line_Idx) = ishft(i2word, -7)

   !--- Ch3a/b flag - note difference in position nesdis (due to byte swapping)
   i2word  = Buffer_AVHRR(14)
   if (i2word < 0) then
      i2word = i2word + 256
   end if
   Ascend(Line_Idx) = ishft(i2word, -7)

   !--- Ch3a/b flag - note difference in position and meaning from nesdis
   i2word  = Buffer_AVHRR(13)
   if (i2word < 0) then
      i2word = i2word + 256
   end if
   onebyte = ishft(ishft(i2word, 6),-6)

   Ch3a_On_AVHRR(Line_Idx) = onebyte
   if (onebyte == 0) then 
      Ch3a_On_AVHRR(Line_Idx) = 1
   endif
   if (onebyte == 1) then 
      Ch3a_On_AVHRR(Line_Idx) = 0
   endif
   Ch6_On_Pixel_Mask(:,Line_Idx) = Ch3a_On_AVHRR(Line_Idx)

   i2word  = Buffer_AVHRR(28)
   if (i2word < 0) then
      i2word = i2word + 256
   end if
   Resync_AVHRR(Line_Idx) = ishft(ishft(i2word, 6),-7)

   !---------------------------------------------------------------------
   ! unpack space view
   !---------------------------------------------------------------------
   do j = 1,10
    Space_Count_Temp(1,j) = MAKE_I2WORD(Buffer_AVHRR(1161+10*(j-1):1162+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
    Space_Count_Temp(2,j) = MAKE_I2WORD(Buffer_AVHRR(1163+10*(j-1):1164+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
    Space_Count_Temp(3,j) = MAKE_I2WORD(Buffer_AVHRR(1165+10*(j-1):1166+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
    Space_Count_Temp(4,j) = MAKE_I2WORD(Buffer_AVHRR(1167+10*(j-1):1168+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
    Space_Count_Temp(5,j) = MAKE_I2WORD(Buffer_AVHRR(1169+10*(j-1):1170+10*(j-1)),sym%SIGNED,Byte_Swap_1b)
  enddo

  !----------------------------------------------------------------------
  ! unpack calibration coefficients
  !----------------------------------------------------------------------
  !--- note switch in scale factors from nesdis
  !--- in AAPP, the operational and test coefficients are always zero, so
  !--- here we read in the prelaunch

  do Chan_Idx = 1, 3
    vis_slope_1(Chan_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(89+(Chan_Idx-1)*60:92+(Chan_Idx-1)*60), &
                                           sym%SIGNED,Byte_Swap_1b) / (10.0**10)
    vis_intercept_1(Chan_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(93+(Chan_Idx-1)*60:96+(Chan_Idx-1)*60), & 
                                               sym%SIGNED,Byte_Swap_1b) / (10.0**7)
    vis_slope_2(Chan_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(97+(Chan_Idx-1)*60:100+(Chan_Idx-1)*60), &
                                           sym%SIGNED,Byte_Swap_1b) / (10.0**10)
    vis_intercept_2(Chan_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(101+(Chan_Idx-1)*60:104+(Chan_Idx-1)*60), &
                                               sym%SIGNED,Byte_Swap_1b) / (10.0**7)
    vis_intersection(Chan_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(105+(Chan_Idx-1)*60:108+(Chan_Idx-1)*60), &
                                                sym%SIGNED,Byte_Swap_1b)
  enddo

  !--- note switch from nesdis in the order
  do Chan_Idx = 1,3
    ir_coef_3_1b(Chan_Idx+2,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(229+(Chan_Idx-1)*24: &
                          232+(Chan_Idx-1)*24),sym%SIGNED,Byte_Swap_1b)/(10.0**9)
    ir_coef_2_1b(Chan_Idx+2,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(233+(Chan_Idx-1)*24: &
                          236+(Chan_Idx-1)*24),sym%SIGNED,Byte_Swap_1b) / (10.0**6)
    ir_coef_1_1b(Chan_Idx+2,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(237+(Chan_Idx-1)*24: &
                          240+(Chan_Idx-1)*24),sym%SIGNED,Byte_Swap_1b) / (10.0**6)
  enddo

  !----------------------------------------------------------------------
  ! unpack number of Anchor points
  !----------------------------------------------------------------------
  !Num_Loc = Buffer_AVHRR(53)
  Num_Loc = 51

  !----------------------------------------------------------------------
  ! unpack angle data
  !----------------------------------------------------------------------
  do Anchor_Idx = 1, Num_Loc
        Byte_Idx = 329 + (Anchor_Idx-1)*6
        Solzen_Anchor(Anchor_Idx,Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(Byte_Idx:Byte_Idx+1),sym%SIGNED,Byte_Swap_1b)/(10.0**2)
        Satzen_Anchor(Anchor_Idx,Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(Byte_Idx+2:Byte_Idx+3),sym%SIGNED,Byte_Swap_1b)/(10.0**2)
        Relaz_Anchor(Anchor_Idx,Line_Idx) = MAKE_I2WORD(Buffer_AVHRR(Byte_Idx+4:Byte_Idx+5),sym%SIGNED,Byte_Swap_1b)/(10.0**2)
  enddo

  !----------------------------------------------------------------------
  ! unpack earth location
  !----------------------------------------------------------------------
  Anchor_Idx = 1
  do Anchor_Idx = 1,Num_Loc
        Byte_Idx = 641 + (Anchor_Idx-1)*8
        Lat_Anchor_1b(Anchor_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(Byte_Idx:Byte_Idx+3),sym%SIGNED,Byte_Swap_1b)/(10.0**4)
        Lon_Anchor_1b(Anchor_Idx,Line_Idx) = MAKE_I4WORD(Buffer_AVHRR(Byte_Idx+4:Byte_Idx+7),sym%SIGNED,Byte_Swap_1b)/(10.0**4)
  enddo

  !----------------------------------------------------------------------
  ! unpack channel  data
  !----------------------------------------------------------------------

  !--- AAPP has each observation stored as two byte integers - not packed
  First_Data_Byte = 1265

  do Pix_Idx = 1, 2048
    do Chan_Idx = 1,5
       Byte_Idx = First_Data_Byte + (Pix_Idx-1)*10 + (Chan_Idx-1)*2
       Chan_Counts_Avhrr(Chan_Idx,Pix_Idx,Line_Idx) =  &
                         MAKE_I2WORD(Buffer_AVHRR(Byte_Idx:Byte_Idx+1),sym%UNSIGNED,Byte_Swap_1b) 
    enddo
  enddo

end subroutine UNPACK_AVHRR_DATA_RECORD_AAPP

!---------------------------------------------------------
! compute the ascending/descending node flag 
!---------------------------------------------------------
  subroutine CALCULATE_ASC_DES(Line_Idx_Min_Segment,Num_Scans_Read)

        integer, intent(in):: Line_Idx_Min_Segment
        integer, intent(in):: Num_Scans_Read
        real:: diff
        real:: diff_before
        real:: diff_after
        integer:: j

        do j = Line_Idx_Min_Segment+1, Num_Scans_Read - 2
        
         if (Bad_Scan_Flag(j) == sym%YES) then

            if (j > Line_Idx_Min_Segment + 1) then
              Ascend(j) = Ascend(j-1)
            endif

            cycle ! Cycle if there is missing navigation

         endif

         diff = Lat_Anchor_1b(25,j) - Lat_Anchor_1b(25,j+1)

         Ascend(j) = 1
         if (diff <= 0) Ascend(j) = 0

         if (Bad_Scan_Flag(j+1) == sym%YES .and.  &
              Bad_Scan_Flag(j-1) == sym%NO) then
              Ascend(j) = Ascend(j-1)
         endif

         diff_before = Lat_Anchor_1b(25,j-1) - Lat_Anchor_1b(25,j)
         diff_after = Lat_Anchor_1b(25,j+1) - Lat_Anchor_1b(25,j+2)

         ! check to see if current one is a fluke in diff
         if (diff_after >= 0 .and. Ascend(j-1) == 1 .and.  &
              diff <= 0 .and. Ascend(j) == 0)  then
              Ascend(j) = 1
         endif
         if (Ascend(j-1) == 0 .and. diff_after <= 0 .and.  &
             diff >= 0 .and. Ascend(j) == 1)  then
             Ascend(j) = 0
         endif

         enddo

        Ascend(Line_Idx_Min_Segment) = Ascend(Line_Idx_Min_Segment+1)
        Ascend(Num_Scans_Read-1) = Ascend(Num_Scans_Read - 2)
        Ascend(Num_Scans_Read) = Ascend(Num_Scans_Read - 2)

  end subroutine CALCULATE_ASC_DES
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_AVHRR_ARRAYS(dim1,dim2,dim3)
   integer, intent(in):: dim1, dim2, dim3
   allocate(Buffer_AVHRR(dim3))
   allocate(Segment_Buffer_AVHRR(dim2*dim3))
   allocate(Fatal_AVHRR(dim2))
   allocate(Resync_AVHRR(dim2))
   allocate(Clavr_Status_1b_AVHRR(dim2))
   !--- Note Header_BUffer_AVHRR is allocated in avhrr_module.f90
   allocate(Chan_Counts_Avhrr(5,dim1,dim2))
   allocate(Scan_Space_Counts_Avhrr(5,dim2))
   allocate(Scan_BB_Counts_Avhrr(5,dim2))
   allocate(Chan_Counts_Avhrr_Sg(5,dim1,dim2))
   allocate(IR_Coef_1_1b(3:5,dim2))
   allocate(IR_Coef_2_1b(3:5,dim2))
   allocate(IR_Coef_3_1b(3:5,dim2))
   allocate(IR_Linear_Slope_1b(3:5,dim2))
   allocate(IR_Linear_Intercept_1b(3:5,dim2))
   allocate(IR_Linear_Slope_New(3:5,dim2))
   allocate(IR_Linear_Intercept_New(3:5,dim2))
   allocate(Vis_Slope_1(1:3,dim2))
   allocate(Vis_Slope_2(1:3,dim2))
   allocate(Vis_Intercept_1(1:3,dim2))
   allocate(Vis_Intercept_2(1:3,dim2))
   allocate(Vis_Intersection(1:3,dim2))
end subroutine CREATE_AVHRR_ARRAYS
subroutine RESET_AVHRR_ARRAYS()
   Buffer_AVHRR = 0
   Segment_Buffer_AVHRR = 0
   Fatal_AVHRR = 0
   Resync_AVHRR = 0
   Clavr_Status_1b_AVHRR = Missing_Value_Int1
   Chan_Counts_Avhrr = Missing_Value_Int2
   Scan_Space_Counts_Avhrr = Missing_Value_Int2
   Scan_BB_Counts_Avhrr = Missing_Value_Int2
   Chan_Counts_Avhrr_Sg = Missing_Value_Int2
   IR_Coef_1_1b = Missing_Value_Real4
   IR_Coef_2_1b = Missing_Value_Real4
   IR_Coef_3_1b = Missing_Value_Real4
   IR_Linear_Slope_1b = Missing_Value_Real4
   IR_Linear_Intercept_1b = Missing_Value_Real4
   IR_Linear_Slope_New = Missing_Value_Real4
   IR_Linear_Intercept_New = Missing_Value_Real4
   Vis_Slope_1 = Missing_Value_Real4
   Vis_Slope_2 = Missing_Value_Real4
   Vis_Intercept_1 = Missing_Value_Real4
   Vis_Intercept_2 = Missing_Value_Real4
   Vis_Intersection = Missing_Value_Real4
end subroutine RESET_AVHRR_ARRAYS
subroutine DESTROY_AVHRR_ARRAYS()
   deallocate(Buffer_AVHRR)
   deallocate(Segment_Buffer_AVHRR)
   deallocate(Fatal_AVHRR)
   deallocate(Resync_AVHRR)
   deallocate(Clavr_Status_1b_AVHRR)
   deallocate(Chan_Counts_Avhrr)
   deallocate(Scan_Space_Counts_Avhrr)
   deallocate(Scan_BB_Counts_Avhrr)
   deallocate(Chan_Counts_Avhrr_Sg)
   deallocate(IR_Coef_1_1b)
   deallocate(IR_Coef_2_1b)
   deallocate(IR_Coef_3_1b)
   deallocate(IR_Linear_Slope_1b)
   deallocate(IR_Linear_Intercept_1b)
   deallocate(IR_Linear_Slope_New)
   deallocate(IR_Linear_Intercept_New)
   deallocate(Vis_Slope_1)
   deallocate(Vis_Slope_2)
   deallocate(Vis_Intercept_1)
   deallocate(Vis_Intercept_2)
   deallocate(Vis_Intersection)
end subroutine DESTROY_AVHRR_ARRAYS
!-----------------------------------------------------------
! this converts the dual gain reflectance counts to single
! gain for the AVHRR_KLM_Flag sensors
!-----------------------------------------------------------
subroutine CONVERT_AVHRR_COUNTS_SINGLE_GAIN(AVHRR_KLM_Flag,j1,j2)
 integer, intent(in):: j1,j2
 integer, intent(in):: AVHRR_KLM_Flag
 integer:: i, j


 if (AVHRR_FLAG == sym%NO) return

 if (Chan_On_Flag_Default(1) == sym%NO) return
 if (Chan_On_Flag_Default(2) == sym%NO) return
 if (Chan_On_Flag_Default(6) == sym%NO) return

 if (AVHRR_KLM_Flag == sym%NO) then

  Chan_Counts_Avhrr_Sg(:,:,j1:j2) = Chan_Counts_Avhrr(:,:,j1:j2)
  do j = j1,j1+j2-1
    Chan_Counts_Avhrr_Sg(1,:,j) = Chan_Counts_Avhrr_Sg(1,:,j) - Scan_Space_Counts_Avhrr(1,j)
    Chan_Counts_Avhrr_Sg(2,:,j) = Chan_Counts_Avhrr_Sg(2,:,j) - Scan_Space_Counts_Avhrr(2,j)
    if (Ch3a_On_AVHRR(j) == sym%YES) then
      Chan_Counts_Avhrr_Sg(3,:,j) = Chan_Counts_Avhrr_Sg(3,:,j) - Scan_Space_Counts_Avhrr(3,j)
    else
      Chan_Counts_Avhrr_Sg(3,:,j) = Missing_Value_Int2
    endif
  enddo

 else

  do j = j1, j1+j2-1

   !--- check for bad scans (note Bad_Pixel_Mask does yet exist)
   if (Bad_Scan_Flag(j) == sym%YES) then
        cycle
   endif

   do i = 1, Num_Pix

     !--- channel 1
     if (Chan_On_Flag_Default(1) == sym%YES) then
      if (Chan_Counts_Avhrr(1,i,j) < Ch1_Switch_Count) then
         Chan_Counts_Avhrr_Sg(1,i,j) = Scan_Space_Counts_Avhrr(1,j) + 0.5*(Chan_Counts_Avhrr(1,i,j) - Scan_Space_Counts_Avhrr(1,j))
       else
         Chan_Counts_Avhrr_Sg(1,i,j) = (Scan_Space_Counts_Avhrr(1,j) + 0.5*(Ch1_Switch_Count - Scan_Space_Counts_Avhrr(1,j))) + &
                                 3.0*(Chan_Counts_Avhrr(1,i,j) - Ch1_Switch_Count)/2.0
       endif
       Chan_Counts_Avhrr_Sg(1,i,j) = Chan_Counts_Avhrr_Sg(1,i,j) - Scan_Space_Counts_Avhrr(1,j)
     endif

     !--- channel 2
     if (Chan_On_Flag_Default(2) == sym%YES) then
      if (Chan_Counts_Avhrr(2,i,j) < Ch2_Switch_Count) then
       Chan_Counts_Avhrr_Sg(2,i,j) = Scan_Space_Counts_Avhrr(2,j) + 0.5*(Chan_Counts_Avhrr(2,i,j) - Scan_Space_Counts_Avhrr(2,j))
      else
       Chan_Counts_Avhrr_Sg(2,i,j) = (Scan_Space_Counts_Avhrr(2,j) + 0.5*(Ch2_Switch_Count - Scan_Space_Counts_Avhrr(2,j))) + &
                                 3.0*(Chan_Counts_Avhrr(2,i,j) - Ch2_Switch_Count)/2.0
      endif
      Chan_Counts_Avhrr_Sg(2,i,j) = Chan_Counts_Avhrr_Sg(2,i,j) - Scan_Space_Counts_Avhrr(1,j)
    endif

    !--- channel 3a
    if (Chan_On_Flag_Default(6) == sym%YES) then

     if (Ch3a_On_AVHRR(j) == sym%YES) then

        if (Chan_Counts_Avhrr(3,i,j) < Ch3a_Switch_Count) then
           Chan_Counts_Avhrr_Sg(3,i,j) = Scan_Space_Counts_Avhrr(3,j) + 0.25*(Chan_Counts_Avhrr(3,i,j) - Scan_Space_Counts_Avhrr(3,j))
        else
           Chan_Counts_Avhrr_Sg(3,i,j) = (Scan_Space_Counts_Avhrr(3,j) + 0.25*(Ch3a_Switch_Count - Scan_Space_Counts_Avhrr(3,j))) + &
                                 1.75*(Chan_Counts_Avhrr(3,i,j) - Ch3a_Switch_Count)
        endif

        Chan_Counts_Avhrr_Sg(3,i,j) = Chan_Counts_Avhrr_Sg(3,i,j) - Scan_Space_Counts_Avhrr(1,j)

     else

        Chan_Counts_Avhrr_Sg(3,i,j) = Missing_Value_Int2

     endif

    endif

    enddo
   enddo

 endif

 Ch1_Counts = Chan_Counts_Avhrr_Sg(1,:,:)
 Ch2_Counts = Chan_Counts_Avhrr_Sg(2,:,:)
 Ch6_Counts = Chan_Counts_Avhrr_Sg(3,:,:)

end subroutine CONVERT_AVHRR_COUNTS_SINGLE_GAIN

!------------------------------------------- ! end of module !----------------------------------------------------------- end 
end module AVHRR_MODULE
