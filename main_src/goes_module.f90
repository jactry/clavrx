!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: goes_module.f90 (src)
!       GOES_MODULE (program)
!
! PURPOSE: This module handles the processing of GOES Level-1b files in PATMOS-x
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  CIRA provided many of the core GOES navigation routines
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
!  This code subject to US Govt Copyright Regulations
!--------------------------------------------------------------------------------------
module GOES_MODULE
use CONSTANTS
use PIXEL_COMMON
use CALIBRATION_CONSTANTS
use PLANCK
use NUMERICAL_ROUTINES
use FILE_UTILITY
use VIEWING_GEOMETRY_MODULE
implicit none

public:: GET_GOES_HEADERS, &
         READ_GOES, &
         READ_GOES_SNDR, &
         GET_IMAGE_FROM_AREAFILE, &
         COMPUTE_SATELLITE_ANGLES, &
         LMODEL, &
         DETERMINE_DARK_COMPOSITE_NAME, &
         ASSIGN_GOES_SAT_ID_NUM_INTERNAL, &
         ASSIGN_GOES_SNDR_ID_NUM_INTERNAL, &
         READ_GOES_SNDR_INSTR_CONSTANTS, &
         READ_GOES_INSTR_CONSTANTS, &
         READ_DARK_COMPOSITE_COUNTS, &
         CALIBRATE_GOES_DARK_COMPOSITE, &
         POST_PROCESS_GOES_DARK_COMPOSITE, &
         DARK_COMPOSITE_CLOUD_MASK

private:: GET_GOES_NAVIGATION, &
          CALIBRATE_GOES_REFLECTANCE, &
          READ_IMGR_RP, &
          READ_IMGR_SIN, &
          READ_IMGR_MON, &
          FIX_VAR, &
          FIX_VAR2, &
          LPOINT, &
          COMP_ES, &
          COMP_LP, &
          INST2E, &
          GPOINT, &
          TIME50
private

 integer(kind=int4), public, parameter:: Goes_Xstride = 2    ! goes is oversampled by 50% in x
 integer(kind=int4), public, parameter:: Goes_Sndr_Xstride = 1
 integer(kind=int4), private, parameter:: Num_4km_Scans_Goes_Fd = 2704 
 integer(kind=int4), private, parameter:: Num_4km_Elem_Goes_Fd = 5200 
 integer(kind=int4), private, parameter:: time_for_fd_Scan_goes =  1560000 !milliseconds
 real, private, save:: Scan_rate    !scan rate in millsec / line
 character(len=120), save, public:: Dark_Comp_Data_Dir_Temp
 integer(kind=int4), private, parameter:: Goes_Imager_Byte_Shift = -5
 integer(kind=int4), private, parameter:: Goes_Sounder_Byte_Shift = 0  !I don't know this
 real(kind=real4), private:: GEO_ALTITUDE = 35786.0 !km

!-----------------------------------------------------------------------------
! define derived data types used for holding GVAR parameters
!-----------------------------------------------------------------------------
!
!C	Imager ONA repeat sinusoid T.
!C

        type, public:: IMGR_SIN 
          integer :: mag_Sinu
          integer :: phase_ang_Sinu
        end type IMGR_SIN
!C
!C	Imager repeat monomial T.
!C
        type, public:: IMGR_MON
          integer :: order_appl_Sinu
          integer :: order_mono_Sinu
          integer :: mag_mono_Sinu
          integer :: phase_ang_Sinu
          integer :: ang_from_epoch
        end type IMGR_MON
!C
!C	Imager repeat T.
!C
        type, public:: IMGR_RP 
           integer :: exp_mag
           integer :: exp_Time_const
           integer :: mean_att_ang_const
           integer :: num_Sinu_per_angle
           type(IMGR_SIN), dimension(15):: sinusoid
           integer :: num_mono_Sinu
           type(IMGR_MON), dimension(4) :: monomial
        end type IMGR_RP
!
!     Define McIDAS area structre
!
type, public:: AREA_STRUCT
   integer ::Area_Status           ! Area status
!  integer(kind=int4) :: Swap_Bytes !whether or not bytes need to be swapped
   integer ::Version_Num           ! Area version number
   integer ::Sat_Id_Num            ! Satellite ID (SSS)
   integer ::Img_Date              ! Year and Julian day (YYDDD)
   integer ::Img_Time              ! Time of image (HHMMSS)
   integer ::North_Bound           ! Upper left line in Sat coords (Y-Coord)
   integer ::West_Vis_Pixel        ! Upper left element in Sat coords (X-Coord)
   integer ::Z_Coor                ! Upper z-coord (Z-Coord)
   integer ::Num_Line              ! Number of lines in image (Y-SIZE)
   integer ::Num_Elem              ! Number of elememts in image (X-SIZE)
   integer ::Bytes_Per_Pixel       ! Number of bytes per data element
   integer ::Line_Res              ! Line resolution (Y-RES) 
   integer ::Elem_Res              ! Element resolution (X-RES)
   integer ::Num_Chan              ! Number of bands (Z-RES)
   integer ::Num_Byte_Ln_Prefix    ! Number of bytes in line prefix (multiple of 4)
   integer ::Proj_Num              ! Project number
   integer ::Creation_Date         ! Creation date (YYDDD)
   integer ::Creation_Time         ! Creation time (HHMMSS)
   integer ::Sndr_Filter_Map       ! Filter map for soundings
   integer ::img_id_num            ! Image ID number
   integer, dimension(4) :: id     ! Reserved for radar appications
   character(len=32):: comment     ! 32 char comments
   integer ::pri_key_calib         ! Calibration colicil (area number)
   integer ::pri_key_nav           ! Primary navigation codicil  (data)
   integer ::sec_key_nav           ! Secondary navigation codicil (nav)
   integer ::val_code              ! Validity code
   integer, dimension(8) :: pdl    ! PDL in packed-byte format
   integer ::band8                 ! Where band-8 came from
   integer ::act_Img_Date          ! Actual image start day (YYDDD)
   integer ::act_Img_Time          ! Actual image start time (HHMMSS)
   integer ::act_Start_Scan        ! Actual start scan
   integer ::len_prefix_doc        ! Length of prefix documentation (DOC)
   integer ::len_prefix_calib      ! Length of prefix calibration (CAL)
   integer ::len_prefix_lev        ! Length of prefix level (LEV)
   character(len=4)::src_Type      ! Source type
   character(len=4)::calib_Type    ! Calibration type
   integer ::avg_or_Sample         ! Data was averaged (0) or sampled (1)
   integer ::poes_Signal           ! LAC, GAC, HRPT
   integer ::poes_up_down          ! POES ascending/descending
   integer ::cal_offset            ! needed for SEVIRI
   character(len=4) :: orig_Src_Type   ! Original source type of data
   integer, dimension(7) ::reserved    ! Reserved (6 calib pointer)
end type AREA_STRUCT
!
!	Define GVAR Navigation block
!
        type, public:: GVAR_NAV
          character(len=4):: nav_Type
          character(len=4):: IMC_Status
          integer, dimension(12) ::  spare1
          integer ::    stat
          integer ::    ref_long
          integer ::    ref_rad_dist
          integer ::    ref_lat
          integer ::    ref_orb_yaw
          integer ::    ref_att_roll
          integer ::    ref_att_pitch
          integer ::    ref_att_yaw
          integer(kind=int1), dimension(8)::  epoch_Time  ! BCD_TIME
          integer ::    start_Time
          integer ::    IMC_corr_roll
          integer ::    IMC_corr_pitch
          integer ::    IMC_corr_yaw
          integer, dimension(13) ::    ref_long_change
          integer, dimension(11) ::    ref_rad_dist_change
          integer, dimension(9) ::    sine_lat
          integer, dimension(9) ::    sine_orb_yaw
          integer ::    solar_rate
          integer ::    exp_Start_Time
          type(IMGR_RP) :: roll_att
          integer, dimension(10) ::     spare2
          character(len=4):: more1
          character(len=4):: gvar1
          type(IMGR_RP) :: pitch_att
          type(IMGR_RP) :: yaw_att
          integer, dimension(16) ::    spare3
          character(len=4):: more2
          character(len=4):: gvar2
          type(IMGR_RP) :: roll_misalgn
          type(IMGR_RP) :: pitch_misalgn
          integer ::    Img_Date
          integer ::    Img_Time
          integer ::    instr
          integer, dimension(9) ::    spare4
          integer ::    NS_CYL
          integer ::    EW_CYL
          integer ::    NS_INC
          integer ::    EW_INC
          character(len=4):: more3
          character(len=4):: gvar3
          integer, dimension(126) ::     spare5
          character(len=4):: more4
          character(len=4):: gvar4
          integer, dimension(127) ::     spare6

          !--Needed for MTSAT
          integer(kind=int4), dimension(672,4) :: MAP
          integer(kind=int4), dimension(25,25,4) :: JSMT
          character(len=1), dimension(3200) :: COBAT
	  real(kind=real4), dimension(4) :: RESLIN,RESELM,RLIC,RELMFC,SENSSU, &
                                            RLINE,RELMNT
          real(kind=real4), dimension(3) :: VMIS
          real(kind=real4), dimension(8) :: RINF
          real(kind=real4), dimension(3,3) :: ELMIS
          real(kind=real8) :: DSPIN,DTIMS,DSCT
          real(kind=real8), dimension(10,33) :: ATIT
          real(kind=real8), dimension(35,8) :: ORBT1
		  real(kind=real4) :: sublon, sublat

          !--Needed for HRIT MTSAT
          real(KIND(0.0d0)) :: sub_lon
          integer(kind=int4) :: CFAC, LFAC, COFF, LOFF

        end type GVAR_NAV
!
!     Calibration table. 
!    
      type, public:: MC_IMGR_CAL_VIS_DET_T
        integer(kind=int4), dimension(8) :: det
      end type MC_IMGR_CAL_VIS_DET_T

      type, public:: MC_IMGR_CAL_IR_DET_T
        integer(kind=int4) :: det0
        integer(kind=int4) :: det2
        integer(kind=int4) :: det4
        integer(kind=int4) :: det6
      end type MC_IMGR_CAL_IR_DET_T

      type, public:: MC_IMGR_CAL_ALL_IR_DET_T 
        type(MC_IMGR_CAL_IR_DET_T):: PRI
        type(MC_IMGR_CAL_IR_DET_T):: SEC
      end type MC_IMGR_CAL_ALL_IR_DET_T

      type, public:: CALIB_STRUCT
          type(MC_IMGR_CAL_VIS_DET_T):: VIS_BIAS
          type(MC_IMGR_CAL_VIS_DET_T):: VIS_GAIN1
          type(MC_IMGR_CAL_VIS_DET_T):: VIS_GAIN2
          integer(kind=int4) :: vis_albedo
          type (MC_IMGR_CAL_ALL_IR_DET_T):: IR_SCL_BIAS
          type (MC_IMGR_CAL_ALL_IR_DET_T):: IR_SCL_GAIN1
          integer(kind=int4), dimension(8) :: unknown
          integer(kind=int4), dimension(79) :: z
      end type CALIB_STRUCT

contains

!--------------------------------------------------------------------
! assign internal sat id's and const file names for goes
!--------------------------------------------------------------------
subroutine ASSIGN_GOES_SAT_ID_NUM_INTERNAL(Mcidas_Id_Num, Elem_Res)
    integer(kind=int4), intent(in):: Mcidas_Id_Num
    integer(kind=int4), intent(in):: Elem_Res

    Goes_Mop_Flag = sym%NO
    Goes_1km_Flag = sym%NO
    Goes_Sndr_Flag = sym%NO ! Make sure sounder is turned off
    if (Mcidas_Id_Num == 70)   then
        Sc_Id_WMO = 252
        Instr_Const_file = 'goes_08_instr.dat'
        Algo_Const_file = 'goes_08_algo.dat'
        Goes_Mop_Flag = sym%NO
        Platform_Name_Attribute = 'GOES-8'
        Sensor_Name_Attribute = 'IMAGER'
    endif
    if (Mcidas_Id_Num == 72)   then
        Instr_Const_file = 'goes_09_instr.dat'
        Algo_Const_file = 'goes_09_algo.dat'
        Goes_Mop_Flag = sym%NO
        Platform_Name_Attribute = 'GOES-9'
        Sensor_Name_Attribute = 'IMAGER'
    endif
    if (Mcidas_Id_Num == 74)   then
        Sc_Id_WMO = 254
        Instr_Const_file = 'goes_10_instr.dat'
        Algo_Const_file = 'goes_10_algo.dat'
        Goes_Mop_Flag = sym%NO
        Platform_Name_Attribute = 'GOES-10'
        Sensor_Name_Attribute = 'IMAGER'
    endif
   if (Mcidas_Id_Num == 76)   then
        Sc_Id_WMO = 255
        Instr_Const_file = 'goes_11_instr.dat'
        Algo_Const_file = 'goes_11_algo.dat'
        Goes_Mop_Flag = sym%NO
        Platform_Name_Attribute = 'GOES-11'
        Sensor_Name_Attribute = 'IMAGER'
    endif
    if (Mcidas_Id_Num == 78)   then
        Sc_Id_WMO = 256
        Instr_Const_file = 'goes_12_instr.dat'
        Algo_Const_file = 'goes_12_algo.dat'
        Goes_Mop_Flag = sym%YES
        Platform_Name_Attribute = 'GOES-12'
        Sensor_Name_Attribute = 'IMAGER'
    endif
    if (Mcidas_Id_Num == 180)   then
        Sc_Id_WMO = 257
        Instr_Const_file = 'goes_13_instr.dat'
        Algo_Const_file = 'goes_13_algo.dat'
        Goes_Mop_Flag = sym%YES
        Platform_Name_Attribute = 'GOES-13'
        Sensor_Name_Attribute = 'IMAGER'
    endif
    if (Mcidas_Id_Num == 182)   then
        Sc_Id_WMO = 258
        Instr_Const_file = 'goes_14_instr.dat'
        Algo_Const_file = 'goes_14_algo.dat'
        Goes_Mop_Flag = sym%YES
        Platform_Name_Attribute = 'GOES-14'
        Sensor_Name_Attribute = 'IMAGER'
    endif
    if (Mcidas_Id_Num == 184)   then
        Sc_Id_WMO = 259
        Instr_Const_file = 'goes_15_instr.dat'
        Algo_Const_file = 'goes_15_algo.dat'
        Goes_Mop_Flag = sym%YES
        Platform_Name_Attribute = 'GOES-15'
        Sensor_Name_Attribute = 'IMAGER'
    endif

   Instr_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Instr_Const_file)
   Algo_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Algo_Const_file)

   if (Goes_Mop_Flag == sym%YES) then
           Chan_On_Flag_Default(32) = sym%NO
   else
           Chan_On_Flag_Default(33) = sym%NO
   endif


   !---- Set flag for processing 1 km GOES Imager data
   if (Goes_Sndr_Flag == sym%NO .and. Elem_Res == 1) then 
      Goes_1km_Flag = sym%YES
   endif

end subroutine ASSIGN_GOES_SAT_ID_NUM_INTERNAL
!--------------------------------------------------------------------
! assign internal sat id's and const file names for goes sounder
!--------------------------------------------------------------------
subroutine ASSIGN_GOES_SNDR_ID_NUM_INTERNAL(Mcidas_Id_Num)
    integer(kind=int4), intent(in):: Mcidas_Id_Num

    Goes_Sndr_Flag = sym%YES

    if (Mcidas_Id_Num == 71)   then
        Sc_Id_WMO = 252
        Instr_Const_file = 'goes_8_Sndr_instr.dat'
        Algo_Const_file = 'goes_8_Sndr_algo.dat'
        Platform_Name_Attribute = 'GOES-8'
        Sensor_Name_Attribute = 'SOUNDER'
    endif
    if (Mcidas_Id_Num == 73)   then
        Sc_Id_WMO = 253
        Instr_Const_file = 'goes_9_Sndr_instr.dat'
        Algo_Const_file = 'goes_9_Sndr_algo.dat'
        Platform_Name_Attribute = 'GOES-9'
        Sensor_Name_Attribute = 'SOUNDER'
    endif
    if (Mcidas_Id_Num == 75)   then
        Sc_Id_WMO = 254
        Instr_Const_file = 'goes_10_Sndr_instr.dat'
        Algo_Const_file = 'goes_10_Sndr_algo.dat'
        Platform_Name_Attribute = 'GOES-10'
        Sensor_Name_Attribute = 'SOUNDER'
    endif
    if (Mcidas_Id_Num == 77)   then
        Sc_Id_WMO = 255
        Instr_Const_file = 'goes_11_Sndr_instr.dat'
        Algo_Const_file = 'goes_11_Sndr_algo.dat'
        Platform_Name_Attribute = 'GOES-11'
        Sensor_Name_Attribute = 'SOUNDER'
    endif
    if (Mcidas_Id_Num == 79)   then
        Sc_Id_WMO = 256
        Instr_Const_file = 'goes_12_Sndr_instr.dat'
        Algo_Const_file = 'goes_12_Sndr_algo.dat'
        Platform_Name_Attribute = 'GOES-12'
        Sensor_Name_Attribute = 'SOUNDER'
    endif
    if (Mcidas_Id_Num == 181)   then
        Sc_Id_WMO = 257
        Instr_Const_file = 'goes_13_Sndr_instr.dat'
        Algo_Const_file = 'goes_13_Sndr_algo.dat'
        Platform_Name_Attribute = 'GOES-13'
        Sensor_Name_Attribute = 'SOUNDER'
    endif
    if (Mcidas_Id_Num == 183)   then
        Sc_Id_WMO = 258
        Instr_Const_file = 'goes_14_Sndr_instr.dat'
        Algo_Const_file = 'goes_14_Sndr_algo.dat'
        Platform_Name_Attribute = 'GOES-14'
        Sensor_Name_Attribute = 'SOUNDER'
    endif
    if (Mcidas_Id_Num == 185)   then
        Sc_Id_WMO = 259
        Instr_Const_file = 'goes_15_Sndr_instr.dat'
        Algo_Const_file = 'goes_15_Sndr_algo.dat'
        Platform_Name_Attribute = 'GOES-15'
        Sensor_Name_Attribute = 'SOUNDER'
    endif

   Instr_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Instr_Const_file)
   Algo_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Algo_Const_file)

end subroutine ASSIGN_GOES_SNDR_ID_NUM_INTERNAL

!----------------------------------------------------------------
! read the goes constants into memory
!-----------------------------------------------------------------
subroutine READ_GOES_INSTR_CONSTANTS(Instr_Const_file)
 character(len=*), intent(in):: Instr_Const_file
 integer:: ios0, erstat
 integer:: Instr_Const_lun

 Instr_Const_lun = GET_LUN()

 open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

 erstat = 0
 if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, "Error opening GOES constants file, ios0 = ", ios0
    stop 19
 endif

  read(unit=Instr_Const_lun,fmt="(a3)") sat_name
  read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
  read(unit=Instr_Const_lun,fmt=*) ew_Ch20
  read(unit=Instr_Const_lun,fmt=*) launch_date
  read(unit=Instr_Const_lun,fmt=*) Ch1_Dark_Count
  read(unit=Instr_Const_lun,fmt=*) Ch1_Gain_Low_0,Ch1_Degrad_Low_1, Ch1_Degrad_Low_2
  read(unit=Instr_Const_lun,fmt=*) a1_20, a2_20, nu_20
  read(unit=Instr_Const_lun,fmt=*) a1_27, a2_27, nu_27
  read(unit=Instr_Const_lun,fmt=*) a1_31, a2_31, nu_31
  read(unit=Instr_Const_lun,fmt=*) a1_32, a2_32, nu_32
  read(unit=Instr_Const_lun,fmt=*) a1_33, a2_33, nu_33
  read(unit=Instr_Const_lun,fmt=*) goes_Ch2_Thermal_Slope,goes_Ch2_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch3_Thermal_Slope,goes_Ch3_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch4_Thermal_Slope,goes_Ch4_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch5_Thermal_Slope,goes_Ch5_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch6_Thermal_Slope,goes_Ch6_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) b1_day_mask,b2_day_mask,b3_day_mask,b4_day_mask
  close(unit=Instr_Const_lun)

  !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / ew_Ch20

end subroutine READ_GOES_INSTR_CONSTANTS

!----------------------------------------------------------------
! read the goes sounder constants into memory
!-----------------------------------------------------------------
subroutine READ_GOES_SNDR_INSTR_CONSTANTS(Instr_Const_file)
 character(len=*), intent(in):: Instr_Const_file
 integer:: ios0, erstat
 integer:: Instr_Const_lun

 Instr_Const_lun = GET_LUN()

 open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

 erstat = 0
 if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, "Error opening GOES constants file, ios0 = ", ios0
    stop 19
 endif

  read(unit=Instr_Const_lun,fmt="(a3)") sat_name
  read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
  read(unit=Instr_Const_lun,fmt=*) ew_Ch20
  read(unit=Instr_Const_lun,fmt=*) launch_date
  read(unit=Instr_Const_lun,fmt=*) Ch1_Dark_Count
  read(unit=Instr_Const_lun,fmt=*) Ch1_Gain_Low_0,Ch1_Degrad_Low_1, Ch1_Degrad_Low_2
  read(unit=Instr_Const_lun,fmt=*) a1_20, a2_20, nu_20
  read(unit=Instr_Const_lun,fmt=*) a1_21, a2_21, nu_21
  read(unit=Instr_Const_lun,fmt=*) a1_23, a2_23, nu_23
  read(unit=Instr_Const_lun,fmt=*) a1_24, a2_24, nu_24
  read(unit=Instr_Const_lun,fmt=*) a1_25, a2_25, nu_25
  read(unit=Instr_Const_lun,fmt=*) a1_27, a2_27, nu_27
  read(unit=Instr_Const_lun,fmt=*) a1_28, a2_28, nu_28
  read(unit=Instr_Const_lun,fmt=*) a1_30, a2_30, nu_30
  read(unit=Instr_Const_lun,fmt=*) a1_31, a2_31, nu_31
  read(unit=Instr_Const_lun,fmt=*) a1_32, a2_32, nu_32
  read(unit=Instr_Const_lun,fmt=*) a1_33, a2_33, nu_33
  read(unit=Instr_Const_lun,fmt=*) a1_34, a2_34, nu_34
  read(unit=Instr_Const_lun,fmt=*) a1_35, a2_35, nu_35
  read(unit=Instr_Const_lun,fmt=*) a1_36, a2_36, nu_36
  read(unit=Instr_Const_lun,fmt=*) goes_Ch2_Thermal_Slope,goes_Ch2_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch3_Thermal_Slope,goes_Ch3_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch4_Thermal_Slope,goes_Ch4_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch5_Thermal_Slope,goes_Ch5_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch7_Thermal_Slope,goes_Ch7_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch8_Thermal_Slope,goes_Ch8_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch9_Thermal_Slope,goes_Ch9_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch10_Thermal_Slope,goes_Ch10_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch12_Thermal_Slope,goes_Ch12_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch13_Thermal_Slope,goes_Ch13_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch14_Thermal_Slope,goes_Ch14_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch16_Thermal_Slope,goes_Ch16_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch17_Thermal_Slope,goes_Ch17_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) goes_Ch18_Thermal_Slope,goes_Ch18_Thermal_Intercept
  read(unit=Instr_Const_lun,fmt=*) b1_day_mask,b2_day_mask,b3_day_mask,b4_day_mask
  close(unit=Instr_Const_lun)

  !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / ew_Ch20

end subroutine READ_GOES_SNDR_INSTR_CONSTANTS


!-------------------------------------------------------------------------------
! public routine to read data from an AREA file for one segment into memory
!-------------------------------------------------------------------------------
subroutine READ_GOES(Segment_Number,Channel_1_Filename, &
                     jday, image_Time_ms, Time_Since_Launch, &
                     AREAstr,NAVstr)

   integer(kind=int4), intent(in):: Segment_Number
   character(len=*), intent(in):: Channel_1_Filename
   type (AREA_STRUCT), intent(in) :: AREAstr
   type (GVAR_NAV), intent(in)    :: NAVstr
   integer(kind=int2), intent(in):: jday
   integer(kind=int4), intent(in):: image_Time_ms
   real(kind=real4), intent(in):: Time_Since_Launch

   character(len=120):: Channel_X_Filename
   character(len=120):: Channel_X_Filename_Full
   character(len=120):: Channel_X_Filename_Full_uncompressed
   character(len=180):: System_String

   integer:: ipos
   integer:: ilen
   integer:: Elem_Idx
   integer:: Line_Idx
   integer:: Ichan_Goes
   integer:: Ichan_Modis
   real(kind=real4):: image_Time_Hours
   integer(kind=int4):: image_jday
   integer(kind=int4):: First_Line_In_Segment
   character(len=1):: Ichan_Goes_String

   !--- assume Channel_1_file name has a unique "_1_" in the name. 
   !--- determine indices needed to replace that string
   ipos = index(Channel_1_Filename, "_1_")
   ilen = len(Channel_1_Filename)


!-------------------------------------------------------------------------------
! uncompress (only on first segment and channel 1 is already done during header read)
!-----------------------------------------------------------------------------
   if (Segment_Number == 1) then

       do Ichan_Goes = 2,6

       if (Ichan_Goes == 2) Ichan_Modis = 20
       if (Ichan_Goes == 3) Ichan_Modis = 27
       if (Ichan_Goes == 4) Ichan_Modis = 31
       if (Ichan_Goes == 5) Ichan_Modis = 32
       if (Ichan_Goes == 6) Ichan_Modis = 33

       write(Ichan_Goes_String,fmt="(I1)") Ichan_Goes

       if (Chan_On_Flag_Default(Ichan_Modis) == sym%YES) then

          Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_"//Ichan_Goes_String//"_" // &
                            Channel_1_Filename(ipos+3:ilen)

          if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
          else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
          endif

          Channel_X_Filename_Full_uncompressed = trim(Dir_1b)//trim(Channel_X_Filename)
          if (L1b_Gzip == sym%YES) then
              System_String = "gunzip -c "//trim(Channel_X_Filename_Full_uncompressed)//".gz"// &
                                " > "//trim(Channel_X_Filename_Full)
              call system(System_String)

              Number_of_Temporary_Files = Number_of_Temporary_Files + 1
              Temporary_File_Name(Number_of_Temporary_Files) = trim(Channel_X_Filename)

          endif
          if (l1b_bzip2 == sym%YES) then
              System_String = "bunzip2 -c "//trim(Channel_X_Filename_Full_uncompressed)//".bz2"// &
                                " > "//trim(Channel_X_Filename_Full)
              call system(System_String)

              Number_of_Temporary_Files = Number_of_Temporary_Files + 1
              Temporary_File_Name(Number_of_Temporary_Files) = trim(Channel_X_Filename)

          endif

      endif

    enddo

   endif


   !---   read channel 1 (GOES channel 1)
   if (Chan_On_Flag_Default(1) == sym%YES) then

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_1_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_1_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Imager_Byte_Shift, &
                                    AREAstr, Goes_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

     !--- calibrate
     call CALIBRATE_GOES_REFLECTANCE(Two_Byte_Temp,Ch1_Dark_Count,Time_Since_Launch,ch(1)%Ref_Toa)

     !--- store ch1 counts for support of PATMOS-x calibration studies
     Ch1_Counts = Two_Byte_Temp

   endif

   !---   read channel 20 (GOES channel 2)
   if (Chan_On_Flag_Default(20) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_2_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Imager_Byte_Shift, &
                                    AREAstr, Goes_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)
      ch(20)%Rad_Toa = (Two_Byte_Temp - Goes_Ch2_Thermal_Intercept)/Goes_Ch2_Thermal_Slope

   endif

   !---   read channel 27 (GOES channel 3)
   if (Chan_On_Flag_Default(27) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_3_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Imager_Byte_Shift, &
                                    AREAstr, Goes_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)
      ch(27)%Rad_Toa = (Two_Byte_Temp - Goes_Ch3_Thermal_Intercept)/Goes_Ch3_Thermal_Slope

   endif

   !---   read channel 31 (GOES channel 4)
   if (Chan_On_Flag_Default(31) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_4_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Imager_Byte_Shift, &
                                    AREAstr, Goes_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)
        ch(31)%Rad_Toa = (Two_Byte_Temp - Goes_Ch4_Thermal_Intercept)/Goes_Ch4_Thermal_Slope

   endif

   !---   read channel 32 (GOES channel 5)
   if (Chan_On_Flag_Default(32) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_5_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Imager_Byte_Shift, &
                                    AREAstr, Goes_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(32)%Rad_Toa = (Two_Byte_Temp - Goes_Ch5_Thermal_Intercept)/Goes_Ch5_Thermal_Slope

   endif

   !---   read channel 33 (GOES channel 6)
   if (Chan_On_Flag_Default(33) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_6_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Imager_Byte_Shift, &
                                    AREAstr, Goes_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(33)%Rad_Toa = (Two_Byte_Temp - Goes_Ch6_Thermal_Intercept)/Goes_Ch6_Thermal_Slope

   endif

   !--------------------------------------------------------------------------
   ! Compute IR Brightness Temperature
   !--------------------------------------------------------------------------
   do Elem_Idx = 1,num_pix
     do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1

        if (Chan_On_Flag_Default(20) == sym%YES) then
           if (ch(20)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(20)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(20,ch(20)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(27) == sym%YES) then
           if (ch(27)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(27)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(27,ch(27)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(31) == sym%YES) then
           if (ch(31)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(31)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(31,ch(31)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(32) == sym%YES) then
           if (ch(32)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(32)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(32,ch(32)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(33) == sym%YES) then
           if (ch(33)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(33)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(33,ch(33)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif
     enddo
   enddo

   !------------------------------------------------------------------------------
   ! Goes Navigation
   !------------------------------------------------------------------------------
   call GET_GOES_NAVIGATION(Segment_Number, Num_Scans_Per_Segment, Num_Scans_Read, NAVstr, AREAstr, Goes_Xstride)

   !------------------------------------------------------------------------------
   ! Goes Angles
   !------------------------------------------------------------------------------
   image_jday = jday
   image_Time_Hours = image_Time_ms / 60.0 / 60.0 / 1000.0
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1
     do Elem_Idx = 1,num_pix
        call POSSOL(image_jday,image_Time_Hours, &
                    Lon_1b(Elem_Idx,Line_Idx),Lat_1b(Elem_Idx,Line_Idx), &
                    solzen(Elem_Idx,Line_Idx),solaz(Elem_Idx,Line_Idx))
     enddo
     call COMPUTE_SATELLITE_ANGLES(goes_Sub_Satellite_longitude,  &
                      goes_Sub_Satellite_latitude, Line_Idx)
   enddo

   !--- scan number and time
   First_Line_In_Segment = (Segment_Number-1)*Num_Scans_Per_Segment
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1
     Scan_number(Line_Idx) = First_Line_In_Segment + Line_Idx
         ! - if normal (header-based) scan_Time detection failed use==>  (AW 2012/02/23)
         if (Scan_Time(Line_Idx) == 0) then
           Scan_Time(Line_Idx) = image_Time_ms + (Scan_number(Line_Idx)-1) * Scan_rate
         endif
   enddo

   !--- ascending node
   Elem_Idx = num_pix/2
   do Line_Idx = Line_Idx_Min_Segment+1, Line_Idx_Min_Segment + Num_Scans_Read - 1
     ascend(Line_Idx) = 0
     if (Lat_1b(Elem_Idx,Line_Idx) < Lat_1b(Elem_Idx,Line_Idx-1)) then
       ascend(Line_Idx) = 1
     endif
   enddo
   ascend(Line_Idx_Min_Segment) = ascend(Line_Idx_Min_Segment+1)


end subroutine READ_GOES

!-------------------------------------------------------------------------------
! public routine to read data from an AREA file for one segment into memory
!-------------------------------------------------------------------------------
subroutine READ_GOES_SNDR(Segment_Number,Channel_1_Filename, &
                     jday, image_Time_ms, Time_Since_Launch, &
                     AREAstr,NAVstr)

   integer(kind=int4), intent(in):: Segment_Number
   character(len=*), intent(in):: Channel_1_Filename
   type (AREA_STRUCT), intent(in) :: AREAstr
   type (GVAR_NAV), intent(in)    :: NAVstr
   integer(kind=int2), intent(in):: jday
   integer(kind=int4), intent(in):: image_Time_ms
   real(kind=real4), intent(in):: Time_Since_Launch

   character(len=120):: Channel_X_Filename
   character(len=120):: Channel_X_Filename_Full
   character(len=120):: Channel_X_Filename_Full_uncompressed
   character(len=180):: System_String

   integer:: ipos
   integer:: ilen
   integer:: Elem_Idx
   integer:: Line_Idx
   integer:: Ichan_Goes
   integer:: Ichan_Modis
   real(kind=real4):: image_Time_Hours
   integer(kind=int4):: image_jday
   integer(kind=int4):: First_Line_In_Segment
   character(len=2):: Ichan_Goes_String

   !--- assume Channel_1_file name has a unique "_1_" in the name. 
   !--- determine indices needed to replace that string
   ipos = index(Channel_1_Filename, "_1_")
   ilen = len(Channel_1_Filename)


!-------------------------------------------------------------------------------
! uncompress (only on first segment and channel 1 is already done during header read)
!-----------------------------------------------------------------------------
   if (Segment_Number == 1) then

       do Ichan_Goes = 2,18

       if (ichan_goes == 2) ichan_modis = 36
       if (ichan_goes == 3) ichan_modis = 35
       if (ichan_goes == 4) ichan_modis = 34
       if (ichan_goes == 5) ichan_modis = 33
       if (ichan_goes == 7) ichan_modis = 32
       if (ichan_goes == 8) ichan_modis = 31
       if (ichan_goes == 9) ichan_modis = 30
       if (ichan_goes == 10) ichan_modis = 28
       if (ichan_goes == 12) ichan_modis = 27
       if (ichan_goes == 13) ichan_modis = 25
       if (ichan_goes == 14) ichan_modis = 24
       if (ichan_goes == 16) ichan_modis = 23
       if (ichan_goes == 17) ichan_modis = 21
       if (ichan_goes == 18) ichan_modis = 20

       if (ichan_goes < 10) then
          write(ichan_goes_string,fmt="(I1.1)") ichan_goes
       else
          write(ichan_goes_string,fmt="(I2.2)") ichan_goes
       endif

       if (Chan_On_Flag_Default(Ichan_Modis) == sym%YES) then

          Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_"//trim(Ichan_Goes_String)//"_" // &
                            Channel_1_Filename(ipos+3:ilen)

          if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
          else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
          endif

          Channel_X_Filename_Full_uncompressed = trim(Dir_1b)//trim(Channel_X_Filename)
          if (L1b_Gzip == sym%YES) then
              System_String = "gunzip -c "//trim(Channel_X_Filename_Full_uncompressed)//".gz"// &
                                " > "//trim(Channel_X_Filename_Full)
              call system(System_String)

              Number_of_Temporary_Files = Number_of_Temporary_Files + 1
              Temporary_File_Name(Number_of_Temporary_Files) = trim(Channel_X_Filename)

          endif
          if (l1b_bzip2 == sym%YES) then
              System_String = "bunzip2 -c "//trim(Channel_X_Filename_Full_uncompressed)//".bz2"// &
                                " > "//trim(Channel_X_Filename_Full)
              call system(System_String)

              Number_of_Temporary_Files = Number_of_Temporary_Files + 1
              Temporary_File_Name(Number_of_Temporary_Files) = trim(Channel_X_Filename)

          endif

      endif

    enddo

   endif


   !---   read channel 36 (GOES Sounder channel 2)
   if (Chan_On_Flag_Default(36) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_2_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(36)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(36)%Rad_Toa = (ch(36)%Rad_Toa - Goes_Ch2_Thermal_Intercept)/Goes_Ch2_Thermal_Slope

   endif


   !---   read channel 35 (GOES Sounder channel 3)
   if (Chan_On_Flag_Default(35) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_3_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(35)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(35)%Rad_Toa = (ch(35)%Rad_Toa - Goes_Ch3_Thermal_Intercept)/Goes_Ch3_Thermal_Slope

   endif

   !---   read channel 34 (GOES Sounder channel 4)
   if (Chan_On_Flag_Default(34) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_4_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(34)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(34)%Rad_Toa = (ch(34)%Rad_Toa - Goes_Ch4_Thermal_Intercept)/Goes_Ch4_Thermal_Slope

   endif
   
   !---   read channel 33 (GOES Sounder channel 5)
   if (Chan_On_Flag_Default(33) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_5_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(33)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(33)%Rad_Toa = (ch(33)%Rad_Toa - Goes_Ch5_Thermal_Intercept)/Goes_Ch5_Thermal_Slope

   endif


   !---   read channel 32 (GOES Sounder channel 7)
   if (Chan_On_Flag_Default(32) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_7_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(32)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(32)%Rad_Toa = (ch(32)%Rad_Toa - Goes_Ch7_Thermal_Intercept)/Goes_Ch7_Thermal_Slope

   endif
   
   
      !---   read channel 31 (GOES Sounder channel 8)
   if (Chan_On_Flag_Default(31) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_8_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(31)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(31)%Rad_Toa = (ch(31)%Rad_Toa - Goes_Ch8_Thermal_Intercept)/Goes_Ch8_Thermal_Slope
      
   endif
   
   
      !---   read channel 30 (GOES Sounder channel 9)
   if (Chan_On_Flag_Default(30) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_9_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(30)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(30)%Rad_Toa = (ch(30)%Rad_Toa - Goes_Ch9_Thermal_Intercept)/Goes_Ch9_Thermal_Slope

   endif
   
      !---   read channel 28 (GOES Sounder channel 10)
   if (Chan_On_Flag_Default(28) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_10_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(28)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(28)%Rad_Toa = (ch(28)%Rad_Toa - Goes_Ch10_Thermal_Intercept)/Goes_Ch10_Thermal_Slope

   endif
   
      !---   read channel 27 (GOES Sounder channel 12)
   if (Chan_On_Flag_Default(27) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_12_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(27)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(27)%Rad_Toa = (ch(27)%Rad_Toa - Goes_Ch12_Thermal_Intercept)/Goes_Ch12_Thermal_Slope

   endif
   
      !---   read channel 25 (GOES Sounder channel 13)
   if (Chan_On_Flag_Default(25) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_13_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)
       
      ch(25)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(25)%Rad_Toa = (ch(25)%Rad_Toa - Goes_Ch13_Thermal_Intercept)/Goes_Ch13_Thermal_Slope

   endif
   
   
      !---   read channel 24 (GOES Sounder channel 14)
   if (Chan_On_Flag_Default(24) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_14_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(24)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(24)%Rad_Toa = (ch(24)%Rad_Toa - Goes_Ch14_Thermal_Intercept)/Goes_Ch14_Thermal_Slope

   endif
   
      !---   read channel 23 (GOES Sounder channel 16)
   if (Chan_On_Flag_Default(23) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_16_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(23)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(23)%Rad_Toa = (ch(23)%Rad_Toa - Goes_Ch16_Thermal_Intercept)/Goes_Ch16_Thermal_Slope

   endif
   
      !---   read channel 21 (GOES Sounder channel 17)
   if (Chan_On_Flag_Default(21) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_17_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)
      ch(21)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(21)%Rad_Toa = (ch(21)%Rad_Toa - Goes_Ch17_Thermal_Intercept)/Goes_Ch17_Thermal_Slope

   endif
   
   
      !---   read channel 20 (GOES Sounder channel 18)
   if (Chan_On_Flag_Default(20) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_18_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      ch(20)%Rad_Toa= UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(20)%Rad_Toa = (ch(20)%Rad_Toa - Goes_Ch18_Thermal_Intercept)/Goes_Ch18_Thermal_Slope

   endif


   
      !---   read channel 2 (GOES Sounder channel 19)
   if (Chan_On_Flag_Default(1) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_19_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (L1b_Gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Goes_Sounder_Byte_Shift, AREAstr, &
                                    Goes_Sndr_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      Ch1_Counts = UNSIGNED_TO_REAL4(Two_Byte_Temp)
      ch(1)%Ref_Toa = Missing_Value_Real4

   endif




   !--------------------------------------------------------------------------
   ! Compute IR Brightness Temperature
   !--------------------------------------------------------------------------
   do Elem_Idx = 1,num_pix
     do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1
        if (Chan_On_Flag_Default(36) == sym%YES) then
           if (ch(36)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(36)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(36,ch(36)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(35) == sym%YES) then
           if (ch(35)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(35)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(35,ch(35)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(34) == sym%YES) then
           if (ch(34)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(34)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(34,ch(34)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif


        if (Chan_On_Flag_Default(33) == sym%YES) then
           if (ch(33)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(33)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(33,ch(33)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif


        if (Chan_On_Flag_Default(32) == sym%YES) then
           if (ch(32)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(32)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(32,ch(32)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif


        if (Chan_On_Flag_Default(31) == sym%YES) then
           if (ch(31)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(31)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(31,ch(31)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif


        if (Chan_On_Flag_Default(30) == sym%YES) then
           if (ch(30)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(30)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(30,ch(30)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif


        if (Chan_On_Flag_Default(28) == sym%YES) then
           if (ch(28)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(28)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(28,ch(28)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(27) == sym%YES) then
           if (ch(27)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(27)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(27,ch(27)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(25) == sym%YES) then
            if (ch(25)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
             ch(25)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(25,ch(25)%Rad_Toa(Elem_Idx,Line_Idx))
            endif
        endif

        if (Chan_On_Flag_Default(24) == sym%YES) then
           if (ch(24)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(24)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(24,ch(24)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(23) == sym%YES) then
            if (ch(23)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
             ch(23)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(23,ch(23)%Rad_Toa(Elem_Idx,Line_Idx))
            endif
        endif

        if (Chan_On_Flag_Default(21) == sym%YES) then
           if (ch(21)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(21)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(21,ch(21)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(20) == sym%YES) then
           if (ch(20)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(20)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(20,ch(20)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

     enddo
   enddo
   
   !------------------------------------------------------------------------------
   ! Goes Navigation
   !------------------------------------------------------------------------------
   call GET_GOES_NAVIGATION(Segment_Number, Num_Scans_Per_Segment, Num_Scans_Read, NAVstr, AREAstr, Goes_Sndr_Xstride)

   !------------------------------------------------------------------------------
   ! Goes Angles
   !------------------------------------------------------------------------------
   image_jday = jday
   image_Time_Hours = image_Time_ms / 60.0 / 60.0 / 1000.0
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1
     do Elem_Idx = 1,num_pix
        call POSSOL(image_jday,image_Time_Hours, &
                    Lon_1b(Elem_Idx,Line_Idx),Lat_1b(Elem_Idx,Line_Idx), &
                    solzen(Elem_Idx,Line_Idx),solaz(Elem_Idx,Line_Idx))
     enddo
     call COMPUTE_SATELLITE_ANGLES(goes_Sub_Satellite_longitude,  &
                      goes_Sub_Satellite_latitude, Line_Idx)
   enddo
   
   !--- scan number and time
   First_Line_In_Segment = (Segment_Number-1)*Num_Scans_Per_Segment
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1
     Scan_number(Line_Idx) = First_Line_In_Segment + Line_Idx
         ! - For now use image time for GOES Sounder
         if (Scan_Time(Line_Idx) == 0) then
           Scan_Time(Line_Idx) = image_Time_ms  !+ (Scan_number(Line_Idx)-1) * Scan_rate
         endif
   enddo


end subroutine READ_GOES_SNDR

!-----------------------------------------------------------------------------------------------------
! routine to read to AREA and NAVIGATION headers
!------------------------------------------------------------------------------------------------------
 subroutine GET_GOES_HEADERS(filename,AREAstr,NAVStr)
   character(len=*), intent(in):: filename
   type(AREA_STRUCT), intent(out):: AREAstr
   type(GVAR_NAV), intent(out):: NAVstr
   integer:: i
   integer:: recnum
   integer:: nav_offset
   integer:: a
   character(len=4):: in
   integer:: Num_Elements_This_image
   integer:: num_Scans_This_image

   open(unit=1,file=trim(filename),form="unformatted",access="direct",recl=4,status="old",action="read")

        read (unit=1, rec = 1) AREAstr%area_Status
        read (unit=1, rec = 2) AREAstr%Version_Num
        if (AREAstr%Version_Num .ne. 4) then
!          print *, "Area file cannot be read"  ! byte swapping may cause this
!          AREAstr%swap_bytes = 1
           return
        endif
        read (unit=1, rec = 3) AREAstr%Sat_Id_Num
        read (unit=1, rec = 4) AREAstr%Img_Date
        read (unit=1, rec = 5) AREAstr%Img_Time
        read (unit=1, rec = 6) AREAstr%North_Bound
        read (unit=1, rec = 7) AREAstr%West_Vis_Pixel
        read (unit=1, rec = 8) AREAstr%z_coor
        read (unit=1, rec = 9) AREAstr%Num_Line
        read (unit=1, rec = 10) AREAstr%Num_Elem
        read (unit=1, rec = 11) AREAstr%Bytes_Per_Pixel
        read (unit=1, rec = 12) AREAstr%Line_Res
        read (unit=1, rec = 13) AREAstr%Elem_Res
        read (unit=1, rec = 14) AREAstr%num_chan
        read (unit=1, rec = 15) AREAstr%num_byte_ln_prefix
        read (unit=1, rec = 16) AREAstr%proj_num
        read (unit=1, rec = 17) AREAstr%creation_date
        read (unit=1, rec = 18) AREAstr%creation_Time
        read (unit=1, rec = 19) AREAstr%sndr_filter_map
        read (unit=1, rec = 20) AREAstr%img_id_num
        do i=1,4
           read (unit=1, rec = i + 20) AREAstr%id(i)
        enddo
        do i=1,8
           read (unit=1, rec = i + 24) in
           AREAstr%comment = AREAstr%comment // in
        enddo
        read (unit=1, rec = 33) AREAstr%pri_key_calib
        read (unit=1, rec = 34) AREAstr%pri_key_nav
        read (unit=1, rec = 35) AREAstr%sec_key_nav
        read (unit=1, rec = 36) AREAstr%val_code
        do i=1,8
           read (unit=1, rec=i+36) AREAstr%pdl(i)
        enddo
        read (unit=1, rec = 45) AREAstr%band8
        read (unit=1, rec = 46) AREAstr%act_Img_Date
        read (unit=1, rec = 47) AREAstr%act_Img_Time
        read (unit=1, rec = 48) AREAstr%act_Start_Scan
        read (unit=1, rec = 49) AREAstr%len_prefix_doc
        read (unit=1, rec = 50) AREAstr%len_prefix_calib
        read (unit=1, rec = 51) AREAstr%len_prefix_lev
        do i=1,1
           read (unit=1, rec = i + 51) in
           AREAstr%src_Type = AREAstr%src_Type//in
        enddo
        do i=1,1
           read (unit=1, rec = i + 52) in
           AREAstr%calib_Type = AREAstr%calib_Type//in
        enddo
        read (unit=1, rec = 54) AREAstr%avg_or_Sample
                read (unit=1, rec = 55) AREAstr%poes_Signal
        read (unit=1, rec= 56) AREAstr%poes_up_down
        do i=1,4
           read (unit=1, rec = i + 56) in
           AREAstr%orig_Src_Type = AREAstr%orig_Src_Type//in
        enddo
        do i=1,4
           read (unit=1, rec = i + 60) AREAstr%reserved(i)
        enddo
        AREAstr%cal_offset = AREAstr%reserved(3)

        !------------------------------------------------------------
        ! read in navigation data
        !------------------------------------------------------------
        nav_offset = AREAstr%sec_key_nav/4
        read (unit=1, rec = nav_offset + 1) NAVstr%nav_Type
        read (unit=1, rec = nav_offset + 2) NAVstr%IMC_Status
        read (unit=1, rec = nav_offset + 3) NAVstr%spare1(1)
        read (unit=1, rec = nav_offset + 4) NAVstr%spare1(2)
        read (unit=1, rec = nav_offset + 5) NAVstr%stat
        read (unit=1, rec = nav_offset + 6) NAVstr%ref_long
        read (unit=1, rec = nav_offset + 7) NAVstr%ref_rad_dist
        read (unit=1, rec = nav_offset + 8) NAVstr%ref_lat
        read (unit=1, rec = nav_offset + 9) NAVstr%ref_orb_yaw
        read (unit=1, rec = nav_offset + 10) NAVstr%ref_att_roll
        read (unit=1, rec = nav_offset + 11) NAVstr%ref_att_pitch
        read (unit=1, rec = nav_offset + 12) NAVstr%ref_att_yaw
        read (unit=1, rec = nav_offset + 13) a !NAVstr%epoch_Time(1:4)
        read (unit=1, rec = nav_offset + 14) a !NAVstr%epoch_Time(5:8)
        read (unit=1, rec = nav_offset + 15) NAVstr%start_Time
        read (unit=1, rec = nav_offset + 16) NAVstr%imc_corr_roll
        read (unit=1, rec = nav_offset + 17) NAVstr%imc_corr_pitch
        read (unit=1, rec = nav_offset + 18) NAVstr%imc_corr_yaw
        do i=1,13
           read (unit=1, rec = i + nav_offset + 18) NAVstr%ref_long_change(i)
        enddo
        do i=1,11
           read (unit=1, rec = i + nav_offset + 31) NAVstr%ref_rad_dist_change(i)
        enddo
        do i=1,9
           read (unit=1, rec = i + nav_offset + 42) NAVstr%sine_lat(i)
        enddo
        do i=1,9
           read (unit=1, rec = i + nav_offset + 51) NAVstr%sine_orb_yaw(i)
        enddo
        read (unit=1, rec = nav_offset + 61) NAVstr%solar_rate
        read (unit=1, rec = nav_offset + 62) NAVstr%exp_Start_Time
        recnum = nav_offset + 62
        call READ_IMGR_RP(recnum, NAVstr%roll_att)
        do i=1,10
           recnum = recnum + 1
           read (unit=1, rec = recnum) NAVstr%spare2(i)
        enddo
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%more1
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%gvar1
        call READ_IMGR_RP(recnum, NAVstr%pitch_att)
        call READ_IMGR_RP(recnum, NAVstr%yaw_att)
        do i=1,16
           recnum = recnum + 1
           read (unit=1, rec = recnum) NAVstr%spare3(i)
        enddo
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%more2
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%gvar2
        call READ_IMGR_RP(recnum, NAVstr%roll_misalgn)
        call READ_IMGR_RP(recnum, NAVstr%pitch_misalgn)
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%Img_Date
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%Img_Time
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%instr
        do i=1,9
           recnum = recnum + 1
          read (unit=1, rec = recnum) NAVstr%spare4(i)
        enddo
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%ns_cyl
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%ew_cyl
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%ns_inc
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%ew_inc
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%more3
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%gvar3
        do i=1,126
           recnum = recnum + 1
           read (unit=1, rec = recnum) NAVstr%spare5(i)
        enddo
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%more4
        recnum = recnum + 1
        read (unit=1, rec = recnum) NAVstr%gvar4
        do i=1,127
           recnum = recnum + 1
           read (unit=1, rec = recnum) NAVstr%spare6(i)
        enddo


   close(unit=1)

   !--- compute scan rate for future use
   Num_Elements_This_image =  int(AREAstr%Num_Elem / Goes_Xstride) + 1
   num_Scans_This_image = AREAstr%Num_Line
   Scan_rate = real((Num_Elements_This_image)/ real(Num_4km_Elem_Goes_Fd/Goes_Xstride)) * &
               real((num_Scans_This_image) / real(Num_4km_Scans_Goes_Fd)) * &
               real(time_for_fd_Scan_goes) / real(num_Scans_This_image)

 ! print *, "elements = ", Num_Elements_This_image, Num_4km_Elem_Goes_Fd/Goes_Xstride
 ! print *, "lines = ", num_Scans_This_image, Num_4km_Scans_Goes_Fd
 ! print *, "time = ", time_for_fd_Scan_goes
 ! print *, "scan rate = ", Scan_rate
 end subroutine GET_GOES_HEADERS

!---------------------------------------------------------------------------
!
! words = unit of storage in these files.
! AREAstr%Bytes_Per_Pixel = number of bytes in a pixel or word
! Words_In_Prefix = number of words in header on each scan line
! Words_Per_Line = number of words in a line (similar to record length)
! First_Line_In_Segment = first data line for this segment
!---------------------------------------------------------------------------
 subroutine GET_IMAGE_FROM_AREAFILE(filename,Byte_Shift, &
                                    AREAstr, XStride, &
                                    Segment_Number, &
                                    Num_Lines_Per_Segment, &
                                    Num_Lines_Read, image)

 character(len=*), intent(in):: filename 
 integer(kind=int4), intent(in):: Byte_Shift 
 type (AREA_STRUCT), intent(in) :: AREAstr
 integer(kind=int4), intent(in):: Xstride
 integer(kind=int4), intent(in):: Segment_Number
 integer(kind=int4), intent(in):: Num_Lines_Per_Segment
 integer(kind=int4), intent(out):: Num_Lines_Read
 integer(kind=int2), dimension(:,:), intent(out):: image
 integer(kind=int4):: bytes_per_pixel !needed for HiRID
 integer(kind=int4):: bytes_per_line
 integer(kind=int4):: num_byte_ln_prefix
 integer(kind=int4):: Words_In_Prefix
 integer(kind=int4):: Words_Per_Line
 integer(kind=int4):: First_Line_In_Segment
 integer(kind=int4):: last_line_in_Segment
 integer(kind=int4):: First_Byte_In_Segment
 integer(kind=int4):: Number_Of_Words_In_Segment
 integer(kind=int4):: Number_Of_Words_Read
 integer(kind=int4):: Word_Start
 integer(kind=int4):: Word_Start_Prefix
 integer(kind=int4):: Word_End
 integer(kind=int4):: Line_Idx
 integer(kind=int2), dimension(:), allocatable:: Word_Buffer,imgbuf 
 integer(kind=int1), dimension(:), allocatable:: Word_Buffer_I1
 integer:: Nwords
 INTEGER(kind=int4), DIMENSION(64) :: i4buf_temp
 integer(kind=int4):: dummy
  
 ! get number of bytes per pixel and num bytes per line for current file
 ! this is needed because the 0.64 and other channels have different values
 ! in the MTSAT HIRID format.
  
 call mreadf_int(trim(filename)//CHAR(0),0,4,64,dummy,i4buf_temp)
 bytes_per_pixel = i4buf_temp(11)
 num_byte_ln_prefix = i4buf_temp(15)

 image = 0

 bytes_per_line = num_byte_ln_prefix + (AREAstr%Num_Elem*Bytes_Per_Pixel) 

 Words_In_Prefix = num_byte_ln_prefix / Bytes_Per_Pixel

 Words_Per_Line = Words_In_Prefix + AREAstr%Num_Elem 

 First_Line_In_Segment = (Segment_Number-1)*Num_Lines_Per_Segment + 1

 last_line_in_Segment = min(AREAstr%Num_Line,Segment_Number*Num_Lines_Per_Segment)

 First_Byte_In_Segment = AREAstr%pri_key_nav + &
                         Bytes_Per_Pixel*(First_Line_In_Segment-1) * Words_Per_Line + &
                         Bytes_Per_Pixel 

 Number_Of_Words_In_Segment = Words_Per_Line * Num_Lines_Per_Segment

 allocate(Word_Buffer(Number_Of_Words_In_Segment), imgbuf(Number_Of_Words_In_Segment))

 !--- Account for different Bytes_Per_Pixel value (Some MTSAT data is 1 byte per pixel) 
 select case (bytes_per_pixel)

   case(1)

     allocate(Word_Buffer_I1(Number_Of_Words_In_Segment))

     call mreadf_int(filename//CHAR(0), &
                 First_Byte_In_Segment,  &
                 bytes_per_pixel, &
                 Number_Of_Words_In_Segment, &
                 Number_Of_Words_Read,Word_Buffer_I1)

     Word_Buffer = int(Word_Buffer_I1,kind=int2)

     !Since 1 byte values are always signed in Fortran, add 256 to the negative values
     where (Word_Buffer < 0)
      Word_Buffer = Word_Buffer + 256
     end where
  
   case(2)

     call mreadf_int(filename//CHAR(0), &
                 First_Byte_In_Segment,  &
                 bytes_per_pixel, &
                 Number_Of_Words_In_Segment, &
                 Number_Of_Words_Read,Word_Buffer)

   case default

     print *, EXE_PROMPT, "Unsupported Bytes_Per_Pixel Value in GET_IMAGE_FROM AREAFILE, stopping"
     stop

 end select

 !--- update number of scans read
 Num_Lines_Read = Number_Of_Words_Read / Words_Per_Line

 do Line_Idx = 1, Num_Lines_Read
     Word_Start = (Words_In_Prefix + AREAstr%Num_Elem)*(Line_Idx-1) + Words_In_Prefix + 1
!    Word_End = min(Word_Start + AREAstr%Num_Elem,Number_Of_Words_In_Segment)
     Word_End = min(Word_Start + (AREAstr%Num_Elem-1),Number_Of_Words_In_Segment)
     Nwords = int(Word_End - Word_Start)/Xstride + 1
     image(1:Nwords,Line_Idx) = ishft(Word_Buffer(Word_Start:Word_End:Xstride),Byte_Shift)

    Word_Start_Prefix = (Words_In_Prefix + AREAstr%Num_Elem)*(Line_Idx-1) + 1

    if (allocated(Scan_Time)) then
      if (Words_In_Prefix .LT. 128) then
        Scan_Time(Line_Idx) = 0
      else
        call PRINT_PREFIX(Word_Buffer(Word_Start_Prefix:Word_Start_Prefix+Words_In_Prefix), &
                        Scan_Time(Line_Idx))
      end if
    endif

 enddo

 !--- deallocate allocated arrays
 deallocate(Word_Buffer, imgbuf)
 if (allocated(Word_Buffer_I1)) deallocate(Word_Buffer_I1)


 end subroutine GET_IMAGE_FROM_AREAFILE

!------------------------------------------------------------------------------
! Goes Navigation
!------------------------------------------------------------------------------
subroutine  GET_GOES_NAVIGATION(Segment_Number, Num_Lines_Per_Segment, &
                                Num_Lines_Read,NAVstr,AREAstr, Xstride)

 type (GVAR_NAV), intent(in) :: NAVstr
 type (AREA_STRUCT), intent(in) ::AREAstr
 integer(kind=int4), intent(in):: Segment_Number
 integer(kind=int4), intent(in):: Num_Lines_Per_Segment
 integer(kind=int4), intent(in):: Num_Lines_Read
 integer(kind=int4), intent(in):: Xstride

 real (kind=real8) :: Line, Elem, Elev, Scan
 real (kind=real8) :: Dlat, Dlon
 integer(kind=int4):: Ierr
 integer(kind=int4):: Elem_Idx_Temp
 integer(kind=int4):: Line_Idx_Temp
 integer(kind=int4):: Elem_Idx
 integer(kind=int4):: Line_Idx
 integer(kind=int4):: First_Line_In_Segment

 !--- initialize
 Lat_1b = Missing_Value_Real4
 Lon_1b = Missing_Value_Real4
 Space_Mask = sym%SPACE

 !-- compute first line in data space for this segment
 First_Line_In_Segment = (Segment_Number-1)*Num_Lines_Per_Segment + 1

 !-- compute first line in satellite space
 Line_Idx_Temp = 1 + (First_Line_In_Segment - 1)*AREAstr%Line_Res

 !--- loop over all lines in segment
 do Line_Idx = 1,Num_Lines_Read

!     Line = real(AREAstr%North_Bound) + real(Line_Idx_Temp - 1) + &
!            real(AREAstr%Line_Res)/2.0

      if(Goes_Sndr_Flag == sym%YES) then 
         line = (real(AREAstr%north_bound) + real(Line_Idx_Temp - 1) + 9.0) / 10.0
      else
         line = real(AREAstr%north_bound) + real(Line_Idx_Temp - 1) + &
              real(AREAstr%line_res)/2.0
      endif

      !--- loop over all elements
      do Elem_Idx = 1, Num_Pix

         !-- element index in satellite space
         Elem_Idx_Temp = (Elem_Idx - 1)*AREAstr%Elem_Res*Xstride + 1

!         Elem = real(AREAstr%West_Vis_Pixel) + real(Elem_Idx_Temp - 1) + &
!                real(AREAstr%Elem_Res*(Xstride))/2.0

        ! For sounder, convert McIDAS elem to GVAR elem. 
        ! Copied from McIDAS code, nvxgvar.dlm, v1.16
        if(Goes_Sndr_Flag == sym%YES) then 
           elem = (real(AREAstr%west_vis_pixel) + real(Elem_Idx_Temp - 1) + 9.0) / 10.0  
        else
           elem = real(AREAstr%west_vis_pixel) + real(Elem_Idx_Temp - 1) + &
                real(AREAstr%elem_res*(xstride))/2.0
        end if
        
        
         CALL COMP_ES(NAVstr, Line, Elem, Elev, Scan)

         !--- Convert angles to lat/lon
         CALL LPOINT(NAVstr%instr, 1, Elev, Scan, Dlat, Dlon, Ierr)

         Space_Mask(Elem_Idx,Line_Idx) = sym%SPACE

         if (Ierr == 0) then
          Lat_1b(Elem_Idx,Line_Idx) = Dlat / DTOR
          Lon_1b(Elem_Idx,Line_Idx) = Dlon / DTOR
          Space_Mask(Elem_Idx,Line_Idx) = sym%NO_SPACE
         endif

      end do

      !--- increment line index in satellite space
      Line_Idx_Temp = Line_Idx_Temp + AREAstr%Line_Res

    end do

end subroutine GET_GOES_NAVIGATION


!-----------------------------------------------------------------------
 subroutine READ_IMGR_RP(recnum, struct)
        integer, intent(inout):: recnum
        type(IMGR_RP), intent(out):: struct
        integer:: i

        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%exp_mag
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%exp_Time_const
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%mean_att_ang_const
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%num_Sinu_per_angle
        do i=1,15
           call READ_IMGR_SIN(recnum, struct%sinusoid(i))
        enddo
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%num_mono_Sinu
        do i=1,4
           call READ_IMGR_MON(recnum, struct%monomial(i))
        enddo
 end subroutine READ_IMGR_RP

!-----------------------------------------------------------------------
  subroutine READ_IMGR_SIN(recnum, struct)
    integer, intent(inout):: recnum
    type(imgr_Sin), intent(out):: struct
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%mag_Sinu
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%phase_ang_Sinu
  end subroutine READ_IMGR_SIN
!-----------------------------------------------------------------------
  subroutine READ_IMGR_MON(recnum, struct)
    integer, intent(inout):: recnum
     type(imgr_mon), intent(out):: struct
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%order_appl_Sinu
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%order_mono_Sinu
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%mag_mono_Sinu
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%phase_ang_Sinu
        recnum = recnum + 1
        read (unit=1, rec = recnum) struct%ang_from_epoch
   end subroutine READ_IMGR_MON

!-----------------------------------------------------------------------
        subroutine FIX_VAR(var)
        integer, intent(inout):: var

        var=ishftc(ishftc(ishftc(var,-8,32),-8,24),-8,16)

        end subroutine FIX_VAR
!-----------------------------------------------------------------------
        subroutine FIX_VAR2(var)
!       integer(kind=int2), intent(inout):: var
        integer*2, intent(inout):: var

        var=ishftc(var,-8,16)

        end subroutine FIX_VAR2
!======================================================================
! BEGIN CIRA GOES IMAGER ROUTINES
!======================================================================
!========================================================================
!     I N S T 2 E
!========================================================================
      subroutine INST2E( R, P, Y, A, AT )
!
!     AUTHOR:         KELLY DEAN
!
!     CREATED:        October 1994
!   
!     DEVELOPED FOR:  CIRA/COLORAdo STATE UNIVERSITY
!
!     PURPOSE:        
!	Procedure INST2E accepts the single precision roll, pitch and yaw
!    angles of an instrument and returns the double precision instrument
!    to earth coordinates transformation matrix.
!
!     REVISION:       0.0
!
!     REFERENCES:
!                     OTHER DOCUMENTS
!
!     COMMENTS:
!               Adapted from Igor Levine program for Integral Systems, Inc.
!
!     ARGUMENTS:
!       NAME:   type:       PURPOSE:                        IN/OUT:
!	 R      real*8	    Roll angle (rad)		     IN
!	 P      real*8      Pitch angle (rad)		     IN
!	 Y      real*8	    Yaw angle (rad)		     IN
!	 A      real*8	    Spacecraft to ECEF coordinates  OUT
!				 transformation matrix
!	 AT     real*8	    Instrument to ECEF coordinates  OUT
!				 transformation matrix
!    
!     VARIABLES:
!       NAME:    PURPOSE:
!       ****************  integer    *****************
!	I    Indices
!	J    Indices
!       ****************  real*8     *****************
!	RPY  Instrument to body coordinates transformation matrix
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     VARIABLE DECLARATION SECTION:
      real*8  A(3,3),AT(3,3),R,RPY(3,3),P,Y
      integer*4   I,J
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
!     Compute instrument to body coordinates transformation matrix
!     by using a small angle approximation of  trigonometric function
!     of the roll, pitch and yaw.
!
      RPY(1,1) = 1.0D0 - 0.5D0 * ( P * P + Y * Y )
      RPY(1,2) = -Y
      RPY(1,3) = P
      RPY(2,1) = Y + P * R
      RPY(2,2) = 1.0D0 - 0.5D0 * ( Y * Y + R * R )
      RPY(2,3) = -R
      RPY(3,1) = -P + R * Y
      RPY(3,2) = R + P * Y
      RPY(3,3) = 1.0D0 - 0.5D0 * ( P * P + R * R )
!
!     Multiplication of matrices A and RPY
!
      do I = 1,3
       do J = 1,3
        AT(I,J) = A(I,1)*RPY(1,J)+A(I,2)*RPY(2,J)+A(I,3)*RPY(3,J)
       ENDDO
      ENDDO
!
!
      return
      end subroutine INST2E
!======================================================================
!     T I M E 5 0
!======================================================================
      real*8 FUNCTION TIME50(btim)
!
!     AUTHOR:         Garrett Campbell and Kelly Dean
!
!     CREATED:        October 1994
!   
!     DEVELOPED FOR:  CIRA/COLORAdo STATE UNIVERSITY
!
!     PURPOSE:        
!	Function TIME50 will take the epoch time from the GVAR NAVstr%and
!      convert it to Minutes from January 1, 1950.  NOTE - Epoch time in
!      the NAVstr%is not in the same format as other BCD times.C
!
!     REVISION:       0.0
!
!     ARGUMENTS:
!       NAME:   type:       PURPOSE:                          IN/OUT:
!	btim     BYTE        Binary coded data (BCD) time      IN
!    
!     FUNCTIONS:
!       NAME:   type:       PURPOSE:                        LIBRARY:
!       MOD     integer     Returns a remainder              Intrinsic
!
!       NAME:    PURPOSE:
!       ****************  integer    *****************
!       day_100   Part of day extracted from BCD
!       day_10    Part of day extracted from BCD
!       day_1     Part of day extracted from BCD
!       Hour_10   Part of Hour extracted from BCD
!       Hour_1    Part of Hour extracted from BCD
!       min_10    Part of Minute extracted from BCD
!       min_1     Part of Minute extracted from BCD
!       NY        YEAR
!       ND        DAY OF YEAR
!       NH        HOUR
!       NM        MINUTE
!       ibt
!       J         Loop control variable
!       year_1000 Part of year extracted from BCD
!       year_100  Part of year extracted from BCD
!       year_10   Part of year extracted from BCD
!       year_1    Part of year extracted from BCD
!       ****************  real*8     *****************
!       S        SECONDS - Double precision
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     VARIABLE DECLARATION SECTION:
!
      byte btim(8),bt
      integer NY,ND,NH,NM,J
      integer year_1000, year_100, year_10, year_1
      integer day_100, day_10, day_1, Hour_10, Hour_1, min_10, min_1
      integer sec_10, sec_1, msec_100, msec_10, msec_1
      integer ibt
      real flywheel
      real*8 S
!
!     Equivalence DECLARATION SECTION:
!
      equivalence (bt,ibt)
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
!    Extract the Binary Coded Time into separate year, Julian day, Hour
!    Minutes, seconds.
!
!     bt = btim(1)
      day_1   = ibt/16
      Hour_10 = mod(ibt,16)
      bt = btim(2)
      day_100  = mod(ibt/16,8)
      day_10   = mod(ibt,16)
      flywheel = mod(ibt,2)
      bt = btim(3)
      year_10   = ibt/16
      year_1   = mod(ibt,16)
      bt = btim(4)
      year_1000 = ibt/16
      year_100  = mod(ibt,16)
      bt = btim(5)
      msec_10 = ibt/16
      msec_1  = mod(ibt,16)
      bt = btim(6)
      sec_1    = ibt/16
      msec_100 = mod(ibt,16)
      bt = btim(7)
      min_1  = ibt/16
      sec_10 = mod(ibt,16)
      bt = btim(8)
      Hour_1 = ibt/16
      min_10 = mod(ibt,16)
!
!     Make the year, Julian day, Hour, Minute, and seconds.
!
      ny = year_1000 * 1000 + year_100 * 100 + year_10 * 10 + year_1
      nd = day_100 * 100 + day_10 * 10 + day_1 
      nh = Hour_10 * 10 + Hour_1
      nm = min_10 * 10 + min_1
      s  = sec_10 * 10.0D0 + sec_1 +                               &
          msec_100 * 0.1D0 + msec_10 * 0.01D0 + msec_1 * 0.001D0

!
!     HERE WE CONVERT integer YEAR AND DAY OF YEAR TO NUMBER OF                 
!     DAYS FROM 0 HOUR UT, 1950 JAN. 1.0                                        
!     THIS CONVERTION IS BASED ON AN ALGORITHM BY FLIEGEL AND VAN               
!     FLANDERN, COMM. OF ACM, VOL.11, NO. 10, OCT. 1968 (P.657)                 
!
      j = nd + 1461 * (ny + 4799) / 4 - 3 *     &
         ( ( ny + 4899 ) / 100 ) / 4 - 2465022
!
!    Compute time in Minutes from January 1.0, 1950 as double precision.
!
      TIME50 = j * 1440.0D0 + nh * 60.0D0 + nm + s / 60.0D0
!
!
      return
      end FUNCTION TIME50
!=======================================================================
!    L M O D E L
!=======================================================================
      subroutine LMODEL( T, TU, NAVstr, RLAT, RLON )
!
!     AUTHOR:         KELLY DEAN
!
!     CREATED:        October 1994
!   
!     DEVELOPED FOR:  CIRA/COLORAdo STATE UNIVERSITY
!
!     PURPOSE:        
!	Procedure LModel accepts an input time and a set of O&A parameters
!	and computes position of the satellite, the attitude angles and
!	attitudes misalignment and the instrument to earth fixed coordinates
!	transformation matrix.
!
!	This procedure computes the position of the satellite and the
!	attitude of the imager or sounder. The calculations are based
!	on the Oats orbit and attitude model represented by the O&A
!	parameter set in NAVstr%
!
!     REVISION:       0.0
!
!     REFERENCES:
!                 Part of this code was adapted from Igor Levine work
!            for Integal System, Inc.
!
!     ARGUMENTS:
!       NAME:   type:       PURPOSE:                                IN/OUT:
!	T        real*8      Input time from Jan 1, 1950 (Minutes)   IN
!	TU       real*8      Epoch time from Jan 1, 1950 (Minutes)   IN
!	RLAT     real*8	     Subsatellite Geodetic latitude (rad)    OUT
!	RLON     real*8	     Subsatellite Geodetic Longitude (rad)   OUT
!    
!     subroutineS:
!       NAME:   PURPOSE:                                    LIBRARY:
!        INST2E   Computes instrument to earth coordinates      GVARnav
!
!     FUNCTIONS:
!       NAME:   type:       PURPOSE:                        LIBRARY:
!        DATAN    real*8      Arc tangent (double precision)        Intrinsic
!        DATAN2   real*8      Arc Tangent (double precision)        Intrinsic
!        DCOS     real*8      Cosine (double precision)             Intrinsic
!        DSIN     real*8      Sine (double precision)               Intrinsic
!        DTAN     real*8      tagent (double precision)             Intrinsic
!        GATT     real*8      Compute attitude and misalignment angle  GVARnav
!
!     VARIABLES:
!       NAME:    PURPOSE:
!       ****************  real*8       *****************
!       XS     NORMALIZED S/C POSITION IN ECEF COORDINATES              
!       BT    ECEF TO INSTRUMENT COORDINATES TRANSFORMATION            
!       Q3    USED IN subroutine LPOINT                                
!       PITCH PITCH ANGLES OF INSTRUMENT (RAD)            
!       ROLL  ROLL ANGLES OF INSTRUMENT (RAD)            
!       YAW   YAW ANGLES OF INSTRUMENT (RAD)             
!       PMA  PITCH MISALIGNMENTS OF INSTRUMENT (RAD)         
!       RMA  ROLL MISALIGNMENTS OF INSTRUMENT (RAD)         
!	R    Normalized satellite distance (km)
!	TS   Time from EPOCH (Minutes)
!	B    Spacecraft to earth fixed coordinates transmation matrix
!	TE   Exponential time delay from EPOCH (Minutes)
!	PHI  Subsatellite geocentric latitude (rad)
!	DR   Radial distance from the nominal (km)
!	PSI  Orbital yaw (rad)
!	LAM  IMC longitude (rad)
!	U    Argument of latitude (rad)
!	SU   DSIN(U)
!	CU   DCOS(U)                                             
!	SI   Sine of the orbit inclination
!	CI   Cosine of the orbit inclination
!	SLAT Sine of geocentric latitude
!	ASC  Longitude of the ascending node (rad)
!	SA   Sine of ASC
!	CA   Cosine of ASC
!	SYAW Sine of the orbit yaw
!	WA   Solar orbit angle (rad)
!	W    Orbit angle (rad)
!	SW   DSIN(W)
!	CW   DCOS(W)
!	S2W  DSIN(2*W)
!	C2W  DCOS(2*W)
!	SW1  DSIN(0.927*W)
!	CW1  DCOS(0.927*W)
!	SW3  Sine of 1.9268*W
!	CW3  Cosine of 1.9268*W
!	DLAT Change in sine of geocentric latitude
!	DYAW Change in sine of orbit yaw
!	A1   Work area
!	A2   Work area
!	XS   S/C position in ECEF coordinates
!
!     COMMON BLOCKS:
!       NAME:     CONTENTS:
!       ELCOMM     Instrument position and attitude variables and 
!                   transformation matrix
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     CONSTANT DECLARATION SECTION:
!
!     include 'kdf.inc'
!
!     VARIABLE DECLARATION SECTION:
!
      integer IMCstatus
      real*8 T, TU
!     real*8 REC(336)
      real*8 PI, DEG, RAD, NOMORB, AE, FER, AEBE2, AEBE3, AEBE4
      real*8 RLAT, RLON, R, TS, TE, PHI, DR, PSI, LAM, U, SU, CU
      real*8 SI, CI, SLAT, ASC, SA, CA, SYAW, WA, W, SW, CW, S2W, C2W
!     real*8 SW1, CW1, SW3, CW3, DLAT, DYAW, A1, A2
      real*8 SW1, CW1, SW3, CW3, DLAT, DYAW
      real*8 B(3,3), BT(3,3), XS(3)
      real*8 Q3, PITCH, ROLL, YAW
      real*8 PMA, RMA
      type(GVAR_NAV):: NAVstr
!
!     FUNCTION DECLARATION SECTION:
!
      real*8 DATAN,DATAN2,DCOS,DSIN,DTAN!,GATT
!
!     COMMON BLOCKS:
!
      COMMON /ELCOMM/ xs, bt, q3, pitch, roll, yaw, pma, rma 
!
!     INITIALIZATIONS: (Description mathematical and earth-related constants)
!
!      PI     = 3.141592653589793D0
      PI     = 4.0*DATAN(1.0D0)
      DEG    = 180D0 / PI
      RAD    = PI / 180D0  ! Degrees to radians conversion (PI/180)
      NOMORB = 42164.365D0 !  Nominal radial distance of satellite (km)
      AE     = 6378.137D0  !  Earth equatorial radius (km)
      FER    = 1.0D0 - ( 6356.7533D0 / AE )  ! Earth flattening coefficient 
      AEBE2  = 1.0D0 / (1.0D0 - FER )**2 
      AEBE3  = AEBE2 - 1. 
      AEBE4  = ( 1.0D0 - FER )**4-1.
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
!     Determine the IMC status
! 
      IMCstatus = IBITS(NAVstr%stat,6,1)  
!
!     Assign referenec values to the subsatellite longitude and
!     latitude, the radial distance and the orbit yaw.
!                                                                               
      LAM = NAVstr%ref_long * 1.0D-7
      DR  = NAVstr%ref_rad_dist
      PHI = NAVstr%ref_lat
      PSI = NAVstr%ref_orb_yaw
!
!     Assign reference values to the attitudes and misalignments
!
      ROLL  = NAVstr%ref_att_roll
      PITCH = NAVstr%ref_att_pitch
      YAW   = NAVstr%ref_att_yaw
      RMA   = 0.0
      PMA   = 0.0
!
!     if IMC_active is OFF, compute changes in the satellite orbit
!
      if ( IMCstatus .NE. 0 ) then
!	PRINT *, ' IMC turned off............'
!
!       Compute time since EPOCH (Minutes)
!
        TS = T - TU
!
!       Compute orbite angle and the related trigonometric functions.
!       earth rotational rate (.729115E-4 rad/sec).
!
        W   = 0.729115e-4 * 60.0D0 * TS
        SW  = DSIN(W)
        CW  = DCOS(W)
        SW1 = DSIN(0.927D0*W)
        CW1 = DCOS(0.927D0*W)
        S2W = DSIN(2.0D0*W)
        C2W = DCOS(2.0D0*W)
        SW3 = DSIN(1.9268D0*W)
        CW3 = DCOS(1.9268D0*W)
!
!     Computes change in the IMC_active longitude from the reference.
!
        LAM = LAM + ( NAVstr%ref_long_change(1) * 1.0D-7 ) +             &
      	    ( ( NAVstr%ref_long_change(2) * 1.0D-7 ) +             &
      	      ( NAVstr%ref_long_change(3) * 1.0D-7 ) * W ) * W +             &
              ( NAVstr%ref_long_change(10) * 1.0D-7 ) * SW1 +             &
      	      ( NAVstr%ref_long_change(11) * 1.0D-7 ) * CW1 +             &
      	    ( ( NAVstr%ref_long_change(4) * 1.0D-7 )  * SW  +             &
      	      ( NAVstr%ref_long_change(5) * 1.0D-7 ) * CW +             &
              ( NAVstr%ref_long_change(6) * 1.0D-7 ) * S2W +             &
      	      ( NAVstr%ref_long_change(7) * 1.0D-7 ) * C2W +             &
      	      ( NAVstr%ref_long_change(8) * 1.0D-7 ) * SW3 +             &
      	      ( NAVstr%ref_long_change(9) * 1.0D-7 ) * CW3 +             &
              W * ( ( NAVstr%ref_long_change(12) * 1.0D-7 ) * SW +             &
      	      ( NAVstr%ref_long_change(13) * 1.0D-7 ) * CW ) ) * 2.0D0
!
!       Computes change in radial distance from the reference (km)
!
        DR = DR + ( NAVstr%ref_rad_dist_change(1) * 1.0D-7 ) +             &
      		  ( NAVstr%ref_rad_dist_change(2) * 1.0D-7 ) * CW  +             &
      		  ( NAVstr%ref_rad_dist_change(3) * 1.0D-7 ) * SW  +             &
           	  ( NAVstr%ref_rad_dist_change(4) * 1.0D-7 ) * C2W +             &
      		  ( NAVstr%ref_rad_dist_change(5) * 1.0D-7 ) * S2W +             &
      		  ( NAVstr%ref_rad_dist_change(6) * 1.0D-7 ) * CW3 +             &
      		  ( NAVstr%ref_rad_dist_change(7) * 1.0D-7 ) * SW3 +             &
      		  ( NAVstr%ref_rad_dist_change(8) * 1.0D-7 ) * CW1 +             &
      		  ( NAVstr%ref_rad_dist_change(9) * 1.0D-7 ) * SW1 +             &
      		  W * ( ( NAVstr%ref_rad_dist_change(10) * 1.0D-7 ) * CW +             &
      		  ( NAVstr%ref_rad_dist_change(11) * 1.0D-7 ) * SW )
!
!       Computes the sine of the change in the geocentric latitude.
!
        DLAT = ( NAVstr%sine_lat(1) * 1.0D-7 ) +             &
      	       ( NAVstr%sine_lat(2) * 1.0D-7 ) * CW  +             &
      	       ( NAVstr%sine_lat(3) * 1.0D-7 ) * SW  +             &
      	       ( NAVstr%sine_lat(4) * 1.0D-7 ) * C2W +             &
      	       ( NAVstr%sine_lat(5) * 1.0D-7 ) * S2W +             &
      	       W * ( ( NAVstr%sine_lat(6) * 1.0D-7 ) * CW +             &
      	       ( NAVstr%sine_lat(7) * 1.0D-7 ) * SW ) +             &
      	       ( NAVstr%sine_lat(8) * 1.0D-7 ) * CW1 +             &
      	       ( NAVstr%sine_lat(9) * 1.0D-7 ) * SW1
!
!	Computes geocentric latitude by using an expansion for arcsine.
!
        PHI = PHI + DLAT * ( 1.0D0 + DLAT * DLAT / 6.0D0 )
!
!	Computes sine of the change in the orbit yaw.
!
        DYAW = ( NAVstr%sine_orb_yaw(1) * 1.0D-7 ) +             &
      	       ( NAVstr%sine_orb_yaw(2) * 1.0D-7 ) * SW  +             &
      	       ( NAVstr%sine_orb_yaw(3) * 1.0D-7 ) * CW  +             &
      	       ( NAVstr%sine_orb_yaw(4) * 1.0D-7 ) * S2W +             &
      	       ( NAVstr%sine_orb_yaw(5) * 1.0D-7 ) * C2W +             &
      	       W * ( ( NAVstr%sine_orb_yaw(6) * 1.0D-7 ) * SW +             &
      	       ( NAVstr%sine_orb_yaw(7) * 1.0D-7 ) * CW ) +             &
      	       ( NAVstr%sine_orb_yaw(8) * 1.0D-7 ) * SW1 +             &
      	       ( NAVstr%sine_orb_yaw(9) * 1.0D-7 ) * CW1
!
!	Computes the orbit yaw by using an expansion for arcsine.
!
        PSI = PSI + DYAW * ( 1.0D0 + DYAW * DYAW / 6.0D0 )
!      else
!C         WRITE(6,*) ' IMC is turned on .......... >',IMCstatus
      endif
!
!     Conversion of the IMC_active longitude and orbit yaw to the subsatellite
!     longitude and the orbit inclination (REF: GOES-PCC-TM-2473). Inputs
!     required for earth location and gridding 
!
      SLAT = DSIN(PHI)
      SYAW = DSIN(PSI)

      SI = SLAT**2 + SYAW**2
      CI = DSQRT(1.0D0 - SI )
      SI = DSQRT(SI)
      if ( SYAW .NE. 0.0D0 ) then
        U = DATAN2(SLAT,SYAW)
      else if (SLAT .GT. 0.0D0 ) then
        U = 1.570796D0
      else if (SLAT .LT. 0.0D0 ) then
        U = 4.712389D0
      else 
        U = LAM
      endif
!
      SU = DSIN(U)
      CU = DCOS(U)
!
!     Computes longitude of the ascending node.
!
      ASC = LAM - U
      SA = DSIN(ASC)
      CA = DCOS(ASC)
!
!     Computes the subsatellite geographic latitude (rad)
!
      RLAT = DATAN(AEBE2*DTAN(PHI))
!
!     Computes the subsatellite geographic longitude (rad)
!
      RLON = ASC + DATAN2(CI*SU,CU)
!
!     Computes the spacecraft to earth fixed coordinates transformation matrix.
!
!         (VECTOR IN ECEF COORDINATES) = B * (VECTOR IN S/C COORDINATES)
!
      B(1,2) = -SA * SI
      B(2,2) = CA * SI
      B(3,2) = -CI
      B(1,3) = -CA * CU + SA * SU * CI
      B(2,3) = -SA * CU - CA * SU * CI
      B(3,3) = -SLAT
      B(1,1) = -CA * SU - SA * CU * CI
      B(2,1) = -SA * SU + CA * CU * CI
      B(3,1) = CU * SI
!
!     Computes the normalized spacecraft position vector in earth fixed
!     coordinates - XS.
!
      R = (NOMORB + DR) / AE
      XS(1) = -B(1,3) * R
      XS(2) = -B(2,3) * R
      XS(3) = -B(3,3) * R
!
!     Precomputes Q3 ( Used in LPoint ).
!
      Q3 = XS(1)**2 + XS(2)**2 + AEBE2 * XS(3)**2 - 1.0D0
!
!     Computes the attitudes and misalignments if IMC_active is OFF
!
      if ( IMCstatus .NE. 0 ) then
!	PRINT *, ' IMC turned off............'
!
!     Computes the solar orbit angle
!                                                                               
         WA = ( NAVstr%solar_rate * 1.0D-7 ) * TS
!
!     Computes the difference between current time, TS, and the
!     exponential time. Note that both times are since EPOCH.
!
         TE = TS - ( NAVstr%exp_Start_Time * 1.0D-7 )
!
!     Computes ROLL + ROLL Misalignment
!
         ROLL = ROLL + GATT(NAVstr%roll_att,WA,TE)
!
!     Computes Pitch + Pitch Misalignment
!                                                                               
         PITCH = PITCH + GATT(NAVstr%pitch_att,WA,TE)
!
!     Computes YAW
!
         YAW = YAW + GATT(NAVstr%yaw_att,WA,TE)
!
!     Computes roll misalignment
!                                                                               
         RMA = GATT(NAVstr%roll_misalgn,WA,TE)
!
!     Computes pitch misalignment
!
         PMA = GATT(NAVstr%pitch_misalgn,WA,TE)
!
!     Apply the Earth Sensor compensation if needed.
!
         if ( TS .GE. ( NAVstr%start_Time * 1.0D-2 ) ) then
                ROLL = ROLL + NAVstr%IMC_corr_roll * 1.0D-7
                PITCH = PITCH + NAVstr%IMC_corr_pitch * 1.0D-7
                YAW = YAW + NAVstr%IMC_corr_yaw * 1.0D-7
         endif
!      else
!C         WRITE(6,*) ' IMC is turned on .......... >',IMCstatus
      endif
!
!     Computes the instrument to earth fixed coordinates transformation
!     matrix - BT
!
      CALL INST2E( ROLL, PITCH, YAW, B, BT)
!
      return

      CONTAINS
 
!=======================================================================
!     G A T T
!=======================================================================
      real*8 FUNCTION GATT( IMGR_REP, WA, TE)
!
!     AUTHOR:         KELLY DEAN
!
!     CREATED:        October 1994
!   
!     DEVELOPED FOR:  CIRA/COLORAdo STATE UNIVERSITY
!
!     PURPOSE:        
!	   This function computes an attitude/misalignment angle from
!	 a given subset of the O&A parameters.
!
!     REVISION:       0.0
!
!     ARGUMENTS:
!       NAME:   type:       PURPOSE:                        IN/OUT:
!        IMGR_REP  STRUCTURE
!        TE        real*8     Input exponential time          IN 
!                             delay from epoch (Minutes)
!	 WA	   real*8     Input solar orbit angle (rad)   IN
!    
!     FUNCTIONS:
!       NAME:   type:       PURPOSE:                        LIBRARY:
!        DCOS    real*8      Cosine ( Double precision )      INTRINSIC
!
!     VARIABLES:
!       NAME:    PURPOSE:
!       ****************  integer    *****************
!       l      Loop control variable
!       m      Temporary variable for order of monomial sinusoids
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     INCLUDE DECLARATION SECTION:
!
!     include 'kdf.inc'
!
!     VARIABLE DECLARATION SECTION:
!
      integer l,m
      real*8  TE, WA
      type (IMGR_RP) :: IMGR_REP
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
      GATT = IMGR_REP%mean_att_ang_const * 1.0D-7
!
!	Computes the exponential term.
!
      if ( TE .GE. 0.0D0 ) then 
       GATT = GATT + (IMGR_REP%exp_mag * 1.0D-7 ) *     &
            EXP(-te / ( IMGR_REP%exp_Time_const * 1.0D-2 ))
      endif
!
!	Calculation of sinusoids.
!
      do l = 1, IMGR_REP%num_Sinu_per_angle
         GATT = GATT + ( IMGR_REP%sinusoid(l)%mag_Sinu * 1.0D-7 ) *    &
               DCOS(wa * l +                                           &
               ( IMGR_REP%sinusoid(l)%phase_ang_Sinu * 1.0D-7 ) )
      ENDDO
!
!	Computes monomial sinusoids.
!
      do l = 1, IMGR_REP%num_mono_Sinu
          m = IMGR_REP%monomial(l)%order_mono_Sinu
          GATT = GATT + (IMGR_REP%monomial(l)%mag_mono_Sinu * 1.0D-7) *    &
            (wa-(IMGR_REP%monomial(l)%ang_from_epoch * 1.0D-7) )**m *    &
            DCOS( IMGR_REP%monomial(l)%order_appl_Sinu * wa +     &
            ( IMGR_REP%monomial(l)%phase_ang_Sinu * 1.0D-7 ) )
       ENDDO
!
!
       return
       end FUNCTION GATT

      end subroutine LMODEL
!***********************************************************************
!***********************************************************************
!**
!** INTEGRAL SYSTEMS, INC.
!**
!***********************************************************************
!**
!** PROJECT : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT
!** SYSTEM : EARTH LOCATION USERS GUIDE
!** ROUTINE : LPOINT
!** SOURCE : LPOINT.FOR
!** LOAD NAME : ANY
!** PROGRAMMER: IGOR LEVINE
!**
!** VER. DATA BY COMMENT
!** ---- -------- --- ---------------------------------------------
!** 1 01/09/89 IL INITIAL CREATION
!** 2 06/02/89 IL COORDINATE AXES CHANGED ACCORDING TO
!** FORD'S DEFINITION IN SDAIP, DRL504-01
!** 3 12/01/93 IL IMPLEMENTED NEW FORMULAE FOR SCAN ANGLE
!** CORRECTIONS DUE TO MISALIGNMENTS
!**
!***********************************************************************
!**
!** THIS SUBROUTINE CONVERTS THE INSTRUMENT ELEVATION AND SCAN
!** ANGLES TO THE RELATED GEOGRAPHIC LATITUDE AND LONGITUDE.
!**
!***********************************************************************
!**
!** CALLED BY : ANY
!** COMMONS MODIFIED: NONE
!** INPUTS : NONE
!** OUTPUTS : NONE
!** ROUTINES CALLED : NONE
!**
!***********************************************************************
!***********************************************************************
      SUBROUTINE LPOINT(INSTR,FLIP_FLG,ALPHA0,ZETA0,RLAT,RLON,IERR)
      IMPLICIT NONE
!
! CALLING PARAMETERS
!
      INTEGER*4 INSTR
! INSTRUMENT CODE (1=IMAGER,2=SOUNDER)
      INTEGER*4 FLIP_FLG
! S/C ORIENTATION FLAG (1=NORMAL,-1=INVERTED)
      REAL*8 ALPHA0
! ELEVATION ANGLE (RAD)
      REAL*8 ZETA0
! SCAN ANGLE (RAD)
      REAL*8 RLAT
! LATITUDE IN RADIANS (OUTPUT)
      REAL*8 RLON
! LONGITUDE IN RADIANS (OUTPUT)
      INTEGER IERR
! OUTPUT STATUS; 0 - POINT ON THE EARTH
! 1 - INSTRUMENT POINTS OFF EARTH
!
! LOCAL VARIABLES
!
      REAL*8 G1(3)
! POINTING VECTOR IN EARTH-CENTERED COORDINATES
      REAL*8 H
! SLANT DISTANCE TO THE EARTH POINT (KM)
      REAL*8 Q1,Q2,D
! WORK SPACE
      REAL*4 G(3)
! POINTING VECTOR IN INSTRUMENT COORDINATES
      REAL*4 U(3)
! COORDINATES OF THE EARTH POINT (KM)
      REAL*4 SA,CA,DA,DZ,D1,CZ,SZ,FF,DOFF,ALPHA,ZETA
! WORK SPACE
!
! INCLUDE FILES
!
!     INCLUDE 'ELCONS.INC'
!***********************************************************************
!***********************************************************************
!**
!** INTEGRAL SYSTEMS, INC.
!**
!***********************************************************************
!**
!** PROJECT : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT
!** SYSTEM : EARTH LOCATION USERS GUIDE
!** NAME : ELCONS
!** TYPE : DATA AREA
!** SOURCE : ELCONS.INC
!**
!** VER. DATA BY COMMENT
!** ---- -------- ---------------- --------------------------------
!** A 01/09/89 I. LEVINE INITIAL CREATION
!**
!***********************************************************************
!**
!** DESCRIPTION
!** MATHEMATICAL AND EARTH-RELATED CONSTANTS
!**
!***********************************************************************
!***********************************************************************
!
      REAL*8 PI
      PARAMETER (PI=3.141592653589793D0)
      REAL*8 DEG
      PARAMETER (DEG=180.D0/PI)
      REAL*8 RAD
      PARAMETER (RAD=PI/180.D0)
! DEGREES TO RADIANS CONVERSION PI/180
      REAL*8 NOMORB
      PARAMETER (NOMORB=42164.365D0)
! NOMINAL RADIAL DISTANCE OF SATELLITE (km)
      REAL*8 AE
      PARAMETER (AE=6378.137D0)
! EARTH EQUATORIAL RADIUS (km)
      REAL*8 FER
      PARAMETER (FER=1.D0/298.25D0)
! EARTH FLATTENING COEFFICIENT = 1-(BE/AE)
      REAL*4 AEBE2
      PARAMETER (AEBE2=1.D0/(1.D0-FER)**2)
      REAL*4 AEBE3
      PARAMETER (AEBE3=AEBE2-1.)
      REAL*4 AEBE4
      PARAMETER (AEBE4=(1.D0-FER)**4-1.)
!
!***********************************************************************
!     INCLUDE 'ELCOMM.INC'
!***********************************************************************
!***********************************************************************
!**
!** INTEGRAL SYSTEMS, INC.
!**
!***********************************************************************
!**
!** PROJECT : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT
!** SYSTEM : EARTH LOCATION USERS GUIDE
!** NAME : ELCOMM
!** TYPE : DATA AREA
!** SOURCE : ELCOMM.INC
!**
!** VER. DATA BY COMMENT
!** ---- -------- ---------------- --------------------------------
!** A 01/09/89 I. LEVINE INITIAL CREATION
!**
!***********************************************************************
!**
!** DESCRIPTION
!** INSTRUMENT POSITION AND ATTITUDE VARIABLES AND TRANSFORMATION
!** MATRIX
!**
!***********************************************************************
!***********************************************************************
!
! COMMON VARIABLES
!
      REAL*8 XS(3)
! NORMALIZED S/C POSITION IN ECEF COORDINATES
      REAL*8 BT(3,3)
! ECEF TO INSTRUMENT COORDINATES TRANSFORMATION
      REAL*8 Q3
! USED IN SUBROUTINE LPOINT
      REAL*8 PITCH,ROLL,YAW
! PITCH,ROLL,YAW ANGLES OF INSTRUMENT (RAD)
      REAL*8 PMA,RMA
! PITCH,ROLL MISALIGNMENTS OF INSTRUMENT (RAD)
      COMMON /ELCOMM/ XS,BT,Q3,PITCH,ROLL,YAW,PMA,RMA
!***********************************************************************
!     INCLUDE 'INSTCOMM.INC'
!***********************************************************************
!***********************************************************************
!**
!** INTEGRAL SYSTEMS, INC.
!**
!***********************************************************************
!**
!** PROJECT : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT
!** SYSTEM : EARTH LOCATION USERS GUIDE
!** NAME : INSTCOMM
!** TYPE : DATA AREA
!** SOURCE : INSTCOMM.INC
!**
!** VER. DATA BY COMMENT
!** ---- -------- ---------------- --------------------------------
!** A 02/16/89 I. LEVINE INITIAL CREATION
!**
!***********************************************************************
!**
!** DESCRIPTION
!** COMMON AREA FOR INSTRUMENT-RELATED CONTROL PARAMETERS
!**
!***********************************************************************
!***********************************************************************
!
! VARIABLES
! CONSTANTS NEEDED TO PERFORM TRANSFORMATIONS BETWEEN THE
! LATITUDE/LONGITUDE, LINE/PIXEL AND INSTRUMENT CYCLES/INCREMENTS
! COORDINATES.
!
      INTEGER*4 INCMAX(2)
! NUMBER OF INCREMENTS PER CYCLE
      REAL*4 ELVMAX(2)
! BOUNDS IN ELEVATION (RADIANS)
      REAL*4 SCNMAX(2)
! BOUNDS IN SCAN ANGLE (RADIANS)
      REAL*4 ELVINCR(2)
! CHANGE IN ELEVATION ANGLE PER INCREMENT (RAD)
      REAL*4 SCNINCR(2)
! CHANGE IN SCAN ANGLE PER INCREMENT (RADIANS)
      REAL*4 ELVLN(2)
! ELEVATION ANGLE PER DETECTOR LINE (RADIANS)
      REAL*4 SCNPX(2)
! SCAN ANGLE PER PIXEL (RADIANS)
      REAL*4 EWNOM(2)
! EW CENTER OF INSTRUMENT
      REAL*4 NSNOM(2)
! NS CENTER OF INSTRUMENT
      COMMON /INSTCOMM/ INCMAX,ELVMAX,SCNMAX,    &
            ELVINCR,SCNINCR,ELVLN,SCNPX,EWNOM,NSNOM

!***********************************************************************
!***********************************************************************
      IERR=1
!
! COMPUTE SIGN OF MISALIGNMENT CORRECTIONS AND ORIGIN OFFSET CORRECTIONS
!
      FF = FLIP_FLG
      IF (INSTR == 2) FF = - FF
            DOFF = SCNMAX(INSTR) - EWNOM(INSTR)
!
! ADD THE NEW SECOND ORDER ORIGIN OFFSET CORRECTION
!
      ALPHA = ALPHA0 - ALPHA0*ZETA0*DOFF
      ZETA = ZETA0 + 0.5*ALPHA0*ALPHA0*DOFF

!
! COMPUTES TRIGONOMETRIC FUNCTIONS OF THE SCAN AND ELEVATION
! ANGLES CORRECTED FOR THE ROLL AND PITCH MISALIGNMENTS
!
      CA=COS(ALPHA)
      SA=SIN(ALPHA)
      CZ=COS(ZETA)
      DA=ALPHA-PMA*SA*(FF/CZ+TAN(ZETA))-RMA*(1.-CA/CZ)
      DZ=ZETA+FF*RMA*SA
!
! COMPUTES POINTING VECTOR IN INSTRUMENT COORDINATES
!
      CZ=COS(DZ)
      G(1)=SIN(DZ)
      G(2)=-CZ*SIN(DA)
      G(3)=CZ*COS(DA)
!
! TRANSFORMS THE POINTING VECTOR TO EARTH-FIXED COORDINATES
!
      G1(1)=BT(1,1)*G(1)+BT(1,2)*G(2)+BT(1,3)*G(3)
      G1(2)=BT(2,1)*G(1)+BT(2,2)*G(2)+BT(2,3)*G(3)
      G1(3)=BT(3,1)*G(1)+BT(3,2)*G(2)+BT(3,3)*G(3)
!
! COMPUTES COEFFICIENTS AND SOLVES A QUADRATIC EQUATION TO
! FIND THE INTERSECT OF THE POINTING VECTOR WITH THE EARTH
! SURFACE
!
      Q1=G1(1)**2+G1(2)**2+AEBE2*G1(3)**2
      Q2=XS(1)*G1(1)+XS(2)*G1(2)+AEBE2*XS(3)*G1(3)
      D=Q2*Q2-Q1*Q3
      IF (DABS(D).LT.1.D-9) D=0.
!
! IF THE DISCIMINANTE OF THE EQUATION, D, IS NEGATIVE, THE
! INSTRUMENT POINTS OFF THE EARTH
!
      IF (D.LT.0) THEN
            RLAT=999999.
            RLON=999999.
      RETURN
      END IF
      D=DSQRT(D)
!
! SLANT DISTANCE FROM THE SATELLITE TO THE EARTH POINT
!
      H=-(Q2+D)/Q1
!
! CARTESIAN COORDINATES OF THE EARTH POINT
!
      U(1)=XS(1)+H*G1(1)
      U(2)=XS(2)+H*G1(2)
      U(3)=XS(3)+H*G1(3)
!
! SINUS OF GEOCENTRIC LATITUDE
!
      D1=U(3)/SQRT(U(1)**2+U(2)**2+U(3)**2)
!
! GEOGRAPHIC (GEODETIC) COORDINATES OF THE POINT
!
      RLAT=ATAN(AEBE2*D1/SQRT(1.-D1*D1))
      RLON=ATAN2(U(2),U(1))
      IERR=0
      RETURN
      END SUBROUTINE LPOINT
!***********************************************************************        
!***********************************************************************        
!**                                                                             
!**   INTEGRAL SYSTEMS, INC.                                                    
!**                                                                             
!***********************************************************************        
!**                                                                             
!**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
!**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
!**   ROUTINE   : GPOINT                                                        
!**   SOURCE    : F.GPOINT                                                      
!**   LOAD NAME : ANY                                                           
!**   PROGRAMMER: IGOR LEVINE                                                   
!**                                                                             
!**   VER.    DATA    BY   COMMENT                                              
!**   ----  --------  ---  ---------------------------------------------        
!**   A     12/10/87  IL   INITIAL CREATION                                     
!**   A     06/10/88  IL   REPLACED ASIN WITH ATAN TO SAVE TIME                 
!**   A     06/02/89  IL   COORDINATE AXES CHANGED ACCORDING TO                 
!**                        FORD'S DEFINITION IN SDAIP, DRL 504-01               
!**                                                                             
!***********************************************************************        
!**                                                                             
!**   THIS subroutine CONVERTS GEOGRAPHIC LATITUDE AND LONGITUDE                
!**   TO THE RELATED ELEVATION AND SCAN ANGLES.                                 
!**                                                                             
!***********************************************************************        
!**                                                                             
!**   CALLED BY       : ANY                                                     
!**   COMMONS MODIFIED: NONE                                                    
!**   INPUTS          : NONE                                                    
!**   OUTPUTS         : NONE                                                    
!**   ROUTINES CALLED : NONE                                                    
!**                                                                             
!***********************************************************************        
!***********************************************************************        
      subroutine GPOINT(RLAT,RLON,ALF,GAM,IERR)                                 
!                                                                               
!     CALLING PARAMETERS                                                        
!                                                                               
      real*8   RLAT	! GEOGRAPHIC LATITUDE IN RADIANS (INPUT)            
      real*8   RLON	! GEOGRAPHIC LONGITUDE IN RADIANS (INPUT)           
      real*8   ALF	! ELEVATION ANGLE IN RADIANS (OUTPUT)               
      real*8   GAM	! SCAN ANGLE IN RADIANS (OUTPUT)                    
      integer IERR	! OUTPUT STATUS
!                             0 - SUCCESSFUL COMPLETION,         
!                             1 - POINT WITH GIVEN LAT/LON IS INVISIBLE         
!
!     LOCAL VARIABLES                                                           
!
      real*8 F(3)	! POINTING VECTOR IN EARTH CENTERED COORDINATES
      real*8 FT(3)	! POINTING VECTOR IN INSTRUMENT COORDINATES
      real*8 U(3)	! COORDINATES OF THE EARTH POINT (KM)
      real*8 SING,SLAT,W1,W2  ! WORK SPACE
!                                                                               
!     INCLUDE FILES                                                             
!                                                                               
      real*8 PI                                                                 
           PARAMETER (PI=3.141592653589793D0)                                   
      real*8 DEG                                                                
           PARAMETER (DEG=180.D0/PI)                                            
      real*8 RAD                                                                
           PARAMETER (RAD=PI/180.D0)                                            
!                    DEGREES TO RADIANS CONVERSION PI/180                       
      real*8 NOMORB                                                             
           PARAMETER (NOMORB=42164.365D0)                                       
!                    NOMINAL RADIAL DISTANCE OF SATELLITE (km)                  
      real*8 AE                                                                 
           PARAMETER (AE=6378.137D0)                                            
!                    EARTH EQUATORIAL RADIUS (km)                               
      real*8 FER                                                                
           PARAMETER (FER=1.D0-(6356.7533D0/AE))                                
!                    EARTH FLATTENING COEFFICIENT = 1-(BE/AE)                   
      real*4 AEBE2                                                              
           PARAMETER (AEBE2=1.D0/(1.D0-FER)**2)                                 
      real*4 AEBE3                                                              
           PARAMETER (AEBE3=AEBE2-1.)                                           
      real*4 AEBE4                                                              
           PARAMETER (AEBE4=(1.D0-FER)**4-1.)
      real*8 XS(3)                                                              
!                      NORMALIZED S/C POSITION IN ECEF COORDINATES              
      real*8 BT(3,3)                                                            
!                      ECEF TO INSTRUMENT COORDINATES TRANSFORMATION            
      real*8  Q3                                                                
!                      USED IN subroutine LPOINT                                
      real*8 PITCH,ROLL,YAW                                                     
!                          PITCH,ROLL,YAW ANGLES OF INSTRUMENT (RAD)            
      real*8 PMA,RMA                                                            
!                          PITCH,ROLL MISALIGNMENTS OF INSTRUMENT (RAD)         
         COMMON /ELCOMM/ XS,BT,Q3,PITCH,ROLL,YAW,PMA,RMA
!***********************************************************************        
!                                                                               
!     COMPUTES SINUS OF GEOGRAPHIC (GEODETIC) LATITUDE                          
!                                                                               
      SING=DSIN(RLAT)
      W1=AEBE4*SING*SING                                                        
!                                                                               
!     SINUS OF THE GEOCENTRIC LATITUDE                                          
!                                                                               
      SLAT=((0.375D0*W1-0.5D0)*W1+1.0D0)*SING/AEBE2
!                                                                               
!     COMPUTES LOCAL EARTH RADIUS AT SPECIFIED POINT                            
!                                                                               
      W2=SLAT*SLAT                                                              
      W1=AEBE3*W2                                                               
      W1=(0.375D0*W1-0.5D0)*W1+1.D0
!                                                                               
!     COMPUTES CARTESIAN COORDINATES OF THE POINT                               
!                                                                               
      U(3)=SLAT*W1                                                              
      W2=W1*DSQRT(1.0D0-W2)
      U(1)=W2*DCOS(RLON)
      U(2)=W2*DSIN(RLON)
!                                                                               
!     POINTING VECTOR FROM SATELLITE TO THE EARTH POINT                         
!                                                                               
      F(1)=U(1)-XS(1)                                                           
      F(2)=U(2)-XS(2)                                                           
      F(3)=U(3)-XS(3)                                                           
      W2=U(1)*SNGL(F(1))+U(2)*SNGL(F(2))+           & 
         U(3)*SNGL(F(3))*AEBE2                                                  
!                                                                               
!     VERIFIES VISIBILITY OF THE POINT                                          
!                                                                               
      if (W2.GT.0.0D0) then
!                               INVISIBLE POINT ON THE EARTH                    
                   IERR=1                                                       
                   ALF=99999.0D0
                   GAM=99999.0D0
                   return                                                       
       end if                                                                   
!                                                                               
!     CONVERTS POINTING VECTOR TO INSTRUMENT COORDINATES                        
!                                                                               
      FT(1)=BT(1,1)*F(1)+BT(2,1)*F(2)+BT(3,1)*F(3)                              
      FT(2)=BT(1,2)*F(1)+BT(2,2)*F(2)+BT(3,2)*F(3)                              
      FT(3)=BT(1,3)*F(1)+BT(2,3)*F(2)+BT(3,3)*F(3)                              
!                                                                               
!     CONVERTS POINTING VECTOR TO SCAN AND ELEVATION ANGLES AND                 
!     CORRECTS FOR THE ROLL AND PITCH MISALIGNMENTS                             
!                                                                               
      GAM=ATAN(FT(1)/SQRT(FT(2)**2+FT(3)**2))                                   
      ALF=-DATAN(FT(2)/FT(3))
      W1=DSIN(ALF)
      W2=DCOS(ALF)
      ALF=ALF+RMA*(1.0D0-W2)+PMA*W1*(1.0D0+DTAN(GAM))
      GAM=GAM-RMA*W1                                                            
      IERR=0                                                                    
      return                                                                    
      end  subroutine GPOINT                                                                      
!=======================================================================
!     C O M P _ E S
!=======================================================================
      subroutine COMP_ES( NAVstr, line, pixel, elev, scan )
!
!     AUTHOR:         KELLY DEAN
!
!     CREATED:        January 1995
!   
!     DEVELOPED FOR:  CIRA/COLORADO STATE UNIVERSITY
!
!     PURPOSE:        
!		Compute the elevation and scan angles related to the
!	satellite line and pixel numbers.
!
!     REVISION:       1.0
!
!     ARGUMENTS:
!    NAME:    type:       PURPOSE:                        IN/OUT:
!	 NAVstr   Structure   Navigation information         IN
! 	 line     real*8      Satellite line number          IN
! 	 pixel    real*8      Satellite pixel number         IN
! 	 elev     real*8      Elevation angle (rad)          OUT
! 	 scan     real*8      Scan angle (rad)               OUT
!    
!     CONSTANTS:
!       NAME:    PURPOSE:
!       ****************  real*8     *****************
!	elvln   Elevation angle per detector line (rad)
!	elvmax  Bounds in elevation
!	scnmax  Bounds in scan angle
!	scnpx   Scan angle per pixel (rad)
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     INCLUDE DECLARATION SECTION:
!
!     include 'kdf.inc'
!
!     CONSTANT DECLARATION SECTION:
!
      integer incmax(2) /6136,2805/
      real*8  elvmax(2) / 0.2208960D0, 0.22089375D0/
      real*8  elvln(2)  /28.0D-6,    280.0D-6/
      real*8  elvinc(2) / 8.0D-6,     17.5D-6/
      real*8  scnmax(2) / 0.245440D0,  0.2454375D0 /
      real*8  scnpx(2)  /16.0D-6,    280.0D-6/
      real*8  scninc(2) /16.0D-6,     35.0D-6/
!
!     VARIABLE DECLARATION SECTION:
!
      type(GVAR_NAV):: NAVstr
      real*8 line, pixel, elev, scan
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
      if ( NAVstr%instr == 1 ) then 
!	Recompute elevation and scan biases based on user inputs of
!	cycles and increments obtained from GVAR.
       elvmax(NAVstr%instr) = ( NAVstr%ns_cyl *      &
                               incmax(NAVstr%instr) +      &
                               NAVstr%ns_inc ) *      &
                               elvinc(NAVSTR%instr)
       scnmax(NAVstr%instr) = ( NAVstr%ew_cyl *     &
                                incmax(NAVstr%instr) +     & 
                                NAVstr%ew_inc ) *      &
                                scninc(NAVstr%instr)
!      Compute elevation angle (rad)
       elev = elvmax(NAVstr%instr) + (4.50 - line) * elvln(NAVstr%instr) 
!      Compute scan angle (rad)
       scan = (pixel - 1.0) * scnpx(NAVstr%instr) - scnmax(NAVstr%instr) 
      else IF( NAVstr%instr == 2 ) then
!	Recompute elevation and scan biases based on user inputs of
!	cycles and increments obtained from GVAR.
       elvmax(NAVstr%instr)  = ( (9 - NAVstr%ns_cyl) *    &
                                incmax(NAVstr%instr) -     &
                                NAVstr%ns_inc ) *     &
                                elvinc(NAVstr%instr)
       scnmax(NAVstr%instr) = ( NAVstr%ew_cyl *     &
                                incmax(NAVstr%instr) + & 
                                NAVstr%ew_inc ) *  &
                                scninc(NAVstr%instr)
!      Compute elevation angle (rad)
       elev = elvmax(NAVstr%instr) + (2.50 - line) * elvln(NAVstr%instr) 
!      Compute scan angle (rad)
       scan = (pixel - 1.0)*scnpx(NAVstr%instr)-scnmax(NAVstr%instr)
      else
!      Unknown instrument.....
       elev = 0.0D0
       scan = 0.0D0
      endif
!
!
      return
      end subroutine COMP_ES
!======================================================================
!     C O M P _ L P
!======================================================================
      subroutine COMP_LP( NAVstr, ELEV, SCAN, RL, RP )
!
!     AUTHOR:         KELLY DEAN
!
!     CREATED:        January 1995
!   
!     DEVELOPED FOR:  CIRA/COLORAdo STATE UNIVERSITY
!
!     PURPOSE:        
!	  Subroutine COMP_LP converts elevation and scan angles to the
!       fractional line and pixel numbers.
!
!     REVISION:       0.0
!
!     ARGUMENTS:
!       NAME:   type:       PURPOSE:                        IN/OUT:
!	 NAVstr   Structure   Navigation information         IN
!	 ELEV	  real*8     Elevation angle (rad)           IN
!	 SCAN	  real*8:    Scan angle (rad)                IN
!	 RL	  real*8     Line Number                     OUT
!	 RP	  real*8     Pixel Number                    OUT
!
!     CONSTANTS:
!       NAME:    PURPOSE:
!       ****************  real*8       *****************
!  	elvln    Elevation angle per detector line (rad)
! 	elvmax   Bounds in elevation
!	scnmax   Bounds in scan angle
!	scnpx    Scan angle per pixel (rad)
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     INCLUDE DECLARATION SECTION:
!
!     include 'kdf.inc'
!
!     CONSTANT DECLARATION SECTION:
!
      integer incmax(2) /6136,2805/
      real*8  elvmax(2) / 0.2208960D0, 0.22089375D0/
      real*8  elvln(2)  /28.0D-6,    280.0D-6/
      real*8  elvinc(2) / 8.0D-6,     17.5D-6/
      real*8  scnmax(2) / 0.245440D0,  0.2454375D0 /
      real*8  scnpx(2)  /16.0D-6,    280.0D-6/
      real*8  scninc(2) /16.0D-6,     35.0D-6/
!
!     VARIABLE DECLARATION SECTION:
!
      type(GVAR_NAV):: NAVstr
      real*8 ELEV, SCAN, RL, RP
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
      if ( NAVstr%instr == 1 ) then
!	Recompute elevation and scan biases based on user inputs of
!	cycles and increments obtained from GVAR.
       elvmax(NAVstr%instr) = ( NAVstr%ns_cyl *      &
                               incmax(NAVstr%instr) +      &
                               NAVstr%ns_inc ) *      &
                               elvinc(NAVSTR%instr)
       scnmax(NAVstr%instr) = ( NAVstr%ew_cyl *     &
                                incmax(NAVstr%instr) +      &
                                NAVstr%ew_inc ) *           &
                                scninc(NAVstr%instr)
!       Compute fractional line number.
      RL = ( ELVMAX(NAVstr%instr) - ELEV ) / ELVLN(NAVstr%instr) 
      RL = RL + 4.5D0
!       Compute fractional pixel number.
       RP = ( SCNMAX(NAVstr%instr) + SCAN ) / SCNPX(NAVstr%instr)+1.0D0 
      else if ( NAVstr%instr == 2 ) then
!       Recompute elevation and scan biases based on user inputs of
!       cycles and increments obtained from GVAR.
        elvmax(NAVstr%instr)  = ( (9 - NAVstr%ns_cyl) *    &
                                incmax(NAVstr%instr) -     &
                                NAVstr%ns_inc ) *     &
                                elvinc(NAVstr%instr)   
        scnmax(NAVstr%instr) = ( NAVstr%ew_cyl *     &
                                incmax(NAVstr%instr) +    & 
                                NAVstr%ew_inc ) *     &
                                scninc(NAVstr%instr)
!       Compute fractional line number.
        RL = ( ELVMAX(NAVstr%instr) - ELEV ) / ELVLN(NAVstr%instr) 
        RL = RL + 2.5D0
!       Compute fractional pixel number.
        RP = ( SCNMAX(NAVstr%instr)+SCAN ) / SCNPX(NAVstr%instr)+1.0D0 
      else
!      Unknown instrument.....
       RL = 0.0D0
       RP = 0.0D0
      endif
!
      return
      end subroutine COMP_LP
!======================================================================
! End of CIRA GOES IMAGER ROUTINES
!======================================================================
!--------------------------------------------------------------
! Subroutine to make geostationary satellite azimuth field
!   
!     glon = longitude of sub-satellite point (positive for western hem)
!     glat = latitude of sub-satellite point (positive for northern hem)
!
!     zenith  = satellite zenith view angle 
!     azimuth = satellite azimuth angle clockwise from north
!--------------------------------------------------------------
 subroutine COMPUTE_SATELLITE_ANGLES(glon, glat, Line_Idx)

!  arguments
   real (kind=real8),  intent(in)  :: glon, glat
   integer(kind=int4), intent(in):: Line_Idx

!  Local variables
   integer :: Elem_Idx, n
   real:: Satlon
   real:: Satlat

   n = size(lon(:,Line_Idx))

   Satlon = glon
   Satlat = glat

   do Elem_Idx = 1, n


    if (Space_Mask(Elem_Idx,Line_Idx) == sym%NO) then

!     Satzen(Elem_Idx,Line_Idx) = COMPUTE_SENSOR_ZENITH_GEO(satlon,satlat, &
!                                Lon_1b(Elem_Idx,Line_Idx),Lat_1b(Elem_Idx,Line_Idx))

      Satzen(Elem_Idx,Line_Idx) = SENSOR_ZENITH(GEO_ALTITUDE,satlon,satlat, &
                                 Lon_1b(Elem_Idx,Line_Idx),Lat_1b(Elem_Idx,Line_Idx))
      
      Sataz(Elem_Idx,Line_Idx) = SENSOR_AZIMUTH(Satlon,Satlat, &
                                                Lon_1b(Elem_Idx,Line_Idx),Lat_1b(Elem_Idx,Line_Idx))

      Relaz(Elem_Idx,Line_Idx) = RELATIVE_AZIMUTH(Solaz(Elem_Idx,Line_Idx), Sataz(Elem_Idx,Line_Idx))

      Glintzen(Elem_Idx,Line_Idx) = GLINT_ANGLE(Solzen(Elem_Idx,Line_Idx), &
                                      Satzen(Elem_Idx,Line_Idx),Relaz(Elem_Idx,Line_Idx)) 

      Scatangle(Elem_Idx,Line_Idx) = SCATTERING_ANGLE(Solzen(Elem_Idx,Line_Idx), &
                                      Satzen(Elem_Idx,Line_Idx),Relaz(Elem_Idx,Line_Idx)) 

    else

     Solzen(Elem_Idx,Line_Idx) = Missing_Value_Real4
     Sataz(Elem_Idx,Line_Idx) = Missing_Value_Real4
     Satzen(Elem_Idx,Line_Idx) = Missing_Value_Real4
     Relaz(Elem_Idx,Line_Idx) = Missing_Value_Real4
     Glintzen(Elem_Idx,Line_Idx) = Missing_Value_Real4
     Scatangle(Elem_Idx,Line_Idx) = Missing_Value_Real4

    endif

   enddo

 end subroutine COMPUTE_SATELLITE_ANGLES
!-----------------------------------------------------------------------------
! Logic copied from MCIDAS gvar_pfx.for, version unknown (last edit 23Dec2006) 
! Print the data block prefix (time)
!
! buf = an int2 vector of the areafile words on a line including the prefix
!-----------------------------------------------------------------------------
subroutine PRINT_PREFIX(buf, ms_Time)

  integer (kind=int2), DIMENSION(:), INTENT(in) :: buf
  integer (kind=int4), INTENT(out) :: ms_Time
  real (kind=real4):: frac_Hours
  integer ITIMES(16)
  integer year,dayofyr,Hour,min,sec,msec
  real (kind=real4) :: Minute_Time
  integer (kind=int2), dimension(128) :: buf2
  logical*1 :: ldoc(256)

!  integer, parameter :: LOC = 4  ! imager value, from mcidas code
   integer, parameter :: LOC = -1   ! value experimentally determined to work for CLAVRx. Shifted over 2 words, hence -1
!  integer, parameter :: LOC = 0  ! sounder value, from mcidas code
 
  equivalence(buf2, ldoc)
  
  buf2(:) = buf(:128)

  CALL UNPKTIME (ldoc ,ITIMES,9+LOC)
  
  year = itimes(1)*1000 + itimes(2)*100 + itimes(3)*10 + itimes(4)
  dayofyr = itimes(5)*100 + itimes(6)*10 + itimes(7)
  Hour = itimes(8)*10 + itimes(9)
  min  = itimes(10)*10 + itimes(11)
  sec  = itimes(12)*10 + itimes(13)
  msec = itimes(14)*100 + itimes(15)*10 + itimes(16)
  
! PRINT STATEMENT KEPT for verification
!      write (6,6005) year,dayofyr,Hour,min,sec,msec
!6005  FORMAT ('time: year',i5,' day of year',i4,' hh:mm:ss ',i2,1h:,i2,1h:,i2,' msec',i4)

  !Calculate miliseconds since midnight of current day
   ms_Time = (((Hour * 60 * 60) + (min * 60) + sec) * 1000) + msec
 
  !Calculates fractions of an Hour since midnight of current day
  Minute_Time = min + (((msec / 1000.0) + sec) / 60.0 )
  frac_Hours = Hour + (Minute_Time / 60.0)


end subroutine PRINT_PREFIX


! FIXME: provisional routine -- from Mcidas, extracts line date and time from data buffer
! Minimally modified to compile with gfortran options 
subroutine UNPKTIME (LDOC,ITIMES,LOC)
!C LOC is 13 for sounder,  11 for imager
!CCC      integer IDOC(256)                    ! This array is NOT USED -- 10/24/2003  JPN
  integer ITIMES(16)

  integer mask1, mask2
  integer mask3
  logical lword
!  logical*1 ldoc(256)
  logical*1, intent(in) :: ldoc(256)

  integer LOC       
  integer I,K,iword 

  equivalence (lword,iword)
  data mask1 /z'0000000F'/, mask2 /z'000000F0'/
  data mask3 /z'000000FF'/
  !C
  K = 0
  !      do 1 I=1,8
  do I=1, 8
     lword = ldoc(LOC+i)
     iword = iand(iword,mask3)
     K = K + 1
     ITIMES(k) = iand(iword,mask2)
     ITIMES(k) = ISHFT(ITIMES(k),-4)
     K = K + 1
     ITIMES(k) = iand(iword,mask1)
     !1        CONTINUE
  end do
  return
end subroutine UNPKTIME

!----------------------------------------------------------------------
! Determine the file name of this dark composite name
!----------------------------------------------------------------------
subroutine DETERMINE_DARK_COMPOSITE_NAME(AREAstr)

 type(AREA_STRUCT), intent(in):: AREAstr
 character(len=9):: Goes_Name
 character(len=4):: Year_String
 character(len=3):: Jday_String
 character(len=4):: Time_String
 logical:: Does_File_Exist
 integer:: Hour
 integer:: Minute
 integer:: itime
 integer(kind=int4), parameter :: MAX_LATENCY = 5 !including current day
 integer:: year
 integer:: year_in
 integer:: day_of_year_in
 integer:: doy
 integer:: doy_idx

 Dark_Composite_Name = "no_file"

 !--- only do this for geostationary imager data
 if (Goes_Flag /= sym%NO .AND.  &
     COMS_Flag == sym%NO .AND. &
     Mtsat_Flag == sym%NO .AND. &
     Seviri_Flag /= sym%NO ) then
     return
 endif

 select case (Sc_Id_WMO)
    case(252)
         Goes_Name = "goes08"
    case(253)
         Goes_Name = "goes09"
    case(254)
         Goes_Name = "goes10"
    case(255)
         Goes_Name = "goes11"
    case(256)
         Goes_Name = "goes12"
    case(257)
         Goes_Name = "goes13"
    case(258)
         Goes_Name = "goes14"
    case(259)
         Goes_Name = "goes15"
    case(55)
         Goes_Name = "met8"
    case(56) 
         Goes_Name = "met9"
    case(57)
         Goes_Name = "met10"
    case(171) 
         Goes_Name = "mtsat-1r"
    case(172) 
         Goes_Name = "mtsat-2"
    case(810) 
         Goes_Name = "coms-1"
    case default
         return
 end select

 Hour = AREAstr%Img_Time/10000
 Minute = (AREAstr%Img_Time - Hour*10000)/100
 Itime = Hour*100 + Minute
 write (Time_String,fmt="(I4.4)") Itime

 !----------------------------------------------------------------------
 ! loop through today and previous days to find the most current file
 ! if no file can be found, the file name is set to "no_file"
 !----------------------------------------------------------------------
 year_in = Start_Year
 day_of_year_in = Start_Day

 Dark_Composite_Name = "no_file"

 do doy_idx = 0, MAX_LATENCY - 1
    doy = day_of_year_in - doy_idx
    year = year_in
    ileap = leap_year_fct(year)
    if (doy < 1) then
           year = year - 1
           ileap = leap_year_fct(year)
           doy = (365 + ileap) + doy
    endif

    write (Year_String,fmt="(I4.4)") year
    write (Jday_String,fmt="(I3.3)") doy
 
    Dark_Composite_Name = trim(Goes_Name)//"_"//Year_String//"_"// &
                           Jday_String//"_"//Time_String// &
                           "_drk_ch1_pix.dat"

    Dark_Comp_Data_Dir_Temp = trim(Dark_Comp_Data_Dir) // &
                              trim(Goes_Name)//"/"//Year_String//"/"


    !--- test for existence - assume uncompressed
    Dark_Composite_Name = trim(Goes_Name)//"_"//Year_String//"_"// &
                           Jday_String//"_"//Time_String// &
                           "_drk_ch1_pix.dat"

    Does_File_Exist = file_exists(trim(Dark_Comp_Data_Dir_Temp)// &
                                  trim(Dark_Composite_Name))

    !-- if found, exit loop
    if (Does_File_Exist .neqv. .false.) then
        exit
    endif

    !--- test for existence - assume gzip compressed
    Dark_Composite_Name = trim(Dark_Composite_Name)//".gz"

    Does_File_Exist = file_exists(trim(Dark_Comp_Data_Dir_Temp)// &
                                  trim(Dark_Composite_Name))

    !-- if found, exit loop
    if (Does_File_Exist .neqv. .false.) then
        exit
    endif

 end do

end subroutine DETERMINE_DARK_COMPOSITE_NAME

!-----------------------------------------------------------------------------------------------------
! public routine to read data from an Ch1 Dark Composite area file for one segment into memory
!-----------------------------------------------------------------------------------------------------
subroutine READ_DARK_COMPOSITE_COUNTS(Segment_Number,Xstride,Dark_Composite_Filename, &
                                      AREAstr_Image,Ch1_Dark_Composite_Counts)

   integer(kind=int4), intent(in):: Segment_Number
   integer(kind=int4), intent(in):: Xstride
   character(len=*), intent(in):: Dark_Composite_Filename
   type (AREA_STRUCT), intent(in) :: AREAstr_Image
   integer(kind=int2), dimension(:,:), intent(out):: Ch1_Dark_Composite_Counts
   integer:: Io_Status

   type (AREA_STRUCT), save :: AREAstr_Dark
   integer:: Dark_Lun_Header
   integer, save:: Dark_Lun_Data
   integer, save:: Element_Offset
   integer, save:: Line_Offset
   integer:: Store_Time_Flag
   integer:: Rec_Num
   integer:: First_Rec_Num
   integer:: Last_Rec_Num
   integer:: Num_Recs_To_Read
   integer:: Num_Elements
   integer, save:: Num_Elements_Dark
   integer:: First_Line_In_Segment
   integer:: Last_Line_In_Segment
   integer:: Num_Lines_Per_Segment
   integer:: Num_Lines_Read
   integer:: Big_Count
   integer:: Line_Idx

   integer(kind=int2), dimension(:), allocatable, save:: Dark_Comp_Counts_Temp
   integer(kind=int2), dimension(:), allocatable, save:: Dark_Comp_Counts

   integer:: Elem_Idx
   integer:: Temp_Idx 
   character(len=200):: System_String
   character(len=200):: Dark_Name_Full
   character(len=3):: Dark_Ext

   !--- aliases
   Num_Elements = Num_Pix
   Num_Lines_Per_Segment = Num_Scans_Per_Segment
   Num_Lines_Read = Num_Scans_Read

   !--- initialize
   Big_Count = 999
   Io_Status = 0
   Ch1_Dark_Composite_Counts = 0

   !--- check to see if a dark composite file exists for this image
   if (trim(Dark_Composite_Filename) == "no_file" .and. Segment_Number == 1) then
          print *, "No Dark Composite Available for this Image"
          return
   endif

   !--- do not store time from the dark composite
   Store_Time_Flag = sym%NO

   !--- on the first segment, open the file and read the Area and Nav headers
   !--- and compute the offsets between the image and the dark composite
   if (Segment_Number == 1) then 


        Dark_Name_Full = trim(Dark_Comp_Data_Dir_Temp)//trim(Dark_Composite_Name)

        !--- check for compressed files
        Dark_Ext = Dark_Composite_Name(len_trim(Dark_Composite_Name)-2:len_trim(Dark_Composite_Name))

        if (Dark_Ext == ".gz") then
          System_String = "gunzip -c "//trim(Dark_Comp_Data_Dir_Temp)// &
                        trim(Dark_Composite_Name)//".gz"// &
                        " > "//trim(Temporary_Data_Dir)//trim(Dark_Composite_Name)
          call system(System_String)
 
          Number_of_Temporary_Files = Number_of_Temporary_Files + 1
          Temporary_File_Name(Number_of_Temporary_Files) = trim(Dark_Composite_Name)

          Dark_Name_Full = trim(Temporary_Data_Dir)//trim(Dark_Composite_Name)

        endif

        !--- open for header read
        Dark_Lun_Header = Get_Lun()
        open(unit=Dark_Lun_Header,file=trim(Dark_Name_Full), &
             form="unformatted", access="direct",recl=4, &
             status="old",action="read", iostat=Io_Status)

        read (unit=Dark_Lun_Header, rec = 1, iostat=Io_Status) AREAstr_Dark%area_Status
        if (Io_Status /= 0) then
           Dark_Composite_Name = "no_file"
           print *, "Dark Composite Area file cannot be read"  ! byte swapping may cause this
           return
        endif
        read (unit=Dark_Lun_Header, rec = 2) AREAstr_Dark%Version_Num
        if (AREAstr_Dark%Version_Num .ne. 4) then
           print *, "Dark Composite Area file cannot be read due to version number filter"  ! byte swapping may cause this
           Dark_Composite_Name = "no_file"
           return
        endif
        read (unit=Dark_Lun_Header, rec = 3) AREAstr_Dark%Sat_Id_Num
        read (unit=Dark_Lun_Header, rec = 4) AREAstr_Dark%Img_Date
        read (unit=Dark_Lun_Header, rec = 5) AREAstr_Dark%Img_Time
        read (unit=Dark_Lun_Header, rec = 6) AREAstr_Dark%North_Bound        !beg_Sc
        read (unit=Dark_Lun_Header, rec = 7) AREAstr_Dark%West_Vis_Pixel     !beg_el
        read (unit=Dark_Lun_Header, rec = 8) AREAstr_Dark%z_coor
        read (unit=Dark_Lun_Header, rec = 9) AREAstr_Dark%Num_Line
        read (unit=Dark_Lun_Header, rec = 10) AREAstr_Dark%Num_Elem
        read (unit=Dark_Lun_Header, rec = 11) AREAstr_Dark%Bytes_Per_Pixel
        read (unit=Dark_Lun_Header, rec = 12) AREAstr_Dark%Line_Res
        read (unit=Dark_Lun_Header, rec = 13) AREAstr_Dark%Elem_Res

        close(unit=Dark_Lun_Header)

        !--------------------------------------------------------------
        ! check to see if dark image is consistent with image data
        !--------------------------------------------------------------
        if ((AREAstr_Image%Num_Elem /= AREAstr_Dark%Num_Elem) .or.    &
            (AREAstr_Image%Num_Line /= AREAstr_Dark%Num_Line)) then
            print *, EXE_PROMPT, "Dark Composite File Is Incompatible, Ignoring It"
            Dark_Composite_Name = "no_file"
            return
        endif

        Num_Elements_Dark = AREAstr_Dark%Num_Elem
        allocate(Dark_Comp_Counts_Temp(Num_Elements_Dark))
        allocate(Dark_Comp_Counts(Num_Elements_Dark))

        !--- compute offsets

        !--- east-west - note this is an absolute value
        Element_Offset = 0
        if (AREAstr_Dark%West_Vis_Pixel /= AREAstr_Image%West_Vis_Pixel) then
                Element_Offset = abs(AREAstr_Dark%West_Vis_Pixel  -    &
                                     AREAstr_Image%West_Vis_Pixel) / &
                                     AREAstr_Image%Elem_Res
        endif

        !--- north-south
        Line_Offset = 0
        if (AREAstr_Dark%North_Bound /= AREAstr_Image%North_Bound) then
                Line_Offset =  (AREAstr_Dark%North_Bound  -    &
                                AREAstr_Image%North_Bound) / &
                                AREAstr_Image%Line_Res
        endif

        !print *, "Shifts for this Dark Comp = ", Element_Offset, Line_Offset

        !-- open file for data read - use int2
        Dark_Lun_Data=GET_LUN()
        open(unit=Dark_Lun_Header,file=trim(Dark_Name_Full), &
             form="unformatted", access="direct",recl=2*Num_Elements_Dark, &
             status="old",action="read",iostat=Io_Status)

        if (Io_Status /= 0) then
            print *, EXE_PROMPT, "Dark Composite File Could Not be Reopened"
            Dark_Composite_Name = "no_file"
            return
        endif
       
   endif 

   if (Dark_Composite_Name /= "no_file") then

    !--- for this segment, compute the words to read in (accounting for offset)
    First_Line_In_Segment = (Segment_Number-1)*Num_Lines_Per_Segment + 1

    First_Line_In_Segment = First_Line_In_Segment + Line_Offset    !check this!!!!
 
    Last_Line_In_Segment = min(AREAstr_Image%Num_Line,Segment_Number*Num_Lines_Per_Segment)


    !--- determine records to read - this accounts for north-south shift
    First_Rec_Num = 1  + First_Line_In_Segment + Line_Offset
    Last_Rec_Num = First_Rec_Num + Num_Scans_Read - 1
    First_Rec_Num = max(2,First_Rec_Num)
    Last_Rec_Num = min(Last_Rec_Num,Num_Scans)
    Num_Recs_To_Read = Last_Rec_Num - First_Rec_Num + 1

    Line_Idx = 1
    do Rec_Num = First_Rec_Num, Last_Rec_Num

     !--- read data
     read(unit=Dark_Lun_Data,rec=Rec_Num,iostat=Io_Status) Dark_Comp_Counts_Temp
     if (Io_Status /= 0) then
          if (Segment_Number == 1) print *, EXE_PROMPT, "Dark Composite File Could Not be Read"
          Dark_Composite_Name = "no_file"
          return
     endif

     !-- apply east-west shift
     Dark_Comp_Counts = Dark_Comp_Counts_Temp

     if (AREAstr_Dark%West_Vis_Pixel < AREAstr_Image%West_Vis_Pixel) then
       Dark_Comp_Counts(1:Num_Elements_Dark - Element_Offset) = Dark_Comp_Counts_Temp(Element_Offset+1:Num_Elements_Dark) 
       Dark_Comp_Counts(Num_Elements_Dark-Element_Offset+1:Num_Elements_Dark) = Big_Count
     endif

     if (AREAstr_Dark%West_Vis_Pixel > AREAstr_Image%West_Vis_Pixel) then
       Dark_Comp_Counts(1:Element_Offset) = Big_Count
       Dark_Comp_Counts(Element_Offset+1:Num_Elements_Dark) = Dark_Comp_Counts_Temp(1:Num_Elements_Dark-Element_Offset) 
       Dark_Comp_Counts(Num_Elements_Dark-Element_Offset+1:Num_Elements_Dark) = Big_Count
     endif

     !-- apply x offset
     do Elem_Idx = 1, Num_Elements
        Temp_Idx = min(Num_Elements_Dark,1 + Xstride*(Elem_Idx-1))
        Ch1_Dark_Composite_Counts(Elem_Idx,Line_Idx) = Dark_Comp_Counts(Temp_Idx)
     enddo

     !-- increment Line_Idx
     Line_Idx = Line_Idx + 1

    enddo

   endif
  
   !--- close file
   if (Segment_Number == Num_Segments) then
      close(unit=Dark_Lun_Data)
      if (allocated(Dark_Comp_Counts_Temp)) deallocate(Dark_Comp_Counts_Temp)
      if (allocated(Dark_Comp_Counts)) deallocate(Dark_Comp_Counts)

      !--- delete uncompressed file
      System_String = "rm "//trim(Temporary_Data_Dir)//trim(Dark_Composite_Name)
      call system(System_String)

   endif

end subroutine READ_DARK_COMPOSITE_COUNTS
!----------------------------------------------------------------------------------------------
! GOES-specific calibration of the dark-sky ch1 composite counts
!----------------------------------------------------------------------------------------------
subroutine CALIBRATE_GOES_DARK_COMPOSITE(Ch1_Counts_Composite,Time_Since_Launch,Ref_Ch1_Dark)

integer(kind=int2), dimension(:,:), intent(in):: Ch1_Counts_Composite
real(kind=real4), dimension(:,:),  intent(out):: Ref_Ch1_Dark
real(kind=real4), intent(in):: Time_Since_Launch

   call CALIBRATE_GOES_REFLECTANCE(Ch1_Counts_Composite,Ch1_Dark_Count,Time_Since_Launch,Ref_Ch1_Dark)

end subroutine CALIBRATE_GOES_DARK_COMPOSITE
!----------------------------------------------------------------------------------------------
! Apply PATMOS-x reflectance coefficients to GOES counts
!----------------------------------------------------------------------------------------------
subroutine CALIBRATE_GOES_REFLECTANCE(Counts,Dark_Count,Time_Since_Launch,Reflectance)

integer(kind=int2), dimension(:,:), intent(in):: Counts
real(kind=real4), intent(in):: Dark_Count
real(kind=real4), intent(in):: Time_Since_Launch
real(kind=real4), dimension(:,:), intent(out):: Reflectance
real:: Ch1_Gain

   !--- determine the gain
   Ch1_Gain = Ch1_Gain_Low_0*(100.0+Ch1_Degrad_Low_1*Time_Since_Launch + &
              Ch1_Degrad_Low_2*Time_Since_Launch**2)/100.0

   !--- calibrate the counts (this is an un-normalized reflectance)
   Reflectance = Missing_Value_Real4
   where(Counts /= Missing_Value_Int2) 
       Reflectance = Ch1_Gain*(Counts - Dark_Count)
   end where

end subroutine CALIBRATE_GOES_REFLECTANCE
!----------------------------------------------------------------------------------------------
! post process the dark-sky ch1 composite counts
! goal is to improve accuracy in perennially cloudy regions
!
! Ref_Ch1_Dark - The 0.65 micron reflectance from the dark composite
! Land_Class_Local - local copy of global Land array.  (7 = Deep Ocean)
! N_Box = pixel radius of area to search for minimum
!----------------------------------------------------------------------------------------------
subroutine POST_PROCESS_GOES_DARK_COMPOSITE(Ref_Ch1_Dark, Land_Class_Local)

   real(kind=real4), dimension(:,:),  intent(inout):: Ref_Ch1_Dark
   integer(kind=int1), dimension(:,:),  intent(in):: Land_Class_Local
   integer, parameter:: N_Box = 1
!  integer, parameter:: N_Box = 5
   integer, parameter:: Land_Class_Allowed = 7
   integer:: Elem_Idx
   integer:: Line_Idx
   integer:: Num_Elements
   integer:: Num_Lines
   integer:: Elem_Idx_1
   integer:: Elem_Idx_2
   integer:: Line_Idx_1
   integer:: Line_Idx_2
   real(kind=real4):: Z_Median
   real(kind=real4):: Z_Mean
   real(kind=real4):: Z_Std_Median
  
   Num_Elements = size(Ref_Ch1_Dark(:,1)) 
   Num_Lines = size(Ref_Ch1_Dark(1,:)) 

   !--- copy to a global scratch array
   Temp_Pix_Array = Ref_Ch1_Dark

element_loop:   Do Elem_Idx = 1, Num_Elements

line_loop:  DO Line_Idx = 1, Num_Lines

   !--- check for valid data
   if (Space_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle 
   if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle 

   ELem_Idx_1 = max(1,min(Num_Elements,Elem_Idx - N_box))
   ELem_Idx_2 = max(1,min(Num_Elements,Elem_Idx + N_box))
   Line_Idx_1 = max(1,min(Num_Lines,Line_Idx - N_box))
   Line_Idx_2 = max(1,min(Num_Lines,Line_Idx + N_box))

!  if (Land_Class_Local(Elem_Idx,Line_Idx) /= Land_Class_Allowed) cycle 
!  Ref_Ch1_Dark_Min = minval(Temp_Pix_Array(Elem_Idx_1:Elem_Idx_2,Line_Idx_1:Line_Idx_2))
!  if (Ref_Ch1_Dark_Min /= Missing_Value_Real4) then
!     Ref_Ch1_Dark(Elem_Idx,Line_Idx) = Ref_Ch1_Dark_Min
!  endif 


   !--- compute median
    call COMPUTE_MEDIAN(Ref_Ch1_Dark(Elem_Idx_1:Elem_Idx_2,Line_Idx_1:Line_Idx_2), &
                        Bad_Pixel_Mask(Elem_Idx_1:Elem_Idx_2,Line_Idx_1:Line_Idx_2), &
                        Z_Median,Z_Mean,Z_Std_Median)

   if (Z_Median /= Missing_Value_Real4) then
      Ref_Ch1_Dark(Elem_Idx,Line_Idx) = Z_Median
   endif 

end do line_loop

end do element_loop

   !--- wipe clean to a global scratch array
   Temp_Pix_Array =  0.0

end subroutine POST_PROCESS_GOES_DARK_COMPOSITE
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
subroutine DARK_COMPOSITE_CLOUD_MASK(Cloud_Mask)

   integer(kind=int1), dimension(:,:),  intent(inout):: Cloud_Mask
   integer:: Elem_Idx
   integer:: Line_Idx
   integer:: Num_Elements
   integer:: Num_Lines
   real, parameter:: Solzen_Max_Threshold = 60.0
   real, parameter:: Ref_Delta_Cloud = 10.0
   real, parameter:: Ref_Delta_Clear = 5.0
  
   Num_Elements = size(Cloud_Mask(:,1)) 
   Num_Lines = size(Cloud_Mask(1,:)) 
 
element_loop:   Do Elem_Idx = 1, Num_Elements

line_loop:  DO Line_Idx = 1, Num_Lines

   !--- check for valid data
   if (Space_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle 
   if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle 
   if (Solzen(Elem_Idx,Line_Idx) > Solzen_Max_Threshold) cycle 
   if (Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) cycle 
   if (ch(1)%Ref_Toa(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle 
   if (Ref_Ch1_Dark_Composite(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle 

   !---  if clear or prob clear, look for cloud
   if (Cloud_Mask(Elem_Idx,Line_Idx) == sym%CLEAR .or.      &
       Cloud_Mask(Elem_Idx,Line_Idx) == sym%PROB_CLEAR) then

       if (ch(1)%Ref_Toa(Elem_Idx,Line_Idx) - Ref_Ch1_Dark_Composite(Elem_Idx,Line_Idx) > Ref_Delta_Cloud) then
          Cloud_Mask(Elem_Idx,Line_Idx) = sym%PROB_CLOUDY
       endif

   endif

   !---  if cloudy or prob cloudy, look for clear
   if (Cloud_Mask(Elem_Idx,Line_Idx) == sym%CLOUDY .or.      &
       Cloud_Mask(Elem_Idx,Line_Idx) == sym%PROB_CLOUDY) then

       if (ch(1)%Ref_Toa(Elem_Idx,Line_Idx) - Ref_Ch1_Dark_Composite(Elem_Idx,Line_Idx) < Ref_Delta_Clear) then
          Cloud_Mask(Elem_Idx,Line_Idx) = sym%PROB_CLEAR
       endif

   endif


end do line_loop

end do element_loop

end subroutine DARK_COMPOSITE_CLOUD_MASK



! convert integer to real, but interpret it as an unsigned integer
! This is needed for the sounder
elemental function UNSIGNED_TO_REAL4(i) result (r)
  integer(kind=int2), intent(in) :: i
  real(kind=real4) :: r

  INTEGER(kind=int4), PARAMETER :: INT2_SIGN_CORRECTION_OFFSET = 65536
  
 
  if(i < 0) then
     r = real(i,kind=real4) + INT2_SIGN_CORRECTION_OFFSET
  else 
     r = real(i,kind=real4)
  end if
end function UNSIGNED_TO_REAL4





!
!-- end of module
!
end module GOES_MODULE
