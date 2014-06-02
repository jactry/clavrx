! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: sensor_module.f90 (src)
!       SENSOR_MODULE (program)
!
! PURPOSE: This module houses routines that apply to multiple sensors
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
! Routines in this module and their purpose:
!
!  the routines are called i process_ckavrx.f90 in this order
!   file_loop
!      detect_sensor_from_file
!      SET_FILE_DIMENSIONS(
!      SET_SENSOR_CONSTANTS 
!      READ_INSTR_CONSTANTS()
!      READ_ALGO_CONSTANTS()
!      READ_LEVEL1B_DATA
!       ...
!       ... 
! end loop
!--------------------------------------------------------------------------------------
module SENSOR_MODULE
 use PIXEL_COMMON
 use CALIBRATION_CONSTANTS
 use ALGORITHM_CONSTANTS
 use CONSTANTS
 use FILE_UTILITY
 use AVHRR_MODULE
 use GOES_MODULE
 use MODIS_MODULE
 use FY2_MODULE
 use COMS_MODULE
 use IFF_CLAVRX_BRIDGE , only : &
     READ_IFF_DATA &
     , READ_IFF_VIIRS_INSTR_CONSTANTS &
     , READ_IFF_DATE_TIME &
     , GET_IFF_DIMS_BRIDGE
 use MTSAT_MODULE
 use SEVIRI_MODULE
#ifdef HDF5LIBS
   use VIIRS_CLAVRX_BRIDGE , only : &
       READ_VIIRS_DATE_TIME &
       , READ_VIIRS_DATA &
       , GET_NUMBER_OF_SCANS_FROM_VIIRS &
       , READ_VIIRS_INSTR_CONSTANTS
#endif

  implicit none

  public:: SET_SENSOR_CONSTANTS
  public:: READ_INSTR_CONSTANTS
  public:: READ_ALGO_CONSTANTS
  public:: DETECT_SENSOR_FROM_FILE
  public:: SET_FILE_DIMENSIONS
  public:: READ_LEVEL1B_DATA 

  character(24), parameter, private :: MOD_PROMPT = " SENSOR_MODULE: "
  character(38) :: orbit_identifier
  
  contains
!==============================================================================
!
!==============================================================================
subroutine SET_SENSOR_CONSTANTS(AREAstr)

   TYPE (AREA_STRUCT), intent(in) :: AREAstr

   INTEGER:: Idx
   INTEGER(kind=int4):: Start_Year_Tmp
   INTEGER(kind=int4):: Start_Day_Tmp
   INTEGER(kind=int4):: End_Year_Tmp
   INTEGER(kind=int4):: End_Day_Tmp
   INTEGER(kind=int4):: Start_Time_Tmp
   INTEGER(kind=int4):: End_Time_Tmp
   INTEGER(kind=int4):: Orbit_Number_Tmp
   INTEGER(kind=int4):: Hour
   INTEGER(kind=int4):: Minute
   INTEGER(kind=int4):: Second
   INTEGER(kind=int4):: Avhrr_Number

if (Avhrr_Flag == sym%YES) then

   !-----------------------------------------------------------------
   !-- determine if this data is from the AVHRR/1 series
   !-----------------------------------------------------------------
   call DETERMINE_AVHRR_1(Start_Year,AVHRR_KLM_Flag,Avhrr_1_Flag)

   !-------------------------------------------------------------------
   !-- define a sc_Id that is unique value for each satellite
   !-------------------------------------------------------------------
   if (Avhrr_Flag == sym%YES) then
      call ASSIGN_AVHRR_SAT_ID_NUM_INTERNAL(Sc_Id,Avhrr_Number,Sc_Id_Char)
   endif

   !--------------------------------------------------------------
   !--- make an INTEGER Orbit_Number from proc_block_Id
   !--------------------------------------------------------------
    Orbit_Number = 0
    do Idx = 1,len(Proc_Block_Id)
      Orbit_Number = Orbit_Number + (ichar(Proc_Block_Id(Idx:Idx))-48)* &
       (10**(len(Proc_Block_Id) - Idx))
    enddo
 
   !--------------------------------------------------------------
   !--- based on spacecraft id, set up constants
   !--------------------------------------------------------------
   if (Avhrr_Number < 10) then
    Instr_Const_File = trim(Ancil_Data_Dir)//"avhrr_data/"// &
              "avhrr_"//char(Avhrr_Number+48)//"_instr.dat"
    Algo_Const_File = trim(Ancil_Data_Dir)//"avhrr_data/"// &
              "avhrr_"//char(Avhrr_Number+48)//"_algo.dat"
   else
    Instr_Const_File = trim(Ancil_Data_Dir)//"avhrr_data/"// &
             "avhrr_1"//char(Avhrr_Number - 10 + 48)//"_instr.dat" 
    Algo_Const_File = trim(Ancil_Data_Dir)//"avhrr_data/"// &
             "avhrr_1"//char(Avhrr_Number - 10 + 48)//"_algo.dat" 
   endif

endif

!----------------------------------------------
! for Modis, take time from file name
!----------------------------------------------
if (Modis_Flag == sym%YES) then

print *, "calling modis time attr"
  CALL READ_MODIS_TIME_ATTR(trim(dir_1b), trim(File_1b), Start_Year, &
                            Start_Day, Start_Time, End_Year, End_Day, &
                            End_Time)
  
print *, "returned from  modis time attr"
  if (Modis_Aqua_Flag == sym%YES .OR. Modis_Aqua_Mac_Flag == sym%YES) then
         Sc_Id_WMO = 784
         Instr_Const_File = 'modis_aqua_instr.dat'
         Algo_Const_File = 'modis_aqua_algo.dat'
         Platform_Name_Attribute = 'AQUA'
         Sensor_Name_Attribute = 'MODIS'
  endif
  if (Modis_Terra_Flag == sym%YES) then
         Sc_Id_WMO = 783
         Instr_Const_File = 'modis_terra_instr.dat'
         Algo_Const_File = 'modis_terra_algo.dat'
         Platform_Name_Attribute = 'TERRA'
         Sensor_Name_Attribute = 'MODIS'
  endif
  Instr_Const_File = trim(Ancil_Data_Dir)//"avhrr_data/"//trim(Instr_Const_File)
  Algo_Const_File = trim(Ancil_Data_Dir)//"avhrr_data/"//trim(Algo_Const_File)

! if (Modis_Terra_Flag == sym%YES) print *, EXE_PROMPT, "Satellite = Terra"
! if (Modis_Aqua_Flag == sym%YES) print *, EXE_PROMPT, "Satellite = Aqua"
! if (Modis_Aqua_Mac_Flag == sym%YES) print *, EXE_PROMPT, "Satellite = Aqua for ATRAIN subregion"

endif

!----------------------------------------------
! for VIIRS take time from file's name
!----------------------------------------------
if (Viirs_Flag == sym%YES) then
#ifdef HDF5LIBS
   call READ_VIIRS_DATE_TIME(trim(File_1b),Start_Year_Tmp,Start_Day_Tmp,Start_Time_Tmp, &
                             End_Time_Tmp,Orbit_Number_Tmp ,orbit_identifier, end_year_tmp , end_day_tmp)
   Start_Year = Start_Year_tmp
   End_Year = End_Year_tmp
   Start_Day = Start_Day_tmp
   End_Day = End_Day_tmp
 
   Start_Time = Start_Time_tmp
   End_Time = End_Time_tmp
   Orbit_Number = Orbit_Number_tmp
   Sc_Id_WMO = 224
   Instr_Const_File = 'viirs_npp_instr.dat'
   Algo_Const_File = 'viirs_npp_algo.dat'
   Platform_Name_Attribute = 'NPP'
   Sensor_Name_Attribute = 'VIIRS'
!  Sensor_Name_Attribute = 'NPP : VIIRS'
   
   Instr_Const_File = trim(Ancil_Data_Dir)//"avhrr_data/"//trim(Instr_Const_File)
   Algo_Const_File = trim(Ancil_Data_Dir)//"avhrr_data/"//trim(Algo_Const_File)

#else
        PRINT *, "No HDF5 libraries installed, stopping"
        stop
#endif

endif 

!----------------------------------------------
! for IFF take time and set some constants
! could be VIRRS or MODIS sensor
!----------------------------------------------
if (Iff_Viirs_Flag == sym%YES .or. Iff_Modis_Flag == sym%YES) then
   call READ_IFF_DATE_TIME(trim(File_1b),Start_Year_tmp,Start_Day_tmp,Start_Time_tmp, &
                      End_Year_tmp,End_Day_tmp,End_Time_tmp)
   Start_Year = Start_Year_tmp
   End_Year = End_Year_tmp
   Start_Day = Start_Day_tmp
   End_Day = End_Day_tmp
   Start_Time = Start_Time_tmp
   End_Time = End_Time_tmp
!   Orbit_Number = Orbit_Number_tmp
   if (Iff_Viirs_Flag == sym%YES) then
      Sc_Id_WMO = 224
      Instr_Const_File = 'iff_viirs_npp_instr.dat'
      Algo_Const_File = 'viirs_npp_algo.dat'
      Platform_Name_Attribute = 'NPP'
      Sensor_Name_Attribute = 'VIIRS'
   elseif (Iff_Modis_Flag == sym%YES) then
      Sc_Id_WMO = 784
      Instr_Const_File = 'modis_aqua_instr.dat'
      Algo_Const_File = 'modis_aqua_algo.dat'
      Platform_Name_Attribute = 'AQUA'      
      Sensor_Name_Attribute = 'MODIS'      
   endif

   Instr_Const_File = trim(Ancil_Data_Dir)//"avhrr_data/"//trim(Instr_Const_File)
   Algo_Const_File = trim(Ancil_Data_Dir)//"avhrr_data/"//trim(Algo_Const_File)

endif

!----------------------------------------------
! for GOES, MTSAT and MSG, take time from AREAstr
!----------------------------------------------
if (Goes_Flag == sym%YES .OR. Goes_Sndr_Flag == sym%YES .OR. &
    Seviri_Flag == sym%YES .OR. &
    FY2_Flag == sym%YES .OR. Mtsat_Flag==sym%YES .OR. Coms_Flag == sym%YES) then
    Start_Year = 1900 + int(AREAstr%act_img_Date/1000)
    End_Year = Start_Year
    Start_Day = AREAstr%act_img_Date - (Start_Year-1900)*1000
    End_Day = Start_Day
    hour = AREAstr%act_img_Time/10000 
    minute = (AREAstr%act_img_Time - hour*10000)/100
    second = (AREAstr%act_img_Time - hour*10000 - minute*100)/100
    Start_Time = ((hour * 60 + minute) * 60 + second) * 1000 !millisec
    End_Time = Start_Time
    
    !--- assign id numbers
    if (Goes_Flag == sym%YES) call ASSIGN_GOES_SAT_ID_NUM_INTERNAL(AREAstr%Sat_Id_Num,AREAstr%Elem_Res)

    if (Goes_Sndr_Flag == sym%YES) call ASSIGN_GOES_SNDR_ID_NUM_INTERNAL(AREAstr%Sat_Id_Num)
    
    if (Seviri_Flag == sym%YES) call ASSIGN_MSG_SAT_ID_NUM_INTERNAL(AREAstr%Sat_Id_Num)

    if (Mtsat_Flag == sym%YES) call ASSIGN_MTSAT_SAT_ID_NUM_INTERNAL(AREAstr%Sat_Id_Num)

    if (FY2_Flag == sym%YES) call ASSIGN_FY_SAT_ID_NUM_INTERNAL(AREAstr%Sat_Id_Num)

    if (Coms_Flag == sym%YES) call ASSIGN_COMS_ID_NUM_INTERNAL(AREAstr%Sat_Id_Num)

endif

!--- screen output
print *, EXE_PROMPT, "Satellite : Sensor = ", Platform_Name_Attribute, ' : ',Sensor_Name_Attribute
print *,EXE_PROMPT, "Spacecraft WMO number = ", Sc_Id_WMO
print *,EXE_PROMPT, "Start Day of Year= ", Start_Day
print *,EXE_PROMPT, "Start Time of Day = ", Start_Time /1000.0/ 60.0 / 60.0
print *,EXE_PROMPT, "Number of Scans = ", Num_Scans

if (Avhrr_Flag == sym%YES) then
    print *,EXE_PROMPT, "Sensor = AVHRR"
    print *,EXE_PROMPT, "spacecraft number = ", Avhrr_Number
    print *,EXE_PROMPT, "data type = ", Data_Type
    print *,EXE_PROMPT, "level 1b version = ", Ver_1b
    print *,EXE_PROMPT, "Gac flag = ", AVHRR_GAC_Flag
    print *,EXE_PROMPT, "AVHRR_KLM_Flag flag = ", AVHRR_KLM_Flag
    print *,EXE_PROMPT, "AAPP_Flag flag = ", AVHRR_AAPP_Flag
    print *,EXE_PROMPT, "avhrr/1 flag = ", Avhrr_1_Flag
endif

if (GOES_Flag == sym%YES .or. MTSAT_Flag == sym%YES .or.  &
    SEVIRI_Flag == sym%YES .or. FY2_FLAG == sym%YES .or. &
    COMS_Flag == sym%YES) then
    print *,EXE_PROMPT, "Pixel Resolution (km) = ", AREAstr%Elem_Res
endif


!--- set Resolution_KM for global attribute
Sensor_Resolution_KM = -999.0

if (Goes_Flag == sym%YES) then
     Sensor_Resolution_KM = 4.0
     if (GOES_1km_Flag == sym%YES) Sensor_Resolution_KM = 1.0
endif
if (MODIS_Flag == sym%YES) Sensor_Resolution_KM = 1.0
if (VIIRS_Flag == sym%YES) Sensor_Resolution_KM = 0.75
if (AVHRR_Flag == sym%YES) then
     Sensor_Resolution_KM = 1.1
     if (AVHRR_GAC_Flag == sym%YES) then
       Sensor_Resolution_KM = 4.0
     endif
endif
if (SEVIRI_Flag == sym%YES) Sensor_Resolution_KM = 3.0
if (COMS_Flag == sym%YES) Sensor_Resolution_KM = 4.0
if (FY2_Flag == sym%YES) Sensor_Resolution_KM = 4.0
if (MTSAT_Flag == sym%YES) Sensor_Resolution_KM = 4.0
if (GOES_Sndr_Flag == sym%YES) Sensor_Resolution_KM = 10.0

end subroutine SET_SENSOR_CONSTANTS

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
subroutine READ_INSTR_CONSTANTS()

   if (Avhrr_Flag == sym%YES) call READ_AVHRR_INSTR_CONSTANTS(trim(Instr_Const_File))

   if (Modis_Flag == sym%YES .or. Iff_Modis_Flag == sym%YES) call READ_MODIS_INSTR_CONSTANTS(trim(Instr_Const_File))

   if (Goes_Flag == sym%YES) call READ_GOES_INSTR_CONSTANTS(trim(Instr_Const_File))

   if (Goes_Sndr_Flag == sym%YES) call READ_GOES_SNDR_INSTR_CONSTANTS(trim(Instr_Const_File))

   if (Seviri_Flag == sym%YES) call READ_MSG_INSTR_CONSTANTS(trim(Instr_Const_File))

   if (Mtsat_Flag == sym%YES) call READ_MTSAT_INSTR_CONSTANTS(trim(Instr_Const_File))

   if (FY2_Flag == sym%YES) call READ_FY_INSTR_CONSTANTS(trim(Instr_Const_File))

   if (Coms_Flag == sym%YES) call READ_COMS_INSTR_CONSTANTS(trim(Instr_Const_File))

   if (Viirs_Flag == sym%YES) then

#ifdef HDF5LIBS !--- read in Viirs Instrument Constants from appropriate file
     call READ_VIIRS_INSTR_CONSTANTS(trim(Instr_Const_File))
#else
       print *, "No HDF5 library installed, stopping"
       stop
#endif
   endif

   if (iff_Viirs_Flag == sym%YES) call READ_IFF_VIIRS_INSTR_CONSTANTS(trim(Instr_Const_File))

end subroutine READ_INSTR_CONSTANTS
!----------------------------------------------------------------
! read the avhrr algorithm constants into memory
!-----------------------------------------------------------------
SUBROUTINE READ_ALGO_CONSTANTS()
 integer:: ios
 integer:: Algo_Lun

 Algo_Lun = GET_LUN()

 open(unit=Algo_Lun,file=trim(Algo_Const_File),status="old",position="rewind",action="read",iostat = ios)

  IF (ios /= 0) THEN
    print *, "Error opening algorithm constant file, iostat = ", ios
    stop
  ENDIF

  !nlSst
  read(unit=Algo_Lun,fmt=*)
! read(unit=Algo_Lun,fmt=*) a1_day,a2_day,a3_day,a4_day

  !triple Sst
  read(unit=Algo_Lun,fmt=*)
! read(unit=Algo_Lun,fmt=*)
! a1_night,a2_night,a3_night,a4_night,a5_night,a6_night

  !-- Olr
  read(unit=Algo_Lun,fmt=*) Win_0,Win_1,Win_2,Win_3
  read(unit=Algo_Lun,fmt=*) Olr_0,Olr_1,Olr_2,Olr_3,Olr_4

  close(unit=Algo_Lun)

END SUBROUTINE READ_ALGO_CONSTANTS

!--------------------------------------------------------------------------------------------------
!  Apply various tests to determine from which sensor this data comes
!--------------------------------------------------------------------------------------------------
subroutine DETECT_SENSOR_FROM_FILE(File_1b_Full,File_1b_Temp,AREAstr,NAVstr,Ierror)

  CHARACTER(len=*), intent(in) :: File_1b_Full
  CHARACTER(len=*), intent(in) :: File_1b_Temp
  TYPE (AREA_STRUCT), intent(out) :: AREAstr
  TYPE (GVAR_NAV), intent(out)    :: NAVstr
  INTEGER(kind=int4):: Ierror

  INTEGER(kind=int4):: Itest_Aqua_5km
  INTEGER(kind=int4):: Itest_Terra_5km
  INTEGER(kind=int4):: Itest_Aqua_1km
  INTEGER(kind=int4):: Itest_Terra_1km
  INTEGER(kind=int4):: Itest_Aqua_Mac
  INTEGER(kind=int4):: Itest_Aqua_CSPP
  INTEGER(kind=int4):: Itest_Terra_CSPP

  INTEGER :: itest_viirs
  INTEGER :: itest_iff_viirs
  INTEGER :: itest_iff_modis


  Ierror = sym%NO
  Goes_Flag = sym%NO
  Goes_Sndr_Flag = sym%NO
  Modis_Flag = sym%NO
  Avhrr_Flag = sym%NO
  Seviri_Flag = sym%NO
  FY2_Flag = sym%NO
  Mtsat_Flag = sym%NO
  Coms_Flag = sym%NO
  Modis_Aqua_Flag = sym%NO
  Modis_Terra_Flag = sym%NO
  Modis_1km_Flag = sym%NO
  Modis_5km_Flag = sym%NO
  Modis_Aqua_Mac_Flag = sym%NO
  Viirs_Flag = sym%NO
  Iff_Viirs_Flag = sym%NO
  Iff_Modis_Flag = sym%NO

  !-------------------------------------------------------------------------
  !--- First, test for a modis file
  !-------------------------------------------------------------------------
  Itest_Aqua_5km = index(File_1b_Temp, 'MYD02SSH')
  Itest_Aqua_1km = index(File_1b_Temp, 'MYD021KM')
  Itest_Terra_5km = index(File_1b_Temp, 'MOD02SSH')
  Itest_Terra_1km = index(File_1b_Temp, 'MOD021KM')
  Itest_Aqua_Mac = index(File_1b_Temp, 'MAC02')
  Itest_Aqua_CSPP = index(File_1b_Temp, 'a1.')
  Itest_Terra_CSPP = index(File_1b_Temp, 't1.')

  if (Itest_Aqua_1km == 1 .or. Itest_Aqua_5km == 1 .or. &
      Itest_Terra_1km == 1 .or. Itest_Terra_5km == 1 .or. &
      Itest_Aqua_Mac == 1 .or. Itest_Aqua_CSPP == 1 .or. &
      Itest_Terra_CSPP == 1) then 
       Modis_Flag = sym%YES
  endif

  if (Modis_Flag == sym%YES) then
          if (Itest_Aqua_5km == 1 .or. Itest_Aqua_1km == 1) Modis_Aqua_Flag = sym%YES

          if (Itest_Terra_5km == 1 .or. Itest_Terra_1km == 1) Modis_Terra_Flag = sym%YES

          if (Itest_Terra_5km == 1 .or. Itest_Aqua_5km == 1) Modis_5km_Flag = sym%YES

          if (Itest_Terra_1km == 1 .or. Itest_Aqua_1km == 1 .or.&
              Itest_Aqua_Mac ==1 .or. Itest_Aqua_CSPP == 1 .or. &
              Itest_Terra_CSPP == 1) Modis_1km_Flag = sym%YES

          if (Itest_Aqua_CSPP == 1 .or. &
              Itest_Terra_CSPP == 1) Modis_CSPP_Flag = sym%YES

          if (Itest_Aqua_Mac ==1 ) Modis_Aqua_Mac_Flag = sym%YES
  endif
  
  !-- for 1 km MODIS, determine name of separate geolocation file
  if (Modis_1km_Flag == sym%YES) then
     call DETERMINE_MODIS_GEOLOCATION_FILE(File_1b_Temp,Dir_Nav_in,Modis_Geo_Name)
          
     if (trim(Modis_Geo_Name) == "no_file") then
        Ierror = sym%YES
     endif
      
     if (Cloud_Mask_Aux_Flag /= sym%No_AUX_CLOUD_MASK) then
        call DETERMINE_MODIS_CLOUD_MASK_FILE(File_1b_Temp,Dir_1b,Modis_Cloud_Mask_Name)
        if (trim(Modis_Cloud_Mask_Name) == "no_file") then
          Ierror = sym%YES
        endif
     endif
  endif
 
  
  ! Set AQUA/Terra flags here for CSPP files
  
  if (Itest_Aqua_CSPP == 1) Modis_Aqua_Flag = sym%YES 
  if (Itest_Terra_CSPP == 1) Modis_Terra_Flag = sym%YES
    
  !-------------------------------------------------------------------------------------------------
  !--- Second, if not MODIS, test if this is a valid  McIdas Areafile.
  !--- If it is, assume it is GOES unless it's satellite id indicates otherwise
  !-------------------------------------------------------------------------------------------------
  if (Modis_Flag == sym%NO) then                                !begin MODIS test

    call GET_GOES_HEADERS(File_1b_Full, AREAstr, NAVstr)

    if (AREAstr%version_Num == 4) then                          !begin valid Areafile test

          Seviri_Flag = sym%NO
          FY2_Flag = sym%NO
          Goes_Flag = sym%YES
          Mtsat_Flag = sym%NO
          Coms_Flag = sym%NO
         
          select case(AREAstr%Sat_Id_Num)
 
          !test for SEVIRI
          case (51:53)
              Seviri_Flag = sym%YES
              Goes_Flag = sym%NO
              Goes_Sub_Satellite_Latitude = 0.0
              if (AREAstr%Sat_Id_Num == 51 ) Goes_Sub_Satellite_Longitude = -3.477996     ! Longitude of actual Sub-Satellite Point for Met-8
              if (AREAstr%Sat_Id_Num == 52 ) Goes_Sub_Satellite_Longitude = -0.159799     ! Longitude of actual Sub-Satellite Point for Met-9
              if (AREAstr%Sat_Id_Num == 53 ) Goes_Sub_Satellite_Longitude = -0.159799     ! Longitude of actual Sub-Satellite Point for Met-10 (?? TBD ??)
              AREAstr%cal_Offset = AREAstr%reserved(3)
                    
          !test for MTSAT
          case (84,85)
              Goes_Flag = sym%NO
              Mtsat_Flag = sym%YES
              !This is needed to determine type of navigation
              !as Nav coefficents specific to MTSAT (and FY2)
              CALL READ_NAVIGATION_BLOCK_MTSAT_FY(trim(File_1b_Full), AREAstr, NAVstr)
              Goes_Sub_Satellite_Latitude = NAVstr%sublat
              Goes_Sub_Satellite_longitude = NAVstr%sublon
          
!WCS3 - FY2-D AREA files have the subsat lat/lon flipped. Fix here
!        Fix with McIDAS-X is fixed.
          case (36, 37)
              Goes_Flag = sym%NO
              FY2_Flag = sym%YES
              !This is needed to determine type of navigation
              !as Nav coefficents specific to FY2D/E. They are stored in
              ! the same manner as MTSAT, hence using the same routine
              CALL READ_NAVIGATION_BLOCK_MTSAT_FY(trim(File_1b_Full), AREAstr, NAVstr)
              Goes_Sub_Satellite_Latitude = NAVstr%sublat
              Goes_Sub_Satellite_Longitude = NAVstr%sublon


          !test for COMS
          case (250)
              Goes_Flag = sym%NO
              Coms_Flag=sym%YES
              !This is needed to determine type of navigation
              !as Nav coefficents specific to COMS
              CALL READ_NAVIGATION_BLOCK_COMS(trim(File_1b_Full), AREAstr,NAVstr)
              Goes_Sub_Satellite_Latitude = NAVstr%sublat
              Goes_Sub_Satellite_Longitude = NAVstr%sublon
                                          
          !Test for the Goes Sounder
          case (71,73,75,77,79,181,183,185)
              Goes_Flag = sym%NO
              Goes_Sndr_Flag = sym%YES
          end select
          
          if (Goes_Flag == sym%YES .or. Goes_Sndr_Flag == sym%YES) THEN
                CALL LMODEL(Goes_Input_Time,  &
                      Goes_Epoch_Time, &
                      NAVstr, &
                      Goes_Sub_Satellite_Latitude, &
                      Goes_Sub_Satellite_Longitude)
                      Goes_Sub_Satellite_Longitude = Goes_Sub_Satellite_Longitude / Dtor
                      Goes_Sub_Satellite_Latitude = Goes_Sub_Satellite_Latitude  / Dtor
          endif

    endif                                                        !end Valid Areafile test
  endif                                                          !end MODIS test

  !--------------------------------------------------------------------------------------
  !--- Third, Test for a IDPS VIIRS data set 
  !--------------------------------------------------------------------------------------
  itest_viirs = index(File_1b_Temp, 'GMTCO')
  if (itest_viirs == 1) then
      Viirs_Flag = sym%YES
  end if

  !--------------------------------------------------------------------------------------
  !--- Fourth, Test for a PEATE Intermediate File Format (IFF) VIIRS or MODIS data set 
  !--------------------------------------------------------------------------------------
  itest_iff_viirs = index(File_1b_Temp, 'IFFSDR_npp')
  itest_iff_modis = index(File_1b_Temp, 'IFFSDR_aqua')
  if (itest_iff_viirs == 1) then
      Iff_Viirs_Flag = sym%YES
  elseif (itest_iff_modis == 1) then
      Iff_Modis_Flag = sym%YES
  endif

  !--------------------------------------------------------------------------------------
  !--- Fifth, if none of the above, set to AVHRR
  !--------------------------------------------------------------------------------------
  if (Modis_Flag == sym%NO .and.  &
      Goes_Flag == sym%NO .and. Goes_Sndr_Flag == sym%NO .and. &
      Seviri_Flag == sym%NO .and. Mtsat_Flag == sym%NO .and. &
      Viirs_Flag == sym%NO .and. Iff_Viirs_Flag == sym%NO .and. &
      Iff_Modis_Flag == sym%NO .and. FY2_Flag == sym%NO .and. &
      Coms_Flag == sym%NO) then
          AVHRR_Flag = sym%YES
  endif

  end subroutine DETECT_SENSOR_FROM_FILE

!---------------------------------------------------------------------------------------------
! Determine the number of elements (Num_Pix) and Number of Scans (Num_Scans)
! expected.  Also, 
!---------------------------------------------------------------------------------------------
  subroutine SET_FILE_DIMENSIONS(File_1b_Full,AREAstr,Nrec_Avhrr_Header,Nword_Clavr,Nword_Clavr_Start)

  CHARACTER(len=*), intent(in) :: File_1b_Full
  TYPE (AREA_STRUCT), intent(in) :: AREAstr
  INTEGER:: Erstat
  INTEGER(kind=int4), intent(out):: Nword_Clavr
  INTEGER(kind=int4), intent(out):: Nword_Clavr_Start
  INTEGER(kind=int4), intent(out):: Nrec_Avhrr_Header
  INTEGER(kind=int4):: Ierror_Viirs_Nscans

  character(len=150):: Dir_File

   if ((Modis_1km_Flag == sym%YES) .or. (Modis_5km_Flag == sym%YES) )then
        CALL READ_MODIS_SIZE_ATTR(trim(File_1b_Full),Num_Pix,Num_Scans)
   endif
   
   if (Modis_Aqua_Mac_Flag == sym%YES) then
      Num_Pix =  11
      Num_Scans = 2030
   endif
   
   

   if (Viirs_Flag == sym%YES) then
      
      Num_Pix = 3200
      Dir_File = trim(Dir_1b) // trim(File_1b)
#ifdef HDF5LIBS
      call GET_NUMBER_OF_SCANS_FROM_VIIRS(trim(Dir_File),Num_Scans,Ierror_Viirs_Nscans)

      ! If error reading, then go to next file
      if (ierror_viirs_nscans /= 0) then
              print *, EXE_PROMPT, "Error reading VIIRS, skipping this file"
              return      !AKH  - WHAT SHOULD HAPPEN HERE?
      endif

      ! Check if VIIRS Number of scans is regular (48) and calculate Number of y pixels
      if (Num_Scans == 48) then
        Num_Scans = Num_Scans * 16      !16pix per Scan
      elseif (Num_Scans == 47) then
        Num_Scans = (Num_Scans+1) * 16
      else
        PRINT *,"!!!!!!!!!!!! Abnormal number of Scan_Lines, skipping this file !!!!!!!!!!!!", Num_Scans
        return      !AKH  - WHAT SHOULD HAPPEN HERE?
      endif

#else
       print *, "No HDF5 library installed. VIIRS unable to process. Stopping"
       stop
#endif

   endif
      
!   if (Iff_Viirs_Flag == sym%YES) then
!      Num_Pix = 3200
!      Num_Scans = 2704
!   elseif (Iff_Modis_Flag == sym%YES) then
!      Num_Pix = 1354
!      Num_Scans = 2020
   !--- if an IFF, call routine to determine dimensions from latitude sds
   if (Iff_Viirs_Flag == sym%YES .or. IFF_Modis_Flag == sym%YES) then
      call GET_IFF_DIMS_BRIDGE(trim(Dir_1b)//trim(File_1b),Num_Pix,Num_Scans)
   endif
     
   if (Avhrr_Flag == sym%YES) then

   !-------------------------------------------------------
   ! Determine the type of level 1b file
   !-------------------------------------------------------
   call DETERMINE_AVHRR_FILE_TYPE(trim(File_1b_Full),AVHRR_GAC_FLAG,AVHRR_KLM_Flag,AVHRR_AAPP_Flag, &
                                  Ver_1b,Data_Type,Byte_Swap_1b)

   !--- based on AVHRR_KLM_Flag and Gac flags, check for option consistency
      if ((AVHRR_KLM_Flag == sym%NO) .and. (Cloud_Mask_Aux_Flag /= sym%NO_AUX_CLOUD_MASK)) then
        Erstat = 2
        print *, EXE_PROMPT,  &
           "Error: using cloud mask from 1bx option not "// & 
           "available for pre-AVHRR_KLM_Flag data, stopping"
        stop 2
      endif
      if ((AVHRR_KLM_Flag == sym%NO) .and. (BX_File_Flag == sym%YES)) then
        Erstat = 3
        print *,EXE_PROMPT, "Error: 1bx option not available "// & 
                            "for pre-AVHRR_KLM_Flag data, stopping"
        stop 3
     endif

    !-------------------------------------------------------------------
    !-- based on file type (AVHRR_KLM_Flag and Gac), determine parameters needed
    !-- to read in header and data records for this orbit
    !------------------------------------------------------------------- 
    call DEFINE_1B_DATA(AVHRR_GAC_Flag,AVHRR_KLM_Flag,AVHRR_AAPP_Flag,Nrec_Avhrr_Header,Nword_Clavr_Start,Nword_Clavr)

    !-------------------------------------------------------------------
    !-- read in header
    !-------------------------------------------------------------------
    call READ_AVHRR_LEVEL1B_HEADER(trim(File_1b_Full))


   endif

   !------------------------------------------------------------------------
   !  if GOES, SEVIRI and MTSAT, use elements of AREAstr to determine filesize
   !------------------------------------------------------------------------
   if (Goes_Flag == sym%YES) then
     Num_Pix =  int(AREAstr%Num_Elem / Goes_Xstride)
     Num_Scans = AREAstr%Num_Line
     L1b_Rec_Length = AREAstr%Num_Byte_Ln_Prefix +  &
                     (AREAstr%Num_Elem*AREAstr%Bytes_Per_Pixel)
   endif

   if (Goes_Sndr_Flag == sym%YES) then
     Num_Pix =  int(AREAstr%Num_Elem / Goes_Sndr_Xstride)
     Num_Scans = AREAstr%Num_Line
     L1b_Rec_Length = AREAstr%Num_Byte_Ln_Prefix +  &
                     (AREAstr%Num_Elem*AREAstr%Bytes_Per_Pixel)
   endif

   if (Seviri_Flag == sym%YES) then
     Num_Pix =  int(AREAstr%Num_Elem)
     Num_Scans = AREAstr%Num_Line
   endif

   if (Mtsat_Flag == sym%YES .OR. FY2_Flag == sym%YES .OR. Coms_Flag == sym%YES) then
     Num_Pix =  int(AREAstr%Num_Elem)
     Num_Scans = AREAstr%Num_Line
   endif

   end subroutine SET_FILE_DIMENSIONS

!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
 subroutine READ_LEVEL1B_DATA(File_1b_Full,Segment_Number,Time_Since_Launch,AREAstr,NAVstr,Nrec_Avhrr_Header,Ierror_Level1b)

     character(len=*), intent(in):: File_1b_Full
     integer, intent(in):: Segment_Number
     integer, intent(in):: Nrec_Avhrr_Header
     TYPE (AREA_STRUCT), intent(in) :: AREAstr
     TYPE (GVAR_NAV), intent(in)    :: NAVstr
     real, intent(in):: Time_Since_Launch
     integer, intent(out):: Ierror_Level1b

     Ierror_Level1b = 0
     Cloud_Mask_Aux_Read_Flag = sym%NO

     if (Modis_Flag == sym%YES) then
       call READ_MODIS(Segment_Number,Ierror_Level1b)
       if (Ierror_Level1b /= 0) return
     endif

     if (Goes_Flag == sym%YES) then
       call READ_GOES(Segment_Number,File_1b, &
                     Start_Day, Start_Time, &
                     Time_Since_Launch, &
                     AREAstr,NAVstr)

      if (Chan_On_Flag_Default(1)==sym%YES) then
       call READ_DARK_COMPOSITE_COUNTS(Segment_Number, Goes_Xstride, &
                     Dark_Composite_Name,AREAstr,Two_Byte_Temp) 
       call CALIBRATE_GOES_DARK_COMPOSITE(Two_Byte_Temp,Time_Since_Launch,Ref_Ch1_Dark_Composite)
      endif

     endif

     if (Goes_Sndr_Flag == sym%YES) then
       call READ_GOES_SNDR(Segment_Number,File_1b, &
                     Start_Day, Start_Time, &
                     Time_Since_Launch, &
                     AREAstr,NAVstr)
     endif

    !--------  MSG/SEVIRI
    if (Seviri_Flag == sym%YES) then
       call READ_SEVIRI(Segment_Number,File_1b, &
                     Start_Day, Start_Time, &
                     Time_Since_Launch, &
                     AREAstr,NAVstr)
      call READ_DARK_COMPOSITE_COUNTS(Segment_Number,Seviri_Xstride, &
                     Dark_Composite_Name,AREAstr,Two_Byte_Temp) 
      call CALIBRATE_SEVIRI_DARK_COMPOSITE(Two_Byte_Temp,Ref_Ch1_Dark_Composite)
    endif


    !--- MTSAT
    if (Mtsat_Flag == sym%YES) then
      call READ_MTSAT(Segment_Number,File_1b, &
                     Start_Day, Start_Time, &
                     Time_Since_Launch, &
                     AREAstr,NAVstr)
      call READ_DARK_COMPOSITE_COUNTS(Segment_Number,Mtsat_Xstride, &
                     Dark_Composite_Name,AREAstr,Two_Byte_Temp) 
      call CALIBRATE_MTSAT_DARK_COMPOSITE(Two_Byte_Temp,Ref_Ch1_Dark_Composite)
    endif

    !--- FY2
    if (FY2_Flag == sym%YES) then
      call READ_FY(Segment_Number,File_1b, &
                     Start_Day, Start_Time, &
                     Time_Since_Launch, &
                     AREAstr,NAVstr)
    endif

    !--- COMS
    if (Coms_Flag == sym%YES) then
      call READ_COMS(Segment_Number,File_1b, &
                     Start_Day, Start_Time, &
                     Time_Since_Launch, &
                     AREAstr,NAVstr)
    endif

      
    if (Avhrr_Flag == sym%YES) then
      call READ_AVHRR_LEVEL1B_DATA(trim(File_1b_Full), &
              AVHRR_KLM_Flag,AVHRR_AAPP_Flag,Therm_Cal_1b,&
              Time_Since_Launch,Nrec_Avhrr_Header,Segment_Number)
    endif

    if(Viirs_Flag == sym%Yes) then
#ifdef HDF5LIBS
      call READ_VIIRS_DATA ( segment_number , trim(File_1b),Ierror_Level1b)
      
      ! If error reading, then go to next file
      if (Ierror_Level1b /= 0) return

  
#else
       print *, "No HDF5 library installed, stopping"
       stop
#endif
    endif

    if(Iff_Viirs_Flag == sym%Yes .or. Iff_Modis_Flag == sym%Yes) then
       IFF_File = trim(File_1b)
        call READ_IFF_DATA (Segment_Number, trim(File_1b_Full), Ierror_Level1b)

      ! If error reading, then go to next file
      if (Ierror_Level1b /= 0) return

    endif

 end subroutine READ_LEVEL1B_DATA


end module SENSOR_MODULE
