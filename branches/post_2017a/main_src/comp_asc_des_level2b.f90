! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: comp_asc_des_level2b.f90 (src)
!       COMPILE_ASC_DES_LEVEL2B (program)
!
! PURPOSE: A main code that generates one of the executables in the CLAVR-x processing
!          system
!
! DESCRIPTION: this code takes the level2 files created for each orbit
!              for one day from one satellite and writes separate level2b 
!              files for the ascending and descending nodes. This code runs
!              after clavrxorb has processed the level-1b files for one day.
!              This program reads input (directories and filenames) from a file
!              called comp_asc_des_level2b_input.  Currently, there are no
!              command line arguments.  
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!
! (c) This code is copyrighted by the author and all NOAA restrictions apply
!
! Revision History:
!     May 2009: created
!     June 2010: Added options to use geo and nav files for speed   
!     August 2010: Added FCDR attributes
!     September 2010: Added descriptive global attributes
!     February 2011: Added support for GOES satellites 
!     November 2013:  Added time dimension for NCDC compliance
!
! Format and meaning of command line arguments
!    1: node option, values are "asc", "des", "geo" or "zen"
!    2: spatial sampling option, values are "near" or "rand"
!    3: orbital overlap option, values are "nadir_overlap" or "random_overlap"
! 
! Format of Required Input File:
!    line 1: directory of level2 files (input)
!    line 2: directory of level2b files (output)
!    line 3: year
!    line 4: julian day
!    line 5: time and time window in hours (only for zen node, ignored for others)
!    line 6: satellite number (i.e. 18 = NOAA-18)
!    line 7: longitudes (west, east, spacing)
!    line 8: latitudes (south, north, spacing)
!    line 9+: level2 file name (one per line)
!
!
! Reference:  Level2 and Level2b file contents are given on
!             http://cimss.ssec.wisc.edu/patmosx
!
! Dependencies:  (The following are names of modules)
!    CONSTANTS
!    HDF
!    HDF_PARAMS
!    NUMERICAL_ROUTINES
!    SCALING_PARAMETERS
!    LEVEL2B_ROUTINES
!    FILE_UTILITY
!
! Calling Sequence:  
!    comp_asc_des_level2 node
!
!    where:   
!          node = asc or des or geo 
!
! Public Routines Within This Program:
!     None
!
!--------------------------------------------------------------------------------------
program COMPILE_ASC_DES_LEVEL2B
 use CONSTANTS
 use HDF
 use HDF_PARAMS
 use NUMERICAL_ROUTINES
 use SCALING_PARAMETERS
 use LEVEL2B_ROUTINES 
 use FILE_UTILITY

 implicit none

 integer, parameter:: N_Files_Max = 1200
 integer(kind=int4):: N_Command_Line_Args
 integer, parameter:: Source_Length_Max = 1000

 integer(kind=int4):: First_Valid_Input
 integer(kind=int4):: Sd_Id_Input
 integer(kind=int4):: Sd_Id_Output

 integer(kind=int4), dimension(2):: Sds_Dims_Output_XY
 integer(kind=int4), dimension(1):: Sds_Dims_Output_X
 integer(kind=int4), dimension(1):: Sds_Dims_Output_Y
 integer(kind=int4), dimension(1):: Sds_Output_Start_X
 integer(kind=int4), dimension(1):: Sds_Output_Stride_X
 integer(kind=int4), dimension(1):: Sds_Output_Start_Y
 integer(kind=int4), dimension(1):: Sds_Output_Stride_Y
 integer(kind=int4), dimension(2):: Sds_Output_Start_XY
 integer(kind=int4), dimension(2):: Sds_Output_Stride_XY
 integer(kind=int4), dimension(2):: Sds_Dims
 logical, dimension(:,:), allocatable:: Mask_Output
 real(kind=real4), dimension(:), allocatable:: Input_Array_1d
 real(kind=real4), dimension(:), allocatable:: Output_Array_1d

 real(kind=real4), dimension(:,:), allocatable:: Scaled_Sds_Data_Input
 real(kind=real4), dimension(:,:), allocatable:: Scaled_Sds_Data_Output
 real(kind=real4), dimension(:), allocatable:: Scaled_Lat_Output_1d
 real(kind=real4), dimension(:), allocatable:: Scaled_Lon_Output_1d
 real(kind=real4), dimension(:,:), allocatable:: Unscaled_Sds_Data_Output
 real(kind=real4), dimension(:,:,:), allocatable:: Scaled_Sds_Data_Output_Full


 character(len=1020),dimension(N_Files_Max):: File_Input
 character(len=1020):: File_Output
 character(len=1020):: Dir_In
 character(len=1020):: Dir_Out
 character(len=1020):: Root_Name
 character(len=1020):: Temp_Name
 integer:: Comp_Asc_Des_Input_Lun
 integer:: Ielem
 integer:: Iline
 integer:: Ilon
 integer:: Ilat
 integer:: Ifile
 integer:: Isds
 integer:: Ipoint
 integer:: Itime

 integer:: N_Files
 integer:: Asc_Des_Node
 character(len=12):: Sensor_Name_Output
 integer(kind=int2):: Sc_Id_Output
 integer(kind=int4):: Sc_Id_Input
 integer:: Istatus_Sum
 integer:: ios
 integer(kind=int4):: Jday
 integer(kind=int4):: Year
 integer(kind=int4):: time
 real(kind=real4):: hour
 real(kind=real4)::  minute
 integer(kind=int4):: Time_Window
 real(kind=real4):: Start_Time_Window
 real(kind=real4):: End_Time_Window
 real(kind=real4):: Start_Date_Window
 real(kind=real4):: End_Date_Window
 real(kind=real4):: Start_Time_Input
 real(kind=real4):: End_Time_Input
 real(kind=real4):: Start_Date_Input
 real(kind=real4):: End_Date_Input
 integer(kind=int2):: Start_Day_Input
 integer(kind=int2):: End_Day_Input
 integer(kind=int2):: Start_Year_Input
 integer(kind=int2):: End_Year_Input
 integer(kind=int4):: Num_Elements_Input
 integer(kind=int4):: Num_Lines_Input
 integer(kind=int4):: Num_Sds_Input
 integer(kind=int4):: Num_Global_Attrs
 character(len=40), dimension(2):: Sds_2d_dim_string

 real(kind=real4):: Lon_West
 real(kind=real4):: Lon_East
 real(kind=real4):: dlon_Output
 real(kind=real4):: Lat_North
 real(kind=real4):: Lat_South
 real(kind=real4):: dlat_Output
 integer(kind=int4):: Ntime_Output
 integer(kind=int4):: Nlon_Output
 integer(kind=int4):: Nlat_Output
 integer(kind=int4):: num_points
 integer(kind=int4):: Dateline_Flag
 real(kind=real4), dimension(:), allocatable:: lon_Output_1d
 real(kind=real4), dimension(:), allocatable:: lat_Output_1d
 integer:: ipos

 type(Sds_Struct), dimension(0:300) :: sds
 type(Sds_Struct) :: Sds_Lat
 type(Sds_Struct) :: Sds_Lon
 type(Sds_Struct) :: Sds_Time
 type(Sds_Struct) :: Sds_Scan_Line_Number
 type(Sds_Struct) :: Sds_Scan_Element_Number
 type(Sds_Struct) :: Sds_Asc_Des_Flag
 type(Sds_Struct) :: Sds_Temp
 type(Sds_Struct) :: Sds_Lat_Out
 type(Sds_Struct) :: Sds_Lon_Out

 character(len=4):: Year_String
 character(len=3):: Jday_String
 character(len=2):: Month_String
 character(len=2):: Day_String
 character(len=3):: Node_String
 character(len=14):: Spatial_Option_String
 character(len=14):: Overlap_Option_String
 character(len=10):: Geo_Option_String
 character(len=4):: Time_String
 character(len=2):: Time_Win_String
 character(len=36)::  Date_Created_String

 integer(kind=int4), dimension(:), allocatable:: Input_Update_Index
 integer(kind=int4), dimension(:), allocatable:: Output_Update_Index

 integer(kind=int1), dimension(:,:), allocatable:: Update_Output_Mask
 real(kind=real4), dimension(:,:), allocatable:: Overlap_Random_Output

 real(kind=real4), dimension(:), allocatable:: time_1d_Input
 real(kind=real4), dimension(:,:), allocatable:: Bad_Pixel_Mask_Input
 real(kind=real4), dimension(:), allocatable:: Asc_Des_Flag_Input
 real(kind=real4), dimension(:), allocatable:: Scan_Line_Number_1d_Input

 real(kind=real4), dimension(:,:), allocatable:: Date_Input
 real(kind=real4), dimension(:,:), allocatable:: Time_Input
 real(kind=real4), dimension(:,:), allocatable:: Scan_Line_Number_Input

 real(kind=real4), dimension(:,:), allocatable:: Solzen_Input
 real(kind=real4), dimension(:,:), allocatable:: Senzen_Input
 real(kind=real4), dimension(:,:), allocatable:: Senzen_Output
 real(kind=real4), dimension(:,:), allocatable:: Solzen_Output
 integer(kind=int4), dimension(:,:), allocatable:: Temp_Mask_Output
 real(kind=real4), dimension(:,:), allocatable:: Time_Output
 real(kind=real4), dimension(:,:), allocatable:: Scan_Line_Number_Output
 real(kind=real4), dimension(:,:), allocatable:: Scan_Element_Number_Output
 real(kind=real4), dimension(:,:), allocatable:: Bad_Pixel_Mask_Output
 integer(kind=int2), dimension(:,:), allocatable:: Element_Number_Output

 real:: Processing_Time_Start_Hours
 real:: Processing_Time_End_Hours
 INTEGER, parameter:: num_Segment_Time_points=2
 REAL(kind=REAL4), dimension(num_Segment_Time_Points):: Segment_Time_Point_Seconds
 REAL(kind=REAL4) :: Start_Time_Point_Hours
 REAL(kind=REAL4) :: End_Time_Point_Hours
 integer:: sfstart
 integer:: sffattr
 integer:: sfrnatt
 integer:: sfsnatt
 integer:: sfscatt
 integer:: sffinfo
 integer:: sfginfo
 integer:: sfselect

 integer:: sfendacc
 integer:: sfend

 real(kind=real4):: Domain_Data_Fraction       !percentage of domain with data
 real(kind=real4):: Domain_Valid_Data_Fraction !percentage of domain with data
 integer(kind=int4):: Count_Total              !number of domain points with data
 integer(kind=int4):: Count_Valid              !number of domain points with valid data
 real(kind=real4):: Domain_Day_Fraction        !percentage of domain with daytime data
 integer(kind=int4):: Count_Day                !number of domain points with daytime data
 real(kind=real4):: Domain_Cloud_Fraction      !percentage of domain with cloudy data
 real(kind=real4):: Domain_Ice_Cloud_Fraction  !percentage of domain with ice-phase cloudy data
 integer(kind=int4):: Count_Cloud              !number of domain points with cloudy data
 real(kind=real4):: Domain_Mean_Day_Reflectance_0_65um_nom
 real(kind=real4):: Domain_Mean_Day_Reflectance_0_86um_nom
 real(kind=real4):: Domain_Mean_Day_Reflectance_1_60um_nom
 real(kind=real4):: Domain_Mean_Temperature_3_75um_nom
 real(kind=real4):: Domain_Mean_Temperature_11_0um_nom
 real(kind=real4):: Domain_Mean_Temperature_12_0um_nom
 real(kind=real4):: Domain_Mean_Acha_Success_Fraction
 real(kind=real4):: Domain_Mean_Dcomp_Success_Fraction

 character(len=22), parameter:: LOCAL_EXE_PROMPT = "COMP_ASC_DES_LEVEL2B: "

 real(kind=real4), parameter:: Senzen_Thresh = 5.0
 integer(kind=int4):: Spatial_Flag
 integer(kind=int4):: Overlap_Flag

 integer(kind=int4):: Temp_Update_Flag

 real(kind=real4):: xrand

! Some variables for CF compliance required for NCDC delivery
 integer, parameter:: Epoch_Date = 2440588     !Julian date for January 1st, 1970 (used for time calc)
 integer:: Days_Since_Epoch    !Days since January 1st, 1970
 character(len=100) :: Time_Coverage_Start
 character(len=100) :: Time_Coverage_End
 character(len=100) :: L1b_Input
 character(len=Source_Length_Max) :: Source_String
 integer:: Source_Length_Start
 integer:: Source_Length_End
 integer:: L1b_Length
 
 character(len=100) :: Id_String
 character(len=10), dimension(3):: ctime
 integer:: L0,N0,I0,J0,K0
 integer:: Geo_2d_Flag

 character(len=100):: Title_String
 character(len=100):: Calibration_String
 character(len=100):: Product_Version_String
 character(len=100):: Status_String
 character(len=100):: Institution_String
 character(len=100):: Program_String
 character(len=300):: Summary_String
 character(len=200):: Variable_String
 character(len=500):: Keywords_String
 character(len=200):: Keywords_Vocabulary_String
 character(len=100):: Time_Coverage_Resolution_String
 character(len=100):: Metadata_Link_String
 character(len=100):: Spatial_Resolution_String

 INCLUDE 'version.inc'

!----------------------------------------------------------------------
! Begin Executable Code
!----------------------------------------------------------------------
 Processing_Time_Start_Hours = COMPUTE_TIME_HOURS()
 Segment_Time_Point_Seconds = 0.0

 Spatial_Option_String = 'near'
 Spatial_Flag = 0

 Overlap_Option_String = 'nadir_overlap'   !'random_overlap'
 Overlap_Flag = 0

 Geo_Option_String = '1d_lat_lon' 
 Geo_2d_Flag = 0

!----------------------------------------------------------------------
! command line arguments
!----------------------------------------------------------------------
 N_Command_Line_Args = iargc()

 if (N_Command_Line_Args == 0) then
     print *, "comp_asc_des_level2b called with no command line arguments, printing help information"
     print *, " "
     print *, "three command line arguments are expected, last two are optional"
     print *, "argument 1 = node selection (asc, des, geo or zen)"
     print *, "argument 2 = spatial sampling option (near or rand) default is near"
     print *, "argument 3 = orbital overlap sampling option (nadir_overlap or randon_overlap) default is nadir_overlap"
     print *, "argument 4 = optional argument, 1d_lat_lon or 2d_lat_lon"
     print *, " "
     print *, "other options are read in from a file names comp_asc_des_level2b_input"
     print *, "line 1: directory of level2 files (input)"
     print *, "line 2: directory of level2b files (output)"
     print *, "line 3: year:"
     print *, "line 4: julian day"
     print *, "line 5: time and time window in hours (only for zen node, ignored for others)"
     print *, "line 6: satellite name"
     print *, "line 7: longitudes (west, east, spacing)"
     print *, "line 8: latitudes (south, north, spacing)"
     print *, "line 9+: level2 file name (one per line)"
     stop
 endif

 Geo_Option_String = '1d_lat_lon'
 if (N_Command_Line_Args == 1) then 
     call getarg(1, Node_String)
 elseif (N_Command_Line_Args == 2) then 
     call getarg(1, Node_String)
     call getarg(2, Spatial_Option_String)
 elseif (N_Command_Line_Args == 3) then 
     call getarg(1, Node_String)
     call getarg(2, Spatial_Option_String)
     call getarg(3, Overlap_Option_String)
 elseif (N_Command_Line_Args == 4) then 
     call getarg(1, Node_String)
     call getarg(2, Spatial_Option_String)
     call getarg(3, Overlap_Option_String)
     call getarg(4, Geo_Option_String)
 else
     print *, LOCAL_EXE_PROMPT, "COMP_ASC_DES_LEVEL2b ERROR:: Unexpected Number of Command Line Arguments"
 endif
 if (Node_String /= "asc" .and. Node_String /= "des" .and. Node_String /= "geo" .and. Node_String /= "zen") then
     print *, LOCAL_EXE_PROMPT, "COMP_ASC_DES_LEVEL2b ERROR:: First node argument is not asc, des, geo or zen"
     stop 
 endif
 if (trim(Spatial_Option_String) /= "near" .and. trim(Spatial_Option_String) /= "rand") then
     print *, LOCAL_EXE_PROMPT, "COMP_ASC_DES_LEVEL2b ERROR:: Second node argument is not near or rand"
     stop 
 endif
 if (trim(Overlap_Option_String) /= "nadir_overlap" .and. trim(Overlap_Option_String) /= "random_overlap") then
     print *, LOCAL_EXE_PROMPT, "COMP_ASC_DES_LEVEL2b ERROR:: Third node argument is not nadir_overlap or random_overlap"
     print *, trim(Overlap_Option_String)
     stop 
 endif
 if (trim(Geo_Option_String) /= "1d_lat_lon" .and. trim(Geo_Option_String) /= "2d_lat_lon") then
     print *, LOCAL_EXE_PROMPT, "COMP_ASC_DES_LEVEL2b ERROR:: Fourth argument is not 1d_lat_lon or 2d_lat_lon"
     print *, trim(Overlap_Option_String)
     stop 
 endif

 if (trim(Spatial_Option_String) == 'rand') then
    Spatial_Flag = 1
 endif

 if (trim(Overlap_Option_String) == 'random_overlap') then
    Overlap_Flag = 1
 endif

 if (Spatial_Flag == 1 .or. Overlap_Flag == 1) then
    call INIT_RANDOM_SEED()
 endif

 if (trim(Geo_Option_String) == '2d_lat_lon') then
    Geo_2d_Flag = 1
 endif

 Sds_2d_Dim_String(1) = "longitude_index"
 Sds_2d_Dim_String(2) = "latitude_index"

!-----------------------------------------------------
! open input file
!----------------------------------------------------
Comp_Asc_Des_Input_Lun = GET_LUN()
open(unit=Comp_Asc_Des_Input_Lun,file="comp_asc_des_level2b_input",status="old",action="read",position="rewind",iostat=ios)
if (ios /= 0) then
  print *, LOCAL_EXE_PROMPT, "error opening comp_asc_dec control file, ios = ", ios
  stop 1
else
  read(unit=Comp_Asc_Des_Input_Lun,fmt="(A)") dir_in
  read(unit=Comp_Asc_Des_Input_Lun,fmt="(A)") dir_out
  read(unit=Comp_Asc_Des_Input_Lun,fmt=*) Year
  read(unit=Comp_Asc_Des_Input_Lun,fmt=*) Jday
  read(unit=Comp_Asc_Des_Input_Lun,fmt=*) Time, Time_Window
  read(unit=Comp_Asc_Des_Input_Lun,fmt=*) Sensor_Name_Output
  read(unit=Comp_Asc_Des_Input_Lun,fmt=*) Lon_West, Lon_East, Dlon_Output
  read(unit=Comp_Asc_Des_Input_Lun,fmt=*) Lat_South, Lat_North, Dlat_Output
endif


!--- calculate days since January 1st, 1970
Days_Since_Epoch = 1 - 32075 + 1461*(Year+4800+(1-14)/12)/4 +  &
                   367*(1-2-((1-14)/12)*12)/12 - 3*((Year+4900+(1-14)/12)/100)/4 + &
                   Jday - Epoch_Date - 1

!---COMPUTES THE GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY)
!   GIVEN THE JULIAN DATE (JD).
L0= Days_Since_Epoch+Epoch_Date+68569
N0= 4*L0/146097
L0= L0-(146097*N0+3)/4
I0= 4000*(L0+1)/1461001
L0= L0-1461*I0/4+31
J0= 80*L0/2447
K0= L0-2447*J0/80       ! DAY
L0= J0/11               
J0= J0+2-12*L0          ! MONTH
I0= 100*(N0-49)+I0+L0   ! YEAR

! Timestamp string format in accordance with the ISO-8601 standard.
write (Year_String,  '(I4.4)') Year
write (Month_String,  '(I2.2)') J0
write (Day_String,  '(I2.2)') K0
Time_Coverage_Start = Year_String//"-"//Month_String//"-"//Day_String&
                //"T"//"00:00:00Z"
Time_Coverage_End = Year_String//"-"//Month_String//"-"//Day_String&
                //"T"//"23:59:59Z"
! CACM 1968 11(10):657, LETTER TO THE EDITOR BY HENRY F. FLIEGEL AND THOMAS C. VAN FLANDERN.

!--- read in file names
do ifile = 1, N_Files_Max
    read(unit=Comp_Asc_Des_Input_Lun,fmt=*,iostat=ios) File_Input(ifile)
    if (ios /= 0) then
       exit
    endif
end do
N_Files = ifile - 1
if (N_Files == 0) then
   print *, LOCAL_EXE_PROMPT, "No files to process, stopping"
   stop
endif

!------------------------------------------------------------------------
! determine if this grid spans the dateline
!------------------------------------------------------------------------
Dateline_Flag = sym%NO
if (Lon_West > 0 .and. Lon_East < 0.0) Dateline_Flag = sym%YES

!------------------------------------------------------------------------
! compute map parameters based on input
!------------------------------------------------------------------------
Ntime_Output = 1
if (Dateline_Flag == sym%NO) then
  Nlon_Output = nint((Lon_East - Lon_West) / Dlon_Output) + 1
  Lon_East = Lon_West + dlon_Output * (Nlon_Output-1)
else
  Nlon_Output = nint( ((360.0+Lon_East) - Lon_West) / dlon_Output) + 1
  Lon_East = Lon_West + dlon_Output * (Nlon_Output-1) - 360.0
endif

Nlat_Output = nint((Lat_North - Lat_South) / Dlat_Output) + 1
Lat_North = Lat_South + dlat_Output * (Nlat_Output-1)

if (Nlat_Output <= 0 .or. Nlon_Output <= 0) then
   print *, LOCAL_EXE_PROMPT, "Error in map parameters, stopping"
   stop
endif

!--- make string needed for output filename
!write (Year_String,  '(I4.4)') Year
write (Jday_string,  '(I3.3)') Jday
write (Time_String,  '(I4.4)') time
write (time_win_string,  '(I2.2)') time_window

if (Node_String == 'zen') then

hour = real(time/100)
minute = real(time - hour*100.0)
Start_Time_Window = hour + (minute-time_window)/60.0
if (Start_Time_Window < 0) then
   Start_Time_Window = 24.0 + Start_Time_Window
   Start_Date_Window = Jday - 1 + Start_Time_Window / 24.0
else
   Start_Date_Window = Jday + Start_Time_Window / 24.0
endif

End_Time_Window = hour + (minute+time_window)/60.0
if (End_Time_Window >= 24.0) then
   End_Time_Window = End_Time_Window - 24.0
   End_Date_Window = Jday + 1 + End_Time_Window / 24.0
else
   End_Date_Window = Jday + End_Time_Window / 24.0
endif

else

    Start_Date_Window = Jday
    End_Date_Window = Jday + 0.999999

endif

if (Node_String == "geo") then
  root_name = trim(File_Input(1))
  ipos = index(root_name,'level2.hdf')
  Time_String = root_name(ipos-5:ipos-1)
endif
if (Node_String == "zen") then
  Sensor_Name_Output = "zen"
endif

call date_and_time(ctime(1), ctime(2), ctime(3))
Date_Created_String = ctime(1)(1:4)//"-"//ctime(1)(5:6)//"-"//ctime(1)(7:8) &
                      //"T"//ctime(2)(1:2)//":"//ctime(2)(3:4)//":"//ctime(2)(5:6)//"Z"

Id_String = "patmosx_"//trim(Product_Version_String)//"_"//trim(Sensor_Name_Output)//"_"// &
            trim(Node_String)//"_d"//trim(Year_String)//trim(Month_String)//trim(Day_String)//"_c"// &
            trim(ctime(1)(1:8))//".nc"

call COMPUTE_WMO_ID_KNOWING_SENSOR_NAME(Sensor_Name_Output,Sc_Id_Output)

!--- store output dimensions
Sds_Dims_Output_xy = (/Nlon_Output,Nlat_Output/)
Sds_Dims_Output_x = (/Nlon_Output/)
Sds_Dims_Output_y = (/Nlat_Output/)
Sds_Output_Start_X = (/0/)
Sds_Output_Stride_X = (/1/)
Sds_Output_Start_Y = (/0/)
Sds_Output_Stride_Y = (/1/)
Sds_Output_Start_XY = (/0,0/)
Sds_Output_Stride_XY = (/1,1/)

    !--- set node integer value
    if (Node_String == "asc") then
     Asc_Des_Node = 0
    endif
    if (Node_String == "des") then
     Asc_Des_Node = 1
    endif
    if (Node_String == "geo") then
     Asc_Des_Node = 1
    endif

    !--------------------------------------------------------------------
    ! create name for composite
    !--------------------------------------------------------------------

    !--- construct output file name for this node
    if (Node_String == "geo") then
     File_Output = "patmosx_"//trim(Sensor_Name_Output)//"_"//trim(Node_String)//"_"// &
                   trim(Year_String)//"_"//trim(Jday_string)//"_"//&
                   trim(Time_String)//".level2b.hdf"
    elseif (Node_String == "zen") then
     File_Output = "patmosx_"//trim(Node_String)//"_"// &
                   trim(Year_String)//"_"//trim(Jday_string)//"_"//&
                   trim(Time_String)//".level2b.hdf"
    else
     File_Output = "patmosx_"//trim(Sensor_Name_Output)//"_"//trim(Node_String)//"_"// &
                   trim(Year_String)//"_"//trim(Jday_string)// &
                   ".level2b.hdf"
    endif

    print *, LOCAL_EXE_PROMPT, "Creating ", trim(File_Output)

    !--- open file for reading and writing
    Sd_Id_Output = sfstart(trim(dir_out)//trim(File_Output),DFACC_CREATE)

    if (Sd_Id_Output <  0) then
       print *, LOCAL_EXE_PROMPT, "output HDF file creation failed. Exiting...", Asc_Des_Node
       stop 38
    endif

    !--------------------------------------------------------------------------------
    ! allocate and initialize arrays for composite
    !--------------------------------------------------------------------------------
     allocate(Time_Output(Nlon_Output,Nlat_Output))
     allocate(Scan_Line_Number_Output(Nlon_Output,Nlat_Output))
     allocate(Scan_Element_Number_Output(Nlon_Output,Nlat_Output))
     allocate(Senzen_Output(Nlon_Output,Nlat_Output))
     allocate(Solzen_Output(Nlon_Output,Nlat_Output))
     allocate(Bad_Pixel_Mask_Output(Nlon_Output,Nlat_Output))
     allocate(Scaled_Sds_Data_Output(Nlon_Output,Nlat_Output))
     allocate(Unscaled_Sds_Data_Output(Nlon_Output,Nlat_Output))
     allocate(Lon_Output(Nlon_Output,Nlat_Output))
     allocate(Lat_Output(Nlon_Output,Nlat_Output))
     allocate(Ielem_Output(Nlon_Output,Nlat_Output))
     allocate(Iline_Output(Nlon_Output,Nlat_Output))
     allocate(Update_Output_Mask(Nlon_Output,Nlat_Output))
     allocate(Overlap_Random_Output(Nlon_Output,Nlat_Output))
     allocate(Lon_Output_1d(Nlon_Output))
     allocate(Lat_Output_1d(Nlat_Output))
     allocate(Scaled_Lon_Output_1d(Nlon_Output))
     allocate(Scaled_Lat_Output_1d(Nlat_Output))
     allocate(Element_Number_Output(Nlon_Output,Nlat_Output))
     allocate(Temp_Mask_Output(Nlon_Output,Nlat_Output))

     allocate(Output_Update_Index(Nlon_Output*Nlat_Output))
     allocate(Output_Array_1d(Nlon_Output*Nlat_Output))

     Output_Array_1d = Missing_Value_Real4
     Time_Output = Missing_Value_Real4
     Scan_Line_Number_Output = missing_value_int4
     Scan_Element_Number_Output = missing_value_int4
     Senzen_Output = Missing_Value_Real4
     Solzen_Output = Missing_Value_Real4
     Bad_Pixel_Mask_Output = Missing_Value_Real4
     Scaled_Sds_Data_Output = Missing_Value_Real4
     Unscaled_Sds_Data_Output = Missing_Value_Real4
     Lon_Output = Missing_Value_Real4
     Lat_Output = Missing_Value_Real4
     Ielem_Output = Missing_Value_Int4
     Iline_Output = Missing_Value_Int4
     Update_Output_Mask = sym%NO
     Overlap_Random_Output = Missing_Value_Real4
     Output_Update_Index = 0
     Element_Number_Output = missing_value_int2
     Temp_Mask_Output = 0
     Source_Length_End = 1

     First_Valid_Input = sym%YES

    !---------------------------------------------------------------------------------
    ! loop through orbit files
    !---------------------------------------------------------------------------------
    file_loop: do ifile = 1, n_files

          print *, LOCAL_EXE_PROMPT, "processing file ", ifile, " of ", n_files, " ", trim(File_Input(ifile))

          Istatus_Sum = 0

          !---- determine if this file should be skipped
          Sd_Id_Input = sfstart(trim(dir_in)//trim(File_Input(ifile)), DFACC_READ)

          if (Sd_Id_Input == FAIL) then
             print *, LOCAL_EXE_PROMPT, "Unable to open: ", trim(File_Input(ifile))
             cycle
          endif

          !--- determine number of sds's in these files (assume the same)
          Istatus_Sum = sffinfo(Sd_Id_Input,Num_Sds_Input,Num_Global_Attrs)

          if (Num_Sds_Input <= 1) then
             print *, LOCAL_EXE_PROMPT, "Empty File:  ", trim(File_Input(ifile))
             cycle
          endif

          !--- read attributes of this file
          Istatus_Sum = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"NUMBER_OF_ELEMENTS"),Num_Elements_Input) + Istatus_Sum
          Istatus_Sum = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"NUMBER_OF_SCANS_LEVEL2"),Num_Lines_Input) + Istatus_Sum
          Istatus_Sum = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"START_YEAR"),Start_Year_Input) + Istatus_Sum
          Istatus_Sum = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"END_YEAR"),End_Year_Input) + Istatus_Sum
          Istatus_Sum = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"START_DAY"),Start_Day_Input) + Istatus_Sum
          Istatus_Sum = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"END_DAY"),End_Day_Input) + Istatus_Sum
          Istatus_Sum = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"START_TIME"),Start_Time_Input) + Istatus_Sum
          Istatus_Sum = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"END_TIME"),End_Time_Input) + Istatus_Sum
          Istatus_Sum = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"WMO_SATELLITE_CODE"),Sc_Id_Input) + Istatus_Sum
          Istatus_Sum = sfrnatt(Sd_Id_Input,sffattr(Sd_Id_Input,"L1B"),L1b_Input) + Istatus_Sum


          Start_Date_Input = Start_Day_Input + Start_Time_Input / 24.0
          End_Date_Input = End_Day_Input + End_Time_Input / 24.0

          !--- based on attributes, see if this file should analyzed further
          if ((Start_Year_Input /= Year) .and. (End_Year_Input /= Year)) then
             print *, "level-2 file year outside bounds for level-2b, skipping file ", &
             Start_Year_Input, End_Year_Input, Year
             cycle
          endif
          if ((Node_String /= "zen") .and. (Start_Day_Input /= Jday) .and. (End_Day_Input /= Jday)) then
             print *, "for nodes other than zen, level-2 file doy outside bounds for level-2b, skipping file ",  &
                      Start_Day_Input, End_Day_Input, Jday
             cycle
          endif
          if ((Node_String /= "zen") .and. (Sc_Id_Output /= Sc_Id_Input)) then
             print *, "for nodes other than zen, level-2 Sat Id outside bounds for level-2b, skipping file ", &
                      Sc_Id_Input, Sc_Id_Output
             cycle
          endif
          if ((Node_String == "zen") .and. (Start_Date_Input < Start_Date_Window)) then
             print *, "for zen node, level-2 start date outside bounds for level-2b, skipping file ", &
                      Start_Date_Input, Start_Date_Window
             cycle
          endif
          if ((Node_String == "zen") .and. (End_Date_Input > End_Date_Window)) then
             print *, "for zen node, level-2 end date outside bounds for level-2b, skipping file ",  &
                      End_Date_Input, End_Date_Window
             cycle
          endif

          !--- pass along Level1b file for a source global variable
          L1b_Length = index(L1b_Input,achar(0)) - 1
          Source_Length_Start = Source_Length_End
          Source_Length_End = Source_Length_Start+L1b_Length+2
          if (Source_Length_End <= Source_Length_Max) then 
           Source_String(Source_Length_Start:Source_Length_End) = L1b_Input(1:L1b_Length)//"; "
          endif

          !--- copy global attributes of first valid file to output file
          if (First_Valid_Input == sym%YES) then

             call COPY_GLOBAL_ATTRIBUTES(Sd_Id_Input,Sd_Id_Output)

             !--- add in attributes that define the output grid
             Istatus_Sum = sfsnatt(Sd_Id_Output,"geospatial_lon_number",DFNT_INT32,1,Nlon_Output) + Istatus_Sum
             Istatus_Sum = sfsnatt(Sd_Id_Output,"geospatial_lat_number",DFNT_INT32,1,Nlat_Output) + Istatus_Sum
             Istatus_Sum = sfsnatt(Sd_Id_Output,"geospatial_lon_spacing",DFNT_FLOAT32,1,dlon_Output) + Istatus_Sum
             Istatus_Sum = sfsnatt(Sd_Id_Output,"geospatial_lat_spacing",DFNT_FLOAT32,1,dlat_Output) + Istatus_Sum
             Istatus_Sum = sfsnatt(Sd_Id_Output,"geospatial_lat_units",DFNT_CHAR8,len_trim("degrees_north"),"degrees_north") + Istatus_Sum
             Istatus_Sum = sfsnatt(Sd_Id_Output,"geospatial_lon_units",DFNT_CHAR8,len_trim("degrees_east"),"degrees_east") + Istatus_Sum

             !------------------------------------------------------------------
             ! add CF-compliant global attributes required for NCDC delivery
             !------------------------------------------------------------------
             Istatus_Sum = sfsnatt(Sd_Id_Output,"geospatial_lon_min",DFNT_FLOAT32,1,Lon_West) + Istatus_Sum
             Istatus_Sum = sfsnatt(Sd_Id_Output,"geospatial_lon_max",DFNT_FLOAT32,1,Lon_East) + Istatus_Sum
             Istatus_Sum = sfsnatt(Sd_Id_Output,"geospatial_lat_max",DFNT_FLOAT32,1,Lat_North) + Istatus_Sum
             Istatus_Sum = sfsnatt(Sd_Id_Output,"geospatial_lat_min",DFNT_FLOAT32,1,Lat_South) + Istatus_Sum
             Istatus_Sum = sfscatt(Sd_Id_Output,"cdm_data_type", DFNT_CHAR8, len_trim("Grid"),trim("Grid")) + Istatus_Sum
             Istatus_Sum = sfsnatt(Sd_Id_Output,"time_coverage_start",DFNT_CHAR8, &
                           len_trim(Time_Coverage_Start),trim(Time_Coverage_Start)) + Istatus_Sum
             Istatus_Sum = sfsnatt(Sd_Id_Output,"time_coverage_end",DFNT_CHAR8, &
                           len_trim(Time_Coverage_End),Time_Coverage_End) + Istatus_Sum
             Istatus_Sum = sfscatt(Sd_Id_Output,"date_created", DFNT_CHAR8, &
                           len_trim(Date_Created_String), trim(Date_Created_String)) + Istatus_Sum
             Istatus_Sum = sfscatt(Sd_Id_Output,"id", DFNT_CHAR8, &
                           len_trim(Id_String),trim(Id_String)) + Istatus_Sum

             
             Istatus_Sum = sfscatt(Sd_Id_Output,"DATA_NODE",DFNT_CHAR8, &
                           len_trim(Node_String),trim(Node_String)) + Istatus_Sum

             if (Node_String /= "zen") then
                Istatus_Sum = sfsnatt(Sd_Id_Output,"COMPOSITE_TIME_HOURS",DFNT_CHAR8, &
                           len_trim(Time_String),Time_String) + Istatus_Sum
                Istatus_Sum = sfscatt(Sd_Id_Output,"COMPOSITE_TIME_WINDOW_HOURS",DFNT_CHAR8, &
                           len_trim(Time_Win_String),Time_Win_String) + Istatus_Sum
             endif

             if (Overlap_Option_String /= "nadir_overlap") then
                Istatus_Sum = sfscatt(Sd_Id_Output,"SENSOR_ZENITH_SAMPLING_METHOD",DFNT_CHAR8, &
                              len_trim("RANDOM"),"RANDOM") + Istatus_Sum
             endif

             if (Overlap_Option_String /= "random_overlap") then
                Istatus_Sum = sfscatt(Sd_Id_Output,"SENSOR_ZENITH_SAMPLING_METHOD",DFNT_CHAR8, &
                              len_trim("MOST_NADIR"),"MOST_NADIR") + Istatus_Sum
             endif


             if (Spatial_Option_String /= "near") then
                Istatus_Sum = sfscatt(Sd_Id_Output,"SPATIAL_SAMPLING_METHOD",DFNT_CHAR8, &
                              len_trim("RANDOM"),"RANDOM") + Istatus_Sum
             endif

             if (Spatial_Option_String /= "rand") then
                Istatus_Sum = sfscatt(Sd_Id_Output,"SPATIAL_SAMPLING_METHOD",DFNT_CHAR8, &
                              len_trim("NEAREST_NEIGHBOR"),"NEAREST_NEIGHBOR") + Istatus_Sum
             endif

          endif

          !--- allocate memory for input data
          allocate(Scaled_Sds_Data_Input(Num_Elements_Input,Num_Lines_Input))
          allocate(Lon_Input(Num_Elements_Input,Num_Lines_Input))
          allocate(Lat_Input(Num_Elements_Input,Num_Lines_Input))
          allocate(Senzen_Input(Num_Elements_Input,Num_Lines_Input))
          allocate(Solzen_Input(Num_Elements_Input,Num_Lines_Input))
          allocate(Gap_Pixel_Mask_Input(Num_Elements_Input,Num_Lines_Input))

          Sds_lon%Variable_Name = "longitude"
          call READ_SDS(Sd_Id_Input,Scaled_Sds_Data_Input,Sds_lon)
          call UNSCALE_SDS(Sds_lon,Scaled_Sds_Data_Input, Lon_Input)

          Sds_Lat%Variable_Name = "latitude"
          call READ_SDS(Sd_Id_Input,Scaled_Sds_Data_Input,Sds_Lat)
          call UNSCALE_SDS(Sds_Lat,Scaled_Sds_Data_Input, Lat_Input)

          !--- read gap pixel mask (goes does not have one)
          Sds_Temp%Variable_Name = "gap_pixel_mask"
          Sds_Temp%Rank = 3
          call READ_SDS(Sd_Id_Input,Gap_Pixel_Mask_Input,Sds_Temp)
          Start_Time_point_Hours = COMPUTE_TIME_HOURS()

          !--- compute remapping arrays for this orbit
          call REGRID( Lat_South,        &
                       Dlat_Output,      &
                       Nlat_Output,      &
                       Lon_West,         &
                       Dlon_Output,      &
                       Nlon_Output,      &
                       Spatial_Flag)

         End_Time_point_Hours = COMPUTE_TIME_HOURS()

         !--- update time summation for level-1b processing
         Segment_Time_Point_Seconds(1) =  Segment_Time_Point_Seconds(1) + &
              60.0*60.0*(End_Time_point_Hours - Start_Time_Point_Hours)

             
         !--- if first valid file, write output lon and lat vectors to output file
         if (First_Valid_Input == sym%YES) then

              !--- longitude
              Sds_Lon_Out = Sds_Lon
              if (Geo_2d_Flag == 0) then    !first row only

                Sds_Lon_Out%Rank = 1
                Lon_Output_1d = Lon_Output(:,1)
                call CORRECT_SDS_STRINGS (Sds_Lon_Out) 
                call DEFINE_SDS_RANK1(Sd_Id_Output, Sds_Dims_Output_X,  &
                                Sds_Dims_Output_X,Sds_Lon_Out)
                call SCALE_SDS(Sds_Lon_Out,Lon_Output_1d,Scaled_Lon_Output_1d)
                call WRITE_SDS(Scaled_Lon_Output_1d,Sds_Lon_Out)
                Istatus_Sum = sfendacc(Sds_Lon_Out%Id_Output) + Istatus_Sum

              else

                Sds_Lon_Out%Rank = 2
                call CORRECT_SDS_STRINGS(Sds_Lon_Out) 
                call DEFINE_SDS_RANK2(Sd_Id_Output, Sds_Dims_Output_XY,  &
                                Sds_Dims_Output_XY,Sds_Lon_Out)
                call SCALE_SDS(Sds_Lon_Out,Lon_Output,Scaled_Sds_Data_Output)
                call WRITE_SDS(Scaled_Sds_Data_Output,Sds_Lon_Out)
                Istatus_Sum = sfendacc(Sds_Lon_Out%Id_Output) + Istatus_Sum

              endif

              !--- latitude 
              Sds_Lat_Out = Sds_Lat
              if (Geo_2d_Flag == 0) then   !first column only

                Sds_Lat_Out%Rank = 1
                Lat_Output_1d = Lat_Output(1,:)
                call CORRECT_SDS_STRINGS (Sds_Lat_Out) 
                call DEFINE_SDS_RANK1(Sd_Id_Output, Sds_Dims_Output_y,  &
                                Sds_Dims_Output_y,Sds_Lat_Out)
                call SCALE_SDS(Sds_Lat_Out,Lat_Output_1d,Scaled_Lat_Output_1d)
                call WRITE_SDS(Scaled_Lat_Output_1d,Sds_Lat_Out)
                Istatus_Sum = sfendacc(Sds_Lat_Out%Id_Output) + Istatus_Sum

              else

                Sds_Lat_Out%Rank = 2
                call CORRECT_SDS_STRINGS (Sds_Lat_Out) 
                call DEFINE_SDS_RANK2(Sd_Id_Output, Sds_Dims_Output_XY,  &
                                Sds_Dims_Output_XY,Sds_Lat_Out)
                call SCALE_SDS(Sds_Lat_Out,Lat_Output,Scaled_Sds_Data_Output)
                call WRITE_SDS(Scaled_Sds_Data_Output,Sds_Lat_Out)
                Istatus_Sum = sfendacc(Sds_Lat_Out%Id_Output) + Istatus_Sum

              endif

              !--- add extra attributes
              Temp_Name = "Y"
              Istatus_Sum = sfsnatt(Sds_Lat_Out%Id_Output, "axis", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum
              Temp_Name = "X"
              Istatus_Sum = sfsnatt(Sds_Lon_Out%Id_Output, "axis", DFNT_CHAR8, len_trim(Temp_Name), trim(Temp_Name)) + Istatus_Sum

         endif

         Num_Points = Num_Elements_Input * Num_Lines_Input
         allocate(Input_Update_Index(Num_Points))
         allocate(Input_Array_1d(Num_Points))

         !--- allocate memory for output on First_Valid_Input
         if (First_Valid_Input == sym%YES) then
            allocate(Scaled_Sds_Data_Output_Full(Nlon_Output, Nlat_Output, Num_Sds_Input))
         endif

         !--------------------------------------------------------------
         !--- bad pixel mask (assumed no scaling applied)
         !--------------------------------------------------------------
         Sds_Temp%Variable_Name = "bad_pixel_mask"
         allocate(Bad_Pixel_Mask_Input(Num_Elements_Input,Num_Lines_Input))
         call READ_SDS(Sd_Id_Input,Bad_Pixel_Mask_Input,Sds_Temp)

         !--------------------------------------------------------------
         !--- ascending / descending flag (assumed no scaling applied)
         !--------------------------------------------------------------
         allocate(Asc_Des_Flag_Input(Num_Lines_Input))
         Sds_Asc_Des_Flag%Rank = 1
         Sds_Asc_Des_Flag%Variable_Name = "asc_des_flag"
         if (Node_String /= "geo") then
             call READ_SDS(Sd_Id_Input,Asc_Des_Flag_Input,Sds_Asc_Des_Flag)
         else
             Asc_Des_Flag_Input(:) = -1 
         endif

         !--------------------------------------------------------------
         !--- read time (assumed no scaling applied)
         !--------------------------------------------------------------
         Sds_Time%Variable_Name = "scan_line_time"
         allocate(time_1d_Input(Num_Lines_Input))
         call READ_SDS(Sd_Id_Input,time_1d_Input,Sds_Time)
         Sds_Time%Rank = 3

         !---- extrapolate scan-line time to each pixel
         allocate(Time_Input(Num_Elements_Input,Num_Lines_Input))
         allocate(Date_Input(Num_Elements_Input,Num_Lines_Input))
         do Iline = 1, Num_Lines_Input
              Time_Input(:,Iline) = Time_1d_Input(Iline)
         enddo
         deallocate(Time_1d_Input)

         !--- convert to a fractional day
         ! Note, this needs to change if we fix geostationary time attributes
         !---
         Date_Input = Start_Day_Input + Time_Input / 24.0
!        if (Start_Day_Input /= End_Day_Input) then 
         if (Node_String == 'zen') then 
            where(Time_Input < Start_Time_Input)
               Date_Input = Start_Day_Input + 1 + Time_Input/24.0
            endwhere
         endif

         !--- define output or read in output values
         if (First_Valid_Input == sym%YES) then   !define sds in output
                 call CORRECT_SDS_STRINGS (Sds_Time) 
                 call DEFINE_SDS_RANK2(Sd_Id_Output, &
                                 Sds_dims_Output_xy, &
                                 Sds_dims_Output_xy, &
                                 Sds_Time)
                 Time_Output(:,:) = Missing_Value_Real4
         endif

         !--------------------------------------------------------------
         !--- scan line number
         !--------------------------------------------------------------
         Sds_scan_line_number%Variable_Name = "scan_line_number"
         allocate(Scan_Line_Number_1d_Input(Num_Lines_Input))

         call READ_SDS(Sd_Id_Input,Scan_Line_Number_1d_Input,Sds_scan_line_number)

         !-- modify level-2 sds parameters for level-2b
         Sds_Scan_Line_Number%Rank = 3
         Sds_Scan_Element_Number%Long_Name = "scan line index of the "// &
                          "pixel chosen for inclusion in level-2b"  

         !---- extrapolate scan-line number to each pixel
         allocate(Scan_Line_Number_Input(Num_Elements_Input,Num_Lines_Input))
         do Iline = 1, Num_Lines_Input
              Scan_Line_Number_Input(:,Iline) = Scan_Line_Number_1d_Input(Iline)
         enddo
         deallocate(Scan_Line_Number_1d_Input)

         !--- define output or read in output values
         if (First_Valid_Input == sym%YES) then   !define sds in output
                 call CORRECT_SDS_STRINGS (Sds_Scan_Line_Number ) 
                 call DEFINE_SDS_RANK2(Sd_Id_Output,Sds_Dims_Output_xy, &
                            Sds_dims_Output_xy,Sds_scan_line_number)
                 Scan_Line_Number_Output = Missing_Value_Int4
         endif

         !--------------------------------------------------------------
         !--- scan element number (this is not read from level-2)
         !--------------------------------------------------------------
         Sds_scan_element_number%Variable_Name = "scan_element_number"
         Sds_scan_element_number%Data_Type = DFNT_INT32
         Sds_scan_element_number%Rank = 3
         Sds_scan_element_number%Scaling_Type = 0
         Sds_scan_element_number%Units = "none"
         Sds_scan_element_number%Standard_Name = "not specified"
         Sds_scan_element_number%Long_Name = "scan element index of the "// &
                          "pixel chosen for inclusion in level-2b"  
         Sds_Scan_Element_Number%Unscaled_Missing = Missing_Value_Real4

         !--- define output or read in output values
         if (First_Valid_Input == sym%YES) then   !define sds in output
             call DEFINE_SDS_RANK2(Sd_Id_Output,Sds_dims_Output_xy, &
                 Sds_dims_Output_xy,Sds_Scan_Element_Number)
                 Scan_Element_Number_Output = missing_value_int4
         endif

         !--------------------------------------------------------------
         !--- sensor zenith
         !--------------------------------------------------------------
         Sds_Temp%Variable_Name = "sensor_zenith_angle"
         call READ_SDS(Sd_Id_Input,Scaled_Sds_Data_Input,Sds_Temp)
         call UNSCALE_SDS(Sds_Temp,Scaled_Sds_Data_Input, Senzen_Input)

         !--------------------------------------------------------------
         !--- solar zenith
         !--------------------------------------------------------------
         Sds_Temp%Variable_Name = "solar_zenith_angle"
         call READ_SDS(Sd_Id_Input,Scaled_Sds_Data_Input,Sds_Temp)
         call UNSCALE_SDS(Sds_Temp,Scaled_Sds_Data_Input,Solzen_Input)
         
         !-------------------------------------------------------------
         !  determine which output grid_points are to be updated by
         !  based on time, date and sensor zenith from the input file
         !-------------------------------------------------------------
         Update_Output_Mask = sym%NO

         Itime = 1

         do Ilon = 1, Nlon_Output
           do Ilat = 1, Nlat_Output

              Temp_Update_Flag = sym%NO

              Ielem = Ielem_Output(Ilon,Ilat)
              Iline = Iline_Output(Ilon,Ilat)

              if ((Ielem > 0) .and. (Iline > 0)) then  !valid Ielem, Iline check

               if (Bad_Pixel_Mask_Input(Ielem,Iline) == sym%NO) then

                !--- check node
                if (Node_String == "zen" .or.  &
                    Node_String == "geo" .or.  &
                    Asc_Des_Flag_Input(Iline) == Asc_Des_Node) then

                 !-- valid data check
                 if (Senzen_Input(Ielem,Iline) /= Missing_Value_Real4) then

                  if (Senzen_Output(Ilon,Ilat) == Missing_Value_Real4)  then
                     Temp_Update_Flag = sym%YES
                  endif
        
                  !--- nadir check
                  if ((Overlap_Flag == 0) .and.  &
                      (Senzen_Output(Ilon,Ilat) /= Missing_Value_Real4) .and. &
                      (Senzen_Output(Ilon,Ilat) - Senzen_Input(Ielem,Iline) > Senzen_Thresh) .and. &
                      (Senzen_Input(Ielem,Iline) < Senzen_Output(Ilon,Ilat))) then
                      Temp_Update_Flag = sym%YES
                  endif

                  !--- random check
                  if (Overlap_Flag == 1) then
                     call random_number(xrand)
                     if ((Overlap_Random_Output(Ilon,Ilat) == Missing_Value_Real4) .or. &
                         (xrand < Overlap_Random_Output(Ilon,Ilat))) then
                        Temp_Update_Flag = sym%YES  
                        Overlap_Random_Output(Ilon,Ilat) = xrand
                     endif
                  endif

                  !--- check time to exclude data from wrong days for orbits that span midnight
                  if (Temp_Update_Flag == sym%YES) then
                      if ((Date_Input(Ielem,Iline) >= Start_Date_Window) .and.  &
                          (Date_Input(Ielem,Iline) <= End_Date_Window)) then
                        Update_Output_Mask(Ilon,Ilat) = sym%YES
                      endif
                  endif

                  !--- update senzen output array
                  if (Update_Output_Mask(Ilon,Ilat) == sym%YES) then
!print *, "updating ", Ilon, Ilat, ielem, iline, Time_Input(ielem, iline)
                    Senzen_Output(Ilon,Ilat) = Senzen_Input(ielem,iline)
                    Time_Output(Ilon,Ilat) = Time_Input(ielem,iline)
                  endif

                endif  !check for valid data
               endif   !check Asc_Des_Node
              endif    !check for bad pixel
             endif    !check for valid Ielem and Iline
            end do
          end do
      
         !----- compute 1d input to output index mapping 
         ipoint = 0
         do Ilon = 1, Nlon_Output
           do Ilat = 1, Nlat_Output
                if (Update_Output_Mask(Ilon,Ilat) == sym%YES) then
                  ipoint = ipoint + 1
                  ipoint = max(0,min(Ipoint,Num_Points))
                  Ielem = Ielem_Output(Ilon,Ilat)
                  Iline = Iline_Output(Ilon,Ilat)
                  Input_Update_Index(Ipoint) = Ielem + (Iline-1)*Num_Elements_Input
                  Output_Update_Index(Ipoint) = Ilon + (Ilat-1)*Nlon_Output
                  Scan_Element_Number_Output(Ilon,Ilat) = Ielem
                  Scan_Line_Number_Output(Ilon,Ilat) = Iline
                  Solzen_Output(Ilon,Ilat) = Solzen_Input(Ielem,Iline)
                  Bad_Pixel_Mask_Output(Ilon,Ilat) = Bad_Pixel_Mask_Input(Ielem,Iline)
                endif
           enddo
          enddo

         deallocate(Senzen_Input)
         deallocate(Solzen_Input)
         deallocate(Time_Input)
         deallocate(Date_Input)
         deallocate(Scan_Line_Number_Input)
         deallocate(Lon_Input)
         deallocate(Lat_Input)
         deallocate(Bad_Pixel_Mask_Input)
         deallocate(Gap_Pixel_Mask_Input)

         !---------------------------------------------------------------------
         ! loop through SDS's and store the information about them
         ! note, this has to be zero based indices
         !---------------------------------------------------------------------
         !--- make a mask for selecting pixels - this will be same for all sds'
         allocate(Mask_Output(Nlon_Output,Nlat_Output))
         Mask_Output = .false.
         where(Update_Output_Mask == sym%YES .and. Ielem_Output > 0 .and. Iline_Output > 0)
            Mask_Output = .true.
         endwhere

         if (First_Valid_Input == sym%YES) then
           Sds_loop_1: do Isds = 0, Num_Sds_Input-1
            
            sds%Id_Input = sfselect(Sd_Id_Input,Isds) 

            Istatus_Sum = sfginfo(sds(Isds)%id_Input,   &
                              sds(Isds)%Variable_Name,  &
                              sds(Isds)%Rank,           &
                              Sds_Dims,                 &
                              sds(Isds)%Data_Type,      &
                              sds(Isds)%Num_Attrs) + Istatus_Sum
           enddo Sds_loop_1
         endif

         !---------------------------------------------------------------------
         ! loop through SDS's and store the information about them
         !---------------------------------------------------------------------
         Sds_loop_2: do Isds = 0, Num_Sds_Input-1

            !--- do not composite lat or lon or other array that are already processed
            if (trim(sds(Isds)%Variable_Name) == trim(Sds_scan_line_number%Variable_Name)) cycle    
            if (trim(sds(Isds)%Variable_Name) == trim(Sds_lon%Variable_Name)) cycle    
            if (trim(sds(Isds)%Variable_Name) == trim(Sds_Lat%Variable_Name)) cycle    
            if (trim(sds(Isds)%Variable_Name) == trim(Sds_Time%Variable_Name)) cycle    
            if (trim(sds(Isds)%Variable_Name) == trim(Sds_Asc_Des_Flag%Variable_Name)) cycle    

            !--- this 3-d packed sds is not supported yet in this code
            if (trim(sds(Isds)%Variable_Name) == 'cloud_mask_test_packed_results') cycle    

            !--- bad_scan_line_flag is 1d and should not be processed
            if (trim(sds(Isds)%Variable_Name) == 'bad_scan_line_flag') cycle    
       
            !-----------------------------------------------------------
            !--- read from input input variable
            !-----------------------------------------------------------
            call READ_SDS(Sd_Id_Input,Scaled_Sds_Data_Input,Sds(Isds))

            !-----------------------------------------------------------
            !--- read from or define  output variable
            !-----------------------------------------------------------
            if (First_Valid_Input == sym%YES) then 
                        
                call CORRECT_SDS_STRINGS (Sds(Isds)) 
                call DEFINE_SDS_RANK2(Sd_Id_Output, Sds_Dims_Output_xy,  &
                                Sds_Dims_Output_xy,Sds(Isds))

                !--- initialize output to missing
                if (sds(Isds)%data_type == DFNT_INT8) then
                    Scaled_Sds_Data_Output = Missing_Value_Int1
                    Scaled_Sds_Data_Output_Full(:,:,Isds) = Missing_Value_Int1
                elseif (sds(Isds)%data_type == DFNT_INT16) then
                    Scaled_Sds_Data_Output = missing_value_int2
                    Scaled_Sds_Data_Output_Full(:,:,Isds) = Missing_Value_Int2
                elseif (sds(Isds)%data_type == DFNT_INT32) then
                    Scaled_Sds_Data_Output = missing_value_int4
                    Scaled_Sds_Data_Output_Full(:,:,Isds) = Missing_Value_Int4
                elseif (sds(Isds)%data_type == DFNT_FLOAT32) then
                    Scaled_Sds_Data_Output = Missing_Value_Real4
                    Scaled_Sds_Data_Output_Full(:,:,Isds) = Missing_Value_Real4
                endif

            endif

            !-----------------------------------------------------------
            !--- update composite
            !--- switching between 2d to 3d arrays necessitated 
            !----by seg faults on ifort
            !-----------------------------------------------------------
            !if (ipoint > 0) then
            ! Input_Array_1d = reshape(Scaled_Sds_Data_Input, (/Num_Elements_Input * Num_Lines_Input/))  
            ! Scaled_Sds_Data_Output = Scaled_Sds_Data_Output_Full(:,:,Isds)
            ! Output_Array_1d = reshape(Scaled_Sds_Data_Output, (/Nlon_Output*Nlat_Output/))  
            ! Output_Array_1d(Output_Update_Index(1:ipoint)) = Input_Array_1d(Input_Update_Index(1:ipoint))
            ! Scaled_Sds_Data_Output = reshape(Output_Array_1d, (/Nlon_Output,Nlat_Output/))
            ! Scaled_Sds_Data_Output_Full(:,:,Isds) = Scaled_Sds_Data_Output
            !endif

            !-----------------------------------------------------------
            !--- update composite
            !-----------------------------------------------------------
            if (ipoint > 0) then
             do Ilon = 1, Nlon_Output
              do Ilat = 1, Nlat_Output
                if (Mask_Output(Ilon,Ilat)) then
                 Scaled_Sds_Data_Output_Full(Ilon,Ilat,Isds) =  &
                           Scaled_Sds_Data_Input(Ielem_Output(Ilon,Ilat),Iline_Output(Ilon,Ilat))
               endif
              enddo
             enddo
            endif

         end do Sds_loop_2
         if (allocated(Mask_Output)) deallocate(Mask_Output)

         !----------------------------------------------------------------------
         ! deallocate allcoated input arrays
         !----------------------------------------------------------------------
         deallocate(Asc_Des_Flag_Input)
         deallocate(Scaled_Sds_Data_Input)
         deallocate(Input_Update_Index)
         deallocate(Input_Array_1d)
         

         !----- close input hdf file
         Istatus_Sum = sfend(Sd_Id_Input) + Istatus_Sum 
      
         First_Valid_Input = sym%NO

    end do file_loop
    !----------------------------------------------------------------------------
    ! after last file, write output if at least one file was successfully read in
    !----------------------------------------------------------------------------

    if (First_Valid_Input == sym%NO) then

            !---- write time to output files
            call SCALE_SDS(Sds_Time,Time_Output,Scaled_Sds_Data_Output)
            call WRITE_SDS(Scaled_Sds_Data_Output,Sds_Time)

            !---- write scan line_number to output files
            call SCALE_SDS(Sds_scan_line_number,Scan_Line_Number_Output,Scaled_Sds_Data_Output)
            call WRITE_SDS(Scaled_Sds_Data_Output,Sds_scan_line_number)

            !---- write scan element number to output files
            call SCALE_SDS(Sds_scan_element_number,Scan_Element_Number_Output,Scaled_Sds_Data_Output)
            call WRITE_SDS(Scaled_Sds_Data_Output,Sds_scan_element_number)

            !---- write level1b source files to global attribute
            Istatus_Sum = sfscatt(Sd_Id_Output,"source", DFNT_CHAR8, &
                          Source_Length_End,Source_String(1:Source_Length_End)) + Istatus_Sum

            Sds_loop_3: do Isds = 0, Num_Sds_Input-1

               !--- do not send these to output
               if (trim(sds(Isds)%Variable_Name) == trim(Sds_scan_line_number%Variable_Name)) cycle    
               if (trim(sds(Isds)%Variable_Name) == trim(Sds_lon%Variable_Name)) cycle    
               if (trim(sds(Isds)%Variable_Name) == trim(Sds_Lat%Variable_Name)) cycle    
               if (trim(sds(Isds)%Variable_Name) == trim(Sds_Time%Variable_Name)) cycle    
               if (trim(sds(Isds)%Variable_Name) == trim(Sds_Asc_Des_Flag%Variable_Name)) cycle    

               !--- write to output
               Scaled_Sds_Data_Output = Scaled_Sds_Data_Output_Full(:,:,Isds)
               call WRITE_SDS(Scaled_Sds_Data_Output,Sds(Isds))

            enddo Sds_loop_3

    else
           print *, LOCAL_EXE_PROMPT, "No data composited"
           !----- close output hdf file
           Istatus_Sum = sfend(Sd_Id_Output) + Istatus_Sum
           call system( 'rm '//trim(dir_out)//trim(File_Output)) 
           print *, LOCAL_EXE_PROMPT, "Level-2b File Deleted, stopping"
           stop
    endif





    !------------------------------------------------------------------
    ! Read certain fields back in and compute required attributes
    !------------------------------------------------------------------

    !--- Domain Mean Ch1 Reflectance
    Sds_Temp%Variable_Name = "refl_0_65um_nom"
    call READ_SDS(Sd_Id_Output,Scaled_Sds_Data_Output,Sds_Temp)
    if ( Sds_Temp%Data_Type >= 0 ) then
       call UNSCALE_SDS(Sds_Temp,Scaled_Sds_Data_Output, Unscaled_Sds_Data_Output)
       Temp_Mask_Output = 0
       where (Solzen_Output /= Missing_Value_Real4 .and. &
           Solzen_Output < 80.0 .and. &
           Unscaled_Sds_Data_Output /= Missing_Value_Real4)
         Temp_Mask_Output = 1
       end where
       Count_Valid = sum(Temp_Mask_Output)
       Domain_Mean_Day_Reflectance_0_65um_nom = Missing_Value_Real4
       if (Count_Valid > 0) then
          Domain_Mean_Day_Reflectance_0_65um_nom =  &
            sum(Temp_Mask_Output*Unscaled_Sds_Data_Output)/Count_Valid
       endif
    endif

    !--- Domain Mean Ch2 Reflectance (goes does not have)
    Sds_Temp%Variable_Name = "refl_0_86um_nom"
    call READ_SDS(Sd_Id_Output,Scaled_Sds_Data_Output,Sds_Temp)
    if ( Sds_Temp%Data_Type >= 0 ) then
       call UNSCALE_SDS(Sds_Temp,Scaled_Sds_Data_Output, Unscaled_Sds_Data_Output)
       Temp_Mask_Output = 0
       where (Solzen_Output /= Missing_Value_Real4 .and. &
              Solzen_Output < 80.0 .and. &
              Unscaled_Sds_Data_Output /= Missing_Value_Real4)
           Temp_Mask_Output = 1
       end where
       Count_Valid = sum(Temp_Mask_Output)
       Domain_Mean_Day_Reflectance_0_86um_nom = Missing_Value_Real4
       if (Count_Valid > 0) then
         Domain_Mean_Day_Reflectance_0_86um_nom =  &
               sum(Temp_Mask_Output*Unscaled_Sds_Data_Output)/Count_Valid
       endif
    endif

    !--- Domain Mean Ch6 Reflectance (goes does not have)
    Sds_Temp%Variable_Name = "refl_1_60um_nom"
    call READ_SDS(Sd_Id_Output,Scaled_Sds_Data_Output,Sds_Temp)
    if ( Sds_Temp%Data_Type >= 0 ) then
       call UNSCALE_SDS(Sds_Temp,Scaled_Sds_Data_Output, Unscaled_Sds_Data_Output)
       Temp_Mask_Output = 0
       where (Solzen_Output /= Missing_Value_Real4 .and. &
            Solzen_Output < 80.0 .and. &
            Unscaled_Sds_Data_Output /= Missing_Value_Real4)
       Temp_Mask_Output = 1
       end where
       Count_Valid = sum(Temp_Mask_Output)
       Domain_Mean_Day_Reflectance_1_60um_nom = Missing_Value_Real4
       if (Count_Valid > 0) then
          Domain_Mean_Day_Reflectance_1_60um_nom =  &
               sum(Temp_Mask_Output*Unscaled_Sds_Data_Output)/Count_Valid
       endif
    endif

    !--- Domain Mean Ch20 Temperature
    Sds_Temp%Variable_Name = "temp_3_75um_nom"
    call READ_SDS(Sd_Id_Output,Scaled_Sds_Data_Output,Sds_Temp)
    call UNSCALE_SDS(Sds_Temp,Scaled_Sds_Data_Output, Unscaled_Sds_Data_Output)
    Temp_Mask_Output = 0
    where (Unscaled_Sds_Data_Output /= Missing_Value_Real4)
        Temp_Mask_Output = 1
    end where
    Count_Valid = sum(Temp_Mask_Output)
    Domain_Mean_Temperature_3_75um_nom = Missing_Value_Real4
    if (Count_Valid > 0) then
       Domain_Mean_Temperature_3_75um_nom =  &
            sum(Temp_Mask_Output*Unscaled_Sds_Data_Output)/Count_Valid
    endif

    !--- Domain Mean Ch31 Temperature
    Sds_Temp%Variable_Name = "temp_11_0um_nom"
    call READ_SDS(Sd_Id_Output,Scaled_Sds_Data_Output,Sds_Temp)
    call UNSCALE_SDS(Sds_Temp,Scaled_Sds_Data_Output, Unscaled_Sds_Data_Output)
    Temp_Mask_Output = 0
    where (Unscaled_Sds_Data_Output /= Missing_Value_Real4)
        Temp_Mask_Output = 1
    end where
    Count_Valid = sum(Temp_Mask_Output)
    Domain_Mean_Temperature_11_0um_nom = Missing_Value_Real4
    if (Count_Valid > 0) then
       Domain_Mean_Temperature_11_0um_nom =  &
            sum(Temp_Mask_Output*Unscaled_Sds_Data_Output)/Count_Valid
    endif
 
    !--- Domain Mean Ch32 Temperature
    Sds_Temp%Variable_Name = "temp_12_0um_nom"
    call READ_SDS(Sd_Id_Output,Scaled_Sds_Data_Output,Sds_Temp)
    call UNSCALE_SDS(Sds_Temp,Scaled_Sds_Data_Output, Unscaled_Sds_Data_Output)
    Temp_Mask_Output = 0
    where (Unscaled_Sds_Data_Output /= Missing_Value_Real4)
        Temp_Mask_Output = 1
    end where
    Count_Valid = sum(Temp_Mask_Output)
    Domain_Mean_Temperature_12_0um_nom = Missing_Value_Real4
    if (Count_Valid > 0) then
       Domain_Mean_Temperature_12_0um_nom =  &
            sum(Temp_Mask_Output*Unscaled_Sds_Data_Output)/Count_Valid
    endif

    !--- Domain Mean Cloud Amount
    Sds_Temp%Variable_Name = "cloud_probability"
    call READ_SDS(Sd_Id_Output,Scaled_Sds_Data_Output,Sds_Temp)
    call UNSCALE_SDS(Sds_Temp,Scaled_Sds_Data_Output, Unscaled_Sds_Data_Output)
    Temp_Mask_Output = 0
    where (Unscaled_Sds_Data_Output /= Missing_Value_Real4)
        Temp_Mask_Output = 1
    end where
    Count_Valid = sum(Temp_Mask_Output)
    Count_Cloud = count(Unscaled_Sds_Data_Output >= 0.5)
    Domain_Cloud_Fraction = Missing_Value_Real4
    if (Count_Valid > 0) then
       Domain_Cloud_Fraction = real(Count_Cloud) / real(Count_Valid)
    endif

    !--- Ice Cloud
    Sds_Temp%Variable_Name = "cloud_type"
    call READ_SDS(Sd_Id_Output,Scaled_Sds_Data_Output,Sds_Temp)
    call UNSCALE_SDS(Sds_Temp,Scaled_Sds_Data_Output, Unscaled_Sds_Data_Output)
    Temp_Mask_Output = 0
    where (Unscaled_Sds_Data_Output /= real(Missing_Value_Int1))
        Temp_Mask_Output = 1
    end where
    Count_Valid = sum(Temp_Mask_Output)
    Count_Cloud = count(Unscaled_Sds_Data_Output >= 6  .and. Unscaled_Sds_Data_Output <= 9)
    Domain_Ice_Cloud_Fraction = Missing_Value_Real4
    if (Count_Valid > 0) then
       Domain_Ice_Cloud_Fraction = real(Count_Cloud) / real(Count_Valid)
    endif

    !--- ACHA Success
    Domain_Mean_Acha_Success_Fraction = Missing_Value_Real4
    Sds_Temp%Variable_Name = "acha_quality"
    call READ_SDS(Sd_Id_Output,Scaled_Sds_Data_Output,Sds_Temp)
    if ( Sds_Temp%Data_Type >= 0 ) then
       !--- tried count
       Count_Total = count(btest(int(Scaled_Sds_Data_Output),0))

       !--- success count
       Count_Valid = count((btest(int(Scaled_Sds_Data_Output),1)) .and.  &
                           (btest(int(Scaled_Sds_Data_Output),2)) .and. &
                           (btest(int(Scaled_Sds_Data_Output),3)))

       if (Count_Total > 0) then
          Domain_Mean_Acha_Success_Fraction = float(Count_Valid)/float(Count_Total)
       endif

    endif

    !--- DCOMP Success
    Domain_Mean_Dcomp_Success_Fraction = Missing_Value_Real4
    Sds_Temp%Variable_Name = "dcomp_quality"
    call READ_SDS(Sd_Id_Output,Scaled_Sds_Data_Output,Sds_Temp)
    if ( Sds_Temp%Data_Type >= 0 ) then

       !--- tried count
       Count_Total = count(btest(int(Scaled_Sds_Data_Output),0))

       !--- success count
       Count_Valid = count( btest(int(Scaled_Sds_Data_Output),0) .and. &
                            (.not. btest(int(Scaled_Sds_Data_Output),1)) .and.  &
                            (.not. btest(int(Scaled_Sds_Data_Output),2)) )

       if (Count_Total > 0) then
          Domain_Mean_Dcomp_Success_Fraction = float(Count_Valid)/float(Count_Total)
       endif

    endif

    !--- compute some metrics to include as attributes
    Start_Time_point_Hours = COMPUTE_TIME_HOURS()

    Count_Total = count(Bad_Pixel_Mask_Output /= Missing_Value_Real4)

    Domain_Data_Fraction = real(Count_Total) /  &
                           real(Nlon_Output * Nlat_Output)

    Count_Valid = count(Bad_Pixel_Mask_Output == sym%NO)

    Domain_Valid_Data_Fraction = real(Count_Valid) /  &
                           real(Nlon_Output * Nlat_Output)

    Count_Day = count(Solzen_Output /= Missing_Value_Real4 .and. &
                      Solzen_Output < 90.0)

    Domain_Day_Fraction = real(Count_Day) / real(Nlon_Output * Nlat_Output)


    !--- write these attributes
    Istatus_Sum = sfsnatt(Sd_Id_Output,"FRACTION_WITH_DATA", &
                      DFNT_FLOAT32,1,Domain_Data_Fraction) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"FRACTION_WITH_VALID_DATA", &
                      DFNT_FLOAT32,1,Domain_Valid_Data_Fraction) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"FRACTION_WITH_DAYTIME_DATA", &
                      DFNT_FLOAT32,1,Domain_Day_Fraction) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"FRACTION_WITH_CLOUDY_DATA", &
                      DFNT_FLOAT32,1,Domain_Cloud_Fraction) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"FRACTION_WITH_ICE_CLOUDY_DATA", &
                      DFNT_FLOAT32,1,Domain_Ice_Cloud_Fraction) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"MEAN_DAY_REFLECTANCE_0_65UM_NOM", &
                      DFNT_FLOAT32,1,Domain_Mean_Day_Reflectance_0_65um_nom) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"MEAN_DAY_REFLECTANCE_0_86UM_NOM", &
                      DFNT_FLOAT32,1,Domain_Mean_Day_Reflectance_0_86um_nom) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"MEAN_DAY_REFLECTANCE_1_60UM_NOM", &
                      DFNT_FLOAT32,1,Domain_Mean_Day_Reflectance_1_60um_nom) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"MEAN_TEMPERATURE_3_750UM_NOM", &
                      DFNT_FLOAT32,1,Domain_Mean_Temperature_3_75um_nom) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"MEAN_TEMPERATURE_11_0UM_NOM", &
                      DFNT_FLOAT32,1,Domain_Mean_Temperature_11_0um_nom) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"MEAN_TEMPERATURE_12_0UM_NOM", &
                      DFNT_FLOAT32,1,Domain_Mean_Temperature_12_0um_nom) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"ACHA_SUCCESS_FRACTION", &
                      DFNT_FLOAT32,1,Domain_Mean_Acha_Success_Fraction) + Istatus_Sum

    Istatus_Sum = sfsnatt(Sd_Id_Output,"DCOMP_SUCCESS_FRACTION", &
                      DFNT_FLOAT32,1,Domain_Mean_Dcomp_Success_Fraction) + Istatus_Sum

    End_Time_point_Hours = COMPUTE_TIME_HOURS()

    !--- update time summation for level-1b processing
    Segment_Time_Point_Seconds(2) =  Segment_Time_Point_Seconds(2) + &
              60.0*60.0*(End_Time_point_Hours - Start_Time_Point_Hours)

    !--- deallocate output fields
    if (allocated(Time_Output)) deallocate(Time_Output)
    if (allocated(Scan_Line_Number_Output)) deallocate(Scan_Line_Number_Output)
    if (allocated(Scan_Element_Number_Output)) deallocate(Scan_Element_Number_Output)
    if (allocated(Senzen_Output)) deallocate(Senzen_Output)
    if (allocated(Solzen_Output)) deallocate(Solzen_Output)
    if (allocated(Scaled_Sds_Data_Output)) deallocate(Scaled_Sds_Data_Output)
    if (allocated(Unscaled_Sds_Data_Output)) deallocate(Unscaled_Sds_Data_Output)
    if (allocated(Lon_Output)) deallocate(Lon_Output)
    if (allocated(Lat_Output)) deallocate(Lat_Output)
    if (allocated(Ielem_Output)) deallocate(Ielem_Output)
    if (allocated(Iline_Output)) deallocate(Iline_Output)
    if (allocated(Update_Output_Mask)) deallocate(Update_Output_Mask)
    if (allocated(Overlap_Random_Output)) deallocate(Overlap_Random_Output)
    if (allocated(Output_Update_Index)) deallocate(Output_Update_Index)
    if (allocated(Output_Array_1d)) deallocate(Output_Array_1d)
    if (allocated(Lon_Output_1d)) deallocate(Lon_Output_1d)
    if (allocated(Lon_Output_1d)) deallocate(Lat_Output_1d)
    if (allocated(Scaled_Lon_Output_1d)) deallocate(Scaled_Lon_Output_1d)
    if (allocated(Scaled_Lat_Output_1d)) deallocate(Scaled_Lat_Output_1d)
    if (allocated(Scaled_Sds_Data_Output_Full)) deallocate(Scaled_Sds_Data_Output_Full)
    if (allocated(Element_Number_Output)) deallocate(Element_Number_Output)
    if (allocated(Bad_Pixel_Mask_Output)) deallocate(Bad_Pixel_Mask_Output)


    !----- close output hdf file
    Istatus_Sum = sfend(Sd_Id_Output) + Istatus_Sum

    !--- notify if accumulated hdf call istatus is non-zero
    if (Istatus_Sum /= 0) then
           print *, LOCAL_EXE_PROMPT, "Non-zero accumulated error codes on hdf calls"
    endif

    !------------------------------------------------------------------------------
    ! open files and compute some metrics
    !------------------------------------------------------------------------------

    !--- open file for reading and writing
    Sd_Id_Output = sfstart(trim(dir_out)//trim(File_Output),DFACC_RDWR)

    !----- close output hdf file
    Istatus_Sum = sfend(Sd_Id_Output) + Istatus_Sum


    Processing_Time_End_Hours = COMPUTE_TIME_HOURS()

    !--- diagnostic screen output
    print *, EXE_PROMPT, "Time for Regridding (sec) = ", Segment_Time_Point_Seconds(1)
    print *, EXE_PROMPT, "Time for Global Attributes (sec) = ", Segment_Time_Point_Seconds(2)
    print *, LOCAL_EXE_PROMPT, 'Total Processing Time (minutes) = ',  &
          60.0*(Processing_Time_End_Hours - Processing_Time_Start_Hours)

!----------------------------------------------------------------------
! End Executabe Code
!----------------------------------------------------------------------
   contains
   
   ! - corrects bad string charaters
   subroutine CORRECT_SDS_STRINGS( sds_data)
        implicit none  
    type(Sds_Struct) , intent(inout) :: Sds_data
   
      sds_data%units = TRIM_C ( sds_data%units, len (sds_data%units) )
      sds_data%long_name = TRIM_C ( sds_data%long_name, len (sds_data%long_name) )
      sds_data%standard_name = TRIM_C ( sds_data%standard_name, len (sds_data%standard_name) )
   end subroutine CORRECT_SDS_STRINGS
   
   ! - removes char(0) characters
   function TRIM_C ( string, length) result (corr_string)
      implicit none
      character ( len = *)  , intent(in) :: string
      integer , intent(in) :: length
      character( len = length)  :: corr_string
      integer :: pos
   
      pos = INDEX(string ,CHAR(0) )
  
      if (pos > 1) corr_string = TRIM(string( 1: pos -1 ) )
      if (pos == 1) corr_string = 'not specified'
      if (pos == 0 ) corr_string = string
   
   end function TRIM_C
   
end program COMPILE_ASC_DES_LEVEL2B
