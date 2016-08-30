!$Id: viirs_nasa_read_module.f90 1541 2016-02-26 23:18:41Z dbotambekov $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version
! 5.4
!
! NAME: viirs_nasa_read_module.f90
!
! PURPOSE: VIIRS NASA (NetCDF4) read tool
!
! DESCRIPTION:  This module deals with reading VIIRS NASA data
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
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
! HISTORY:   created      May 2016 (Denis B.)
!--------------------------------------------------------------------------------------

module VIIRS_NASA_READ_MODULE 

  use HDF5
  use READH5DATASET 
  use FILE_TOOLS, only: &
        FILE_SEARCH &
      , GETLUN
  use PIXEL_COMMON, only: &
        Sensor &
      , Gap_Pixel_Mask
  use CONSTANTS, only: &
        Real4 &
      , Real8 &
      , Int8 &
      , Int4 &
      , Int2 &
      , Int1 &
      , Sym &
      , Exe_Prompt &
      , Missing_Value_Real4 &
      , Missing_Value_Int4 

  implicit none

  public :: READ_VIIRS_NASA_DATA
  public :: READ_VIIRS_NASA_DATE_TIME
  public :: READ_NUMBER_OF_SCANS_VIIRS_NASA
  private :: DETERMINE_VIIRS_NASA_FILE
  private :: JULIAN
  private :: CONVERT_VIIRS_NASA_RADIANCE
  private :: CONVERT_RAD_2_SOL_REF_DNB

  character(len=18), parameter:: VIIRS_NASA_PROMPT="VIIRS_NASA_MODULE:"
 
  ! - bowtie gaps values
  integer, parameter :: Ny_Pattern = 16
  integer, parameter :: Nx_Pattern = 3200
  integer(kind=int4), dimension(:,:), allocatable, public :: Gap_Line_Idx_Pattern
  integer(kind=int1), dimension(:,:), allocatable, public :: Gap_Pixel_Mask_Pattern

      contains

!--------------------------------------------------------------------
! read viirs nasa data
!--------------------------------------------------------------------
subroutine READ_VIIRS_NASA_DATA (Segment_Number, VGEOM_File, Error_Out)

   use Pixel_Common , only : &
        Image &
      , Sensor &
      , Geo &
      , Nav &
      , Scan_Number &
      , Ancil_Data_Dir &
      , Cloud_Mask_Aux_Flag &
      , Cloud_Mask_Aux_Read_Flag &
      , Cld_Mask_Aux &
      , Cld_Type_Aux &
      , Cld_Phase_Aux &
      , Scan_Time &
      , Ch &
      , Ref_ChI1 &
      , Ref_ChI2 &
      , Ref_ChI3 &
      , Bt_ChI4 &
      , Bt_ChI5 &
      , Ref_Min_ChI1 &
      , Ref_Max_ChI1 &
      , Ref_Mean_ChI1 &
      , Ref_Uni_ChI1 &
      , Ref_Min_ChI2 &
      , Ref_Max_ChI2 &
      , Ref_Mean_ChI2 &
      , Ref_Uni_ChI2 &
      , Ref_Min_ChI3 &
      , Ref_Max_ChI3 &
      , Ref_Mean_ChI3 &
      , Ref_Uni_ChI3 &
      , Bt_Min_ChI4 &
      , Bt_Max_ChI4 &
      , Bt_Mean_ChI4 &
      , Bt_Uni_ChI4 &
      , Bt_Min_ChI5 &
      , Bt_Max_ChI5 &
      , Bt_Mean_ChI5 &
      , Bt_Uni_ChI5

   use PLANCK

   use VIEWING_GEOMETRY_MODULE, only: &
        GLINT_ANGLE &
      , SCATTERING_ANGLE &
      , RELATIVE_AZIMUTH

   use CALIBRATION_CONSTANTS, only: &
        Planck_Nu

      integer(kind=int4), intent(in):: Segment_Number
      character(len=*), intent(in):: VGEOM_File
      integer(kind=int4), intent(out):: Error_Out

      integer(kind=int4) :: Nx_Start , Nx_End , Ny_Start , Ny_End , N_Seg_Lines
      integer(kind=int4) :: I_Geo
      integer(kind=int4) :: Mband_Start
      integer(kind=int4) :: I_Mband
      integer(kind=int4) :: k , i
      integer(kind=int4) :: Lun
      integer(kind=int4) :: Io_Err_Stat
      integer(kind=int4), dimension(2) :: Dim_Seg
      integer(kind=int4), dimension(2) :: Dim_Dnb_Seg
      integer(kind=int4), dimension(2) :: Offset_Mband
      integer(kind=int4), dimension(3200) :: D2m_Idx
      integer(kind=int4), dimension(:), allocatable :: Time_Sec_Day
      integer(kind=int4), dimension(7) :: Scaled_Geo = [0, 0, 0, 1, 1, 1, 1]
      integer(kind=int4), dimension(16) :: Modis_Chn_List
      real(kind=real4) :: Fill_Value
      real(kind=real4) :: Scale_Factor , Add_Offset
      real(kind=real8), parameter :: Sec_Per_Day = 86400.
      real(kind=real8), dimension(:), pointer :: R1d_Buffer
      real(kind=real4), dimension(:,:), pointer :: R2d_Buffer
      integer(kind=int4), dimension(:,:), pointer :: I2d_Buffer
      character(len=1020) :: File_2_Read
      character(len=1020) :: File_Dnb_Idx
      character(len=2) :: Band_Num_Str
      character(len=5) :: Day_Night_Flag
      character(len=100), dimension ( 7 ) :: Setname_Geo_List = (/ character (len =300) :: &
                           'geolocation_data/latitude             ' & ! 1
                         , 'geolocation_data/longitude            ' & ! 2
                         , 'scan_line_attributes/scan_start_time  ' & ! 3
                         , 'geolocation_data/sensor_azimuth       ' & ! 4
                         , 'geolocation_data/sensor_zenith        ' & ! 5
                         , 'geolocation_data/solar_azimuth        ' & ! 6
                         , 'geolocation_data/solar_zenith         ' & ! 7
                                                    /)
      character(len=100) :: Setname_Name
      logical, dimension(16) :: Is_Mband_On
      logical, dimension(16) :: Is_Mband_Read

      Error_Out = 0

      ! --- calculate start and end segment limits
      Nx_Start = 1
      Ny_Start = (Segment_Number - 1) * Image%Number_Of_Lines_Per_Segment + 1
      N_Seg_Lines = min (Ny_Start + Image%Number_Of_Lines_Per_Segment -1 , &
                                Image%Number_Of_Lines) - Ny_Start  + 1
      Nx_End = Nx_Start + Image%Number_Of_Elements - 1
      Ny_End = Ny_Start + N_Seg_Lines - 1
      Offset_Mband = [Nx_Start - 1, Ny_Start - 1]
      Dim_Seg = [Nx_End - Nx_Start + 1, Ny_End - Ny_Start + 1]
            
      !print *, trim(Image%Level1b_Path)//trim(VGEOM_File)

      ! --- read geo file data
      ! loop over geo data
      do I_Geo = 1, 7
        if ( I_Geo == 3 ) then
           call H5READDATASET (trim(Image%Level1b_Path)//trim(VGEOM_File), &
                      trim(Setname_Geo_List(I_Geo)), R1d_Buffer)
        else
           call H5READDATASET (trim(Image%Level1b_Path)//trim(VGEOM_File), &
                      trim(Setname_Geo_List(I_Geo)), Offset_Mband, Dim_Seg, R2d_Buffer)
           if (Scaled_Geo(I_Geo) .eq. 1) then
              call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(VGEOM_File), &
                      trim(Setname_Geo_List(I_Geo)//'/scale_factor'), Scale_Factor)
              call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(VGEOM_File), &
                      trim(Setname_Geo_List(I_Geo)//'/add_offset'), Add_Offset)
              call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(VGEOM_File), &
                      trim(Setname_Geo_List(I_Geo)//'/_FillValue'), Fill_Value)
              where ( R2d_Buffer .ne. Fill_Value)
                 R2d_Buffer = (R2d_Buffer * Scale_Factor) + Add_Offset
              endwhere
           endif

        endif

        ! save read data to output in correct format
        select case (I_Geo)
         case(1)
            Nav % Lat_1b (:,1:N_Seg_Lines) = R2d_Buffer
         case(2)
            Nav % Lon_1b (:,1:N_Seg_Lines) = R2d_Buffer
         case(3)
            allocate (Time_Sec_Day (size(R1d_Buffer)))
            Time_Sec_Day = int((mod(R1d_Buffer, Sec_Per_Day)) * 1000)

            ! make data missing of missing scan time
            where (R1d_Buffer < 0)
               Time_Sec_Day = Missing_Value_Int4
            endwhere
            Scan_Time (1:N_Seg_Lines) = (/(Time_Sec_Day((k - 1) / 16 + 1), &
                                            k = Ny_Start, Ny_End)/)
            deallocate ( Time_Sec_Day )
         case(4)
            Geo % Sataz (:,1:N_Seg_Lines) = R2d_Buffer
         case(5)
            Geo % Satzen (:,1:N_Seg_Lines) = R2d_Buffer
         case(6)
            Geo % Solaz (:,1:N_Seg_Lines) = R2d_Buffer
         case(7)
            Geo % Solzen (:,1:N_Seg_Lines) = R2d_Buffer
         end select

         if (I_Geo /= 3)  deallocate ( R2d_Buffer)
         if (I_Geo == 3) deallocate ( R1d_Buffer)
      enddo

      ! --- calculate ascending/descending
      Nav % Ascend (1:N_Seg_Lines) = 0
      do i = 1, N_Seg_Lines - 1
         if (Nav%Lat_1b (Dim_Seg(1) / 2, i + 1) <= Nav%Lat_1b(Dim_Seg(1) / 2, i)) &
                    Nav % Ascend( i ) = 1
      end do
      ! --- fix for the last line
      if ( Nav%Lat_1b(Nx_End / 2, N_Seg_Lines) <= &
           Nav%Lat_1b(Nx_End / 2, N_Seg_Lines - 1) ) &
                Nav % Ascend (N_Seg_Lines) = 1

      ! --- calculate relative azimuth
      Geo % Relaz = RELATIVE_AZIMUTH (Geo%Solaz, Geo%Sataz)

      !--- compute the glint zenith angle
      Geo % Glintzen = GLINT_ANGLE (Geo%Solzen, Geo%Satzen, Geo%Relaz )

      !--- compute the scattering angle
      Geo % Scatangle = SCATTERING_ANGLE (Geo%Solzen, Geo%Satzen, Geo%Relaz)


      ! --- read m-band data
      ! find m-band file
      call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                    'VL1BM',File_2_Read)   
      ! --- if no file quit
      if (trim(File_2_Read) .eq. 'no_file') then
         Error_Out = 1
         return
      endif

      ! --- calculate bowtie pixels
      call COMPUTE_VIIRS_BOWTIE_GAP_PATTERN()

      ! --- read global attribute Day_Night_Flag
      !!!!! Note: if Day_Night_Flag = 'Night' NO M01 - M06, M09, and M11 bands  !!!!!
      !!!!! If Day_Night_Flag = 'Day' or 'Mixed' then ALL M-bands exist         !!!!!
      call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                                       'day_night_flag', Day_Night_Flag)
      Mband_Start = 1
      if (trim(Day_Night_Flag) .eq. 'Night') Mband_Start = 7

      ! --- mapping modis to viirs
      !                 041 044 048 055 068  074  085 124 138  160 225 375  405  855  108  120
      !                 M1  M2  M3  M4  M5   M6   M7  M8  M9   M10 M11 M12  M13  M14  M15  M16
      Modis_Chn_List = [ 8 , 9 , 3 , 4 , 1 , 15 , 2 , 5 , 26 , 6 , 7 , 20 , 22 , 29 , 31 , 32 ]
      Is_Mband_On = Sensor%Chan_On_Flag_Default (Modis_Chn_List) == sym%YES

      ! --- loop over m-band channels
      do I_Mband = Mband_Start , 16

         ! - one more filter for missing night bands
         if (trim(Day_Night_Flag) .eq. 'Night' .and. &
                (I_Mband .eq. 9 .or. I_Mband .eq. 11)) cycle

         ! - check if channel is on 
         Is_Mband_Read(I_Mband) = .false.
         if (.not. Is_Mband_On(I_Mband)) cycle
          
         ! - make string band number
         write (Band_Num_Str, '(i0.2)' )  I_Mband
         Setname_Name = 'observation_data/M'//Band_Num_Str

         ! - read data
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                             Setname_Name, Offset_Mband, Dim_Seg, I2d_Buffer)

         ! - read attribute to unscale
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/scale_factor', Scale_Factor)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/add_offset', Add_Offset)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/_FillValue', Fill_Value)

         ! - separate reflectance and radiance data
         if (I_Mband .le. 11) then
            ! - initialize output
            Ch (Modis_Chn_List(I_Mband)) % Ref_Toa (:, 1:N_Seg_Lines) = Missing_Value_Real4

            ! - unscale
            where ( I2d_Buffer .ne. Fill_Value)
               Ch (Modis_Chn_List(I_Mband)) % Ref_Toa (:, 1:N_Seg_Lines) = &
                                 ((I2d_Buffer * Scale_Factor) + Add_Offset) * 100.
            endwhere

            ! - fill the gaps
            call FILL_VIIRS_BOWTIE_GAPS ( Ny_Start, N_Seg_Lines, Ch (Modis_Chn_List(I_Mband)) % Ref_Toa )

         else
            ! - initialize output
            Ch (Modis_Chn_List(I_Mband)) % Rad_Toa (:, 1:N_Seg_Lines) = Missing_Value_Real4

            ! - unscale
            where ( I2d_Buffer .ne. Fill_Value)
               Ch (Modis_Chn_List(I_Mband)) % Rad_Toa (:, 1:N_Seg_Lines) = &
                                 ((I2d_Buffer * Scale_Factor) + Add_Offset)
            endwhere

            ! - fill the gaps
            call FILL_VIIRS_BOWTIE_GAPS ( Ny_Start, N_Seg_Lines, Ch (Modis_Chn_List(I_Mband)) % Rad_Toa )

            ! - convert radiance to noaa format
            call CONVERT_VIIRS_NASA_RADIANCE (Ch(Modis_Chn_List(I_Mband))%Rad_Toa(:, 1:N_Seg_Lines), &
                                Planck_Nu(Modis_Chn_List(I_Mband)), Missing_Value_Real4) 

            ! - calculate bt
            call COMPUTE_BT_ARRAY (Ch(Modis_Chn_List(I_Mband))%Bt_Toa (:,1:N_Seg_Lines), &
                                   Ch(Modis_Chn_List(I_Mband))%Rad_Toa (:,1:N_Seg_Lines), &
                                   Modis_Chn_List(I_Mband), Missing_Value_Real4)

         endif

         deallocate (I2d_Buffer)
         Is_Mband_Read(I_Mband) = .true.

      enddo  ! - m-bands loop


      ! --- read dnb data (use do loop to quit if no file)
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES) then
      do k = 1, 1

         ! - find dnb geo file
         call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                    'VGEOD',File_2_Read)
         ! - if no file quit
         if (trim(File_2_Read) .eq. 'no_file') then
            print *,'Error: No VGEOD File Found, Skipping'
            cycle
         endif

         ! - mapping file (maps from dnb to M-bands resolution)
         File_Dnb_Idx = trim(Ancil_Data_Dir)//'static/viirs/dnb2m_indx.txt'
         Lun = GETLUN()
         Dim_Dnb_Seg(1) = 4064
         Dim_Dnb_Seg(2) = Dim_Seg(2)

         open (unit = Lun , file=trim ( File_Dnb_Idx) , status="old",action="read")
         read (unit = Lun , fmt=* , iostat = Io_Err_Stat) D2m_Idx
         close (unit = Lun)

         ! - read dnb geo data
         Setname_Name = 'geolocation_data/lunar_azimuth'
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                        Setname_Name, Offset_Mband, Dim_Dnb_Seg, I2d_Buffer)

         ! - read attribute to unscale
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/scale_factor', Scale_Factor)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/add_offset', Add_Offset)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/_FillValue', Fill_Value)

         ! - unscale and remap
         do i = 1, Dim_Seg(1)
            Geo % Lunaz (i, 1:N_Seg_Lines) = (I2d_Buffer(D2m_Idx(i), :) &
                                    * Scale_Factor) + Add_Offset
         end do
         deallocate ( I2d_Buffer )

         Setname_Name = 'geolocation_data/lunar_zenith'
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                        Setname_Name, Offset_Mband, Dim_Dnb_Seg, I2d_Buffer)

         ! - read attribute to unscale
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/scale_factor', Scale_Factor)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/add_offset', Add_Offset)
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                              Setname_Name//'/_FillValue', Fill_Value)

         ! - unscale and remap
         do i = 1, Dim_Seg(1)
            Geo % Lunzen (i, 1:N_Seg_Lines) = (I2d_Buffer(D2m_Idx(i), :) &
                                    * Scale_Factor) + Add_Offset
         end do

         ! - set missing
         where ( I2d_Buffer .eq. Fill_Value)
            Geo % Lunaz (:, 1:N_Seg_Lines) = Missing_Value_Real4
            Geo % Lunzen (:, 1:N_Seg_Lines) = Missing_Value_Real4
         endwhere
         deallocate ( I2d_Buffer )

         ! - read moon phase angle and moon illumination fraction
         call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(File_2_Read), &
                      'moon_phase_angle', Geo % Moon_Phase_Angle)
         call H5READATTRIBUTE (trim(Image%Level1b_Path)//trim(File_2_Read), &                                                                                   
                      'moon_illumination_fraction', Geo % Moon_Illum_Frac)

         ! - compute lunar relative azimuth
         Geo % Lunrelaz (:, 1:N_Seg_Lines) = RELATIVE_AZIMUTH ( &
                   Geo % Lunaz (:, 1:N_Seg_Lines), &
                   Geo % Sataz (:, 1:N_Seg_Lines) )

         ! - compute lunar glint zenith
         Geo % Glintzen_Lunar (:, 1:N_Seg_Lines) = GLINT_ANGLE ( &
                   Geo % Lunzen (:, 1:N_Seg_Lines), &
                   Geo % Satzen (:, 1:N_Seg_Lines), &
                   Geo % Lunrelaz (:, 1:N_Seg_Lines) )

         ! - compute lunar scattering angle
         Geo % Scatangle_Lunar (:, 1:N_Seg_Lines) = SCATTERING_ANGLE ( &
                   Geo % Lunzen (:, 1:N_Seg_Lines), &
                   Geo % Satzen (:, 1:N_Seg_Lines), &
                   Geo % Lunrelaz (:, 1:N_Seg_Lines) )

         ! - find dnb file
         call DETERMINE_VIIRS_NASA_FILE(trim(Image%Level1b_Path),trim(VGEOM_File), &
                                    'VL1BD',File_2_Read)
         ! - if no file quit
         if (trim(File_2_Read) .eq. 'no_file') then
            print *,'Error: No VL1BD File Found, Skipping'
            cycle
         endif

         ! - read dnb data
         Setname_Name = 'observation_data/DNB_observations'
         call H5READDATASET (trim(Image%Level1b_Path)//trim(File_2_Read), &
                        Setname_Name, Offset_Mband, Dim_Dnb_Seg, R2d_Buffer)

         ! - remap
         do i = 1, Dim_Seg(1)
            Ch (44) % Rad_Toa (i, 1:N_Seg_Lines) = R2d_Buffer(D2m_Idx(i), :)
         end do
         deallocate ( R2d_Buffer )

         ! - convert radiance to reflectance
         call CONVERT_RAD_2_SOL_REF_DNB (Ch (44) % Rad_Toa, Geo % Solzen, &
                      Image % Start_Doy, Missing_Value_Real4, Ch (44) % Ref_Toa)

      enddo ! one time loop
      endif ! dnb on


      ! --- global variables which have to be set
      Image%Number_Of_Lines_Read_This_Segment = N_Seg_Lines
      do i = 1, Image%Number_Of_Lines_Per_Segment
         Scan_Number(i) = Ny_Start + i - 1
      end do

end subroutine READ_VIIRS_NASA_DATA

!--------------------------------------------------------------------
! read viirs nasa time
!--------------------------------------------------------------------
subroutine READ_VIIRS_NASA_DATE_TIME (Path, Infile, Year , Doy , Start_Time &
                , End_Time , Orbit , End_Year, End_Doy )

      character(len=*), intent(in) :: Path
      character(len=*), intent(in) :: Infile
      integer, intent(out) :: Year
      integer, intent(out) :: Doy    !day of year
      integer, intent(out) :: Start_Time  !millisec
      integer, intent(out) :: End_Time    !millisec
      integer, intent(out) :: Orbit
      integer , intent(out) :: End_Year
      integer, intent(out) :: End_Doy    !day of year

      character(len=35) :: Tmp_String 
      integer :: Month
      integer :: Day
      integer :: Start_Hour
      integer :: Start_Minute
      integer :: Start_Sec

      integer :: Days_Of_Year
      integer :: End_Hour
      integer :: End_Minute
      integer :: End_Sec


      ! --- read time and date from the attributes
      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'OrbitNumber', Orbit)

      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'time_coverage_start', Tmp_String)
      !0        1         2
      !123456789012345678901234567890
      !2016-04-20T12:00:00.000Z
      read(Tmp_String(1:4), fmt="(I4)") Year
      read(Tmp_String(6:7), fmt="(I2)") Month
      read(Tmp_String(9:10), fmt="(I2)") Day
      read(Tmp_String(12:13), fmt="(I2)") Start_Hour
      read(Tmp_String(15:16), fmt="(I2)") Start_Minute
      read(Tmp_String(18:19), fmt="(I2)") Start_Sec

      !--- compute start day of year
      call JULIAN ( Day, Month, Year, Doy )
      
      call H5ReadAttribute( trim(Path)//trim(Infile), &
           'time_coverage_end', Tmp_String)
      read(Tmp_String(1:4), fmt="(I4)") End_Year
      read(Tmp_String(6:7), fmt="(I2)") Month
      read(Tmp_String(9:10), fmt="(I2)") Day
      read(Tmp_String(12:13), fmt="(I2)") End_Hour
      read(Tmp_String(15:16), fmt="(I2)") End_Minute
      read(Tmp_String(18:19), fmt="(I2)") End_Sec

      !--- compute end day of year
      call JULIAN ( Day, Month, End_Year, End_Doy )

      ! --- Calculate start and end time
      Start_Time = ((Start_Hour * 60 + Start_Minute) * 60 + Start_Sec) * 1000
      End_Time = ((End_Hour * 60 + End_Minute) * 60 + End_Sec) * 1000

end subroutine READ_VIIRS_NASA_DATE_TIME

!--------------------------------------------------------------------
! read viirs nasa number of scans
!--------------------------------------------------------------------
subroutine READ_NUMBER_OF_SCANS_VIIRS_NASA (Infile, Number_Of_Viirs_Lines, &
                 Error_Out)

      character(len=*), intent(in) :: Infile
      integer(kind=int4), intent(out) :: Error_Out
      integer(kind=int4), intent(out) :: Number_of_Viirs_Lines

      character(len=100) :: Setname
      integer, dimension(:), pointer :: Test
      integer, dimension(1) :: Dims
      integer :: Status


      Error_Out = 0
      Setname = 'scan_line_attributes/scan_quality'

      call H5ReadDataset( Infile, trim(Setname), Test )
      Dims = shape(Test)

      ! --- error 
      if (Dims(1) .lt. 1) then
         Number_Of_Viirs_Lines = -999
         Error_Out = 1
         Return
      endif
         
      Number_Of_Viirs_Lines = Dims(1) * 16

end subroutine READ_NUMBER_OF_SCANS_VIIRS_NASA

!--------------------------------------------------------------------
! find viirs nasa files
!--------------------------------------------------------------------
subroutine DETERMINE_VIIRS_NASA_FILE(Path_In,File_In,File_Type_In,File_Out)

      character(len=*), intent(in):: Path_In
      character(len=*), intent(in):: File_In
      character(len=*), intent(in):: File_Type_In
      character(len=*), intent(out):: File_Out

      character(len=500):: Search_String
      character(len=1020), dimension(:), pointer:: Files
      integer(kind=int4):: Num_Files


      !0        1         2         3         4         5
      !12345678901234567890123456789012345678901234567890
      !VGEOM_snpp_d20160420_t120000_c20160420175142.nc
      !VL1BM_snpp_d20160420_t120000_c20160420172929.nc

      Search_String = trim(File_Type_In)//trim(File_In(6:28))//'*.nc'

      Files => FILE_SEARCH(trim(Path_In),trim(Search_String),count=Num_Files)

      if (Num_Files == 0) then
         print *, EXE_PROMPT, VIIRS_NASA_PROMPT, "No "//trim(File_Type_In)//" File Found"
         File_Out = "no_file"
         return
      endif

      if (Num_Files > 1) then
         print *, EXE_PROMPT, VIIRS_NASA_PROMPT, "Multiple "//trim(File_Type_In)//" Files Found"
!         File_Out = "no_file"
      endif

      File_Out = Files(1)
!      print *, EXE_PROMPT, VIIRS_NASA_PROMPT, trim(File_Type_In)//" File Found, ",trim(File_Out)

      Files => null()

end subroutine DETERMINE_VIIRS_NASA_FILE

!-------------------------------------------------
! subroutine JULIAN(iday,imonth,iyear,jday)
! compute julian day
! input:
!         iday - integer day
!         imonth - integer month
!         iyear - integer year (2 or four digits)
!         
! output : jday - julian day
!--------------------------------------------------
subroutine JULIAN(iday,imonth,iyear,jday)

        integer, intent(in)::  iday,imonth,iyear
        integer, intent(out):: jday
        integer::  j
        integer, dimension(12)::  jmonth

        jmonth = reshape ((/31,28,31,30,31,30,31,31,30,31,30,31/),(/12/))

        jday = iday
        if (modulo(iyear,4) == 0) then
            jmonth(2)=29
        endif

        do j = 1,imonth-1
           jday = jday + jmonth(j)
        end do

end subroutine JULIAN

!--------------------------------------------------
! Function Name: CONVERT_VIIRS_RADIANCE
!
! Function:
!    Convert to units of the VIIRS radiance values from the that used
!    in the IDPS level-1b to that expected by CLAVR-x
!
! Description: 
!   
! Calling Sequence: rad_new =
! convert_viirs_radiance(rad_old,nu,missing_value)
!   
!
! Inputs:
!   rad_old = radiance in units of W/m^2/micron/str (2d array)
!   nu = channels equivalent width in units of cm^-1
!   missing_value = value assigned to missing radiance values
!
! Outputs: 
!   rad_new = radiance in units of mW/m^2/cm^-1/str (2d array)
!
! Dependencies:
!
! Restrictions:  None
!
! Reference: algebraic manipulation of Planck Equation
! ---------------------------------------------------------------------------------------
subroutine CONVERT_VIIRS_NASA_RADIANCE(Radiance,Nu,Missing_Value)
      real (kind=real4), dimension(:,:), intent(inout):: Radiance
      real (kind=real4), intent(in):: Nu
      real (kind=real4), intent(in):: Missing_Value

      where(Radiance /= Missing_Value)
         Radiance = Radiance * (((10000.0 / Nu )**2) / 10.0)
      end where

      return

end subroutine CONVERT_VIIRS_NASA_RADIANCE

!--------------------------------------------------
subroutine CONVERT_RAD_2_SOL_REF_DNB( Radiance, Solzen , Day_Of_Year , &
                                      Missing_Value , Reflectance )

      real (kind=real4), dimension(:,:), intent(in) :: Radiance
      real (kind=real4), dimension(:,:), intent(in) :: Solzen
      integer(kind=int2), intent(in) :: Day_Of_Year
      real (kind=real4), intent(in) :: Missing_Value
      real (kind=real4), dimension(:,:), intent(out) :: Reflectance
      real (kind=real4), parameter :: Fo = 0.044217282   ! dnb solar energy in w/cm2 ( Source? )
      real , parameter :: PI = 3.14159265359
      real , parameter :: DTOR = PI / 180.
      real :: Sun_Earth_Distance

      Sun_Earth_Distance = 1.0 - 0.016729 * cos( 0.9856 * ( Day_Of_Year - 4.0 ) * DTOR )

      Reflectance = Missing_Value

      where(Radiance /= Missing_Value)
         Reflectance = 100. * (Pi * Radiance * Sun_Earth_Distance ** 2) / &
                            (cos(Solzen * DTOR) * Fo)
      end where

end subroutine CONVERT_RAD_2_SOL_REF_DNB

!-------------------------------------------------------------------------------------
! subroutine to compute the bowtie gap pattern.  
! This pattern repeats every 48 scans 
! all VIIRS files should be integer multiples of this pattern
!
! Gap_Pixel_Mask_Pattern = a binary mask that identifies these bowtie gaps 
! Gap_Line_Idx = line index for each pixel in pattern including gap pixels
!
! A description of the pattern
!
!--------  line type 3
!----      line type 4
! 12 lines without gaps
!----      line type 1
!--------  line type 2
!--------  line type 3
!----      line type 4
! 12 lines without gaps
!----      line type 1
!--------  line type 2
!--------  line type 3
!----      line type 4
! 12 lines without gaps
!----      line type 1
!--------  line type 2
!
!  line types 1 and 4 have 640 missing pixels
!  line types 2 and 3 have 1008 missing pixels
! 
!-------------------------------------------------------------------------------------
   subroutine COMPUTE_VIIRS_BOWTIE_GAP_PATTERN()

      integer (kind=int4), dimension(Ny_Pattern):: Line_Type
      integer (kind=int4), parameter:: Ngap_1 = 640  !1280
      integer (kind=int4), parameter:: Ngap_2 = 1008 !2016
      integer (kind=int4), parameter:: Ngap_3 = 1008 !2016
      integer (kind=int4), parameter:: Ngap_4 = 640  !1280
      integer (kind=int4):: Iline
      integer (kind=int4):: i1
      integer (kind=int4):: i2

      !--- define the line patterns as described above
      Line_Type = &
      (/3,4,0,0,0,0,0,0,0,0,0,0,0,0,1,2/)


      if (.not. allocated(Gap_Line_Idx_Pattern)) &
           allocate(Gap_Line_Idx_Pattern(Nx_Pattern,Ny_Pattern))
      if (.not. allocated(Gap_Pixel_Mask_Pattern)) &
           allocate(Gap_Pixel_Mask_Pattern(Nx_Pattern,Ny_Pattern))

      do Iline = 1 , Ny_Pattern

         Gap_Line_Idx_Pattern( : , Iline ) = -999
         Gap_Pixel_Mask_Pattern( : , Iline ) = 0

         if (Line_Type(Iline) == 1) then

            i1 = 1
            i2 = Ngap_1
            Gap_Line_Idx_Pattern( i1 : i2 , Iline) = Iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline) = 1

            i1 = Nx_Pattern - Ngap_1 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline ) = 1
         end if                                                                                                                                                 
                                                                                                                                                                
         if (Line_Type(Iline) == 2) then
            i1 = 1
            i2 = Ngap_1
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline - 2
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline ) = 1

            i1 = Ngap_1 + 1
            i2 = Ngap_2
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline ) = 1

            i1 = Nx_Pattern - Ngap_1 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline - 2
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline ) = 1

            i1 = Nx_Pattern - Ngap_2 + 1
            i2 = i1  + (Ngap_2 - Ngap_1)
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline - 1
            Gap_Pixel_Mask_Pattern( i1 : i2 , Iline ) = 1
         end if

         if (Line_Type(Iline) == 3) then
            i1 = 1
            i2 = Ngap_4
            Gap_Line_Idx_Pattern(i1:i2,Iline) = Iline + 2
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1

            i1 = Ngap_4 + 1
            i2 = Ngap_3
            Gap_Line_Idx_Pattern(i1:i2,Iline) = Iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1

            i1 = Nx_Pattern - Ngap_4 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern(i1:i2,Iline) = Iline + 2
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1

            i1 = Nx_Pattern - Ngap_3 + 1
            i2 = i1 + (Ngap_3 - Ngap_4)
            Gap_Line_Idx_Pattern(i1:i2,Iline) = Iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1
         end if

         if (Line_Type(Iline) ==  4) then
            i1 = 1
            i2 = Ngap_4
            Gap_Line_Idx_Pattern(i1:i2,Iline) = Iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1

            i1 = Nx_Pattern - Ngap_4 + 1
            i2 = Nx_Pattern
            Gap_Line_Idx_Pattern( i1 : i2 , Iline ) = Iline + 1
            Gap_Pixel_Mask_Pattern(i1:i2,Iline) = 1
         end if

      end do

   end subroutine COMPUTE_VIIRS_BOWTIE_GAP_PATTERN

!------------------------------------------------------------------------------
! this routine uses the bowtie gap pattern, applies it to an arbitrary 
! segment of data and fills in the observations with the closest valid data
!------------------------------------------------------------------------------
    subroutine FILL_VIIRS_BOWTIE_GAPS ( Line_Start, Number_of_Lines  , Variable )
      integer(kind=int4), intent(in) :: Line_Start
      integer(kind=int4), intent(in) :: Number_of_Lines
      real(kind=real4), dimension(:,:), intent (inout) :: Variable

      integer(kind=int4) :: Line_Offset
      integer(kind=int4) :: Line_In_Pattern
      integer(kind=int4) :: Line_in_Segment
      integer(kind=int4) :: Line_Idx
      integer(kind=int4) :: Elem_Idx
      integer(kind=int4) :: Num_Pix = 3200
      integer(kind=int4), dimension(:,:), allocatable :: Gap_Line_Idx

      allocate (Gap_Line_Idx(3200, Number_Of_Lines))
      Gap_Line_Idx = Missing_Value_Int4

      do Line_In_Segment = 1,  Number_Of_Lines

         Line_In_Pattern = mod (Line_Start - 1 + Line_In_Segment , Ny_Pattern)
         if (Line_In_Pattern == 0) Line_In_Pattern = Ny_Pattern

         Line_Offset = Line_In_Segment - Line_In_Pattern

         if (Gap_Line_Idx_Pattern(1,Line_In_Pattern) /= Missing_Value_Int4) &
             Gap_Line_Idx(:, Line_In_Segment) = Gap_Line_Idx_Pattern(:,Line_In_Pattern) + Line_Offset

         where ( Gap_Line_Idx(:,Line_in_Segment) <= 0 )
            Gap_Line_Idx(:,Line_in_Segment) = 1
         end where

         where ( Gap_Line_Idx(:,Line_in_Segment) > Number_of_Lines )
            Gap_Line_Idx(:,Line_in_Segment) = Number_of_Lines
         end where

         ! - write gap mask to output
         Gap_Pixel_Mask(:,Line_in_Segment) = Gap_Pixel_Mask_Pattern(:,Line_in_Pattern)

         ! - fill gaps
         do Elem_Idx = 1, Num_Pix
            if (Gap_Pixel_Mask(Elem_Idx,Line_in_Segment) == 1) then
               Line_Idx = Gap_Line_Idx(Elem_Idx,Line_in_Segment)

               Variable(Elem_Idx, Line_In_Segment) = Variable(Elem_Idx,Line_Idx)

            end if
         end do
      end do

      deallocate (Gap_Line_Idx)

   end subroutine FILL_VIIRS_BOWTIE_GAPS

!====================================================================

end module VIIRS_NASA_READ_MODULE 


