!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: insat3_module.f90 (src)
!       INSAT3_module (program)
!
! PURPOSE: This module contains all the subroutines needed to perform navigation and
!          calibration for INSAT3D.
!
! DESCRIPTION: This module assumes band separated images only.  It also assumes
!              that the calibration type is INDI.  Calibration calculations are
!              taken from McIDAS-X 2016.1, kbxindi.dlm.
!
! NOTE: Calibration could change to a coefficient based model.  However, right
!       now those values are no good.
!
!       Navigation is of type GEOS.
!
! AUTHORS:
!  Steve Wanzong, CIMSS, stevew@ssec.wisc.edu
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
!--------------------------------------------------------------------------------------

module INSAT3_MODULE

use CONSTANTS
use PIXEL_COMMON
use CALIBRATION_CONSTANTS
use PLANCK
use CGMS_NAV
use NUMERICAL_ROUTINES
use GOES_MODULE
use FILE_UTILITY
use VIEWING_GEOMETRY_MODULE

implicit none
public:: READ_INSAT3
public:: READ_NAVIGATION_BLOCK_INSAT3
public:: READ_INSAT3_INSTR_CONSTANTS
         
private :: INSAT3_RADIANCE_BT, INSAT3_NAVIGATION, INSAT3_Reflectance
 
type (GVAR_NAV), PRIVATE    :: NAVstr_INSAT3_NAV
integer, PARAMETER, PRIVATE :: Nchan_INSAT3= 6
integer, PARAMETER, PRIVATE :: Ndet_INSAT3 = 4
integer, PARAMETER, PRIVATE :: Ntable_INSAT3 = 1024

integer, PRIVATE :: Nref_Table_INSAT3
integer, PRIVATE :: Nbt_Table_INSAT3
CHARACTER(len=4), PRIVATE:: Calib_Type

!real (kind=real4), dimension(Ntable_INSAT3), PRIVATE  :: Ref_Table
real (kind=real4), dimension(2,Ntable_INSAT3), PRIVATE  :: Ref_Table
integer (kind=int4), dimension(nchan_INSAT3,ndet_INSAT3,Ntable_INSAT3), PRIVATE  :: bt_table
integer (kind=int4), dimension(nchan_INSAT3,ndet_INSAT3,Ntable_INSAT3), PRIVATE  :: rad_table

integer(kind=int4), private, parameter:: INSAT3_Xstride = 1
integer(kind=int4), private, parameter:: num_4km_scans_fd = 2816
integer(kind=int4), private, parameter:: num_4km_elem_fd = 2808
integer(kind=int4), private, parameter:: time_for_fd_scan =  1560000 !milliseconds (26min)
real, private, save:: Scan_rate    !scan rate in millsec / line
integer(kind=int4), private, parameter:: IN3_Byte_Shift = 0 !number of bytes to shift for FY2

integer(kind=int4), private :: ptr

integer(kind=int4), dimension(64) :: i4buf

!--- Need to set up individual file names.
character(len=1020):: channel_2_filename
character(len=1020):: channel_3_filename
character(len=1020):: channel_4_filename
character(len=1020):: channel_5_filename
character(len=1020):: channel_6_filename

CONTAINS

!----------------------------------------------------------------
! read the INSAT3 constants into memory
!-----------------------------------------------------------------
subroutine READ_INSAT3_INSTR_CONSTANTS(Instr_Const_file)
  character(len=*), intent(in):: Instr_Const_file
  integer:: ios0, erstat
  integer:: Instr_Const_lun

  Instr_Const_lun = GET_LUN()
  open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

  print *, "opening ", trim(Instr_Const_file)
  erstat = 0
  if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, "Error opening INSAT3 constants file, ios0 = ", ios0
    stop 19
  endif
  read(unit=Instr_Const_lun,fmt="(a7)") sat_name
  read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
  read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
  read(unit=Instr_Const_lun,fmt=*) planck_a1(20), planck_a2(20), planck_nu(20)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(27), planck_a2(27), planck_nu(27)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(31), planck_a2(31), planck_nu(31)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(32), planck_a2(32), planck_nu(32)

  read(unit=Instr_Const_lun,fmt=*) Ch1_Dark_Count
  read(unit=Instr_Const_lun,fmt=*) Ch1_Gain_Low_0,Ch1_Degrad_Low_1, Ch1_Degrad_Low_2
  read(unit=Instr_Const_lun,fmt=*) Launch_Date

  close(unit=Instr_Const_lun)

  !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / ew_Ch20

  !-- hardwire ch1 dark count
  Ch1_Dark_Count = 29

end subroutine READ_INSAT3_INSTR_CONSTANTS

 ! Perform INSAT3 Reflectance and BT calibration
 subroutine READ_INSAT3(segment_number,channel_1_filename, &
                     jday, image_time_ms, &
                     AREAstr,NAVstr_INSAT3)

   integer(kind=int4), intent(in):: segment_number
   character(len=*), intent(in):: channel_1_filename
   type (AREA_STRUCT), intent(in) :: AREAstr
   type (GVAR_NAV), intent(in)    :: NAVstr_INSAT3
   integer(kind=int2), intent(in):: jday
   integer(kind=int4), intent(in):: image_time_ms

   character(len=1020):: channel_x_filename
   character(len=1020):: channel_x_filename_full
   character(len=1020):: channel_x_filename_full_uncompressed
   character(len=180):: System_String
   integer:: ipos
   integer:: ilen
   integer:: ichan_goes
   integer:: ichan_modis
   integer:: INSAT3_file_id
   real(kind=real4):: image_time_hours
   integer(kind=int4):: image_jday
   integer(kind=int4):: first_line_in_segment
   character(len=2):: ichan_goes_string
   integer :: Line_Idx
   integer :: Elem_Idx
   integer:: num_elements_this_image
   integer:: num_scans_this_image
   integer(kind=int2),  parameter::  Chan1_int2 = 1
   integer(kind=int2),  parameter::  Chan2_int2 = 2

   !--- assume channel_1_file name has a unique "_1_" in the name. 
   !--- determine indices needed to replace that string
   ipos = index(channel_1_filename, "_1_")
   ilen = len(channel_1_filename)
    
   first_line_in_segment = (segment_number-1)*Image%Number_Of_Lines_Per_Segment

   !---------------------------------------------------------------------------
   ! INSAT3 Navigation (Do Navigation and Solar angles first)
   !---------------------------------------------------------------------------
   
   call INSAT3_NAVIGATION(1,first_line_in_segment,&
                              Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment,1,&
                              AREAstr,NAVstr_INSAT3)
   
   if (segment_number == 1) then

     image_jday = jday
     image_time_hours = image_time_ms / 60.0 / 60.0 / 1000.0

     !--- compute scan rate for future use
     num_elements_this_image =  int(AREAstr%num_elem / INSAT3_Xstride) + 1
     num_scans_this_image = AREAstr%num_line
     Scan_Rate = real((num_elements_this_image)/               &
       real(num_4km_elem_fd/INSAT3_Xstride)) * &
       real((num_scans_this_image) / real(num_4km_scans_fd)) * &
       real(time_for_fd_scan) / real(num_scans_this_image)

     !--- Need to loop over all channels as the calibration block is unique for
     !--- each band.
     do ichan_goes = 1,6

       if (ichan_goes == 1) ichan_modis = 1
       if (ichan_goes == 2) ichan_modis = 6
       if (ichan_goes == 3) ichan_modis = 20
       if (ichan_goes == 4) ichan_modis = 27
       if (ichan_goes == 5) ichan_modis = 31
       if (ichan_goes == 6) ichan_modis = 32
       
       write(ichan_goes_string,fmt="(I1.1)") ichan_goes
       if(ichan_goes > 9) write(ichan_goes_string,fmt="(I2.2)") ichan_goes

       if (Sensor%Chan_On_Flag_Default(ichan_modis) == sym%YES) then

         channel_x_filename = channel_1_filename(1:ipos-1) // "_"//trim(ichan_goes_string)//"_" // &
                            channel_1_filename(ipos+3:ilen)

         if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
           channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
         else
           channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
         endif

         channel_x_filename_full_uncompressed = trim(Image%Level1b_Path)//trim(channel_x_filename)

         if (l1b_gzip == sym%YES) then
           System_String = "gunzip -c "//trim(channel_x_filename_full_uncompressed)//".gz"// &
                                " > "//trim(channel_x_filename_full)
                                
           call system(System_String)

           Number_of_Temporary_Files = Number_of_Temporary_Files + 1
           Temporary_File_Name(Number_of_Temporary_Files) = trim(channel_x_filename)

         endif
         if (l1b_bzip2 == sym%YES) then
           System_String = "bunzip2 -c "//trim(channel_x_filename_full_uncompressed)//".bz2"// &
                              " > "//trim(channel_x_filename_full)
           call system(System_String)

           Number_of_Temporary_Files = Number_of_Temporary_Files + 1
           Temporary_File_Name(Number_of_Temporary_Files) = trim(channel_x_filename)
         endif

       endif

     enddo
     
     !--- On first segment, reflectance, BT and rad tables.  For INSAT3D, we will
     !--- need to read each AREA file separately, as the calibration block is
     !--- unique for each AREA file.

     !--- Band 1 - 0.65 um
     INSAT3_file_id = get_lun()   
     print*,"Channel 1 calibration for : ", trim(channel_1_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_1_filename)//CHAR(0), INSAT3_file_id)
     else
        call mread_open(trim(Image%Level1b_Path)//trim(channel_1_filename)//CHAR(0), INSAT3_file_id)
     endif  

     call load_INSAT3_calibration(INSAT3_file_id, AREAstr)
     call mread_close(INSAT3_file_id)

     !--- Band 2 - 1.62 um
     INSAT3_file_id = get_lun()   
     channel_2_filename = channel_1_filename(1:ipos-1) // "_2_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 2 calibrtaion for : ", trim(channel_2_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_2_filename)//CHAR(0), INSAT3_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_2_filename)//CHAR(0), INSAT3_file_id)
     endif  

     call load_INSAT3_calibration(INSAT3_file_id, AREAstr)
     call mread_close(INSAT3_file_id)

     !--- Band 3 - 3.9 um
     INSAT3_file_id = get_lun()   
     channel_3_filename = channel_1_filename(1:ipos-1) // "_3_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 3 calibration for  : ", trim(channel_3_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_3_filename)//CHAR(0), INSAT3_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_3_filename)//CHAR(0), INSAT3_file_id)
     endif  

     call load_INSAT3_calibration(INSAT3_file_id, AREAstr)
     call mread_close(INSAT3_file_id)

     !--- Band 4 - 6.8 um
     INSAT3_file_id = get_lun()   
     channel_4_filename = channel_1_filename(1:ipos-1) // "_4_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 4 calibration for  : ", trim(channel_4_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_4_filename)//CHAR(0), INSAT3_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_4_filename)//CHAR(0), INSAT3_file_id)
     endif  

     call load_INSAT3_calibration(INSAT3_file_id, AREAstr)
     call mread_close(INSAT3_file_id)

     !--- Band 5 - 10.8 um
     INSAT3_file_id = get_lun()   
     channel_5_filename = channel_1_filename(1:ipos-1) // "_5_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 5 calibration for  : ", trim(channel_5_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_5_filename)//CHAR(0), INSAT3_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_5_filename)//CHAR(0), INSAT3_file_id)
     endif  

     call load_INSAT3_calibration(INSAT3_file_id, AREAstr)
     call mread_close(INSAT3_file_id)

     !--- Band 6 - 12.0 um
     INSAT3_file_id = get_lun()   
     channel_6_filename = channel_1_filename(1:ipos-1) // "_6_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 6 calibration for : ", trim(channel_6_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_6_filename)//CHAR(0), INSAT3_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_6_filename)//CHAR(0), INSAT3_file_id)
     endif  

     call load_INSAT3_calibration(INSAT3_file_id, AREAstr)
     call mread_close(INSAT3_file_id)

   endif !--- Segment 1
   
   !---   read channel 1 (INSAT3 channel 1)
   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_1_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_1_filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    IN3_Byte_Shift, &
                                    AREAstr, INSAT3_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)

       !--- Make reflectance table.
       call INSAT3_Reflectance(Two_Byte_Temp,ch(1)%Ref_Toa(:,:),Chan1_int2)

       !--- store ch1 counts for support of PATMOS-x calibration studies
       Ch1_Counts = Two_Byte_Temp

   endif

   !---   read channel 6 (INSAT3 channel 2)
   if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_2_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_2_filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    IN3_Byte_Shift, &
                                    AREAstr, INSAT3_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)

       !--- Make reflectance table.
       call INSAT3_Reflectance(Two_Byte_Temp,ch(6)%Ref_Toa(:,:),Chan2_int2)

       !--- store ch2 counts for support of PATMOS-x calibration studies
       Ch2_Counts = Two_Byte_Temp

   endif
   
   !---   read channel 20 (INSAT3 channel 3)
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_3_" // &
                            channel_1_filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    IN3_Byte_Shift, &
                                    AREAstr, INSAT3_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)

       call INSAT3_RADIANCE_BT(3_int1, Two_Byte_Temp, ch(20)%Rad_Toa, ch(20)%Bt_Toa)

   endif
                    
   
   !---   read channel 27 (INSAT3 channel 4)
   if (Sensor%Chan_On_Flag_Default(27) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_4_" // &
                            channel_1_filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif

      call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    IN3_Byte_Shift, &
                                    AREAstr, INSAT3_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)

      call INSAT3_RADIANCE_BT(4_int1, Two_Byte_Temp, ch(27)%Rad_Toa, ch(27)%Bt_Toa)

   endif
   
   !---   read channel 31 (INSAT3 channel 5)
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_5_" // &
                            channel_1_filename(ipos+3:ilen)
       
       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    IN3_Byte_Shift, &
                                    AREAstr, INSAT3_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)
       
      call INSAT3_RADIANCE_BT(5_int1, Two_Byte_Temp, ch(31)%Rad_Toa, ch(31)%Bt_Toa)

   endif
   
   
   !---   read channel 32 (INSAT3 channel 6)
   if (Sensor%Chan_On_Flag_Default(32) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_6_" // &
                            channel_1_filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    IN3_Byte_Shift, &
                                    AREAstr, INSAT3_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)

      call INSAT3_RADIANCE_BT(6_int1, Two_Byte_Temp, ch(32)%Rad_Toa, ch(32)%Bt_Toa)

   endif
    
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     Scan_Number(Line_Idx) = first_line_in_segment + Line_Idx
     Scan_Time(Line_Idx) = image_time_ms + (Scan_Number(Line_Idx)-1) * Scan_rate
   enddo

!------------------------------------------------------------------------------
! INSAT3 Angles
! NOTE: These were private routines in the GOES module. Suggest they become
!       public with different names, since they are used cross platform  
!------------------------------------------------------------------------------
   image_jday = jday
   image_time_hours = image_time_ms / 60.0 / 60.0 / 1000.0
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     do Elem_Idx = 1,Image%Number_Of_Elements
        call POSSOL(image_jday,image_time_hours, &
                    Nav%Lon_1b(Elem_Idx,Line_Idx),Nav%Lat_1b(Elem_Idx,Line_Idx), &
                    Geo%Solzen(Elem_Idx,Line_Idx),Geo%Solaz(Elem_Idx,Line_Idx))
     enddo
      call COMPUTE_SATELLITE_ANGLES(Sensor%Geo_Sub_Satellite_Longitude,  &
                                    Sensor%Geo_Sub_Satellite_Latitude, Line_Idx)
   enddo

   !--- Ascending node
   Elem_Idx = Image%Number_Of_Elements/2
   do Line_Idx = Line_Idx_Min_Segment+1, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     Nav%Ascend(Line_Idx) = 0
     if (Nav%Lat_1b(Elem_Idx,Line_Idx) < Nav%Lat_1b(Elem_Idx,Line_Idx-1)) then
       Nav%Ascend(Line_Idx) = 1
     endif
   enddo
   Nav%Ascend(Line_Idx_Min_Segment) = Nav%Ascend(Line_Idx_Min_Segment+1)
    
 end subroutine READ_INSAT3
 
 subroutine load_INSAT3_calibration(lun, AREAstr)
  integer(kind=int4), intent(in) :: lun
  type(AREA_STRUCT), intent(in):: AREAstr
  integer(kind=int4), dimension(2048) :: ibuf !kbxindi.dlm
  integer :: nref, nbt, i, offset
  integer(kind=int4) :: band_offset_2, band_offset_14, band_offset_15, &
                        band_offset_9, band_offset_7, dir_offset, band_offset_5, &
                        squared_factor, linear_factor
  real(kind=real4) :: albedo, temperature, radiance, squared_coef, linear_coef
  integer(kind=int4), parameter:: MAXVAL=1023
  integer(kind=int4) :: local_sndr_filter_map

  !--- Store the correct calibration source type in the AREA structure.  It will be
  !--- important if the calibration type changes.
  call mreadf_int_o(lun,0,4,64,i4buf)
  call move_bytes(4,i4buf(1+51:1+51),AREAstr%src_Type,0)
  local_sndr_filter_map = i4buf(19) ! Need band map for current open AREA file.

  if (AREAstr%src_Type .eq. 'INDI') then
    print*,"INDI calibration, proceeding ...", AREAstr%num_chan, AREAstr%sndr_filter_map
  else
    print*,"Unknown INSAT3D calibration type, exiting ..."
    stop
  endif

  !--- Single band images are required.
  if (AREAstr%num_chan > 1) then
    print*,"Multi-band image, exiting ..."
    stop
  endif

  !---Read band map to determine which type of calibration table we are reading.
  if (local_sndr_filter_map == 1) then

    !Band 1 Albedo.
    ptr=0
    call mreadf_int_o(lun,AREAstr%cal_offset+ptr,4,(MAXVAL+1),ibuf)
    do i = 1, MAXVAL+1
      Ref_Table(1,i) = ibuf(i+1)/10.0 ! Match McIDAS-X
    end do
    !--- Read Band 1 Radiances (.65 um rad)
    !--- Don't do anything with this yet.
    ptr=ptr+((MAXVAL+1)*4)
    call mreadf_int_o(lun,AREAstr%cal_offset+ptr,4,(MAXVAL+1),ibuf)
    !---*** This does not capture McIDAS 999.00.  Only to 998.0

  elseif (local_sndr_filter_map == 2) then

    !Band 2 Radiance.
    !--- Don't do anything with this yet.
    ptr=0
    call mreadf_int_o(lun,AREAstr%cal_offset+ptr,4,(MAXVAL+1),ibuf)
    !--- Read Band 2 Albedo (1.62 um Albedo)
    ptr=ptr+((MAXVAL+1)*4)
    call mreadf_int_o(lun,AREAstr%cal_offset+ptr,4,(MAXVAL+1),ibuf)
    do i = 1, MAXVAL+1
      Ref_Table(2,i) = ibuf(i+1)/10.0 ! Match McIDAS-X
    end do

  elseif (local_sndr_filter_map == 4 .OR. local_sndr_filter_map == 8 .OR. &
          local_sndr_filter_map == 16 .OR. local_sndr_filter_map == 32) then

    !Band 3-6 Radiances
    ptr=0
    call mreadf_int_o(lun,AREAstr%cal_offset+ptr,4,(MAXVAL+1),ibuf)
    do i = 1, MAXVAL+1
      select case(local_sndr_filter_map)
        case(4)
          rad_table(3,1,i) = ibuf(i+1)
        case(8)
          rad_table(4,1,i) = ibuf(i+1)
        case(16)
          rad_table(5,1,i) = ibuf(i+1)
        case(32)
          rad_table(6,1,i) = ibuf(i+1)
      end select
    end do
    !--- Read BT's
    ptr=ptr+((MAXVAL+1)*4)
    call mreadf_int_o(lun,AREAstr%cal_offset+ptr,4,(MAXVAL+1),ibuf)
    do i = 1, MAXVAL+1
      select case(local_sndr_filter_map)
        case(4)
          bt_table(3,1,i) = ibuf(i+1)
          temperature = real(bt_table(3,1,i),kind=real4) / 100.0
          radiance = PLANCK_RAD(20,temperature)
          rad_table(3,1,i) = nint(radiance * 1000.)
        case(8)
          bt_table(4,1,i) = ibuf(i+1)
          temperature = real(bt_table(4,1,i),kind=real4) / 100.0
          radiance = PLANCK_RAD(27,temperature)
          rad_table(4,1,i) = nint(radiance * 1000.)
        case(16)
          bt_table(5,1,i) = ibuf(i+1)
          temperature = real(bt_table(5,1,i),kind=real4) / 100.0
          radiance = PLANCK_RAD(31,temperature)
          rad_table(5,1,i) = nint(radiance * 1000.)
        case(32)
          bt_table(6,1,i) = ibuf(i+1)
          temperature = real(bt_table(6,1,i),kind=real4) / 100.0
          radiance = PLANCK_RAD(32,temperature)
          rad_table(6,1,i) = nint(radiance * 1000.)
      end select

    end do

  endif

  !---stwNref = 256
  !---stwNbt = 1024
  !---stwNref_Table_INSAT3 = nref
  !---stwNbt_Table_INSAT3 = nbt
  
 end subroutine load_INSAT3_calibration

 ! Perform INSAT3 Navigation

 subroutine INSAT3_NAVIGATION(xstart,ystart,xsize,ysize,xstride, &
                            AREAstr,NAVstr_INSAT3)
    integer(kind=int4) :: xstart, ystart
    integer(kind=int4) :: xsize, ysize
    integer(kind=int4) :: xstride  
    type (AREA_STRUCT) :: AREAstr
    type (GVAR_NAV), intent(in)    :: NAVstr_INSAT3
    
    integer :: i, j, ii, jj, imode
    real(kind(0.0d0)) :: latitude, longitude
    real(kind=real4) :: height
    integer :: FGF_type = 3 !INSAT uses JMA GEOS navigation, so set type here

    NAVstr_INSAT3_NAV = NAVstr_INSAT3
    
    imode = -1
    height = 0.0   !Used for parallax correction
    
    Nav%Lat = Missing_Value_Real4
    Nav%Lon = Missing_Value_Real4
    
    if (NAVstr_INSAT3%nav_type == 'GEOS') then      
        !HRIT requires actual line and element of being processed.
                
          do j=1, ysize
            
            jj = ystart + (j-1) + (AREAstr%north_bound / real(AREAstr%line_res))
            
            do i=1, xsize
                ii = (i - 1)*(xstride) + xstart ! get element of the image segement
                ii = ii  + (AREAstr%west_vis_pixel / real(AREAstr%elem_res))

                CALL fgf_to_earth(FGF_type,                  &
                                  DBLE(ii),                  &
                                  DBLE(jj),                  &
                                  DBLE(NAVstr_INSAT3%CFAC),    &
                                  DBLE(NAVstr_INSAT3%COFF),    &
                                  DBLE(NAVstr_INSAT3%LFAC),    &
                                  DBLE(NAVstr_INSAT3%LOFF),    &
                                  DBLE(NAVstr_INSAT3%sub_lon), &
                                  longitude,                 &
                                  latitude)

             if (latitude .LE. -999.0) then  ! -999.99 is MSV nav missing value
                    Nav%Lat_1b(i,j) = Missing_Value_Real4
                    Nav%Lon_1b(i,j) = Missing_Value_Real4
                    Space_Mask(i,j) = sym%SPACE
                else
                    Nav%Lat_1b(i,j) = real(latitude,kind=real4)
                    Nav%Lon_1b(i,j) = real(longitude,kind=real4)

                    ! we want 180 to -180, one last check.
                    
                    if (longitude .GT. 180.0 ) then
                        Nav%Lon_1b(i,j) = real(longitude,kind=real4) - 360.0
                    endif
                                        
                    Space_Mask(i,j) = sym%NO_SPACE
                endif
        
            end DO
                        
        end do     
        
    endif
      
 end subroutine INSAT3_NAVIGATION
 
 
!------------------------------------------------------------------
! subroutine to convert INSAT3 counts to radiance and brightness
! temperature
!------------------------------------------------------------------
  
  subroutine INSAT3_RADIANCE_BT(chan_num,INSAT3_Counts, rad2, temp1)

    integer (kind=INT2), dimension(:,:), intent(in):: INSAT3_Counts
    integer (kind=int1), INTENT(in) :: chan_num
    real (kind=real4), dimension(:,:), INTENT(out):: temp1, rad2
    
    integer :: i, j, index

    do j = 1, Image%Number_Of_Lines_Read_This_Segment
      do i = 1, Image%Number_Of_Elements
        if (Space_Mask(i,j) == sym%NO_SPACE) then
          index = int(INSAT3_Counts(i,j),kind=int2) + 1
          rad2(i,j) = real(rad_table(chan_num,1,index),kind=real4)/1000.0
          temp1(i,j) = real(bt_table(chan_num,1,index),kind=real4)/100.0                    
        else
          rad2(i,j) = Missing_Value_Real4
          temp1(i,j) = Missing_Value_Real4
        endif
        if (chan_num == 6) then
        endif
      end DO
    end do    
  
  end subroutine INSAT3_RADIANCE_BT
  
!------------------- INSAT3 NAV BLOC

 subroutine READ_NAVIGATION_BLOCK_INSAT3(filename, AREAstr, NAVstr)
  CHARACTER(len=*), intent(in):: filename
  type(AREA_STRUCT), intent(in):: AREAstr
  type(GVAR_NAV), intent(inout):: NAVstr
 
  integer :: geos_nav
  integer(kind=int4)nav_offset
  integer:: number_of_words_read
  integer(kind=int4), dimension(640) :: i4buf
  
  nav_offset = AREAstr%sec_key_nav
    
  !determine GEOS or other navigation
  geos_nav = sym%NO
  call mreadf_int(trim(filename)//CHAR(0),nav_offset,4,640,&
                    number_of_words_read, i4buf)
!  if (AREAstr%swap_bytes > 0) call swap_bytes4(i4buf,640)
  call move_bytes(4,i4buf(1),NAVstr%nav_type,0)
 !SUBLON stored as SUBLON *10 in McIDAS NAV block
  NAVstr%sub_lon = real(i4buf(6),kind=real4) / 10
  NAVstr%sublon = NAVstr%sub_lon

 ! LOFF, COFF, CFAC, LFAC stored in McIDAS header for 1km data. All
 ! Multipied by 10. Order from nvxmtst.dlm in McIDAS
  NAVstr%LOFF=(i4buf(2) / 10 ) / real(AREAstr%line_res)
  NAVstr%COFF=(i4buf(3) / 10) / real(AREAstr%elem_res)
  NAVstr%LFAC=(i4buf(4) / 10 ) / real(AREAstr%line_res)
  NAVstr%CFAC=(i4buf(5) / 10 ) / real(AREAstr%elem_res)

 end subroutine READ_NAVIGATION_BLOCK_INSAT3

 ! Perform Insat3 Reflectance calculation
 subroutine INSAT3_Reflectance(Insat3_Counts, Alb_Temp,chan_num)

    integer (kind=INT2), dimension(:,:), intent(in):: Insat3_Counts
    integer (kind=INT2), intent(in):: chan_num
    real (kind=real4), dimension(:,:), intent(out):: Alb_Temp

    integer :: i, j
    integer :: index

    DO j=1, Image%Number_Of_Lines_Read_This_Segment
      DO i=1, Image%Number_Of_Elements
        IF (space_mask(i,j) == sym%NO_SPACE .and. Geo%solzen(i,j) < 90.0) THEN

          !Not sure if I need to add 1 here
          index = int(Insat3_Counts(i,j),kind=int2)

          alb_temp(i,j) = Missing_Value_Real4

          IF((index .GT. 0) .AND. (index .LE. 1024)) THEN
            Alb_Temp(i,j) = (real(Ref_Table(chan_num,index),kind=real4))
          ENDIF

        ELSE
         Alb_Temp(i,j) = Missing_Value_Real4
        ENDif
      END DO
    END DO

 end subroutine INSAT3_Reflectance

end module INSAT3_MODULE
