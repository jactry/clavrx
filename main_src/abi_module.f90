!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: abi_module.f90 (src)
!       ABI_MODULE (program)
!
! PURPOSE: This module contains all the subroutines needed to perform navigation and
!          calibration for GOES-16 ABI.
!
! DESCRIPTION: This module assumes band separated images only.  It also assumes
!              that the calibration type is ABIN.  
!
!       Navigation is of type ABIN.
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
! AHI Channel Mapping
!
!  wvl        ahi  modis/clavrx
!
!  0.47       1      3
!  0.64       2      1
!  0.87       3      2
!  1.38       4     26
!  1.61       5      6
!  2.25       6      7
!  3.90       7     20
!  6.19       8     37
!  6.95       9     27
!  7.34      10     28
!  8.50      11     29
!  9.61      12     30
!  10.4      13     38
!  11.2      14     31
!  12.3      15     32
!  13.3      16     33
!
!--------------------------------------------------------------------------------------

module ABI_MODULE

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
public:: READ_ABI
public:: READ_NAVIGATION_BLOCK_ABI
public:: READ_ABI_INSTR_CONSTANTS
         
private :: ABI_RADIANCE_BT, ABI_NAVIGATION, ABI_Reflectance, GET_ABI_IMAGE
 
type (GVAR_NAV), PRIVATE    :: NAVstr_ABI_NAV
integer, PARAMETER, PRIVATE :: Nchan_ABI= 16
integer, PARAMETER, PRIVATE :: Ndet_ABI = 1
integer, PARAMETER, PRIVATE :: Ntable_ABI = 65536

integer, PRIVATE :: Nref_Table_ABI
integer, PRIVATE :: Nbt_Table_ABI
CHARACTER(len=4), PRIVATE:: Calib_Type

!real (kind=real4), dimension(Ntable_ABI), PRIVATE  :: Ref_Table
real (kind=real4), dimension(6,Ntable_ABI), PRIVATE  :: Ref_Table
integer (kind=int4), dimension(16,ndet_ABI,Ntable_ABI), PRIVATE  :: bt_table
integer (kind=int4), dimension(16,ndet_ABI,Ntable_ABI), PRIVATE  :: rad_table

integer(kind=int4), private, parameter:: ABI_Xstride = 1
integer(kind=int4), private, parameter:: num_4km_scans_fd = 5424
integer(kind=int4), private, parameter:: num_4km_elem_fd = 5424
integer(kind=int4), private, parameter:: time_for_fd_scan =  636000 !milliseconds (10.6 min,)
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
character(len=1020):: channel_7_filename
character(len=1020):: channel_8_filename
character(len=1020):: channel_9_filename
character(len=1020):: channel_10_filename
character(len=1020):: channel_11_filename
character(len=1020):: channel_12_filename
character(len=1020):: channel_13_filename
character(len=1020):: channel_14_filename
character(len=1020):: channel_15_filename
character(len=1020):: channel_16_filename

CONTAINS

!----------------------------------------------------------------
! read the ABI constants into memory
!-----------------------------------------------------------------
subroutine READ_ABI_INSTR_CONSTANTS(Instr_Const_file)
  character(len=*), intent(in):: Instr_Const_file
  integer:: ios0, erstat
  integer:: Instr_Const_lun

  Instr_Const_lun = GET_LUN()
  open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

  print *, "opening ", trim(Instr_Const_file)
  erstat = 0
  if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, "Error opening ABI constants file, ios0 = ", ios0
    stop 19
  endif
  read(unit=Instr_Const_lun,fmt="(a7)") sat_name
  read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
  read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
  read(unit=Instr_Const_lun,fmt=*) planck_a1(20), planck_a2(20), planck_nu(20)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(27), planck_a2(27), planck_nu(27)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(28), planck_a2(28), planck_nu(28)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(29), planck_a2(29), planck_nu(29)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(30), planck_a2(30), planck_nu(30)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(31), planck_a2(31), planck_nu(31)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(32), planck_a2(32), planck_nu(32)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(33), planck_a2(33), planck_nu(33)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(37), planck_a2(37), planck_nu(37)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(38), planck_a2(38), planck_nu(38)
  close(unit=Instr_Const_lun)

  !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / ew_Ch20

  !-- hardwire ch1 dark count
  Ch1_Dark_Count = 29

end subroutine READ_ABI_INSTR_CONSTANTS

!-------------------------------------------------------------------------------
! Private routine to read data from an AREA file for one segment into memory
!-------------------------------------------------------------------------------
subroutine READ_ABI(segment_number,channel_1_filename, &
                    jday, image_time_ms, &
                    AREAstr,NAVstr_ABI)

  integer(kind=int4), intent(in):: segment_number
  character(len=*), intent(in):: channel_1_filename
  type (AREA_STRUCT), intent(in) :: AREAstr
  type (GVAR_NAV), intent(in)    :: NAVstr_ABI
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
  integer:: ABI_file_id
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
  ! ABI Navigation (Do Navigation and Solar angles first)
  !---------------------------------------------------------------------------
   
  call ABI_NAVIGATION(1,first_line_in_segment,&
                             Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment,1,&
                             AREAstr,NAVstr_ABI)
   
   if (segment_number == 1) then

     Ref_Table = missing_value_int4

     image_jday = jday
     image_time_hours = image_time_ms / 60.0 / 60.0 / 1000.0

     !--- FIXME - Will depend on which ABI Scan Mode is uses.  This should
     !--- be in the AREA file in an upcoming McIDAS release.

     !--- compute scan rate for future use
     num_elements_this_image =  int(AREAstr%num_elem / ABI_Xstride) + 1
     num_scans_this_image = AREAstr%num_line
     Scan_rate = real((num_elements_this_image)/               &
       real(num_4km_elem_fd/ABI_Xstride)) * &
       real((num_scans_this_image) / real(num_4km_scans_fd)) * &
       real(time_for_fd_scan) / real(num_scans_this_image)
    
     !---stw Debug
     !print*,"num_elements_this_image : ", num_elements_this_image
     !print*,"num_scans_this_image : ", num_scans_this_image
     !print*,"num_4km_elem_fd : ", num_4km_elem_fd
     !print*,"num_4km_scans_fd : ", num_4km_scans_fd
     !print*,"time_for_fd_scan : ", time_for_fd_scan
     !print*,"Scan_rate: ", Scan_rate
     !stop

     !--- Need to loop over all channels as the calibration block is unique for
     !--- each band.
     do ichan_goes = 1,16

       if (ichan_goes == 1)  ichan_modis = 3
       if (ichan_goes == 2)  ichan_modis = 1
       if (ichan_goes == 3)  ichan_modis = 2
       if (ichan_goes == 4)  ichan_modis = 26
       if (ichan_goes == 5)  ichan_modis = 6
       if (ichan_goes == 6)  ichan_modis = 7
       if (ichan_goes == 7)  ichan_modis = 20
       if (ichan_goes == 8)  ichan_modis = 37
       if (ichan_goes == 9)  ichan_modis = 27
       if (ichan_goes == 10) ichan_modis = 28
       if (ichan_goes == 11) ichan_modis = 29
       if (ichan_goes == 12) ichan_modis = 30
       if (ichan_goes == 13) ichan_modis = 38
       if (ichan_goes == 14) ichan_modis = 31
       if (ichan_goes == 15) ichan_modis = 32
       if (ichan_goes == 16) ichan_modis = 33
       
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
     
     !--- On first segment grab the calibration block from the AREA file.
     !--- Need to read each AREA file separately, as the calibration block is
     !--- unique for each AREA file.

     !--- Band 1 - 0.47 um
     ABI_file_id = get_lun()   
     print*,"Channel 1 calibration for : ", trim(channel_1_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_1_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_1_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 3
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 2 - 0.54 um
     ABI_file_id = get_lun()   
     channel_2_filename = channel_1_filename(1:ipos-1) // "_2_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 2 calibration for : ", trim(channel_2_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_2_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_2_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 1
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 3 - 0.87 um
     ABI_file_id = get_lun()   
     channel_3_filename = channel_1_filename(1:ipos-1) // "_3_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 3 calibration for  : ", trim(channel_3_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_3_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_3_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 2
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 4 - 1.38 um
     ABI_file_id = get_lun()   
     channel_4_filename = channel_1_filename(1:ipos-1) // "_4_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 4 calibration for  : ", trim(channel_4_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_4_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_4_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 26
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 5 - 1.61 um
     ABI_file_id = get_lun()   
     channel_5_filename = channel_1_filename(1:ipos-1) // "_5_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 5 calibration for  : ", trim(channel_5_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_5_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_5_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 6
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 6 - 2.25 um
     ABI_file_id = get_lun()   
     channel_6_filename = channel_1_filename(1:ipos-1) // "_6_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 6 calibration for : ", trim(channel_6_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_6_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_6_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 7
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 7 - 3.90 um
     ABI_file_id = get_lun()   
     channel_7_filename = channel_1_filename(1:ipos-1) // "_7_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 7 calibration for : ", trim(channel_7_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_7_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_7_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 20
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 8 - 6.19 um
     ABI_file_id = get_lun()   
     channel_8_filename = channel_1_filename(1:ipos-1) // "_8_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 8 calibration for : ", trim(channel_8_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_8_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_8_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 37
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 9 - 6.95 um
     ABI_file_id = get_lun()   
     channel_9_filename = channel_1_filename(1:ipos-1) // "_9_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 9 calibration for : ", trim(channel_9_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_9_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_9_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 27
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 10 - 7.34 um
     ABI_file_id = get_lun()   
     channel_10_filename = channel_1_filename(1:ipos-1) // "_10_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 10 calibration for : ", trim(channel_10_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_10_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_10_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 28
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 11 - 8.50 um
     ABI_file_id = get_lun()   
     channel_11_filename = channel_1_filename(1:ipos-1) // "_11_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 11 calibration for : ", trim(channel_11_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_11_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_11_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 29
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 12 - 9.61 um
     ABI_file_id = get_lun()   
     channel_12_filename = channel_1_filename(1:ipos-1) // "_12_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 12 calibration for : ", trim(channel_12_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_12_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_12_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 30
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 13 - 10.35 um
     ABI_file_id = get_lun()   
     channel_13_filename = channel_1_filename(1:ipos-1) // "_13_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 13 calibration for : ", trim(channel_13_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_13_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_13_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 38
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 14 - 11.20 um
     ABI_file_id = get_lun()   
     channel_14_filename = channel_1_filename(1:ipos-1) // "_14_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 14 calibration for : ", trim(channel_14_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_14_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_14_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 31
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 15 - 12.30 um
     ABI_file_id = get_lun()   
     channel_15_filename = channel_1_filename(1:ipos-1) // "_15_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 15 calibration for : ", trim(channel_15_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_15_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_15_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 32
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

     !--- Band 16 - 13.30 um
     ABI_file_id = get_lun()   
     channel_16_filename = channel_1_filename(1:ipos-1) // "_16_" // &
                            channel_1_filename(ipos+3:ilen)
     print*,"Channel 16 calibration for : ", trim(channel_16_filename)
     if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
       call mread_open(trim(Temporary_Data_Dir)//trim(channel_16_filename)//CHAR(0), ABI_file_id)
     else
       call mread_open(trim(Image%Level1b_Path)//trim(channel_16_filename)//CHAR(0), ABI_file_id)
     endif  

     ichan_modis = 33
     call load_ABI_calibration(ABI_file_id, AREAstr,ichan_modis)
     call mread_close(ABI_file_id)

   endif !--- Segment 1

   !--- END CALIBRATION BLOCK READ
   
   !--- START READ OF AREA FILE

   !---   read channel 3 (ABI channel 1)
   if (Sensor%Chan_On_Flag_Default(3) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_1_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_1_filename)
       endif

       !---stw Debug
       !print*,"filename full   : ", trim(Channel_X_Filename_Full)
       !print*,"AREA string     : ", AREAstr
       !print*,"Segment Num     : ", Segment_Number
       !print*,"Lines per seg   : ", Image%Number_Of_Lines_Per_Segment
       !print*,"Lines read this seg : ", Image%Number_Of_Lines_Read_This_Segment
       !stop
       !---stw End Debug

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

       !--- Make reflectance table.
       call ABI_Reflectance(Two_Byte_Temp,ch(3)%Ref_Toa(:,:),1)

   endif

   !---   read channel 1 (ABI channel 2)
   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_2_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_2_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

       !--- Make reflectance table.
       call ABI_Reflectance(Two_Byte_Temp,ch(1)%Ref_Toa(:,:),2)

       !--- store ch2 counts for support of PATMOS-x calibration studies
       Ch1_Counts = Two_Byte_Temp

   endif

   !---   read channel 2 (ABI channel 3)
   if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_3_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_3_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

       !--- Make reflectance table.
       call ABI_Reflectance(Two_Byte_Temp,ch(2)%Ref_Toa(:,:),3)

   endif

   !---   read channel 26 (ABI channel 4)
   if (Sensor%Chan_On_Flag_Default(26) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_4_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_4_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

       !--- Make reflectance table.
       call ABI_Reflectance(Two_Byte_Temp,ch(26)%Ref_Toa(:,:),4)

   endif

   !---   read channel 6 (ABI channel 5)
   if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_5_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_5_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

       !--- Make reflectance table.
       call ABI_Reflectance(Two_Byte_Temp,ch(6)%Ref_Toa(:,:),5)

   endif

   !---   read channel 7 (ABI channel 6)
   if (Sensor%Chan_On_Flag_Default(7) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_6_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_6_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

       !--- Make reflectance table.
       call ABI_Reflectance(Two_Byte_Temp,ch(7)%Ref_Toa(:,:),6)

   endif

   !---   read channel 20 (ABI channel 7)
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_7_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_7_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

       call ABI_RADIANCE_BT(7, Two_Byte_Temp, ch(20)%Rad_Toa, ch(20)%Bt_Toa)

   endif

   !---   read channel 37 (ABI channel 8)
   if (Sensor%Chan_On_Flag_Default(37) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_8_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_8_filename)
       endif

      call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

      call ABI_RADIANCE_BT(8, Two_Byte_Temp, ch(37)%Rad_Toa, ch(37)%Bt_Toa)

   endif

   !---   read channel 27 (ABI channel 9)
   if (Sensor%Chan_On_Flag_Default(27) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_9_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_9_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)
       
      call ABI_RADIANCE_BT(9, Two_Byte_Temp, ch(27)%Rad_Toa, ch(27)%Bt_Toa)

   endif

   !---   read channel 28 (ABI channel 10)
   if (Sensor%Chan_On_Flag_Default(28) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_10_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_10_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

      call ABI_RADIANCE_BT(10, Two_Byte_Temp, ch(28)%Rad_Toa, ch(28)%Bt_Toa)

   endif

   !---   read channel 29 (ABI channel 11)
   if (Sensor%Chan_On_Flag_Default(29) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_11_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_11_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

      call ABI_RADIANCE_BT(11, Two_Byte_Temp, ch(29)%Rad_Toa, ch(29)%Bt_Toa)

   endif

   !---   read channel 30 (ABI channel 12)
   if (Sensor%Chan_On_Flag_Default(30) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_12_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_12_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

      call ABI_RADIANCE_BT(12, Two_Byte_Temp, ch(30)%Rad_Toa, ch(30)%Bt_Toa)

   endif

   !---   read channel 38 (ABI channel 13)
   if (Sensor%Chan_On_Flag_Default(38) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_13_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_13_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

      call ABI_RADIANCE_BT(13, Two_Byte_Temp, ch(38)%Rad_Toa, ch(38)%Bt_Toa)

   endif

   !---   read channel 31 (ABI channel 14)
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_14_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_14_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

      call ABI_RADIANCE_BT(14, Two_Byte_Temp, ch(31)%Rad_Toa, ch(31)%Bt_Toa)

   endif

   !---   read channel 32 (ABI channel 15)
   if (Sensor%Chan_On_Flag_Default(32) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_15_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_15_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

      call ABI_RADIANCE_BT(15, Two_Byte_Temp, ch(32)%Rad_Toa, ch(32)%Bt_Toa)

   endif

   !---   read channel 33 (ABI channel 16)
   if (Sensor%Chan_On_Flag_Default(33) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_16_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_16_filename)
       endif

       call GET_ABI_IMAGE(trim(Channel_X_Filename_Full), &
                                    AREAstr, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment, &
                                    Two_Byte_Temp)

      call ABI_RADIANCE_BT(16, Two_Byte_Temp, ch(33)%Rad_Toa, ch(33)%Bt_Toa)

   endif

   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     Scan_Number(Line_Idx) = first_line_in_segment + Line_Idx
     Scan_Time(Line_Idx) = image_time_ms + (Scan_Number(Line_Idx)-1) * Scan_rate
   enddo

!------------------------------------------------------------------------------
! ABI Angles
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
    
 end subroutine READ_ABI
 
 !-------------------------------------------------------------------------------
 ! LOAD ABI CALIBRATION
 ! 
 !  For single band AREA files, the calibration block is distinct for each
 !  channel.  It must be called for all 16 channels to load the calibration
 !  tables.
 !
 !  The calibration block is an (18,16) array.  Second array index is the band.
 !  The below example will be for band 1.
 !
 !   Channel Number     -> (1,1)
 !   Channel wavelength -> (2,1)
 !   Scale Factor       -> (3,1)
 !   Add Offset         -> (4,1)
 !   Radianc to Albedo  -> (12,1)
 !   FK1                -> (14,1)
 !   FK2                -> (15,1)
 !   BC1                -> (16,1)
 !   BC2                -> (17,1)
 !   Invalid Value      -> (18,1)
 !
 ! Current calibration source type is 'ABIN'
 !
 ! Band Number Filter Map for single band AREA files.
 ! From the AREA directory, index 19.
 !  Band 1  -> 1
 !  Bamd 2  -> 2
 !  Band 3  -> 4
 !  Band 4  -> 8
 !  Band 5  -> 16
 !  Band 6  -> 32
 !  Band 7  -> 64
 !  Band 8  -> 128
 !  Band 9  -> 256
 !  Band 10 -> 512
 !  Band 11 -> 1024
 !  Band 12 -> 2048
 !  Band 13 -> 4096
 !  Band 14 -> 8192
 !  Band 15 -> 16384
 !  Band 16 -> 32768
 !
 ! REFERENCE - McIDAS-X 2016.1 kbxabin.dlm
 !
 !-------------------------------------------------------------------------------
 subroutine load_ABI_calibration(lun, AREAstr, ichan_modis)
   integer(kind=int4), intent(in) :: lun, ichan_modis
   type(AREA_STRUCT), intent(in):: AREAstr
   integer :: i
   integer(kind=int4) :: local_sndr_filter_map
   integer, parameter :: calb_size = 289
   integer(kind=int4), parameter :: MAXVAL = 65536
   integer(kind=int4), dimension(18,16) :: ibuf
   integer(kind=int4) :: band
   integer(kind=int4) :: iCNT
   integer(kind=int4) :: calb_bandNo
   integer(kind=int4) :: calb_errorCount
   integer(kind=int4) :: calb_invalidValue
   real(kind=real8) :: dRAD
   real(kind=real4) :: rRAD
   real(kind=real8) :: dALB
   real(kind=real4) :: rALB
   real(kind=real4) :: rTEMP
   real(kind=real8) :: calb_gainCnt2rad
   real(kind=real8) :: calb_cnstCnt2rad
   real(kind=real8) :: calb_rad2albedo

   !--- Store the correct calibration source type in the AREA structure.  It will be
   !--- important if the calibration type changes.
   call mreadf_int_o(lun,0,4,64,i4buf)
   call move_bytes(4,i4buf(1+51:1+51),AREAstr%src_Type,0)
   local_sndr_filter_map = i4buf(19) ! Need band map for current open AREA file.

   if (AREAstr%src_Type .eq. 'ABIN') then
     print*,"ABIN calibration, proceeding ...", AREAstr%num_chan, local_sndr_filter_map, ichan_modis
   else
     print*,"Unknown ABI calibration type, exiting ..."
     stop
   endif

   !--- Single band images are required.
   if (AREAstr%num_chan > 1) then
     print*,"Multi-band image, exiting ..."
     stop
   endif

   !--- Initialize these variables
   calb_errorCount = -999
   calb_invalidValue = -999

   !--- Extract calibration block per image.
   select case(local_sndr_filter_map)
     !--- Bands 1-6 are albedo/reflectance.
     case(1,2,4,8,16,32)

       !--- Set band number
       if (local_sndr_filter_map == 1)  band=1
       if (local_sndr_filter_map == 2)  band=2
       if (local_sndr_filter_map == 4)  band=3
       if (local_sndr_filter_map == 8)  band=4
       if (local_sndr_filter_map == 16) band=5
       if (local_sndr_filter_map == 32) band=6

       !--- Read the calibration block.
       call mreadf_int_o(lun,AREAstr%cal_offset+4,4,calb_size,ibuf)

       !--- Fill variables with calibration block information.
       calb_bandNo       = ibuf(1,band)
       calb_gainCnt2rad  = DBLE(ibuf(3,band)) * 0.0000001D0
       calb_cnstCnt2rad  = DBLE(ibuf(4,band)) * 0.0000001D0
       calb_rad2albedo   = DBLE(ibuf(12,band)) * 0.0000001D0
       calb_errorCount   = ibuf(18,band)      
       calb_invalidValue = ibuf(18,band)

       !--- Loop through MAXVAL RAW counts to build tables.
       do i = 1, MAXVAL

         !--- Count value is table index.
         iCNT = i - 1

         !--- Convert RAW to RAD
         if (iCNT == calb_errorCount) then
           dRAD = calb_invalidValue
         else
           dRAD = dble(iCNT)*calb_gainCnt2rad+calb_cnstCnt2rad
           if (dRAD <= 0.0D0 ) then
             dRAD = calb_invalidValue
           endif
         endif

         !--- Convert RAD to ALB.
         if (dRAD == calb_invalidValue) then 
           dALB = 0
           rALB = 0.0
         else
           dALB = calb_rad2albedo*dRAD
           rALB = REAL(dALB) * 100.0   ! To match McIDAS.
         endif

         !--- Fill lookup table with albedo.
         Ref_Table(band,i) = rAlb

       enddo

     !--- Bands 7-16 are thermal or mixed.
     case(64,128,256,512,1024,2048,4096,8192,16384,32768)
     
       !--- Set band number
       if (local_sndr_filter_map == 64)    band=7
       if (local_sndr_filter_map == 128)   band=8
       if (local_sndr_filter_map == 256)   band=9
       if (local_sndr_filter_map == 512)   band=10
       if (local_sndr_filter_map == 1024)  band=11
       if (local_sndr_filter_map == 2048)  band=12
       if (local_sndr_filter_map == 4096)  band=13
       if (local_sndr_filter_map == 8192)  band=14
       if (local_sndr_filter_map == 16384) band=15
       if (local_sndr_filter_map == 32768) band=16

       !--- Read the calibration block.
       call mreadf_int_o(lun,AREAstr%cal_offset+4,4,calb_size,ibuf)
    
       !--- Fill variables with calibration block information.
       calb_bandNo       = ibuf(1,band)
       calb_gainCnt2rad  = DBLE(ibuf(3,band)) * 0.0000001D0
       calb_cnstCnt2rad  = DBLE(ibuf(4,band)) * 0.0000001D0

       !--- Loop through MAXVAL RAW counts to build tables.
       do i = 1, MAXVAL

         !--- Count value is table index.
         iCNT = i - 1

         !--- Convert RAW to RAD
         if (iCNT == calb_errorCount) then
           dRAD = calb_invalidValue
         else
           dRAD = dble(iCNT)*calb_gainCnt2rad+calb_cnstCnt2rad
           if (dRAD <= 0.0D0 ) then
             dRAD = calb_invalidValue
           endif
         endif

         !--- Convert RAD to TEMP.
         if (dRAD == calb_invalidValue) then
           dRAD =  calb_invalidValue
           rTEMP = calb_invalidValue
         else
           rRAD = REAL(dRAD)
           rTEMP = PLANCK_TEMP_FAST(ichan_modis,rRAD)
         endif

         !--- Fill look up tables with RAD and TEMP. Scale to integers.
         if (dRAD == calb_invalidValue) then
           bt_table(band,1,i) = NINT(Missing_Value_Real4)
           rad_table(band,1,i) = NINT(Missing_Value_Real4)
         else
           bt_table(band,1,i) = NINT(rTEMP * 100.0)
           rad_table(band,1,i) = NINT(dRAD * 1000.0)
         endif

       enddo

   end select

 end subroutine load_ABI_calibration

 ! Perform ABI Navigation

 subroutine ABI_NAVIGATION(xstart,ystart,xsize,ysize,xstride, &
                            AREAstr,NAVstr_ABI)
    integer(kind=int4) :: xstart, ystart
    integer(kind=int4) :: xsize, ysize
    integer(kind=int4) :: xstride  
    type (AREA_STRUCT) :: AREAstr
    type (GVAR_NAV), intent(in)    :: NAVstr_ABI
    
    integer :: i, j, ii, jj, imode
    real(kind(0.0d0)) :: latitude, longitude
    real(kind=real4) :: height
    integer :: FGF_type = 1 ! GOES-16

    NAVstr_ABI_NAV = NAVstr_ABI
    
    imode = -1
    height = 0.0   !Used for parallax correction
    
    Nav%Lat = Missing_Value_Real4
    Nav%Lon = Missing_Value_Real4

    !---stw Debug
    !print*,"nav type : ", NAVstr_ABI%nav_type
    !print*,"ysize : ", ysize
    !print*,"xsize : ", xsize
    !print*,"ystart : ", ystart
    !print*,"north_bound : ", AREAstr%north_bound
    !print*,"west vis pixel : ", AREAstr%west_vis_pixel
    !print*,"line_res : ", AREAstr%line_res
    !print*,"AREA String : ", AREAstr
    !print*,"CFAC : ", NAVstr_ABI%abiCFAC
    !print*,"COFF : ", NAVstr_ABI%abiCOFF
    !print*,"LFAC : ", NAVstr_ABI%abiLFAC
    !print*,"LOFF : ", NAVstr_ABI%abiLOFF
    !print*,"sub lon : ", NAVstr_ABI%sub_lon
    !print*,"xstride : ", xstride
    !stop
    !---stw End Debug
    
    if (NAVstr_ABI%nav_type == 'ABIN') then      
                
      do j=1, ysize
            
        !--- For 2 km full disk, we need navigation transformations in 1 km
        !--- space, as Band 1 is the active AREA file at this point.  However,
        !--- lcor/ecor are in 1/2 km space.
        !---stw ORIGINAL jj = ystart + (j-1)*AREAstr%line_res
        jj = ystart*AREAstr%line_res + (j-1)*AREAstr%line_res
        jj = ((jj + AREAstr%north_bound + (NAVstr_ABI%BRES-1)))/NAVstr_ABI%BRES
            
        !---stw Debug
        !print*,"north_bound : ", AREAstr%north_bound
        !print*,"ystart : ", ystart
        !print*,"j : ", j
        !print*,"BRES : ", NAVstr_ABI%BRES
        !print*,"jj : ", jj
        !print*,"xsize : ", xsize
        !print*,"elem res : ", AREAstr%elem_res
        !stop
        !---stw End Debug
            
        do i=1, xsize

          !ii = (xstart + (i-1)*xstride + AREAstr%west_vis_pixel + (NAVstr_ABI%BRES-1))/NAVstr_ABI%BRES
          ii = xstart + (i-1)*xstride*AREAstr%elem_res
          ii = ((ii + AREAstr%west_vis_pixel + (NAVstr_ABI%BRES-1)))/NAVstr_ABI%BRES

          !---stw Debug
          !print*,"xstart : ", xstart
          !print*,"xstride : ", xstride
          !print*,"i : ", i
          !print*,"ii : ", ii
          !print*,"BRES : ", NAVstr_ABI%BRES
          !print*,"west_vis_pixel : ", AREAstr%west_vis_pixel
          !print*,"elem res : ", AREAstr%elem_res
          !stop
          !---stw End Debug

          CALL fgf_to_earth(FGF_type,                  &
                                  DBLE(ii),                  &
                                  DBLE(jj),                  &
                                  DBLE(NAVstr_ABI%abiCFAC),    &
                                  DBLE(NAVstr_ABI%abiCOFF),    &
                                  DBLE(NAVstr_ABI%abiLFAC),    &
                                  DBLE(NAVstr_ABI%abiLOFF),    &
                                  DBLE(NAVstr_ABI%sub_lon), &
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
      
 end subroutine ABI_NAVIGATION
 
!------------------------------------------------------------------
! subroutine to convert ABI counts to radiance and brightness
! temperature
!------------------------------------------------------------------
  
  subroutine ABI_RADIANCE_BT(chan_num, ABI_Counts, rad2, temp1)

    integer (kind=INT2), dimension(:,:), intent(in):: ABI_Counts
    !integer (kind=int2), INTENT(in) :: chan_num
    integer (kind=int4), INTENT(in) :: chan_num
    real (kind=real4), dimension(:,:), INTENT(out):: temp1, rad2
    
    integer :: i, j, index

    do j = 1, Image%Number_Of_Lines_Read_This_Segment
      do i = 1, Image%Number_Of_Elements
        if (Space_Mask(i,j) == sym%NO_SPACE) then
          index = int(ABI_Counts(i,j),kind=int2)
          rad2(i,j) = real(rad_table(chan_num,1,index+1),kind=real4)/1000.0
          temp1(i,j) = real(bt_table(chan_num,1,index+1),kind=real4)/100.0                    
        else
          rad2(i,j) = Missing_Value_Real4
          temp1(i,j) = Missing_Value_Real4
        endif
      end do
    end do    
  
  end subroutine ABI_RADIANCE_BT
  
!------------------- ABI NAV BLOC ------------------------------
!--- From McIDAS-X 2016.2 nvxabin.dlm
!--- Navigation Block:
!---   1: ABIN -> Navigation module string name.
!---   2: LOFF
!---   3: COFF
!---   4: LFAC
!---   5: CFAC
!---   6: Satellite SubPoint Longitude
!---   7: Base Image Resolution.
!---
!--- The call from sensor_module.f90 is for Band 1.  The data is
!--- at 2km. The navigation module will take care of the necessary
!--- transformations to account for 2km data. 
!---  LOFF = 15185800
!---  COFF = -15185800
!---  LFAC = -2800
!---  CFAC = 2800
!-------------------------------------------------------------

 subroutine READ_NAVIGATION_BLOCK_ABI(filename, AREAstr, NAVstr)
  CHARACTER(len=*), intent(in):: filename
  type(AREA_STRUCT), intent(in):: AREAstr
  type(GVAR_NAV), intent(inout):: NAVstr
 
  integer :: geos_nav
  integer(kind=int4)nav_offset
  integer:: number_of_words_read
  integer(kind=int4), dimension(640) :: i4buf

  nav_offset = AREAstr%sec_key_nav

  ! ABI navigation is names ABIN, but will use the GEOS transforms
  ! from geos_transform_pix.c

  ! Read the navigation block into the buffer.
  geos_nav = sym%NO
  call mreadf_int(trim(filename)//CHAR(0),nav_offset,4,640,&
                    number_of_words_read, i4buf)
  call move_bytes(4,i4buf(1),NAVstr%nav_type,0)

  ! LOFF, COFF, CFAC, LFAC stored in McIDAS header for 2km data,
  ! assuming that we are basing the navigation off of Band 1.
  NAVstr%abiLOFF=real(i4buf(2),kind=real8) / 100000000.0
  NAVstr%abiCOFF=real(i4buf(3),kind=real8) / 100000000.0
  NAVstr%abiLFAC=real(i4buf(4),kind=real8) / 100000000.0
  NAVstr%abiCFAC=real(i4buf(5),kind=real8) / 100000000.0
  NAVstr%sub_lon = real(i4buf(6),kind=real4) / 10       ! SUBLON stored as SUBLON *10 in McIDAS NAV block
  NAVstr%sublon = NAVstr%sub_lon
  NAVstr%BRES=i4buf(7)

  !---stw Debug
  !print*,"a LOFF : ", NAVstr%abiLOFF
  !print*,"a COFF : ", NAVstr%abiCOFF
  !print*,"a LFAC : ", NAVstr%abiLFAC
  !print*,"a CFAC : ", NAVstr%abiCFAC
  !print*,"a sub_lon : ", NAVstr%sublon
  !stop
  !---stw End Debug

 end subroutine READ_NAVIGATION_BLOCK_ABI

 ! Perform ABI Reflectance calculation.  These vaules have been
 ! stored in an array, previous to this call.
 subroutine ABI_Reflectance(ABI_Counts, Alb_Temp,chan_num)

   integer (kind=INT2), dimension(:,:), intent(in):: ABI_Counts
   !integer (kind=INT2), intent(in):: chan_num
   integer (kind=INT4), intent(in):: chan_num
   real (kind=real4), dimension(:,:), intent(out):: Alb_Temp

   integer :: i, j
   integer :: index

   !---stw Debug
   !print*,"Channel Number : ", chan_num
   !stop
   !---stw End Debug

   DO j=1, Image%Number_Of_Lines_Read_This_Segment
     DO i=1, Image%Number_Of_Elements

       IF (space_mask(i,j) == sym%NO_SPACE .and. Geo%solzen(i,j) < 90.0) THEN

         index = int(ABI_Counts(i,j),kind=int2)

         alb_temp(i,j) = Missing_Value_Real4

         IF((index .GT. 0) .AND. (index .LE. Ntable_ABI)) THEN
           Alb_Temp(i,j) = (real(Ref_Table(chan_num,index+1),kind=real4))
         ENDIF

       ELSE
         Alb_Temp(i,j) = Missing_Value_Real4
       ENDIF

     END DO

   END DO

 end subroutine ABI_Reflectance

 subroutine GET_ABI_IMAGE(filename,AREAstr, &
                                    segment_number, &
                                    num_lines_per_segment, &
                                    num_lines_read, image)

   character(len=*), intent(in):: filename
   type (AREA_STRUCT), intent(in) :: AREAstr
   integer(kind=int4), intent(in):: segment_number
   integer(kind=int4), intent(in):: num_lines_per_segment
   integer(kind=int4), intent(out):: num_lines_read
   integer(kind=int2), dimension(:,:), intent(out):: image

   integer(kind=int4):: bytes_per_line
   integer(kind=int4):: words_in_prefix
   integer(kind=int4):: words_per_line
   integer(kind=int4):: first_line_in_segment
   integer(kind=int4):: last_line_in_segment
   integer(kind=int4):: first_byte_in_segment
   integer(kind=int4):: number_of_words_in_segment
   integer(kind=int4):: bytes_per_pixel
   integer(kind=int4):: num_byte_ln_prefix
   integer(kind=int4):: pri_key_nav
   integer(kind=int4):: number_of_words_read
   integer(kind=int4):: word_start
   integer(kind=int4):: dummy
   integer(kind=int4):: word_end
   integer(kind=int4):: bytemove
   integer(kind=int4):: Line_Idx
   integer(kind=int2), dimension(:), allocatable:: word_buffer
   integer (kind=int1), dimension(:), allocatable :: buffer1
   integer(kind=int4), dimension(64) :: i4buf_temp
   integer:: nwords

   bytemove = 0 !Not sure what this needs to be for ABI?
   image = 0

   ! get number of bytes per pixel and num bytes per line for current file
   ! this is needed because the 0.64 and other channels have different values
   ! in the COMS HIRID format.

   call mreadf_int(trim(filename)//CHAR(0),0,4,64,dummy,i4buf_temp)
   bytes_per_pixel = i4buf_temp(11)
   num_byte_ln_prefix = i4buf_temp(15)
   pri_key_nav = AREAstr%pri_key_nav
   bytes_per_line = num_byte_ln_prefix + (AREAstr%num_elem*bytes_per_pixel)
   words_in_prefix = num_byte_ln_prefix / bytes_per_pixel
   words_per_line = words_in_prefix + AREAstr%num_elem
   first_line_in_segment = (segment_number-1)*num_lines_per_segment + 1
   last_line_in_segment = min(AREAstr%num_line,segment_number*num_lines_per_segment)
   first_byte_in_segment = pri_key_nav + &
                         bytes_per_pixel*(first_line_in_segment-1) * words_per_line + &
                         bytes_per_pixel
   number_of_words_in_segment = words_per_line * num_lines_per_segment

   !---stw Debug
   !print*,"bytes per pixel : ", bytes_per_pixel
   !print*,"bytes in prefix : ", num_byte_ln_prefix
   !print*,"pri key nav : ", pri_key_nav
   !print*,"bytes per line : ", bytes_per_line
   !print*,"words in prefix : ", words_in_prefix
   !print*,"words per line : ", words_per_line
   !print*,"1st line in seg : ", first_line_in_segment
   !print*,"last line in seg : ", last_line_in_segment
   !print*,"1st byte in seg : ", first_byte_in_segment
   !print*,"num words in seg : ", number_of_words_in_segment
   !stop
   !---stw End Debug

   allocate(word_buffer(number_of_words_in_segment))

   if ( bytes_per_pixel == 2) then
     call mreadf_int(filename//CHAR(0), &
                 first_byte_in_segment,  &
                 bytes_per_pixel, &
                 number_of_words_in_segment, &
                 number_of_words_read,word_buffer)
   endif

   !---stw Debug
   !print*,"two_byte_temp : ", word_buffer
   !---stw End Debug

   if( bytes_per_pixel == 1) then
     allocate(buffer1(number_of_words_in_segment))

     call mreadf_int(filename//CHAR(0), &
                 first_byte_in_segment,  &
                 bytes_per_pixel, &
                 number_of_words_in_segment, &
                 number_of_words_read,buffer1)

     word_buffer = int(buffer1,kind=int2)
     !Since 1 byte values are always signed in Fortran, I need to add 256 to the
     !negative values
     where (word_buffer < 0)
       word_buffer = word_buffer + 256
     end where
   endif

   !--- update number of scans read
   num_lines_read = number_of_words_read / words_per_line

   do Line_Idx = 1, num_lines_read
     word_start = (words_in_prefix + AREAstr%num_elem)*(Line_Idx-1) + words_in_prefix + 1
     word_end = min(word_start + AREAstr%num_elem - 1,number_of_words_in_segment)
     nwords = int(word_end - word_start)/ABI_Xstride + 1
     image(1:nwords,Line_Idx) = ishft(word_buffer(word_start:word_end:1),bytemove)
   enddo

   deallocate(word_buffer)
   if (allocated(buffer1)) deallocate(buffer1)

 end subroutine GET_ABI_IMAGE

end module ABI_MODULE
