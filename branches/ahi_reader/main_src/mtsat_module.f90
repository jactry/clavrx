!// $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: mtsat_module.f90 (src)
!       MTSAT_module (program)
!
! PURPOSE: This module contains all the subroutines needed to perform navigation
!          and calibration for MTSAT, bot HiRID and HRIT
!
! DESCRIPTION: 
!
! AUTHORS:
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
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
module MTSAT_MODULE

use CONSTANTS
use PIXEL_COMMON
use CALIBRATION_CONSTANTS
use PLANCK
use NUMERICAL_ROUTINES
use CGMS_NAV
use GOES_MODULE
use FILE_UTILITY
use VIEWING_GEOMETRY_MODULE

implicit none
public :: READ_MTSAT
public :: READ_NAVIGATION_BLOCK_MTSAT_FY
public :: CALIBRATE_MTSAT_DARK_COMPOSITE
public :: READ_MTSAT_INSTR_CONSTANTS
         
private :: MTSAT_RADIANCE_BT
private :: MTSAT_REFLECTANCE_PRELAUNCH
private :: MTSAT_REFLECTANCE_GSICS
private :: MTSAT_NAVIGATION
private :: MGIVSR 
 

 type (GVAR_NAV), PRIVATE    :: NAVstr_MTSAT_NAV
 integer, PARAMETER, PRIVATE :: nchan_mtsat= 5
 integer, PARAMETER, PRIVATE :: ndet_mtsat = 4
 integer, PARAMETER, PRIVATE :: ntable_mtsat = 1024

 integer, PRIVATE :: nref_table_mtsat
 integer, PRIVATE :: nbt_table_mtsat
 character(len=4), SAVE, PRIVATE:: calib_type

 integer (kind=int4), dimension(nchan_mtsat,ndet_mtsat,ntable_mtsat), PRIVATE  :: ref_table
 integer (kind=int4), dimension(nchan_mtsat,ndet_mtsat,ntable_mtsat), PRIVATE  :: bt_table
 integer (kind=int4), dimension(nchan_mtsat,ndet_mtsat,ntable_mtsat), PRIVATE  :: rad_table

 integer(kind=int4), public, parameter:: Mtsat_Xstride = 1
 integer(kind=int4), private, parameter:: Num_4km_Scans_Fd = 3712
 integer(kind=int4), private, parameter:: Num_4km_elem_fd = 3712
 integer(kind=int4), private, parameter:: Time_For_Fd_Scan =  1560000 !milliseconds (26min)
 real, private, save:: Scan_rate    !scan rate in millsec / line
 integer(kind=int4), private, parameter:: Mtsat_Byte_Shift = 0

 contains

!----------------------------------------------------------------
! read the MTSAT constants into memory
!-----------------------------------------------------------------
subroutine READ_MTSAT_INSTR_CONSTANTS(Instr_Const_file)
 character(len=*), intent(in):: Instr_Const_file
 integer:: ios0, erstat
 integer:: Instr_Const_lun

 Instr_Const_lun = GET_LUN()

 open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

 print *, "opening ", trim(Instr_Const_file)
 erstat = 0
 if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, "Error opening MTSAT constants file, ios0 = ", ios0
    stop 19
 endif
  read(unit=Instr_Const_lun,fmt="(a3)") sat_name
  read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
  read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
  read(unit=Instr_Const_lun,fmt=*) planck_a1(20), planck_a2(20),planck_nu(20)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(27), planck_a2(27),planck_nu(27)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(31), planck_a2(31),planck_nu(31)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(32), planck_a2(32),planck_nu(32)

  read(unit=Instr_Const_lun,fmt=*) Ch1_Dark_Count
  read(unit=Instr_Const_lun,fmt=*) Ch1_Gain_Low_0,Ch1_Degrad_Low_1, Ch1_Degrad_Low_2
  read(unit=Instr_Const_lun,fmt=*) Launch_Date

  close(unit=Instr_Const_lun)

  !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / ew_Ch20

end subroutine READ_MTSAT_INSTR_CONSTANTS

! Perform MTSAT Reflectance and BT calibration
subroutine READ_MTSAT(segment_number,Channel_1_Filename, &
                     jday, image_time_ms, Time_Since_Launch, &
                     AREAstr,NAVstr_MTSAT)

   integer(kind=int4), intent(in):: segment_number
   character(len=*), intent(in):: Channel_1_Filename
   TYPE (AREA_STRUCT), intent(in) :: AREAstr
   TYPE (GVAR_NAV), intent(in)    :: NAVstr_MTSAT
   integer(kind=int2), intent(in):: jday
   integer(kind=int4), intent(in):: image_time_ms
   real(kind=real4), intent(in):: Time_Since_Launch

   character(len=120):: Channel_X_Filename
   character(len=120):: Channel_X_Filename_Full
   character(len=120):: Channel_X_Filename_Full_uncompressed
   character(len=180):: System_String
   integer:: ipos
   integer:: ilen
   integer:: ielem
   integer:: iline
   integer:: Chan_Idx_Mtsat
   integer:: Chan_Idx_Modis
   integer:: mtsat_file_id
   real(kind=real4):: image_time_hours
   integer(kind=int4):: image_jday
   integer(kind=int4):: first_line_in_segment
   character(len=2):: Chan_Idx_Mtsat_String
   integer:: Num_Elements_This_Image
   integer:: Num_Scans_This_Image
   integer:: Line_Idx

   

   !--- assume Channel_1_file name has a unique "_1_" in the name. 
   !--- determine indices needed to replace that string
   ipos = index(Channel_1_Filename, "_1_")
   ilen = len(Channel_1_Filename)
    
   first_line_in_segment = (segment_number-1)*Image%Number_Of_Lines_Per_Segment

   !---------------------------------------------------------------------------
   ! MTSAT Navigation (Do Navigation and Solar angles first)
   !---------------------------------------------------------------------------
   
   call mtsat_navigation(1,first_line_in_segment,&
                              Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment,1,&
                              AREAstr,NAVstr_MTSAT)
   
   if (segment_number == 1) then

      image_jday = jday
      image_time_hours = image_time_ms / 60.0 / 60.0 / 1000.0

      !--- compute scan rate for future use
      Num_Elements_This_Image =  int(AREAstr%num_elem / MTSAT_Xstride) + 1
      Num_Scans_This_Image = AREAstr%num_line
      Scan_Rate = real((Num_Elements_This_Image)/               &
                  real(Num_4km_elem_fd/MTSAT_Xstride)) * &
                  real((Num_Scans_This_Image) / real(Num_4km_Scans_Fd)) * &
                  real(Time_For_Fd_Scan) / real(Num_Scans_This_Image)

       do Chan_Idx_Mtsat = 2,5

       if (Chan_Idx_Mtsat == 5) Chan_Idx_Modis = 20
       if (Chan_Idx_Mtsat == 4) Chan_Idx_Modis = 27
       if (Chan_Idx_Mtsat == 2) Chan_Idx_Modis = 31
       if (Chan_Idx_Mtsat == 3) Chan_Idx_Modis = 32
       
       write(Chan_Idx_Mtsat_String,fmt="(I1.1)") Chan_Idx_Mtsat
       if(Chan_Idx_Mtsat > 9) write(Chan_Idx_Mtsat_String,fmt="(I2.2)") Chan_Idx_Mtsat

       if (Sensor%Chan_On_Flag_Default(Chan_Idx_Modis) == sym%YES) then

          Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_"//trim(Chan_Idx_Mtsat_String)//"_" // &
                            Channel_1_Filename(ipos+3:ilen)
          
          if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
          else
               Channel_X_Filename_Full = trim(Image%Level1b_Path)//trim(Channel_X_Filename)
          endif

          Channel_X_Filename_Full_uncompressed = trim(Image%Level1b_Path)//trim(Channel_X_Filename)
          if (l1b_gzip == sym%YES) then
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

    ! On first segment, reflectance, BT and rad tables
    ! On first segment, get slope/offset information from McIDAS Header
    mtsat_file_id = get_lun()   
    if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
        call mread_open(trim(Temporary_Data_Dir)//trim(Channel_1_Filename)//CHAR(0), mtsat_file_id)
    else
      call mread_open(trim(Image%Level1b_Path)//trim(Channel_1_Filename)//CHAR(0), mtsat_file_id)
    endif  

    call load_mtsat_calibration(mtsat_file_id, AREAstr)
    call mread_close(mtsat_file_id)


   endif
   
   
   !---   read channel 1 (MTSAT channel 1)
   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_1_Filename)
       else
               Channel_X_Filename_Full = trim(Image%Level1b_Path)//trim(Channel_1_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Mtsat_Byte_Shift, &
                                    AREAstr, Mtsat_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)

        call MTSAT_REFLECTANCE_GSICS(Two_Byte_Temp,Time_Since_Launch,ch(1)%Ref_Toa(:,:))
! old   call MTSAT_REFLECTANCE_GSICS(Two_Byte_Temp,Time_Since_Launch,Ref_Ch1(:,:))

        Ch1_Counts = Two_Byte_Temp

   endif
   
   !---   read channel 20 (MTSAT channel 5)
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_5_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Image%Level1b_Path)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Mtsat_Byte_Shift, &
                                    AREAstr, Mtsat_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)

      call MTSAT_RADIANCE_BT(5_int1, Two_Byte_Temp, ch(20)%Rad_Toa, ch(20)%Bt_Toa)
!cspp call MTSAT_RADIANCE_BT(5_int1, Two_Byte_Temp, Rad_Ch20, Bt_Ch20)

   endif
                    
   
   !---   read channel 27 (MTSAT channel 4)
   if (Sensor%Chan_On_Flag_Default(27) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_4_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Image%Level1b_Path)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Mtsat_Byte_Shift, &
                                    AREAstr, Mtsat_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)

      call MTSAT_RADIANCE_BT(4_int1, Two_Byte_Temp, ch(27)%Rad_Toa, ch(27)%Bt_Toa)
!cspp call MTSAT_RADIANCE_BT(4_int1, Two_Byte_Temp, Rad_Ch27, Bt_Ch27)

   endif


   
   !---   read channel 31 (MTSAT channel 2)
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_2_" // &
                            Channel_1_Filename(ipos+3:ilen)
       
       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Image%Level1b_Path)//trim(Channel_X_Filename)
       endif
       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Mtsat_Byte_Shift, &
                                    AREAstr, Mtsat_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)

      call MTSAT_RADIANCE_BT(2_int1, Two_Byte_Temp, ch(31)%Rad_Toa, ch(31)%Bt_Toa)
!cspp call MTSAT_RADIANCE_BT(2_int1, Two_Byte_Temp, Rad_Ch31, Bt_Ch31)

   endif
   
   
   !---   read channel 32 (MTSAT channel 3)
   if (Sensor%Chan_On_Flag_Default(32) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_3_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Image%Level1b_Path)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Mtsat_Byte_Shift, &
                                    AREAstr, Mtsat_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)

      call MTSAT_RADIANCE_BT(3_int1, Two_Byte_Temp, ch(32)%Rad_Toa, ch(32)%Bt_Toa)
!cspp call MTSAT_RADIANCE_BT(3_int1, Two_Byte_Temp, Rad_Ch32, Bt_Ch32)

   endif
    
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     Scan_Number(Line_Idx) = first_line_in_segment + Line_Idx
     Scan_Time(Line_Idx) = image_time_ms + (Scan_Number(Line_Idx)-1) * Scan_rate
   enddo

!------------------------------------------------------------------------------
! MTSAT Angles
! NOTE: These were private routines in the GOES module. Suggest they become
!       public with different names, since they are used cross platform  
!------------------------------------------------------------------------------
   image_jday = jday
   image_time_hours = image_time_ms / 60.0 / 60.0 / 1000.0
   do iline = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     do ielem = 1,Image%Number_Of_Elements
        call POSSOL(image_jday,image_time_hours, &
                    Nav%Lon_1b(ielem,iline),Nav%Lat_1b(ielem,iline), &
                    Geo%Solzen(ielem,iline),Geo%Solaz(ielem,iline))
     enddo
     call COMPUTE_SATELLITE_ANGLES(Sensor%Geo_Sub_Satellite_Longitude,  &
                                   Sensor%Geo_Sub_Satellite_Latitude, iline)                      
   enddo

   !--- ascending node
   ielem = Image%Number_Of_Elements/2
   do iline = Line_Idx_Min_Segment+1, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     Nav%Ascend(iline) = 0
     if (Nav%Lat_1b(ielem,iline) < Nav%Lat_1b(ielem,iline-1)) then
       Nav%Ascend(iline) = 1
     endif
   enddo
   Nav%Ascend(Line_Idx_Min_Segment) = Nav%Ascend(Line_Idx_Min_Segment+1)

end subroutine READ_MTSAT
 
subroutine LOAD_MTSAT_CALIBRATION(lun, AREAstr)
  integer(kind=int4), intent(in) :: lun
  type(AREA_STRUCT), intent(in):: AREAstr
  integer(kind=int4), dimension(6528) :: ibuf
  character(len=25) :: cbuf
  integer :: nref, nbt, i, j, offset
  integer(kind=int4) :: band_offset_2, band_offset_14, band_offset_15, &
                        band_offset_9, band_offset_7, dir_offset
  real(kind=real4) :: albedo, temperature, radiance
  real(kind=real4), dimension(5)  :: a_mtsat, b_mtsat, nu_mtsat

  call mreadf_int_o(lun,AREAstr%cal_offset,4,6528,ibuf)
  !if (AREAstr%swap_bytes > 0) call swap_bytes4(ibuf,6528)
  
  dir_offset = ibuf(4)
  band_offset_2 = ibuf(6)
  band_offset_14 = ibuf(8)
  band_offset_15 = ibuf(10)
  band_offset_9 = ibuf(12)
  band_offset_7 = ibuf(14)  
  
  nref = 256
  nbt = 1024
  nref_table_mtsat = nref
  nbt_table_mtsat = nbt
  
  ! We need to get the vis calibration type. This is in 5th
  ! 4 character long block of the cal block (as per line 284, kbxmtst.dlm, v1.5)
  ! WCS3 (3/29/2010)
  
  call mreadf_int_o(lun,AREAstr%cal_offset,1,25,cbuf)
  calib_type = cbuf(17:20)

  !---------------------------------------------------------------------
  ! Load the visible channel calibration table for HiRID/MVIS calibration -WCS3
  !---------------------------------------------------------------------
  
  do j = 1,4
    do i=1,256,4
      offset = band_offset_2/4 + (j-1) * 64
      albedo = real(ibuf(offset+i/4+1),kind=real4) / 10000.
      
      ref_table(2,j,i)   = nint(albedo * 100.) 
      ref_table(2,j,i+1) = nint(albedo * 100.)
      ref_table(2,j,i+2) = nint(albedo * 100.)
      ref_table(2,j,i+3) = nint(albedo * 100.)
    end do
  end do

  !---------------------------------------------------------------------
  ! Load the IR channel calibration table.
  !---------------------------------------------------------------------
  
  !--- first set a, b, nu in form that is needed for table. Outside this routine, these
  !    coefficents are not needed.
  
  if (AREAstr%sat_id_num  == 84) then
  ! Updated with values from GSICS activities at JMA http://mscweb.kishou.go.jp/monitoring/gsics/ir/techinfo.htm
        nu_mtsat(5) = 2652.9316   !  3.9 micron
        nu_mtsat(4) = 1482.2068   !  6.8 micron
        nu_mtsat(2) = 926.6118   ! 10.8 micron
        nu_mtsat(3) = 833.1675   ! 12.0 micron
  
        a_mtsat(5) = 0.9969755    !  3.9 micron
        a_mtsat(4) = 0.9991187    !  6.8 micron
        a_mtsat(2) = 0.9987587   ! 10.8 micron
        a_mtsat(3) = 0.9992525   ! 12.0 micron
  
        b_mtsat(5) =  2.3473427   !  3.9 micron
        b_mtsat(4) =  0.3785336   !  6.8 micron
        b_mtsat(2) = 0.3592380   ! 10.8 micron
        b_mtsat(3) = 0.1968675   ! 12.0 micron
  
  
  endif
  


  if (AREAstr%sat_id_num  == 85) then
  ! Coeff provided by JMA to WCS3. Will be on GSICS/JMA website at some point
  ! http://mscweb.kishou.go.jp/monitoring/gsics/ir/techinfo.htm

        nu_mtsat(5) = 2684.1181    !  3.9 micron
        nu_mtsat(4) = 1476.6898    !  6.8 micron
        nu_mtsat(2) =  926.4627   ! 10.8 micron
        nu_mtsat(3) =  835.6672   ! 12.0 micron
  
        a_mtsat(5) = 0.9967825    !  3.9 micron
        a_mtsat(4) = 0.9991492    !  6.8 micron
        a_mtsat(2) = 0.9991676   ! 10.8 micron
        a_mtsat(3) = 0.9987568   ! 12.0 micron
  
        b_mtsat(5) = 2.4635230    !  3.9 micron
        b_mtsat(4) = 0.3645235    !  6.8 micron
        b_mtsat(2) = 0.3597581   ! 10.8 micron
        b_mtsat(3) = 0.2195110   ! 12.0 micron

  endif
  
  do i=1, 1024
    !  3.9 micron
    offset = band_offset_7/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(5,1,i) = nint(temperature * 100.)
    radiance = c1*(nu_mtsat(5)**3)/(exp((c2*nu_mtsat(5))/ &
               (a_mtsat(5)*temperature+b_mtsat(5)))-1.0)
    rad_table(5,1,i) = nint(radiance * 1000.)
    
    !  6.8 micron
    offset = band_offset_9/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(4,1,i) = nint(temperature * 100.)
    radiance = c1*(nu_mtsat(4)**3)/(exp((c2*nu_mtsat(4))/ &
               (a_mtsat(4)*temperature+b_mtsat(4)))-1.0)
    rad_table(4,1,i) = nint(radiance * 1000.)
    
    ! 10.8 micron
    offset = band_offset_14/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(2,1,i) = nint(temperature * 100.)
    radiance = c1*(nu_mtsat(2)**3)/(exp((c2*nu_mtsat(2))/ &
               (a_mtsat(2)*temperature+b_mtsat(2)))-1.0)
    rad_table(2,1,i) = nint(radiance * 1000.)
        
    ! 12.0 micron
    offset = band_offset_15/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(3,1,i) = nint(temperature * 100.)
    radiance = c1*(nu_mtsat(3)**3)/(exp((c2*nu_mtsat(3))/ &
               (a_mtsat(3)*temperature+b_mtsat(3)))-1.0)
    rad_table(3,1,i) = nint(radiance * 1000.) 
  end do
 
end subroutine LOAD_MTSAT_CALIBRATION

!----------------------------------------------------------------------
! Perform MTSAT Reflectance calculation using fits to JMA GSICS tables
!----------------------------------------------------------------------
subroutine MTSAT_REFLECTANCE_GSICS(Mtsat_Counts, Time_Temp_Since_Launch, Alb_Temp)
    integer (kind=INT2), dimension(:,:), intent(in):: Mtsat_Counts
    real (kind=real4), intent(in):: Time_Temp_Since_Launch
    real (kind=real4), dimension(:,:), intent(out):: Alb_Temp

    integer :: index
    integer:: i, j

    Ch1_Gain_Low = Ch1_Gain_Low_0*(100.0+Ch1_Degrad_Low_1*Time_Temp_Since_Launch + &
                                  Ch1_Degrad_Low_2*Time_Temp_Since_Launch**2)/100.0

    Alb_Temp = Ch1_Gain_Low * ( Mtsat_Counts - Ch1_Dark_Count)
    where (Space_Mask == sym%YES)
     Alb_Temp = Missing_Value_Real4
    endwhere
    
!   do j = 1,Image%Number_Of_Lines_Read_This_Segment
!     do i = 1,Image%Number_Of_Elements
!       if (Space_Mask(i,j) == sym%NO) then
!          Alb_Temp(i,j) = Ch1_Gain_Low * ( Mtsat_Counts(i,j) - Ch1_Dark_Count)
!       endif
!     enddo
!   enddo

end subroutine MTSAT_REFLECTANCE_GSICS

!----------------------------------------------------------------------
! Perform MTSAT Reflectance calculation using fits to prelaunch tables
!----------------------------------------------------------------------
subroutine MTSAT_REFLECTANCE_PRELAUNCH(Mtsat_Counts, alb_temp)
                              
    integer (kind=INT2), dimension(:,:), intent(in):: Mtsat_Counts
    real (kind=real4), dimension(:,:), intent(out):: alb_temp

    integer :: index
    integer:: i, j
    
    
    !---- WCS3 ------!
    ! Reflectance linear fit of reflectance table for HRIT files is taken from JMA Excel table
    ! of calibration coefficents.
    ! http://mscweb.kishou.go.jp/operation/calibration/mt1r/HRIT/mtsat1r_hrit_calibration_table.xls
    !Currently the same for MTSAT-1r and MTSAT-2. McIDAS currently doesn't calculate this.
    !---- WCS3 ------!

    do j=1, Image%Number_Of_Lines_Read_This_Segment
      do i=1, Image%Number_Of_Elements
        if (Space_Mask(i,j) == sym%NO) then
          !  Because MTSAT has two vis calibration type, we need to 
          !  which route it needs to go. calibration.f90 was edited to 
          !  gather the vis calibration type, which is stored in thesat_info
          !  structure. Now we USE it here, based off of how to do it via
          !  kbxmtst.dlm (v1.5) from McIDAS - WCS3 - 3/29/2010
                    
          if ( calib_type == 'MVSH' .OR. &
               calib_type == 'HSVM') then

            !---- WCS3 ------!
            ! JMA has a simple conversion table (for reflectance
            ! factor, so 0 to 1) at http://mscweb.kishou.go.jp/operation/calibration/mt1r/HRIT/mt1r_hrit.htm
            ! So you need to multiply by 100. In addition, JMA doesn't divide by solzen or SED.
            
            ! The digital count goes from 0 - 1023, but linear fit assumes indicies go from 1
            ! to 1024. So, you need to add 1 to get the right index. See above website for reference
            !---- WCS3 ------!
            
            index = int(Mtsat_Counts(i,j),kind=int2) + 1.0
                        !---- WCS3 ------!
            !linear fit based off of JMA conversion table
            !y = 0.000978494701531316*x - 0.00197851402698724
            !y = reflectance factor
            !x = index (1 to 1024)

            !MTSAT-2 = = 0.000978494701531316*L14 - 0.00197851402698724
            !---- WCS3 ------!
                               
!AKH - AIT Standards say these numbers should be in variables or parameters

            Alb_Temp(i,j) = 0.000978494701531316*(index) - 0.00197851402698724
                        
            !Recall that reflectance factor goes from 0 to 1 and is not corrected by
            ! cossolzen and SED. Correction with SED and Cossolzen done elsewhere - WCS3
            
            Alb_Temp(i,j) = (Alb_Temp(i,j) * 100.0) 
                              
          endif
          
          if (calib_type == 'MVIS' .OR. &
               calib_type == 'SIVM') then
          
            !Not sure if I need to add 1 here
            index = int(Mtsat_Counts(i,j),kind=int2)

            if((index .GT. 0) .AND. (index .LE. 1024)) then
            
                alb_temp(i,j) = &
             (real(ref_table(2,1,index),kind=real4) /  100.0) 
     
            endif                        
            
          endif
        else

            alb_temp(i,j) = Missing_Value_Real4

        endif      
      end DO
    end DO
        
end subroutine MTSAT_REFLECTANCE_PRELAUNCH

 ! Perform MTSAT Navigation

 subroutine mtsat_navigation(xstart,ystart,xsize,ysize,xstride, &
                            AREAstr,NAVstr_MTSAT)
    integer(kind=int4) :: xstart, ystart
    integer(kind=int4) :: xsize, ysize
    integer(kind=int4) :: xstride  
    type (AREA_STRUCT) :: AREAstr
    TYPE (GVAR_NAV), intent(in)    :: NAVstr_MTSAT
    
    integer :: i, j, ii, jj, ierr, imode
    real(kind(0.0d0)) :: latitude, longitude
    real(kind=real4) :: elem, line, height
    real(kind=real4) :: dlon, dlat
    real(kind=real8) :: mjd
    real(kind=real4), dimension(8) :: angles
    integer :: FGF_TYPE = 3 !MTSAT uses JMA GEOS navigation, so set type here

    NAVstr_MTSAT_NAV = NAVstr_MTSAT
    
    imode = -1
    height = 0.0    !Used for parallax correction
    
    Nav%Lat = Missing_Value_Real4
    Nav%Lon = Missing_Value_Real4
       
    if (NAVstr_MTSAT%nav_type == 'GMSX') then
    
        jj = 1 + (ystart - 1)*AREAstr%line_res
    
        do j=1, ysize
            line = real(AREAstr%north_bound) + real(jj - 1) + &
                    real(AREAstr%line_res)/2.0
               
            do i=1, xsize
                ii = ((i+(xstart-1)) - 1)*(AREAstr%elem_res*(xstride)) + 1
        
                elem = real(AREAstr%west_vis_pixel) + real(ii - 1) + &
	                   real(AREAstr%elem_res*(xstride))/2.0
                   
                call MGIVSR(imode,elem,line,dlon,dlat,height,&
                            angles,mjd,ierr)
                
                Space_Mask(i,j) = sym%YES
            
                if (ierr == 0) then
                     Space_Mask(i,j) = sym%NO
                     Nav%Lat_1b(i,j) = dlat
                     Nav%Lon_1b(i,j) = dlon
                endif
            end DO
            
            jj = jj + AREAstr%line_res
        end DO
                        
    endif

    if (NAVstr_MTSAT%nav_type == 'GEOS') then      
        !HRIT requires actual line and element of being processed.
        ! Unlike MSG, MTSAT requires no switching to different corrdinates.
                
          do j=1, ysize
!            jj = ystart + (j-1)
            jj = (AREAstr%north_bound / AREAstr%line_res) + ystart + (j-1)
    
               
            do i=1, xsize
!                ii = (i - 1)*(xstride) + xstart	 ! get element of the image segement
                ii = (AREAstr%west_vis_pixel / AREAstr%elem_res) + (i - 1)*(xstride) + xstart   ! get element of the image segement


                !call pixcoord2geocoord_mtsat(ii, jj, NAVstr_MTSAT%LOFF, NAVstr_MTSAT%COFF, &
                !     latitude,longitude, NAVstr_MTSAT%sub_lon, NAVstr_MTSAT%CFAC, NAVstr_MTSAT%LFAC)

                ! Move to C-routine navigation below.
                ! again, use common algorithm for CGMS navigation
                !call pixcoord2geocoord_cgms(ii,                  &
                !                            jj,                  &
                !                            NAVstr_MTSAT%LOFF,   &
                !                            NAVstr_MTSAT%COFF,   &
                !                            NAVstr_MTSAT%LFAC,   &
                !                            NAVstr_MTSAT%CFAC,   &
                !                            1,             &
                !                            NAVstr_MTSAT%sub_lon, &
                !                            latitude,            &
                !                            longitude)

                 ! Putting in hooks for updated navigation code
                 ! to be edited after 1/12

                 call fgf_to_earth(FGF_TYPE,                  &
                                   DBLE(ii),                  &
                                   DBLE(jj),                  &
                                   DBLE(NAVstr_MTSAT%CFAC),   &
                                   DBLE(NAVstr_MTSAT%COFF),   &
                                   DBLE(NAVstr_MTSAT%LFAC),   &
                                   DBLE(NAVstr_MTSAT%LOFF),   &
                                   DBLE(NAVstr_MTSAT%sub_lon), &
                                   longitude,            &
                                   latitude)

             if (latitude .LE. -999.0) then  ! -999.99 is MSV nav missing value
                    Nav%Lat_1b(i,j) = Missing_Value_Real4
                    Nav%Lon_1b(i,j) = Missing_Value_Real4
                    Space_Mask(i,j) = sym%YES
                else
                    Nav%Lat_1b(i,j) = real(latitude,kind=real4)
                    Nav%Lon_1b(i,j) = real(longitude,kind=real4)
                    
                    ! Because JMA sets their longitudes from 0 to 360, and
                    ! we want 180 to -180, one last check.
                    
                    if (longitude .GT. 180.0 ) then
                        Nav%Lon_1b(i,j) = real(longitude,kind=real4) - 360.0
                    endif
                                        
                    Space_Mask(i,j) = sym%NO

                endif

        
            end DO
                        
        end do     
        
    endif
      
 end subroutine mtsat_navigation
 
 
!------------------------------------------------------------------
! subroutine to convert MTSAT counts to radiance and brightness
! temperature
!------------------------------------------------------------------
  
  subroutine MTSAT_RADIANCE_BT(chan_num,Mtsat_Counts, rad2, temp1)

    integer (kind=INT2), dimension(:,:), intent(in):: Mtsat_Counts
    integer (kind=int1), intent(in) :: chan_num
    real (kind=real4), dimension(:,:), intent(out):: temp1, rad2
    
    integer :: i, j, index
                                   
    do j = 1, Image%Number_Of_Lines_Read_This_Segment
      do i = 1, Image%Number_Of_Elements
       
       index = int(Mtsat_Counts(i,j),kind=int2) + 1
       
       if ((Space_Mask(i,j) == sym%NO) .AND. &
           (index .LE. 1024) .AND. (index .GE. 1)) then 
       !only do valid counts
          rad2(i,j) = real(rad_table(chan_num,1,index),kind=real4)/1000.0
          temp1(i,j) = real(bt_table(chan_num,1,index),kind=real4)/100.0                    
       else
          rad2(i,j) = Missing_Value_Real4
          temp1(i,j) = Missing_Value_Real4
       endif
      end DO
    end do    
  
  end subroutine MTSAT_RADIANCE_BT
  


!---- MTSAT HiRID Navigation


subroutine MGIVSR(IMODE,RPIX,RLIN,RLON,RLAT,RHGT, &
                  RINF,DSCT,IRTN)
 !
 !***********************************************************************
 !
 !  THIS PROGRAM CONVERTS GEOGRAPHICAL CO-ORDINATES (LATITUDE,LONGITUDE,
 !  HEIGHT) TO VISSR IMAGE CO-ORDINATES (LINE,PIXEL) AND VICE VERSA.
 !
 !  THIS PROGRAM IS PROVIDED BY THE METEOROLOGICAL SATELLITE CENTER OF
 !  THE JAPAN METEOROLOGICAL AGENCY TO USERS OF GMS DATA.
 !
 !                                              MSC TECH. NOTE NO.23
 !                                                          JMA/MSC 1991
 !
 !***********************************************************************
 !***********************************************************************
 !***********************************************************************
 !***********************************************************************
 !
 !***********************************************************************
 !            I/O  TYPE
 !   IMODE     I   I*4   CONVERSION MODE & IMAGE kind
 !                         IMAGE kind
 !                               GMS-4 GMS-5  MTSAT
 !                          1,-1  VIS   VIS    VIS
 !                          2,-2  IR    IR1  IR1,IR4
 !                          3,-3  --    IR2    IR2
 !                          4,-4  --    WV     WV
 !                         CONVERSION MODE
 !                           1 TO 4    (LAT,LON,HGT)=>(LINE,PIXEL)
 !                           -1 TO -4  (LAT,LON    )<=(LINE,PIXEL)
 !   RPIX     I/O  R*4   PIXEL OF POINT
 !   RLIN     I/O  R*4   LINE OF POINT
 !   RLON     I/O  R*4   LONGITUDE OF POINT (DEGREES,EAST:+,WEST:-)
 !   RLAT     I/O  R*4   LATITUDE OF POINT (DEGREES,NORTH:+,SOUTH:-)
 !   RHGT      I   R*4   HEIGHT OF POINT (METER)
 !   RINF(8)   O   R*4   (1) SATELLITE ZENITH DISTANCE (DEGREES)
 !                       (2) SATELLITE AZIMUTH ANGLE (DEGREES)
 !                       (3) SUN ZENITH DISTANCE (DEGREES)
 !                       (4) SUN AZIMUTH ANGLE (DEGREES)
 !                       (5) SATELLTE-SUN DIPARTURE ANGLE (DEGREES)
 !                       (6) SATELLTE DISTANCE (METER)
 !                       (7) SUN DISTANCE (KIRO-METER)
 !                       (8) SUN GRINT ANGLE (DEGREES)
 !   DSCT      O   R*8   SCAN TIME (MJD)
 !   IRTN      O   I*4   return CODE (0=O.K.)
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 !   COMMON /MMAP1/MAP(672,4)
 !
 !     1. COORDINATE TRANSFORMATION PARAMETERS SEGMENT
 !                                                  MAP(1,1)-MAP(672,1)
 !     2. ATTITUDE PREDICTION DATA SEGMENT          MAP(1,2)-MAP(672,2)
 !     3. ORBIT PREDICTION DATA 1 SEGMENT           MAP(1,3)-MAP(672,3)
 !     4. ORBIT PREDICTION DATA 2 SEGMENT           MAP(1,4)-MAP(672,4)
 !***********************************************************************
 !
 !!!!!!!!!!!!!!!!!! DEFINITION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !      COMMON /MMAP1/MAP
 !
       integer, intent(in) :: IMODE
       real(kind=real4), intent(inout) :: RPIX, RLIN, RLON, RLAT
       real(kind=real4), intent(in) :: RHGT
       real(kind=real4), dimension(8), intent(out) :: RINF
       real(kind=real8), intent(out) :: DSCT
       integer, intent(out) :: IRTN
       
       integer :: LMODE
       real(kind=real8) :: WKCOS, WKSIN
    
!       real*4     RPIX,RLIN,RLON,RLAT,RHGT,RINF(8)
 !     integer*4  MAP(672,4)
 !
       real*4     EPS,RI0,RI,RJ,RSTEP,RSAMP,RFCL,RFCP,SENS,RFTL,RFTP
       real*4     RESLIN(4),RESELM(4),RLIC(4),RELMFC(4),SENSSU(4), &
                  VMIS(3),ELMIS(3,3),RLINE(4),RELMNT(4)
       real*8     BC,BETA,BS,CDR,CRD,DD,DDA,DDB,DDC,DEF,DK,DK1,DK2, &
                  DLAT,DLON,DPAI,DSPIN,DTIMS,EA,EE,EF,EN,HPAI,PC,PI,PS, &
                  QC,QS,RTIM,TF,TL,TP, &
                  SAT(3),SL(3),SLV(3),SP(3),SS(3),STN1(3),STN2(3), &
                  SX(3),SY(3),SW1(3),SW2(3),SW3(3)
       real*8     DSATZ,DSATA,DSUNZ,DSUNA,DSSDA,DSATD,SUNM,SDIS, &
                  DLATN,DLONN,STN3(3),DSUNG
 !
 !!!!!!!!!!!!!!!!!! EQUIVALENCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !      EQUIVALENCE (MAP( 5,1),DTIMS),    (MAP( 7,1),RESLIN(1))
 !      EQUIVALENCE (MAP(11,1),RESELM(1)),(MAP(15,1),RLIC(1))
 !      EQUIVALENCE (MAP(19,1),RELMFC(1)),(MAP(27,1),SENSSU(1))
 !      EQUIVALENCE (MAP(31,1),RLINE(1)), (MAP(35,1),RELMNT(1))
 !      EQUIVALENCE (MAP(39,1),VMIS(1)),  (MAP(42,1),ELMIS)
 !      EQUIVALENCE (MAP(131,1),DSPIN)
 !
 !***********************************************************************
 !
       DTIMS = NAVstr_MTSAT_NAV%DTIMS
       RESLIN = NAVstr_MTSAT_NAV%RESLIN
       RESELM = NAVstr_MTSAT_NAV%RESELM
       RLIC = NAVstr_MTSAT_NAV%RLIC
       RELMFC = NAVstr_MTSAT_NAV%RELMFC
       SENSSU = NAVstr_MTSAT_NAV%SENSSU
       RLINE = NAVstr_MTSAT_NAV%RLINE
       RELMNT = NAVstr_MTSAT_NAV%RELMNT
       VMIS = NAVstr_MTSAT_NAV%VMIS
       ELMIS = NAVstr_MTSAT_NAV%ELMIS
       DSPIN = NAVstr_MTSAT_NAV%DSPIN     
       
       PI    =  3.141592653D0
       CDR   =  PI/180.D0
       CRD   =  180.D0/PI
       HPAI  =  PI/2.D0
       DPAI  =  PI*2.D0
       EA    =  6378136.D0
       EF    =  1.D0/298.257D0
       EPS   =  1.0
 !!!!!!!!!!!!!!!!!! PARAMETER CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IRTN  =  0
       if(ABS(IMODE).GT.4)  IRTN=1
       if(ABS(RLAT).GT.90. .AND. IMODE.GT.0)  IRTN=2
       if(IRTN.NE.0)  return
 !!!!!!!!!!!!!!!!!! VISSR FRAME INFORMATION SET !!!!!!!!!!!!!!!!!!!!!!!!!
       LMODE   = ABS(IMODE)                                         ![3.1]
       RSTEP   = RESLIN(LMODE)
       RSAMP   = RESELM(LMODE)
       RFCL    = RLIC(LMODE)
       RFCP    = RELMFC(LMODE)
       SENS    = SENSSU(LMODE)
       RFTL    = RLINE(LMODE)+0.5
       RFTP    = RELMNT(LMODE)+0.5
 !!!!!!!!!!!!!!!!!! TRANSFORMATION (GEOGRAPHICAL=>VISSR) !!!!!!!!!!!!!!!!
       if( IMODE.GT.0 .AND. IMODE.LT.5 )  then
         DLAT    = DBLE(RLAT)*CDR                                   ![3.2]
         DLON    = DBLE(RLON)*CDR
         EE      = 2.D0*EF-EF*EF
         EN      = EA/DSQRT(1.D0-EE*DSIN(DLAT)*DSIN(DLAT))
         STN1(1) = (EN+DBLE(RHGT))*DCOS(DLAT)*DCOS(DLON)
         STN1(2) = (EN+DBLE(RHGT))*DCOS(DLAT)*DSIN(DLON)
         STN1(3) = (EN*(1.D0-EE)+DBLE(RHGT))*DSIN(DLAT)
 !
         RI0     = RFCL-ATAN(SIN(SNGL(DLAT))/(6.610689-COS(SNGL(DLAT)))) &
                   /RSTEP
         RTIM    = DTIMS+DBLE(RI0/SENS/1440.)/DSPIN                 ![3.3]
 !
   100   CONTINUE
         call  MGI100(RTIM,CDR,SAT,SP,SS,BETA)                      ![3.4]
 !-----------------------------------------------------------------------
         call  MGI220(SP,SS,SW1)                                    ![3.5]
         call  MGI220(SW1,SP,SW2)
         BC      = DCOS(BETA)
         BS      = DSIN(BETA)
         SW3(1)  = SW1(1)*BS+SW2(1)*BC
         SW3(2)  = SW1(2)*BS+SW2(2)*BC
         SW3(3)  = SW1(3)*BS+SW2(3)*BC
         call  MGI200(SW3,SX)
         call  MGI220(SP,SX,SY)
         SLV(1)  = STN1(1)-SAT(1)                                  ![3.6]
         SLV(2)  = STN1(2)-SAT(2)
         SLV(3)  = STN1(3)-SAT(3)
         call  MGI200(SLV,SL)                                      ![3.7]
         call  MGI210(SP,SL,SW2)
         call  MGI210(SY,SW2,SW3)
         call  MGI230(SY,SW2,TP)
         TF      = SP(1)*SW3(1)+SP(2)*SW3(2)+SP(3)*SW3(3)
         if(TF.LT.0.D0)  TP=-TP
         call  MGI230(SP,SL,TL)
 !
         RI      = SNGL(HPAI-TL)/RSTEP+RFCL-VMIS(2)/RSTEP
         RJ      = SNGL(TP)/RSAMP+RFCP &
                  +VMIS(3)/RSAMP-SNGL(HPAI-TL)*TAN(VMIS(1))/RSAMP
 !
         if(ABS(RI-RI0).GE.EPS)  then                              ![3.8]
           RTIM  = DBLE(AINT((RI-1.)/SENS)+RJ*RSAMP/SNGL(DPAI))/ &
                   (DSPIN*1440.D0)+DTIMS
           RI0   = RI
           GO TO  100
         endif
         RLIN    = RI
         RPIX    = RJ
         DSCT    = RTIM
         if(RLIN.LT.0 .OR. RLIN.GT.RFTL)  IRTN=4
         if(RPIX.LT.0 .OR. RPIX.GT.RFTP)  IRTN=5
 !
 !!!!!!!!!!!!!!!!!! TRANSFORMATION (VISSR=>GEOGRAPHICAL) !!!!!!!!!!!!!!!!
       elseif(IMODE.LT.0 .AND. IMODE.GT.-5)  then
 !
         RTIM    = DBLE(AINT((RLIN-1.)/SENS)+RPIX*RSAMP/SNGL(DPAI))/ &
                   (DSPIN*1440.D0)+DTIMS                           ![3.9]
         call  MGI100(RTIM,CDR,SAT,SP,SS,BETA)                     ![3.10]
         call  MGI220(SP,SS,SW1)                                   ![3.11]
         call  MGI220(SW1,SP,SW2)
         BC      = DCOS(BETA)
         BS      = DSIN(BETA)
         SW3(1)  = SW1(1)*BS+SW2(1)*BC
         SW3(2)  = SW1(2)*BS+SW2(2)*BC
         SW3(3)  = SW1(3)*BS+SW2(3)*BC
         call  MGI200(SW3,SX)
         call  MGI220(SP,SX,SY)
         PC      = DCOS(DBLE(RSTEP*(RLIN-RFCL)))                   ![3.12]
         PS      = DSIN(DBLE(RSTEP*(RLIN-RFCL)))
         QC      = DCOS(DBLE(RSAMP*(RPIX-RFCP)))
         QS      = DSIN(DBLE(RSAMP*(RPIX-RFCP)))
         SW1(1)  = DBLE(ELMIS(1,1))*PC+DBLE(ELMIS(1,3))*PS
         SW1(2)  = DBLE(ELMIS(2,1))*PC+DBLE(ELMIS(2,3))*PS
         SW1(3)  = DBLE(ELMIS(3,1))*PC+DBLE(ELMIS(3,3))*PS
         SW2(1)  = QC*SW1(1)-QS*SW1(2)
         SW2(2)  = QS*SW1(1)+QC*SW1(2)
         SW2(3)  = SW1(3)
         SW3(1)  = SX(1)*SW2(1)+SY(1)*SW2(2)+SP(1)*SW2(3)          ![3.13]
         SW3(2)  = SX(2)*SW2(1)+SY(2)*SW2(2)+SP(2)*SW2(3)
         SW3(3)  = SX(3)*SW2(1)+SY(3)*SW2(2)+SP(3)*SW2(3)
         call  MGI200(SW3,SL)                                      ![3.14]
         DEF     = (1.D0-EF)*(1.D0-EF)
         DDA     = DEF*(SL(1)*SL(1)+SL(2)*SL(2))+SL(3)*SL(3)
         DDB     = DEF*(SAT(1)*SL(1)+SAT(2)*SL(2))+SAT(3)*SL(3)
         DDC     = DEF*(SAT(1)*SAT(1)+SAT(2)*SAT(2)-EA*EA)+SAT(3)*SAT(3)
         DD      = DDB*DDB-DDA*DDC
         if(DD.GE.0.D0 .AND. DDA.NE.0.D0)  then
           DK1     = (-DDB+DSQRT(DD))/DDA
           DK2     = (-DDB-DSQRT(DD))/DDA
         else
           IRTN    = 6
           GO TO  9000
         endif
         if(DABS(DK1).LE.DABS(DK2))  then
           DK    = DK1
         else
           DK    = DK2
         endif
         STN1(1) = SAT(1)+DK*SL(1)
         STN1(2) = SAT(2)+DK*SL(2)
         STN1(3) = SAT(3)+DK*SL(3)
         DLAT    = DATAN(STN1(3)/(DEF*DSQRT(STN1(1)*STN1(1)+ &
                   STN1(2)*STN1(2))))                              ![3.15]
         if(STN1(1).NE.0.D0)  then
           DLON  = DATAN(STN1(2)/STN1(1))
           if(STN1(1).LT.0.D0 .AND. STN1(2).GE.0.D0)  DLON=DLON+PI
           if(STN1(1).LT.0.D0 .AND. STN1(2).LT.0.D0)  DLON=DLON-PI
         else
           if(STN1(2).GT.0.D0)  then
             DLON=HPAI
           else
             DLON=-HPAI
           endif
         endif
         RLAT    = SNGL(DLAT*CRD)
         RLON    = SNGL(DLON*CRD)
         DSCT    = RTIM
       endif
 !
 !!!!!!!!!!!!!!!!!! TRANSFORMATION (ZENITH/AZIMUTH) !!!!!!!!!!!!!!![3.16]
       STN2(1)   = DCOS(DLAT)*DCOS(DLON)                           ![3.17]
       STN2(2)   = DCOS(DLAT)*DSIN(DLON)
       STN2(3)   = DSIN(DLAT)
       SLV(1)    = SAT(1)-STN1(1)                                  ![3.18]
       SLV(2)    = SAT(2)-STN1(2)
       SLV(3)    = SAT(3)-STN1(3)
       call  MGI200(SLV,SL)
 !
       call  MGI230(STN2,SL,DSATZ)                                 ![3.19]
       if(DSATZ.GT.HPAI)  IRTN = 7
 !
       SUNM    = 315.253D0+0.985600D0*RTIM                         ![3.20]
       SUNM    = DMOD(SUNM,360.D0)*CDR
       SDIS    = (1.00014D0-0.01672D0*DCOS(SUNM)-0.00014*DCOS(2.D0* &
                 SUNM))*1.49597870D8
 !
       if(DLAT.GE.0.D0) then                                       ![3.21]
         DLATN   = HPAI-DLAT
         DLONN = DLON-PI
         if(DLONN.LE.-PI)  DLONN=DLONN+DPAI
       else
         DLATN   = HPAI+DLAT
         DLONN   = DLON
       endif
       STN3(1) = DCOS(DLATN)*DCOS(DLONN)
       STN3(2) = DCOS(DLATN)*DSIN(DLONN)
       STN3(3) = DSIN(DLATN)
       SW1(1)  = SLV(1)+SS(1)*SDIS*1.D3                            ![3.22]
       SW1(2)  = SLV(2)+SS(2)*SDIS*1.D3
       SW1(3)  = SLV(3)+SS(3)*SDIS*1.D3
       call  MGI200(SW1,SW2)                                       ![3.23]
       call  MGI230(STN2,SW2,DSUNZ)
       call  MGI230(SL,SW2,DSSDA)                                  ![3.24]
       call  MGI240(SL,STN2,STN3,DPAI,DSATA)                       ![3.25]
       call  MGI240(SW2,STN2,STN3,DPAI,DSUNA)                      ![3.26]
       DSATD   = DSQRT(SLV(1)*SLV(1)+SLV(2)*SLV(2)+SLV(3)*SLV(3))  ![3.27]
 !
 !
       call  MGI200(STN1,SL)                                       ![3.28]
       call  MGI230(SW2,SL,DSUNG)
       call  MGI220(SL, SW2,SW3)
       call  MGI220(SW3,SL, SW1)
       WKCOS=DCOS(DSUNG)
       WKSIN=DSIN(DSUNG)
       SW2(1)=WKCOS*SL(1)-WKSIN*SW1(1)
       SW2(2)=WKCOS*SL(2)-WKSIN*SW1(2)
       SW2(3)=WKCOS*SL(3)-WKSIN*SW1(3)
       call  MGI230(SW2,SLV,DSUNG)
 !
       RINF(6) = SNGL(DSATD)
       RINF(7) = SNGL(SDIS)
       RINF(1) = SNGL(DSATZ*CRD)
       RINF(2) = SNGL(DSATA*CRD)
       RINF(3) = SNGL(DSUNZ*CRD)
       RINF(4) = SNGL(DSUNA*CRD)
       RINF(5) = SNGL(DSSDA*CRD)
       RINF(8) = SNGL(DSUNG*CRD)
 !!!!!!!!!!!!!!!!!!! STOP/end !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  9000 CONTINUE
       return
end subroutine MGIVSR

subroutine  MGI100(RTIM,CDR,SAT,SP,SS,BETA)
!       COMMON /MMAP1/MAP
       real*8    ATTALP,ATTDEL,BETA,CDR,DELT,RTIM,SITAGT,SUNALP,SUNDEL, &
                 WKCOS,WKSIN
       real*8    ATT1(3),ATT2(3),ATT3(3),NPA(3,3), &
                 SAT(3),SP(3),SS(3),ORBT2(35,8)
 !     integer*4 MAP(672,4)
       integer :: I
 !
 !      EQUIVALENCE (MAP(13,3),ORBT1(1,1))
 !      EQUIVALENCE (MAP(13,2),ATIT(1,1))
 !
       do 1000 I=1,7
         if(RTIM.GE.NAVstr_MTSAT_NAV%ORBT1(1,I).AND.RTIM.LT.NAVstr_MTSAT_NAV%ORBT1(1,I+1))  then
           call  MGI110 &
                (I,RTIM,CDR,NAVstr_MTSAT_NAV%ORBT1,ORBT2,SAT,SITAGT,SUNALP,SUNDEL,NPA)
           GO TO  1200
         endif
  1000 CONTINUE
  1200 CONTINUE
 !
       do 3000 I=1,9
         if(RTIM.GE.NAVstr_MTSAT_NAV%ATIT(1,I) .AND. RTIM.LT.NAVstr_MTSAT_NAV%ATIT(1,I+1))  then
           DELT = (RTIM-NAVstr_MTSAT_NAV%ATIT(1,I))/(NAVstr_MTSAT_NAV%ATIT(1,I+1)-NAVstr_MTSAT_NAV%ATIT(1,I))
           ATTALP = NAVstr_MTSAT_NAV%ATIT(3,I)+(NAVstr_MTSAT_NAV%ATIT(3,I+1)-NAVstr_MTSAT_NAV%ATIT(3,I))*DELT
           ATTDEL = NAVstr_MTSAT_NAV%ATIT(4,I)+(NAVstr_MTSAT_NAV%ATIT(4,I+1)-NAVstr_MTSAT_NAV%ATIT(4,I))*DELT
           BETA   = NAVstr_MTSAT_NAV%ATIT(5,I)+(NAVstr_MTSAT_NAV%ATIT(5,I+1)-NAVstr_MTSAT_NAV%ATIT(5,I))*DELT
           if( (NAVstr_MTSAT_NAV%ATIT(5,I+1)-NAVstr_MTSAT_NAV%ATIT(5,I)).GT.0.D0 ) &
             BETA   = NAVstr_MTSAT_NAV%ATIT(5,I)+(NAVstr_MTSAT_NAV%ATIT(5,I+1)-NAVstr_MTSAT_NAV%ATIT(5,I)-360.D0*CDR)*DELT
           GO TO  3001
         endif
  3000 CONTINUE
  3001 CONTINUE
 !
       WKCOS    = DCOS(ATTDEL)
       ATT1(1)  = DSIN(ATTDEL)
       ATT1(2)  = WKCOS       *(-DSIN(ATTALP))
       ATT1(3)  = WKCOS       *DCOS(ATTALP)
       ATT2(1)  = NPA(1,1)*ATT1(1)+NPA(1,2)*ATT1(2)+NPA(1,3)*ATT1(3)
       ATT2(2)  = NPA(2,1)*ATT1(1)+NPA(2,2)*ATT1(2)+NPA(2,3)*ATT1(3)
       ATT2(3)  = NPA(3,1)*ATT1(1)+NPA(3,2)*ATT1(2)+NPA(3,3)*ATT1(3)
       WKSIN    = DSIN(SITAGT)
       WKCOS    = DCOS(SITAGT)
       ATT3(1)  = WKCOS*ATT2(1)+WKSIN*ATT2(2)
       ATT3(2)  =-WKSIN*ATT2(1)+WKCOS*ATT2(2)
       ATT3(3)  = ATT2(3)
       call  MGI200(ATT3,SP)
 !
       WKCOS    = DCOS(SUNDEL)
       SS(1)    = WKCOS       *DCOS(SUNALP)
       SS(2)    = WKCOS       *DSIN(SUNALP)
       SS(3)    = DSIN(SUNDEL)
 !
       return
end subroutine MGI100

subroutine MGI110(I,RTIM,CDR,ORBTA,ORBTB,SAT,SITAGT,SUNALP,SUNDEL,NPA)
       real*8    CDR,SAT(3),RTIM,ORBTA(35,8),ORBTB(35,8)
       real*8    SITAGT,SUNDEL,SUNALP,NPA(3,3),DELT
       integer*4 I
       if(I.NE.8)  then
         DELT=(RTIM-ORBTA(1,I))/(ORBTA(1,I+1)-ORBTA(1,I))
         SAT(1)   = ORBTA( 9,I)+(ORBTA( 9,I+1)-ORBTA( 9,I))*DELT
         SAT(2)   = ORBTA(10,I)+(ORBTA(10,I+1)-ORBTA(10,I))*DELT
         SAT(3)   = ORBTA(11,I)+(ORBTA(11,I+1)-ORBTA(11,I))*DELT
         SITAGT   =(ORBTA(15,I)+(ORBTA(15,I+1)-ORBTA(15,I))*DELT)*CDR
         if( (ORBTA(15,I+1)-ORBTA(15,I)).LT.0.D0 ) &
           SITAGT   =(ORBTA(15,I)+(ORBTA(15,I+1)-ORBTA(15,I)+360.D0) &
                     *DELT)*CDR
         SUNALP   =(ORBTA(18,I)+(ORBTA(18,I+1)-ORBTA(18,I))*DELT)*CDR
         if( (ORBTA(18,I+1)-ORBTA(18,I)).GT.0.D0 ) &
           SUNALP   =(ORBTA(18,I)+(ORBTA(18,I+1)-ORBTA(18,I)-360.D0) &
                     *DELT)*CDR
         SUNDEL   =(ORBTA(19,I)+(ORBTA(19,I+1)-ORBTA(19,I))*DELT)*CDR
         NPA(1,1) = ORBTA(20,I)
         NPA(2,1) = ORBTA(21,I)
         NPA(3,1) = ORBTA(22,I)
         NPA(1,2) = ORBTA(23,I)
         NPA(2,2) = ORBTA(24,I)
         NPA(3,2) = ORBTA(25,I)
         NPA(1,3) = ORBTA(26,I)
         NPA(2,3) = ORBTA(27,I)
         NPA(3,3) = ORBTA(28,I)
       endif
       return
end subroutine MGI110

subroutine MGI200(VECT,VECTU)
       real*8  VECT(3),VECTU(3),RV1,RV2
       RV1=VECT(1)*VECT(1)+VECT(2)*VECT(2)+VECT(3)*VECT(3)
       if(RV1 == 0.D0)  return
       RV2=DSQRT(RV1)
       VECTU(1)=VECT(1)/RV2
       VECTU(2)=VECT(2)/RV2
       VECTU(3)=VECT(3)/RV2
       return
end subroutine MGI200

subroutine MGI210(VA,VB,VC)
       real*8  VA(3),VB(3),VC(3)
       VC(1)= VA(2)*VB(3)-VA(3)*VB(2)
       VC(2)= VA(3)*VB(1)-VA(1)*VB(3)
       VC(3)= VA(1)*VB(2)-VA(2)*VB(1)
       return
end subroutine MGI210

subroutine MGI220(VA,VB,VD)
       real*8  VA(3),VB(3),VC(3),VD(3)
       VC(1)= VA(2)*VB(3)-VA(3)*VB(2)
       VC(2)= VA(3)*VB(1)-VA(1)*VB(3)
       VC(3)= VA(1)*VB(2)-VA(2)*VB(1)
       call  MGI200(VC,VD)
       return
end subroutine MGI220

subroutine MGI230(VA,VB,ASITA)
       real*8  VA(3),VB(3),ASITA,AS1,AS2
       AS1= VA(1)*VB(1)+VA(2)*VB(2)+VA(3)*VB(3)
       AS2=(VA(1)*VA(1)+VA(2)*VA(2)+VA(3)*VA(3))* &
           (VB(1)*VB(1)+VB(2)*VB(2)+VB(3)*VB(3))
       if(AS2 == 0.D0)  return
       ASITA=DACOS(AS1/DSQRT(AS2))
       return
end subroutine MGI230

subroutine MGI240(VA,VH,VN,DPAI,AZI)
       real*8  VA(3),VH(3),VN(3),VB(3),VC(3),VD(3),DPAI,AZI,DNAI
       call  MGI220(VN,VH,VB)
       call  MGI220(VA,VH,VC)
       call  MGI230(VB,VC,AZI)
       call  MGI220(VB,VC,VD)
       DNAI = VD(1)*VH(1)+VD(2)*VH(2)+VD(3)*VH(3)
       if(DNAI.GT.0.D0)  AZI=DPAI-AZI
       return
end subroutine MGI240


!------------------- MTSAT NAV BLOC

 subroutine READ_NAVIGATION_BLOCK_MTSAT_FY(filename, AREAstr, NAVstr)
  character(len=*), intent(in):: filename
  TYPE(AREA_STRUCT), intent(in):: AREAstr
  TYPE(GVAR_NAV), intent(inout):: NAVstr
 
  character(len=1), dimension(3200) :: CBUF
  integer :: i, j, geos_nav
  integer(kind=int4)nav_offset
  real(kind=real4) :: R4DMY
  real(kind=real8) :: R8DMY  
! real(kind=real8) :: LOFF,COFF, LFAC, CFAC  
  integer:: number_of_words_read
  integer(kind=int4), dimension(640) :: i4buf
  
  nav_offset = AREAstr%sec_key_nav
    
  !determine GEOS or other navigation
  geos_nav = sym%NO
  call mreadf_int(trim(filename)//CHAR(0),nav_offset,4,640,&
                    number_of_words_read, i4buf)
!  if (AREAstr%swap_bytes > 0) call swap_bytes4(i4buf,640)
  call move_bytes(4,i4buf(1),NAVstr%nav_type,0)
  
  if (NAVstr%nav_type == 'GEOS') then
        !SUBLON stored as SUBLON *10 in McIDAS NAV block
        NAVstr%sub_lon = real(i4buf(6),kind=real4) / 10
        NAVstr%sublon = NAVstr%sub_lon

        ! LOFF, COFF, CFAC, LFAC stored in McIDAS header for 1km data. All
        ! Multipied by 10. Order from nvxmtst.dlm in McIDAS
        NAVstr%LOFF=(i4buf(2) / 10 ) / real(AREAstr%line_res)
        NAVstr%COFF=(i4buf(3) / 10) / real(AREAstr%elem_res)
        NAVstr%LFAC=(i4buf(4) / 10 ) / real(AREAstr%line_res)
        NAVstr%CFAC=(i4buf(5) / 10 ) / real(AREAstr%elem_res)
     
  endif
    
  if (NAVstr%nav_type == 'GMSX') then
    call mreadf_int(trim(filename)//CHAR(0),nav_offset,1,3200,&
                    number_of_words_read,NAVstr%COBAT)
        
    CBUF(1:504) = NAVstr%COBAT(5:508)
    CBUF(505:1012) = NAVstr%COBAT(513:1020)
    CBUF(1013:1520) = NAVstr%COBAT(1025:1532)
    CBUF(1521:2028) = NAVstr%COBAT(1537:2044)
    CBUF(2029:2536) = NAVstr%COBAT(2049:2556)
    
    call SV0100( 6, 8, CBUF( 1:   6), R4DMY, NAVstr%DTIMS )
    call SV0100( 4, 8, CBUF( 7:  10), NAVstr%RESLIN(1), R8DMY )
    call SV0100( 4, 8, CBUF( 11: 14), NAVstr%RESLIN(2), R8DMY )
    call SV0100( 4, 8, CBUF( 11: 14), NAVstr%RESLIN(3), R8DMY )
    call SV0100( 4, 8, CBUF( 11: 14), NAVstr%RESLIN(4), R8DMY )
    call SV0100( 4,10, CBUF( 15: 18), NAVstr%RESELM(1), R8DMY )
    call SV0100( 4,10, CBUF( 19: 22), NAVstr%RESELM(2), R8DMY )
    call SV0100( 4,10, CBUF( 19: 22), NAVstr%RESELM(3), R8DMY )
    call SV0100( 4,10, CBUF( 19: 22), NAVstr%RESELM(4), R8DMY )
    call SV0100( 4, 4, CBUF( 23: 26), NAVstr%RLIC(1), R8DMY )
    call SV0100( 4, 4, CBUF( 27: 30), NAVstr%RLIC(2), R8DMY )
    call SV0100( 4, 4, CBUF(111:114), NAVstr%RLIC(3), R8DMY )
    call SV0100( 4, 4, CBUF(115:118), NAVstr%RLIC(4), R8DMY )
    call SV0100( 4, 4, CBUF( 31: 34), NAVstr%RELMFC(1), R8DMY )
    call SV0100( 4, 4, CBUF( 35: 38), NAVstr%RELMFC(2), R8DMY )
    call SV0100( 4, 4, CBUF(119:122), NAVstr%RELMFC(3), R8DMY )
    call SV0100( 4, 4, CBUF(123:126), NAVstr%RELMFC(4), R8DMY )
    call SV0100( 4, 0, CBUF( 39: 42), NAVstr%SENSSU(1), R8DMY )
    call SV0100( 4, 0, CBUF( 43: 46), NAVstr%SENSSU(2), R8DMY )
    call SV0100( 4, 0, CBUF( 43: 46), NAVstr%SENSSU(3), R8DMY )
    call SV0100( 4, 0, CBUF( 43: 46), NAVstr%SENSSU(4), R8DMY )
    call SV0100( 4, 0, CBUF( 47: 50), NAVstr%RLINE(1) , R8DMY )
    call SV0100( 4, 0, CBUF( 51: 54), NAVstr%RLINE(2) , R8DMY )
    call SV0100( 4, 0, CBUF( 51: 54), NAVstr%RLINE(3) , R8DMY )
    call SV0100( 4, 0, CBUF( 51: 54), NAVstr%RLINE(4) , R8DMY )
    call SV0100( 4, 0, CBUF( 55: 58), NAVstr%RELMNT(1), R8DMY )
    call SV0100( 4, 0, CBUF( 59: 62), NAVstr%RELMNT(2), R8DMY )
    call SV0100( 4, 0, CBUF( 59: 62), NAVstr%RELMNT(3), R8DMY )
    call SV0100( 4, 0, CBUF( 59: 62), NAVstr%RELMNT(4), R8DMY )
    call SV0100( 4,10, CBUF( 63: 66), NAVstr%VMIS(1), R8DMY )
    call SV0100( 4,10, CBUF( 67: 70), NAVstr%VMIS(2), R8DMY )
    call SV0100( 4,10, CBUF( 71: 74), NAVstr%VMIS(3), R8DMY )
    call SV0100( 4, 7, CBUF( 75: 78), NAVstr%ELMIS(1,1), R8DMY )
    call SV0100( 4,10, CBUF( 79: 82), NAVstr%ELMIS(2,1), R8DMY )
    call SV0100( 4,10, CBUF( 83: 86), NAVstr%ELMIS(3,1), R8DMY )
    call SV0100( 4,10, CBUF( 87: 90), NAVstr%ELMIS(1,2), R8DMY )
    call SV0100( 4, 7, CBUF( 91: 94), NAVstr%ELMIS(2,2), R8DMY )
    call SV0100( 4,10, CBUF( 95: 98), NAVstr%ELMIS(3,2), R8DMY )
    call SV0100( 4,10, CBUF( 99:102), NAVstr%ELMIS(1,3), R8DMY )
    call SV0100( 4,10, CBUF(103:106), NAVstr%ELMIS(2,3), R8DMY )
    call SV0100( 4, 7, CBUF(107:110), NAVstr%ELMIS(3,3), R8DMY )
    call SV0100( 6, 8, CBUF(241:246), R4DMY, NAVstr%DSPIN )
    call SV0100( 6, 6, CBUF(199:204), NAVstr%sublon, R8DMY )
    call SV0100( 6, 6, CBUF(205:210), NAVstr%sublat, R8DMY )
    
    do i=1, 10
        j = (i-1)*48 + 256
        call SV0100(6, 8,CBUF( 1+j: 6+j),R4DMY,NAVstr%ATIT(1,i))
        call SV0100(6, 8,CBUF(13+j:18+j),R4DMY,NAVstr%ATIT(3,i))
        call SV0100(6,11,CBUF(19+j:24+j),R4DMY,NAVstr%ATIT(4,i))
        call SV0100(6, 8,CBUF(25+j:30+j),R4DMY,NAVstr%ATIT(5,i))
        call SV0100(6, 8,CBUF(31+j:36+j),R4DMY,NAVstr%ATIT(6,i))
    end DO
    
    do i=1, 8
        j = (i-1)*200 + (256 + 10*48)
        call SV0100(6, 8,CBUF(1+j:6+j),R4DMY,NAVstr%ORBT1( 1,i))
        call SV0100(6, 6,CBUF(49+j:54+j),R4DMY,NAVstr%ORBT1( 9,i))
        call SV0100(6, 6,CBUF(55+j:60+j),R4DMY,NAVstr%ORBT1(10,i))
        call SV0100(6, 6,CBUF(61+j:66+j),R4DMY,NAVstr%ORBT1(11,i))
        call SV0100(6, 8,CBUF(85+j:90+j),R4DMY,NAVstr%ORBT1(15,i))
        call SV0100(6, 8,CBUF(103+j:108+j),R4DMY,NAVstr%ORBT1(18,i))
        call SV0100(6, 8,CBUF(109+j:114+j),R4DMY,NAVstr%ORBT1(19,i))
        call SV0100(6,12,CBUF(129+j:134+j),R4DMY,NAVstr%ORBT1(20,i))
        call SV0100(6,14,CBUF(135+j:140+j),R4DMY,NAVstr%ORBT1(21,i))
        call SV0100(6,14,CBUF(141+j:146+j),R4DMY,NAVstr%ORBT1(22,i))
        call SV0100(6,14,CBUF(147+j:152+j),R4DMY,NAVstr%ORBT1(23,i))
        call SV0100(6,12,CBUF(153+j:158+j),R4DMY,NAVstr%ORBT1(24,i))
        call SV0100(6,16,CBUF(159+j:164+j),R4DMY,NAVstr%ORBT1(25,i))
        call SV0100(6,12,CBUF(165+j:170+j),R4DMY,NAVstr%ORBT1(26,i))
        call SV0100(6,16,CBUF(171+j:176+j),R4DMY,NAVstr%ORBT1(27,i))
        call SV0100(6,12,CBUF(177+j:182+j),R4DMY,NAVstr%ORBT1(28,i))
    end DO
    NAVstr%sub_lon = NAVstr%sublon    
  endif

 end subroutine READ_NAVIGATION_BLOCK_MTSAT_FY


 subroutine  SV0100( IWORD, IPOS, C, R4DAT, R8DAT )
 !---- ------------------------------------------------------------------
 !     TYPE CONVERT ROUTINE ( R-TYPE )
 !---- ------------------------------------------------------------------
       integer*4  IWORD,IPOS,IDATA1
       character  C(*)*1
       real*4     R4DAT
       real*8     R8DAT
       R4DAT = 0.0
       R8DAT = 0.D0
       if( IWORD == 4 )  then
         IDATA1 = ICHAR( C(1)(1:1) )/128
         R8DAT  = DFLOAT( MOD(ICHAR(C(1)(1:1)),128) )*2.D0**(8*3)+ &
                  DFLOAT( ICHAR(C(2)(1:1)) )*2.D0**(8*2)+ &
                  DFLOAT( ICHAR(C(3)(1:1)) )*2.D0**(8*1)+ &
                  DFLOAT( ICHAR(C(4)(1:1)) )
         R8DAT  = R8DAT/10.D0**IPOS
         if( IDATA1 == 1 )  R8DAT = -R8DAT
         R4DAT  = SNGL( R8DAT )
       elseif( IWORD == 6 )  then
         IDATA1 = ICHAR( C(1)(1:1) )/128
         R8DAT  = DFLOAT( MOD(ICHAR(C(1)(1:1)),128) )*2.D0**(8*5)+ &
                  DFLOAT( ICHAR(C(2)(1:1)) )*2.D0**(8*4)+ &
                  DFLOAT( ICHAR(C(3)(1:1)) )*2.D0**(8*3)+ &
                  DFLOAT( ICHAR(C(4)(1:1)) )*2.D0**(8*2)+ &
                  DFLOAT( ICHAR(C(5)(1:1)) )*2.D0**(8*1)+ &
                  DFLOAT( ICHAR(C(6)(1:1)) )
         R8DAT  = R8DAT/10.D0**IPOS
         if( IDATA1 == 1 )  R8DAT = -R8DAT
         R4DAT  = SNGL( R8DAT )
       endif
       return
end subroutine SV0100

!---------------------------------------------------------------------------------------------
! Calibrate the Ch1 Dark Counts using the Mtsat Reflectance Calibration Function
!---------------------------------------------------------------------------------------------
subroutine CALIBRATE_MTSAT_DARK_COMPOSITE(Dark_Comp_Counts,Ref_Ch1_Dark)

integer(kind=int2), dimension(:,:), intent(in):: Dark_Comp_Counts
real(kind=real4), dimension(:,:),  intent(out):: Ref_Ch1_Dark

call MTSAT_REFLECTANCE_PRELAUNCH(Dark_Comp_Counts,Ref_Ch1_Dark)

end subroutine CALIBRATE_MTSAT_DARK_COMPOSITE

!
!--- end of module
!
end module MTSAT_MODULE
