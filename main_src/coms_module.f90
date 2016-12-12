!// $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: coms_module.f90 (src)
!       COMS_module (program)
!
! PURPOSE: This module contains all the subroutines needed to perform navigation and
!          calibration for COMS, both HiRID and HRIT
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

module COMS_MODULE

use CX_CONSTANTS_MOD
use PIXEL_COMMON
use CALIBRATION_CONSTANTS
use PLANCK
use NUMERICAL_TOOLS_MOD
use CGMS_NAV
use GOES_MODULE
use FILE_TOOLS, only: &
   get_lun
use VIEWING_GEOMETRY_MODULE
   
   use CX_SSEC_AREAFILE_MOD
 implicit none
   private
 public :: READ_COMS
 public:: READ_NAVIGATION_BLOCK_COMS
 public:: GET_COMS_IMAGE
 public:: READ_COMS_INSTR_CONSTANTS
         
 private :: COMS_RADIANCE_BT, &
            COMS_REFLECTANCE_PRELAUNCH, &
            COMS_NAVIGATION
 
 type (GVAR_NAV), PRIVATE    :: NAVstr_COMS_NAV
 integer, PARAMETER, PRIVATE :: Nchan_COMS= 5
 integer, PARAMETER, PRIVATE :: Ndet_COMS = 4
 integer, PARAMETER, PRIVATE :: Ntable_COMS = 1024

 integer, PRIVATE :: Nref_Table_COMS
 integer, PRIVATE :: Nbt_Table_COMS
 CHARACTER(len=4), PRIVATE:: Calib_Type

 !---stw integer (kind=int4), dimension(nchan_COMS,ndet_COMS,Ntable_COMS), PRIVATE  :: Ref_Table
 real (kind=real4), dimension(Ntable_COMS), PRIVATE  :: Ref_Table
 integer (kind=int4), dimension(nchan_COMS,ndet_COMS,Ntable_COMS), PRIVATE  :: bt_table
 integer (kind=int4), dimension(nchan_COMS,ndet_COMS,Ntable_COMS), PRIVATE  :: rad_table

 integer(kind=int4), private, parameter:: COMS_Xstride = 1
 !---stw Full disk COMS scans are 2750x2750
 !---stw integer(kind=int4), private, parameter:: num_4km_scans_fd = 3712
 !---stw integer(kind=int4), private, parameter:: num_4km_elem_fd = 3712
 integer(kind=int4), private, parameter:: num_4km_scans_fd = 2750
 integer(kind=int4), private, parameter:: num_4km_elem_fd = 2750
 integer(kind=int4), private, parameter:: time_for_fd_scan =  1560000 !milliseconds (26min)
 real, private, save:: Scan_rate    !scan rate in millsec / line

 integer(kind=int4), private :: ptr

 CONTAINS
!----------------------------------------------------------------
! read the COMS constants into memory
!-----------------------------------------------------------------
subroutine READ_COMS_INSTR_CONSTANTS(Instr_Const_file)
 character(len=*), intent(in):: Instr_Const_file
 integer:: ios0, erstat
 integer:: Instr_Const_lun

 Instr_Const_lun = GET_LUN()
 open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

 print *, "opening ", trim(Instr_Const_file)
 erstat = 0
 if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, "Error opening COMS constants file, ios0 = ", ios0
    stop 19
 endif
  read(unit=Instr_Const_lun,fmt="(a3)") sat_name
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


end subroutine READ_COMS_INSTR_CONSTANTS

 ! Perform COMS Reflectance and BT calibration
 subroutine READ_COMS(segment_number,channel_1_filename, &
                     jday, image_time_ms, &
                     AREAstr,NAVstr_COMS)

   integer(kind=int4), intent(in):: segment_number
   character(len=*), intent(in):: channel_1_filename
   type (area_header_type), intent(in) :: AREAstr
   type (GVAR_NAV), intent(in)    :: NAVstr_COMS
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
   integer:: COMS_file_id
   real(kind=real4):: image_time_hours
   integer(kind=int4):: image_jday
   integer(kind=int4):: first_line_in_segment
   character(len=2):: ichan_goes_string
   integer :: Line_Idx
   integer :: Elem_Idx
   integer:: num_elements_this_image
   integer:: num_scans_this_image

   

   !--- assume channel_1_file name has a unique "_1_" in the name. 
   !--- determine indices needed to replace that string
   ipos = index(channel_1_filename, "_1_")
   ilen = len(channel_1_filename)
    
   first_line_in_segment = (segment_number-1)*Image%Number_Of_Lines_Per_Segment

   !---------------------------------------------------------------------------
   ! COMS Navigation (Do Navigation and Solar angles first)
   !---------------------------------------------------------------------------
   
   call COMS_NAVIGATION(1,first_line_in_segment,&
                              Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment,1,&
                              AREAstr,NAVstr_COMS)
   
   if (segment_number == 1) then

      image_jday = jday
      image_time_hours = image_time_ms / 60.0 / 60.0 / 1000.0

      !--- compute scan rate for future use
      num_elements_this_image =  int(AREAstr%num_elem / COMS_Xstride) + 1
      num_scans_this_image = AREAstr%num_line
      Scan_Rate = real((num_elements_this_image)/               &
               real(num_4km_elem_fd/COMS_Xstride)) * &
               real((num_scans_this_image) / real(num_4km_scans_fd)) * &
               real(time_for_fd_scan) / real(num_scans_this_image)

       do ichan_goes = 2,5

       if (ichan_goes == 2) ichan_modis = 20
       if (ichan_goes == 3) ichan_modis = 27
       if (ichan_goes == 4) ichan_modis = 31
       if (ichan_goes == 5) ichan_modis = 32
       
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

    ! On first segment, reflectance, BT and rad tables
    ! On first segment, get slope/offset information from McIDAS Header
    COMS_file_id = get_lun()   
    if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) then
        call mread_open(trim(Temporary_Data_Dir)//trim(channel_1_filename)//CHAR(0), COMS_file_id)
    else
        call mread_open(trim(Image%Level1b_Path)//trim(channel_1_filename)//CHAR(0), COMS_file_id)
    endif  

    call load_COMS_calibration(COMS_file_id, AREAstr)
    call mread_close(COMS_file_id)


   endif
   
   
   !---   read channel 1 (COMS channel 1)
   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_1_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_1_filename)
       endif

       call GET_COMS_IMAGE(trim(channel_x_filename_full), &
                                    AREAstr, &
                                    segment_number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp)
                                    
        call COMS_REFLECTANCE_PRELAUNCH(Two_Byte_Temp,ch(1)%Ref_Toa(:,:))
!cspp   call COMS_REFLECTANCE_PRELAUNCH(Two_Byte_Temp,Ref_Ch1)

        Ch1_Counts = Two_Byte_Temp

   endif
   
   !---   read channel 20 (COMS channel 2)
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_2_" // &
                            channel_1_filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif

       call GET_COMS_IMAGE(trim(channel_x_filename_full), &
                                    AREAstr, &
                                    segment_number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp)
      call COMS_RADIANCE_BT(2_int1, Two_Byte_Temp, ch(20)%Rad_Toa, ch(20)%Bt_Toa)
!cspp call COMS_RADIANCE_BT(2_int1, Two_Byte_Temp, Rad_Ch20, Bt_Ch20)

   endif
                    
   
   !---   read channel 27 (COMS channel 3)
   if (Sensor%Chan_On_Flag_Default(27) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_3_" // &
                            channel_1_filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif

       call GET_COMS_IMAGE(trim(channel_x_filename_full), &
                                    AREAstr, &
                                    segment_number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp)
      call COMS_RADIANCE_BT(3_int1, Two_Byte_Temp, ch(27)%Rad_Toa, ch(27)%Bt_Toa)
!cspp call COMS_RADIANCE_BT(3_int1, Two_Byte_Temp, Rad_Ch27, Bt_Ch27)

   endif


   
   !---   read channel 31 (COMS channel 4)
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_4_" // &
                            channel_1_filename(ipos+3:ilen)
       
       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif
       
       call GET_COMS_IMAGE(trim(channel_x_filename_full), &
                                    AREAstr, &
                                    segment_number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp)
      call COMS_RADIANCE_BT(4_int1, Two_Byte_Temp, ch(31)%Rad_Toa, ch(31)%Bt_Toa)
!cspp call COMS_RADIANCE_BT(4_int1, Two_Byte_Temp, Rad_Ch31, Bt_Ch31)

   endif
   
   
   !---   read channel 32 (COMS channel 5)
   if (Sensor%Chan_On_Flag_Default(32) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_5_" // &
                            channel_1_filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif

       call GET_COMS_IMAGE(trim(channel_x_filename_full), &
                                    AREAstr, &
                                    segment_number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp)
      call COMS_RADIANCE_BT(5_int1, Two_Byte_Temp, ch(32)%Rad_Toa, ch(32)%Bt_Toa)
!cspp call COMS_RADIANCE_BT(5_int1, Two_Byte_Temp, Rad_Ch32, Bt_Ch32)

   endif
    
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     Scan_Number(Line_Idx) = first_line_in_segment + Line_Idx
     Scan_Time(Line_Idx) = image_time_ms + (Scan_Number(Line_Idx)-1) * Scan_rate
   enddo

!------------------------------------------------------------------------------
! COMS Angles
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
    
 end subroutine READ_COMS
 
 subroutine load_COMS_calibration(lun, AREAstr)
  integer(kind=int4), intent(in) :: lun
  type(area_header_type), intent(in):: AREAstr
  integer(kind=int4), dimension(6528) :: ibuf
  CHARACTER(len=25) :: cbuf
  integer :: nref, nbt, i, offset
  integer(kind=int4) :: band_offset_2, band_offset_14, band_offset_15, &
                        band_offset_9, band_offset_7, dir_offset
  real(kind=real4) :: albedo, temperature, radiance

  call mreadf_int_o(lun,AREAstr%cal_offset,4,6528,ibuf)
  !if (AREAstr%swap_bytes > 0) call swap_bytes4(ibuf,6528)
  
  dir_offset = ibuf(4)
  band_offset_2 = ibuf(6)
  band_offset_7 = ibuf(8)
  band_offset_9 = ibuf(10)
  band_offset_14 = ibuf(12)
  band_offset_15 = ibuf(14)  
  
  Nref = 256
  Nbt = 1024
  Nref_Table_COMS = nref
  Nbt_Table_COMS = nbt
  
  ! We need to get the vis calibration type. This is in 5th
  ! 4 character long block of the cal block (as per line 284, kbxmtst.dlm, v1.5)
  ! WCS3 (3/29/2010)
  
  call mreadf_int_o(lun,AREAstr%cal_offset,1,25,cbuf)
  calib_type = cbuf(17:20)

  !---------------------------------------------------------------------
  ! Load the visible channel calibration from the McIDAS AREA file
  !---------------------------------------------------------------------
  
  ptr = Band_Offset_2/4

  do i = 1, 1024

    !---stw At this point, it still needs a divide by 10 to
    !---stw to match McIDAS-X.

    Albedo = real(ibuf(ptr),kind=real4) / 100.0

    Ref_Table(i) = Albedo / 10.0

    ptr = ptr+1

  end do

  !---------------------------------------------------------------------
  ! Load the IR channel calibration table.
  !---------------------------------------------------------------------
  
  do i=1, 1024
    !  3.8 micron
    offset = band_offset_7/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(2,1,i) = nint(temperature * 100.)
    radiance = PLANCK_RAD(20, temperature)
    rad_table(2,1,i) = nint(radiance * 1000.)
    
    !  6.8 micron
    offset = band_offset_9/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(3,1,i) = nint(temperature * 100.)
    radiance = PLANCK_RAD(27, temperature)
    rad_table(3,1,i) = nint(radiance * 1000.)
    
    ! 10.8 micron
    offset = band_offset_14/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(4,1,i) = nint(temperature * 100.)
    radiance = PLANCK_RAD(31, temperature)
    rad_table(4,1,i) = nint(radiance * 1000.)
        
    ! 12.0 micron
    offset = band_offset_15/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(5,1,i) = nint(temperature * 100.)
    radiance = PLANCK_RAD(32, temperature)
    rad_table(5,1,i) = nint(radiance * 1000.) 
  end do
 
 end subroutine load_COMS_calibration


 ! Perform COMS Reflectance calculation

 subroutine COMS_REFLECTANCE_PRELAUNCH(COMS_Counts, alb_temp)
                              
    integer (kind=INT2), dimension(:,:), intent(in):: COMS_Counts
    real (kind=real4), dimension(:,:), intent(out):: alb_temp

    integer :: i, j
    integer :: index
    
    do j=1, Image%Number_Of_Lines_Read_This_Segment
      do i=1, Image%Number_Of_Elements
        if (Space_Mask(i,j) == sym%NO_SPACE) then
          
            !Not sure if I need to add 1 here.  Seems to be necessary
            !to match McIDAS-X.
            index = int(COMS_Counts(i,j),kind=int2) + 1
            
            ! until Mc-X is fixed, will use MTSAT HRIT scaling - WCS3
            !alb_temp(i,j) = 0.000978494701531316*(index) - 0.00197851402698724
            !alb_temp(i,j) = (alb_temp(i,j) * 100.0) 

            Alb_Temp(i,j) = Ref_Table(index)
  
            else

             Alb_Temp(i,j) = Missing_Value_Real4
            
            endif      
      end DO
    end DO
                
 end subroutine COMS_REFLECTANCE_PRELAUNCH

 ! Perform COMS Navigation

 subroutine COMS_NAVIGATION(xstart,ystart,xsize,ysize,xstride, &
                            AREAstr,NAVstr_COMS)
    integer(kind=int4) :: xstart, ystart
    integer(kind=int4) :: xsize, ysize
    integer(kind=int4) :: xstride  
    type (area_header_type) :: AREAstr
    type (GVAR_NAV), intent(in)    :: NAVstr_COMS
    
    integer :: i, j, ii, jj, imode
    real(kind(0.0d0)) :: latitude, longitude
    real(kind=real4) :: height
    integer :: FGF_type = 3 !COMS uses JMA GEOS navigation, so set type here

    NAVstr_COMS_NAV = NAVstr_COMS
    
    imode = -1
    height = 0.0   !Used for parallax correction
    
    Nav%Lat = Missing_Value_Real4
    Nav%Lon = Missing_Value_Real4
    
    if (NAVstr_COMS%nav_type == 'GEOS') then      
        !HRIT requires actual line and element of being processed.
        ! Unlike MSG, COMS requires no switching to different corrdinates.
                
          do j=1, ysize
            
            jj = ystart + (j-1) + (AREAstr%north_bound / real(AREAstr%line_res))
            
            
            !print *, AREAstr%north_bound / real(AREAstr%line_res)
            !stop
               
            do i=1, xsize
                ii = (i - 1)*(xstride) + xstart ! get element of the image segement
                ii = ii  + (AREAstr%west_vis_pixel / real(AREAstr%elem_res))

                ! again, use common algorithm for CGMS navigation
                !call pixcoord2geocoord_cgms(ii,                  &
                !                            jj,                  &
                !                            NAVstr_COMS%LOFF,   &
                !                            NAVstr_COMS%COFF,   & 
                !                            NAVstr_COMS%LFAC,   &
                !                            NAVstr_COMS%CFAC,   &
                !                            1,             &
                !                            NAVstr_COMS%sub_lon, &
                !                            latitude,            &
                !                            longitude)

                 ! Putting in hooks for updated navigation code
                 ! to be edited after 1/12

                CALL fgf_to_earth(FGF_type,                  &
                                  DBLE(ii),                  &
                                  DBLE(jj),                  &
                                  DBLE(NAVstr_COMS%CFAC),    &
                                  DBLE(NAVstr_COMS%COFF),    &
                                  DBLE(NAVstr_COMS%LFAC),    &
                                  DBLE(NAVstr_COMS%LOFF),    &
                                  DBLE(NAVstr_COMS%sub_lon), &
                                  longitude,                 &
                                  latitude)

             if (latitude .LE. -999.0) then  ! -999.99 is MSV nav missing value
                    Nav%Lat_1b(i,j) = Missing_Value_Real4
                    Nav%Lon_1b(i,j) = Missing_Value_Real4
                    Space_Mask(i,j) = sym%SPACE
                else
                    Nav%Lat_1b(i,j) = real(latitude,kind=real4)
                    Nav%Lon_1b(i,j) = real(longitude,kind=real4)
                    
                    ! BecaUSE JMA sets their longitudes from 0 to 360, and
                    ! we want 180 to -180, one last check.
                    
                    if (longitude .GT. 180.0 ) then
                        Nav%Lon_1b(i,j) = real(longitude,kind=real4) - 360.0
                    endif
                                        
                    Space_Mask(i,j) = sym%NO_SPACE
                endif

        
            end DO
                        
        end do     
        
    endif
      
 end subroutine COMS_NAVIGATION
 
 
!------------------------------------------------------------------
! subroutine to convert COMS counts to radiance and brightness
! temperature
!------------------------------------------------------------------
  
  subroutine COMS_RADIANCE_BT(chan_num,COMS_Counts, rad2, temp1)

    integer (kind=INT2), dimension(:,:), intent(in):: COMS_Counts
    integer (kind=int1), INTENT(in) :: chan_num
    real (kind=real4), dimension(:,:), INTENT(out):: temp1, rad2
    
    integer :: i, j, index
                                   
    do j = 1, Image%Number_Of_Lines_Read_This_Segment
      do i = 1, Image%Number_Of_Elements
        if (Space_Mask(i,j) == sym%NO_SPACE) then
          index = int(COMS_Counts(i,j),kind=int2) + 1
          rad2(i,j) = real(rad_table(chan_num,1,index),kind=real4)/1000.0
          temp1(i,j) = real(bt_table(chan_num,1,index),kind=real4)/100.0                    
        else
          rad2(i,j) = Missing_Value_Real4
          temp1(i,j) = Missing_Value_Real4
        endif
      end DO
    end do    
  
  end subroutine COMS_RADIANCE_BT
  
!------------------- COMS NAV BLOC

 subroutine READ_NAVIGATION_BLOCK_COMS(filename, AREAstr, NAVstr)
  CHARACTER(len=*), intent(in):: filename
  type(area_header_type), intent(in):: AREAstr
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

 end subroutine READ_NAVIGATION_BLOCK_COMS

!----------------------------------------------
! Copied subroutines
!----------------------------------------------
subroutine GET_COMS_IMAGE(filename,AREAstr, &
                                    segment_number, &
                                    num_lines_per_segment, &
                                    num_lines_read, image)

 character(len=*), intent(in):: filename 
 type (area_header_type), intent(in) :: AREAstr
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

 bytemove = 0 !COMS IS 0 byte offset

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
 
 allocate(word_buffer(number_of_words_in_segment))

 if( bytes_per_pixel == 2) then
    call mreadf_int(filename//CHAR(0), &
                 first_byte_in_segment,  &
                 bytes_per_pixel, &
                 number_of_words_in_segment, &
                 number_of_words_read,word_buffer)
 endif

 if( bytes_per_pixel == 1) then
    allocate(buffer1(number_of_words_in_segment))
 
    call mreadf_int(filename//CHAR(0), &
                 first_byte_in_segment,  &
                 bytes_per_pixel, &
                 number_of_words_in_segment, &
                 number_of_words_read,buffer1)
                 
    word_buffer = int(buffer1,kind=int2)
    !Since 1 byte values are always signed in Fortran, I need to add 256 to the negative values
    where (word_buffer < 0)
      word_buffer = word_buffer + 256
    end where
 endif

 !--- update number of scans read
 num_lines_read = number_of_words_read / words_per_line

 do Line_Idx = 1, num_lines_read
    word_start = (words_in_prefix + AREAstr%num_elem)*(Line_Idx-1) + words_in_prefix + 1
    word_end = min(word_start + AREAstr%num_elem - 1,number_of_words_in_segment)
    nwords = int(word_end - word_start)/COMS_Xstride + 1
    image(1:nwords,Line_Idx) = ishft(word_buffer(word_start:word_end:1),bytemove)       
 enddo

 deallocate(word_buffer)
 if (allocated(buffer1)) deallocate(buffer1)
 
 end subroutine GET_COMS_IMAGE

end module COMS_MODULE
