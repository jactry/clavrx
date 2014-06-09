!// $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: coms_module.f90 (src)
!       COMS_MODULE (program)
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

MODULE COMS_MODULE

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
 public :: READ_COMS
 public:: READ_NAVIGATION_BLOCK_COMS
 public:: GET_COMS_IMAGE
 public:: READ_COMS_INSTR_CONSTANTS
 public:: ASSIGN_COMS_ID_NUM_INTERNAL
         
 private :: COMS_RADIANCE_BT,COMS_Reflectance, &
            COMS_navigation
 

 TYPE (GVAR_NAV), PRIVATE    :: NAVstr_COMS_NAV
 integer, PARAMETER, PRIVATE :: nchan_COMS= 5
 INTEGER, PARAMETER, PRIVATE :: ndet_COMS = 4
 INTEGER, PARAMETER, PRIVATE :: ntable_COMS = 1024

 INTEGER, PRIVATE :: nref_table_COMS
 INTEGER, PRIVATE :: nbt_table_COMS
 CHARACTER(len=4), PRIVATE:: calib_type

 !---stw INTEGER (kind=int4), dimension(nchan_COMS,ndet_COMS,ntable_COMS), PRIVATE  :: ref_table
 REAL (kind=real4), dimension(ntable_COMS), PRIVATE  :: ref_table
 INTEGER (kind=int4), dimension(nchan_COMS,ndet_COMS,ntable_COMS), PRIVATE  :: bt_table
 INTEGER (kind=int4), dimension(nchan_COMS,ndet_COMS,ntable_COMS), PRIVATE  :: rad_table

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
!--------------------------------------------------------------------
! assign internal sat id's and const file names for COMS
!--------------------------------------------------------------------
subroutine ASSIGN_COMS_ID_NUM_INTERNAL(Mcidas_Id_Num)
    integer(kind=int4), intent(in):: Mcidas_Id_Num

    !COMS-1
    if (Mcidas_Id_Num == 250)   then
        Sc_Id_WMO = 810
        Instr_Const_file = 'coms1_instr.dat'
        Algo_Const_file = 'coms1_algo.dat'
        Platform_Name_Attribute = 'COMS-1'
        Sensor_Name_Attribute = 'IMAGER'
    endif

   Instr_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Instr_Const_file)
   Algo_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Algo_Const_file)

end subroutine ASSIGN_COMS_ID_NUM_INTERNAL
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
  read(unit=Instr_Const_lun,fmt=*) a1_20, a2_20,nu_20
  read(unit=Instr_Const_lun,fmt=*) a1_27, a2_27,nu_27
  read(unit=Instr_Const_lun,fmt=*) a1_31, a2_31,nu_31
  read(unit=Instr_Const_lun,fmt=*) a1_32, a2_32,nu_32

  read(unit=Instr_Const_lun,fmt=*) b1_day_mask,b2_day_mask,b3_day_mask,b4_day_mask
  close(unit=Instr_Const_lun)

  !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / ew_Ch20

  !-- hardwire ch1 dark count
  Ch1_Dark_Count = 29


end subroutine READ_COMS_INSTR_CONSTANTS

 ! Perform COMS Reflectance and BT calibration
 SUBROUTINE READ_COMS(segment_number,channel_1_filename, &
                     jday, image_time_ms, Time_Since_Launch, &
                     AREAstr,NAVstr_COMS)

   integer(kind=int4), intent(in):: segment_number
   character(len=*), intent(in):: channel_1_filename
   TYPE (AREA_STRUCT), intent(in) :: AREAstr
   TYPE (GVAR_NAV), intent(in)    :: NAVstr_COMS
   integer(kind=int2), intent(in):: jday
   integer(kind=int4), intent(in):: image_time_ms
   real(kind=real4), intent(in):: Time_Since_Launch

   character(len=120):: channel_x_filename
   character(len=120):: channel_x_filename_full
   character(len=120):: channel_x_filename_full_uncompressed
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
   INTEGER :: Line_Idx
   INTEGER :: Elem_Idx
   integer:: num_elements_this_image
   integer:: num_scans_this_image

   

   !--- assume channel_1_file name has a unique "_1_" in the name. 
   !--- determine indices needed to replace that string
   ipos = index(channel_1_filename, "_1_")
   ilen = len(channel_1_filename)
    
   first_line_in_segment = (segment_number-1)*num_scans_per_segment

   !---------------------------------------------------------------------------
   ! COMS Navigation (Do Navigation and Solar angles first)
   !---------------------------------------------------------------------------
   
   call COMS_navigation(1,first_line_in_segment,&
                              num_pix,Num_Scans_Per_Segment,1,&
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

       if (Chan_On_Flag_Default(ichan_modis) == sym%YES) then

          channel_x_filename = channel_1_filename(1:ipos-1) // "_"//trim(ichan_goes_string)//"_" // &
                            channel_1_filename(ipos+3:ilen)
          
          if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
          else
               channel_x_filename_full = trim(Dir_1b)//trim(channel_x_filename)
          endif

          channel_x_filename_full_uncompressed = trim(Dir_1b)//trim(channel_x_filename)
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
    IF (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) THEN
        CALL mread_open(trim(Temporary_Data_Dir)//trim(channel_1_filename)//CHAR(0), COMS_file_id)
    ELSE
        CALL mread_open(trim(Dir_1b)//trim(channel_1_filename)//CHAR(0), COMS_file_id)
    ENDIF  

    CALL load_COMS_calibration(COMS_file_id, AREAstr)
    CALL mread_close(COMS_file_id)


   endif
   
   
   !---   read channel 1 (COMS channel 1)
   if (Chan_On_Flag_Default(1) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_1_filename)
       else
               channel_x_filename_full = trim(Dir_1b)//trim(channel_1_filename)
       endif

       call GET_COMS_IMAGE(trim(channel_x_filename_full), &
                                    AREAstr, &
                                    segment_number, &
                                    num_scans_per_segment, &
                                    num_scans_read,   &
                                    Two_Byte_Temp)
                                    
        call COMS_Reflectance(Two_Byte_Temp,ch(1)%Ref_Toa(:,:))

        Ch1_Counts = Two_Byte_Temp

   endif
   
   !---   read channel 20 (COMS channel 2)
   if (Chan_On_Flag_Default(20) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_2_" // &
                            channel_1_filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Dir_1b)//trim(channel_x_filename)
       endif

       call GET_COMS_IMAGE(trim(channel_x_filename_full), &
                                    AREAstr, &
                                    segment_number, &
                                    num_scans_per_segment, &
                                    num_scans_read,   &
                                    Two_Byte_Temp)
      call COMS_RADIANCE_BT(2_int1, Two_Byte_Temp, ch(20)%Rad_Toa, ch(20)%Bt_Toa)

   endif
                    
   
   !---   read channel 27 (COMS channel 3)
   if (Chan_On_Flag_Default(27) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_3_" // &
                            channel_1_filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Dir_1b)//trim(channel_x_filename)
       endif

       call GET_COMS_IMAGE(trim(channel_x_filename_full), &
                                    AREAstr, &
                                    segment_number, &
                                    num_scans_per_segment, &
                                    num_scans_read,   &
                                    Two_Byte_Temp)
      call COMS_RADIANCE_BT(3_int1, Two_Byte_Temp, ch(27)%Rad_Toa, ch(27)%Bt_Toa)

   endif


   
   !---   read channel 31 (COMS channel 4)
   if (Chan_On_Flag_Default(31) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_4_" // &
                            channel_1_filename(ipos+3:ilen)
       
       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Dir_1b)//trim(channel_x_filename)
       endif
       
       call GET_COMS_IMAGE(trim(channel_x_filename_full), &
                                    AREAstr, &
                                    segment_number, &
                                    num_scans_per_segment, &
                                    num_scans_read,   &
                                    Two_Byte_Temp)
      call COMS_RADIANCE_BT(4_int1, Two_Byte_Temp, ch(31)%Rad_Toa, ch(31)%Bt_Toa)

   endif
   
   
   !---   read channel 32 (COMS channel 5)
   if (Chan_On_Flag_Default(32) == sym%YES) then

       channel_x_filename = channel_1_filename(1:ipos-1) // "_5_" // &
                            channel_1_filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Dir_1b)//trim(channel_x_filename)
       endif

       call GET_COMS_IMAGE(trim(channel_x_filename_full), &
                                    AREAstr, &
                                    segment_number, &
                                    num_scans_per_segment, &
                                    num_scans_read,   &
                                    Two_Byte_Temp)
      call COMS_RADIANCE_BT(5_int1, Two_Byte_Temp, ch(32)%Rad_Toa, ch(32)%Bt_Toa)

   endif
    
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + num_scans_read - 1
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
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + num_scans_read - 1
     do Elem_Idx = 1,num_pix
        call POSSOL(image_jday,image_time_hours, &
                    Lon_1b(Elem_Idx,Line_Idx),Lat_1b(Elem_Idx,Line_Idx), &
                    Solzen(Elem_Idx,Line_Idx),Solaz(Elem_Idx,Line_Idx))
     enddo
      call COMPUTE_SATELLITE_ANGLES(goes_sub_satellite_longitude,  &
                      goes_sub_satellite_latitude, Line_Idx)
   enddo

   !--- ascending node
   Elem_Idx = num_pix/2
   do Line_Idx = Line_Idx_Min_Segment+1, Line_Idx_Min_Segment + num_scans_read - 1
     ascend(Line_Idx) = 0
     if (lat_1b(Elem_Idx,Line_Idx) < lat_1b(Elem_Idx,Line_Idx-1)) then
       ascend(Line_Idx) = 1
     endif
   enddo
   ascend(Line_Idx_Min_Segment) = ascend(Line_Idx_Min_Segment+1)

    
 END SUBROUTINE READ_COMS
 
 
 
 SUBROUTINE load_COMS_calibration(lun, AREAstr)
  INTEGER(kind=int4), intent(in) :: lun
  type(AREA_STRUCT), intent(in):: AREAstr
  INTEGER(kind=int4), dimension(6528) :: ibuf
  CHARACTER(len=25) :: cbuf
  INTEGER :: nref, nbt, i, offset
  INTEGER(kind=int4) :: band_offset_2, band_offset_14, band_offset_15, &
                        band_offset_9, band_offset_7, dir_offset
  REAL(kind=real4) :: albedo, temperature, radiance

  call mreadf_int_o(lun,AREAstr%cal_offset,4,6528,ibuf)
  !if (AREAstr%swap_bytes > 0) call swap_bytes4(ibuf,6528)
  
  dir_offset = ibuf(4)
  band_offset_2 = ibuf(6)
  band_offset_7 = ibuf(8)
  band_offset_9 = ibuf(10)
  band_offset_14 = ibuf(12)
  band_offset_15 = ibuf(14)  
  
  nref = 256
  nbt = 1024
  nref_table_COMS = nref
  nbt_table_COMS = nbt
  
  ! We need to get the vis calibration type. This is in 5th
  ! 4 character long block of the cal block (as per line 284, kbxmtst.dlm, v1.5)
  ! WCS3 (3/29/2010)
  
  call mreadf_int_o(lun,AREAstr%cal_offset,1,25,cbuf)
  calib_type = cbuf(17:20)

  !---------------------------------------------------------------------
  ! Load the visible channel calibration from the McIDAS AREA file
  !---------------------------------------------------------------------
  
  ptr = band_offset_2/4

  do i = 1, 1024

    !---stw At this point, it still needs a divide by 10 to
    !---stw to match McIDAS-X.

    albedo = real(ibuf(ptr),kind=real4) / 100.0
    !print*,"albedo : ", i,albedo

    ref_table(i) = albedo / 10.0
    !print *, i, ref_table(i)

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
 
 END SUBROUTINE load_COMS_calibration


 ! Perform COMS Reflectance calculation

 SUBROUTINE COMS_Reflectance(COMS_Counts, alb_temp)
                              
    INTEGER (kind=INT2), dimension(:,:), intent(in):: COMS_Counts
    REAL (KIND=real4), dimension(:,:), intent(out):: alb_temp

    INTEGER :: i, j
    INTEGER :: index
    
    DO j=1, num_scans_read
      DO i=1, num_pix
        IF (Space_Mask(i,j) == sym%NO_SPACE) THEN
          
            !Not sure if I need to add 1 here.  Seems to be necessary
            !to match McIDAS-X.
            index = int(COMS_Counts(i,j),KIND=int2) + 1
            
            ! until Mc-X is fixed, will use MTSAT HRIT scaling - WCS3
            !alb_temp(i,j) = 0.000978494701531316*(index) - 0.00197851402698724
            !alb_temp(i,j) = (alb_temp(i,j) * 100.0) 

            alb_temp(i,j) = ref_table(index)
  
            ELSE

             alb_temp(i,j) = Missing_Value_Real4
            
            ENDIF      
      END DO
    END DO
                
 END SUBROUTINE COMS_Reflectance

 ! Perform COMS Navigation

 SUBROUTINE COMS_navigation(xstart,ystart,xsize,ysize,xstride, &
                            AREAstr,NAVstr_COMS)
    INTEGER(KIND=int4) :: xstart, ystart
    INTEGER(KIND=int4) :: xsize, ysize
    INTEGER(KIND=int4) :: xstride  
    type (AREA_STRUCT) :: AREAstr
    TYPE (GVAR_NAV), intent(in)    :: NAVstr_COMS
    
    INTEGER :: i, j, ii, jj, imode
    REAL(KIND(0.0d0)) :: latitude, longitude
    REAL(KIND=REAL4) :: height

    NAVstr_COMS_NAV = NAVstr_COMS
    
    imode = -1
    height = 0.0   !Used for parallax correction
    
    lat = Missing_Value_Real4
    lon = Missing_Value_Real4
    
    IF (NAVstr_COMS%nav_type == 'GEOS') THEN      
        !HRIT requires actual line and element of being processed.
        ! Unlike MSG, COMS requires no switching to different corrdinates.
                
          DO j=1, ysize
            
            jj = ystart + (j-1) + (AREAstr%north_bound / real(AREAstr%line_res))
            
            
            !print *, AREAstr%north_bound / real(AREAstr%line_res)
            !stop
               
            DO i=1, xsize
                ii = (i - 1)*(xstride) + xstart ! get element of the image segement
                ii = ii  + (AREAstr%west_vis_pixel / real(AREAstr%elem_res))

                ! again, use common algorithm for CGMS navigation
                CALL pixcoord2geocoord_cgms(ii,                  &
                                            jj,                  &
                                            NAVstr_COMS%LOFF,   &
                                            NAVstr_COMS%COFF,   & 
                                            NAVstr_COMS%LFAC,   &
                                            NAVstr_COMS%CFAC,   &
                                            1,             &
                                            NAVstr_COMS%sub_lon, &
                                            latitude,            &
                                            longitude)
                                          
             IF (latitude .LE. -999.0) THEN  ! -999.99 is MSV nav missing value
                    Lat_1b(i,j) = Missing_Value_Real4
                    Lon_1b(i,j) = Missing_Value_Real4
                    Space_Mask(i,j) = sym%SPACE
                ELSE
                    Lat_1b(i,j) = REAL(latitude,kind=REAL4)
                    Lon_1b(i,j) = REAL(longitude,kind=REAL4)
                    
                    ! BecaUSE JMA sets their longitudes from 0 to 360, and
                    ! we want 180 to -180, one last check.
                    
                    IF (longitude .GT. 180.0 ) THEN
                        Lon_1b(i,j) = REAL(longitude,kind=REAL4) - 360.0
                    ENDIF
                                        
                    Space_Mask(i,j) = sym%NO_SPACE
                ENDIF

        
            END DO
                        
        END DO     
        
    ENDIF
      
 END SUBROUTINE COMS_navigation
 
 
!------------------------------------------------------------------
! SUBROUTINE to convert COMS counts to radiance and brightness
! temperature
!------------------------------------------------------------------
  
  SUBROUTINE COMS_RADIANCE_BT(chan_num,COMS_Counts, rad2, temp1)

    INTEGER (kind=INT2), dimension(:,:), intent(in):: COMS_Counts
    INTEGER (kind=int1), INTENT(in) :: chan_num
    REAL (kind=real4), DIMENSION(:,:), INTENT(out):: temp1, rad2
    
    INTEGER :: i, j, index
                                   
    DO j = 1, num_scans_read
      DO i = 1, num_pix
        IF (Space_Mask(i,j) == sym%NO_SPACE) THEN
          index = int(COMS_Counts(i,j),KIND=int2) + 1

          rad2(i,j) = REAL(rad_table(chan_num,1,index),KIND=REAL4)/1000.0
          temp1(i,j) = REAL(bt_table(chan_num,1,index),KIND=REAL4)/100.0                    
        ELSE
          rad2(i,j) = Missing_Value_Real4
          temp1(i,j) = Missing_Value_Real4
        ENDIF
      END DO
    END DO    
  
  END SUBROUTINE COMS_RADIANCE_BT
  
!------------------- COMS NAV BLOC

 SUBROUTINE READ_NAVIGATION_BLOCK_COMS(filename, AREAstr, NAVstr)
  CHARACTER(len=*), intent(in):: filename
  TYPE(AREA_STRUCT), intent(in):: AREAstr
  TYPE(GVAR_NAV), intent(inout):: NAVstr
 
  INTEGER :: geos_nav
  INTEGER(kind=int4)nav_offset
  integer:: number_of_words_read
  INTEGER(kind=int4), DIMENSION(640) :: i4buf
  
  nav_offset = AREAstr%sec_key_nav
    
  !determine GEOS or other navigation
  geos_nav = sym%NO
  CALL mreadf_int(trim(filename)//CHAR(0),nav_offset,4,640,&
                    number_of_words_read, i4buf)
!  IF (AREAstr%swap_bytes > 0) CALL swap_bytes4(i4buf,640)
  CALL move_bytes(4,i4buf(1),NAVstr%nav_type,0)
 !SUBLON stored as SUBLON *10 in McIDAS NAV block
  NAVstr%sub_lon = REAL(i4buf(6),kind=real4) / 10
  NAVstr%sublon = NAVstr%sub_lon

	 ! LOFF, COFF, CFAC, LFAC stored in McIDAS header for 1km data. All
	 ! Multipied by 10. Order from nvxmtst.dlm in McIDAS
  NAVstr%LOFF=(i4buf(2) / 10 ) / REAL(AREAstr%line_res)
  NAVstr%COFF=(i4buf(3) / 10) / REAL(AREAstr%elem_res)
  NAVstr%LFAC=(i4buf(4) / 10 ) / REAL(AREAstr%line_res)
  NAVstr%CFAC=(i4buf(5) / 10 ) / REAL(AREAstr%elem_res)

 END SUBROUTINE READ_NAVIGATION_BLOCK_COMS

!----------------------------------------------
! Copied subroutines
!----------------------------------------------


!COMS IS TRICKY. Need to keep seperate because of wierdness in HiRID bytes per pixel

subroutine GET_COMS_IMAGE(filename,AREAstr, &
                                    segment_number, &
                                    num_lines_per_segment, &
                                    num_lines_read, image)

 character(len=*), intent(in):: filename 
 TYPE (AREA_STRUCT), intent(in) :: AREAstr
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
 INTEGER (kind=int1), dimension(:), allocatable :: buffer1
 INTEGER(kind=int4), DIMENSION(64) :: i4buf_temp
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

END MODULE COMS_MODULE
