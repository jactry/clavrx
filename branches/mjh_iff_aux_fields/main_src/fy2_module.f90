!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: fy2_module.f90 (src)
!       fy2_module (program)
!
! PURPOSE: This module contains all the subroutines needed to perform navigation and
!          calibration for FY-2D/E
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
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
MODULE FY2_MODULE



use CONSTANTS
use PIXEL_COMMON
use CALIBRATION_CONSTANTS
use PLANCK
use NUMERICAL_ROUTINES
use GOES_MODULE
use FILE_UTILITY
use VIEWING_GEOMETRY_MODULE
  
 implicit none
 public :: FY_navigation, READ_FY
 public:: READ_FY_INSTR_CONSTANTS
 private :: FY_RADIANCE_BT,FY_Reflectance, MGIVSR


 TYPE (GVAR_NAV), PRIVATE    :: NAVstr_FY_NAV
 integer, PARAMETER, PRIVATE :: nchan_FY= 5
 integer, PARAMETER, PRIVATE :: ndet_FY = 4
 integer, PARAMETER, PRIVATE :: ntable_FY = 1024

 integer, PRIVATE :: nref_table_FY
 integer, PRIVATE :: nbt_table_FY
 CHARACTER(len=4), PRIVATE:: calib_type

 integer (kind=int4), dimension(nchan_FY,ndet_FY,ntable_FY), PRIVATE  :: ref_table
 integer (kind=int4), dimension(nchan_FY,ndet_FY,ntable_FY), PRIVATE  :: bt_table
 integer (kind=int4), dimension(nchan_FY,ndet_FY,ntable_FY), PRIVATE  :: rad_table

 integer(kind=int4), private, parameter:: FY_Xstride = 1
 integer(kind=int4), private, parameter:: num_4km_scans_fd = 3712
 integer(kind=int4), private, parameter:: num_4km_elem_fd = 3712
 integer(kind=int4), private, parameter:: time_for_fd_scan =  1560000 !milliseconds (26min)
 real, private, save:: Scan_rate    !scan rate in millsec / line
 integer(kind=int4), private, parameter:: FY2_Byte_Shift = 0 !number of bytes to shift for FY2
 
 CONTAINS
!----------------------------------------------------------------
! read the FY2 constants into memory
!-----------------------------------------------------------------
subroutine READ_FY_INSTR_CONSTANTS(Instr_Const_File)
 character(len=*), intent(in):: Instr_Const_File
 integer:: ios0, erstat
 integer:: Instr_Const_lun

 Instr_Const_lun = GET_LUN()

 open(unit=Instr_Const_lun,file=trim(Instr_Const_File),status="old",position="rewind",action="read",iostat=ios0)

 print *, "opening ", trim(Instr_Const_File)
 erstat = 0
 if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, "Error opening FY constants file, ios0 = ", ios0
    stop 19
 endif
  read(unit=Instr_Const_lun,fmt="(a3)") sat_name
  read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
  read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
  read(unit=Instr_Const_lun,fmt=*) planck_a1(20), planck_a2(20), planck_nu(20)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(27), planck_a2(27), planck_nu(27)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(31), planck_a2(31), planck_nu(31)
  read(unit=Instr_Const_lun,fmt=*) planck_a1(32), planck_a2(32), planck_nu(32)
  close(unit=Instr_Const_lun)

  !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / ew_Ch20

  !-- hardwire ch1 dark count
  Ch1_Dark_Count = 29

end subroutine READ_FY_INSTR_CONSTANTS
!============================================================================
! Perform fy Reflectance and BT calibration
!============================================================================
    subroutine READ_FY(segment_number,channel_1_filename, &
                     jday, image_time_ms, &
                     AREAstr,NAVstr_FY)

   integer(kind=int4), intent(in):: segment_number
   character(len=*), intent(in):: channel_1_filename
   TYPE (AREA_STRUCT), intent(in) :: AREAstr
   TYPE (GVAR_NAV), intent(in)    :: NAVstr_FY
   integer(kind=int2), intent(in):: jday
   integer(kind=int4), intent(in):: image_time_ms

   character(len=120):: channel_x_filename
   character(len=120):: channel_x_filename_full
   character(len=120):: channel_x_filename_full_uncompressed
   character(len=180):: System_String
   integer:: ipos
   integer:: ilen
   integer:: Elem_Idx
   integer:: Line_Idx
   integer:: ichan_goes
   integer:: ichan_modis
   integer:: fy_file_id
   real(kind=real4):: image_time_hours
   integer(kind=int4):: image_jday
   integer(kind=int4):: first_line_in_segment
   character(len=2):: ichan_goes_string
   integer:: num_elements_this_image
   integer:: num_scans_this_image

   !--- assume channel_1_file name has a unique "_1_" in the name. 
   !--- determine indices needed to replace that string
   ipos = index(channel_1_filename, "_1_")
   ilen = len(channel_1_filename)
    
   first_line_in_segment = (segment_number-1)*Image%Number_Of_Lines_Per_Segment

   !---------------------------------------------------------------------------
   ! FY Navigation (Do Navigation first)
   !---------------------------------------------------------------------------
   
   call FY_navigation(1,first_line_in_segment,&
                              Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment,1,&
                              AREAstr,NAVstr_FY)

   if (segment_number == 1) then

      image_jday = jday
      image_time_hours = image_time_ms / 60.0 / 60.0 / 1000.0

      !--- compute scan rate for future use *NOT CORRECT FOR FY RIGHT NOW*
      num_elements_this_image =  int(AREAstr%num_elem / FY_Xstride) + 1
      num_scans_this_image = AREAstr%num_line
      Scan_Rate = real((num_elements_this_image)/               &
               real(num_4km_elem_fd/FY_Xstride)) * &
               real((num_scans_this_image) / real(num_4km_scans_fd)) * &
               real(time_for_fd_scan) / real(num_scans_this_image)

       do ichan_goes = 2,5

       if (ichan_goes == 5) ichan_modis = 20
       if (ichan_goes == 4) ichan_modis = 27
       if (ichan_goes == 2) ichan_modis = 31
       if (ichan_goes == 3) ichan_modis = 32
       
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
    
    ! On first segment, reflectance, BT and rad tables from McIDAS Header
    fy_file_id = get_lun()   
    if(l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) THEN
      call MREAD_OPEN(trim(Temporary_Data_Dir)//trim(channel_1_filename)//CHAR(0), fy_file_id)
    else 
      call MREAD_OPEN(trim(Image%Level1b_Path)//trim(channel_1_filename)//CHAR(0), fy_file_id)
    endif  

    call LOAD_FY_CALIBRATION(fy_file_id, AREAstr)
    call MREAD_CLOSE(fy_file_id)

   endif
  
        
    IF(Sensor%Chan_On_Flag_Default(1) == sym%YES) THEN

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_1_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_1_filename)
       endif

        call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    FY2_Byte_Shift, &
                                    AREAstr, FY_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)
     
       
       call FY_Reflectance(Two_Byte_Temp,ch(1)%Ref_Toa(:,:))
       
       !--- store ch1 counts for support of PATMOS-x calibration studies
       Ch1_Counts = Two_Byte_Temp
       
    
    endif
    
        
    IF(Sensor%Chan_On_Flag_Default(20) == sym%YES) THEN

       channel_x_filename = channel_1_filename(1:ipos-1) // "_5_" // &
                            channel_1_filename(ipos+3:ilen)
       
       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif

       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    FY2_Byte_Shift, &
                                    AREAstr, FY_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)
     
       
       
       
      call FY_RADIANCE_BT(5_int1, Two_Byte_Temp, ch(20)%Rad_Toa, ch(20)%Bt_Toa)
        
    endif



    if (Sensor%Chan_On_Flag_Default(27) == sym%YES) THEN

       channel_x_filename = channel_1_filename(1:ipos-1) // "_4_" // &
                            channel_1_filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif

        call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    FY2_Byte_Shift, &
                                    AREAstr, FY_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)
     
       
      call FY_RADIANCE_BT(4_int1, Two_Byte_Temp, ch(27)%Rad_Toa, ch(27)%Bt_Toa)
           
    ENDif    
    
    IF(Sensor%Chan_On_Flag_Default(31) == sym%YES) THEN

       channel_x_filename = channel_1_filename(1:ipos-1) // "_2_" // &
                            channel_1_filename(ipos+3:ilen)
       
       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif

        call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    FY2_Byte_Shift, &
                                    AREAstr, FY_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)
     
       
         call FY_RADIANCE_BT(2_int1, Two_Byte_Temp, ch(31)%Rad_Toa, ch(31)%Bt_Toa)
         
    endif
    


    IF(Sensor%Chan_On_Flag_Default(32) == sym%YES) THEN
           

       channel_x_filename = channel_1_filename(1:ipos-1) // "_3_" // &
                            channel_1_filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               channel_x_filename_full = trim(Temporary_Data_Dir)//trim(channel_x_filename)
       else
               channel_x_filename_full = trim(Image%Level1b_Path)//trim(channel_x_filename)
       endif

        call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    FY2_Byte_Shift, &
                                    AREAstr, FY_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)
     
         call FY_RADIANCE_BT(3_int1, Two_Byte_Temp, ch(32)%Rad_Toa, ch(32)%Bt_Toa)
    
    endif



!to here

    
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     Scan_Number(Line_Idx) = first_line_in_segment + Line_Idx
     Scan_Time(Line_Idx) = image_time_ms + (Scan_Number(Line_Idx)-1) * Scan_rate
   enddo

!------------------------------------------------------------------------------
! FY Angles
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
                    
        !--- because FY2 has solar contamination, we will set the
        !    Solar_Contamination_Mask here as well for night pixels
        
        Solar_Contamination_Mask(Elem_Idx,Line_Idx) = sym%NO
        
!        IF ( (Geo%Solzen(Elem_Idx,Line_Idx) .GT. 90.) .AND. &
!             (ch(1)%Ref_Toa(Elem_Idx,Line_Idx) .GT. 0)) THEN
        IF ( (Geo%Solzen(Elem_Idx,Line_Idx) .GT. 90.)) THEN
             
             
             Solar_Contamination_Mask(Elem_Idx,Line_Idx) = sym%YES
             
        ENDIF
        
        
         
     enddo
      call COMPUTE_SATELLITE_ANGLES(Sensor%Geo_Sub_Satellite_Longitude,  &
                                    Sensor%Geo_Sub_Satellite_Latitude, Line_Idx)
   enddo
   
      
   !--- ascending node
   Elem_Idx = Image%Number_Of_Elements/2
   do Line_Idx = Line_Idx_Min_Segment+1, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     Nav%Ascend(Line_Idx) = 0
     if (Nav%Lat_1b(Elem_Idx,Line_Idx) < Nav%Lat_1b(Elem_Idx,Line_Idx-1)) then
       Nav%Ascend(Line_Idx) = 1
     endif
   enddo
   Nav%Ascend(Line_Idx_Min_Segment) = Nav%Ascend(Line_Idx_Min_Segment+1)
    
   
    
 end subroutine READ_FY


 ! Perform fy Reflectance calculation
 subroutine FY_Reflectance(FY_Counts, Alb_Temp)

    integer (kind=INT2), dimension(:,:), intent(in):: FY_Counts
    real (kind=real4), dimension(:,:), intent(out):: Alb_Temp

    integer :: i, j
    integer :: index

    DO j=1, Image%Number_Of_Lines_Read_This_Segment
      DO i=1, Image%Number_Of_Elements
        if (space_mask(i,j) == sym%NO_SPACE .and. Geo%solzen(i,j) < 90.0) THEN
         
            !Not sure if I need to add 1 here
            index = int(FY_Counts(i,j),kind=int2)

            alb_temp(i,j) = Missing_Value_Real4

            IF((index .GT. 0) .AND. (index .LE. 1024)) THEN
                        
                Alb_Temp(i,j) = &
                    (real(ref_table(2,1,index),kind=real4) / &
                     100.)
            ENDif                        
 
        ELSE
            Alb_Temp(i,j) = Missing_Value_Real4
        ENDif      
      END DO
    END DO

 end subroutine FY_Reflectance

 ! Perform FY-2 Navigation

 subroutine FY_navigation(xstart,ystart,xsize,ysize,xstride, &
                            AREAstr,NAVstr_FY)
    integer(kind=int4) :: xstart, ystart
    integer(kind=int4) :: xsize, ysize
    integer(kind=int4) :: xstride  
    type (AREA_STRUCT) :: AREAstr
    type (GVAR_NAV), intent(in)    :: NAVstr_FY
    
    integer :: i, j, ii, jj, ierr, imode
    real(kind=real4) :: elem, line, height
    real(kind=real4) :: dlon, dlat
    real(kind=real8) :: mjd
    real(kind=real4), dimension(8) :: angles

    NAVstr_FY_NAV = NAVstr_FY
    
    imode = -1
    height = 0.0     !Used for parallax correction
    
    Nav%Lat = Missing_Value_Real4
    Nav%Lon = Missing_Value_Real4
    
               
    if (NAVstr_FY%nav_type == 'GMSX') THEN
    
        jj = 1 + (ystart)*AREAstr%line_res
            
        DO j=1, ysize
            line = real(AREAstr%north_bound) + real(jj - 1) + &
                    real(AREAstr%line_res)/2.0
               
            DO i=1, xsize
                ii = ((i+(xstart-1)) - 1)*(AREAstr%elem_res*(xstride)) + 1
        
                elem = real(AREAstr%west_vis_pixel) + real(ii - 1) + &
	                   real(AREAstr%elem_res*(xstride))/2.0
                   
                call MGIVSR(imode,elem,line,dlon,dlat,height,&
                            angles,mjd,ierr)
                
                Space_Mask(i,j) = sym%SPACE
            
                if (ierr == 0) THEN
                     Space_Mask(i,j) = sym%NO_SPACE
                     Nav%Lat_1b(i,j) = dlat
                     Nav%Lon_1b(i,j) = dlon
                endif
                
            enddo
            
            jj = jj + AREAstr%line_res
        enddo
        
    endif
    
      
 end subroutine FY_NAVIGATION

!------------------------------------------------------------------
! subroutine to convert fy counts to radiance and brightness
! temperature
!------------------------------------------------------------------
  subroutine FY_RADIANCE_BT(Chan_Num,FY_Counts, rad2, temp1)

    integer (kind=INT2), dimension(:,:), intent(in):: FY_Counts
    integer (kind=int1), intent(in) :: chan_num
    real (kind=real4), DIMENSION(:,:), intent(out):: temp1, rad2
    
    integer :: i, j, index

    DO j = 1, Image%Number_Of_Lines_Read_This_Segment
      DO i = 1, Image%Number_Of_Elements
        index = int(FY_Counts(i,j),kind=int2) + 1
        
        if (space_mask(i,j) == sym%NO_SPACE  .AND. &
           (index <= 1024) .AND. (index >= 1)) THEN 
           
         rad2(i,j) = real(rad_table(chan_num,1,index),kind=real4)/1000.0
         temp1(i,j) = real(bt_table(chan_num,1,index),kind=real4)/100.0
        ELSE
         rad2(i,j) = Missing_Value_Real4
         temp1(i,j) = Missing_Value_Real4
        endif
      END DO
    END DO
    
  
  end subroutine FY_RADIANCE_BT
  


!---------------------------------------------------------------------
! Subroutine that loads the FY2 calibration table into the
! spacecraft structure.
!---------------------------------------------------------------------

subroutine load_fy_calibration(lun, AREAstr)
  integer(kind=int4), intent(in) :: lun
  type(AREA_STRUCT), intent(in):: AREAstr
  integer(kind=int4), dimension(6528) :: ibuf
  integer :: nref, nbt, i, j, offset
  integer(kind=int4) :: band_offset_2, band_offset_14, band_offset_15, &
                        band_offset_9, band_offset_7, dir_offset, avoid_warning
  real(kind=real4) :: albedo, temperature, radiance
  real(kind=real4), dimension(5)  :: a_fy, b_fy, nu_fy
    
  avoid_warning = lun
  
  ! Constants taken from:
  ! http://fengyunuds.cma.gov.cn/FYCV_EN/PublishInfo/PublishPage.aspx?SuperiorID=1&ChildCatalogID=%2010082-1
  ! GSICS calibrations
  
  if (AREAstr%sat_id_num == 36) THEN
  
    nu_fy(5) = 2601.9680
    nu_fy(4) = 1429.0728
    nu_fy(2) = 923.4366
    nu_fy(3) = 839.4041
  
    a_fy(5) = 0.9814
    a_fy(4) = 0.9874
    a_fy(2) = 0.9979
    a_fy(3) = 0.9980
  
    b_fy(5) = 2.98396
    b_fy(4) = 1.5094
    b_fy(2) = 0.3865
    b_fy(3) = 0.3710
  endif

  if (AREAstr%sat_id_num  == 37) THEN
  
    nu_fy(5) = 2568.2084
    nu_fy(4) = 1436.5964
    nu_fy(2) = 923.0511
    nu_fy(3) = 820.0376
  
    a_fy(5) = 0.9815
    a_fy(4) = 0.9883
    a_fy(2) = 0.9981
    a_fy(3) = 0.9986
  
    b_fy(5) = 2.9366
    b_fy(4) = 1.3981
    b_fy(2) = 0.3609
    b_fy(3) = 0.2661
  endif
     
  !---------------------------------------------------------------------
  ! Read from the calibration block.
  !---------------------------------------------------------------------
  
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
  nref_table_fy = nref
  nbt_table_fy = nbt
    
  !---------------------------------------------------------------------
  ! Load the visible channel calibration table for HiRID/MVIS calibration
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
  do i=1, 1024
    !  3.9 micron
    offset = band_offset_7/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(5,1,i) = nint(temperature * 100.)
       
    radiance = c1*(nu_fy(5)**3)/(exp((c2*nu_fy(5))/ &
               (a_fy(5)*temperature+b_fy(5)))-1.0)
    rad_table(5,1,i) = nint(radiance * 1000.)
    
    !  6.8 micron
    offset = band_offset_9/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(4,1,i) = nint(temperature * 100.)
    radiance = c1*(nu_fy(4)**3)/(exp((c2*nu_fy(4))/ &
               (a_fy(4)*temperature+b_fy(4)))-1.0)
    rad_table(4,1,i) = nint(radiance * 1000.)
    
    ! 10.8 micron
    offset = band_offset_14/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(2,1,i) = nint(temperature * 100.)
    radiance = c1*(nu_fy(2)**3)/(exp((c2*nu_fy(2))/ &
               (a_fy(2)*temperature+b_fy(2)))-1.0)
    rad_table(2,1,i) = nint(radiance * 1000.)
        
    ! 12.0 micron
    offset = band_offset_15/4
    temperature = real(ibuf(offset + i),kind=real4) / 1000.
    bt_table(3,1,i) = nint(temperature * 100.)
    radiance = c1*(nu_fy(3)**3)/(exp((c2*nu_fy(3))/ &
               (a_fy(3)*temperature+b_fy(3)))-1.0)
    rad_table(3,1,i) = nint(radiance * 1000.) 
  end do
 
end subroutine load_fy_calibration


!FY2 Navigation


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
!      integer*4  MAP(672,4)
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
       DTIMS = NAVstr_FY_NAV%DTIMS
       RESLIN = NAVstr_FY_NAV%RESLIN
       RESELM = NAVstr_FY_NAV%RESELM
       RLIC = NAVstr_FY_NAV%RLIC
       RELMFC = NAVstr_FY_NAV%RELMFC
       SENSSU = NAVstr_FY_NAV%SENSSU
       RLINE = NAVstr_FY_NAV%RLINE
       RELMNT = NAVstr_FY_NAV%RELMNT
       VMIS = NAVstr_FY_NAV%VMIS
       ELMIS = NAVstr_FY_NAV%ELMIS
       DSPIN = NAVstr_FY_NAV%DSPIN     
       
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
       IF(ABS(IMODE).GT.4)  IRTN=1
       IF(ABS(RLAT).GT.90. .AND. IMODE.GT.0)  IRTN=2
       IF(IRTN.NE.0)  return
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
       IF( IMODE.GT.0 .AND. IMODE.LT.5 )  THEN
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
         IF(TF.LT.0.D0)  TP=-TP
         call  MGI230(SP,SL,TL)
 !
         RI      = SNGL(HPAI-TL)/RSTEP+RFCL-VMIS(2)/RSTEP
         RJ      = SNGL(TP)/RSAMP+RFCP &
                  +VMIS(3)/RSAMP-SNGL(HPAI-TL)*TAN(VMIS(1))/RSAMP
 !
         IF(ABS(RI-RI0).GE.EPS)  THEN                              ![3.8]
           RTIM  = DBLE(AINT((RI-1.)/SENS)+RJ*RSAMP/SNGL(DPAI))/ &
                   (DSPIN*1440.D0)+DTIMS
           RI0   = RI
           GO TO  100
         endif
         RLIN    = RI
         RPIX    = RJ
         DSCT    = RTIM
         IF(RLIN.LT.0 .OR. RLIN.GT.RFTL)  IRTN=4
         IF(RPIX.LT.0 .OR. RPIX.GT.RFTP)  IRTN=5
 !
 !!!!!!!!!!!!!!!!!! TRANSFORMATION (VISSR=>GEOGRAPHICAL) !!!!!!!!!!!!!!!!
       ELSEIF(IMODE.LT.0 .AND. IMODE.GT.-5)  THEN
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
         IF(DD.GE.0.D0 .AND. DDA.NE.0.D0)  THEN
           DK1     = (-DDB+DSQRT(DD))/DDA
           DK2     = (-DDB-DSQRT(DD))/DDA
         ELSE
           IRTN    = 6
           GO TO  9000
         endif
         IF(DABS(DK1).LE.DABS(DK2))  THEN
           DK    = DK1
         ELSE
           DK    = DK2
         endif
         STN1(1) = SAT(1)+DK*SL(1)
         STN1(2) = SAT(2)+DK*SL(2)
         STN1(3) = SAT(3)+DK*SL(3)
         DLAT    = DATAN(STN1(3)/(DEF*DSQRT(STN1(1)*STN1(1)+ &
                   STN1(2)*STN1(2))))                              ![3.15]
         IF(STN1(1).NE.0.D0)  THEN
           DLON  = DATAN(STN1(2)/STN1(1))
           IF(STN1(1).LT.0.D0 .AND. STN1(2).GE.0.D0)  DLON=DLON+PI
           IF(STN1(1).LT.0.D0 .AND. STN1(2).LT.0.D0)  DLON=DLON-PI
         ELSE
           IF(STN1(2).GT.0.D0)  THEN
             DLON=HPAI
           ELSE
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
       IF(DSATZ.GT.HPAI)  IRTN = 7
 !
       SUNM    = 315.253D0+0.985600D0*RTIM                         ![3.20]
       SUNM    = DMOD(SUNM,360.D0)*CDR
       SDIS    = (1.00014D0-0.01672D0*DCOS(SUNM)-0.00014*DCOS(2.D0* &
                 SUNM))*1.49597870D8
 !
       IF(DLAT.GE.0.D0) THEN                                       ![3.21]
         DLATN   = HPAI-DLAT
         DLONN = DLON-PI
         IF(DLONN.LE.-PI)  DLONN=DLONN+DPAI
       ELSE
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
 !!!!!!!!!!!!!!!!!!! STOP/END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  9000 CONTINUE
       return
end subroutine MGIVSR

subroutine  MGI100(RTIM,CDR,SAT,SP,SS,BETA)
!       COMMON /MMAP1/MAP
       real*8    ATTALP,ATTDEL,BETA,CDR,DELT,RTIM,SITAGT,SUNALP,SUNDEL, &
                 WKCOS,WKSIN
       real*8    ATT1(3),ATT2(3),ATT3(3),NPA(3,3), &
                 SAT(3),SP(3),SS(3)
 !     integer*4 MAP(672,4)
       integer :: I
 !
 !      EQUIVALENCE (MAP(13,3),ORBT1(1,1))
 !      EQUIVALENCE (MAP(13,2),ATIT(1,1))
 !
       DO 1000 I=1,7
         IF(RTIM.GE.NAVstr_FY_NAV%ORBT1(1,I).AND.RTIM.LT.NAVstr_FY_NAV%ORBT1(1,I+1))  THEN
!          call  MGI110 &
!               (I,RTIM,CDR,NAVstr_FY_NAV%ORBT1,ORBT2,SAT,SITAGT,SUNALP,SUNDEL,NPA)
           call  MGI110 &
                (I,RTIM,CDR,NAVstr_FY_NAV%ORBT1,SAT,SITAGT,SUNALP,SUNDEL,NPA)
           GO TO  1200
         endif
  1000 CONTINUE
  1200 CONTINUE
 !
       DO 3000 I=1,9
         IF(RTIM.GE.NAVstr_FY_NAV%ATIT(1,I) .AND. RTIM.LT.NAVstr_FY_NAV%ATIT(1,I+1))  THEN
           DELT = (RTIM-NAVstr_FY_NAV%ATIT(1,I))/(NAVstr_FY_NAV%ATIT(1,I+1)-NAVstr_FY_NAV%ATIT(1,I))
           ATTALP = NAVstr_FY_NAV%ATIT(3,I)+(NAVstr_FY_NAV%ATIT(3,I+1)-NAVstr_FY_NAV%ATIT(3,I))*DELT
           ATTDEL = NAVstr_FY_NAV%ATIT(4,I)+(NAVstr_FY_NAV%ATIT(4,I+1)-NAVstr_FY_NAV%ATIT(4,I))*DELT
           BETA   = NAVstr_FY_NAV%ATIT(5,I)+(NAVstr_FY_NAV%ATIT(5,I+1)-NAVstr_FY_NAV%ATIT(5,I))*DELT
           IF( (NAVstr_FY_NAV%ATIT(5,I+1)-NAVstr_FY_NAV%ATIT(5,I)).GT.0.D0 ) &
             BETA   = NAVstr_FY_NAV%ATIT(5,I)+(NAVstr_FY_NAV%ATIT(5,I+1)-NAVstr_FY_NAV%ATIT(5,I)-360.D0*CDR)*DELT
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

!subroutine MGI110(I,RTIM,CDR,ORBTA,ORBTB,SAT,SITAGT,SUNALP,SUNDEL,NPA)
!      real*8    CDR,SAT(3),RTIM,ORBTA(35,8),ORBTB(35,8)
subroutine MGI110(I,RTIM,CDR,ORBTA,SAT,SITAGT,SUNALP,SUNDEL,NPA)
       real*8    CDR,SAT(3),RTIM,ORBTA(35,8)
       real*8    SITAGT,SUNDEL,SUNALP,NPA(3,3),DELT
       integer*4 I
       IF(I.NE.8)  THEN
         DELT=(RTIM-ORBTA(1,I))/(ORBTA(1,I+1)-ORBTA(1,I))
         SAT(1)   = ORBTA( 9,I)+(ORBTA( 9,I+1)-ORBTA( 9,I))*DELT
         SAT(2)   = ORBTA(10,I)+(ORBTA(10,I+1)-ORBTA(10,I))*DELT
         SAT(3)   = ORBTA(11,I)+(ORBTA(11,I+1)-ORBTA(11,I))*DELT
         SITAGT   =(ORBTA(15,I)+(ORBTA(15,I+1)-ORBTA(15,I))*DELT)*CDR
         IF( (ORBTA(15,I+1)-ORBTA(15,I)).LT.0.D0 ) &
           SITAGT   =(ORBTA(15,I)+(ORBTA(15,I+1)-ORBTA(15,I)+360.D0) &
                     *DELT)*CDR
         SUNALP   =(ORBTA(18,I)+(ORBTA(18,I+1)-ORBTA(18,I))*DELT)*CDR
         IF( (ORBTA(18,I+1)-ORBTA(18,I)).GT.0.D0 ) &
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
       IF(RV1 == 0.D0)  return
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
       IF(AS2 == 0.D0)  return
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
       IF(DNAI.GT.0.D0)  AZI=DPAI-AZI
       return
end subroutine MGI240





END MODULE FY2_MODULE
