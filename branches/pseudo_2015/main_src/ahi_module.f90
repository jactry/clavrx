!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: AHI_module.f90 (src)
!       AHI_MODULE (program)
!
! PURPOSE: 
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
! AHI Channel Mapping
!
! wvl       ahi  modis/clavrx 
!
! 0.47       1      3   
! 0.51       2      4
! 0.64       3      1  
! 0.86       4      2  
! 1.6        5      6
! 2.2        6     r7
! 3.9        7     20
! 6.2        8     37
! 6.9        9     27
! 7.3       10     28
! 8.6       11     29
! 9.6       12     30
! 10.4      13     38
! 11.2      14     31
! 12.3      15     32 
! 13.3      16     33
!
!--------------------------------------------------------------------------------------
module AHI_MODULE

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

 private
 private:: READ_AHI_INSTR_CONSTANTS
 public:: READ_AHI
 public:: READ_AHI_CHANNEL
 public:: READ_AHI_INSTR_CONSTANTS
 public :: ASSIGN_AHI_SAT_ID_NUM_INTERNAL


 integer(kind=int2), parameter, private:: fill_value = 32767
 real(kind=real4), parameter, private:: missing_value = -999.0
 character(len=13), parameter:: MODULE_PROMPT="AHI_MODULE:"
 type (GVAR_NAV), PRIVATE    :: NAVstr_AHI_NAV
 integer, PARAMETER, PRIVATE :: nchan_ahi= 16

 integer(kind=int4), private, parameter:: Num_2km_Scans_Fd = 3712
 integer(kind=int4), private, parameter:: Num_2km_elem_fd = 3712
 integer, dimension(16), parameter, private:: Chan_Idx=[1,2,6,20,27,28,29,30,31,32,33]


contains


!--------------------------------------------------------------------
! assign internal sat id's and const file names for AHI
!--------------------------------------------------------------------
subroutine ASSIGN_AHI_SAT_ID_NUM_INTERNAL(Mcidas_Id_Num)
    integer(kind=int4), intent(in):: Mcidas_Id_Num



    if (Mcidas_Id_Num == 220)   then ! from MUG
        Sc_Id_WMO = 173
        Instr_Const_file = 'Himawari8_instr.dat'
        Algo_Const_file = 'Himawari8_algo.dat'
        Platform_Name_Attribute = 'Himawari-8'
        Sensor_Name_Attribute = 'IMAGER'
    endif

    if (Mcidas_Id_Num == 222)   then  ! Fudged for now
        Sc_Id_WMO = 174
        Instr_Const_file = 'Himawari9_instr.dat'
        Algo_Const_file = 'Himawari9_algo.dat'
        Platform_Name_Attribute = 'Himawari-8'
        Sensor_Name_Attribute = 'IMAGER'
    endif

   Instr_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Instr_Const_file)
   Algo_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Algo_Const_file)

end subroutine ASSIGN_AHI_SAT_ID_NUM_INTERNAL

!----------------------------------------------------------------
! read the AHI constants into memory
!-----------------------------------------------------------------
subroutine READ_AHI_INSTR_CONSTANTS(Instr_Const_file)
 character(len=*), intent(in):: Instr_Const_file
 integer:: ios0, erstat
 integer:: Instr_Const_lun

 Instr_Const_lun = GET_LUN()

 open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

 print *, EXE_PROMPT, MODULE_PROMPT, " Opening ", trim(Instr_Const_file)
 erstat = 0
 if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, MODULE_PROMPT, "Error opening AHI constants file, ios0 = ", ios0
    stop 19
 endif

  read(unit=Instr_Const_lun,fmt="(a3)") sat_name
  read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
  read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
  read(unit=Instr_Const_lun,fmt=*) a1_20, a2_20,nu_20 ! Band 7
  !Note AHI has a 6.2 (Band 8), but MODIS doesn't have one
  read(unit=Instr_Const_lun,fmt=*) a1_27, a2_27,nu_27 !Band 9
  read(unit=Instr_Const_lun,fmt=*) a1_28, a2_28,nu_28 !Band 10
  read(unit=Instr_Const_lun,fmt=*) a1_29, a2_29,nu_29 !Band 11
  read(unit=Instr_Const_lun,fmt=*) a1_30, a2_30,nu_30 !Band 12
  !NOTE AHI as a 10.4 (Band 13), but MODIS doesn't have one
  read(unit=Instr_Const_lun,fmt=*) a1_31, a2_31,nu_31 !Band 14
  read(unit=Instr_Const_lun,fmt=*) a1_32, a2_32,nu_32 !Band 15
  read(unit=Instr_Const_lun,fmt=*) a1_33, a2_33,nu_33 !Band 16
  read(unit=Instr_Const_lun,fmt=*) b1_day_mask,b2_day_mask,b3_day_mask,b4_day_mask
  close(unit=Instr_Const_lun)

  !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

end subroutine READ_AHI_INSTR_CONSTANTS


 ! Perform AHI Reflectance and BT calibration
 subroutine READ_AHI(segment_number,Channel_1_Filename, &
                     jday, image_time_ms, Time_Since_Launch, &
                     AREAstr,NAVstr_AHI)

   integer(kind=int4), intent(in):: segment_number
   character(len=*), intent(in):: Channel_1_Filename
   type (AREA_STRUCT), intent(in) :: AREAstr
   type (GVAR_NAV), intent(in)    :: NAVstr_AHI
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
   integer:: Chan_Idx_AHI
   integer:: Chan_Idx_Modis
   integer:: AHI_file_id
   real(kind=real4):: image_time_hours
   integer(kind=int4):: image_jday
   integer(kind=int4):: first_line_in_segment
   character(len=2):: Chan_Idx_AHI_String
   integer:: Num_Elements_This_Image
   integer:: Num_Scans_This_Image
   integer:: Line_Idx

   

   !--- assume Channel_1_file name has a unique "_1_" in the name. 
   !--- determine indices needed to replace that string
   ipos = index(Channel_1_Filename, "_1_")
   ilen = len(Channel_1_Filename)
    
   first_line_in_segment = (segment_number-1)*num_scans_per_segment

   !---------------------------------------------------------------------------
   ! AHI Navigation (Do Navigation and Solar angles first)
   !---------------------------------------------------------------------------
   
   call AHI_NAVIGATON(1,first_line_in_segment,&
                      num_pix,Num_Scans_Per_Segment,1,&
                      AREAstr,NAVstr_AHI)
   
   if (segment_number == 1) then

      image_jday = jday
      image_time_hours = image_time_ms / 60.0 / 60.0 / 1000.0

      !--- compute scan rate for future use
!      Num_Elements_This_Image =  int(AREAstr%num_elem / AHI_Xstride) + 1
!      Num_Scans_This_Image = AREAstr%num_line
!      Scan_Rate = real((Num_Elements_This_Image)/               &
!                  real(Num_4km_elem_fd/AHI_Xstride)) * &
!                  real((Num_Scans_This_Image) / real(Num_4km_Scans_Fd)) * &
!                  real(Time_For_Fd_Scan) / real(Num_Scans_This_Image)  


      call (trim(Image%LeveL1b_Path)//trim(Channel_1_Filename)//CHAR(0), Severi_File_Id)
   
   endif
  
   !---   read channel 3 (AHI channel 1)
   if (Chan_On_Flag_Default(3) == sym%YES) then

       call READ_AHI_CHANNEL(1, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)

   endif
   

   !---   read channel 12 (AHI channel 2)
   if (Chan_On_Flag_Default(3) == sym%YES) then

       call READ_AHI_CHANNEL(2, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)

   endif
   
   !---   read channel 1 (AHI channel 3)
   if (Chan_On_Flag_Default(1) == sym%YES) then

       call READ_AHI_CHANNEL(3, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)

   endif

   !---   read channel 2 (AHI channel 4)
   if (Chan_On_Flag_Default(1) == sym%YES) then

       call READ_AHI_CHANNEL(4, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)

   endif


   !---   read channel 6 (AHI channel 5)
   if (Chan_On_Flag_Default(6) == sym%YES) then

       call READ_AHI_CHANNEL(5, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)

   endif

   !---   read channel 7 (AHI channel 6)
   if (Chan_On_Flag_Default(1) == sym%YES) then

       call READ_AHI_CHANNEL(6, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)

   endif


   
   !---   read channel 20 (AHI channel 7)
   if (Chan_On_Flag_Default(20) == sym%YES) then

       call READ_AHI_CHANNEL(7, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)

   endif
                    
   
   !---   read channel 27 (AHI channel 9)
   if (Chan_On_Flag_Default(27) == sym%YES) then

       call READ_AHI_CHANNEL(9, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)
   endif

   !---   read channel 28 (AHI channel 10)
   if (Chan_On_Flag_Default(28) == sym%YES) then

       call READ_AHI_CHANNEL(10, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)
   endif


   !---   read channel 29 (AHI channel 11)
   if (Chan_On_Flag_Default(29) == sym%YES) then

       call READ_AHI_CHANNEL(11, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)
   endif


   !---   read channel 30 (AHI channel 12)
   if (Chan_On_Flag_Default(30) == sym%YES) then

       call READ_AHI_CHANNEL(12, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)
   endif

   
   !---   read channel 31 (AHI channel 14)
   if (Chan_On_Flag_Default(31) == sym%YES) then
      
       call READ_AHI_CHANNEL(14, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)
   endif
   
   
   !---   read channel 32 (AHI channel 15)
   if (Chan_On_Flag_Default(32) == sym%YES) then

       call READ_AHI_CHANNEL(15, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)

   endif
    
   !---   read channel 33 (AHI channel 16)
   if (Chan_On_Flag_Default(33) == sym%YES) then

       call READ_AHI_CHANNEL(16, AHI_Xstride, Segment_Number, &
                             Num_Scans_Per_Segment, Num_Scans_Read,   &
                             Two_Byte_Temp)

   endif
    
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + num_scans_read - 1
     Scan_Number(Line_Idx) = first_line_in_segment + Line_Idx
     Scan_Time(Line_Idx) = image_time_ms + (Scan_Number(Line_Idx)-1) * Scan_rate
   enddo

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

        if (Chan_On_Flag_Default(28) == sym%YES) then
           if (ch(28)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(28)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(28,ch(28)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(29) == sym%YES) then
           if (ch(29)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(29)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(29,ch(29)%Rad_Toa(Elem_Idx,Line_Idx))
           endif
        endif

        if (Chan_On_Flag_Default(30) == sym%YES) then
           if (ch(30)%Rad_Toa(Elem_Idx,Line_Idx) > 0.0) then
            ch(30)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(30,ch(30)%Rad_Toa(Elem_Idx,Line_Idx))
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
! AHI Angles
! NOTE: These were private routines in the GOES module. Suggest they become
!       public with different names, since they are used cross platform  
!------------------------------------------------------------------------------
   image_jday = jday
   image_time_hours = image_time_ms / 60.0 / 60.0 / 1000.0
   do iline = Line_Idx_Min_Segment, Line_Idx_Min_Segment + num_scans_read - 1
     do ielem = 1,num_pix
        call POSSOL(image_jday,image_time_hours, &
                    Lon_1b(ielem,iline),Lat_1b(ielem,iline), &
                    Solzen(ielem,iline),Solaz(ielem,iline))
     enddo
     call COMPUTE_SATELLITE_ANGLES(goes_sub_satellite_longitude,  &
                                   goes_sub_satellite_latitude, iline)                      
   enddo
    
 end subroutine READ_AHI
 

 ! Perform AHI Navigation

 subroutine AHI_NAVIGATION(xstart,ystart,xsize,ysize,xstride, &
                           AREAstr,NAVstr_AHI)
    integer(kind=int4) :: xstart, ystart
    integer(kind=int4) :: xsize, ysize
    integer(kind=int4) :: xstride  
    type (AREA_STRUCT) :: AREAstr
    type (GVAR_NAV), intent(in)    :: NAVstr_AHI
    
    integer :: i, j, ii, jj, ierr, imode
    real(kind(0.0d0)) :: latitude, longitude
    real(kind=real4) :: elem, line, height
    real(kind=real4) :: dlon, dlat
    real(kind=real8) :: mjd
    real(kind=real4), dimension(8) :: angles
    integer :: FGF_TYPE = 3 !AHI uses JMA GEOS navigation, so set type here

    NAVstr_AHI_NAV = NAVstr_AHI
    
    imode = -1
    height = 0.0    !Used for parallax correction
    
    lat = Missing_Value_Real4
    lon = Missing_Value_Real4
          
    !HRIT requires actual line and element of being processed.
    !Unlike MSG, AHI requires no switching to different corrdinates.
                
    do j = 1, ysize

        jj = ystart + (j-1)
        
        do i = 1, xsize

            ii = (i - 1) + xstart      ! get element of the image segement
                
            ! again, use common algorithm for CGMS navigation
                 CALL fgf_to_earth(FGF_TYPE,                  &
                                   DBLE(ii),                  &
                                   DBLE(jj),                  &
                                   DBLE(NAVstr_AHI%CFAC),   &
                                   DBLE(NAVstr_AHI%COFF),   &
                                   DBLE(NAVstr_AHI%LFAC),   &
                                   DBLE(NAVstr_AHI%LOFF),   & 
                                   DBLE(NAVstr_AHI%sub_lon), &
                                   longitude,            &
                                   latitude)
                                          
            if (latitude .LE. -999.0) THEN  ! -999.99 is MSV nav missing value
                    Lat_1b(i,j) = Missing_Value_Real4
                    Lon_1b(i,j) = Missing_Value_Real4
                    Space_Mask(i,j) = sym%SPACE
            else
                    Lat_1b(i,j) = real(latitude,kind=real4)
                    Lon_1b(i,j) = real(longitude,kind=real4)
                    
                    ! BecaUSE JMA sets their longitudes from 0 to 360, and
                    ! we want 180 to -180, one last check.
                    
                    if (longitude .GT. 180.0 ) THEN
                        Lon_1b(i,j) = real(longitude,kind=real4) - 360.0
                    endif
                                        
                    Space_Mask(i,j) = sym%NO_SPACE
            endif

        
            enddo
                        
        enddo     
        
    endif
      
 end subroutine AHI_navigation

!======================================================================
subroutine READ_AHI_CHANNEL(Channel, Segment_Number, &
                    Num_Scans_Per_Segment, &
                    Num_Scans_Read, XStride &
                    image)

! character(len=*), intent(in):: filename !May not need depending on APIs
 integer(kind=int4), intent(in):: Channel
 integer(kind=int4), intent(in):: Xstride
 integer(kind=int4), intent(in):: Segment_Number
 integer(kind=int4), intent(in):: Num_Lines_Per_Segment
 integer(kind=int4), intent(out):: Num_Lines_Read
 integer(kind=int2), dimension(:,:), intent(out):: image

!---- hooks into AHI reader here. Basically it'll provide 


end subroutine READ_AHI_CHANNEL

!----------------------------------------------------------------
! read the MSG constants into memory
!-----------------------------------------------------------------
subroutine READ_AHI_INSTR_CONSTANTS(Instr_Const_file)
  character(len=*), intent(in):: Instr_Const_file
  integer:: ios0, erstat
  integer:: Instr_Const_lun
 
  Instr_Const_lun = GET_LUN()

  open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

  print *, "opening ", trim(Instr_Const_file)
  erstat = 0
  if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, "Error opening AHI constants file, ios0 = ", ios0
    stop 19
  endif
  read(unit=Instr_Const_lun,fmt="(a3)") Sat_Name
  read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
  read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
  read(unit=Instr_Const_lun,fmt=*) a1_20, a2_20,nu_20
  read(unit=Instr_Const_lun,fmt=*) a1_27, a2_27,nu_27
  read(unit=Instr_Const_lun,fmt=*) a1_28, a2_28,nu_28
  read(unit=Instr_Const_lun,fmt=*) a1_29, a2_29,nu_29
  read(unit=Instr_Const_lun,fmt=*) a1_30, a2_30,nu_30
  read(unit=Instr_Const_lun,fmt=*) a1_31, a2_31,nu_31
  read(unit=Instr_Const_lun,fmt=*) a1_32, a2_32,nu_32
  read(unit=Instr_Const_lun,fmt=*) a1_33, a2_33,nu_33
  read(unit=Instr_Const_lun,fmt=*) a1_43, a2_43,nu_43     
  read(unit=Instr_Const_lun,fmt=*) a1_44, a2_44,nu_44
  read(unit=Instr_Const_lun,fmt=*) b1_day_mask,b2_day_mask,b3_day_mask,b4_day_mask
  close(unit=Instr_Const_lun)

  !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

  !--- hardwire dark counts   (TBD)
  Ch1_Dark_Count = 29
  Ch2_Dark_Count = 29

end subroutine READ_AHI_INSTR_CONSTANTS

end module AHI_MODULE
