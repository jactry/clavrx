!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: seviri_module.f90 (src)
!       SEVIRI_MODULE (program)
!
! PURPOSE:
!
! DESCRIPTION: 
!
! AUTHORS:
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE public
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED public USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
!--------------------------------------------------------------------------------------

module SEVIRI_MODULE
   
   use CALIBRATION_CONSTANTS,only: &
    Planck_nu &
    , Planck_A1 &
    , Planck_A2 &
    , sat_name &
    , Solar_Ch20 &
    , EW_Ch20 &
    , Solar_ch20_nu &
    , Ch1_dark_count &
    , Ch2_dark_count
   
   use CGMS_NAV,only:
   
   use CX_CONSTANTS_MOD
   
   use FILE_TOOLS,only: &
      get_lun
   
   use GOES_MODULE,only: &
      gvar_nav &
      , COMPUTE_SATELLITE_ANGLES &
      , GET_IMAGE_FROM_AREAFILE
   
    use NUMERICAL_TOOLS_MOD
   
   use PIXEL_COMMON
   
   use PLANCK
  
   use VIEWING_GEOMETRY_MODULE
   
   use CX_SSEC_AREAFILE_MOD

   implicit none
   private

   public::  READ_SEVIRI,  &
          CALIBRATE_SEVIRI_DARK_COMPOSITE, &
          READ_NAVIGATION_BLOCK_SEVIRI
   public :: READ_MSG_INSTR_CONSTANTS        
private:: GET_SEVIRI_NAVIGATION,  &
          LOAD_SEVIRI_CAL_AREA,  &
          MSG_RAD_BT

! SUB_LON_MSG needs to be read from the nav block as it can be different
! from 0.0 when looking at Met-8 Indian Ocean data.
!---------- info needed for navigation   
!---stwreal(KIND(0.0d0)), parameter, private ::  SUB_LON_MSG = 0.0     ! Longitude of Projection Sub-Satellite Point in degrees
!real, parameter, public ::  SUB_LON_MET8 = -3.477996     ! Longitude of actual Sub-Satellite Point for Met-8
!real, parameter, public ::  SUB_LON_MET9 = -0.159799     ! Longitude of actual Sub-Satellite Point for Met-9
integer, parameter, private ::  CFAC = -781648343  
integer, parameter, private ::  LFAC = -781648343 
integer, parameter, private ::  COFF = 1856
integer, parameter, private ::  LOFF = 1856
  
!---------- info needed for count -> rad calibration
real, private, dimension(11) ::  Slope_Sev
real, private, dimension(11) ::  Offset_Sev
 
!------- Since we are going to do Rad -> BT internal, we'll need a, b, nu 
real, private, dimension(4)  ::  Solar_Const_Sev


integer(kind=int4), public, parameter:: Seviri_Xstride = 1
integer(kind=int4), private, parameter:: Num_3km_scans_fd = 3712
integer(kind=int4), private, parameter:: Num_3km_Elem_fd = 3712
integer(kind=int4), private, parameter:: time_for_fd_scan = 900000 !milliseconds (15min)
real, private, save:: Scan_rate    !scan rate in millsec / line
integer(kind=int4), private, parameter:: Seviri_Byte_Shift = 0
integer, dimension(11), parameter, private:: Chan_Idx=[1,2,6,20,27,28,29,30,31,32,33]

contains

!--------------------------------------------------------------------
! assign internal sat id's and const file names for MSG
!--------------------------------------------------------------------
subroutine ASSIGN_MSG_SAT_ID_NUM_INTERNAL(Mcidas_Id_Num)
    integer(kind=int4), intent(in):: Mcidas_Id_Num

    !--- Met-08
    Sensor%Sensor_Name = 'SEVIRI'
    if (Mcidas_Id_Num == 51)   then
        Sensor%WMO_Id = 55
        Sensor%Instr_Const_File = 'met8_instr.dat'
        Sensor%Algo_Const_File = 'met8_algo.dat'
        Sensor%Platform_Name = 'Meteosat-8'
    endif
    !--- Met-09
    if (Mcidas_Id_Num == 52)   then
        Sensor%WMO_Id = 56
        Sensor%Instr_Const_File = 'met9_instr.dat'
        Sensor%Algo_Const_File = 'met9_algo.dat'
        Sensor%Platform_Name = 'Meteosat-9'
    endif
    !--- Met-10
    if (Mcidas_Id_Num == 53)   then
        Sensor%WMO_Id = 57
        Sensor%Instr_Const_File = 'met10_instr.dat'
        Sensor%Algo_Const_File = 'met10_algo.dat'
        Sensor%Platform_Name = 'Meteosat-10'
    endif
    !--- Met-11
    if (Mcidas_Id_Num == 54)   then
        Sensor%WMO_Id = 58
        Sensor%Instr_Const_File = 'met11_instr.dat'
        Sensor%Algo_Const_File = 'met11_algo.dat'
        Sensor%Platform_Name = 'Meteosat-11'
    endif

end subroutine ASSIGN_MSG_SAT_ID_NUM_INTERNAL
!----------------------------------------------------------------
! read the MSG constants into memory
!-----------------------------------------------------------------
subroutine READ_MSG_INSTR_CONSTANTS(Instr_Const_file)
 character(len=*), intent(in):: Instr_Const_file
 integer:: ios0, erstat
 integer:: Instr_Const_lun

 Instr_Const_lun = GET_LUN()

 open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)

 print *, "opening ", trim(Instr_Const_file)
 erstat = 0
 if (ios0 /= 0) then
    erstat = 19
    print *, EXE_PROMPT, "Error opening SEVIRI constants file, ios0 = ", ios0
    stop 19
 endif
  read(unit=Instr_Const_lun,fmt="(a3)") sat_name
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
  close(unit=Instr_Const_lun)

  !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

  !--- hardwire dark counts
  Ch1_Dark_Count = 29
  Ch2_Dark_Count = 29

end subroutine READ_MSG_INSTR_CONSTANTS
!-------------------------------------------------------------------------------
! public routine to read data from an AREA file for one segment into memory
!-------------------------------------------------------------------------------
subroutine READ_SEVIRI(Segment_Number,Channel_1_Filename, &
                       Day_Of_Year, Image_Time_Ms, &
                       AREAstr)

   integer(kind=int4), intent(in):: Segment_Number
   character(len=*), intent(in):: Channel_1_Filename
   TYPE (area_header_type), intent(in) :: AREAstr
   integer(kind=int2), intent(in):: Day_Of_Year
   integer(kind=int4), intent(in):: Image_Time_Ms

   character(len=1020):: Channel_X_Filename
   character(len=1020):: Channel_X_Filename_Full
   character(len=1020):: Channel_X_Filename_Full_uncompressed
   character(len=1020):: System_String

   integer:: ipos
   integer:: ilen
   integer:: Elem_Idx
   integer:: Line_Idx
   integer:: Chan_Idx_MSG
   integer:: Severi_File_Id
   real(kind=real4), save:: image_time_hours
   integer(kind=int4), save:: Image_Day_Of_Year
   integer(kind=int4):: First_Line_In_Segment
   character(len=2):: Chan_Idx_MSG_String
   integer :: Num_Elements_this_image
   integer :: Num_Scans_This_Image
  

   !--- assume Channel_1_file name has a unique "_1_" in the name. 
   !--- determine indices needed to replace that string
   ipos = index(Channel_1_Filename, "_1_")
   ilen = len(Channel_1_Filename)
    
   First_Line_In_Segment = (Segment_Number-1)*Image%Number_Of_Lines_Per_Segment
   
   !---------------------------------------------------------------------------
   ! SEVIRI Navigation (Do Navigation first)
   !---------------------------------------------------------------------------
   call GET_SEVIRI_NAVIGATION(1,First_Line_In_Segment,&
                              Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment,1,&
                              AREAstr)

   !--------------------------------------------------------------------------------------
   ! uncompress (only on first segment and channel 1 is already done during header read)
   !--------------------------------------------------------------------------------------
   if (Segment_Number == 1) then  !first segment check

      Image_Day_Of_Year = Day_Of_Year
      image_time_hours = Image_Time_Ms / 60.0 / 60.0 / 1000.0

      !--- compute scan rate for future use
      Num_Elements_this_image =  int(AREAstr%Num_Elem / Seviri_Xstride) + 1
      Num_scans_this_image = AREAstr%Num_Line
      Scan_Rate = real((Num_Elements_this_image)/               &
               real(Num_3km_Elem_fd/Seviri_Xstride)) * &
               real((Num_scans_this_image) / real(Num_3km_scans_fd)) * &
               real(time_for_fd_scan) / real(Num_scans_this_image)

       !-------------------------------------------------------------------------------
       !decompress
       !-------------------------------------------------------------------------------
       do Chan_Idx_MSG = 2,11

       if (Chan_Idx_MSG < 10) then
          write(Chan_Idx_MSG_String,fmt="(I1.1)") Chan_Idx_MSG
       else
          write(Chan_Idx_MSG_String,fmt="(I2.2)") Chan_Idx_MSG
       endif

       if (Sensor%Chan_On_Flag_Default(Chan_Idx(Chan_Idx_Msg)) ) then

          Channel_X_Filename = Channel_1_Filename(1:ipos-1) //  &
                               "_"//trim(Chan_Idx_MSG_String)//"_" // &
                               Channel_1_Filename(ipos+3:ilen)

          if (L1b_gzip  .or. L1b_bzip2) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
          else
               Channel_X_Filename_Full = trim(Image%LeveL1b_Path)//trim(Channel_X_Filename)
          endif

          Channel_X_Filename_Full_uncompressed = trim(Image%LeveL1b_Path)//trim(Channel_X_Filename)
          if (L1b_gzip ) then
              System_String = "gunzip -c "//trim(Channel_X_Filename_Full_uncompressed)//".gz"// &
                                " > "//trim(Channel_X_Filename_Full)
              call system(System_String)

              Number_of_Temporary_Files = Number_of_Temporary_Files + 1
              Temporary_File_Name(Number_of_Temporary_Files) = trim(Channel_X_Filename)

          endif
          if (L1b_bzip2 ) then
              System_String = "bunzip2 -c "//trim(Channel_X_Filename_Full_uncompressed)//".bz2"// &
                                " > "//trim(Channel_X_Filename_Full)
              call system(System_String)

              Number_of_Temporary_Files = Number_of_Temporary_Files + 1
              Temporary_File_Name(Number_of_Temporary_Files) = trim(Channel_X_Filename)
          endif

      endif

    enddo
    
    ! On first segment, get slope/offset information from McIDAS Header
    Severi_File_Id = get_lun()   

    if (L1b_gzip  .OR. L1b_bzip2 ) then
      call MREAD_OPEN(trim(Temporary_Data_Dir)//trim(Channel_1_Filename)//CHAR(0), Severi_File_Id)
    else
      call MREAD_OPEN(trim(Image%LeveL1b_Path)//trim(Channel_1_Filename)//CHAR(0), Severi_File_Id)
    endif
    call LOAD_SEVIRI_CAL_AREA(Severi_File_Id, AREAstr)
    call MREAD_CLOSE(Severi_File_Id)

   endif  !end of check of first segment
     
   !---------------------------------------------------------------------------------------------------
   !  Read and Calibrate
   !---------------------------------------------------------------------------------------------------
   read_cal_loop: do Chan_Idx_MSG = 1,11

       if (Chan_Idx_MSG < 10) then
          write(Chan_Idx_MSG_String,fmt="(I1.1)") Chan_Idx_MSG
       else
          write(Chan_Idx_MSG_String,fmt="(I2.2)") Chan_Idx_MSG
       endif

       if (Sensor%Chan_On_Flag_Default(Chan_Idx(Chan_Idx_MSG))) then

         Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_"//trim(Chan_Idx_MSG_String)//"_" // &
                            Channel_1_Filename(ipos+3:ilen)

         if (L1b_gzip  .or. L1b_bzip2 ) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
         else
               Channel_X_Filename_Full = trim(Image%LeveL1b_Path)//trim(Channel_X_Filename)
         endif

         call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Image%Number_Of_Lines_Per_Segment, &
                                    Image%Number_Of_Lines_Read_This_Segment,   &
                                    Two_Byte_Temp, &
                                    Goes_Scan_Line_Flag)

         !--- solar reflectance cal
         if (Chan_Idx_MSG <= 3) then
           ch(Chan_Idx(Chan_Idx_MSG))%Ref_Toa =  Missing_Value_Real4
           where(Space_Mask == sym%NO .and. Geo%Solzen < 90.0)
              ch(Chan_Idx(Chan_Idx_MSG))%Ref_Toa = 100.0 * ((Slope_Sev(Chan_Idx_MSG) * Two_Byte_Temp) + &
                     Offset_Sev(Chan_Idx_MSG))/ Solar_Const_Sev(Chan_Idx_MSG)
           endwhere
         endif

         !--- thermal cal
         if (Chan_Idx_MSG >= 4) then
           call MSG_RAD_BT(Chan_Idx_MSG, Chan_Idx(Chan_Idx_MSG), Two_Byte_Temp,  &
                           ch(Chan_Idx(Chan_Idx_MSG))%Bt_Toa, ch(Chan_Idx(Chan_Idx_MSG))%Rad_Toa)
         endif

         !--- save these to CLAVR-x global variables
         if (Chan_Idx(Chan_Idx_MSG) == 1) Ch1_Counts = Two_Byte_Temp
         if (Chan_Idx(Chan_Idx_MSG) == 2) Ch2_Counts = Two_Byte_Temp
         if (Chan_Idx(Chan_Idx_MSG) == 3) Ch6_Counts = Two_Byte_Temp

         !--- apply KNMI suggested solar calibration tweaks
         !ch(1)%Ref_Toa(Elem_Idx,Line_Idx) = ch(1)%Ref_Toa(Elem_Idx,Line_Idx) / 0.92 
         !ch(6)%Ref_Toa(Elem_Idx,Line_Idx) = ch(6)%Ref_Toa(Elem_Idx,Line_Idx) / 1.035 

      endif

   enddo read_cal_loop

   !------------------------------------------------------------------------------
   ! Compute SEVIRI Angles
   !------------------------------------------------------------------------------
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     do Elem_Idx = 1,Image%Number_Of_Elements
        call POSSOL(Image_Day_Of_Year,image_time_hours, &
                    Nav%Lon_1b(Elem_Idx,Line_Idx),Nav%Lat_1b(Elem_Idx,Line_Idx), &
                    Geo%Solzen(Elem_Idx,Line_Idx),Geo%solaz(Elem_Idx,Line_Idx))
     enddo

      call COMPUTE_SATELLITE_ANGLES(Sensor%Geo_Sub_Satellite_Longitude,  &
                                    Sensor%Geo_Sub_Satellite_Latitude, Line_Idx)
   enddo

   !--- compute scantime and scan angles
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
     Scan_Number(Line_Idx) = First_Line_In_Segment + Line_Idx
     Scan_Time(Line_Idx) = Image_Time_Ms + (Scan_number(Line_Idx)-1) * Scan_Rate
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
   
end subroutine READ_SEVIRI
 
!=========================================================================================== 
!
!=========================================================================================== 
subroutine MSG_RAD_BT(Chan_Num, Chan_Num_Ref,Sev_Counts, Brit_Temp_Out, Rad_Out)
!This subroutine takes a radiance and converts it to a temperature
!ONLY TO BE USED WITH SEVIRI

    integer (kind=real4), intent(in):: Chan_Num
    integer (kind=real4), intent(in):: Chan_Num_Ref
    integer (kind=INT2), dimension(:,:), intent(in):: Sev_Counts
    real (kind=real4), dimension(:,:), intent(out):: Brit_Temp_Out, Rad_Out
    real (KIND=real4)::  Rad_Temp
    integer :: Line_Idx, Elem_Idx

    Rad_Out = Missing_Value_real4
    Brit_Temp_Out = Missing_Value_real4
    Rad_Temp = Missing_Value_real4
        
    do Line_Idx=1, Image%Number_Of_Lines_Per_Segment
      do Elem_Idx=1, Image%Number_Of_Elements
            
        if (Space_Mask(Elem_Idx,Line_Idx) == sym%SPACE) then
            cycle
        endif
        
        
        Rad_Temp = ((Slope_Sev(Chan_Num) * Sev_Counts(Elem_Idx,Line_Idx)) + &
                     Offset_Sev(Chan_Num))
 
        if (Rad_Temp > 0.0 ) then
                
            Rad_Out(Elem_Idx,Line_Idx) = Rad_Temp
            
            Brit_Temp_Out(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(Chan_Num_Ref,Rad_Temp)

        endif 
   
      enddo
    enddo
    
     
  end subroutine MSG_RAD_BT 
!-------------------------------------------------------------------------------
! routine gets Lat/lon for SEVIRI
!-------------------------------------------------------------------------------
subroutine GET_SEVIRI_NAVIGATION(xstart,ystart,xsize,ysize,xstride,AREAstr)

    integer(kind=int4) :: xstart, ystart
    integer(kind=int4) :: xsize, ysize
    integer(kind=int4) :: xstride
    TYPE (area_header_type), intent(in) ::AREAstr
    type (GVAR_NAV) :: NAVstr
    
    integer :: i, j, elem, line
    real(kind=real8) :: dlat, dlon
    integer ::  CFAC_MSG = -781648343 
    integer ::  LFAC_MSG = -781648343
    integer ::  COFF_MSG = 1856
    integer ::  LOFF_MSG = 1856
    integer :: FGF_TYPE = 2 !MSG uses EUMETSAT GEOS navigation, so set type here
    real(KIND(0.0d0)) ::  SUB_LON_MSG      ! Longitude of Projection Sub-Satellite Point in degrees

    ! Read the satellite sub point longitude from the AREA file.
    call READ_NAVIGATION_BLOCK_SEVIRI(trim(Image%Level1b_Full_Name), AREAstr, NAVstr)
    SUB_LON_MSG = NAVstr%sublon

    do j=1, ysize
      line = (ystart - 1) + j

      ! convert to eumetsat coordinate space if necessary
      if (AREAstr%north_bound == 1) line = AREAstr%Num_Line - line 
      
      do i=1, xsize
         elem = (i - 1)*(xstride) + xstart
         ! convert to eumetsat coordinate space if necessary
         if (AREAstr%west_vis_pixel == 1) elem = AREAstr%Num_Elem - elem + 1 

         !call pixcoord2geocoord_cgms(elem,         &
         !                            line,         &
         !                            COFF_MSG,     &
         !                            LOFF_MSG,     & 
         !                            LFAC_MSG,     &
         !                            CFAC_MSG,     &
         !                            0,       &
         !                            SUB_LON_MSG,  &
         !                            dlat,         &
         !                            dlon)

          ! Putting in hooks for updated navigation code
          ! to be edited after 1/12

          CALL fgf_to_earth(FGF_TYPE,        &
                            DBLE(elem),      &
                            DBLE(line),      &
                            DBLE(CFAC_MSG),  &
                            DBLE(COFF_MSG),  &
                            DBLE(LFAC_MSG),  &
                            DBLE(LOFF_MSG),  &
                            SUB_LON_MSG,     &
                            dlon,            &
                            dlat)

         if (dlat == -999.0) then  ! -999.0 is MSG nav missing value
            Nav%Lat_1b(i,j) = Missing_Value_Real4
            Nav%Lon_1b(i,j) = Missing_Value_Real4
            Space_Mask(i,j) = sym%SPACE
         else
            Nav%Lat_1b(i,j) = real(dlat,kind=real4)
            Nav%Lon_1b(i,j) = real(dlon,kind=real4)
            Space_Mask(i,j) = sym%NO_SPACE
         endif
         
      end do
   end do
  
end subroutine GET_SEVIRI_NAVIGATION

!--------------------------------------------------------------------------------
! AKH This needs to go away!
!--------------------------------------------------------------------------------
subroutine LOAD_SEVIRI_CAL_AREA(lun, AREAstr)
  integer(kind=int4), intent(in) :: lun
  type(area_header_type), intent(in):: AREAstr
  character(len=1252) :: cbuf
  character(len=104) :: cout
  integer :: bandoffset, band, avoid_warning
  real(kind=real8) :: c1w3,c2w,alpha,beta,gain,offset
  integer, parameter :: nbands = 11
  
  avoid_warning = lun  
  
  !Solar constants
  Solar_Const_Sev(1) = 20.76
  Solar_Const_Sev(2) = 23.24
  Solar_Const_Sev(3) = 19.85
  Solar_Const_Sev(4) = 4.92
  
  
  !---------------------------------------------------------------------
  ! Read from the calibration block.  Logic supplied by D. Santek.
  !---------------------------------------------------------------------

  call MREADF_INT_O(lun,AREAstr%Cal_Offset,1,1252,cbuf)

  ! read nu out of header, along with slope/offset incase of shift
  do band=1, nbands
    bandoffset = (band-1)*104 + 5
    cout(1:104) = cbuf(bandoffset:bandoffset+103)
    read(cout,'(6E17.10)') c1w3,c2w,alpha,beta,gain,offset
    Slope_Sev(band) = gain
    Offset_Sev(band) = offset
  enddo
  
end subroutine LOAD_SEVIRI_CAL_AREA
 
!============================================================================================
!  Calibrate the Ch1 Dark-Sky Composite Counts
!============================================================================================
subroutine CALIBRATE_SEVIRI_DARK_COMPOSITE(Ch1_Counts_Composite,Ref_Ch1_Dark)

integer(kind=int2), dimension(:,:), intent(in):: Ch1_Counts_Composite
real(kind=real4), dimension(:,:),  intent(out):: Ref_Ch1_Dark

        !---  initialzie
        Ref_Ch1_Dark = Missing_Value_Real4

        !--- calibrate the counts (this is an un-normalized reflectance)
        where(Ch1_Counts_Composite > 0)
          Ref_Ch1_Dark = 100.0 * ((Slope_Sev(1) * Ch1_Counts_Composite) + Offset_Sev(1))/Solar_Const_Sev(1)    
        end where

end subroutine CALIBRATE_SEVIRI_DARK_COMPOSITE

   !
   !
   !
   subroutine READ_NAVIGATION_BLOCK_SEVIRI(filename, AREAstr, NAVstr)

      CHARACTER(len=*), intent(in):: filename
      type(area_header_type), intent(in):: AREAstr
      type(GVAR_NAV), intent(inout):: NAVstr

      integer(kind=int4) :: nav_offset
      integer:: number_of_words_read
      integer(kind=int4), dimension(640) :: i4buf

      ! Navigation block offset is read from the McIDAS AREA file.
      nav_offset = AREAstr%sec_key_nav

      ! Read the navigation block.
      call mreadf_int(trim(filename)//CHAR(0),nav_offset,4,640,&
                    number_of_words_read, i4buf)

      ! Isolate the navigation block.
      call move_bytes(4,i4buf(1),NAVstr%nav_type,0)

      ! Extract the satellite sub longitude point, and convert to positive east.
      NAVstr%sub_lon = real(i4buf(6),kind=real4) / 10000 * real(-1.0)
      NAVstr%sublon = NAVstr%sub_lon

   end subroutine READ_NAVIGATION_BLOCK_SEVIRI

end module SEVIRI_MODULE
