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

module SEVIRI_MODULE
use CONSTANTS
use PIXEL_COMMON
use CALIBRATION_CONSTANTS
use PLANCK
use NUMERICAL_ROUTINES
use GOES_MODULE
use CGMS_NAV
use FILE_UTILITY
use VIEWING_GEOMETRY_MODULE

implicit none


public::  READ_SEVIRI,  &
          CALIBRATE_SEVIRI_DARK_COMPOSITE

         
private:: GET_SEVIRI_NAVIGATION,  &
          LOAD_SEVIRI_CAL_AREA,  &
          MSG_RAD_BT

!---------- info needed for navigation   
real(KIND(0.0d0)), PARAMETER, PRIVATE ::  SUB_LON_MSG = 0.0     ! Longitude of Projection Sub-Satellite Point in degrees
real, PARAMETER, PUBLIC ::  SUB_LON_MET8 = -3.477996     ! Longitude of actual Sub-Satellite Point for Met-8
real, PARAMETER, PUBLIC ::  SUB_LON_MET9 = -0.159799     ! Longitude of actual Sub-Satellite Point for Met-9
integer, PARAMETER, PRIVATE ::  CFAC = -781648343  
integer, PARAMETER, PRIVATE ::  LFAC = -781648343 
integer, PARAMETER, PRIVATE ::  COFF = 1856
integer, PARAMETER, PRIVATE ::  LOFF = 1856
  
!---------- info needed for count -> rad calibration
real, PRIVATE,dimension(11) ::  Slope_Sev
real, PRIVATE,dimension(11) ::  Offset_Sev
 
!------- Since we are going to do Rad -> BT internal, we'll need a, b, nu 
real, PRIVATE,dimension(4)  ::  Solar_Const_Sev
real, PRIVATE,dimension(11) ::  alpha_sev
real, PRIVATE,dimension(11) ::  beta_sev
real, PRIVATE,dimension(11) ::  nu_sev


integer(kind=int4), public, parameter:: Seviri_Xstride = 1
integer(kind=int4), private, parameter:: Num_3km_scans_fd = 3712
integer(kind=int4), private, parameter:: Num_3km_Elem_fd = 3712
integer(kind=int4), private, parameter:: time_for_fd_scan = 900000 !milliseconds (15min)
real, private, save:: Scan_rate    !scan rate in millsec / line
integer(kind=int4), private, parameter:: Seviri_Byte_Shift = 0

contains

!--------------------------------------------------------------------
! assign internal sat id's and const file names for MSG
!--------------------------------------------------------------------
subroutine ASSIGN_MSG_SAT_ID_NUM_INTERNAL(Mcidas_Id_Num)
    integer(kind=int4), intent(in):: Mcidas_Id_Num

    !--- Met-08
    Sensor_Name_Attribute = 'SEVIRI'
    if (Mcidas_Id_Num == 51)   then
        Sc_Id_WMO = 55
        Instr_Const_file = 'met8_instr.dat'
        Algo_Const_file = 'met8_algo.dat'
        Platform_Name_Attribute = 'Meteosat-8'
    endif
    !--- Met-09
    if (Mcidas_Id_Num == 52)   then
        Sc_Id_WMO = 56
        Instr_Const_file = 'met9_instr.dat'
        Algo_Const_file = 'met9_algo.dat'
        Goes_Mop_Flag = sym%NO
        Platform_Name_Attribute = 'Meteosat-9'
    endif
    !--- Met-10
    if (Mcidas_Id_Num == 53)   then
        Sc_Id_WMO = 57
        Instr_Const_file = 'met10_instr.dat'
        Algo_Const_file = 'met10_algo.dat'
        Platform_Name_Attribute = 'Meteosat-10'
    endif
    !--- Met-11
    if (Mcidas_Id_Num == 54)   then
        Sc_Id_WMO = 58
        Instr_Const_file = 'met10_instr.dat'
        Algo_Const_file = 'met10_algo.dat'
        Platform_Name_Attribute = 'Meteosat-11'
    endif

   Instr_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Instr_Const_file)
   Algo_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Algo_Const_file)

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
  read(unit=Instr_Const_lun,fmt=*) a1_20, a2_20,nu_20
  read(unit=Instr_Const_lun,fmt=*) a1_27, a2_27,nu_27
  read(unit=Instr_Const_lun,fmt=*) a1_28, a2_28,nu_28
  read(unit=Instr_Const_lun,fmt=*) a1_29, a2_29,nu_29
  read(unit=Instr_Const_lun,fmt=*) a1_30, a2_30,nu_30
  read(unit=Instr_Const_lun,fmt=*) a1_31, a2_31,nu_31
  read(unit=Instr_Const_lun,fmt=*) a1_32, a2_32,nu_32
  read(unit=Instr_Const_lun,fmt=*) a1_33, a2_33,nu_33
  read(unit=Instr_Const_lun,fmt=*) b1_day_mask,b2_day_mask,b3_day_mask,b4_day_mask
  close(unit=Instr_Const_lun)

  !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
  Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / ew_Ch20

  !--- hardwire dark counts
  Ch1_Dark_Count = 29
  Ch2_Dark_Count = 29

end subroutine READ_MSG_INSTR_CONSTANTS
!-------------------------------------------------------------------------------
! public routine to read data from an AREA file for one segment into memory
!-------------------------------------------------------------------------------
subroutine READ_SEVIRI(Segment_Number,Channel_1_Filename, &
                     jday, image_time_ms, Time_Since_Launch, &
                     AREAstr,NAVstr)

   integer(kind=int4), intent(in):: Segment_Number
   character(len=*), intent(in):: Channel_1_Filename
   TYPE (AREA_STRUCT), intent(in) :: AREAstr
   TYPE (GVAR_NAV), intent(in)    :: NAVstr
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
   integer:: ichan_goes
   integer:: ichan_modis
   integer:: seviri_file_id
   real(kind=real4), save:: image_time_hours
   integer(kind=int4), save:: image_jday
   integer(kind=int4):: first_Line_in_segment
   character(len=2):: Ichan_Goes_String
   integer :: Line_Idx
   integer :: Elem_Idx
   integer :: Num_Elements_this_image
   integer :: Num_scans_this_image
  

   !--- assume Channel_1_file name has a unique "_1_" in the name. 
   !--- determine indices needed to replace that string
   ipos = index(Channel_1_Filename, "_1_")
   ilen = len(Channel_1_Filename)
    
   first_Line_in_segment = (Segment_Number-1)*Num_scans_per_segment
   
   !---------------------------------------------------------------------------
   ! SEVIRI Navigation (Do Navigation first)
   !---------------------------------------------------------------------------
   call GET_SEVIRI_NAVIGATION(1,first_Line_in_segment,&
                              Num_pix,Num_Scans_Per_Segment,1,&
                              AREAstr)

!-------------------------------------------------------------------------------
! uncompress (only on first segment and channel 1 is already done during header read)
!   This is the same basic logic as goes_routines.f90
!-------------------------------------------------------------------------------

   if (Segment_Number == 1) then

      image_jday = jday
      image_time_hours = image_time_ms / 60.0 / 60.0 / 1000.0

      !--- compute scan rate for future use
      Num_Elements_this_image =  int(AREAstr%Num_Elem / Seviri_Xstride) + 1
      Num_scans_this_image = AREAstr%Num_Line
      Scan_Rate = real((Num_Elements_this_image)/               &
               real(Num_3km_Elem_fd/Seviri_Xstride)) * &
               real((Num_scans_this_image) / real(Num_3km_scans_fd)) * &
               real(time_for_fd_scan) / real(Num_scans_this_image)

       do ichan_goes = 2,11

       if (ichan_goes == 2) ichan_modis = 2
       if (ichan_goes == 3) ichan_modis = 6
       if (ichan_goes == 4) ichan_modis = 20
       if (ichan_goes == 5) ichan_modis = 27
       if (ichan_goes == 6) ichan_modis = 28
       if (ichan_goes == 7) ichan_modis = 29
       if (ichan_goes == 8) ichan_modis = 30
       if (ichan_goes == 9) ichan_modis = 31
       if (ichan_goes == 10) ichan_modis = 32
       if (ichan_goes == 11) ichan_modis = 33

       if (ichan_goes < 10) then
          write(Ichan_Goes_String,fmt="(I1.1)") ichan_goes
       else
          write(Ichan_Goes_String,fmt="(I2.2)") ichan_goes
       endif

       if (Chan_On_Flag_Default(ichan_modis) == sym%YES) then

          Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_"//trim(Ichan_Goes_String)//"_" // &
                            Channel_1_Filename(ipos+3:ilen)

          if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
          else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
          endif

          Channel_X_Filename_Full_uncompressed = trim(Dir_1b)//trim(Channel_X_Filename)
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
    
    
    ! On first segment, get slope/offset information from McIDAS Header
    seviri_file_id = get_lun()   

    if (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) THEN
      CALL mread_open(trim(Temporary_Data_Dir)//trim(Channel_1_Filename)//CHAR(0), seviri_file_id)
    ELSE
      CALL mread_open(trim(Dir_1b)//trim(Channel_1_Filename)//CHAR(0), seviri_file_id)
    ENDif  
    CALL LOAD_SEVIRI_CAL_AREA(seviri_file_id, AREAstr)
    CALL mread_close(seviri_file_id)
   endif
     

   !---   read channel 1 (MSG channel 1)
   if (Chan_On_Flag_Default(1) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_1_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_1_Filename)
       endif
      
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

       Ch1_Counts = Two_Byte_Temp


       !----- How to do Ref % utilizing EUMETSAT values. Solar constant is 20.76 for 0.64
       DO Line_Idx=1, Num_Scans_Per_Segment
           DO Elem_Idx=1, Num_pix
           
           ch(1)%Ref_Toa(Elem_Idx,Line_Idx) =  Missing_Value_Real4

           if ( (Space_Mask(Elem_Idx,Line_Idx) == sym%SPACE) .OR. &
                (Solzen(Elem_Idx,Line_Idx) >= 90.0) ) THEN                
                cycle
           endif
        
           
           !Note - Normalization occurrs in NORMALIZE_REFLECTANCES
           ch(1)%Ref_Toa(Elem_Idx,Line_Idx) = 100.0 * ((Slope_Sev(1) * &
                         Two_Byte_Temp(Elem_Idx,Line_Idx)) + &
                         Offset_Sev(1))/ Solar_Const_Sev(1)
           
           
           
           end DO 
       end DO 
   endif
   


   if (Chan_On_Flag_Default(2) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_2_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif
       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

       Ch2_Counts = Two_Byte_Temp
                                    
       DO Line_Idx=1, Num_Scans_Per_Segment
           DO Elem_Idx=1, Num_pix
           
           ch(2)%Ref_Toa(Elem_Idx,Line_Idx) =  Missing_Value_Real4

           if ( (Space_Mask(Elem_Idx,Line_Idx) == sym%SPACE) .OR. &
                (Solzen(Elem_Idx,Line_Idx) >= 90.0) ) THEN                
                cycle
           endif
           
           !Note - Normalization occurrs in NORMALIZE_REFLECTANCES
           ch(2)%Ref_Toa(Elem_Idx,Line_Idx) = 100.0 * ((Slope_Sev(2) * &
                         Two_Byte_Temp(Elem_Idx,Line_Idx)) + &
                         Offset_Sev(2))/ Solar_Const_Sev(2)
           end DO 
       end DO 

   endif

!--- 1.6 micron
   if (Chan_On_Flag_Default(6) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_3_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)
       
!----- How to do Ref % utilizing EUMETSAT values. Solar constant is 19.85 for 1.64

       DO Line_Idx=1, Num_Scans_Per_Segment
           DO Elem_Idx=1, Num_pix
           
           ch(6)%Ref_Toa(Elem_Idx,Line_Idx) =  Missing_Value_Real4

           if ( (Space_Mask(Elem_Idx,Line_Idx) == sym%SPACE) .OR. &
                (Solzen(Elem_Idx,Line_Idx) .GE. 90.0) ) THEN                
                cycle
           endif
        
           
           !Note - Normalization occurrs in NORMALIZE_REFLECTANCES
           ch(6)%Ref_Toa(Elem_Idx,Line_Idx) = 100.0 * ((Slope_Sev(3) * &
                         Two_Byte_Temp(Elem_Idx,Line_Idx)) + &
                         Offset_Sev(3))/ Solar_Const_Sev(3)
           end DO 
       end DO 
   endif

   
! 3.9 micron
   if (Chan_On_Flag_Default(20) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_4_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif
       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      call MSG_RAD_BT(4, 20, Two_Byte_Temp, ch(20)%Bt_Toa, ch(20)%Rad_Toa)


   endif
   
! 6.2 micron
   if (Chan_On_Flag_Default(27) == sym%YES) then

        Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_5_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif
       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

       call MSG_RAD_BT(5, 27, Two_Byte_Temp, ch(27)%Bt_Toa, ch(27)%Rad_Toa)

   endif

! 7.4 micron
   if (Chan_On_Flag_Default(28) == sym%YES) then

        Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_6_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif
       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

       call MSG_RAD_BT(6, 28, Two_Byte_Temp, ch(28)%Bt_Toa, ch(28)%Rad_Toa)

   endif

! 8.5 micron
   if (Chan_On_Flag_Default(29) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_7_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif
       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

       call MSG_RAD_BT(7, 29, Two_Byte_Temp, ch(29)%Bt_Toa, ch(29)%Rad_Toa)

   endif

!9.7 micron
   if (Chan_On_Flag_Default(30) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_8_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif
       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      call MSG_RAD_BT(8, 30, Two_Byte_Temp, ch(30)%Bt_Toa, ch(30)%Rad_Toa)

   endif


!11 micron
   if (Chan_On_Flag_Default(31) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_9_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif
       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)
       
      call MSG_RAD_BT(9, 31, Two_Byte_Temp, ch(31)%Bt_Toa, ch(31)%Rad_Toa)

   endif
   
!12 micron
   if (Chan_On_Flag_Default(32) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_10_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif
       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      call MSG_RAD_BT(10, 32, Two_Byte_Temp, ch(32)%Bt_Toa, ch(32)%Rad_Toa)

   endif

   
!13.3 micron
   if (Chan_On_Flag_Default(33) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_11_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif
       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Seviri_Byte_Shift, &
                                    AREAstr, Seviri_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      call MSG_RAD_BT(11, 33, Two_Byte_Temp, ch(33)%Bt_Toa, ch(33)%Rad_Toa)
      
   endif  

!------------------------------------------------------------------------------
! SEVIRI Angles
! NOTE: These were private routines in the GOES module. Suggest they become
!       public with different names, since they are used cross platform  
!------------------------------------------------------------------------------
   do iline = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1
     do ielem = 1,Num_pix
        call POSSOL(image_jday,image_time_hours, &
                    Lon_1b(ielem,iline),Lat_1b(ielem,iline), &
                    Solzen(ielem,iline),solaz(ielem,iline))
     enddo

      call COMPUTE_SATELLITE_ANGLES(goes_sub_satellite_longitude,  &
                                    goes_sub_satellite_latitude, iline)
   enddo

   !--- compute scantime and scan angles
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1
     Scan_Number(Line_Idx) = first_Line_in_segment + Line_Idx
     Scan_Time(Line_Idx) = image_time_ms + (Scan_number(Line_Idx)-1) * Scan_Rate
   enddo

   !--- ascending node
   ielem = Num_pix/2
   do iline = Line_Idx_Min_Segment+1, Line_Idx_Min_Segment + Num_Scans_Read - 1
     ascend(iline) = 0
     if (lat_1b(ielem,iline) < lat_1b(ielem,iline-1)) then
       ascend(iline) = 1
     endif
   enddo
   ascend(Line_Idx_Min_Segment) = ascend(Line_Idx_Min_Segment+1)

   
end subroutine READ_SEVIRI
 
 
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
        
    DO Line_Idx=1, Num_Scans_Per_Segment
      DO Elem_Idx=1, Num_pix
            
        if (Space_Mask(Elem_Idx,Line_Idx) == sym%SPACE) THEN
            cycle
        endif
        
        
        Rad_Temp = ((Slope_Sev(Chan_Num) * Sev_Counts(Elem_Idx,Line_Idx)) + &
                     Offset_Sev(Chan_Num))
 
        if (Rad_Temp > 0.0 ) THEN
                
            Rad_Out(Elem_Idx,Line_Idx) = Rad_Temp
            
            !--- this uses EUMETSAT's numbers
!           Brit_Temp_Out(Elem_Idx,Line_Idx) =(((c2*nu_sev(Chan_Num)) / &
!                                     (log( 1 + ((c1*(nu_sev(Chan_Num)**3)) /&
!                                      Rad_Temp)))) - beta_sev(Chan_Num)) / &
!                                      alpha_sev(Chan_Num)

            !--- this uses AKH derived numbers
            Brit_Temp_Out(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(Chan_Num_Ref,Rad_Temp)

        endif 
   
      end DO
    end DO
    
     
  end subroutine MSG_RAD_BT 
 !-------------------------------------------------------------------------------
! routine gets Lat/lon for SEVIRI
!-------------------------------------------------------------------------------
subroutine GET_SEVIRI_NAVIGATION(xstart,ystart,xsize,ysize,xstride,AREAstr)

    integer(kind=int4) :: xstart, ystart
    integer(kind=int4) :: xsize, ysize
    integer(kind=int4) :: xstride
    TYPE (AREA_STRUCT), intent(in) ::AREAstr
    
    integer :: i, j, elem, line
    real(kind=real8) :: dlat, dlon
    integer ::  CFAC_MSG = -781648343 
    integer ::  LFAC_MSG = -781648343
    integer ::  COFF_MSG = 1856
    integer ::  LOFF_MSG = 1856
        
       
    do j=1, ysize
      line = (ystart - 1) + j

      ! convert to eumetsat coordinate space if necessary
      if (AREAstr%north_bound == 1) line = AREAstr%Num_Line - line 
      
      do i=1, xsize
         elem = (i - 1)*(xstride) + xstart
         ! convert to eumetsat coordinate space if necessary
         if (AREAstr%west_vis_pixel == 1) elem = AREAstr%Num_Elem - elem + 1 

         CALL pixcoord2geocoord_cgms(elem,         &
                                     line,         &
                                     COFF_MSG,     &
                                     LOFF_MSG,     & 
                                     LFAC_MSG,     &
                                     CFAC_MSG,     &
                                     0,       &
                                     SUB_LON_MSG,  &
                                     dlat,         &
                                     dlon)

                  
         if (dlat == -999.0) then  ! -999.0 is MSG nav missing value
            Lat_1b(i,j) = Missing_Value_Real4
            Lon_1b(i,j) = Missing_Value_Real4
            Space_Mask(i,j) = sym%SPACE
         else
            Lat_1b(i,j) = real(dlat,kind=real4)
            Lon_1b(i,j) = real(dlon,kind=real4)
            Space_Mask(i,j) = sym%NO_SPACE
         endif
         
      end do
   end do
  
end subroutine GET_SEVIRI_NAVIGATION

subroutine LOAD_SEVIRI_CAL_AREA(lun, AREAstr)
  integer(kind=int4), intent(in) :: lun
  type(AREA_STRUCT), intent(in):: AREAstr
  character(len=1252) :: cbuf
  character(len=104) :: cout
  integer :: bandoffset, band, avoid_warning
  real(kind=real8) :: c1w3,c2w,alpha,beta,gain,offset
  integer, parameter :: nbands = 11
  integer, parameter, dimension(nbands) :: arr_index = (/2,3,5,7,8,10,11,12,14,15,16/)
  
  avoid_warning = lun  
  
!  c2 = C_2*1.0e02_real4
  
!Solar constants
  Solar_Const_Sev(1) = 20.76
  Solar_Const_Sev(2) = 23.24
  Solar_Const_Sev(3) = 19.85
  Solar_Const_Sev(4) = 4.92
  
  
! internal plack coeff
  if (Sc_Id_WMO == 55) then !MET-8
        alpha_sev(4) = 0.9959
        alpha_sev(5) = 0.9963
        alpha_sev(6) = 0.9991
        alpha_sev(7) = 0.9996
        alpha_sev(8) = 0.9999
        alpha_sev(9) = 0.9983
        alpha_sev(10) = 0.9988
        alpha_sev(11) = 0.9981

        beta_sev(4) = 3.471
        beta_sev(5) = 2.219
        beta_sev(6) = 0.485
        beta_sev(7) = 0.181
        beta_sev(8) = 0.060
        beta_sev(9) = 0.627
        beta_sev(10) = 0.397
        beta_sev(11) = 0.576
    
  
  endif  
    

!Met-9
  if (Sc_Id_WMO == 56) then !MET-8
        alpha_sev(4) = 0.9954
        alpha_sev(5) = 0.9963
        alpha_sev(6) = 0.9991
        alpha_sev(7) = 0.9996
        alpha_sev(8) = 0.9999
        alpha_sev(9) = 0.9983
        alpha_sev(10) = 0.9988
        alpha_sev(11) = 0.9981
  
        beta_sev(4) = 3.438
        beta_sev(5) = 2.185
        beta_sev(6) = 0.470
        beta_sev(7) = 0.179
        beta_sev(8) = 0.056
        beta_sev(9) = 0.640
        beta_sev(10) = 0.408
        beta_sev(11) = 0.561
  endif
  

  !---------------------------------------------------------------------
  ! Read from the calibration block.  Logic supplied by D. Santek.
  ! mreadf version as implmented in GEOCAT, which is much cleaner than 
  ! original GSIP code
  !---------------------------------------------------------------------

  call mreadf_int_o(lun,AREAstr%cal_offset,1,1252,cbuf)

  ! read nu out of header, along with slope/offset incase of shift
  do band=1, nbands
    bandoffset = (band-1)*104 + 5
    cout(1:104) = cbuf(bandoffset:bandoffset+103)
    read(cout,'(6E17.10)') c1w3,c2w,alpha,beta,gain,offset
    nu_sev(band) = c2w/c2
    !a(arr_index(band),1) = alpha
    !b(arr_index(band),1) = beta
    Slope_Sev(band) = gain
    Offset_Sev(band) = offset

!    print *,'band=',arr_index(band),'c1w2=',c1w3,'nu=',c2w/c2,'a=',alpha,'b=',beta,'gain=',gain,'offset=',offset
!   print *,'band=',arr_index(band),'gain=',slope(arr_index(band),1),'offset=',&
!  offset(arr_index(band),1)
    
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

end module SEVIRI_MODULE
