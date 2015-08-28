!// $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: mtsat_module.f90 (src)
!       MTSAT_MODULE (program)
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
MODULE MTSAT_MODULE

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
 public :: READ_MTSAT_DARK_COMPOSITE
 public :: CALIBRATE_MTSAT_DARK_COMPOSITE
 public :: READ_MTSAT_INSTR_CONSTANTS
 public :: ASSIGN_MTSAT_SAT_ID_NUM_INTERNAL
         
 private :: MTSAT_RADIANCE_BT,MTSAT_REFLECTANCE, &
            mtsat_navigation, MGIVSR 
 

 TYPE (GVAR_NAV), PRIVATE    :: NAVstr_MTSAT_NAV
 integer, PARAMETER, PRIVATE :: nchan_mtsat= 5
 INTEGER, PARAMETER, PRIVATE :: ndet_mtsat = 4
 INTEGER, PARAMETER, PRIVATE :: ntable_mtsat = 1024

 INTEGER, PRIVATE :: nref_table_mtsat
 INTEGER, PRIVATE :: nbt_table_mtsat
 CHARACTER(len=4), SAVE, PRIVATE:: calib_type

 INTEGER (kind=int4), dimension(nchan_mtsat,ndet_mtsat,ntable_mtsat), PRIVATE  :: ref_table
 INTEGER (kind=int4), dimension(nchan_mtsat,ndet_mtsat,ntable_mtsat), PRIVATE  :: bt_table
 INTEGER (kind=int4), dimension(nchan_mtsat,ndet_mtsat,ntable_mtsat), PRIVATE  :: rad_table

 integer(kind=int4), public, parameter:: Mtsat_Xstride = 1
 integer(kind=int4), private, parameter:: Num_4km_Scans_Fd = 3712
 integer(kind=int4), private, parameter:: Num_4km_elem_fd = 3712
 integer(kind=int4), private, parameter:: Time_For_Fd_Scan =  1560000 !milliseconds (26min)
 real, private, save:: Scan_rate    !scan rate in millsec / line
 integer(kind=int4), private, parameter:: Mtsat_Byte_Shift = 0


 CONTAINS
!--------------------------------------------------------------------
! assign internal sat id's and const file names for FY2
!--------------------------------------------------------------------
subroutine ASSIGN_MTSAT_SAT_ID_NUM_INTERNAL(Mcidas_Id_Num)
    integer(kind=int4), intent(in):: Mcidas_Id_Num

    if (Mcidas_Id_Num == 84)   then
        Sc_Id_Internal = 34
        Sc_Id_WMO = 171
        Instr_Const_file = 'mtsat1r_instr.dat'
        Algo_Const_file = 'mtsat1r_algo.dat'
        Platform_Name_Attribute = 'MTSAT-1R'
        Sensor_Name_Attribute = 'IMAGER'
!       Sensor_Name_Attribute = 'MTSAT-1R : Imager'
    endif
    if (Mcidas_Id_Num == 85)   then
        Sc_Id_Internal = 35
        Sc_Id_WMO = 172
        Instr_Const_file = 'mtsat2_instr.dat'
        Algo_Const_file = 'mtsat2_algo.dat'
        Goes_Mop_Flag = sym%NO
        Platform_Name_Attribute = 'MTSAT-2'
        Sensor_Name_Attribute = 'IMAGER'
!       Sensor_Name_Attribute = 'MTSAT-2 : Imager'
    endif

   Instr_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Instr_Const_file)
   Algo_Const_file = trim(ancil_data_dir)//"avhrr_data/"//trim(Algo_Const_file)

end subroutine ASSIGN_MTSAT_SAT_ID_NUM_INTERNAL
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

end subroutine READ_MTSAT_INSTR_CONSTANTS

 ! Perform MTSAT Reflectance and BT calibration
 SUBROUTINE READ_MTSAT(segment_number,Channel_1_Filename, &
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
    
   first_line_in_segment = (segment_number-1)*num_scans_per_segment

   !---------------------------------------------------------------------------
   ! MTSAT Navigation (Do Navigation and Solar angles first)
   !---------------------------------------------------------------------------
   
   call mtsat_navigation(1,first_line_in_segment,&
                              num_pix,Num_Scans_Per_Segment,1,&
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

       if (Chan_On_Flag_Default(Chan_Idx_Modis) == sym%YES) then

          Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_"//trim(Chan_Idx_Mtsat_String)//"_" // &
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

    ! On first segment, reflectance, BT and rad tables
    ! On first segment, get slope/offset information from McIDAS Header
    mtsat_file_id = get_lun()   
    IF (l1b_gzip == sym%YES .OR. l1b_bzip2 == sym%YES) THEN
        CALL mread_open(trim(Temporary_Data_Dir)//trim(Channel_1_Filename)//CHAR(0), mtsat_file_id)
    ELSE
      CALL mread_open(trim(Dir_1b)//trim(Channel_1_Filename)//CHAR(0), mtsat_file_id)
    ENDIF  

    CALL load_mtsat_calibration(mtsat_file_id, AREAstr)
    CALL mread_close(mtsat_file_id)


   endif
   
   
   !---   read channel 1 (MTSAT channel 1)
   if (Chan_On_Flag_Default(1) == sym%YES) then

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_1_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_1_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Mtsat_Byte_Shift, &
                                    AREAstr, Mtsat_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

        call MTSAT_REFLECTANCE(Two_Byte_Temp,ch(1)%Ref_Toa(:,:))

        Ch1_Counts = Two_Byte_Temp

   endif
   
   !---   read channel 20 (MTSAT channel 5)
   if (Chan_On_Flag_Default(20) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_5_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Mtsat_Byte_Shift, &
                                    AREAstr, Mtsat_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      call MTSAT_RADIANCE_BT(5_int1, Two_Byte_Temp, ch(20)%Rad_Toa, ch(20)%Bt_Toa)

   endif
                    
   
   !---   read channel 27 (MTSAT channel 4)
   if (Chan_On_Flag_Default(27) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_4_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Mtsat_Byte_Shift, &
                                    AREAstr, Mtsat_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      call MTSAT_RADIANCE_BT(4_int1, Two_Byte_Temp, ch(27)%Rad_Toa, ch(27)%Bt_Toa)

   endif


   
   !---   read channel 31 (MTSAT channel 2)
   if (Chan_On_Flag_Default(31) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_2_" // &
                            Channel_1_Filename(ipos+3:ilen)
       
       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif
       
       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Mtsat_Byte_Shift, &
                                    AREAstr, Mtsat_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      call MTSAT_RADIANCE_BT(2_int1, Two_Byte_Temp, ch(31)%Rad_Toa, ch(31)%Bt_Toa)

   endif
   
   
   !---   read channel 32 (MTSAT channel 3)
   if (Chan_On_Flag_Default(32) == sym%YES) then

       Channel_X_Filename = Channel_1_Filename(1:ipos-1) // "_3_" // &
                            Channel_1_Filename(ipos+3:ilen)

       if (l1b_gzip == sym%YES .or. l1b_bzip2 == sym%YES) then
               Channel_X_Filename_Full = trim(Temporary_Data_Dir)//trim(Channel_X_Filename)
       else
               Channel_X_Filename_Full = trim(Dir_1b)//trim(Channel_X_Filename)
       endif

       call GET_IMAGE_FROM_AREAFILE(trim(Channel_X_Filename_Full), &
                                    Mtsat_Byte_Shift, &
                                    AREAstr, Mtsat_Xstride, &
                                    Segment_Number, &
                                    Num_Scans_Per_Segment, &
                                    Num_Scans_Read,   &
                                    Two_Byte_Temp)

      call MTSAT_RADIANCE_BT(3_int1, Two_Byte_Temp, ch(32)%Rad_Toa, ch(32)%Bt_Toa)

   endif
    
   do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + num_scans_read - 1
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
   do iline = Line_Idx_Min_Segment, Line_Idx_Min_Segment + num_scans_read - 1
     do ielem = 1,num_pix
        call POSSOL(image_jday,image_time_hours, &
                    Lon_1b(ielem,iline),Lat_1b(ielem,iline), &
                    Solzen(ielem,iline),Solaz(ielem,iline))
     enddo
     call COMPUTE_SATELLITE_ANGLES(goes_sub_satellite_longitude,  &
                                   goes_sub_satellite_latitude, iline)                      
   enddo

   !--- ascending node
   ielem = num_pix/2
   do iline = Line_Idx_Min_Segment+1, Line_Idx_Min_Segment + num_scans_read - 1
     ascend(iline) = 0
     if (lat_1b(ielem,iline) < lat_1b(ielem,iline-1)) then
       ascend(iline) = 1
     endif
   enddo
   ascend(Line_Idx_Min_Segment) = ascend(Line_Idx_Min_Segment+1)

    
 END SUBROUTINE READ_MTSAT
 
 
 
 SUBROUTINE load_mtsat_calibration(lun, AREAstr)
  INTEGER(kind=int4), intent(in) :: lun
  type(AREA_STRUCT), intent(in):: AREAstr
  INTEGER(kind=int4), dimension(6528) :: ibuf
  CHARACTER(len=25) :: cbuf
  INTEGER :: nref, nbt, i, j, offset
  INTEGER(kind=int4) :: band_offset_2, band_offset_14, band_offset_15, &
                        band_offset_9, band_offset_7, dir_offset
  REAL(kind=real4) :: albedo, temperature, radiance
  REAL(kind=real4), dimension(5)  :: a_mtsat, b_mtsat, nu_mtsat

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
  
  IF (AREAstr%sat_id_num  == 84) THEN
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
  
  
  ENDIF
  


  IF (AREAstr%sat_id_num  == 85) THEN
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

  ENDIF
  
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
 
 END SUBROUTINE load_mtsat_calibration


 ! Perform MTSAT Reflectance calculation

 SUBROUTINE MTSAT_Reflectance(Mtsat_Counts, alb_temp)
                              
    INTEGER (kind=INT2), dimension(:,:), intent(in):: Mtsat_Counts
    REAL (KIND=real4), dimension(:,:), intent(out):: alb_temp

    INTEGER :: index
    INTEGER:: i, j
    
    
    !---- WCS3 ------!
    ! Reflectance linear fit of reflectance table for HRIT files is taken from JMA Excel table
    ! of calibration coefficents.
    ! http://mscweb.kishou.go.jp/operation/calibration/mt1r/HRIT/mtsat1r_hrit_calibration_table.xls
    !Currently the same for MTSAT-1r and MTSAT-2. McIDAS currently doesn't calculate this.
    !---- WCS3 ------!

    DO j=1, num_scans_read
      DO i=1, num_pix
        IF (Space_Mask(i,j) == sym%NO_SPACE) THEN
          !  Because MTSAT has two vis calibration type, we need to 
          !  which route it needs to go. calibration.f90 was edited to 
          !  gather the vis calibration type, which is stored in thesat_info
          !  structure. Now we USE it here, based off of how to do it via
          !  kbxmtst.dlm (v1.5) from McIDAS - WCS3 - 3/29/2010
                    
          IF ( calib_type == 'MVSH' .OR. &
               calib_type == 'HSVM') THEN

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

            alb_temp(i,j) = 0.000978494701531316*(index) - 0.00197851402698724
                        
            !Recall that reflectance factor goes from 0 to 1 and is not corrected by
            ! cossolzen and SED. Correction with SED and Cossolzen done elsewhere - WCS3
            
            alb_temp(i,j) = (alb_temp(i,j) * 100.0) 
                              
          ENDIF
          
          IF (calib_type == 'MVIS' .OR. &
               calib_type == 'SIVM') THEN
          
            !Not sure if I need to add 1 here
            index = int(Mtsat_Counts(i,j),KIND=int2)

            IF((index .GT. 0) .AND. (index .LE. 1024)) THEN
            
                alb_temp(i,j) = &
             (REAL(ref_table(2,1,index),KIND=real4) /  100.) 
     
            ENDIF                        
            
          ENDIF
        ELSE

            alb_temp(i,j) = Missing_Value_Real4

        ENDIF      
      END DO
    END DO
        
 END SUBROUTINE MTSAT_Reflectance

 ! Perform MTSAT Navigation

 SUBROUTINE mtsat_navigation(xstart,ystart,xsize,ysize,xstride, &
                            AREAstr,NAVstr_MTSAT)
    INTEGER(KIND=int4) :: xstart, ystart
    INTEGER(KIND=int4) :: xsize, ysize
    INTEGER(KIND=int4) :: xstride  
    type (AREA_STRUCT) :: AREAstr
    TYPE (GVAR_NAV), intent(in)    :: NAVstr_MTSAT
    
    INTEGER :: i, j, ii, jj, ierr, imode
    REAL(KIND(0.0d0)) :: latitude, longitude
    REAL(KIND=REAL4) :: elem, line, height
    REAL(KIND=REAL4) :: dlon, dlat
    REAL(KIND=REAL8) :: mjd
    REAL(KIND=REAL4), dimension(8) :: angles

    NAVstr_MTSAT_NAV = NAVstr_MTSAT
    
    imode = -1
    height = 0.0    !Used for parallax correction
    
    lat = Missing_Value_Real4
    lon = Missing_Value_Real4
       
    IF (NAVstr_MTSAT%nav_type == 'GMSX') THEN
    
        jj = 1 + (ystart - 1)*AREAstr%line_res
    
        DO j=1, ysize
            line = REAL(AREAstr%north_bound) + REAL(jj - 1) + &
                    REAL(AREAstr%line_res)/2.0
               
            DO i=1, xsize
                ii = ((i+(xstart-1)) - 1)*(AREAstr%elem_res*(xstride)) + 1
        
                elem = REAL(AREAstr%west_vis_pixel) + REAL(ii - 1) + &
	                   REAL(AREAstr%elem_res*(xstride))/2.0
                   
                CALL MGIVSR(imode,elem,line,dlon,dlat,height,&
                            angles,mjd,ierr)
                
                Space_Mask(i,j) = sym%SPACE
            
                IF (ierr == 0) THEN
                     Space_Mask(i,j) = sym%NO_SPACE
                     Lat_1b(i,j) = dlat
                     Lon_1b(i,j) = dlon
                ENDIF
            END DO
            
            jj = jj + AREAstr%line_res
        END DO
                        
    ENDIF

    IF (NAVstr_MTSAT%nav_type == 'GEOS') THEN      
        !HRIT requires actual line and element of being processed.
        ! Unlike MSG, MTSAT requires no switching to different corrdinates.
                
          DO j=1, ysize
!            jj = ystart + (j-1)
            jj = (AREAstr%north_bound / AREAstr%line_res) + ystart + (j-1)
    
               
            DO i=1, xsize
!                ii = (i - 1)*(xstride) + xstart	 ! get element of the image segement
                ii = (AREAstr%west_vis_pixel / AREAstr%elem_res) + (i - 1)*(xstride) + xstart	 ! get element of the image segement
                
                
                !CALL pixcoord2geocoord_mtsat(ii, jj, NAVstr_MTSAT%LOFF, NAVstr_MTSAT%COFF, &
                !     latitude,longitude, NAVstr_MTSAT%sub_lon, NAVstr_MTSAT%CFAC, NAVstr_MTSAT%LFAC)


                ! again, use common algorithm for CGMS navigation
                CALL pixcoord2geocoord_cgms(ii,                  &
                                            jj,                  &
                                            NAVstr_MTSAT%LOFF,   &
                                            NAVstr_MTSAT%COFF,   & 
                                            NAVstr_MTSAT%LFAC,   &
                                            NAVstr_MTSAT%CFAC,   &
                                            1,             &
                                            NAVstr_MTSAT%sub_lon, &
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
      
 END SUBROUTINE mtsat_navigation
 
 
!------------------------------------------------------------------
! SUBROUTINE to convert MTSAT counts to radiance and brightness
! temperature
!------------------------------------------------------------------
  
  SUBROUTINE MTSAT_RADIANCE_BT(chan_num,Mtsat_Counts, rad2, temp1)

    INTEGER (kind=INT2), dimension(:,:), intent(in):: Mtsat_Counts
    INTEGER (kind=int1), INTENT(in) :: chan_num
    REAL (kind=real4), DIMENSION(:,:), INTENT(out):: temp1, rad2
    
    INTEGER :: i, j, index
                                   
    DO j = 1, num_scans_read
      DO i = 1, num_pix
       
       index = int(Mtsat_Counts(i,j),KIND=int2) + 1
       
       IF ((Space_Mask(i,j) == sym%NO_SPACE) .AND. &
           (index .LE. 1024) .AND. (index .GE. 1)) THEN 
       !only do valid counts
          rad2(i,j) = REAL(rad_table(chan_num,1,index),KIND=REAL4)/1000.0
          temp1(i,j) = REAL(bt_table(chan_num,1,index),KIND=REAL4)/100.0                    
       ELSE
          rad2(i,j) = Missing_Value_Real4
          temp1(i,j) = Missing_Value_Real4
       ENDIF
      END DO
    END DO    
  
  END SUBROUTINE MTSAT_RADIANCE_BT
  


!---- MTSAT HiRID Navigation


SUBROUTINE MGIVSR(IMODE,RPIX,RLIN,RLON,RLAT,RHGT, &
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
 !   IMODE     I   I*4   CONVERSION MODE & IMAGE KIND
 !                         IMAGE KIND
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
 !   IRTN      O   I*4   RETURN CODE (0=O.K.)
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
       INTEGER, INTENT(in) :: IMODE
       REAL(KIND=REAL4), INTENT(inout) :: RPIX, RLIN, RLON, RLAT
       REAL(KIND=REAL4), INTENT(in) :: RHGT
       REAL(KIND=REAL4), dimension(8), INTENT(out) :: RINF
       REAL(KIND=REAL8), INTENT(out) :: DSCT
       INTEGER, INTENT(out) :: IRTN
       
       INTEGER :: LMODE
       REAL(KIND=REAL8) :: WKCOS, WKSIN
    
!       REAL*4     RPIX,RLIN,RLON,RLAT,RHGT,RINF(8)
 !     INTEGER*4  MAP(672,4)
 !
       REAL*4     EPS,RI0,RI,RJ,RSTEP,RSAMP,RFCL,RFCP,SENS,RFTL,RFTP
       REAL*4     RESLIN(4),RESELM(4),RLIC(4),RELMFC(4),SENSSU(4), &
                  VMIS(3),ELMIS(3,3),RLINE(4),RELMNT(4)
       REAL*8     BC,BETA,BS,CDR,CRD,DD,DDA,DDB,DDC,DEF,DK,DK1,DK2, &
                  DLAT,DLON,DPAI,DSPIN,DTIMS,EA,EE,EF,EN,HPAI,PC,PI,PS, &
                  QC,QS,RTIM,TF,TL,TP, &
                  SAT(3),SL(3),SLV(3),SP(3),SS(3),STN1(3),STN2(3), &
                  SX(3),SY(3),SW1(3),SW2(3),SW3(3)
       REAL*8     DSATZ,DSATA,DSUNZ,DSUNA,DSSDA,DSATD,SUNM,SDIS, &
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
       IF(ABS(IMODE).GT.4)  IRTN=1
       IF(ABS(RLAT).GT.90. .AND. IMODE.GT.0)  IRTN=2
       IF(IRTN.NE.0)  RETURN
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
         CALL  MGI100(RTIM,CDR,SAT,SP,SS,BETA)                      ![3.4]
 !-----------------------------------------------------------------------
         CALL  MGI220(SP,SS,SW1)                                    ![3.5]
         CALL  MGI220(SW1,SP,SW2)
         BC      = DCOS(BETA)
         BS      = DSIN(BETA)
         SW3(1)  = SW1(1)*BS+SW2(1)*BC
         SW3(2)  = SW1(2)*BS+SW2(2)*BC
         SW3(3)  = SW1(3)*BS+SW2(3)*BC
         CALL  MGI200(SW3,SX)
         CALL  MGI220(SP,SX,SY)
         SLV(1)  = STN1(1)-SAT(1)                                  ![3.6]
         SLV(2)  = STN1(2)-SAT(2)
         SLV(3)  = STN1(3)-SAT(3)
         CALL  MGI200(SLV,SL)                                      ![3.7]
         CALL  MGI210(SP,SL,SW2)
         CALL  MGI210(SY,SW2,SW3)
         CALL  MGI230(SY,SW2,TP)
         TF      = SP(1)*SW3(1)+SP(2)*SW3(2)+SP(3)*SW3(3)
         IF(TF.LT.0.D0)  TP=-TP
         CALL  MGI230(SP,SL,TL)
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
         ENDIF
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
         CALL  MGI100(RTIM,CDR,SAT,SP,SS,BETA)                     ![3.10]
         CALL  MGI220(SP,SS,SW1)                                   ![3.11]
         CALL  MGI220(SW1,SP,SW2)
         BC      = DCOS(BETA)
         BS      = DSIN(BETA)
         SW3(1)  = SW1(1)*BS+SW2(1)*BC
         SW3(2)  = SW1(2)*BS+SW2(2)*BC
         SW3(3)  = SW1(3)*BS+SW2(3)*BC
         CALL  MGI200(SW3,SX)
         CALL  MGI220(SP,SX,SY)
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
         CALL  MGI200(SW3,SL)                                      ![3.14]
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
         ENDIF
         IF(DABS(DK1).LE.DABS(DK2))  THEN
           DK    = DK1
         ELSE
           DK    = DK2
         ENDIF
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
           ENDIF
         ENDIF
         RLAT    = SNGL(DLAT*CRD)
         RLON    = SNGL(DLON*CRD)
         DSCT    = RTIM
       ENDIF
 !
 !!!!!!!!!!!!!!!!!! TRANSFORMATION (ZENITH/AZIMUTH) !!!!!!!!!!!!!!![3.16]
       STN2(1)   = DCOS(DLAT)*DCOS(DLON)                           ![3.17]
       STN2(2)   = DCOS(DLAT)*DSIN(DLON)
       STN2(3)   = DSIN(DLAT)
       SLV(1)    = SAT(1)-STN1(1)                                  ![3.18]
       SLV(2)    = SAT(2)-STN1(2)
       SLV(3)    = SAT(3)-STN1(3)
       CALL  MGI200(SLV,SL)
 !
       CALL  MGI230(STN2,SL,DSATZ)                                 ![3.19]
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
       ENDIF
       STN3(1) = DCOS(DLATN)*DCOS(DLONN)
       STN3(2) = DCOS(DLATN)*DSIN(DLONN)
       STN3(3) = DSIN(DLATN)
       SW1(1)  = SLV(1)+SS(1)*SDIS*1.D3                            ![3.22]
       SW1(2)  = SLV(2)+SS(2)*SDIS*1.D3
       SW1(3)  = SLV(3)+SS(3)*SDIS*1.D3
       CALL  MGI200(SW1,SW2)                                       ![3.23]
       CALL  MGI230(STN2,SW2,DSUNZ)
       CALL  MGI230(SL,SW2,DSSDA)                                  ![3.24]
       CALL  MGI240(SL,STN2,STN3,DPAI,DSATA)                       ![3.25]
       CALL  MGI240(SW2,STN2,STN3,DPAI,DSUNA)                      ![3.26]
       DSATD   = DSQRT(SLV(1)*SLV(1)+SLV(2)*SLV(2)+SLV(3)*SLV(3))  ![3.27]
 !
 !
       CALL  MGI200(STN1,SL)                                       ![3.28]
       CALL  MGI230(SW2,SL,DSUNG)
       CALL  MGI220(SL, SW2,SW3)
       CALL  MGI220(SW3,SL, SW1)
       WKCOS=DCOS(DSUNG)
       WKSIN=DSIN(DSUNG)
       SW2(1)=WKCOS*SL(1)-WKSIN*SW1(1)
       SW2(2)=WKCOS*SL(2)-WKSIN*SW1(2)
       SW2(3)=WKCOS*SL(3)-WKSIN*SW1(3)
       CALL  MGI230(SW2,SLV,DSUNG)
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
       RETURN
END SUBROUTINE MGIVSR

SUBROUTINE  MGI100(RTIM,CDR,SAT,SP,SS,BETA)
!       COMMON /MMAP1/MAP
       REAL*8    ATTALP,ATTDEL,BETA,CDR,DELT,RTIM,SITAGT,SUNALP,SUNDEL, &
                 WKCOS,WKSIN
       REAL*8    ATT1(3),ATT2(3),ATT3(3),NPA(3,3), &
                 SAT(3),SP(3),SS(3),ORBT2(35,8)
 !     INTEGER*4 MAP(672,4)
       INTEGER :: I
 !
 !      EQUIVALENCE (MAP(13,3),ORBT1(1,1))
 !      EQUIVALENCE (MAP(13,2),ATIT(1,1))
 !
       DO 1000 I=1,7
         IF(RTIM.GE.NAVstr_MTSAT_NAV%ORBT1(1,I).AND.RTIM.LT.NAVstr_MTSAT_NAV%ORBT1(1,I+1))  THEN
           CALL  MGI110 &
                (I,RTIM,CDR,NAVstr_MTSAT_NAV%ORBT1,ORBT2,SAT,SITAGT,SUNALP,SUNDEL,NPA)
           GO TO  1200
         ENDIF
  1000 CONTINUE
  1200 CONTINUE
 !
       DO 3000 I=1,9
         IF(RTIM.GE.NAVstr_MTSAT_NAV%ATIT(1,I) .AND. RTIM.LT.NAVstr_MTSAT_NAV%ATIT(1,I+1))  THEN
           DELT = (RTIM-NAVstr_MTSAT_NAV%ATIT(1,I))/(NAVstr_MTSAT_NAV%ATIT(1,I+1)-NAVstr_MTSAT_NAV%ATIT(1,I))
           ATTALP = NAVstr_MTSAT_NAV%ATIT(3,I)+(NAVstr_MTSAT_NAV%ATIT(3,I+1)-NAVstr_MTSAT_NAV%ATIT(3,I))*DELT
           ATTDEL = NAVstr_MTSAT_NAV%ATIT(4,I)+(NAVstr_MTSAT_NAV%ATIT(4,I+1)-NAVstr_MTSAT_NAV%ATIT(4,I))*DELT
           BETA   = NAVstr_MTSAT_NAV%ATIT(5,I)+(NAVstr_MTSAT_NAV%ATIT(5,I+1)-NAVstr_MTSAT_NAV%ATIT(5,I))*DELT
           IF( (NAVstr_MTSAT_NAV%ATIT(5,I+1)-NAVstr_MTSAT_NAV%ATIT(5,I)).GT.0.D0 ) &
             BETA   = NAVstr_MTSAT_NAV%ATIT(5,I)+(NAVstr_MTSAT_NAV%ATIT(5,I+1)-NAVstr_MTSAT_NAV%ATIT(5,I)-360.D0*CDR)*DELT
           GO TO  3001
         ENDIF
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
       CALL  MGI200(ATT3,SP)
 !
       WKCOS    = DCOS(SUNDEL)
       SS(1)    = WKCOS       *DCOS(SUNALP)
       SS(2)    = WKCOS       *DSIN(SUNALP)
       SS(3)    = DSIN(SUNDEL)
 !
       RETURN
END SUBROUTINE MGI100

SUBROUTINE MGI110(I,RTIM,CDR,ORBTA,ORBTB,SAT,SITAGT,SUNALP,SUNDEL,NPA)
       REAL*8    CDR,SAT(3),RTIM,ORBTA(35,8),ORBTB(35,8)
       REAL*8    SITAGT,SUNDEL,SUNALP,NPA(3,3),DELT
       INTEGER*4 I
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
       ENDIF
       RETURN
END SUBROUTINE MGI110

SUBROUTINE MGI200(VECT,VECTU)
       REAL*8  VECT(3),VECTU(3),RV1,RV2
       RV1=VECT(1)*VECT(1)+VECT(2)*VECT(2)+VECT(3)*VECT(3)
       IF(RV1 == 0.D0)  RETURN
       RV2=DSQRT(RV1)
       VECTU(1)=VECT(1)/RV2
       VECTU(2)=VECT(2)/RV2
       VECTU(3)=VECT(3)/RV2
       RETURN
END SUBROUTINE MGI200

SUBROUTINE MGI210(VA,VB,VC)
       REAL*8  VA(3),VB(3),VC(3)
       VC(1)= VA(2)*VB(3)-VA(3)*VB(2)
       VC(2)= VA(3)*VB(1)-VA(1)*VB(3)
       VC(3)= VA(1)*VB(2)-VA(2)*VB(1)
       RETURN
END SUBROUTINE MGI210

SUBROUTINE MGI220(VA,VB,VD)
       REAL*8  VA(3),VB(3),VC(3),VD(3)
       VC(1)= VA(2)*VB(3)-VA(3)*VB(2)
       VC(2)= VA(3)*VB(1)-VA(1)*VB(3)
       VC(3)= VA(1)*VB(2)-VA(2)*VB(1)
       CALL  MGI200(VC,VD)
       RETURN
END SUBROUTINE MGI220

SUBROUTINE MGI230(VA,VB,ASITA)
       REAL*8  VA(3),VB(3),ASITA,AS1,AS2
       AS1= VA(1)*VB(1)+VA(2)*VB(2)+VA(3)*VB(3)
       AS2=(VA(1)*VA(1)+VA(2)*VA(2)+VA(3)*VA(3))* &
           (VB(1)*VB(1)+VB(2)*VB(2)+VB(3)*VB(3))
       IF(AS2 == 0.D0)  RETURN
       ASITA=DACOS(AS1/DSQRT(AS2))
       RETURN
END SUBROUTINE MGI230

SUBROUTINE MGI240(VA,VH,VN,DPAI,AZI)
       REAL*8  VA(3),VH(3),VN(3),VB(3),VC(3),VD(3),DPAI,AZI,DNAI
       CALL  MGI220(VN,VH,VB)
       CALL  MGI220(VA,VH,VC)
       CALL  MGI230(VB,VC,AZI)
       CALL  MGI220(VB,VC,VD)
       DNAI = VD(1)*VH(1)+VD(2)*VH(2)+VD(3)*VH(3)
       IF(DNAI.GT.0.D0)  AZI=DPAI-AZI
       RETURN
END SUBROUTINE MGI240


!------------------- MTSAT NAV BLOC

 SUBROUTINE READ_NAVIGATION_BLOCK_MTSAT_FY(filename, AREAstr, NAVstr)
  CHARACTER(len=*), intent(in):: filename
  TYPE(AREA_STRUCT), intent(in):: AREAstr
  TYPE(GVAR_NAV), intent(inout):: NAVstr
 
  CHARACTER(len=1), dimension(3200) :: CBUF
  INTEGER :: i, j, geos_nav
  INTEGER(kind=int4)nav_offset
  REAL(kind=real4) :: R4DMY
  REAL(kind=real8) :: R8DMY  
! REAL(kind=real8) :: LOFF,COFF, LFAC, CFAC  
  integer:: number_of_words_read
  INTEGER(kind=int4), DIMENSION(640) :: i4buf
  
  nav_offset = AREAstr%sec_key_nav
    
  !determine GEOS or other navigation
  geos_nav = sym%NO
  CALL mreadf_int(trim(filename)//CHAR(0),nav_offset,4,640,&
                    number_of_words_read, i4buf)
!  IF (AREAstr%swap_bytes > 0) CALL swap_bytes4(i4buf,640)
  CALL move_bytes(4,i4buf(1),NAVstr%nav_type,0)
  
  IF (NAVstr%nav_type == 'GEOS') THEN
        !SUBLON stored as SUBLON *10 in McIDAS NAV block
        NAVstr%sub_lon = REAL(i4buf(6),kind=real4) / 10
        NAVstr%sublon = NAVstr%sub_lon

        ! LOFF, COFF, CFAC, LFAC stored in McIDAS header for 1km data. All
        ! Multipied by 10. Order from nvxmtst.dlm in McIDAS
        NAVstr%LOFF=(i4buf(2) / 10 ) / REAL(AREAstr%line_res)
        NAVstr%COFF=(i4buf(3) / 10) / REAL(AREAstr%elem_res)
        NAVstr%LFAC=(i4buf(4) / 10 ) / REAL(AREAstr%line_res)
        NAVstr%CFAC=(i4buf(5) / 10 ) / REAL(AREAstr%elem_res)
     
  END IF
    
  IF (NAVstr%nav_type == 'GMSX') THEN
    CALL mreadf_int(trim(filename)//CHAR(0),nav_offset,1,3200,&
                    number_of_words_read,NAVstr%COBAT)
        
    CBUF(1:504) = NAVstr%COBAT(5:508)
    CBUF(505:1012) = NAVstr%COBAT(513:1020)
    CBUF(1013:1520) = NAVstr%COBAT(1025:1532)
    CBUF(1521:2028) = NAVstr%COBAT(1537:2044)
    CBUF(2029:2536) = NAVstr%COBAT(2049:2556)
    
    CALL SV0100( 6, 8, CBUF( 1:   6), R4DMY, NAVstr%DTIMS )
    CALL SV0100( 4, 8, CBUF( 7:  10), NAVstr%RESLIN(1), R8DMY )
    CALL SV0100( 4, 8, CBUF( 11: 14), NAVstr%RESLIN(2), R8DMY )
    CALL SV0100( 4, 8, CBUF( 11: 14), NAVstr%RESLIN(3), R8DMY )
    CALL SV0100( 4, 8, CBUF( 11: 14), NAVstr%RESLIN(4), R8DMY )
    CALL SV0100( 4,10, CBUF( 15: 18), NAVstr%RESELM(1), R8DMY )
    CALL SV0100( 4,10, CBUF( 19: 22), NAVstr%RESELM(2), R8DMY )
    CALL SV0100( 4,10, CBUF( 19: 22), NAVstr%RESELM(3), R8DMY )
    CALL SV0100( 4,10, CBUF( 19: 22), NAVstr%RESELM(4), R8DMY )
    CALL SV0100( 4, 4, CBUF( 23: 26), NAVstr%RLIC(1), R8DMY )
    CALL SV0100( 4, 4, CBUF( 27: 30), NAVstr%RLIC(2), R8DMY )
    CALL SV0100( 4, 4, CBUF(111:114), NAVstr%RLIC(3), R8DMY )
    CALL SV0100( 4, 4, CBUF(115:118), NAVstr%RLIC(4), R8DMY )
    CALL SV0100( 4, 4, CBUF( 31: 34), NAVstr%RELMFC(1), R8DMY )
    CALL SV0100( 4, 4, CBUF( 35: 38), NAVstr%RELMFC(2), R8DMY )
    CALL SV0100( 4, 4, CBUF(119:122), NAVstr%RELMFC(3), R8DMY )
    CALL SV0100( 4, 4, CBUF(123:126), NAVstr%RELMFC(4), R8DMY )
    CALL SV0100( 4, 0, CBUF( 39: 42), NAVstr%SENSSU(1), R8DMY )
    CALL SV0100( 4, 0, CBUF( 43: 46), NAVstr%SENSSU(2), R8DMY )
    CALL SV0100( 4, 0, CBUF( 43: 46), NAVstr%SENSSU(3), R8DMY )
    CALL SV0100( 4, 0, CBUF( 43: 46), NAVstr%SENSSU(4), R8DMY )
    CALL SV0100( 4, 0, CBUF( 47: 50), NAVstr%RLINE(1) , R8DMY )
    CALL SV0100( 4, 0, CBUF( 51: 54), NAVstr%RLINE(2) , R8DMY )
    CALL SV0100( 4, 0, CBUF( 51: 54), NAVstr%RLINE(3) , R8DMY )
    CALL SV0100( 4, 0, CBUF( 51: 54), NAVstr%RLINE(4) , R8DMY )
    CALL SV0100( 4, 0, CBUF( 55: 58), NAVstr%RELMNT(1), R8DMY )
    CALL SV0100( 4, 0, CBUF( 59: 62), NAVstr%RELMNT(2), R8DMY )
    CALL SV0100( 4, 0, CBUF( 59: 62), NAVstr%RELMNT(3), R8DMY )
    CALL SV0100( 4, 0, CBUF( 59: 62), NAVstr%RELMNT(4), R8DMY )
    CALL SV0100( 4,10, CBUF( 63: 66), NAVstr%VMIS(1), R8DMY )
    CALL SV0100( 4,10, CBUF( 67: 70), NAVstr%VMIS(2), R8DMY )
    CALL SV0100( 4,10, CBUF( 71: 74), NAVstr%VMIS(3), R8DMY )
    CALL SV0100( 4, 7, CBUF( 75: 78), NAVstr%ELMIS(1,1), R8DMY )
    CALL SV0100( 4,10, CBUF( 79: 82), NAVstr%ELMIS(2,1), R8DMY )
    CALL SV0100( 4,10, CBUF( 83: 86), NAVstr%ELMIS(3,1), R8DMY )
    CALL SV0100( 4,10, CBUF( 87: 90), NAVstr%ELMIS(1,2), R8DMY )
    CALL SV0100( 4, 7, CBUF( 91: 94), NAVstr%ELMIS(2,2), R8DMY )
    CALL SV0100( 4,10, CBUF( 95: 98), NAVstr%ELMIS(3,2), R8DMY )
    CALL SV0100( 4,10, CBUF( 99:102), NAVstr%ELMIS(1,3), R8DMY )
    CALL SV0100( 4,10, CBUF(103:106), NAVstr%ELMIS(2,3), R8DMY )
    CALL SV0100( 4, 7, CBUF(107:110), NAVstr%ELMIS(3,3), R8DMY )
    CALL SV0100( 6, 8, CBUF(241:246), R4DMY, NAVstr%DSPIN )
    CALL SV0100( 6, 6, CBUF(199:204), NAVstr%sublon, R8DMY )
    CALL SV0100( 6, 6, CBUF(205:210), NAVstr%sublat, R8DMY )
    
    DO i=1, 10
        j = (i-1)*48 + 256
        CALL SV0100(6, 8,CBUF( 1+j: 6+j),R4DMY,NAVstr%ATIT(1,i))
        CALL SV0100(6, 8,CBUF(13+j:18+j),R4DMY,NAVstr%ATIT(3,i))
        CALL SV0100(6,11,CBUF(19+j:24+j),R4DMY,NAVstr%ATIT(4,i))
        CALL SV0100(6, 8,CBUF(25+j:30+j),R4DMY,NAVstr%ATIT(5,i))
        CALL SV0100(6, 8,CBUF(31+j:36+j),R4DMY,NAVstr%ATIT(6,i))
    END DO
    
    DO i=1, 8
        j = (i-1)*200 + (256 + 10*48)
        CALL SV0100(6, 8,CBUF(1+j:6+j),R4DMY,NAVstr%ORBT1( 1,i))
        CALL SV0100(6, 6,CBUF(49+j:54+j),R4DMY,NAVstr%ORBT1( 9,i))
        CALL SV0100(6, 6,CBUF(55+j:60+j),R4DMY,NAVstr%ORBT1(10,i))
        CALL SV0100(6, 6,CBUF(61+j:66+j),R4DMY,NAVstr%ORBT1(11,i))
        CALL SV0100(6, 8,CBUF(85+j:90+j),R4DMY,NAVstr%ORBT1(15,i))
        CALL SV0100(6, 8,CBUF(103+j:108+j),R4DMY,NAVstr%ORBT1(18,i))
        CALL SV0100(6, 8,CBUF(109+j:114+j),R4DMY,NAVstr%ORBT1(19,i))
        CALL SV0100(6,12,CBUF(129+j:134+j),R4DMY,NAVstr%ORBT1(20,i))
        CALL SV0100(6,14,CBUF(135+j:140+j),R4DMY,NAVstr%ORBT1(21,i))
        CALL SV0100(6,14,CBUF(141+j:146+j),R4DMY,NAVstr%ORBT1(22,i))
        CALL SV0100(6,14,CBUF(147+j:152+j),R4DMY,NAVstr%ORBT1(23,i))
        CALL SV0100(6,12,CBUF(153+j:158+j),R4DMY,NAVstr%ORBT1(24,i))
        CALL SV0100(6,16,CBUF(159+j:164+j),R4DMY,NAVstr%ORBT1(25,i))
        CALL SV0100(6,12,CBUF(165+j:170+j),R4DMY,NAVstr%ORBT1(26,i))
        CALL SV0100(6,16,CBUF(171+j:176+j),R4DMY,NAVstr%ORBT1(27,i))
        CALL SV0100(6,12,CBUF(177+j:182+j),R4DMY,NAVstr%ORBT1(28,i))
    END DO
    NAVstr%sub_lon = NAVstr%sublon    
  END IF

 END SUBROUTINE READ_NAVIGATION_BLOCK_MTSAT_FY


 SUBROUTINE  SV0100( IWORD, IPOS, C, R4DAT, R8DAT )
 !---- ------------------------------------------------------------------
 !     TYPE CONVERT ROUTINE ( R-TYPE )
 !---- ------------------------------------------------------------------
       INTEGER*4  IWORD,IPOS,IDATA1
       CHARACTER  C(*)*1
       REAL*4     R4DAT
       REAL*8     R8DAT
       R4DAT = 0.0
       R8DAT = 0.D0
       IF( IWORD == 4 )  THEN
         IDATA1 = ICHAR( C(1)(1:1) )/128
         R8DAT  = DFLOAT( MOD(ICHAR(C(1)(1:1)),128) )*2.D0**(8*3)+ &
                  DFLOAT( ICHAR(C(2)(1:1)) )*2.D0**(8*2)+ &
                  DFLOAT( ICHAR(C(3)(1:1)) )*2.D0**(8*1)+ &
                  DFLOAT( ICHAR(C(4)(1:1)) )
         R8DAT  = R8DAT/10.D0**IPOS
         IF( IDATA1 == 1 )  R8DAT = -R8DAT
         R4DAT  = SNGL( R8DAT )
       ELSEIF( IWORD == 6 )  THEN
         IDATA1 = ICHAR( C(1)(1:1) )/128
         R8DAT  = DFLOAT( MOD(ICHAR(C(1)(1:1)),128) )*2.D0**(8*5)+ &
                  DFLOAT( ICHAR(C(2)(1:1)) )*2.D0**(8*4)+ &
                  DFLOAT( ICHAR(C(3)(1:1)) )*2.D0**(8*3)+ &
                  DFLOAT( ICHAR(C(4)(1:1)) )*2.D0**(8*2)+ &
                  DFLOAT( ICHAR(C(5)(1:1)) )*2.D0**(8*1)+ &
                  DFLOAT( ICHAR(C(6)(1:1)) )
         R8DAT  = R8DAT/10.D0**IPOS
         IF( IDATA1 == 1 )  R8DAT = -R8DAT
         R4DAT  = SNGL( R8DAT )
       END IF
       RETURN
END SUBROUTINE SV0100

!---------------------------------------------------------------------------------------------
! Calibrate the Ch1 Dark Counts using the Mtsat Reflectance Calibration Function
!---------------------------------------------------------------------------------------------
subroutine CALIBRATE_MTSAT_DARK_COMPOSITE(Dark_Comp_Counts,Ref_Ch1_Dark)

integer(kind=int2), dimension(:,:), intent(in):: Dark_Comp_Counts
real(kind=real4), dimension(:,:),  intent(out):: Ref_Ch1_Dark

call MTSAT_REFLECTANCE(Dark_Comp_Counts,Ref_Ch1_Dark)

end subroutine CALIBRATE_MTSAT_DARK_COMPOSITE

!-----------------------------------------------------------------------------------------------------
! public routine to read data from an Ch1 Dark Composite area file for one segment into memory
!-----------------------------------------------------------------------------------------------------
subroutine READ_MTSAT_DARK_COMPOSITE(Segment_Number,Dark_Composite_Filename, &
                                    Time_Since_Launch,AREAstr_image,Ref_Ch1_Dark)

   integer(kind=int4), intent(in):: Segment_Number
   character(len=*), intent(in):: Dark_Composite_Filename
   TYPE (AREA_STRUCT), intent(in) :: AREAstr_image
   real(kind=real4), intent(in):: Time_Since_Launch
   real(kind=real4), dimension(:,:), intent(out):: Ref_Ch1_Dark
   integer:: io_status

   TYPE (AREA_STRUCT), save :: AREAstr_dark
   integer:: Dark_Lun_Header
   integer, save:: Dark_Lun_Data
   integer, save:: Element_Offset
   integer, save:: Line_Offset
   integer:: Store_Time_Flag
   integer:: Rec_Num
   integer:: First_Rec_Num
   integer:: Last_Rec_Num
   integer:: Num_Recs_To_Read
   integer:: Num_Elements
   integer, save:: Num_Elements_Dark
   integer:: First_Line_In_Segment
   integer:: Last_Line_In_Segment
   integer:: Num_Lines_Per_Segment
   integer:: Num_Lines_Read
   integer:: Big_Count
   integer:: iline

   integer(kind=int2), dimension(:), allocatable, save:: Dark_Comp_Counts_Temp
   integer(kind=int2), dimension(:), allocatable, save:: Dark_Comp_Counts

   real, save:: Ch1_Gain
   integer:: ielem
   integer:: ii
   integer :: index ! needed for MTSAT
   character(len=200):: System_String

   !--- aliases
   Num_Elements = Num_Pix
   Num_Lines_Per_Segment = Num_Scans_Per_Segment
   Num_Lines_Read = Num_Scans_Read
   Big_Count = 999
   io_status = 0

   !--- check to see if a dark composite file exists for this image
   if (trim(Dark_Composite_Filename) == "no_file") then
!         print *, "No Dark Composite Available for this Image"
          return
   endif

   !--- do not store time from the dark composite
   Store_Time_Flag = sym%NO


   !--- on the first segment, open the file and read the Area and Nav headers
   !--- and compute the offsets between the image and the dark composite
   if (Segment_Number == 1) then 


        !--- uncompress
        System_String = "gunzip -c "//trim(Dark_Comp_Data_Dir_Temp)// &
                        trim(Dark_Composite_Name)//".gz"// &
                        " > "//trim(Temporary_Data_Dir)//trim(Dark_Composite_Name)
        call system(System_String)

        Number_of_Temporary_Files = Number_of_Temporary_Files + 1
        Temporary_File_Name(Number_of_Temporary_Files) = trim(Dark_Composite_Name)

        !--- open for header read
        Dark_Lun_Header = Get_Lun()
        open(unit=Dark_Lun_Header,file=trim(Temporary_Data_Dir)//trim(Dark_Composite_Name), &
             form="unformatted", access="direct",recl=4, &
             status="old",action="read", iostat=io_status)

        read (unit=Dark_Lun_Header, rec = 1, iostat=io_status) AREAstr_dark%area_status
        if (io_status /= 0) then
           Dark_Composite_Name = "no_file"
           print *, "Dark Composite Area file cannot be read"  ! byte swapping may cause this
           return
        endif
        read (unit=Dark_Lun_Header, rec = 2) AREAstr_dark%version_num
        if (AREAstr_dark%version_num .ne. 4) then
           Dark_Composite_Name = "no_file"
           print *, "Dark Composite Area file cannot be read"  ! byte swapping may cause this
           return
        endif
        read (unit=Dark_Lun_Header, rec = 3) AREAstr_dark%sat_id_num
        read (unit=Dark_Lun_Header, rec = 4) AREAstr_dark%img_date
        read (unit=Dark_Lun_Header, rec = 5) AREAstr_dark%img_time
        read (unit=Dark_Lun_Header, rec = 6) AREAstr_dark%north_bound        !beg_sc
        read (unit=Dark_Lun_Header, rec = 7) AREAstr_dark%west_vis_pixel     !beg_el
        read (unit=Dark_Lun_Header, rec = 8) AREAstr_dark%z_coor
        read (unit=Dark_Lun_Header, rec = 9) AREAstr_dark%num_line
        read (unit=Dark_Lun_Header, rec = 10) AREAstr_dark%num_elem
        read (unit=Dark_Lun_Header, rec = 11) AREAstr_dark%bytes_per_pixel
        read (unit=Dark_Lun_Header, rec = 12) AREAstr_dark%line_res
        read (unit=Dark_Lun_Header, rec = 13) AREAstr_dark%elem_res

        close(unit=Dark_Lun_Header)

        !--------------------------------------------------------------
        ! check to see if dark image is consistent with image data
        !--------------------------------------------------------------
        if ((AREAstr_image%num_elem /= AREAstr_dark%num_elem) .or.    &
            (AREAstr_image%num_line /= AREAstr_dark%num_line)) then
            print *, EXE_PROMPT, "Dark Composite File Is Incompatible, Ignoring It"
            Dark_Composite_Name = "no_file"
            return
        endif


        Num_Elements_Dark = AREAstr_dark%num_elem
        allocate(Dark_Comp_Counts_Temp(Num_Elements_Dark))
        allocate(Dark_Comp_Counts(Num_Elements_Dark))

!print *, "computing x-offsets ", AREAstr_image%west_vis_pixel, AREAstr_dark%west_vis_pixel
!print *, "computing y-offsets ", AREAstr_image%north_bound, AREAstr_dark%north_bound

        !--- compute offsets

        !--- east-west - note this is an absolute value
        Element_Offset = 0
        if (AREAstr_dark%west_vis_pixel /= AREAstr_image%west_vis_pixel) then
                Element_Offset = abs(AREAstr_dark%west_vis_pixel  -    &
                                     AREAstr_image%west_vis_pixel) / &
                                     AREAstr_image%elem_res
        endif

        !--- north-south
        Line_Offset = 0
        if (AREAstr_dark%north_bound /= AREAstr_image%north_bound) then
                Line_Offset =  (AREAstr_dark%north_bound  -    &
                                AREAstr_image%north_bound) / &
                                AREAstr_image%line_res
        endif

        print *, "Shifts for this Dark Comp = ", Element_Offset, Line_Offset

        !--- determine the gain
        Ch1_Gain = Ch1_Gain_Low_0*(100.0+Ch1_Degrad_Low_1*Time_Since_Launch+ &
                   Ch1_Degrad_Low_2*Time_Since_Launch**2)/100.0

        !-- open file for data read - use int2
        Dark_Lun_Data=GET_LUN()
        open(unit=Dark_Lun_Data,file=trim(Temporary_Data_Dir)//trim(Dark_Composite_Name), &
             form="unformatted", access="direct",recl=2*Num_Elements_Dark, &
             status="old",action="read",iostat=io_status)

        if (io_status /= 0) then
            print *, EXE_PROMPT, "Dark Composite File Could Not be Reopened"
            Dark_Composite_Name = "no_file"
            return
        endif
       
   endif 

   !--- for this segment, compute the words to read in (accounting for offset)
   First_Line_In_Segment = (Segment_Number-1)*Num_Lines_Per_Segment + 1

   First_Line_In_Segment = First_Line_In_Segment + Line_Offset    !check this!!!!

   Last_Line_In_Segment = min(AREAstr_image%num_line,Segment_Number*Num_Lines_Per_Segment)


   !--- determine records to read - this accounts for north-south shift
   First_Rec_Num = 1  + First_Line_In_Segment + Line_Offset
   Last_Rec_Num = First_Rec_Num + Num_Scans_Read - 1
   First_Rec_Num = max(2,First_Rec_Num)
   Last_Rec_Num = min(Last_Rec_Num,Num_Scans)
   Num_Recs_To_Read = Last_Rec_Num - First_Rec_Num + 1

!  print *, "Records to Read for this Seg = ", First_Rec_Num, Last_Rec_Num, &
!           num_recs_to_read, Num_Scans_Read
   iline = 1
   do Rec_Num = First_Rec_Num, Last_Rec_Num

     !--- read data
     read(unit=Dark_Lun_Data,rec=Rec_Num,iostat=io_status) Dark_Comp_Counts_Temp
     if (io_status /= 0) then
          print *, EXE_PROMPT, "Dark Composite File Could Not be Read"
          Dark_Composite_Name = "no_file"
          return
     endif
  
     !-- apply east-west shift
     Dark_Comp_Counts = Dark_Comp_Counts_Temp

     if (AREAstr_dark%west_vis_pixel < AREAstr_image%west_vis_pixel) then
       Dark_Comp_Counts(1:Num_Elements_Dark - Element_Offset) = Dark_Comp_Counts_temp(Element_Offset+1:Num_Elements_Dark) 
       Dark_Comp_Counts(Num_Elements_Dark-Element_Offset+1:Num_Elements_Dark) = Big_Count
     endif

     if (AREAstr_dark%west_vis_pixel > AREAstr_image%west_vis_pixel) then
       Dark_Comp_Counts(1:Element_Offset) = Big_Count
       Dark_Comp_Counts(Element_Offset+1:Num_Elements_Dark) = Dark_Comp_Counts_temp(1:Num_Elements_Dark-Element_Offset) 
       Dark_Comp_Counts(Num_Elements_Dark-Element_Offset+1:Num_Elements_Dark) = Big_Count
     endif

     !-- apply x offset
     do ielem = 1, Num_Elements
        ii = min(Num_Elements_Dark,1 + Mtsat_Xstride*(ielem-1))
        if (Dark_Comp_Counts(ii) == Big_Count) then
           Ref_Ch1_Dark(ielem,iline) = Missing_Value_Real4
        else
           !Recall, MTSAT has 2 different flavors of reflectance calculation
           ! See MTSAT module for documentation
                
           IF ( calib_type == 'MVSH' .OR. &
                calib_type == 'HSVM') THEN
                      
                index = int(Dark_Comp_Counts(ii),kind=int2) + 1.0
                Ref_Ch1_Dark(ielem,iline) = 100.0 * &
                                ((0.000978494701531316*index) - &
                                  0.00197851402698724)
            ENDIF
                 
                 
             IF (AREAstr_image%calib_type == 'MVIS' .OR. &
                  AREAstr_image%calib_type == 'SIVM') THEN

                     index = int(Dark_Comp_Counts(ii),KIND=int2)
                     Ref_Ch1_Dark(ielem,iline) = &
                            (REAL(ref_table(2,1,index),KIND=real4) / 100.0) 

             ENDIF           
           
        endif
     enddo

     

     !-- increment iline
     iline = iline + 1

   enddo

   !--- close file
   if (Segment_Number == Num_Segments) then
      close(unit=Dark_Lun_Data)
      deallocate(Dark_Comp_Counts_Temp)
      deallocate(Dark_Comp_Counts)

      !--- delete uncompressed file
      System_String = "rm "//trim(Temporary_Data_Dir)//trim(Dark_Composite_Name)
      call system(System_String)

   endif

end subroutine READ_MTSAT_DARK_COMPOSITE

!
!--- end of module
!
END MODULE MTSAT_MODULE