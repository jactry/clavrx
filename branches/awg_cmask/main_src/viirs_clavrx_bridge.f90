! $Header$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: viirs_clavrx_bridge.f90 (src)
!       viirs_clavrx_bridge (program)
!
! PURPOSE: VIIRS reader bridge to clavr-x
!
! DESCRIPTION: This module deals with the clav-x global value world with the more 
!              hidden viirs data world.  
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
! REVISION HISTORY:   created      March 2013 (AW)
!            4 Dec 2013 -  added mapping table in comments (AW) 
!                       -  change mapping of modis 3 and 4
!
!            12 January 2014 AW: add lunar rela zimuth and scattering angle 
!                                computation for global pixel_common variable
!
! NOTES:
!  VIIRS  MODIS(CLAVRX)  Wavelength
!           mapping
!    M1    -     8    -    0.412
!    M2    -     9    -    0.445
!    M3    -     3    -    0.488
!    M4    -     4    -    0.555
!    M5    -     1    -    0.672
!    M6    -    15    -    0.746
!    M7    -     2    -    0.865
!    M8    -     5    -    1.240
!    M9    -    26    -    1.378
!    M10   -     6    -    1.610
!    M11   -     7    -    2.250
!    M12   -    20    -    3.700
!    M13   -    22    -    4.050
!    M14   -    29    -    8.550
!    M15   -    31    -   10.763
!    M16   -    32    -   12.013
!    I1    -    39    -    0.640
!    I2    -    40    -    0.865
!    I3    -    41    -    1.610
!    I4    -    42    -    3.740
!    I5    -    43    -   11.450
!    DNB   -    44    -    0.700
!
!--------------------------------------------------------------------------------------

module VIIRS_CLAVRX_BRIDGE 


   use Pixel_Common , only : &
        Image &
      , Sensor &
      , Geo &
      , Nav &
      , Scan_Number &
      , Ancil_Data_Dir & 
      , Cloud_Mask_Aux_Flag &
      , Cloud_Mask_Aux_Read_Flag &
      , Cld_Mask_Aux &
      , Cld_Type_Aux &
      , Cld_Phase_Aux &
      , Scan_Time &
      , Gap_Pixel_Mask &
      , Ch &
      , Ref_ChI1 &
      , Ref_ChI2 &
      , Ref_ChI3 &
      , Bt_ChI4 &
      , Bt_ChI5 &
      , Ref_Min_ChI1 &
      , Ref_Max_ChI1 &
      , Ref_Mean_ChI1 &
      , Ref_Uni_ChI1 &
      , Ref_Min_ChI2 &
      , Ref_Max_ChI2 &
      , Ref_Mean_ChI2 &
      , Ref_Uni_ChI2 & 
      , Ref_Min_ChI3 &
      , Ref_Max_ChI3 &
      , Ref_Mean_ChI3 &
      , Ref_Uni_ChI3 & 
      , Bt_Min_ChI4 &
      , Bt_Max_ChI4 &
      , Bt_Mean_ChI4 &
      , Bt_Uni_ChI4 &
      , Bt_Min_ChI5 &
      , Bt_Max_ChI5 &
      , Bt_Mean_ChI5 &
      , Bt_Uni_ChI5 
      

   use constants, only: &
      int4 &
		, sym &
		, Missing_Value_Real4
      
   use clavrx_message_module   
   
contains
   
   subroutine read_viirs_data ( segment_number ,  file_gmtco_base , error_out )
      use viirs_read_mod , only : &
          viirs_data_config &
           , viirs_data_out &
           , get_viirs_data 
   
      use planck
      use viewing_geometry_module , only: &
          glint_angle &
          , scattering_angle  &
          , relative_azimuth

      use calibration_constants, only: &
            Planck_Nu
      
      
      implicit none
      
     
      integer , intent(in) :: segment_number
      character(len=*), intent(in) :: file_gmtco_base
      integer(kind=int4), intent(out) :: error_out
      
      type ( viirs_data_config )  :: v_conf
      type ( viirs_data_out )  :: out
      integer :: modis_chn_list (16)
      integer :: modis_chn_list_iband (5)
      logical :: is_mband_on (16)
      logical :: is_iband_on (5)
      integer :: i_mband , i_iband
      integer :: y_start , c_seg_lines , c_seg_lines_iband
      integer :: i
      integer :: modis_chn
      
      
  
      error_out = 0
      ! - mapping modis to viirs
      !                 041 044 048 055  068  074   085 124 138 160 225  375  405  855  108  120
      !                 M1  M2   M3   M4  M5   M6   M7  M8  M9  M10 M11  M12  M13  M14  M15  M16  
      modis_chn_list = [ 8 , 9 , 3 , 4 , 1 , 15 , 2 , 5 , 26 , 6 , 7 , 20 , 22 , 29 , 31 , 32 ]
      modis_chn_list_iband = [ 39 , 40 , 41 , 42 , 43 ]
      is_mband_on = Sensor%Chan_On_Flag_Default ( modis_chn_list) == sym%YES
      is_iband_on = Sensor%Chan_On_Flag_Default ( modis_chn_list_iband ) == sym%YES
      
      y_start = ( segment_number -1 ) * Image%Number_Of_Lines_Per_Segment + 1
      c_seg_lines = min (  y_start + Image%Number_Of_Lines_Per_Segment -1 , Image%Number_Of_Lines )  - y_start  + 1
      
      ! - configure viirs interface
      v_conf % chan_on_rfl_mband = is_mband_on
      v_conf % chan_on_iband = is_iband_on
      v_conf % chan_on_dnb = Sensor%Chan_On_Flag_Default(44) == sym%YES
      v_conf % viirs_cloud_mask_on = cloud_mask_aux_flag /= sym%NO_AUX_CLOUD_MASK
      v_conf % viirs_cloud_type_on = cloud_mask_aux_flag /= sym%NO_AUX_CLOUD_MASK
      
      v_conf % offset = [ 1 , y_start]
      v_conf % count = [ Image%Number_Of_Elements  , c_seg_lines  ]
      v_conf % dir_1b = trim(Image%Level1b_Path)
      
      v_conf % Ancil_Data_Dir = trim(Ancil_Data_Dir)
      v_conf % file_gmtco_base =  trim(file_gmtco_base)

      v_conf % Nu_List = 0.0
      v_conf % Nu_List(12:16) = [Planck_Nu(20) , Planck_Nu(22) ,  &
                                Planck_Nu(29) , Planck_Nu(31) , Planck_Nu(32)]

      ! - read the data 
      call get_viirs_data ( v_conf, out )
     

      ! - output to clavrx global variables
      ! geo
      nav % lat_1b(:,1:c_seg_lines)    = out % geo % lat
      nav % lon_1b(:,1:c_seg_lines)    = out % geo % lon
      scan_time(1:c_seg_lines)   = out % geo % scan_time
      geo % sataz(:,1:c_seg_lines)     = out % geo % sataz
      geo % satzen(:,1:c_seg_lines)    = out % geo % satzen
      geo % solaz (:,1:c_seg_lines)    = out % geo % solaz 
      geo % solzen (:,1:c_seg_lines)   = out % geo % solzen 
      nav % ascend (1:c_seg_lines)     = out % geo % ascend
      
      geo % moon_phase_angle = out % geo % Moon_Phase_Angle
      ! rel azimuths  - these are all global variables
      call  COMPUTE_RELATIVE_AZIMUTH_VIIRS( geo % solaz , geo % sataz , geo % relaz )

      !--- compute the glint zenith angle
      geo % glintzen = glint_angle( geo % solzen , geo % satzen , geo % relaz )

      !--- compute the scattering angle
      geo % scatangle = scattering_angle( geo % solzen , geo % satzen , geo % relaz )

      ! gap
      gap_pixel_mask( : ,1:c_seg_lines) = 0
      where ( out % gap % mask )
         gap_pixel_mask( : ,1:c_seg_lines) = 1
      end where 
      
      
      ! - m-bands
      do i_mband = 1 , 16
         modis_chn = modis_chn_list (i_mband)
         if ( .not. out % mband ( i_mband ) % is_read ) then
            sensor % chan_on_flag_per_line (modis_chn ,1:c_seg_lines) = sym % no 
            cycle   
         end if
         
         if ( .not. is_mband_on(i_mband) .or. (size(out % mband (i_mband) % ref) < 1 &
              .and. size(out % mband (i_mband) % rad) < 1) ) cycle
         
         if ( i_mband <= 11 ) then
            ch(modis_chn) % Ref_Toa ( : ,1:c_seg_lines)  =  out % mband (i_mband) % ref   
         end if
         if ( i_mband >= 12 ) then
            
            ch ( modis_chn)  % Rad_Toa( : ,1:c_seg_lines) =   out % mband (i_mband) % rad
            call compute_bt_array ( ch(modis_chn)%bt_toa , ch(modis_chn)%rad_toa , modis_chn , missing_value_real4 )
         end if
         
      end do
      
      ! - i-bands
      do  i_iband = 1 , 5
         
            if ( .not. out % file_exists % svi_file_exists (i_iband)) then
                 ! - switch off chan_on in CLAVR-x if file is not there..
               Sensor%Chan_On_Flag_Default ( modis_chn_list_iband ) = sym % NO
               sensor % chan_on_flag_per_line (modis_chn_list_iband (i_iband) ,1:c_seg_lines) = sym % NO
               cycle
            end if
            
            
           if ( .not. out % iband ( i_iband ) % is_read ) then
            sensor % chan_on_flag_per_line (modis_chn_list_iband (i_iband) ,1:c_seg_lines) = sym % no 
          
            cycle   
         end if
         
         if ( .not. is_iband_on(i_iband) .or. (size(out % iband (i_iband) % ref) < 1 &
              .and. size(out % iband (i_iband) % bt) < 1) ) then    
              cycle
         end if
         
         c_seg_lines_iband  = 2 * c_seg_lines
         select case ( i_iband)
         case(1)
            Ref_Chi1( : ,1 : c_seg_lines_iband) = out % iband (i_iband) % ref           
         case(2)
            Ref_Chi2( : ,1 : c_seg_lines_iband) = out % iband (i_iband) % ref 
         case(3)
            Ref_Chi3( : ,1 : c_seg_lines_iband) = out % iband (i_iband) % ref 
         case(4)
            Bt_Chi4( : ,1 : c_seg_lines_iband) = out % iband (i_iband) % bt
         case(5)
            Bt_Chi5( : ,1 : c_seg_lines_iband) = out % iband (i_iband) % bt 
         end select 
      
      end do 
    
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES .and. size(out % dnb_mgrid % rad) > 1) then
         ch(44)%rad_toa( : ,1:c_seg_lines)  = out % dnb_mgrid % rad
         geo % lunzen( : ,1:c_seg_lines) = out % geo % lunzen
         geo % lunaz( : ,1:c_seg_lines) = out % geo % lunaz
         ch(44)%ref_toa( : ,1:c_seg_lines) = out % dnb_mgrid % ref
         geo % lunrelaz( : ,1:c_seg_lines) = Relative_Azimuth ( geo % lunaz( : ,1:c_seg_lines) &
                                                             , geo % sataz( : ,1:c_seg_lines) )
          
         !--- compute the scattering angle
         geo % scatangle_lunar( : ,1:c_seg_lines) = scattering_angle(  geo % lunzen( : ,1:c_seg_lines) &
                                                , geo % satzen( : ,1:c_seg_lines) &
                                                , geo % lunrelaz( : ,1:c_seg_lines) )
         geo % glintzen_lunar( : ,1:c_seg_lines) = glint_angle( geo % lunzen( : ,1:c_seg_lines) &
                                             , geo % satzen( : ,1:c_seg_lines) &
                                             , geo % lunrelaz( : ,1:c_seg_lines) )                                       
      end if
      
      ! -global variables which has to be set
      Image%Number_Of_Lines_Read_This_Segment = c_seg_lines
      do i = 1, Image%Number_Of_Lines_Per_Segment
         scan_number(i) = y_start + i - 1
      end do
      
      !- ascending  (global varaibel )
      nav % ascend = 0  
      do i = 1 , Image%Number_Of_Lines_Read_This_Segment - 1
         if ( nav % lat_1b(Image%Number_Of_Elements / 2 , i + 1) <= nav % lat_1b( Image%Number_Of_Elements / 2 , i ) ) nav % ascend ( i )  = 1
      end do
      ! --- fix for the last line Denis B.
      if ( nav % lat_1b(Image%Number_Of_Elements / 2 , Image%Number_Of_Lines_Read_This_Segment) <= &
           nav % lat_1b(Image%Number_Of_Elements / 2 , Image%Number_Of_Lines_Read_This_Segment - 1) ) &
            nav % ascend ( Image%Number_Of_Lines_Read_This_Segment ) = 1

      
      !---  statistics I-Band on M-band grid
      if ( is_iband_on( 1 ) ) call COMPUTE_IBAND_STATISTICS (Ref_ChI1 , Ref_Min_ChI1 , Ref_Max_ChI1 , Ref_Mean_ChI1 , Ref_Uni_ChI1)
      if ( is_iband_on( 2 ) ) call COMPUTE_IBAND_STATISTICS (Ref_ChI2 , Ref_Min_ChI2 , Ref_Max_ChI2 , Ref_Mean_ChI2 , Ref_Uni_ChI2) 
      if ( is_iband_on( 3 ) ) call COMPUTE_IBAND_STATISTICS (Ref_ChI3 , Ref_Min_ChI3 , Ref_Max_ChI3 , Ref_Mean_ChI3 , Ref_Uni_ChI3) 
      if ( is_iband_on( 4 ) ) call COMPUTE_IBAND_STATISTICS (Bt_ChI4  , Bt_Min_ChI4  , Bt_Max_ChI4  , Bt_Mean_ChI4  , Bt_Uni_ChI4 )
      if ( is_iband_on( 5 ) ) call COMPUTE_IBAND_STATISTICS (Bt_ChI5  , Bt_Min_ChI5  , Bt_Max_ChI5  , Bt_Mean_ChI5  , Bt_Uni_ChI5 )
  

      if ( v_conf % viirs_cloud_mask_on .and. size(out % prd % cld_mask) > 0 ) then
         Cld_Mask_Aux( : ,1 : c_seg_lines ) = out % prd % cld_mask
         Cloud_Mask_Aux_Read_Flag = 1
      else
         Cloud_Mask_Aux_Read_Flag = 0
      end if   

      if ( v_conf % viirs_cloud_type_on .and. size(out % prd % cld_type) > 0 ) then
         cld_type_aux( : ,1 : c_seg_lines) = out % prd % cld_type
         cld_phase_aux( : ,1 : c_seg_lines ) = out % prd % cld_phase
      end if   
      
      call out % dealloc ()
    
   end subroutine READ_VIIRS_DATA

   !----------------------------------------------------------------   
   !  this routine should be at a different place
   !----------------------------------------------------------------
   subroutine  COMPUTE_RELATIVE_AZIMUTH_VIIRS( ang1 , ang2, rel_az_out )
      real , dimension(:,:), intent(in) :: ang1
      real , dimension(:,:), intent(in) :: ang2
      real , dimension(:,:)  :: rel_az_out
   
      rel_az_out = abs (ang1 - ang2 )
      where ( rel_az_out > 180 )
          rel_az_out = 360.0 - rel_az_out
      end where
      rel_az_out = 180.0 - rel_az_out
     
   end subroutine COMPUTE_RELATIVE_AZIMUTH_VIIRS
   
   !----------------------------------------------------------------
   ! - iband has full file dimension of 6400 x1536
   ! - mband 3200 x 768
   !  - output of min_val ... is 3200 768
   !----------------------------------------------------------------
   subroutine COMPUTE_IBAND_STATISTICS ( iband_array , out_min_val, out_max_val , out_mean_val, out_std_val )
      implicit none
      real, dimension(:,:) , intent(in) :: iband_array
      real, dimension(:,:) , intent(out)  :: out_min_val, out_max_val , out_mean_val , out_std_val
      real, dimension(2,2) :: small_iband
      
      integer :: im , jm
      integer , dimension(2) ::  dim_m
      integer :: iband_x0, iband_x1 ,  iband_y0, iband_y1
     
      
      dim_m = shape ( out_min_val )
      
      out_min_val =-999.
      out_max_val =-999.
      out_mean_val =-999.
      
      do im = 1   , dim_m( 1)
        
         iband_x0 = ( im-1) * 2 +1
         iband_x1 = iband_x0 + 1
         do jm =1 , dim_m( 2 )
            iband_y0 = ( jm-1) * 2 +1
            iband_y1 = iband_y0 + 1
            small_iband = iband_array ( iband_x0 :  iband_x1 ,  iband_y0 :  iband_y1 )
            if ( minval ( small_iband ) > 0 ) then 
               out_min_val ( im, jm ) = minval ( small_iband )
               out_max_val ( im, jm ) = maxval ( small_iband )
               out_mean_val ( im, jm ) = sum ( small_iband ) / 4.0
               out_std_val ( im, jm ) = SQRT ( ( SUM ( (small_iband - out_mean_val(im, jm )) **2 ) ) / ( 4. - 1. ) ) 
            end if
         end do
      end do
   
   end subroutine COMPUTE_IBAND_STATISTICS

   !----------------------------------------------------------------
   ! read the VIIRS constants into memory
   !-----------------------------------------------------------------
   subroutine READ_VIIRS_INSTR_CONSTANTS(Instr_Const_file)
      use calibration_constants
      use file_tools , only: getlun
      
      implicit none
 
      character(len=*), intent(in):: Instr_Const_file
      integer:: ios0, erstat
      integer:: Instr_Const_lun

      Instr_Const_lun = GETLUN()

      open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)
      call mesg ("opening "//trim(Instr_Const_file), level = verb_lev % VERBOSE) 
      erstat = 0
      if (ios0 /= 0) then
         erstat = 19
         print *,  " VIIRS_CLAVRX_BRIDGE.f90: Error opening VIIRS constants file, ios0 = ", ios0
         stop 19
      end if

      read(unit=Instr_Const_lun,fmt="(a3)") sat_name
      read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
      read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
      read(unit=Instr_Const_lun,fmt=*) planck_a1(20), planck_a2(20),planck_nu(20)
      read(unit=Instr_Const_lun,fmt=*) planck_a1(22), planck_a2(22),planck_nu(22)
      read(unit=Instr_Const_lun,fmt=*) planck_a1(29), planck_a2(29),planck_nu(29)
      read(unit=Instr_Const_lun,fmt=*) planck_a1(31), planck_a2(31),planck_nu(31)
      read(unit=Instr_Const_lun,fmt=*) planck_a1(32), planck_a2(32),planck_nu(32)
      read(unit=Instr_Const_lun,fmt=*) planck_a1(42), planck_a2(42),planck_nu(42)
      read(unit=Instr_Const_lun,fmt=*) planck_a1(43), planck_a2(43),planck_nu(43)
      close(unit=Instr_Const_lun)
  
      !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
      Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

   end subroutine READ_VIIRS_INSTR_CONSTANTS

   !-----------------------------------------------------------------------------------------
   !  Get information from VIIRS 
   !-----------------------------------------------------------------------------------------
  
   subroutine READ_VIIRS_DATE_TIME ( Path, Infile &
                , Year , Doy , Start_Time , End_Time , Orbit , Orbit_Identifier &
                , End_Year, End_Doy )

      use VIIRS_READ_MOD, only : &                                                                                                                    
            READ_VIIRS_DATE_TIME_ATT
          

      ! Get the date & time from the file's name
      implicit none
      
      character(len=*), intent(in) :: Path
      character(len=*), intent(in) :: Infile
      integer, intent(inout) , optional :: Year
      integer, intent(inout)  , optional:: Doy    !day of year
      integer, intent(inout) , optional :: Start_Time  !millisec
      integer, intent(inout)  , optional:: End_Time    !millisec
      integer, intent(inout)  , optional:: Orbit
      character(38), intent(inout) , optional :: Orbit_Identifier
      integer , intent(inout) , optional :: End_Year
      integer, intent(inout)  , optional:: End_Doy    !day of year
  

      !--- call READ_VIIRS_DATE_TIME_ATT from module
      call READ_VIIRS_DATE_TIME_ATT (Path, Infile &
                , Year , Doy , Start_Time , End_Time , Orbit , Orbit_Identifier &
                , End_Year, End_Doy )

   end subroutine READ_VIIRS_DATE_TIME 

!---------------------------------------------------------------------------------
!  subroutine GET_NUMBER_OF_SCANS_FROM_VIIRS_BRIDGE ( Infile , Number_Of_Viirs_Lines , Error_Out )
!  it's asking to read number of scans in the viirs_read_mod
!---------------------------------------------------------------------------------
   SUBROUTINE GET_NUMBER_OF_SCANS_FROM_VIIRS_BRIDGE ( Infile , Number_Of_Viirs_Lines , Error_Out )
      use viirs_read_mod , only : &
          READ_NUMBER_OF_SCANS_FROM_VIIRS
   
      CHARACTER(Len=*), INTENT(IN) :: Infile  
      INTEGER(kind=int4), INTENT(OUT) :: Error_Out
      INTEGER(KIND=INT4), INTENT(OUT):: Number_of_Viirs_Lines

      error_out = 0
      call READ_NUMBER_OF_SCANS_FROM_VIIRS ( Infile , Number_of_Viirs_Lines , Error_Out )
     
   END SUBROUTINE GET_NUMBER_OF_SCANS_FROM_VIIRS_BRIDGE
   
!-------------------------------------------------
! subroutine JULIAN(iday,imonth,iyear,jday)
! compute julian day
! input:
!         iday - integer day
!         imonth - integer month
!         iyear - integer year (2 or four digits)
!         
! output : jday - julian day
!--------------------------------------------------
   subroutine JULIAN(iday,imonth,iyear,jday)

      !-- Computes julian day (1-365/366)
      integer, intent(in)::  iday,imonth,iyear
      integer, intent(out):: jday
      integer::  j
      integer, dimension(12)::  jmonth

      jmonth = reshape ((/31,28,31,30,31,30,31,31,30,31,30,31/),(/12/))
      
      jday = iday
      if (modulo(iyear,4) == 0) then
         jmonth(2)=29
         ! make CLAVR-x year 2100 ready !!   
         if ( modulo(iyear,100) == 0 .and. modulo(iyear,400) /= 0 ) then
            jmonth (2) = 28
         end if
      endif
      
      

      do j = 1,imonth-1
         jday = jday + jmonth(j)
      end do

   end subroutine JULIAN

end module viirs_clavrx_bridge

