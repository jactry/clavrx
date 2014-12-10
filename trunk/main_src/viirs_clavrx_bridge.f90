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
!    I1    -    37    -    0.640
!    I2    -    38    -    0.865
!    I3    -    39    -    1.610
!    I4    -    40    -    3.740
!    I5    -    41    -   11.450
!    DNB   -    42    -    0.700
!
!--------------------------------------------------------------------------------------

module VIIRS_CLAVRX_BRIDGE 


   use Pixel_Common , only : &
      Chan_On_Flag_Default &
      , Chan_On_Flag &
      , Num_Scans_Per_Segment &
      , Num_Scans &
      , Num_Scans_Read &
      , Num_Pix &
      , Scan_Number &
      , Dir_1b &
      , Ancil_Data_Dir & 
      , Cloud_Mask_Aux_Flag &
      , Cloud_Mask_Aux_Read_Flag &
      , Cld_Mask_Aux &
      , Cld_Type_Aux &
      , Cld_Phase_Aux &
      , Lat_1b &
      , Lon_1b &
      , Scan_Time &
      , Sataz &
      , Satzen &
      , Solaz &
      , Solzen &
      , Ascend &
      , Moon_Phase_Angle &
      , Relaz &
      , Glintzen &
      , Lunzen &
      , Lunaz &
      , Lunrelaz &
      , Scatangle_Lunar &
      , Scatangle &
      , Glintzen_Lunar &
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
      int4
      
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
          Nu_20 &
          , Nu_22 &
          , Nu_29 &
          , Nu_31 &
          , Nu_32
      
      
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
      modis_chn_list_iband = [ 37 , 38 , 39 , 40 , 41 ]
      is_mband_on = Chan_On_Flag_Default ( modis_chn_list) == sym%YES
      is_iband_on = Chan_On_Flag_Default ( modis_chn_list_iband ) == sym%YES
      
      y_start = ( segment_number -1 ) * num_scans_per_segment + 1
      c_seg_lines = min (  y_start + num_scans_per_segment -1 , num_scans )  - y_start  + 1
      
      ! - configure viirs interface
      v_conf % chan_on_rfl_mband = is_mband_on
      v_conf % chan_on_iband = is_iband_on
      v_conf % chan_on_dnb = Chan_On_Flag_Default(42) == sym%YES
      v_conf % viirs_cloud_mask_on = cloud_mask_aux_flag /= sym%NO_AUX_CLOUD_MASK
      v_conf % viirs_cloud_type_on = cloud_mask_aux_flag /= sym%NO_AUX_CLOUD_MASK
      
      v_conf % offset = [ 1 , y_start]
      v_conf % count = [ num_pix  , c_seg_lines  ]
      v_conf % dir_1b = trim(dir_1b)
      
      v_conf % Ancil_Data_Dir = trim(Ancil_Data_Dir)
      v_conf % file_gmtco_base =  trim(file_gmtco_base)

      v_conf % Nu_List = 0.0
      v_conf % Nu_List(12:16) = [Nu_20 , Nu_22 , Nu_29 , Nu_31 , Nu_32]

      ! - read the data 
      call get_viirs_data ( v_conf, out )
     

      ! - output to clavrx global variables
      ! geo
      lat_1b(:,1:c_seg_lines)    = out % geo % lat
      lon_1b(:,1:c_seg_lines)    = out % geo % lon
      scan_time(1:c_seg_lines)   = out % geo % scan_time
      sataz(:,1:c_seg_lines)     = out % geo % sataz
      satzen(:,1:c_seg_lines)    = out % geo % satzen
      solaz (:,1:c_seg_lines)    = out % geo % solaz 
      solzen (:,1:c_seg_lines)   = out % geo % solzen 
      ascend (1:c_seg_lines)     = out % geo % ascend
      
      moon_phase_angle = out % geo % Moon_Phase_Angle
      ! rel azimuths  - these are all global variables
      call  COMPUTE_RELATIVE_AZIMUTH_VIIRS( solaz , sataz , relaz )

      !--- compute the glint zenith angle
      glintzen = glint_angle( solzen , satzen , relaz )

      !--- compute the scattering angle
      scatangle = scattering_angle( solzen , satzen , relaz )

      ! gap
      gap_pixel_mask( : ,1:c_seg_lines) = 0
      where ( out % gap % mask )
         gap_pixel_mask( : ,1:c_seg_lines) = 1
      end where 
      
      
      ! - m-bands
      do i_mband = 1 , 16
         modis_chn = modis_chn_list (i_mband)
         if ( .not. out % mband ( i_mband ) % is_read ) then
            chan_on_flag (modis_chn ,1:c_seg_lines) = sym % no 
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
           if ( .not. out % iband ( i_iband ) % is_read ) then
            chan_on_flag (modis_chn_list_iband (i_iband) ,1:c_seg_lines) = sym % no 
            cycle   
         end if
         if ( .not. is_iband_on(i_iband) .or. (size(out % iband (i_iband) % ref) < 1 &
              .and. size(out % iband (i_iband) % bt) < 1) ) cycle
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
    
      if ( Chan_On_Flag_Default(42) == sym%YES .and. size(out % dnb_mgrid % rad) > 0) then
         ch(42)%rad_toa( : ,1:c_seg_lines)  = out % dnb_mgrid % rad
         lunzen( : ,1:c_seg_lines) = out % geo % lunzen
         lunaz( : ,1:c_seg_lines) = out % geo % lunaz
         ch(42)%ref_toa( : ,1:c_seg_lines) = out % dnb_mgrid % ref
         lunrelaz( : ,1:c_seg_lines) = Relative_Azimuth ( lunaz( : ,1:c_seg_lines) &
                                                       , sataz( : ,1:c_seg_lines) )
          
         !--- compute the scattering angle
         scatangle_lunar( : ,1:c_seg_lines) = scattering_angle(  lunzen( : ,1:c_seg_lines) &
                                                , satzen( : ,1:c_seg_lines) &
                                                , lunrelaz( : ,1:c_seg_lines) )
         glintzen_lunar( : ,1:c_seg_lines) = glint_angle( lunzen( : ,1:c_seg_lines) &
                                             , satzen( : ,1:c_seg_lines) &
                                             , lunrelaz( : ,1:c_seg_lines) )                                       
      end if
      
      ! -global variables which has to be set
       Num_Scans_Read = c_seg_lines
      do i = 1, num_scans_per_segment
         scan_number(i) = y_start + i - 1
      end do
      
      !- ascending  (global varaibel )
      ascend = 0  
      do i = 1 , Num_Scans_Read - 1
         if ( lat_1b(num_pix / 2 , i + 1) <= lat_1b( num_pix / 2 , i ) ) ascend ( i )  = 1
      end do
      
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
         print *, EXE_PROMPT, "Error opening VIIRS constants file, ios0 = ", ios0
         stop 19
      end if

      read(unit=Instr_Const_lun,fmt="(a3)") sat_name
      read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
      read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
      read(unit=Instr_Const_lun,fmt=*) a1_20, a2_20,nu_20
      read(unit=Instr_Const_lun,fmt=*) a1_22, a2_22,nu_22
      read(unit=Instr_Const_lun,fmt=*) a1_29, a2_29,nu_29
      read(unit=Instr_Const_lun,fmt=*) a1_31, a2_31,nu_31
      read(unit=Instr_Const_lun,fmt=*) a1_32, a2_32,nu_32
      read(unit=Instr_Const_lun,fmt=*) a1_40, a2_40,nu_40
      read(unit=Instr_Const_lun,fmt=*) a1_41, a2_41,nu_41
      read(unit=Instr_Const_lun,fmt=*) b1_day_mask,b2_day_mask,b3_day_mask,b4_day_mask
      close(unit=Instr_Const_lun)
  
      !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
      Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

   end subroutine READ_VIIRS_INSTR_CONSTANTS

 !-----------------------------------------------------------------------------------------
   !  Extract time information from VIIRS filename - should explore use of header for this
   !   assumingly called from outside  ...
   !-----------------------------------------------------------------------------------------
  
   subroutine READ_VIIRS_DATE_TIME ( infile &
                , year , doy , start_time , end_time , orbit , orbit_identifier , end_year, end_doy )
      ! Get the date & time from the file's name
      implicit none
      
      character(len=*), intent(in) :: infile   
      integer, intent(out) , optional :: year
      integer, intent(out)  , optional:: doy    !day of year
      integer, intent(out) , optional :: start_time  !millisec
      integer, intent(out)  , optional:: end_time    !millisec
      integer, intent(out)  , optional:: orbit
      character(38), intent(out) , optional :: orbit_identifier
      integer , intent(out) , optional :: end_year
      integer, intent(out)  , optional:: end_doy    !day of year
  
      integer :: month
      integer :: day
      integer :: start_hour
      integer :: start_minute
      integer :: start_sec

      integer :: days_of_year
      integer :: end_hour
      integer :: end_minute
      integer :: end_sec
      integer:: year_loc
      integer:: doy_loc    !day of year
      integer:: end_year_loc
      integer:: end_doy_loc    !day of year
      integer :: start_time_loc  !millisec
      integer:: end_time_loc    !millisec
      integer:: orbit_loc
      character(38):: orbit_identifier_loc
 
  
      !         1	    2	      3	        4	  5	    6	      7	        8
      !12345678901234567890123456789012345678901234567890123456789012345678901234567890
      !GMODO_npp_d20100906_t2110510_e2112156_b00012_c20110707160532497848_noaa_ops.h5

      ! --- Read data from the file name
      read(Infile(12:15), fmt="(I4)") year_loc
      read(Infile(16:17), fmt="(I2)") month
      read(Infile(18:19), fmt="(I2)") day
      read(Infile(22:23), fmt="(I2)") start_hour
      read(Infile(24:25), fmt="(I2)") start_minute
      read(Infile(26:27), fmt="(I2)") start_sec
      read(Infile(31:32), fmt="(I2)") end_hour
      read(Infile(33:34), fmt="(I2)") end_minute
      read(Infile(35:36), fmt="(I2)") end_sec
      read(Infile(40:44), fmt="(I5)") orbit_loc

      !---- store orbit number
      orbit_identifier_loc = infile(7:44)

      !--- compute day of year
      call JULIAN(day,month,year_loc,doy_loc)

      ! --- Calculate start and end time
      start_time_loc = ((start_hour * 60 + start_minute) * 60 + start_sec) * 1000
      end_time_loc = ((end_hour * 60 + end_minute) * 60 + end_sec) * 1000
  
      end_doy_loc = doy_loc
      end_year_loc = year_loc
      if ( end_time_loc <= start_time_loc) then
         end_doy_loc = end_doy_loc + 1
         days_of_year = 365
         if ( modulo(year_loc,4) == 0)  days_of_year = 366
         if ( end_doy_loc > days_of_year) then
            end_doy_loc = 1
            end_year_loc = end_year_loc + 1
         end if  
      end if
      
      if ( present ( year)) year = year_loc
      if ( present ( doy)) doy = doy_loc
      if ( present ( start_time)) start_time = start_time_loc
      if ( present ( end_time)) end_time = end_time_loc
      if ( present ( orbit)) orbit = orbit_loc
      if ( present ( orbit_identifier)) orbit_identifier = orbit_identifier_loc
      if ( present ( end_year)) end_year = end_year_loc
      if ( present ( end_doy)) end_doy = end_doy_loc
      
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
        endif

        do j = 1,imonth-1
           jday = jday + jmonth(j)
        end do

   end subroutine JULIAN

end module viirs_clavrx_bridge

