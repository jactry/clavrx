!$Header$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: iff_clavrx_bridge.f90 (src)
!       IFF_CLAVRX_BRIDGE (program)
!
! PURPOSE: IFF read tool
!
! DESCRIPTION: This bridge controls reading of IFF data
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
! REVISON HISTORY: 
!  created    October 2013 (Denis B)
!--------------------------------------------------------------------------------------
module IFF_CLAVRX_BRIDGE

   use PIXEL_COMMON, only: &
       Iff_Viirs_Flag &
       , Iff_Modis_Flag &
       , Iff_Avhrr_Flag &
       , Iff_Gap_Mask &
       , Cld_Mask_Aux &
       , Cloud_Mask_Aux_Flag &
       , Cloud_Mask_Aux_Read_Flag &
       , Chan_On_Flag_Default &
       , Chan_On_Flag &
       , Num_Scans_Per_Segment &
       , Num_Scans &
       , Num_Scans_Read &
       , Scan_Number &
       , Start_Time &
       , End_Time &
       , Num_Pix &
       , Dir_1b &
       , Ancil_Data_Dir &
       , Lat_1b &
       , Lon_1b &
       , Scan_Time &
       , Sataz &
       , Satzen &
       , Solaz &
       , Solzen &
       , Relaz &
       , Glintzen &
       , Scatangle &
       , Ascend &
       , ch &
       , Bt_375um_Sounder &
       , Bt_11um_Sounder &
       , Bt_12um_Sounder
   use CONSTANTS
   use IFF_MODULE

   implicit none

contains

!----------------------------------------------------------------------
!   Subroutine to get dimentions from IFF 1b file through iff_module (Latitude)
!
   subroutine GET_IFF_DIMS_BRIDGE (File_Name, Nx, Ny)

    implicit none

    character(len=*), intent(in):: File_Name
    integer, intent(out):: Nx
    integer, intent(out):: Ny

    call GET_IFF_DIMS (File_Name, Nx, Ny)

   end subroutine GET_IFF_DIMS_BRIDGE

!----------------------------------------------------------------------
!   Subroutine to start IFF 1b file read
!
   subroutine READ_IFF_DATA (segment_number, iff_file, error_out)
      use PLANCK , only : &
          COMPUTE_BT_ARRAY
      use viewing_geometry_module , only: &
          RELATIVE_AZIMUTH &
          , GLINT_ANGLE &
          , SCATTERING_ANGLE
   
      implicit none

      integer , intent(in) :: segment_number
      character(len=*), intent(in) :: iff_file
      integer(kind=int4), intent(out) :: error_out

      type ( iff_data_config )  :: iff_conf
      type ( iff_data_out )  :: out
      integer :: num_chan = 36
      integer , dimension(36) :: modis_chn_list
      integer :: i_band
      integer :: y_start , c_seg_lines
      integer :: iline
      logical, dimension(36) :: is_band_on


      error_out = 0
      ! - modis ch to read
      modis_chn_list = [  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, &
                         16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, &
                         29, 30, 31, 32, 33, 34, 35, 36 ]
      is_band_on = Chan_On_Flag_Default ( modis_chn_list) == sym%YES

      y_start = ( segment_number - 1 ) * Num_Scans_Per_Segment + 1
      c_seg_lines = min (  y_start + Num_Scans_Per_Segment - 1 , Num_Scans )  &
                         - y_start  + 1

      ! - configure interface
      iff_conf % n_chan = num_chan
      iff_conf % chan_list = modis_chn_list
      iff_conf % chan_on = Chan_On_Flag_Default ( modis_chn_list ) == sym%YES
      iff_conf % iff_cloud_mask_on = Cloud_Mask_Aux_Flag /= sym%NO_AUX_CLOUD_MASK

      iff_conf % offset = [ 1 , y_start]
      iff_conf % count = [ num_pix , c_seg_lines  ]
      iff_conf % dir_1b = trim(dir_1b)

      iff_conf % Ancil_Data_Dir = trim(Ancil_Data_Dir)
      iff_conf % iff_file =  trim(iff_file)

      ! - read the data 
      call GET_IFF_DATA ( iff_conf, out )

      ! - output to clavrx global variables
      ! geo
      Lat_1b(:,1:c_seg_lines)    = out % geo % lat
      Lon_1b(:,1:c_seg_lines)    = out % geo % lon
      Scan_Time(1:c_seg_lines)   = out % geo % scan_time
      Sataz(:,1:c_seg_lines)     = out % geo % sataz
      Satzen(:,1:c_seg_lines)    = out % geo % satzen
      Solaz (:,1:c_seg_lines)    = out % geo % solaz
      Solzen (:,1:c_seg_lines)   = out % geo % solzen

      ! - compute relative azimuth
      Relaz = RELATIVE_AZIMUTH( Solaz , Sataz )

      ! - compute glint angle
      Glintzen = GLINT_ANGLE( Solzen , Satzen , Relaz )

      ! - compute the scattering angle
      Scatangle = SCATTERING_ANGLE( Solzen , Satzen , Relaz )

      !--- redefine solar azimuth to be consistent with avhrr
      where (out % geo % solaz /= Missing_Value_Real4)
         out % geo % solaz = 180.0 - abs(out % geo % solaz)
      end where

      !---- compute asc/des node, scan_number
      Num_Scans_Read = c_seg_lines
      do iline = 1, Num_Scans_Read
        Scan_Number(iline) = y_start + iline - 1
      enddo

      ! - ascending or descending
      Ascend = 0
      do iline = 1 , Num_Scans_Read - 1
         if ( out % geo % lat(Num_Pix / 2 , iline + 1) <= out % geo % lat( Num_Pix / 2 , iline ) ) &
             Ascend( iline )  = 1
      end do

      ! - save all channel data that were read to global
      do i_band = 1 , iff_conf % n_chan
         if ( .not. out % band ( i_band ) % is_read ) then
            Chan_On_Flag (modis_chn_list (i_band) ,1:c_seg_lines) = sym % no
            cycle
         end if

         if ( .not. is_band_on(i_band) .or. (size(out % band (i_band) % ref) < 1 &
              .and. size(out % band (i_band) % rad) < 1) ) cycle

         ! - ref channels
         if (i_band < 20 .or. i_band == 26) then
            ch(i_band)%Ref_Toa( : ,1:c_seg_lines)  = out % band (i_band) % ref

         ! - rad/bt channels
         elseif ((i_band >= 20 .and. i_band < 26) &
            .or. (i_band > 26 .and. i_band <=36)) then

            !!! - ATTN: USED UNASSIGNED CHANELS FOR HIRS (SPECIAL CASE) !!!
            if (Iff_Avhrr_Flag == sym%YES .and. i_band == 21) then ! 3.75um 
               call COMPUTE_BT_ARRAY ( Bt_375um_Sounder , out % band (21) % rad , &
                                       20 , missing_value_real4 )
            elseif (Iff_Avhrr_Flag == sym%YES .and. i_band == 22) then ! 11um
               call COMPUTE_BT_ARRAY ( Bt_11um_Sounder , out % band (22) % rad , &
                                       31 , missing_value_real4 )
            elseif (Iff_Avhrr_Flag == sym%YES .and. i_band == 29) then ! 12um
               call COMPUTE_BT_ARRAY ( Bt_12um_Sounder , out % band (29) % rad , &
                                       32 , missing_value_real4 )
            else ! the rest of channels are normal
               ch(i_band)%Rad_Toa( : ,1:c_seg_lines)  = out % band (i_band) % rad
               call COMPUTE_BT_ARRAY ( ch(i_band)%Bt_Toa , ch(i_band)%Rad_Toa , &
                                       i_band , missing_value_real4 )
            endif
            ! --- make Iff_Gap_Mask out of CRIS channel
            ! --- 0 = data in 13.3, 1=no data in 13.3
            if (Iff_Viirs_Flag == sym%YES .and. i_band == 33) then
               IFF_Gap_Mask = 0
               where (ch(33)%Bt_Toa == missing_value_real4) 
                  IFF_Gap_Mask = 1
               endwhere
            endif

         ! - skip all not listed channels
         else
            cycle
         endif 
      end do

      ! --- check if we need to read cloud mask aux
      if ( iff_conf % iff_cloud_mask_on .and. size(out % prd % cld_mask) > 0 ) then
         cld_mask_aux( : ,1 : c_seg_lines ) = out % prd % cld_mask
         Cloud_Mask_Aux_Read_Flag = 1
      else
         Cloud_Mask_Aux_Read_Flag = 0
      end if


      call out % dealloc ()


   end subroutine READ_IFF_DATA

!----------------------------------------------------------------
! read the IFF VIIRS + CrIS constants into memory
!-----------------------------------------------------------------
   subroutine READ_IFF_VIIRS_INSTR_CONSTANTS(Instr_Const_file)
      use calibration_constants
      use file_tools , only: getlun

      implicit none

      character(len=*), intent(in):: Instr_Const_file
      integer:: ios0, erstat
      integer:: Instr_Const_lun

      Instr_Const_lun = GETLUN()

      open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)
      print *, "opening ", trim(Instr_Const_file)
      erstat = 0
      if (ios0 /= 0) then
         erstat = 19
         print *, EXE_PROMPT, "Error opening IFF VIIRS constants file, ios0 = ", ios0
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
      read(unit=Instr_Const_lun,fmt=*) a1_33, a2_33,nu_33
      read(unit=Instr_Const_lun,fmt=*) a1_34, a2_34,nu_34
      read(unit=Instr_Const_lun,fmt=*) a1_35, a2_35,nu_35
      read(unit=Instr_Const_lun,fmt=*) a1_36, a2_36,nu_36
      read(unit=Instr_Const_lun,fmt=*) b1_day_mask,b2_day_mask,b3_day_mask,b4_day_mask
      close(unit=Instr_Const_lun)

      !-- convert solar flux in channel 20 to mean with units mW/m^2/cm^-1
      Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

   end subroutine READ_IFF_VIIRS_INSTR_CONSTANTS

!---------------------------------------------------------------------- 
! read the IFF AVHRR + HIRS constants into memory
!----------------------------------------------------------------------
   subroutine READ_IFF_AVHRR_INSTR_CONSTANTS(Instr_Const_file)
      use calibration_constants
      use file_tools , only: getlun

      implicit none

      character(len=*), intent(in):: Instr_Const_file
      integer:: ios0, erstat
      integer:: Instr_Const_lun

      Instr_Const_lun = GETLUN()

      open(unit=Instr_Const_lun,file=trim(Instr_Const_file),status="old",position="rewind",action="read",iostat=ios0)
      print *, "opening ", trim(Instr_Const_file)
      erstat = 0
      if (ios0 /= 0) then
         erstat = 19
         print *, EXE_PROMPT, "Error opening IFF AVHRR constants file, ios0 = ", ios0
         stop 19
      end if

      read(unit=Instr_Const_lun,fmt="(a3)") sat_name
      read(unit=Instr_Const_lun,fmt=*) Solar_Ch20
      read(unit=Instr_Const_lun,fmt=*) Ew_Ch20
      read(unit=Instr_Const_lun,fmt=*) a1_20, a2_20, nu_20
      read(unit=Instr_Const_lun,fmt=*) a1_31, a2_31, nu_31
      read(unit=Instr_Const_lun,fmt=*) a1_32, a2_32, nu_32
      read(unit=Instr_Const_lun,fmt=*) a1_36, a2_36, nu_36
      read(unit=Instr_Const_lun,fmt=*) a1_35, a2_35, nu_35
      read(unit=Instr_Const_lun,fmt=*) a1_34, a2_34, nu_34
      read(unit=Instr_Const_lun,fmt=*) a1_33, a2_33, nu_33
      read(unit=Instr_Const_lun,fmt=*) a1_30, a2_30, nu_30
      read(unit=Instr_Const_lun,fmt=*) a1_28, a2_28, nu_28
      read(unit=Instr_Const_lun,fmt=*) a1_27, a2_27, nu_27
      read(unit=Instr_Const_lun,fmt=*) a1_25, a2_25, nu_25
      read(unit=Instr_Const_lun,fmt=*) a1_24, a2_24, nu_24
      read(unit=Instr_Const_lun,fmt=*) a1_23, a2_23, nu_23
      read(unit=Instr_Const_lun,fmt=*) a1_21, a2_21, nu_21
      read(unit=Instr_Const_lun,fmt=*) b1_day_mask, b2_day_mask, b3_day_mask, b4_day_mask
      close(unit=Instr_Const_lun)

      !-- convert solar flux in channel 20 & 21 to mean with units mW/m^2/cm^-1
      Solar_Ch20_Nu = 1000.0 * Solar_Ch20 / Ew_Ch20

   end subroutine READ_IFF_AVHRR_INSTR_CONSTANTS

!----------------------------------------------------------------------

end module IFF_CLAVRX_BRIDGE

