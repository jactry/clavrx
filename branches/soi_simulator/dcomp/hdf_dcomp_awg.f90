! $Id: hdf_dcomp_awg.f90 5 2014-01-21 23:42:33Z awalther $
module HDF_DCOMP_AWG
!--------------------------------------------------------------------------
!  Module name: HDF_DCOMP_AWG
!
!  Description:
!  This module contains routines for hdf read functions for DCOMP algorithm
!
!    A. Walther 03/28/2010
!    andi.walther@ssec.wisc.edu
!
!
!
!--------------------------------------------------------------------------


 
   implicit  none
 
   private
 
   ! some params 
   integer, parameter, private:: int1 = selected_int_kind(1)
   integer, parameter, private:: real4 = selected_real_kind(6,37)
   character(*), parameter, private :: exe_prompt = "hdf_dcomp_awg>> "
   integer, parameter, private:: dfacc_read = 1
 
   integer,parameter,private::flag_yes = 1
   integer,parameter,private::flag_no = 0
 
   integer(kind=int1),parameter , private:: water_phase = 1
   integer(kind=int1),parameter , private:: ice_phase = 2
 
   integer(kind=int1),parameter , private:: num_phase = 2
 
   integer(kind=int1), parameter, private :: vis_chn_number = 1
   integer(kind=int1), parameter, private :: ir_chn_number_3a = 2      !--- avhrr
   integer(kind=int1), parameter, private :: ir_chn_number_3b = 3      !--- avhrr

   integer (kind=int1), parameter, private :: num_sat_zen = 45              !- size of zenith dimension
   integer (kind=int1), parameter, private :: num_sol_zen = 45              !- size of solar zenith dimension
   integer (kind=int1), parameter, private :: num_rel_azi = 45              !- size of relative az. dimension
   integer (kind=int1), parameter, private :: num_cod = 29
   integer (kind=int1), parameter, private :: num_cps = 9
 
   ! -- look-up-tables structure for cloud parameters 
   type :: reflectance_lookup_table
      integer (kind=int1) :: flag                          !- flag whether populated or not
      real (kind=real4), dimension(num_sat_zen) :: sat_zen            !- satellite zenith vector
      real (kind=real4), dimension(num_sol_zen) :: sol_zen            !- solar zenith vector
      real (kind=real4), dimension(num_rel_azi) :: rel_azi            !- relative azimuth vector
      real (kind=real4), dimension(num_cod) :: cod_vec                !- log10 optical depth vector
      real (kind=real4), dimension(num_cps):: cps_vec                 !- log10 effective radius vector
      real (kind=real4), dimension(num_cps,num_cod) :: sph_alb        !- spherical albedo
      real (kind=real4), dimension(num_cps,num_cod,num_sol_zen) :: trans       !- flux transmission
      real (kind=real4), dimension(num_cps,num_cod,num_sol_zen) :: cld_alb     !- cloud albedo
      real (kind=real4), dimension(num_cps,num_cod,num_sol_zen,num_sat_zen,num_rel_azi) :: refl  !- reflectance
   end type reflectance_lookup_table


   type :: phase_reflectance_table
      integer (kind=int1) :: flag                                 !- phase flag 
      type(reflectance_lookup_table),dimension(num_phase) :: phase  !- cloud phase
   end type phase_reflectance_table

   type :: channel_reflectance_tables1
      integer (kind=int1) :: flag                                  !-channel flag
      type(phase_reflectance_table), dimension(3) :: channel
   end type channel_reflectance_tables1
  
   type (channel_reflectance_tables1), save, target:: cld_refl_lut

   !-----------------------------------------------------------------------
   ! structures for emissivity tables
   !-----------------------------------------------------------------------
   type, private :: emissivity_lookup_table
      integer (kind=int1) :: flag
      real (kind=real4), dimension(num_cps,num_cod,num_sat_zen) :: trans 
      real (kind=real4), dimension(num_cps,num_cod,num_sat_zen) :: emiss 
   end type emissivity_lookup_table

   type, private:: phase_emissivity_tables
      integer (kind=int1):: flag
      type(emissivity_lookup_table),dimension(num_phase):: phase
   end type phase_emissivity_tables

   Type (Phase_Emissivity_Tables), private, SAVE, TARGET:: cld_ems_LUT

   !-- ancil data structure
   type :: ancil_data_type
      integer(kind=int1)::flag
      real(kind=real4),dimension(3,3) :: gas_trans  ! -- gas transmission coefficients
      real(kind=real4),dimension(3) :: ozone_trans        ! -- ozone transmission coefficients 
   end type ancil_data_type
 
   type(ancil_data_type),save :: ancil_data
 
   public  :: populate_all_lut
   private :: populate_cloud_lut_single
   private :: populate_ancil
   public  :: get_lut_data

contains
 
!====================================================================
! Function Name: POPULATE_ALL_LUT
!
! Function:     Populates all main LUTs and ancil LUT
!
! Inputs: lut_path - look-up table path
!         sensor_name - sensor name as appears in LUT file name
!         sensor_chn_number - channel number as appears in LUT file name
!         noaa_sat_string - noaa satellite string as appears in LUT file name
!
! Dependencies:
!         POPULATE_CLOUD_LUT_SINGLE(Lut_File,idx_chn,idx_phase)
!         CALL POPULATE_CLOUD_LUT_SINGLE_EMS(Lut_File,idx_phase)
!         POPULATE_ANCIL(ancil_file)
!
! History:
!
! Local Variables:
!         Phase_String - "wat" or "ice"
!         Chan_String - Channel string
!         Lut_File - LUT file string
!         Ancil_File - Ancillary file string
!         Idx_Chn, Idx_Phase, Identifier_old, Identifier_current
!
!====================================================================
SUBROUTINE POPULATE_ALL_LUT( &
      Lut_Path &          ! - input 
    , Sensor_Name &       ! - input
        , Sensor_Chn_Number & ! - input
		, Sensor_Chn_String &
		, geocat_flag )  ! - input

CHARACTER (LEN=*), INTENT(in) ::  lut_path  ! - path there LUT is stored
CHARACTER (LEN=*), INTENT(in):: sensor_name ! - sensor name as appeared in LUT file name
INTEGER, INTENT(in), DIMENSION(:)  :: Sensor_Chn_Number

CHARACTER (LEN=*), INTENT(in), DIMENSION(:),OPTIONAL  :: Sensor_Chn_String            ! - channel number as appeared in LUT file name
INTEGER(KIND=int1) , INTENT(in), DIMENSION(:),OPTIONAL  :: geocat_Flag 
! CHARACTER (LEN=*), INTENT(in)  :: NOAA_sat_string            ! - channel number as appeared in LUT file name

!-locals
CHARACTER (LEN =3 ) :: Phase_String ! - "wat" or "ice"
CHARACTER (LEN =3 ) :: Chan_String  !
CHARACTER (LEN=1024):: Lut_File     !
CHARACTER (LEN=1024):: Ancil_File   !
INTEGER::Idx_Chn
INTEGER::Idx_Phase
CHARACTER(LEN=1024)::Identifier_old
CHARACTER(LEN=1024)::Identifier_current

!- executable

Identifier_current = TRIM(Lut_Path)//TRIM(Sensor_Name)

if ( trim(sensor_name) == 'MODIS-AQUA' .or. trim(sensor_name) == 'MODIS-TERRA' ) then
  Identifier_current = TRIM(Lut_Path)//'MODIS'
end if

if ( trim(sensor_name) == 'MTSAT-1R' .or. trim(sensor_name) == 'MTSAT-2' ) then
  Identifier_current = TRIM(Lut_Path)//'MTSAT'
end if

IF (Identifier_current == Identifier_old) THEN
   PRINT*,TRIM(Identifier_current)//' already populated!'
   RETURN
ELSE

! PRINT "(a,'LUT  DCOMP starts')",EXE_PROMPT

  Loop_Channel: DO idx_chn = 1 , 3 
    !--- create a string to hold channel number
    IF (PRESENT(Sensor_Chn_String)) THEN
	   Chan_String= trim(sensor_chn_String(Idx_Chn)) 
	ELSE
	  WRITE(Chan_String, '(I1.1)') sensor_chn_number(Idx_Chn)
	END IF   
	
    IF (Chan_String .EQV. 'no') CYCLE 
    
	Loop_Phase: DO idx_phase = 1 , 2 
      IF (idx_phase .EQV. 1 ) THEN
            Phase_String = "wat"
            Lut_File = trim(Identifier_current)//"_ch" &
            //TRIM(Chan_String)//"_ref_lut_"//Phase_String//"_cld.hdf"
          ENDIF  
          
      IF (idx_phase .eqv. 2 ) THEN
            Phase_String = "ice"               
        Lut_File = trim(Identifier_current)//"_ch" &
            //TRIM(Chan_String)//"_ref_lut_"//Phase_String//"_cld.hdf"
      
          ENDIF
         
          CALL POPULATE_CLOUD_LUT_SINGLE(Lut_File,idx_chn,idx_phase)    
    END DO Loop_Phase 
  END DO Loop_channel

  Emm_loop_phase: DO idx_phase = 1 ,2
    IF (idx_phase .EQV. 1 ) THEN
	  Phase_String = 'wat'
      Lut_File = trim(Identifier_current)//"_ch" &
            //TRIM(Chan_String)//"_ems_lut_"//Phase_String//"_cld.hdf"
    ENDIF
	
	IF (idx_phase .EQV. 2 ) THEN
	  Phase_string='ice'
	  Lut_File = trim(Identifier_current)//"_ch" &
            //TRIM(Chan_String)//"_ems_lut_"//Phase_String//"_cld.hdf"
	 !-->IF (INDEX(Sensor_Name,'NOAA').OR. INDEX(Sensor_Name,'METOP')) then
         ! IF ((INDEX(Sensor_Name,'NOAA') > 0) .OR. (INDEX(Sensor_Name,'METOP')>0)) then
	   ! Lut_File = TRIM(Lut_Path)//'AVHRR'//"_ch" &
       !     //TRIM(Chan_String)//"_ems_lut_"//Phase_String//"_cld.hdf"
	 ! ENDIF		
	ENDIF
	    
		IF (Chan_String .EQV. 'no') CYCLE
		
		
        CALL POPULATE_CLOUD_LUT_SINGLE_EMS(Lut_File,idx_phase)

  END DO Emm_loop_phase

  IF (PRESENT(geocat_flag)) THEN
    Ancil_file = Identifier_current//'_gas_coeff.hdf'  
    print*,'DCOMP ancil: ', trim(Ancil_file)
    CALL POPULATE_ANCIL(ancil_file)
    Identifier_old = Identifier_current
  END IF
    
ENDIF

END SUBROUTINE POPULATE_ALL_LUT

!====================================================================
! Subroutine Name: POPULATE_CLOUD_LUT_SINGLE(lut_file, Idx_Chn, Idx_Phase) 
!
! Function: populates cloud look-up table  
!
! Inputs:
!    Lut_file - file there LUT is stored
!    Idx_Chn - channel number  in LUT
!    Idx_Phase - phase number  in LUT
!
!====================================================================
 
SUBROUTINE POPULATE_CLOUD_LUT_SINGLE( &
      lut_file &    ! - input
    , Idx_Chn   &      ! - input
        , Idx_Phase )

IMPLICIT NONE

! ---  input
CHARACTER (LEN=*), INTENT(in) ::  Lut_file   ! - file there LUT is stored
INTEGER , INTENT(in)  ::Idx_Chn               ! - channel number  in LUT 
INTEGER , INTENT(in)  ::Idx_Phase             ! - phase number  in LUT 


   
! -  hdf stuff

INTEGER:: Istatus
 
INTEGER:: Sfstart, Sfn2Index,  Sfendacc,Sfend,Sfrdata,Sfselect
INTEGER:: Sd_Id
INTEGER:: Sds_Id
 

INTEGER, PARAMETER:: Sds_Rank_1D = 1
INTEGER, DIMENSION(Sds_Rank_1D):: Sds_Start_1D, Sds_Edge_1D, Sds_Stride_1D

INTEGER, PARAMETER:: Sds_Rank_2D = 2
INTEGER, DIMENSION(Sds_Rank_2D):: Sds_Start_2D, Sds_Edge_2D, Sds_Stride_2D, Sds_Dims_2D

INTEGER, PARAMETER:: Sds_Rank_3D = 3
INTEGER, DIMENSION(Sds_Rank_3D):: Sds_Start_3D, Sds_Edge_3D, Sds_Stride_3D, Sds_Dims_3D

INTEGER, PARAMETER:: Sds_Rank_5D = 5
INTEGER, DIMENSION(Sds_Rank_5D):: Sds_Start_5D, Sds_Edge_5D, Sds_Stride_5D, Sds_Dims_5D


!----------------------------------------------------------------------------------
! Executable Code
!----------------------------------------------------------------------------------
 
 
 

Sd_Id = Sfstart( &
                 TRIM(lut_file), &
                 DFACC_READ) 
  
  
Istatus = 0


!--- Read in values in the table
Istatus = 0

!-- read in vectors
 Sds_Start_1D = 0
 Sds_Stride_1D = 1

!-- sensor zenith angle   
Sds_Edge_1D = Num_Sat_zen
Sds_Id = Sfselect(Sd_Id,Sfn2Index(Sd_Id,"sensor_zenith_angle"))
Istatus = Sfrdata(Sds_Id,Sds_Start_1D, Sds_Stride_1D, Sds_Edge_1D, &
                  Cld_Refl_Lut%Channel(Idx_Chn)%Phase(Idx_Phase)%Sat_Zen) +  &
                  Istatus
Istatus = Sfendacc(Sds_Id) + Istatus
 
!-- solar sensor zenith angle
Sds_edge_1d = Num_sol_zen
Sds_id = sfselect(sd_id,sfn2index(sd_id,"solar_zenith_angle"))
Istatus = sfrdata(sds_id,sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  Cld_Refl_Lut%Channel(Idx_Chn)%Phase(Idx_Phase)%sol_zen) +  &
                  istatus
Istatus = sfendacc(sds_id) + istatus
 
!-- relative azimith angle
 
sds_edge_1d = Num_Rel_Azi
sds_id = sfselect(sd_id,sfn2index(sd_id,"relative_azimuth_angle"))
istatus = sfrdata(sds_id,sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  Cld_Refl_Lut%Channel(Idx_Chn)%Phase(Idx_Phase)%Rel_Azi) +  &
                  istatus
istatus = sfENDacc(sds_id) + istatus


!-- log10 optical depth
 
sds_edge_1d = NUM_COD
sds_id = sfselect(sd_id,sfn2index(sd_id,"log10_optical_depth"))
istatus = sfrdata(sds_id,sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  Cld_Refl_Lut%Channel(Idx_Chn)%Phase(Idx_Phase)%COD_Vec) +  &
                  istatus
istatus = sfENDacc(sds_id) + istatus
 

!-- log10 effective radius
 
sds_edge_1d = NUM_CPS
sds_id = sfselect(sd_id,sfn2index(sd_id,"log10_eff_radius"))
istatus = sfrdata(sds_id,sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  Cld_Refl_Lut%Channel(Idx_Chn)%Phase(Idx_Phase)%CPS_Vec) +  &
                  istatus
istatus = sfENDacc(sds_id) + istatus
 

!--- 2d arrays
sds_start_2d = 0
sds_stride_2d = 1


!--- spherical albedo
sds_dims_2d(1) = NUM_CPS   ! NUM_CPS
sds_dims_2d(2) = NUM_COD  ! NUM_COD 
sds_edge_2d(1) = NUM_CPS
sds_edge_2d(2) = NUM_COD
sds_id = sfselect(sd_id,sfn2index(sd_id,"spherical_albedo"))
istatus = sfrdata(sds_id,sds_start_2d, sds_stride_2d, sds_edge_2d, &
                           Cld_Refl_Lut%Channel(Idx_Chn)%Phase(Idx_Phase)%sph_alb) +  &
                           istatus
istatus = sfENDacc(sds_id) + istatus
 
 
!--- 3d arrays
sds_start_3d = 0
sds_stride_3d = 1
 
!--- albedo
sds_dims_3d(1) = NUM_CPS
sds_dims_3d(2) = NUM_COD
sds_dims_3d(3) = num_sol_zen
sds_edge_3d(1) = NUM_CPS
sds_edge_3d(2) = NUM_COD
sds_edge_3d(3) = num_sol_zen
sds_id = sfselect(sd_id,sfn2index(sd_id,"albedo"))
istatus = sfrdata(sds_id,sds_start_3d, sds_stride_3d, sds_edge_3d, &
                           Cld_Refl_Lut%Channel(Idx_Chn)%Phase(Idx_Phase)%Cld_Alb) +  &
                           istatus
istatus = sfENDacc(sds_id) + istatus
 

 
!-- transmission
sds_dims_3d(1) = NUM_CPS
sds_dims_3d(2) = NUM_COD
sds_dims_3d(3) = num_sol_zen
sds_edge_3d(1) = NUM_CPS
sds_edge_3d(2) = NUM_COD
sds_edge_3d(3) = num_sol_zen
sds_id = sfselect(sd_id,sfn2index(sd_id,"transmission"))
istatus = sfrdata(sds_id,sds_start_3d, sds_stride_3d, sds_edge_3d, &
                           Cld_Refl_Lut%Channel(Idx_Chn)%Phase(Idx_Phase)%trans) +  &
                           istatus
istatus = sfendacc(sds_id) + istatus
 

!--- 5d arrays
sds_start_5d = 0
sds_stride_5d = 1


!--- reflectance
sds_dims_5d(1) = NUM_CPS
sds_dims_5d(2) = NUM_COD
sds_dims_5d(3) = num_sol_zen
sds_dims_5d(4) = num_sat_zen
sds_dims_5d(5) = Num_Rel_Azi

sds_edge_5d(1) = NUM_CPS
sds_edge_5d(2) = NUM_COD
sds_edge_5d(3) = num_sol_zen
sds_edge_5d(4) = num_sat_zen
sds_edge_5d(5) = Num_Rel_Azi

sds_id = sfselect(sd_id,sfn2index(sd_id,"reflectance"))
istatus = sfrdata(sds_id,sds_start_5d, sds_stride_5d, sds_edge_5d, &
                           Cld_Refl_Lut%Channel(Idx_Chn)%Phase(Idx_Phase)%refl) +  &
                           istatus
istatus = sfENDacc(sds_id) + istatus
 
 
IF (istatus /= 0) THEN
  PRINT "(a,'Error reading sds data from cloud lut files, stopping')",EXE_PROMPT
   print*,trim(lut_file)
  STOP
END IF

!--- Close the lookup table file
istatus = sfend(sd_id) 
IF (istatus /= 0) THEN
  PRINT "(a,'Error closing cloud lut files ')",EXE_PROMPT
  print*,trim(lut_file)
END IF
   
Cld_Refl_Lut%Channel(Idx_Chn)%Phase(Idx_Phase)%Flag = FLAG_YES

END SUBROUTINE POPULATE_CLOUD_LUT_SINGLE          





!====================================================================
! Subroutine Name: POPULATE_CLOUD_LUT_SINGLE_EMS(lut_file,Idx_Phase) 
!
! Function: populate cloud look-up table 
!
! Inputs: 
!     Lut_file - file there LUT is stored
!     Idx_Phase - phase number  in LUT
!
!====================================================================
 
   subroutine populate_cloud_lut_single_ems( &
                lut_file &    ! - input
              , idx_phase )

      implicit none

      ! ---  input
      character (len=*), intent(in) ::  lut_file   ! - file there lut is stored
      integer , intent(in)  :: idx_phase             ! - phase number  in lut 

      ! -  hdf stuff

      integer:: istatus
      integer:: sfstart, sfn2index,  sfendacc,sfend,sfrdata,sfselect
      integer:: sd_id
      integer:: sds_id
      integer, parameter:: sds_rank_3d = 3
      integer, dimension(sds_rank_3d):: sds_start_3d, sds_edge_3d, sds_stride_3d, sds_dims_3d

      
      ! executable code
      
      sd_id = sfstart( &
                 trim(lut_file), &
                 dfacc_read) 
     
      istatus = 0

      !--- 3d arrays
      sds_start_3d = 0
      sds_stride_3d = 1
 
      !--- emissivity
      sds_dims_3d(1) = num_cps
      sds_dims_3d(2) = num_cod
      sds_dims_3d(3) = num_sat_zen
      sds_edge_3d(1) = num_cps
      sds_edge_3d(2) = num_cod
      sds_edge_3d(3) = num_sat_zen
      sds_id = sfselect(sd_id,sfn2index(sd_id,"cloud_emissivity"))
      istatus = sfrdata(sds_id,sds_start_3d, sds_stride_3d, sds_edge_3d, &
                           cld_ems_lut%phase(idx_phase)%emiss) +  &
                           istatus
      istatus = sfendacc(sds_id) + istatus
 
      !-- transmission
      sds_dims_3d(1) = num_cps
      sds_dims_3d(2) = num_cod
      sds_dims_3d(3) = num_sol_zen
      sds_edge_3d(1) = num_cps
      sds_edge_3d(2) = num_cod
      sds_edge_3d(3) = num_sol_zen
      sds_id = sfselect(sd_id,sfn2index(sd_id,"cloud_transmission"))
      istatus = sfrdata(sds_id,sds_start_3d, sds_stride_3d, sds_edge_3d, &
                           cld_ems_lut%phase(idx_phase)%trans) +  &
                           istatus
      istatus = sfendacc(sds_id) + istatus
  
      !--- close the lookup table file
      istatus =  sfend (sd_id) 
      if (istatus /= 0) then
         print "(a,'error reading ems lut files')",exe_prompt
      end if
   
      cld_ems_lut%phase(idx_phase)%flag = flag_yes

   end subroutine populate_cloud_lut_single_ems      

 
!====================================================================
! SUBROUTINE Name: HDF_DCOMP_READ_ANCIL
!
! Function:     Reads ancillary data from HDF file
!
! Description: 
!
! Calling Sequence:
!
!      CALL HDF_DCOMP_READ_ANCIL(ancil_file,ozone,gas)
!
! Inputs: 
!     Ancil_file:  string which contains the ancillary HDF file
!
! Output:
!
!     Ozone : 3 ozone coefficients needed for transmission calculations in VIS channel
!     Ozone : 3 x 3 ozone coefficients needed for transmission calculations in VIS(0.6)  and 1.6 and 3.9 channels
!    
!
! History:    created April 2010
!
!  Global variables: 
!
! Local Variables:
!    
!     ozone_coeff
!     gas_coeff
!
!====================================================================
!

SUBROUTINE POPULATE_ANCIL(ancil_file)

!- input 
CHARACTER (LEN=*), INTENT(in):: Ancil_file

!- locals
REAL(KIND=real4), DIMENSION(3)::ozone_coeff
REAL(KIND=real4), DIMENSION(3,3)::gas_coeff

!-- hdf
INTEGER:: Sfstart,Sfn2Index,  Sfendacc ,Sfrdata,Sfselect
INTEGER:: Sd_Id
INTEGER:: Sds_Id
INTEGER:: Istatus

INTEGER, PARAMETER:: Sds_Rank_1D = 1

INTEGER, PARAMETER:: Sds_Rank_2D = 2
INTEGER, DIMENSION(Sds_Rank_2D):: Sds_Start_2D, Sds_Edge_2D, Sds_Stride_2D

!----------------------------------------------------------------------------------
! Executable Code
!----------------------------------------------------------------------------------

Sd_Id = Sfstart( &
                                TRIM(ancil_file), &
                                DFACC_READ) 

Sds_Id = Sfselect(Sd_Id,Sfn2Index(Sd_Id,"Ozone_transmission"))
Istatus = Sfrdata(Sds_Id,0, 1, 3, ozone_coeff) 

sds_start_2d = 0
sds_stride_2d = 1

sds_edge_2d(1) = 3
sds_edge_2d(2) = 3

Sds_Id = Sfselect(Sd_Id,Sfn2Index(Sd_Id,"Gas_transmission"))
Istatus = Sfrdata(Sds_Id,sds_start_2d, sds_stride_2d, sds_edge_2d, gas_coeff)
istatus = sfENDacc(sds_id)

!--- initialize 
Ancil_Data%Gas_Trans = 0.0
Ancil_Data%Ozone_Trans = 0.0

!--- populate 
Ancil_Data%Gas_Trans  =  Gas_Coeff
Ancil_Data%Ozone_Trans  =  Ozone_Coeff
        
Ancil_Data%Flag = FLAG_YES
IF (istatus /= 0) THEN
  PRINT "(a,'Error reading sds data from ANCIL files, stopping')",EXE_PROMPT
   print*,trim(ancil_file)
  STOP
END IF
 END SUBROUTINE POPULATE_ANCIL
 
 !====================================================================
! SUBROUTINE Name: GET_LUT_DATA
!
! Function:      Interface to look-up-table data
!
! Description:    handles data access to LUT data. HDF read routines 
!
! Calling Sequence:
!
! Inputs: 
!    Idx_Chn     : Channel index 
!    Idx_Phase   : Phase index
!
! Output:
!    Sat_Zen_Vec_LUT - satellite zenith angle vector
!    Sol_Zen_Vec_LUT - solar zenith angle vector
!    Rel_Azi_Vec_LUT - relative azimuth vector
!    COD_Vec_LUT     - log10 optical depth vector
!    CPS_Vec_LUT     - log10 effective radius vector
!    Sph_Alb_2D_LUT  - spherical albedo
!    Trans_3D_LUT    - flux transmission
!    Cld_Alb_3D_LUT  - cloud albedo
!    Refl_5D_LUT - reflectance
!    Trans_Ems_3D_LUT - Transmission ems
!    Ems_3D_LUT - reflectance
!    Num_Sol_Zen_LUT - total # solar zeniths
!    Num_Sat_Zen_LUT - total # satellite zeniths
!    Num_Rel_Azi_LUT - total # relative azimuths
!    Ozone_Coeff
!    Gas_Coeff
!
! History:                   
!
!====================================================================
  
SUBROUTINE GET_LUT_DATA(Idx_Chn,Idx_Phase &   
        , Num_Sol_Zen &
        , Num_Sat_Zen &
        , Num_Rel_Azi &
        , Sat_Zen_Vec_LUT &
        , Sol_Zen_Vec_LUT &
        , Rel_Azi_Vec_LUT &
        , COD_Vec_LUT &
        , CPS_Vec_LUT &
        , Sph_Alb_2D_LUT &
        , Trans_3D_LUT &
        , Cld_Alb_3d_LUT &
        , Refl_5D_LUT &
        , Trans_EMS_3D_LUT &
        , Ems_3D_LUT &
        , Gas_Coeff &
        , Ozone_Coeff  )
        
! ---  input
INTEGER , INTENT(in) :: Idx_chn    ! channel number
INTEGER , INTENT(in) :: Idx_phase  ! phase number

! --- optional output
! - number of parameters 
INTEGER,INTENT(out), OPTIONAL:: Num_Sol_Zen
INTEGER,INTENT(out), OPTIONAL:: Num_Sat_Zen
INTEGER,INTENT(out), OPTIONAL:: Num_Rel_Azi           !- total number relative azimuths
 
! -- data 
REAL (Kind=Real4),INTENT(out),OPTIONAL ,POINTER, DIMENSION(:) ::     Sat_Zen_Vec_LUT   ! - satellite zenith angle vector 
REAL (Kind=Real4),INTENT(out),OPTIONAL ,POINTER, DIMENSION(:) ::     Sol_Zen_Vec_LUT   ! - solar zenith angle vector
REAL (Kind=Real4),INTENT(out),OPTIONAL ,POINTER, DIMENSION(:) ::     Rel_Azi_Vec_LUT   ! - relative azimuth vector
REAL (Kind=Real4),INTENT(out),OPTIONAL ,POINTER, DIMENSION(:) ::     COD_Vec_LUT       ! - log10 optical depth vector
REAL (Kind=Real4),INTENT(out),OPTIONAL ,POINTER, DIMENSION(:) ::     CPS_Vec_LUT       ! - log10 effective radius vector
REAL (Kind=Real4),INTENT(out),OPTIONAL ,POINTER, DIMENSION(:,:) ::   Sph_Alb_2D_LUT    ! - spherical albedo
REAL (Kind=Real4),INTENT(out),OPTIONAL ,POINTER, DIMENSION(:,:,:) :: Trans_3D_LUT      ! - flux transmission
REAL (Kind=Real4),INTENT(out),OPTIONAL ,POINTER, DIMENSION(:,:,:) :: Cld_Alb_3D_LUT    ! - cloud albedo
REAL (Kind=Real4),INTENT(out),OPTIONAL ,POINTER, DIMENSION(:,:,:,:,:) :: Refl_5D_LUT   ! - reflectance
REAL (Kind=Real4),INTENT(out),OPTIONAL ,POINTER, DIMENSION(:,:,:) :: Trans_Ems_3D_LUT   ! - Transmission ems
REAL (Kind=Real4),INTENT(out),OPTIONAL ,POINTER, DIMENSION(:,:,:) :: Ems_3D_LUT   ! - reflectance
  
! - data from ancil file 
REAL(KIND=real4),INTENT(out), OPTIONAL , DIMENSION(3)::Ozone_Coeff
REAL(KIND=real4),INTENT(out), OPTIONAL , DIMENSION(3)::Gas_Coeff
  

INTEGER::Idx_Chn_LUT  ! - substitutes sensor channel number with 1 for VIS or 2 for IR 
   
!-------------------------------------------------------------------------------
! Executable Code
!-------------------------------------------------------------------------------

! - set  table numbers  outdated stuff can be removed after testing
!IF (Idx_chn .EQV. VIS_CHN_NUMBER) Idx_Chn_LUT = 1
!IF (Idx_chn .EQV.  IR_CHN_NUMBER_3A) Idx_Chn_LUT = 2
!IF (Idx_chn .EQV.  IR_CHN_NUMBER_3B) Idx_Chn_LUT = 3

 Idx_chn_LUT = Idx_chn


 ! - populate optional output variables
IF (present(Num_Sol_Zen)) &
  Num_Sol_Zen = SIZE(Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Sol_Zen)
  
IF (present(Num_Sat_Zen)) &
  Num_Sat_Zen = SIZE(Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Sat_Zen)
  
IF (present(Num_Rel_Azi)) &
  Num_Rel_Azi = SIZE(Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Rel_Azi)  

IF (present(Sat_Zen_Vec_LUT)) &
  Sat_Zen_Vec_LUT => Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Sat_Zen  
  
IF (present(Sol_Zen_Vec_LUT)) &
  Sol_Zen_Vec_LUT => Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Sol_Zen 
  
IF (present(Rel_Azi_Vec_LUT)) &
  Rel_Azi_Vec_LUT => Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Rel_Azi 
  
IF (present(CPS_Vec_LUT)) &
  CPS_Vec_LUT =>         Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%CPS_Vec  
  
IF (present(COD_Vec_LUT)) &
  COD_Vec_LUT =>         Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%COD_Vec 
  
IF (present(Sph_Alb_2D_LUT)) &
  Sph_alb_2D_LUT =>  Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Sph_Alb  
  
IF (present(Trans_3D_LUT)) &
  Trans_3D_LUT =>        Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Trans 
  
IF (present(Cld_Alb_3D_LUT)) &
  Cld_Alb_3D_LUT =>  Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Cld_Alb 
  
IF (present(Refl_5D_LUT)) &
  Refl_5D_LUT =>         Cld_Refl_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Refl 
   
IF (present(Trans_ems_3D_LUT)) &
  Trans_ems_3D_LUT =>   Cld_Ems_Lut%Phase(Idx_Phase)%Trans   
 
IF (present(Ems_3d_LUT)) &
   Ems_3D_LUT => Cld_Ems_Lut%Phase(Idx_Phase)%Emiss 

! - populate optional output ancil parameters
IF (present(Gas_coeff)) &
  Gas_coeff = Ancil_Data%Gas_Trans(Idx_Chn_LUT,:)
  
IF (present(Ozone_coeff)) THEN

  select case (Idx_Chn)
   case(1) 
      Ozone_coeff(:) = (/-0.000606266,9.77984e-05,-1.67962e-08/)
   case(2) 
      Ozone_coeff(:) = (/-0.000606266,9.77984e-05,-1.67962e-08/)
   case(5) 
      Ozone_coeff(:) = (/0.000178786,8.88510e-05,3.46661e-09/)
   case(6) 
      Ozone_coeff(:) = (/0.000986278,8.23759e-05,1.07635e-08/)
   case(7) 
      Ozone_coeff(:) = (/0.000638877,8.45589e-05,9.84943e-09/)
   case(8) 
      Ozone_coeff(:) = (/-2.74525e-05,9.02744e-05,2.31578e-09/)
   case(9) 
      Ozone_coeff(:) = (/-0.000856822,9.79630e-05,-1.48645e-08/)
   case(10) 
      Ozone_coeff(:) = (/0.000663486,8.40186e-05,1.07917e-08/)
   case(11) 
      Ozone_coeff(:) = (/0.000897866,8.26629e-05,1.45399e-08/)
   case(12) 
      Ozone_coeff(:) = (/0.000210246,8.78771e-05,4.71260e-09/)
   case(13) 
      Ozone_coeff(:) = (/0.00114584,8.06184e-05,1.62876e-08/)
   case(14) 
      Ozone_coeff(:) = (/-0.000665083,9.51923e-05,-8.21226e-09/)
   case(15)
      Ozone_coeff(:) = (/0.000898769,8.11015e-05,1.69703e-08/)
   case(16) 
      Ozone_coeff(:) = (/-0.000822729,9.65522e-05,-1.17455e-08/)
   case(17)
      Ozone_coeff(:) = (/0.000913279,8.22402e-05,1.27435e-08/)
   case(18)
      Ozone_coeff(:) = (/0.000258731,8.80321e-05,4.91110e-09/)
   case default
      Ozone_coeff(:) = (/0.000258731,8.80321e-05,4.91110e-09/)
 
 end select 

 
ENDIF
END SUBROUTINE GET_LUT_DATA


END MODULE HDF_DCOMP_AWG

