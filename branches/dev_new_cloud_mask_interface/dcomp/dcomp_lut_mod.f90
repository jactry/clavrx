! $Id: dcomp_lut_mod.f90 90 2014-03-31 19:33:53Z awalther $
module dcomp_lut_mod


 
   implicit  none
 
   private
 
   ! some params 
   integer, parameter, private:: int1 = selected_int_kind(1)
   integer, parameter, private:: real4 = selected_real_kind(6,37)
   character(*), parameter, private :: exe_prompt = "dcomp_lut_mod>> "
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
      logical :: is_set                          !- flag whether populated or not
	  logical :: is_alloc
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
      logical :: is_set                                !- phase flag 
      type(reflectance_lookup_table),dimension(num_phase) :: phase  !- cloud phase
   end type phase_reflectance_table

   type :: channel_reflectance_tables1
      logical :: is_set                              !-channel flag
	  logical :: is_alloc
      type(phase_reflectance_table), dimension(:) , allocatable :: channel
   end type channel_reflectance_tables1
  
   type (channel_reflectance_tables1), save, target:: cld_refl_lut

   !-----------------------------------------------------------------------
   ! structures for emissivity tables
   !-----------------------------------------------------------------------
   type, private :: emissivity_lookup_table
      logical :: is_set
	  
      real (kind=real4), dimension(num_cps,num_cod,num_sat_zen) :: trans 
      real (kind=real4), dimension(num_cps,num_cod,num_sat_zen) :: emiss 
   end type emissivity_lookup_table

   type, private:: phase_emissivity_table
      logical :: is_set
      type(emissivity_lookup_table),dimension(num_phase):: phase
   end type phase_emissivity_table
   
   type :: channel_emissivity_tables1
      logical :: is_set                                  !-channel flag
	  logical :: is_alloc
      type(phase_emissivity_table), dimension(:) , allocatable  :: channel
   end type channel_emissivity_tables1


   Type (channel_Emissivity_Tables1), private, SAVE, TARGET:: cld_ems_LUT

   !-- ancil data structure
   type :: ancil_data_type
      integer(kind=int1)::flag
      real(kind=real4),dimension(3,3) :: gas_trans  ! -- gas transmission coefficients
      real(kind=real4),dimension(3) :: ozone_trans        ! -- ozone transmission coefficients 
   end type ancil_data_type
 
   type(ancil_data_type),save :: ancil_data
   
   character ( len = 1024 ) , save :: identifier_old
   character ( len = 1024 ) :: identifier_current
   character ( len = 20 ) :: sensor_current	  
   
   integer , dimension ( : ) , allocatable :: mapped_modis_channels 
 
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
   subroutine populate_all_lut( &
         sensor_name , channels , lut_path )  ! - input

    
      character (len=*), intent(in) :: sensor_name ! - sensor name as appeared in lut file name
      integer, intent(in), dimension(:)  :: channels
      character ( len = 1024 ) , intent ( in ),optional :: lut_path
	  
	  
	!  character (len=255) ::  lut_path  ! - path there lut is stored


      !-locals
	   character ( len = 1024 ) :: lut_path_loc
	   integer :: n_channels 
      character ( len = 3 ) , dimension(2)   :: phase_string = [ 'wat',  'ice' ]
      character ( len = 3 ) , dimension(30) :: chan_string ='no'
	  logical , dimension ( 30 ) :: has_ems_table = .false.
	  logical , dimension ( 30 ) :: has_sol_table = .false.
      character ( len = 255 ) :: lut_file     ! 
      integer :: idx_chn
      integer :: idx_phase
 
	  integer :: i_channel

      !- executable
	  n_channels = size ( channels )
	  sensor_current = trim(sensor_name)
	  
	 
	  lut_path_loc = '/DATA/Ancil_Data/clavrx_ancil_data/luts/cld/'
     
	 
	 if (  present ( lut_path ))  lut_path_loc = trim(lut_path) 
	  
	  ! mapping sensor channel emis yes/no
	  sensor_block: select case ( sensor_name )
	     case ('Meteosat-8','Meteosat-9','Meteosat-10') sensor_block
		    has_sol_table(1:2) = .true.
			has_sol_table(6) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(2) = '2'
			chan_string(6) = '3'
			chan_string(20) = '4'
        case ('NOAA-05','NOAA-06','NOAA-07','NOAA-08','NOAA-09', 'NOAA-10','NOAA-11','NOAA-12', 'NOAA-14','TIROS-N')  sensor_block
		    has_sol_table(1) = .true.
			
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			
			chan_string(20) = '3b'
        
		case ('NOAA-15','NOAA-16', 'NOAA-17','NOAA-18','NOAA-19','METOP-A','METOP-B')  sensor_block
		    has_sol_table(1) = .true.
			has_sol_table(6) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(6) = '3a'
			chan_string(20) = '3b'		
         case ('GOES-08','GOES-09','GOES-10', 'GOES-11' , 'GOES-12' , 'GOES-13',  'GOES-14', 'GOES-15','COMS-1'  )   sensor_block
		    has_sol_table(1) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(20) = '2'
         case ('MODIS-AQUA', 'MODIS-TERRA')    sensor_block
            has_sol_table(1:2) = .true.
			has_sol_table(5:7) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(2) = '2'
			chan_string(5) = '5'
			chan_string(6) = '6'
			chan_string(7) = '7'
			chan_string(20) = '20'
         
         case('GOES-16')  sensor_block
		   has_sol_table(1) = .true.
			has_sol_table(6) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '2'
			chan_string(6) = '5'
			chan_string(20) = '7'
         
         case('ABI') sensor_block
         has_sol_table(1) = .true.
			has_sol_table(6) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '2'
			chan_string(6) = '5'
			chan_string(20) = '7'
            
         case ('AATSR')   sensor_block
		    has_sol_table(1) = .true.
			has_sol_table(6) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(6) = '6'
			chan_string(20) = '20'
             
         case ('VIIRS')   sensor_block
		   has_sol_table(1) = .true.
			has_sol_table(5) = .true.
			has_sol_table(6) = .true.
			has_sol_table(7) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '5'
			chan_string(5) = '8'
			chan_string(6) = '10'
			chan_string(7) = '11'
			chan_string(20) ='12'
         
         case ('MTSAT-1R')   sensor_block
         has_sol_table(1) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(20) = '5'  
         
         case ('MTSAT-2')   sensor_block
         has_sol_table(1) = .true.
			has_sol_table(20) = .true.
			has_ems_table(20) = .true.
			chan_string(1) = '1'
			chan_string(20) = '5'  	
         
         
         
         case default
             print*,'add sensor in dcomp_lut_mod.f90 routine populate...', trim(sensor_name)
                      stop
	  
	  end select sensor_block
	  
	 
	  if ( .not. Cld_Refl_Lut%is_alloc ) allocate (  Cld_Refl_Lut%channel (n_channels) )
	  if ( .not. Cld_Ems_Lut%is_alloc) allocate (  Cld_Ems_Lut%channel (n_channels) )
	  Cld_Refl_Lut%is_alloc = .true.
	  Cld_Ems_Lut%is_alloc = .true.
	       
       identifier_current = trim(lut_path_loc) // trim ( sensor_name )
     
      
      if ( trim(sensor_name) == 'MODIS-AQUA' .or. trim(sensor_name) == 'MODIS-TERRA' ) then
         identifier_current = trim(lut_path_loc)//'MODIS'
      end if
	
	
      if (identifier_current == identifier_old) then
         print*,trim(identifier_current)//' already populated!'
         return
      end if
	  
	 identifier_old = identifier_current
	  
	  if (.not. allocated (mapped_modis_channels) ) allocate (  mapped_modis_channels ( n_channels) )
	  mapped_modis_channels  = channels
	  
	  loop_channel : do i_channel = 1 , n_channels
	     idx_chn = channels ( i_channel ) 
		 
		 if ( .not. has_sol_table ( idx_chn ) )  cycle 
	     loop_phase: do idx_phase = 1 , 2
		    lut_file = trim(identifier_current)//'_ch' &
			         //trim ( chan_string ( idx_chn ) ) &
				     //'_ref_lut_'//phase_string(idx_phase)//'_cld.hdf'
                
            call populate_cloud_lut_single ( lut_file , i_channel , idx_phase ) 
		    
			if ( has_ems_table(idx_chn) ) then
			   lut_file = trim(identifier_current)//'_ch' &
                       //trim(chan_string ( idx_chn) ) &
					   //"_ems_lut_"//phase_string(idx_phase)//'_cld.hdf'
                    
               call populate_cloud_lut_single_ems(lut_file , i_channel , idx_phase)
			end if
 			
		 end do loop_phase
	  end do loop_channel
      

end subroutine populate_all_lut

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
   
Cld_Refl_Lut % channel ( Idx_Chn ) % Phase ( Idx_Phase ) % is_set = .true.

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
			  , idx_chn &	
              , idx_phase )

      implicit none

      ! ---  input
      character (len=*), intent(in) ::  lut_file   ! - file there lut is stored
      integer , intent(in)  :: idx_phase             ! - phase number  in lut 
	  integer , intent(in) :: idx_chn

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
                           cld_ems_lut%Channel(Idx_Chn)%phase(idx_phase)%emiss) +  &
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
                           cld_ems_lut%Channel(Idx_Chn)%phase(idx_phase)%trans) +  &
                           istatus
      istatus = sfendacc(sds_id) + istatus
  
      !--- close the lookup table file
      istatus =  sfend (sd_id) 
      if (istatus /= 0) then
         print "(a,'error reading ems lut files')",exe_prompt
      end if
   
      cld_ems_lut%Channel(Idx_Chn)%phase(idx_phase)%is_set = .true.

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
  
   subroutine get_lut_data(idx_chn,idx_phase &   
        , num_sol_zen &
        , num_sat_zen &
        , num_rel_azi &
        , sat_zen_vec_lut &
        , sol_zen_vec_lut &
        , rel_azi_vec_lut &
        , cod_vec_lut &
        , cps_vec_lut &
        , sph_alb_2d_lut &
        , trans_3d_lut &
        , cld_alb_3d_lut &
        , refl_5d_lut &
        , trans_ems_3d_lut &
        , ems_3d_lut &
        , gas_coeff &
        , ozone_coeff  )
        
      ! ---  input
      integer , intent(in) :: idx_chn    ! channel number modis mapping
      integer , intent(in) :: idx_phase  ! phase number

      ! --- optional output
      ! - number of parameters 
      integer,intent(out), optional:: num_sol_zen
      integer,intent(out), optional:: num_sat_zen
      integer,intent(out), optional:: num_rel_azi           !- total number relative azimuths
 
      ! -- data 
      real (kind=real4),intent(out),optional ,pointer, dimension(:) ::     sat_zen_vec_lut   ! - satellite zenith angle vector 
      real (kind=real4),intent(out),optional ,pointer, dimension(:) ::     sol_zen_vec_lut   ! - solar zenith angle vector
      real (kind=real4),intent(out),optional ,pointer, dimension(:) ::     rel_azi_vec_lut   ! - relative azimuth vector
      real (kind=real4),intent(out),optional ,pointer, dimension(:) ::     cod_vec_lut       ! - log10 optical depth vector 
      real (kind=real4),intent(out),optional ,pointer, dimension(:) ::     cps_vec_lut       ! - log10 effective radius vector
      real (kind=real4),intent(out),optional ,pointer, dimension(:,:) ::   sph_alb_2d_lut    ! - spherical albedo
      real (kind=real4),intent(out),optional ,pointer, dimension(:,:,:) :: trans_3d_lut      ! - flux transmission
      real (kind=real4),intent(out),optional ,pointer, dimension(:,:,:) :: cld_alb_3d_lut    ! - cloud albedo
      real (kind=real4),intent(out),optional ,pointer, dimension(:,:,:,:,:) :: refl_5d_lut   ! - reflectance
      real (kind=real4),intent(out),optional ,pointer, dimension(:,:,:) :: trans_ems_3d_lut   ! - transmission ems
      real (kind=real4),intent(out),optional ,pointer, dimension(:,:,:) :: ems_3d_lut   ! - reflectance
  
      ! - data from ancil file 
      real(kind=real4),intent(out), optional , dimension(3)::ozone_coeff
      real(kind=real4),intent(out), optional , dimension(3)::gas_coeff
  

      integer::idx_chn_lut  ! - substitutes sensor channel number with 1 for vis or 2 for ir 
      integer :: i_channel
!-------------------------------------------------------------------------------
! Executable Code
!-------------------------------------------------------------------------------
    
      idx_chn_lut = -1
	  channel_loop: do i_channel = 1 , size ( mapped_modis_channels )  
	     if ( idx_chn == mapped_modis_channels(i_channel) ) then
		    idx_chn_lut = i_channel
        end if 
      end do channel_loop
	  
	  if ( idx_chn_lut == -1 ) then
	     print*, 'error: channel was not set in populate_lut wrong channel'
	     return
	  end if
	  
      if ( .not. Cld_Refl_Lut % channel ( Idx_Chn_lut ) % Phase ( Idx_Phase ) % is_set ) then
         print*, 'dcomp get_data: this channel was not populated'
		 print*, 'check your settings'
		 print*,'bad MODIS-like channel ', idx_chn , ' for sensor ' , sensor_current
		 stop
      end if 
 
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
  Trans_ems_3D_LUT =>   Cld_Ems_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Trans   
 
IF (present(Ems_3d_LUT)) &
   Ems_3D_LUT => Cld_Ems_Lut%Channel(Idx_Chn_LUT)%Phase(Idx_Phase)%Emiss 

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


END MODULE dcomp_lut_mod

