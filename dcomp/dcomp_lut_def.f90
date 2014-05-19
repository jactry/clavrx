! $Header: /home/repository/cloud_team_dcomp/dcomp_lut_mod.f90,v 1.2 2013/12/05 20:27:16 awalther Exp $
module dcomp_lut_def

   implicit  none
  
   ! some params 
   integer, parameter, private:: int1 = selected_int_kind(1)
   integer, parameter, private:: real4 = selected_real_kind(6,37)
   integer, parameter:: dfacc_read = 1
 
   integer,parameter::flag_yes = 1
   integer,parameter::flag_no = 0
 
!   integer(kind=int1),parameter :: water_phase = 1
!   integer(kind=int1),parameter :: ice_phase = 2
 
   integer(kind=int1),parameter :: num_phase = 2
 
   integer(kind=int1), parameter :: vis_chn_number = 1
   integer(kind=int1), parameter :: ir_chn_number_3a = 2      !--- avhrr
   integer(kind=int1), parameter :: ir_chn_number_3b = 3      !--- avhrr

   integer (kind=int1), parameter :: num_sat_zen = 45              !- size of zenith dimension
   integer (kind=int1), parameter :: num_sol_zen = 45              !- size of solar zenith dimension
   integer (kind=int1), parameter :: num_rel_azi = 45              !- size of relative az. dimension
   integer (kind=int1), parameter :: num_cod = 29
   integer (kind=int1), parameter :: num_cps = 9
 
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
  
   type (channel_reflectance_tables1), save, target, public:: cld_refl_lut

   !-----------------------------------------------------------------------
   ! structures for emissivity tables
   !-----------------------------------------------------------------------
   type :: emissivity_lookup_table
      logical :: is_set
	  
      real (kind=real4), dimension(num_cps,num_cod,num_sat_zen) :: trans 
      real (kind=real4), dimension(num_cps,num_cod,num_sat_zen) :: emiss 
   end type emissivity_lookup_table

   type:: phase_emissivity_table
      logical :: is_set
      type(emissivity_lookup_table),dimension(num_phase):: phase
   end type phase_emissivity_table
   
   type :: channel_emissivity_tables1
      logical :: is_set                                  !-channel flag
	  logical :: is_alloc
      type(phase_emissivity_table), dimension(:) , allocatable  :: channel
   end type channel_emissivity_tables1


   Type (channel_Emissivity_Tables1), SAVE, TARGET, public:: cld_ems_LUT

   !-- ancil data structure
   type :: ancil_data_type
      integer(kind=int1)::flag
      real(kind=real4),dimension(3,3) :: gas_trans  ! -- gas transmission coefficients
      real(kind=real4),dimension(3) :: ozone_trans        ! -- ozone transmission coefficients 
   end type ancil_data_type
 
   type(ancil_data_type), save, public :: ancil_data
   
   character ( len = 1024 ) , save :: identifier_old
   character ( len = 1024 ) :: identifier_current
   character ( len = 20 ) :: sensor_current	  
   
   integer , dimension ( : ) , allocatable, public :: mapped_modis_channels 
 
contains
 


END MODULE dcomp_lut_def

