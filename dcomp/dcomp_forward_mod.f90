! $Id: dcomp_forward_mod.f90 77 2014-02-13 18:55:48Z awalther $
module dcomp_forward_mod
   use dcomp_lut_mod
   use dcomp_tools

   integer, parameter, public:: int1 = selected_int_kind(1)
   integer, parameter, public:: int2 = selected_int_kind(3)
   integer, parameter, public:: int4 = selected_int_kind(8)
   integer, parameter, public:: int8 = selected_int_kind(10)
   integer, parameter, public:: real4 = selected_real_kind(6,37)
   integer, parameter, public:: real8 = selected_real_kind(15,307)
   integer, parameter, public:: ipre = real4
   real , parameter :: PI = 3.14159265
   real , parameter :: DTOR = PI / 180.
   
   type pixel_vec
      real :: sol_zen
	   real :: sat_zen
	   real :: rel_azi
	   real :: ctt
      logical :: is_water_phase
   end type pixel_vec
   
   character ( len = 255 ) :: dcomp_ancil_path
   
   public :: dcomp_forward_computation
   public :: compute_cld_albedo_vis
   integer , save :: sol_pos, sol2_pos , sat_pos , azi_pos
   real , save :: planck_rad
   real , save :: rad_to_refl 
   
contains
   
   subroutine dcomp_forward_computation ( &
                                    state_vec  &   ! -input
                                  , pixel &
                                  , sensor &
                                  , channels &
                                  , alb_sfc &
                                  , air_trans_ac &
                                  , fm_vec &         ! - output
                                  , cld_trns_sol &
                                  , cld_trns_sat &
                                  , cld_sph_alb &
                                  , kernel & 
                                  , rad_abv_cld &
                                  , rad_clear_toc  &
                                  , lut_path &
                                  , cld_albedo_vis  )
                                    
   
      implicit none
      
      real , intent(in) :: state_vec (2)
      type ( pixel_vec ) , intent ( in ) :: pixel
      integer , intent ( in ) :: channels ( :)
      real,  intent ( in ) :: air_trans_ac ( : )   
      character ( len = * ) , intent (in ) :: sensor     
      real, intent(in) ::  alb_sfc (2)  
      real, optional, intent ( in ) ::   rad_abv_cld 
      real, optional, intent ( in ) ::   rad_clear_toc 
      character ( len = 1024 ) , intent ( in ) , optional :: lut_path  
         
      real , intent (out) :: fm_vec ( 2 )
      real , intent ( out) :: kernel( 2 , 2 )
      real , intent ( out ) :: cld_trns_sol( 2 ) 
      real , intent ( out ) :: cld_trns_sat( 2 )
      real , intent ( out ) :: cld_sph_alb( 2 )
      real , intent ( out ) , optional :: cld_albedo_vis
                 
      real , pointer , save :: CPS_Vec_LUT (:)  !- CPS vector in LUT
      real , pointer , save :: COD_Vec_LUT (:)  !- COD vector in LUT
         
      integer , save :: cps_count , cod_count
      integer , save :: sol_zen_count, sat_zen_count, rel_azi_count
      integer , save :: cps_pos, cod_pos
      
      integer :: sol_pos_tmp, sol2_pos_tmp , sat_pos_tmp , azi_pos_tmp
	   real, save :: cps_wgt, cod_wgt
      real, save :: sol_wgt, sol2_wgt, sat_wgt, azi_wgt
         
      real, save :: sol_zen_last, sat_zen_last, rel_azi_last
       real ::  get_rad_refl_factor
        
      real , pointer, save :: Sol_Zen_Vec_LUT(:) !- solar zenith angle vector   
      real , pointer, save :: Sat_Zen_Vec_LUT(:) !- solar zenith angle vector
      real , pointer, save :: Rel_Azi_Vec_LUT(:) !- relative azimuth vector
         
      integer :: phase_num
      real :: fm_nir_terr 
      real :: kernel_nir_terr_cod
      real :: kernel_nir_terr_cps
         
      
      real :: get_planck_radiance_39um
      real, dimension(2) :: trans_two_way
      real :: air_mass_two_way
      real :: air_mass_sat
      real :: trans_abv_cld
         
	   logical , save :: read_angle_val = .false.

      integer :: i_channel , n_channels
      character ( len = 20 ) , save :: sensor_set 
      
      integer :: MODIS_EQUIVALENT_39_CHN = 20
	  
	   ! - executable
	   n_channels = size ( channels )
      ! - populate lut if needed 
	   if ( trim(sensor) /= trim(sensor_set) ) call populate_all_lut ( sensor  , channels , lut_path = lut_path)
	   sensor_set = trim(sensor)
      phase_num = 2
      if ( pixel % is_water_phase ) phase_num = 1
	  
	   if ( .not. read_angle_val ) then
         call Get_Lut_Data ( channels(1) , phase_num &
            & , num_sol_zen = sol_zen_count &
            & , num_sat_zen = sat_zen_count &
            & , num_rel_azi = rel_azi_count &
            & , Sol_Zen_Vec_LUT = Sol_Zen_Vec_LUT &
            & , Sat_Zen_Vec_LUT = Sat_Zen_Vec_LUT &
            & , Rel_Azi_Vec_LUT = Rel_Azi_Vec_LUT )   
         call Get_lut_data ( channels(1), phase_num  &
            & , CPS_Vec_LUT  = CPS_Vec_LUT &
            & , COD_Vec_LUT  = COD_Vec_LUT )
			
		   cps_count = size(CPS_Vec_LUT)
         cod_count = size(COD_Vec_LUT)	
	      read_angle_val = .true.
		 
	   end if
	       
      ! - find angle pos and weight if needed
      if ( pixel % sol_zen /= sol_zen_last &
             &  .or. pixel % sat_zen /= sat_zen_last &
             &  .or. pixel % rel_azi /= rel_azi_last ) then  
         
         call dcomp_Interpolation_Weight ( Sol_Zen_count ,  pixel % Sol_Zen  &
                                     & , Sol_Zen_Vec_LUT , Sol_Wgt , Sol_Pos_tmp , sol_pos)         
         call dcomp_Interpolation_Weight ( Sol_Zen_count ,  pixel % Sat_Zen &
                                     & , Sol_Zen_Vec_LUT , Sol2_Wgt , Sol2_Pos_tmp , sol2_pos )  
         call dcomp_Interpolation_Weight ( Sat_Zen_count ,  pixel % Sat_Zen &
                                     & , Sat_Zen_Vec_LUT , Sat_Wgt , Sat_Pos_tmp , sat_pos)
         call dcomp_Interpolation_Weight ( Rel_Azi_count ,  pixel % Rel_Azi  &
                                     & , Rel_Azi_Vec_LUT , Azi_Wgt , Azi_Pos_tmp  ,azi_pos)
         sol_zen_last = pixel % sol_zen
         sat_zen_last = pixel % sat_zen
         rel_azi_last = pixel % rel_azi
      end if       
         
             
      call Dcomp_Interpolation_Weight (cps_count , state_vec(2) &
                                  ,CPS_Vec_LUT , CPS_Wgt , CPS_Pos)   
								                           
      call Dcomp_Interpolation_Weight (cod_count , state_vec(1) &
                                 ,COD_Vec_LUT , COD_Wgt , COD_Pos )
            
      if ( cps_pos >= ubound ( cps_vec_lut , dim = 1 )) then
         cps_pos = ubound ( cps_vec_lut , dim = 1 ) - 1
         cps_wgt = 1.
      end if
          
      if ( cod_pos >= ubound ( cod_vec_lut , dim = 1 ) ) then
         cod_pos = ubound ( cod_vec_lut , dim = 1) - 1
      end if
                                                                         
      air_mass_two_way = ( 1. / cos (pixel % Sol_zen * DTOR ) + 1./ cos ( pixel % Sat_zen * DTOR ) )
	  
	   do i_channel = 1 , n_channels
        
		   call compute_fm_single_channel (  channels ( i_channel )   , phase_num  &
                                     & , cps_wgt , cod_wgt  &
                                     & , cps_pos , cod_pos, sol_pos , sat_pos , azi_pos & 
                                     & , alb_sfc( i_channel ) &
                                     & , fm_vec(i_channel) , kernel(i_channel,1), kernel(i_channel,2) &
                                     & , cld_trns_sol(i_channel), cld_trns_sat(i_channel) , cld_sph_alb(i_channel) )
         
         trans_two_way = air_trans_ac ** air_mass_two_way  
         fm_vec(i_channel) = fm_vec(i_channel) * trans_two_way(i_channel)
	      kernel(i_channel,1) = kernel(i_channel,1) * trans_two_way(i_channel)
         kernel(i_channel,2) = kernel(i_channel,2) * trans_two_way(i_channel)
		 
		 ! - the only emis channel is channel 20
		   if ( channels ( i_channel ) == MODIS_EQUIVALENT_39_CHN ) then
		      planck_rad = get_planck_radiance_39um ( pixel % ctt  , trim(sensor) )
			
            rad_to_refl = get_rad_refl_factor ( trim ( sensor ) , pixel % sol_zen ) 
			   air_mass_sat =  1./ cos ( pixel % Sat_zen * dtor) 
			   trans_abv_cld = air_trans_ac(i_channel) ** air_mass_sat
            
            call compute_fm_terrestrial ( channels ( i_channel) , phase_num  &
                                     & , cps_wgt , cod_wgt   &
                                     & , cps_pos , cod_pos, sat_pos  & 
                                     & , planck_rad , rad_abv_cld , trans_abv_cld , rad_clear_toc , rad_to_refl &
                                     & , fm_nir_terr , kernel_nir_terr_cod ,  kernel_nir_terr_cps )
                                      
            fm_vec(i_channel) = fm_vec(i_channel) + fm_nir_terr
            kernel ( i_channel, 1) = kernel ( i_channel , 1)   + kernel_nir_terr_cod 
            kernel ( i_channel, 2) = kernel ( i_channel , 2 )  + kernel_nir_terr_cps 
  
         end if
	   end do
      
      
      if ( present ( cld_albedo_vis )  ) then
         cld_albedo_vis =  compute_cld_albedo_vis ( phase_num &
                                     & , cps_wgt , cod_wgt  &
                                     & , cps_pos , cod_pos, sol_pos )        
      end if 
	  	  
        
   end subroutine dcomp_forward_computation
   ! ------
   !
   !  computes forward model for one channel
   !   input is the atmospheric properties (cod and cps) and observation geometry
   !   output is the cloud reflectance and other cloud properties
   !   
   ! -----
   subroutine compute_fm_single_channel  ( chn_idx  , phase_idx  &
               & , cps_wgt , cod_wgt  &
               & , cps_pos , cod_pos, sol_pos , sat_pos , azi_pos & 
               & , alb_sfc &
               & , refl , drefl_dcod, drefl_dcps &
               & , trns_sol, trns_sat , sph_alb )
   
      implicit none
      ! -- input
      integer, intent(in) :: chn_idx             ! - channel index
      integer, intent(in) :: phase_idx           ! - phase index
      integer, intent(in) :: cps_pos     ! - position  in lut
      integer, intent(in) :: cod_pos     ! - position  in lut
      integer, intent(in) :: sol_pos, sat_pos, azi_pos
        
      real, intent(in) :: cps_wgt        ! - weight beween positionn and next index in lut
      real, intent(in) :: cod_wgt        ! - weight beween positionn and next index in lut
      real, intent(in) :: alb_sfc                ! - surface albedo; value between 0 and 1
   
      ! -- output
      real, intent(out):: refl                         ! - reflectance at toa under consideration of surface impact
      real, intent(out):: drefl_dcod           ! - derivation of reflectance to delta cloud optical depth
      real, intent(out):: drefl_dcps           ! - derivation of reflectance to delta cloud particle size

      real, intent(out):: trns_sol             ! - cloud transmission for radiance coming from solar angle
      real, intent(out):: trns_sat                 ! - cloud transmission for radiance coming from satellite angle
      real, intent(out):: sph_alb                   ! - spherical albedo
         
      ! -- local
      real:: dtrans_sol_dcod   ! - derivation transmission solar angle to cod
      real:: dtrans_sol_dcps   ! - derivation transmission solar angle to cps
      real:: dtrans_sat_dcod   ! - derivation transmission satellite angle to cod
      real:: dtrans_sat_dcps   ! - derivation transmission satellite angle to cps
      real:: dsph_alb_dcod     ! - derivation to cod
      real:: dsph_alb_dcps     ! - derivation to cps
      real:: alb_sfc_term      ! - term which accounts for albedo surface 
      real:: cld_refl                  ! - cloud reflectance after interpolation

      ! - data which comes from lut
      real, dimension(:,:), pointer::  sph_alb_2d_lut
      real, dimension(:,:,:), pointer:: trans_3d_lut
      real, dimension(:,:,:,:,:), pointer:: refl_5d_lut
	  
	   real , dimension ( 2 , 2 ) ::  cld_refl_2d , trans_sol_2d , trans_sat_2d , sph_alb_2d
      real :: cod_diff,ref_diff

      !-------------------------------------------------------------------------------
      !   executable
      !-------------------------------------------------------------------------------

      call Get_Lut_Data (chn_idx, phase_idx &
                  , Refl_5D_LUT    = REFL_5D_LUT &
                  , Trans_3D_LUT   = Trans_3D_LUT &
                  , Sph_Alb_2D_LUT = Sph_Alb_2D_LUT)
           
	   cld_refl_2d = REFL_5D_LUT ( CPS_Pos : CPS_Pos + 1 &
                                , COD_Pos : COD_Pos + 1 &
                                , Sol_Pos  &
                                , Sat_Pos  &
                                , Azi_Pos  )
																		
	   trans_sol_2d = Trans_3D_LUT (  CPS_Pos : CPS_Pos+1 &
                                   , COD_Pos : COD_Pos+1 &
                                   , Sol_Pos  )
	  	  
	   trans_sat_2d = Trans_3D_LUT( CPS_Pos : CPS_Pos + 1 &
                                        , COD_Pos : COD_Pos + 1 &
                                        , Sat_Pos )
										
	   sph_alb_2d = Sph_Alb_2D_LUT(CPS_Pos:CPS_Pos+1&
                                        ,COD_Pos:COD_Pos+1)
										
	  
	   ref_diff = 0.2
	   cod_diff = 0.1
	   call interpolate_2d( cld_refl_2d , cps_wgt , cod_wgt , ref_diff , cod_diff , cld_refl , dRefl_dcps      , dRefl_dcod ) 
	   call interpolate_2d( trans_sol_2d, cps_wgt , cod_wgt , ref_diff , cod_diff , trns_sol , dtrans_sol_dcps , dtrans_sol_dcod)
	   call interpolate_2d( trans_sat_2d, cps_wgt , cod_wgt , ref_diff , cod_diff , trns_sat , dTrans_sat_dcod , dTrans_sat_dcps )
	   call interpolate_2d( sph_alb_2d  , cps_wgt , cod_wgt , ref_diff , cod_diff , sph_alb  , dsph_alb_dcod   , dSph_alb_dcps)
	  
	   trns_sol = max ( min ( trns_sol , 1. ) , 0. )
	   trns_sat = max ( min ( trns_sat , 1. ) , 0. )
	   sph_alb = max ( min ( sph_alb , 1. ) , 0. )
	   cld_refl = max ( cld_refl , 0. )
	  
      ! - Account for surface reflection and atmospheric transmission

      ! - The next two terms are described in ATBD eq. 23 
      Alb_Sfc_Term = MAX(0. , Alb_Sfc / (1.0 - Alb_Sfc * sph_alb)) 

      Refl = Cld_Refl + Alb_Sfc_Term * Trns_sol * Trns_sat 
      
      ! - derivations of equation above
      Drefl_Dcod = Drefl_Dcod  &
                 & + Alb_Sfc_Term * Trns_sol  * Dtrans_sat_Dcod &
                 & + Alb_Sfc_Term * Trns_sat  * Dtrans_sol_Dcod &
                 & +  ((Trns_sol * Trns_sat  * Alb_Sfc * Alb_Sfc * Dsph_alb_Dcod) &
                    /(( 1 - Alb_Sfc * sph_alb)**2)) 
                                     
      Drefl_Dcps = Drefl_Dcps  &
                  &  + Alb_Sfc_Term * Trns_Sol  * Dtrans_sat_Dcps &
                  &  + Alb_Sfc_Term * Trns_Sat  * Dtrans_sol_Dcps  &
                  &  +  ((Trns_sol * Trns_Sat  * Alb_Sfc * Alb_Sfc * Dsph_alb_Dcps) &
                                         /(( 1 - Alb_Sfc * Sph_alb)**2 ))
       
         
      NULLIFY (Sph_Alb_2D_Lut)
      NULLIFY (Trans_3D_Lut) 
      NULLIFY (Refl_5D_Lut)  
         
   end subroutine compute_fm_single_channel
   
   ! ---------------------------------------------------
   !     computes the terrestrial amount
   ! --------------------------
   subroutine compute_fm_terrestrial (&
                         chn_idx , phase_idx  &
                         & , cps_wgt , cod_wgt  &
                         & , cps_pos , cod_pos , sat_pos &
                         & , planck_rad_loc ,  rad_clear_ct , trans_abv_cld , rad_clear_toa , rad_to_refl &
                         & , refl_terr , drefl_dcod_terr,drefl_dcps_terr)
                                                 
       
      implicit none
                                                    
      !---input
      integer, intent (in) :: chn_idx             ! - Channel index
      integer, intent (in) :: phase_idx           ! - Phase index
      integer, intent (in) :: cps_pos     ! - Position  in LUT
      integer, intent (in) :: cod_pos     ! - Position  in LUT
      integer, intent (in) :: sat_pos
      
	  
      real, intent (in) :: cps_wgt        ! - Weight beween positionn and next index in LUT
      real, intent (in) :: cod_wgt        ! - Weight beween positionn and next index in LUT
      real, intent (in) :: planck_rad_loc
      real, intent (in) :: rad_clear_ct , trans_abv_cld , rad_clear_toa , rad_to_refl
         
      real, intent(out):: Refl_Terr
      real, intent(out):: Drefl_Dcod_Terr
      real, intent(out):: Drefl_Dcps_Terr

      real, dimension(:,:,:), POINTER:: Trans_3D_Lut
      real, dimension(:,:,:), POINTER:: EMS_3D_Lut

      ! - local sub arrays
     

      real, dimension(2,2):: Trans_sat_2d   ! - Transmission sat 2D with dim CPS and COD array before interpolation
      real, dimension(2,2):: Ems_2d   

      !- clavrx emissivty terms
      real(KIND=real4):: Rad
      real(KIND=real4):: DRad_Dcod
      real(KIND=real4):: DRad_Dcps  

      real(KIND=real4):: Cld_Emiss
      real(KIND=real4):: Dems_Dcod
      real(KIND=real4):: Dems_Dcps
      real(KIND=real4):: Cld_Trans
      real(KIND=real4):: Dtrans_Dcod
      real(KIND=real4):: Dtrans_Dcps
      real :: cod_diff,ref_diff
      !----------------------------------------------------------------------
      ! construct sub-tables from appropriate lookup tables
      !----------------------------------------------------------------------

      call get_lut_data ( Chn_idx , Phase_idx &
                        , EMS_3D_LUT = EMS_3D_LUT &   
                        , TRANS_ems_3d_LUT = TRANS_3D_LUT)
    
	   ems_2d = ems_3d_LUT(CPS_Pos:CPS_Pos+1 &
                                        , COD_Pos:COD_Pos+1 &
                                        , Sat_Pos  )
                            								
      trans_sat_2d = Trans_3D_LUT(CPS_Pos:CPS_Pos+1 &
                                        ,COD_Pos:COD_Pos+1 &
                                        , Sat_Pos )
	 	
		
	   ref_diff = 0.2
	   cod_diff = 0.1
	   call interpolate_2d( ems_2d , cps_wgt , cod_wgt , ref_diff , cod_diff , cld_emiss , dems_dcps      , dems_dcod ) 
	   call interpolate_2d( trans_sat_2d, cps_wgt , cod_wgt , ref_diff , cod_diff , cld_trans , dTrans_dcps , dTrans_dcod )
	 							     
      rad = rad_clear_ct * Cld_Emiss &
                  + trans_abv_cld * cld_emiss * planck_rad_loc &
				  + cld_trans * rad_clear_toa                 
                
      dRad_dcod =  rad_clear_ct * Dems_Dcod & 
               + trans_abv_cld * Planck_Rad_loc * Dems_Dcod  &
               + rad_clear_toa * Dtrans_Dcod

      dRad_dcps =   rad_clear_ct * Dems_Dcps & 
                + trans_abv_cld * Planck_Rad_loc * Dems_Dcps  &
                + rad_clear_toa  *Dtrans_Dcps        
                            
                            
      refl_terr = rad * rad_to_refl
         
      drefl_dcod_terr = dRad_dcod * rad_to_refl
      drefl_dcps_terr = dRad_dcps *  rad_to_refl
 
   end subroutine compute_fm_terrestrial
   
  
   
   
   !
   !
   !-----------------------------------------------------------------------------------------
   function planck_rad_fast_20 ( tmp , db_dt ) result ( rad )
      implicit none
      real, intent ( in) :: tmp
      real, intent ( out), optional :: db_dt
      real :: rad
         
      real :: db_dt_tmp
      integer, parameter:: nplanck = 161
      real, dimension(nplanck) :: B20
      real, parameter :: T_planck_min = 180.0  
      real, parameter :: delta_T_planck = 1.0
      real, dimension(nplanck) :: T_planck
      integer :: l , i
      
      real :: nu_20 , a1_20 , a2_20
      real (kind=real4), parameter :: c1 = 1.191062e-5
      real (kind=real4), parameter :: c2 = 1.4387863
         
      nu_20 = 2562.1383
      a1_20 = -1.661627
      a2_20 = 1.0023207
           
      do i = 1 , nplanck 
         T_planck(i) = T_planck_min + (i-1)*delta_T_planck
         B20(i) = c1*(nu_20**3)/(exp((c2*nu_20)/ &
              ((T_planck(i)-a1_20)/a2_20))-1.0) 
      end do                
     
      l = (tmp - T_planck_min)/delta_T_planck
      l = max(1,min(nplanck-1,l))
         
      dB_dT_tmp = (B20(l+1)-B20(l))/(T_planck(l+1)-T_planck(l))
         
      rad = B20(l) + (tmp - T_planck(l)) * (dB_dT_tmp)
         
      if ( present ( db_dt)) db_dt = db_dt_tmp
    

   end function planck_rad_fast_20

   ! ---------------------------------------------------------------------------------
   !
   !  INTERPOLATE_2D
   !
   !  Linear Interpolation for a 2 x 2 
   ! 
   !  Returns interpolated value of a 2D array with 2 elements for each dimension
   !
   !  INPUT: 
   !     table:      3D array
   !     Wgt_Dim1,Wgt_Dim2 : weights for each dimension
   !
   ! ----------------------------------------------------------------------------------
   subroutine interpolate_2d(table,Wgt_Dim1, Wgt_Dim2 , delta_1 , delta_2, value , dval_d1, dval_d2 ) 
      real, dimension(2,2), intent(in) :: table
      real, intent(in) :: Wgt_Dim1 , Wgt_Dim2
	   real, intent(in) :: delta_1 , delta_2
      
	   real , intent ( out )  :: value
	   real, intent ( out ) :: dval_d1, dval_d2
     

      !- locals
     
      real :: r , s
	  
	   r = wgt_dim1
	   s = wgt_dim2 


      value  =     ( 1.0 - r ) * ( 1.0 - s ) * table( 1 , 1 ) &
               & + r * ( 1.0 - s )          * table( 2 , 1 ) &
               & + (1.0 - r ) * s           * table( 1 , 2 ) &
				   & + r * s                    * table( 2 , 2 )

      dval_d2 =  ( ( 1.0 - r ) * ( table( 1 , 2 ) - table( 1 , 1 ) )   &
               & +  r * ( table( 2 , 2 ) - table( 2 , 1 ) ) ) /  &
               &  delta_2

      dval_d1 =   ( ( 1.0 - s ) * ( table( 2 , 1 ) - table( 1 , 1 ) )  &
                 &    + s * ( table( 2 , 2 ) - table( 1 , 2 ) ) ) /  &
                 &       delta_1 


                  
   end subroutine interpolate_2d  

   ! ------
   !
   ! -----
   function  compute_cld_albedo_vis ( phase_idx &
                                     & , cps_wgt , cod_wgt  &
                                     & , cps_pos , cod_pos, sol_pos   ) result (cld_alb)
   
      implicit none
      ! -- input
      
      integer, intent(in) :: phase_idx           ! - phase index
      integer, intent(in) :: cps_pos     ! - position  in lut
      integer, intent(in) :: cod_pos     ! - position  in lut
      integer, intent(in) :: sol_pos
        
      real, intent(in) :: cps_wgt        ! - weight beween positionn and next index in lut
      real, intent(in) :: cod_wgt        ! - weight beween positionn and next index in lut
     
   
      ! -- output
      real :: cld_alb                         ! - reflectance at toa under consideration of surface impact
      real:: dumm_1 , dumm_2           ! - derivation of reflectance to delta cloud optical depth
     
     

      ! - data which comes from lut
      
      real, dimension(:,:,:), pointer:: cld_alb_3d_lut  
	   real , dimension ( 2 , 2 ) :: cld_alb_2d
      real :: cod_diff,ref_diff

      !-------------------------------------------------------------------------------
      !   executable
      !-------------------------------------------------------------------------------

      call Get_Lut_Data ( 1 , phase_idx &
                  , cld_alb_3d_lut = cld_alb_3d_lut )
           																
	   cld_alb_2d =  cld_alb_3d_lut  (  CPS_Pos : CPS_Pos+1 &
                                     , COD_Pos : COD_Pos+1 &
                                     , Sol_Pos  )
	  	  
		  
	   ref_diff = 0.2
	   cod_diff = 0.1
	   call interpolate_2d( cld_alb_2d , cps_wgt , cod_wgt , ref_diff , cod_diff , cld_alb , dumm_1 , dumm_2 ) 
	 
   end function  compute_cld_albedo_vis
   !------------------------------------------------------------------------------
   !
   !------------------------------------------------------------------------------
   function thick_cloud_cps ( obs_vec, channels, pixel , dcomp_mode ) result ( cps)
      implicit none
   
      real, dimension(2), intent(in) :: obs_vec
      type ( pixel_vec) , intent(in) :: pixel
      integer , dimension ( : ) , intent ( in ) :: channels
      integer , intent ( in ) :: dcomp_mode
      real :: cps
      ! -- local
      integer:: COD_Pos          ! position in LUT
      integer:: CPS_Pos          ! position in LUT
      REAL(KIND=Real4):: CPS_Wgt ! position weights in LUT
      REAL (Kind=Real4), dimension(:),POINTER :: CPS_Vec_LUT  !- CPS vector in LUT
      REAL (Kind=Real4), dimension(:),POINTER :: COD_Vec_LUT  !- COD vector in LUT

      ! - local sub arrays
      REAL, dimension(9):: Refl_LUT_Sub  ! - sub table
      REAL, dimension(9):: rad  ! - sub table
      REAL, dimension(9):: refl_terr  ! - sub table
      REAL, dimension(9):: Trans_zen_LUT_Sub

      ! - data which comes from LUT
      REAL, dimension(:,:), POINTER::  Sph_Alb_2D_Lut
      REAL, dimension(:,:,:), POINTER:: Trans_3D_Lut
      REAL(KIND=real4), dimension(:,:,:,:,:), POINTER:: Refl_5D_Lut
      REAL, dimension(:,:,:), POINTER:: EMS_3D_Lut

      ! - local sub arrays

      REAL, dimension(9):: Ems_LUT_Sub  ! - sub table
    
      REAL, dimension(9):: Refl 						! - Reflectance at TOA under consideration of surface impact
      integer::ii
      integer :: phase_num

      
      cod_pos = 28
      phase_num = 2
      if ( pixel % is_water_phase ) phase_num = 1
      
      call GET_LUT_DATA(channels(2), phase_num &
        , CPS_Vec_LUT  = CPS_Vec_LUT &
        , COD_Vec_LUT  = COD_Vec_LUT &
        , EMS_3D_LUT = EMS_3D_LUT &   
        ,  TRANS_ems_3d_LUT = TRANS_3D_LUT &
         , Refl_5D_LUT    = REFL_5D_LUT)
   
      ! - build sub tables     

      Ems_LUT_Sub   = EMS_3D_LUT(:,COD_Pos,Sat_Pos )
            
      Refl_LUT_Sub      = REFL_5D_LUT(: &
                                ,COD_Pos &
                                ,Sol_Pos &
                                ,Sat_Pos &
                                ,Azi_Pos)
		Trans_zen_LUT_Sub = Trans_3D_LUT(: ,COD_Pos,Sat_Pos )

      !-- toc   
     
      Refl = refl_lut_sub  
      if ( dcomp_mode .eq. 3 ) then      
         rad =  Ems_LUT_sub * Planck_Rad 
				
         refl_terr = rad * Rad_To_Refl

         Refl = refl_lut_sub +  refl_terr
      end if
      
      cps_pos = -999
      cps = -999.
 
      do ii = 1, 8
         if ( refl(ii)-Obs_Vec(2) > 0 .and. refl(ii+1)-Obs_Vec(2) < 0 ) cps_pos = ii
      end do

      if (cps_pos .NE. -999) then
         cps_wgt = (refl(cps_pos) - Obs_Vec(2))/(refl(cps_pos) - refl(cps_pos+1))
         cps=CPS_Vec_LUT(cps_pos) + (cps_wgt * 0.2)
      end if
 
      REFL_5D_LUT=>NULL()
      Sph_Alb_2D_Lut =>NULL()
      Trans_3D_Lut =>NULL()
      return					 
      
      
   end function thick_cloud_cps
 
       
end module dcomp_forward_mod
