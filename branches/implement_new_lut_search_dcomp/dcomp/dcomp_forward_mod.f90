! $Id: dcomp_forward_mod.f90 77 2014-02-13 18:55:48Z awalther $
module dcomp_forward_mod

   use dcomp_data_pool_mod , only : &
       lut_obj , lut_output
   use dcomp_tools, only:

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
      real :: trans_two_way
      real :: air_mass_two_way
      real :: air_mass_sat
      real :: trans_abv_cld
         
	    type ( lut_output) :: lut_data

      integer :: i_channel , n_channels
      character ( len = 20 ) , save :: sensor_set 
      
      integer :: MODIS_EQUIVALENT_39_CHN = 20
	   real :: Alb_Sfc_Term
	   ! - executable
      
      call lut_obj % initialize ( sensor )
      call lut_obj % set_angles ( pixel % sat_zen , pixel % sol_zen, pixel % rel_azi )
      
	   n_channels = size ( channels )
      
	   
      phase_num = 2
      if ( pixel % is_water_phase ) phase_num = 1
	  
	
                                                                         
      air_mass_two_way = ( 1. / cos (pixel % Sol_zen * DTOR ) + 1./ cos ( pixel % Sat_zen * DTOR ) )
	  
	   do i_channel = 1 , n_channels
        
		
         
         call lut_obj % get_data ( channels ( i_channel )  , phase_num , state_vec(1), state_vec(2) , lut_data)
         
         ! -
         Alb_Sfc_Term = MAX(0. , Alb_Sfc( i_channel ) / (1.0 - Alb_Sfc( i_channel ) * lut_data % albsph )) 
         fm_vec(i_channel) = lut_data% Refl + Alb_Sfc_Term * lut_data%Trn_sol * lut_data%Trn_sat 
         trans_two_way = air_trans_ac( i_channel ) ** air_mass_two_way  
         
        
         fm_vec(i_channel) = fm_vec(i_channel) * trans_two_way
         
         
         kernel(i_channel,1) = lut_data%dRefl_dcod &
                     & +  Alb_Sfc_Term * lut_data%Trn_sol *  lut_data % dTrans_sat_dcod &
                     & + Alb_Sfc_Term * lut_data%Trn_sat  * lut_data%Dtrans_sol_Dcod &
                     & +  ((lut_data%Trn_sol * lut_data%Trn_sat  * Alb_Sfc( i_channel ) * Alb_Sfc( i_channel ) * lut_data%Dsph_alb_Dcod) &
                    /(( 1 - Alb_Sfc( i_channel ) * lut_data%albsph)**2))  
                    
                             
         kernel(i_channel,2) = lut_data%dRefl_dcps &
                     & +  Alb_Sfc_Term * lut_data%Trn_sol *  lut_data % dTrans_sat_dcps &
                     & + Alb_Sfc_Term * lut_data%Trn_sat  * lut_data%Dtrans_sol_Dcps &
                     & +  ((lut_data%Trn_sol * lut_data%Trn_sat  * Alb_Sfc( i_channel ) * Alb_Sfc( i_channel ) * lut_data%Dsph_alb_Dcps) &
                    /(( 1 - Alb_Sfc( i_channel ) * lut_data%albsph)**2))  
         
	    kernel(i_channel,1) = kernel(i_channel,1) * trans_two_way
       kernel(i_channel,2) = kernel(i_channel,2) * trans_two_way
        
		 
		 ! - the only emis channel is channel 20
		   if ( channels ( i_channel ) == MODIS_EQUIVALENT_39_CHN ) then
		      planck_rad = get_planck_radiance_39um ( pixel % ctt  , trim(sensor) )
			
            rad_to_refl = get_rad_refl_factor ( trim ( sensor ) , pixel % sol_zen ) 
			   air_mass_sat =  1./ cos ( pixel % Sat_zen * dtor) 
			   trans_abv_cld = air_trans_ac(i_channel) ** air_mass_sat
            
                               
             
            fm_nir_terr =  rad_to_refl * (rad_abv_cld * lut_data % ems &
                  + trans_abv_cld * lut_data %  ems * planck_rad &
				  + lut_data % trn_ems * rad_clear_toc )
            
            kernel_nir_terr_cod = rad_to_refl * (rad_abv_cld * lut_data %  Dems_Dcod & 
               + trans_abv_cld * Planck_Rad * lut_data % Dems_Dcod  &
               + rad_clear_toc *lut_data % Dtrnems_Dcod)
            
            
             
            kernel_nir_terr_cps  =  rad_to_refl * (rad_abv_cld * lut_data % Dems_Dcps & 
                + trans_abv_cld * Planck_Rad * lut_data % Dems_Dcps  &
                + rad_clear_toc  *lut_data %Dtrnems_Dcps )
                       
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

   !   call Get_Lut_Data ( 1 , phase_idx &
    !              , cld_alb_3d_lut = cld_alb_3d_lut )
           																
	 !  cld_alb_2d =  cld_alb_3d_lut  (  CPS_Pos : CPS_Pos+1 &
    !                                 , COD_Pos : COD_Pos+1 &
    !                                 , Sol_Pos  )
	  	  
		  
	 !  ref_diff = 0.2
	 !  cod_diff = 0.1
	 !  call interpolate_2d( cld_alb_2d , cps_wgt , cod_wgt , ref_diff , cod_diff , cld_alb , dumm_1 , dumm_2 ) 
	 
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
      
   !   call GET_LUT_DATA(channels(2), phase_num &
   !     , CPS_Vec_LUT  = CPS_Vec_LUT &
   !     , COD_Vec_LUT  = COD_Vec_LUT &
   !     , EMS_3D_LUT = EMS_3D_LUT &   
   !     ,  TRANS_ems_3d_LUT = TRANS_3D_LUT &
   !      , Refl_5D_LUT    = REFL_5D_LUT)
   
      ! - build sub tables     

   !   Ems_LUT_Sub   = EMS_3D_LUT(:,COD_Pos,Sat_Pos )
            
   !   Refl_LUT_Sub      = REFL_5D_LUT(: &
    !                            ,COD_Pos &
    !                            ,Sol_Pos &
    !                            ,Sat_Pos &
     !                           ,Azi_Pos)
	!	Trans_zen_LUT_Sub = Trans_3D_LUT(: ,COD_Pos,Sat_Pos )

      !-- toc   
     
   !   Refl = refl_lut_sub  
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
