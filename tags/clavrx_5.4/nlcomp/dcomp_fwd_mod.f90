module dcomp_fwd_mod
   use dcomp_lut_mod
   use dcomp_tools
   
   integer, parameter, public:: int1 = selected_int_kind(1)
   integer, parameter, public:: int2 = selected_int_kind(3)
   integer, parameter, public:: int4 = selected_int_kind(8)
   integer, parameter, public:: int8 = selected_int_kind(10)
   integer, parameter, public:: real4 = selected_real_kind(6,37)
   integer, parameter, public:: real8 = selected_real_kind(15,307)
   integer, parameter, public:: ipre = real4

   real , parameter :: Pi = 3.14159265
   real , parameter :: dtor = pi / 180.
   
   type pixel_vec
      real :: sol_zen
      real :: sat_zen
      real :: rel_azi
      real :: ctt
      logical :: is_water_phase
   end type pixel_vec

contains
   subroutine get_forward_model ( &
                              state_vec &
							 , pixel &
							 , sensor &
							 , channels &
							 , alb_sfc &
							 , air_trans_ac &
							 , fwm_vec &
							 , cld_trns_sol &
							 , cld_trns_sat &
							 , cld_sph_alb &
							 , kernel &
							 , rad_abv_cld &
							 , rad_clear_toc )
      
	  implicit none
	  
	  real , dimension ( 2 ) , intent ( in ) :: state_vec
	  type ( pixel_vec ) , intent (in ) :: pixel 
	  character ( len = * ), intent ( in ) :: sensor
	  integer , dimension ( : ) , intent ( in ) :: channels
	  real , dimension ( : ) , intent ( in ) :: alb_sfc
	  real , dimension ( : ) , intent ( in ) :: air_trans_ac
	  real , dimension ( : ) , intent ( out ) :: fwm_vec
	  real , dimension ( : ) , intent ( out ) :: cld_trns_sol
	  real , dimension ( : ) , intent ( out ) :: cld_trns_sat
	  real , dimension ( : ) , intent ( out ) :: cld_sph_alb
	  real , dimension ( : , : ) , intent ( out ) :: kernel
	  real, intent ( in ), optional ::  rad_abv_cld ,  rad_clear_toc 
	   							  
      real, dimension(:,:), POINTER::  Sph_Alb_2D_Lut					 
      integer :: n_channels
      integer :: i
	  
	  logical :: read_angle_val
	  integer :: phase_num
	  integer , save :: sol_zen_count, sat_zen_count, rel_azi_count
	  integer , save :: cps_count , cod_count
	  
	  real, dimension(:),pointer, save :: Sol_Zen_Vec_LUT !- solar zenith angle vector   
      real, dimension(:),pointer, save :: Sat_Zen_Vec_LUT !- solar zenith angle vector
      real, dimension(:),pointer, save :: Rel_Azi_Vec_LUT !- relative azimuth vector
	  real , dimension(:), pointer , save :: CPS_Vec_LUT  !- CPS vector in LUT
      real , dimension(:), pointer , save :: COD_Vec_LUT  !- COD vector in LUT
	  real, save :: sol_zen_last, sat_zen_last, rel_azi_last
	  real :: sol_wgt, sol2_wgt, sat_wgt, azi_wgt
	  integer :: sol_pos_tmp, sol2_pos_tmp , sat_pos_tmp , azi_pos_tmp
	  integer , save :: sol_pos, sol2_pos , sat_pos , azi_pos
	  real, save :: cps_wgt, cod_wgt
	  integer , save :: cps_pos, cod_pos
	  
	  real :: air_mass_two_way
	  real, dimension(30) :: trans_two_way
	  
	  real :: planck_rad
	  real :: get_planck_radiance_39um
	  real :: rad_to_refl , get_rad_refl_factor
	  real :: trans_abv_cld
      real :: air_mass_sat
	  
	  real :: fm_nir_terr , fm_nir_solar
      real :: kernel_nir_solar_cod  , kernel_nir_terr_cod
      real :: kernel_nir_solar_cps  , kernel_nir_terr_cps
	  
	  ! - executable
	  n_channels = size ( channels )
	  
	  call populate_all_lut ( sensor , channels )
	  
	  
	  ! find angle positions in lut
	  phase_num = 2
      if ( pixel % is_water_phase ) phase_num = 1
	  if ( read_angle_val  .eqv. .false. ) then
         call Get_Lut_Data ( 1 , phase_num &
            & , num_sol_zen = sol_zen_count &
            & , num_sat_zen = sat_zen_count &
            & , num_rel_azi = rel_azi_count &
            & , Sol_Zen_Vec_LUT = Sol_Zen_Vec_LUT &
            & , Sat_Zen_Vec_LUT = Sat_Zen_Vec_LUT &
            & , Rel_Azi_Vec_LUT = Rel_Azi_Vec_LUT )   
	     call Get_lut_data ( 1, phase_num  &
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
         
         call dcomp_Interpolation_Weight ( Sol_Zen_count , pixel % Sol_Zen  &
                                     & , Sol_Zen_Vec_LUT , Sol_Wgt , Sol_Pos_tmp , sol_pos)  
									            
         call dcomp_Interpolation_Weight ( Sol_Zen_count , pixel % Sat_Zen &
                                     & , Sol_Zen_Vec_LUT , Sol2_Wgt , Sol2_Pos_tmp , sol2_pos )  
									            
         call dcomp_Interpolation_Weight ( Sat_Zen_count , pixel % Sat_Zen &
                                     & , Sat_Zen_Vec_LUT , Sat_Wgt , Sat_Pos_tmp , sat_pos)
									 
         call dcomp_Interpolation_Weight ( Rel_Azi_count , pixel % Rel_Azi  &
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
                                                                       
	  air_mass_two_way = ( 1. / cos ( pixel % Sol_zen * dtor ) &
	            + 1./ cos ( pixel % sat_zen * dtor) )
	  do i = 1 , n_channels
	     call compute_fm_single_channel (  channels ( i )   , phase_num  &
                                     & , cps_wgt , cod_wgt  &
                                     & , cps_pos , cod_pos, sol_pos , sat_pos , azi_pos & 
                                     & , alb_sfc( i ) &
                                     & , fwm_vec(i) , kernel(i,1), kernel(i,2) &
                                     & , cld_trns_sol(i), cld_trns_sat(i) , cld_sph_alb(i) )
         print*,i
         trans_two_way = air_trans_ac ** air_mass_two_way  
         fwm_vec(i) = fwm_vec(i) * trans_two_way(i)
	     kernel(i,1) = kernel(i,1) * trans_two_way(i)
         kernel(i,2) = kernel(i,2) * trans_two_way(i)
		 
		 if ( channels ( i ) .eq. 20 ) then
		    planck_rad = get_planck_radiance_39um ( pixel % ctt  , trim(sensor) )
			
            rad_to_refl = get_rad_refl_factor ( trim ( sensor ) , pixel % sol_zen ) 
			
            air_mass_sat =  1./ cos ( pixel % Sat_zen * dtor) 
			
            trans_abv_cld = air_trans_ac(i) ** air_mass_sat
            
            call compute_fm_terrestrial ( channels ( i) , phase_num  &
                                     & , cps_wgt , cod_wgt   &
                                     & , cps_pos , cod_pos, sat_pos  & 
                                     & , planck_rad , rad_abv_cld , trans_abv_cld , rad_clear_toc , rad_to_refl &
                                     & , fm_nir_terr , kernel_nir_terr_cod ,  kernel_nir_terr_cps )
                                      
            fwm_vec(i) = fwm_vec(i) + fm_nir_terr
            kernel ( i, 1) = kernel ( i , 1)   + kernel_nir_terr_cod 
            kernel ( i, 2) = kernel ( i , 2 )  + kernel_nir_terr_cps 
		  
		  
		  end if
		  
		  
	  end do
	  
	 
   
   end subroutine get_forward_model

   ! ------
   !
   ! -----
   subroutine compute_fm_single_channel  ( chn_idx  , phase_idx  &
                         & , cps_wgt , cod_wgt  &
                                     & , cps_pos , cod_pos, sol_pos , sat_pos , azi_pos & 
                         & , alb_sfc &
                                     & , refl , drefl_dcod, drefl_dcps &
                                     & , trns_sol, trns_sat , sph_alb )
   
       implicit none
      !---input
      integer, intent(in) :: chn_idx             ! - Channel index
      integer, intent(in) :: phase_idx           ! - Phase index
      integer, intent(in) :: cps_pos     ! - Position  in LUT
      integer, intent(in) :: cod_pos     ! - Position  in LUT
         integer, intent(in) :: sol_pos, sat_pos, azi_pos
        
      real, intent(In) :: cps_wgt        ! - Weight beween positionn and next index in LUT
      real, intent(In) :: cod_wgt        ! - Weight beween positionn and next index in LUT
      real, intent(In) :: alb_Sfc                ! - Surface albedo; value between 0 and 1
   
       !--output
      REAL, intent(Out):: Refl                         ! - Reflectance at TOA under consideration of surface impact
      REAL, intent(Out):: Drefl_Dcod           ! - Derivation of reflectance to delta cloud optical depth
      REAL, intent(Out):: Drefl_Dcps           ! - Derivation of reflectance to delta cloud particle size

      REAL, intent(Out):: Trns_sol             ! - Cloud transmission for radiance coming from solar angle
      REAL, intent(Out):: Trns_sat                 ! - Cloud transmission for radiance coming from satellite angle
      REAL, intent(Out):: sph_alb                   ! - Spherical Albedo
         
          !--local

      REAL:: Dtrans_sol_Dcod   ! - Derivation transmission solar angle to COD
      REAL:: Dtrans_sol_Dcps   ! - Derivation transmission solar angle to CPS
      REAL:: Dtrans_sat_Dcod   ! - Derivation transmission satellite angle to COD
      REAL:: Dtrans_sat_Dcps   ! - Derivation transmission satellite angle to CPS
      REAL:: Dsph_alb_Dcod     ! - Derivation to cod
      REAL:: Dsph_alb_Dcps     ! - Derivation to cps
      REAL:: Alb_Sfc_Term      ! - Term which accounts for albedo surface 
      REAL:: Cld_Refl                  ! - Cloud reflectance after interpolation


     

      ! - data which comes from LUT
      REAL, dimension(:,:), POINTER::  Sph_Alb_2D_Lut
      REAL, dimension(:,:,:), POINTER:: Trans_3D_Lut
      REAL, dimension(:,:,:,:,:), POINTER:: Refl_5D_Lut
	  
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
   subroutine compute_fm_terrestrial (&
                         chn_idx , phase_idx  &
                         & , cps_wgt , cod_wgt  &
                         & , cps_pos , cod_pos , sat_pos &
                         & , planck_rad ,  rad_clear_ct , trans_abv_cld , rad_clear_toa , rad_to_refl &
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
      real, intent (in) :: planck_rad
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
	 								

     
print*,cld_emiss,planck_rad, trans_abv_cld, cld_trans
print*, rad_clear_toa
     
      rad = rad_clear_ct * Cld_Emiss &
                  + trans_abv_cld * cld_emiss * planck_rad &
				  + cld_trans * rad_clear_toa                 
            
      dRad_dcod =  rad_clear_ct * Dems_Dcod & 
               + trans_abv_cld * Planck_Rad * Dems_Dcod  &
               + rad_clear_toa * Dtrans_Dcod

      dRad_dcps =   rad_clear_ct * Dems_Dcps & 
                + trans_abv_cld * Planck_Rad * Dems_Dcps  &
                + rad_clear_toa  *Dtrans_Dcps        
                            
                            
      refl_terr = rad * rad_to_refl
         
      drefl_dcod_terr = dRad_dcod * rad_to_refl
      drefl_dcps_terr = dRad_dcps *  rad_to_refl
    
   end subroutine compute_fm_terrestrial
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
! INTERPOLATE_2D
!
! Linear Interpolation for a 2 x 2 
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

! - executable


!--- spherical albedo
    
      value  =     ( 1.0 - r ) * ( 1.0 - s ) * table( 1 , 1 ) &
	              + r * ( 1.0 - s )          * table( 2 , 1 )  &
                  + (1.0 - r ) * s           * table( 1 , 2 ) &
				  + r * s                    * table( 2 , 2 )

      dval_d2 =  ( ( 1.0 - r ) * ( table( 1 , 2 ) - table( 1 , 1 ) )   &
                    +  r * ( table( 2 , 2 ) - table( 2 , 1 ) ) ) /  &
                       delta_2

      dval_d1 =  ( ( 1.0 - s ) * ( table( 2 , 1 ) - table( 1 , 1 ) )  &
                    + s * ( table( 2 , 2 ) - table( 1 , 2 ) ) ) /  &
                        delta_1 


                  
   end subroutine interpolate_2d  


   


end module dcomp_fwd_mod
