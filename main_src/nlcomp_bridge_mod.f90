! $Header$ 
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: nlcomp_bridge_mod.f90 (src)
!       NLCOMP_BRIDGE_MOD (program)
!
! PURPOSE: 
!
! DESCRIPTION: 
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
! REVISION HISTORY:
!   02/10.2013 : first version
!    01/13/2014 :  This bridge tool calls now array Nlcomp instead of pixel nlcomp (AW)
!--------------------------------------------------------------------------------------
module nlcomp_bridge_mod


   ! -- MODULES USED

   use constants, only: &
        REAL4 , INT4, INT2, INT1 , SYM &
      , MISSING_VALUE_REAL4 , PI &
      , dcomp_version   ! -- this is not a constant!
   
   use rtm_common
!  use rt_utilities, only: &
!        nlevels_rtm &
!       , p_std_rtm
   
   use  clavrx_message_module, only: &
        mesg
   
   use pixel_common, only: &
        sensor &
      , geo &
      , sfc &
      , image &
      , ch  &
      , acha &
      , cld_type &
      , cld_mask &
      , bad_pixel_mask &
      , reff_nlcomp, tau_nlcomp  &
      , tau_dcomp_qf, reff_dcomp_qf &
      , tau_nlcomp_cost , reff_nlcomp_cost &
      , nlcomp_info_flag, nlcomp_quality_flag &
      , cloud_063um_transmission_solar &
      , cloud_063um_transmission_view &
      , cloud_063um_spherical_albedo, cloud_063um_albedo &
      , ancil_data_dir &
      , i_nwp , j_nwp &
      , ch &
      , dcomp_mode  &
      , zen_idx_rtm , solar_rtm
   
   use pixel_common, only: &
        dcomp_diag_1 , dcomp_diag_2
   
   use calibration_constants, only: &
        sun_earth_distance  &      !---- check this is defined in three routines   
      , solar_ch20_nu
   
   use nwp_common, only: &
       t_prof_nwp , z_prof_nwp , tpw_prof_nwp &
      , p_std_nwp , sfc_level_nwp , tropo_level_nwp &
      , inversion_level_nwp, psfc_nwp , ozone_nwp

!  use rt_utilities, only: &
!     rtm

   use planck, only: &
      planck_rad_fast 

   use dcomp_rtm_module
#ifdef HDF5LIBS   

   use nlcomp_interface_def_mod , only: &
         nlcomp_in_type &
       , nlcomp_out_type &
       , alloc_nlcomp &
       , deallocate_nlcompin &
       , deallocate_nlcompout
   implicit none

#endif
   private
   
   logical :: first_call = .true.


   !--- module main subroutine 
   public:: awg_cloud_nlcomp_algorithm

contains

   !----------------------------------------------------------------------
   !
   !---------------------------------------------------------------------- 
   subroutine awg_cloud_nlcomp_algorithm ( iseg_in  )   
 
      implicit none
 
      !--- input
      integer, intent(in),optional:: iseg_in

      integer :: dim_1, dim_2
      integer :: idx_chn
      
      integer :: nlcomp_possible_channels ( 2 ) 
      integer :: i
      
      real , parameter :: PI = 3.1415927
      real, parameter :: DTOR = PI/180.
      
   
     
   

       type(dcomp_rtm_type) :: nlcomp_rtm
      type(nlcomp_in_type)  :: nlcomp_input
      type(nlcomp_out_type) :: nlcomp_output
      
      interface
         subroutine nlcomp_array_loop_sub ( a, b, debug_mode_user)
            import nlcomp_in_type
            import nlcomp_out_type
            type ( nlcomp_in_type ), intent(in) :: a
            type ( nlcomp_out_type ), intent(in) :: b
            integer , optional :: debug_mode_user
         end subroutine
      end interface
     

      ! ----- executable  --------------------------------------------------- !


      if ( first_call ) then
        call mesg ('NL-COMP starts ... ', color=46  ) 
        first_call = .false.
      end if
      
     
      ! - compute DCOMP related RTM 
      call perform_rtm_dcomp ( nlcomp_rtm ) 
      
      dim_1 = Image%Number_Of_Elements
      dim_2 = Image%Number_Of_Lines_Read_This_Segment
      
      
      ! - which channels do we need? possibles are 
      nlcomp_input % is_channel_on = .false.
      
      nlcomp_possible_channels = [  20 , 44 ]
      do i = 1 , size ( nlcomp_possible_channels )   
         if ( sensor % chan_on_flag_default ( nlcomp_possible_channels ( i) ) == 1 ) then
            nlcomp_input % is_channel_on (nlcomp_possible_channels ( i)  )  = .true.
         end if
      end do 
      
      
      
       ! === ALLOCATION
      do idx_chn = 1 , 44
         if ( nlcomp_input % is_channel_on (idx_chn) .eqv. .false.) cycle
        
         call  alloc_nlcomp ( nlcomp_input % refl    (  idx_chn  ) , dim_1,dim_2 ) 
         call  alloc_nlcomp ( nlcomp_input % alb_sfc (  idx_chn  ) ,  dim_1,dim_2 ) 
                 
         if ( idx_chn >= 20 .and. idx_chn /= 40 ) then   
            call  alloc_nlcomp ( nlcomp_input % rad (  idx_chn  ) ,  dim_1,dim_2 )  
            call  alloc_nlcomp ( nlcomp_input % emiss_sfc (  idx_chn  ) ,  dim_1,dim_2 )
            call  alloc_nlcomp ( nlcomp_input % rad_clear_sky_toa ( idx_chn ),  dim_1,dim_2 )
            call  alloc_nlcomp ( nlcomp_input % rad_clear_sky_toc ( idx_chn ), dim_1,dim_2 )
            call alloc_nlcomp (  nlcomp_input % trans_ac_nadir ( idx_chn )   , dim_1,dim_2 )
         end if   
      end do
      
      call  alloc_nlcomp ( nlcomp_input % snow_class,   dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % is_land,      dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % cloud_press,  dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % cloud_temp,   dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % cloud_hgt,    dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % cloud_type,   dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % cloud_mask,   dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % tau_acha,   dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % ozone_nwp,    dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % tpw_ac,       dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % press_sfc,    dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % is_valid,     dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % sol,          dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % sat,          dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % azi,          dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % zen_lunar,          dim_1,dim_2 )
      call  alloc_nlcomp ( nlcomp_input % azi_lunar,          dim_1,dim_2 )
      
      ! -- configure
      
           ! - ancil/lut path
      nlcomp_input % lut_path = trim(ancil_data_dir)//"static//luts/cld/"   
         ! - wmo sensor id
      nlcomp_input % sensor_wmo_id = sensor % wmo_id
      nlcomp_input % sun_earth_dist = sun_earth_distance
            
         ! -  Satellite Data
     
      
     
      if ( nlcomp_input % is_channel_on (20)) nlcomp_input % rad ( 20 ) % d = ch(20)%rad_toa
      if ( nlcomp_input % is_channel_on (44)) nlcomp_input % refl ( 44 ) % d = ch(44)%ref_lunar_toa
      if ( nlcomp_input % is_channel_on (44)) nlcomp_input % rad ( 44 ) % d = ch(44)%rad_toa
      
      nlcomp_input % sat % d = geo % satzen
      nlcomp_input % sol % d = geo % solzen
      nlcomp_input % azi % d = geo % relaz
      nlcomp_input % zen_lunar % d = geo % lunzen
      nlcomp_input % azi_lunar % d = geo % lunrelaz
              
         ! - Cloud products
      nlcomp_input % cloud_press % d = acha % pc
      nlcomp_input % cloud_temp % d  = acha % tc
      nlcomp_input % tau_acha % d    = acha % tau
      nlcomp_input % cloud_mask % d  = cld_mask
      nlcomp_input % cloud_type % d  = cld_type
                  
         ! - Flags
      nlcomp_input % is_land % d = sfc % land_mask == 1 
      nlcomp_input % is_valid % d = bad_pixel_mask /= 1
      
            
         ! - Surface  
      ! - we use for DNB (channel 44) channel 1  
      if ( nlcomp_input % is_channel_on (44 )) nlcomp_input % alb_sfc ( 44 ) % d = ch(1)%sfc_ref_white_sky  
     

      if ( nlcomp_input % is_channel_on (20)) nlcomp_input % alb_sfc ( 20) % d = 100.0*(1.0 - ch(20)%sfc_emiss)    !check this AKH
      if ( nlcomp_input % is_channel_on (20)) nlcomp_input % emiss_sfc ( 20) % d = ch(20)%sfc_emiss
      nlcomp_input % press_sfc % d =  nlcomp_rtm % sfc_nwp
      nlcomp_input % snow_class % d = sfc % snow
            
         ! - Atmospheric contents
         ! ozone column in Dobson
      nlcomp_input % ozone_nwp % d = nlcomp_rtm % ozone_path
         ! Total water Vapour above the cloud
      nlcomp_input % tpw_ac % d = nlcomp_rtm % tpw_ac
           
         ! - RTM coeffeicience to compute transmission in non-IR channels
      nlcomp_input % gas_coeff(1) % d = solar_rtm % tau_h2o_coef(1,:)
      
      nlcomp_input % gas_coeff(6) % d = solar_rtm % tau_h2o_coef(6,:)
      nlcomp_input % gas_coeff(7) % d = solar_rtm % tau_h2o_coef(7,:)
      
         ! -- transmission above the cloud in channel 20 (3.9) from RTM  
      if ( nlcomp_input % is_channel_on (20)) then       
         nlcomp_input % trans_ac_nadir (20) % d = nlcomp_rtm % trans_ir_ac_nadir
         nlcomp_input % rad_clear_sky_toc (20) % d = nlcomp_rtm % rad_clear_sky_toc_ch20
         nlcomp_input % rad_clear_sky_toa (20) % d = nlcomp_rtm % rad_clear_sky_toa_ch20
         ! -- Solar irradiance in channel 20
         nlcomp_input % solar_irradiance ( 20) = solar_ch20_nu
      end if
      
      call nlcomp_array_loop_sub (nlcomp_input,nlcomp_output )
       
     
      
      ! === POPULATE CLAVR-X VARIABLES FROM PIXEL_COMMON
      tau_nlcomp (1:dim_1,1:dim_2)   = nlcomp_output % cod % d
      reff_nlcomp  (1:dim_1,1:dim_2) = nlcomp_output % cps % d
      tau_nlcomp_cost(1:dim_1,1:dim_2) = nlcomp_output % cod_unc % d
      reff_nlcomp_cost(1:dim_1,1:dim_2) = nlcomp_output % ref_unc % d
      nlcomp_quality_flag(1:dim_1,1:dim_2) = nlcomp_output %  quality % d
      nlcomp_info_flag(1:dim_1,1:dim_2) = nlcomp_output % info % d
      
      call deallocate_nlcompout ( nlcomp_output)
  
   
   end subroutine awg_cloud_nlcomp_algorithm

  

end module nlcomp_bridge_mod 

