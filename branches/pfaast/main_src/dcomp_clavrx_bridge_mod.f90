! $Header$ 
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: dcomp_clavrx_bridge_mod.f90 (src)
!       dcomp_clavrx_bridge_mod (program)
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
!   10/21/2013 : bridge to array instead to pixel-based    
!--------------------------------------------------------------------------------------
module dcomp_clavrx_bridge_mod


   ! -- MODULES USED

   use constants, only: &
        REAL4 , INT4, INT2, INT1 , SYM &
      , MISSING_VALUE_REAL4 , PI &
      , dcomp_version   ! -- this is not a constant!
   
   !use rtm_common
   
    use dcomp_interface_types_mod , only: &
         dcomp_in_type &
       , dcomp_out_type &
       , alloc_dcomp &
       , deallocate_dcompin &
       , deallocate_dcompout
   
   use  clavrx_message_module, only: &
        mesg
   
    use pixel_common, only: &
         ch &
       , geo &
       , sfc &
       , sensor &
       , image &
       , acha &
       , cld_type &
       , cld_mask &
       , bad_pixel_mask &
       , lwp_dcomp, reff_dcomp, tau_dcomp, iwp_dcomp &
       , tau_dcomp_qf, reff_dcomp_qf &
       , tau_dcomp_cost , reff_dcomp_cost &
       , dcomp_info_flag, dcomp_quality_flag &
       , cloud_063um_transmission_solar &
       , cloud_063um_transmission_view &
       , cloud_063um_spherical_albedo &
       , cloud_063um_albedo &
       , ancil_data_dir &
       , ch &
       , dcomp_mode &
       , zen_idx_rtm &
       , solar_rtm
   
   use pixel_common, only: &
        dcomp_diag_1 , dcomp_diag_2 , dcomp_diag_3 , dcomp_diag_4 &
         ,  dcomp_diag_toc_rfl1,  dcomp_diag_toc_rfl2 ,dcomp_diag_virt_alb1 &
         ,dcomp_diag_virt_alb2,dcomp_diag_wv1,dcomp_diag_wv2
   
   use calibration_constants, only: &
        sun_earth_distance  &      !---- check, This is defined in three routines   
      , solar_ch20_nu
     
   use dcomp_rtm_module
          
   implicit none

   private
   public :: set_dcomp_version              
   public :: awg_cloud_dcomp_algorithm
   
   character(len = 200) :: DCOMP_RELEASE_VERSION = 'DCOMP version 2_0_0'
  
       
contains

   !----------------------------------------------------------------------
   !  AWG_CLOUD_DCOMP_ALGORITHM
   !    This is the DCOMP bridge from CLAVR-x
   !---------------------------------------------------------------------- 
   subroutine awg_cloud_dcomp_algorithm (  iseg_in , dcomp_run )   
       
      implicit none
 
      !--- input
      integer, intent(in),optional:: iseg_in
      
      ! - output 
      logical , intent(out) :: dcomp_run
      
      type(dcomp_rtm_type) :: dcomp_rtm
      type(dcomp_in_type)  :: dcomp_input
      type(dcomp_out_type) :: dcomp_output
      
      integer :: debug_mode      
      integer :: dim_1, dim_2
      integer :: idx_chn
      
      integer :: dcomp_possible_channels ( 5) 
      integer :: i
      integer :: CHN_VIS, CHN_NIR
      
      interface
         subroutine dcomp_array_loop (a , b , debug_mode_user)
            import dcomp_in_type
            import dcomp_out_type
            type (dcomp_in_type) , intent(in) :: a
            type (dcomp_out_type), intent(out) :: b
            integer , optional :: debug_mode_user
         end subroutine

      end interface
      
      
      ! ----- executable  --------------------------------------------------- !
            
      dcomp_run = .false.
      
      ! - do we need to run dcomp at all? ( night  etc..)
      if ( count ( geo % solzen < 75. .and. geo % solzen >= 0 .and. geo % satzen < 75. ) < 1 ) return
      dcomp_run = .true.
      
      
      if ( iseg_in == 1 ) then
        call mesg ('dcomp start ' ) 
      end if
      
      
      ! - compute DCOMP related RTM 
      call perform_rtm_dcomp ( dcomp_rtm ) 
            
      dim_1 = Image%Number_Of_Elements
      dim_2 = Image%Number_Of_Lines_Read_This_Segment
  
      ! - which channels do we need? possibles are 
      dcomp_input % is_channel_on = .false.
      
      dcomp_possible_channels = [ 1, 5, 6, 7, 20 ]
      do i = 1 , size ( dcomp_possible_channels )   
         if ( sensor % chan_on_flag_default ( dcomp_possible_channels ( i) ) == 1 ) then
            dcomp_input % is_channel_on (dcomp_possible_channels ( i)  )  = .true.
         end if
      end do 
      
      !- check mode
      CHN_VIS = 1
      select case ( dcomp_mode )
         case ( 1 ) 
            CHN_NIR = 6               
         case ( 2 )
            CHN_NIR = 7
         case ( 3 )
            CHN_NIR = 20
         case default
            print*, 'dcomp mode ',dcomp_mode,' not possible'
            return
      end select
      
      if ( dcomp_input % is_channel_on (CHN_NIR) .eqv. .false.) then
         if ( iseg_in == 1 ) then
            print*,'dcomp NIR channel is not set! ==> MODIS equaivalant channel: ', CHN_NIR
            call mesg ( 'all dcomp results are set to missing values', color=41 , level = -1 ) 
         end if
         return
      end if
   
      if ( dcomp_input % is_channel_on (CHN_VIS) .eqv. .false.) then
         if ( iseg_in == 1 ) then
            print*, 'dcomp VIS channel is not set! ==> MODIS equaivalant channel: ', CHN_VIS
            call mesg ( 'all dcomp results are set to missing values', color=41 , level = -1 ) 
         end if
         return
      end if
   
      
      ! === ALLOCATION
      do idx_chn = 1 , 40
         if ( dcomp_input % is_channel_on (idx_chn) .eqv. .false.) cycle
        
         call  alloc_dcomp ( dcomp_input % refl    (  idx_chn  ) , dim_1,dim_2 ) 
         call  alloc_dcomp ( dcomp_input % alb_sfc (  idx_chn  ) ,  dim_1,dim_2 ) 
                 
         if ( idx_chn >= 20 ) then   
            call  alloc_dcomp ( dcomp_input % rad (  idx_chn  ) ,  dim_1,dim_2 )  
            call  alloc_dcomp ( dcomp_input % emiss_sfc (  idx_chn  ) ,  dim_1,dim_2 )
            call  alloc_dcomp ( dcomp_input % rad_clear_sky_toa ( idx_chn ),  dim_1,dim_2 )
            call  alloc_dcomp ( dcomp_input % rad_clear_sky_toc ( idx_chn ), dim_1,dim_2 )
            call alloc_dcomp (  dcomp_input % trans_ac_nadir ( idx_chn )   , dim_1,dim_2 )
         end if   
      end do
      
      call  alloc_dcomp ( dcomp_input % snow_class,   dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % is_land,      dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % cloud_press,  dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % cloud_temp,   dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % cloud_hgt,    dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % cloud_type,   dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % cloud_mask,   dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % tau_acha,   dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % ozone_nwp,    dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % tpw_ac,       dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % press_sfc,    dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % is_valid,     dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % sol,          dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % sat,          dim_1,dim_2 )
      call  alloc_dcomp ( dcomp_input % azi,          dim_1,dim_2 )
      
      ! == CONFIGURE 
      
         ! - dcomp-mode
      dcomp_input % mode = dcomp_mode
         ! - ancil/lut path
      dcomp_input % lut_path = trim(ancil_data_dir)//"/static/luts/cld/"   
         ! - wmo sensor id
      dcomp_input % sensor_wmo_id = sensor % wmo_id
      dcomp_input % sun_earth_dist = sun_earth_distance
            
         ! -  Satellite Data
      dcomp_input % refl ( 1 ) % d = ch(1)%ref_toa
      
      if ( dcomp_input % is_channel_on (5))  dcomp_input % refl ( 5 ) % d = ch(5)%ref_toa
      if ( dcomp_input % is_channel_on (6))  dcomp_input % refl ( 6 ) % d = ch(6)%ref_toa
      if ( dcomp_input % is_channel_on (7))  dcomp_input % refl ( 7 ) % d = ch(7)%ref_toa
      if ( dcomp_input % is_channel_on (20)) dcomp_input % rad ( 20 ) % d = ch(20)%rad_toa
      
      dcomp_input % sat % d = geo % satzen
      dcomp_input % sol % d = geo % solzen
      dcomp_input % azi % d = geo % relaz
              
         ! - Cloud products
      dcomp_input % cloud_press % d = acha % pc
      dcomp_input % cloud_temp % d  = acha % tc
      dcomp_input % tau_acha % d    = acha % tau
      dcomp_input % cloud_mask % d  = cld_mask
      dcomp_input % cloud_type % d  = cld_type
                  
         ! - Flags
      dcomp_input % is_land % d = sfc % land_mask == 1 
      dcomp_input % is_valid % d = bad_pixel_mask /= 1
      
            
         ! - Surface  (AKH - What about Channel 5)
      dcomp_input % alb_sfc ( 1 ) % d = ch(1) % sfc_ref_white_sky
      if ( dcomp_input % is_channel_on (6 )) dcomp_input % alb_sfc ( 6 ) % d = ch(6) % sfc_ref_white_sky
      if ( dcomp_input % is_channel_on (7 )) dcomp_input % alb_sfc ( 7 ) % d = ch(7) % sfc_ref_white_sky
      if ( dcomp_input % is_channel_on (20)) dcomp_input % alb_sfc ( 20) % d = 100.0*(1.0 - ch(20)%sfc_emiss)    
      if ( dcomp_input % is_channel_on (20)) dcomp_input % emiss_sfc ( 20) % d = ch(20)%sfc_emiss
      dcomp_input % press_sfc % d =  dcomp_rtm % sfc_nwp
      dcomp_input % snow_class % d = sfc % snow
            
         ! - Atmospheric contents
         ! ozone column in Dobson
      dcomp_input % ozone_nwp % d = dcomp_rtm % ozone_path
         ! Total water Vapour above the cloud
      dcomp_input % tpw_ac % d = dcomp_rtm % tpw_ac
           
         ! - RTM coeffeicience to compute transmission in non-IR channels
      dcomp_input % gas_coeff(1) % d = solar_rtm % tau_h2o_coef(1,:)
      dcomp_input % gas_coeff(6) % d = solar_rtm % tau_h2o_coef(6,:)
      dcomp_input % gas_coeff(7) % d = solar_rtm % tau_h2o_coef(7,:)
      
         ! -- transmission above the cloud in channel 20 (3.9) from RTM  
      if ( dcomp_input % is_channel_on (20)) then       
         dcomp_input % trans_ac_nadir (20) % d = dcomp_rtm % trans_ir_ac_nadir
         dcomp_input % rad_clear_sky_toc (20) % d = dcomp_rtm % rad_clear_sky_toc_ch20
         dcomp_input % rad_clear_sky_toa (20) % d = dcomp_rtm % rad_clear_sky_toa_ch20
         ! -- Solar irradiance in channel 20
         dcomp_input % solar_irradiance ( 20) = solar_ch20_nu
      end if
      
      call dcomp_rtm % deallocate_it()
      
      ! === THE MAIN CALL of DCOMP ===          
      debug_mode = 1
      call dcomp_array_loop ( dcomp_input , dcomp_output , debug_mode_user = debug_mode)
      
      ! === DEALLOCATION
      call deallocate_dcompin ( dcomp_input )
      
      ! === POPULATE CLAVR-X VARIABLES FROM PIXEL_COMMON
      tau_dcomp (1:dim_1,1:dim_2)   = dcomp_output % cod % d
      reff_dcomp  (1:dim_1,1:dim_2) = dcomp_output % cps % d
      lwp_dcomp (1:dim_1,1:dim_2)   = dcomp_output % lwp % d
      iwp_dcomp (1:dim_1,1:dim_2)   = dcomp_output % iwp % d
      
      tau_dcomp_cost(1:dim_1,1:dim_2)     = dcomp_output % cod_unc % d
      reff_dcomp_cost(1:dim_1,1:dim_2)    = dcomp_output % ref_unc % d
      dcomp_quality_flag(1:dim_1,1:dim_2) = dcomp_output % quality % d
      dcomp_info_flag(1:dim_1,1:dim_2)    = dcomp_output % info % d
    !  
      cloud_063um_transmission_solar(1:dim_1,1:dim_2) = dcomp_output % cld_trn_sol % d
      cloud_063um_transmission_view(1:dim_1,1:dim_2)  = dcomp_output % cld_trn_obs % d
      cloud_063um_albedo(1:dim_1,1:dim_2)             = dcomp_output % cld_alb % d
      cloud_063um_spherical_albedo(1:dim_1,1:dim_2)   = dcomp_output % cld_sph_alb % d
      
      DCOMP_RELEASE_VERSION = dcomp_output % version
      
      call deallocate_dcompout ( dcomp_output)
      
     
   end subroutine awg_cloud_dcomp_algorithm

  
   !---------------------------------------------------------------------------
   ! routine to set the cvs version in a global variable to write to hdf file
   !---------------------------------------------------------------------------
   subroutine set_dcomp_version()
      dcomp_version = DCOMP_RELEASE_VERSION
   end subroutine set_dcomp_version

    
end module dcomp_clavrx_bridge_mod 

