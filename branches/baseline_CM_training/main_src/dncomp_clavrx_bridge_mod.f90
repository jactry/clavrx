! $Header: https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/baseline_CM_training/main_src/dncomp_clavrx_bridge_mod.f90 1757 2016-09-09 20:46:30Z dbotambekov $ 
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
   
    use dncomp_interface_def_mod , only: &
         dncomp_in_type &
       , dncomp_out_type &
       , alloc_dncomp &
       , n_chn 
   
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
       , solar_rtm &
       , tau_nlcomp &
       , reff_nlcomp &
       , tau_nlcomp_cost &
       , reff_nlcomp_cost &
       , nlcomp_quality_flag &
       , nlcomp_info_flag
   
  !!! use pixel_common, only: &
   !      dcomp_diag_2 , dcomp_diag_3 , dcomp_diag_4 &
   !      ,  dcomp_diag_toc_rfl1,  dcomp_diag_toc_rfl2 ,dcomp_diag_virt_alb1 &
   !      ,dcomp_diag_virt_alb2,dcomp_diag_wv1,dcomp_diag_wv2
   
   use calibration_constants, only: &
        sun_earth_distance  &      !---- check, This is defined in three routines   
      , solar_ch20_nu
     
   use dcomp_rtm_module
          
   implicit none

   private
   
   logical :: first_call = .true.
   
   public :: set_dcomp_version              
   public :: awg_cloud_dncomp_algorithm
   
   character(len = 200) :: DCOMP_RELEASE_VERSION = 'DCOMP version 2_0_0'
  
       
contains

   !----------------------------------------------------------------------
   !  AWG_CLOUD_DCOMP_ALGORITHM
   !    This is the DCOMP bridge from CLAVR-x
   !---------------------------------------------------------------------- 
   subroutine awg_cloud_dncomp_algorithm (  iseg_in , nlcomp_mode,  algorithm_started )   
       
      implicit none
 
      !--- input
      integer, intent(in),optional :: iseg_in
      logical, intent(in),optional :: nlcomp_mode
      
      ! - output 
      logical , intent(out) :: algorithm_started
      
      type(dcomp_rtm_type) :: dcomp_rtm
      type(dncomp_in_type)  :: dcomp_input
      type(dncomp_out_type) :: dncomp_output
      
      integer :: debug_mode      
      integer :: dim_1, dim_2
      integer :: idx_chn
      
      logical :: run_nlcomp
      
      integer, allocatable :: possible_channels ( : )
      logical :: chan_on ( N_CHN ) = .false.
      integer :: i
      integer :: CHN_VIS
      integer :: CHN_NIR
      
      interface
         subroutine dcomp_array_loop (a , b , debug_mode_user)
            import dncomp_in_type
            import dncomp_out_type
            type (dncomp_in_type) , intent(in) :: a
            type (dncomp_out_type), intent(out) :: b
            integer , intent(in), optional :: debug_mode_user
         end subroutine

      end interface
      
      
      ! ----- executable  --------------------------------------------------- !
      run_nlcomp = .false.
      if (present(nlcomp_mode)) run_nlcomp = nlcomp_mode
            
      algorithm_started = .false.
      
      ! - do we need to run dcomp at all? ( night  etc..)
      
      if (run_nlcomp) then
         ! add here all conditions which leads to a immediate stop
         if ( count (ch(44)%Ref_Lunar_Toa > 0) < 1 ) return
      else      
         if ( count ( geo % solzen < 75. .and. geo % solzen >= 0 .and. geo % satzen < 75. ) < 1 ) return
      end if   
      algorithm_started = .true.
      
      
      if ( first_call) then
        call mesg ('DNCOMP starts ',color = 46 ) 
        first_call = .false.
      end if
      
      
      ! - compute DCOMP related RTM 
      call perform_rtm_dcomp ( dcomp_rtm ) 
           
      dim_1 = Image%Number_Of_Elements
      dim_2 = Image%Number_Of_Lines_Read_This_Segment
  
      chan_on = .false.
      
      if (run_nlcomp) then
         allocate (possible_channels(2))
         possible_channels =[20,44]
      else
         allocate ( possible_channels(5))
         possible_channels = [ 1, 5, 6, 7, 20 ]
      end if
      
      do i = 1 , size ( possible_channels )   
         if ( sensor % chan_on_flag_default ( possible_channels ( i) ) == 1 ) then
            chan_on (possible_channels ( i)  )  = .true.
         end if
      end do 
      
      ! - here we have to add channels for snow
     
      !-allocate input
      dcomp_input = dncomp_in_type ( dim_1, dim_2, chan_on )
      
      
      if ( .not. run_nlcomp) then      
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
      end if
      
      ! == CONFIGURE 
      
         ! - dcomp-mode
      dcomp_input % mode = dcomp_mode
         ! - ancil/lut path
      dcomp_input % lut_path = trim(ancil_data_dir)//"/static/luts/cld/"   
         ! - wmo sensor id
      dcomp_input % sensor_wmo_id = sensor % wmo_id
      dcomp_input % sun_earth_dist = sun_earth_distance
            
         ! -  Satellite Data
      
      ! - all reflectance channels
      do i = 1, 19 
         if ( dcomp_input % is_channel_on (i)) then
             dcomp_input % refl ( i ) % d = ch(i)%ref_toa
             dcomp_input % alb_sfc ( i ) % d = ch(i) % sfc_ref_white_sky
             dcomp_input % gas_coeff(i) % d = solar_rtm % tau_h2o_coef(i,:)
         end if     
      end do   
      
      if ( dcomp_input % is_channel_on (20)) then
         dcomp_input % rad ( 20 ) % d = ch(20)%rad_toa
         dcomp_input % alb_sfc ( 20) % d = 100.0*(1.0 - ch(20)%sfc_emiss)
         dcomp_input % emiss_sfc ( 20) % d = ch(20)%sfc_emiss
         dcomp_input % trans_ac_nadir (20) % d = dcomp_rtm % trans_ir_ac_nadir
         dcomp_input % rad_clear_sky_toc (20) % d = dcomp_rtm % rad_clear_sky_toc_ch20
         dcomp_input % rad_clear_sky_toa (20) % d = dcomp_rtm % rad_clear_sky_toa_ch20
         ! -- Solar irradiance in channel 20
         dcomp_input % solar_irradiance ( 20) = solar_ch20_nu
          
      end if
      
      do i = 21, 45
         if ( dcomp_input % is_channel_on (i)) dcomp_input % rad ( i ) % d = ch(i)%rad_toa
      end do
      
      
      
      if ( run_nlcomp) then
         if ( dcomp_input % is_channel_on (44)) then
            dcomp_input % refl ( 44 ) % d = ch(44)%ref_lunar_toa
            dcomp_input % alb_sfc ( 44 ) % d = ch(1)%sfc_ref_white_sky
         end if   
         dcomp_input % zen_lunar % d = geo % lunzen
         dcomp_input % azi_lunar % d = geo % lunrelaz
      
      end if
      
      
      dcomp_input % sat % d = geo % satzen
      dcomp_input % sol % d = geo % solzen
      dcomp_input % azi % d = geo % relaz
      
   
              
         ! - Cloud products
      dcomp_input % cloud_press % d = acha % pc
      dcomp_input % cloud_temp % d  = acha % tc
      dcomp_input % tau_acha % d    = acha % tau
      dcomp_input % cloud_mask % d  = cld_mask
!ccm
      dcomp_input % cloud_type % d  = cld_type(1:dim_1,1:dim_2)
!end ccm
!ccm      dcomp_input % cloud_type % d  = cld_type
                  
         ! - Flags
      dcomp_input % is_land % d = sfc % land_mask == 1 
      dcomp_input % is_valid % d = bad_pixel_mask /= 1
      
           
      dcomp_input % press_sfc % d =  dcomp_rtm % sfc_nwp
      dcomp_input % snow_class % d = sfc % snow
            
         ! - Atmospheric contents
         ! ozone column in Dobson
      dcomp_input % ozone_nwp % d = dcomp_rtm % ozone_path
         ! Total water Vapour above the cloud
      dcomp_input % tpw_ac % d = dcomp_rtm % tpw_ac
           


      
      call dcomp_rtm % deallocate_it()
      
      ! === THE MAIN CALL of DCOMP ===  
      
      if (run_nlcomp) then
         call nlcomp_array_loop_sub (dcomp_input,dncomp_output )
      else  
         debug_mode = 0
         call dcomp_input % check_input (debug_mode)
         call dcomp_array_loop ( dcomp_input , dncomp_output , debug_mode_user = debug_mode)
      end if
      ! === POPULATE CLAVR-X VARIABLES FROM PIXEL_COMMON
!ccm


            ! === POPULATE CLAVR-X VARIABLES FROM PIXEL_COMMON
      
      if ( run_nlcomp) then
      
         tau_nlcomp (1:dim_1,1:dim_2)   = dncomp_output % cod % d
         reff_nlcomp  (1:dim_1,1:dim_2) = dncomp_output % cps % d
         tau_nlcomp_cost(1:dim_1,1:dim_2) = dncomp_output % cod_unc % d
         reff_nlcomp_cost(1:dim_1,1:dim_2) = dncomp_output % ref_unc % d
         nlcomp_quality_flag(1:dim_1,1:dim_2) = dncomp_output %  quality % d
         nlcomp_info_flag(1:dim_1,1:dim_2) = dncomp_output % info % d
      
      else
      
         tau_dcomp (1:dim_1,1:dim_2)   = dncomp_output % cod % d(1:dim_1,1:dim_2)
         reff_dcomp  (1:dim_1,1:dim_2) = dncomp_output % cps % d(1:dim_1,1:dim_2)
         lwp_dcomp (1:dim_1,1:dim_2)   = dncomp_output % lwp % d(1:dim_1,1:dim_2)
         iwp_dcomp (1:dim_1,1:dim_2)   = dncomp_output % iwp % d(1:dim_1,1:dim_2)
      
         tau_dcomp_cost(1:dim_1,1:dim_2)     = dncomp_output % cod_unc % d(1:dim_1,1:dim_2)
         reff_dcomp_cost(1:dim_1,1:dim_2)    = dncomp_output % ref_unc % d(1:dim_1,1:dim_2)
         dcomp_quality_flag(1:dim_1,1:dim_2) = dncomp_output % quality % d(1:dim_1,1:dim_2)
         dcomp_info_flag(1:dim_1,1:dim_2)    = dncomp_output % info % d(1:dim_1,1:dim_2)
      
         cloud_063um_transmission_solar(1:dim_1,1:dim_2) = dncomp_output % cld_trn_sol % d(1:dim_1,1:dim_2)
         cloud_063um_transmission_view(1:dim_1,1:dim_2)  = dncomp_output % cld_trn_obs % d(1:dim_1,1:dim_2)
         cloud_063um_albedo(1:dim_1,1:dim_2)             = dncomp_output % cld_alb % d(1:dim_1,1:dim_2)
         cloud_063um_spherical_albedo(1:dim_1,1:dim_2)   = dncomp_output % cld_sph_alb % d(1:dim_1,1:dim_2)
      
         DCOMP_RELEASE_VERSION = dncomp_output % version
      end if
      
      

      

   end subroutine awg_cloud_dncomp_algorithm

  
   !---------------------------------------------------------------------------
   ! routine to set the cvs version in a global variable to write to hdf file
   !---------------------------------------------------------------------------
   subroutine set_dcomp_version()
      dcomp_version = DCOMP_RELEASE_VERSION
   end subroutine set_dcomp_version

    
end module dcomp_clavrx_bridge_mod 
