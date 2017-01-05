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
module dncomp_clavrx_bridge_mod


   ! -- MODULES USED

   use cx_constants_mod, only: &
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
       , tau_dcomp_1, tau_dcomp_2, tau_dcomp_3 &
       , reff_dcomp_1, reff_dcomp_2, reff_dcomp_3 &
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
      
   use calibration_constants, only: &
        sun_earth_distance  &      !---- check, This is defined in three routines   
      , solar_ch20_nu
     
   use dcomp_rtm_module
   
   use cx_array_tools_mod, only: &
      cx_rebin
          
   implicit none

   private
   
 
   public :: set_dcomp_version              
   public :: awg_cloud_dncomp_algorithm
   public :: awg_cloud_dncomp_algorithm_iband
   
   
   logical :: first_call = .true.
   logical :: first_call_iband = .true.
   
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
      
      logical :: run_nlcomp
      
      integer, allocatable :: possible_channels ( : )
      logical :: chan_on ( N_CHN ) = .false.
      integer :: i, i_mode
      integer :: CHN_VIS
      integer :: CHN_NIR
      
      integer :: dcomp_mode_local
      
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
         if ( sensor % chan_on_flag_default ( possible_channels ( i) )  ) then
            chan_on (possible_channels ( i)  )  = .true.
         end if
      end do 
      
     
      ! - all dcomp modes 
      do i_mode = 1, 3 
         dcomp_mode_local = i_mode
         if ( dcomp_mode .ne. i_mode .and. dcomp_mode .ne. 9) cycle
         ! - compute DCOMP related RTM 
      
         !-allocate input
         dcomp_input = dncomp_in_type ( dim_1, dim_2, chan_on )
      
         if ( .not. run_nlcomp) then      
            !- check mode
            CHN_VIS = 1
            select case ( dcomp_mode_local )
               case ( 1 ) 
                  CHN_NIR = 6               
               case ( 2 )
                  CHN_NIR = 7
               case ( 3 )
                  CHN_NIR = 20
               case default
                  print*, 'dcomp mode ',dcomp_mode_local,' not possible'
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
         dcomp_input % mode = dcomp_mode_local
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
         dcomp_input % is_valid % d = .NOT. bad_pixel_mask 
      
           
         dcomp_input % press_sfc % d =  dcomp_rtm % sfc_nwp
         dcomp_input % snow_class % d = sfc % snow
            
         ! - Atmospheric contents
         ! ozone column in Dobson
         dcomp_input % ozone_nwp % d = dcomp_rtm % ozone_path
         ! Total water Vapour above the cloud
         dcomp_input % tpw_ac % d = dcomp_rtm % tpw_ac
  
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
         
            select case (dcomp_mode_local)
         
            case (1)
               tau_dcomp_1 (1:dim_1,1:dim_2)   = dncomp_output % cod % d(1:dim_1,1:dim_2)
               reff_dcomp_1  (1:dim_1,1:dim_2) = dncomp_output % cps % d(1:dim_1,1:dim_2)
            case(2)
         
               tau_dcomp_2 (1:dim_1,1:dim_2)   = dncomp_output % cod % d(1:dim_1,1:dim_2)
               reff_dcomp_2  (1:dim_1,1:dim_2) = dncomp_output % cps % d(1:dim_1,1:dim_2)
         
            case(3)
               tau_dcomp_3 (1:dim_1,1:dim_2)   = dncomp_output % cod % d(1:dim_1,1:dim_2)
               reff_dcomp_3  (1:dim_1,1:dim_2) = dncomp_output % cps % d(1:dim_1,1:dim_2)
            end select
         
         
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
      
      end do
      call dcomp_rtm % deallocate_it()
      

   end subroutine awg_cloud_dncomp_algorithm
   
   
   
   !----------------------------------------------------------------------
   !  AWG_CLOUD_DCOMP_ALGORITHM
   !    This is the DCOMP bridge from CLAVR-x
   !   this is for iband VIIRS 0.6/1.6 algorithm only
   !---------------------------------------------------------------------- 
   subroutine awg_cloud_dncomp_algorithm_iband (  infile, path, algorithm_started )   
      
      use pixel_common, only: &
         Ref_chi1 &
         , ref_chi3   
       
      implicit none
 
      !--- input
		character(len=*) , intent(in) :: infile
      character(len=*) , intent(in) :: path
      ! - output 
      logical , intent(out) :: algorithm_started
      
      type(dcomp_rtm_type) :: dcomp_rtm
      type(dncomp_in_type)  :: dcomp_input
      type(dncomp_out_type) :: dncomp_output
      
      integer :: debug_mode      
      integer :: dim_1, dim_2

      logical :: chan_on ( N_CHN ) = .false.
   
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
      
      
            
      algorithm_started = .false.
      
      ! - do we need to run dcomp at all? ( night  etc..)
    
          
      if ( count ( geo % solzen < 75. .and. geo % solzen >= 0 .and. geo % satzen < 75. ) < 1 ) return
     
      algorithm_started = .true.
      
      
      if ( first_call_iband) then
        call mesg ('DCNOMP Iband starts ',color = 46 ) 
        first_call_iband = .false.
      end if
      
       
      ! - compute DCOMP related RTM 
      call perform_rtm_dcomp ( dcomp_rtm ) 
        
      dim_1 =  2 * Image%Number_Of_Elements
      dim_2 =  2 * Image%Number_Of_Lines_Read_This_Segment
   
      chan_on = .false.
      chan_on ( 1  )  = .true.
      chan_on ( 6  )  = .true.
      

      !-allocate input
      dcomp_input = dncomp_in_type ( dim_1, dim_2, chan_on )
      
      
  
      ! == CONFIGURE 
       
      ! - dcomp-mode
      dcomp_input % mode = 1
      ! - ancil/lut path
      dcomp_input % lut_path = trim(ancil_data_dir)//"/static/luts/cld/"   
      ! - wmo sensor id
      dcomp_input % sensor_wmo_id = sensor % wmo_id
      dcomp_input % sun_earth_dist = sun_earth_distance
            
      ! -  Satellite Data
      dcomp_input % refl ( 1 ) % d = Ref_chi1
      dcomp_input % alb_sfc ( 1 ) % d = cx_rebin(ch(1) % sfc_ref_white_sky(:,1:Image%Number_Of_Lines_Read_This_Segment), dim_1, dim_2)
      dcomp_input % gas_coeff(1) % d = solar_rtm % tau_h2o_coef(1,:)
               
      dcomp_input % refl ( 6 ) % d = ref_chi3
      dcomp_input % alb_sfc ( 6 ) % d = cx_rebin(ch(6) % sfc_ref_white_sky(:,1:Image%Number_Of_Lines_Read_This_Segment), dim_1, dim_2)
      dcomp_input % gas_coeff( 6 ) % d = solar_rtm % tau_h2o_coef(6,:)
      

      dcomp_input % sat % d = cx_rebin(geo % satzen, dim_1, dim_2)
      dcomp_input % sol % d = cx_rebin(geo % solzen, dim_1, dim_2)
      dcomp_input % azi % d = cx_rebin(geo % relaz, dim_1, dim_2)
          
      ! - Cloud products
      dcomp_input % cloud_press % d = cx_rebin(acha % pc, dim_1, dim_2)
      dcomp_input % cloud_temp % d  = cx_rebin(acha % tc, dim_1, dim_2)
      dcomp_input % tau_acha % d    = cx_rebin(acha % tau, dim_1, dim_2)
      dcomp_input % cloud_mask % d  = cx_rebin(cld_mask, dim_1, dim_2)
!ccm
      dcomp_input % cloud_type % d  = cx_rebin(cld_type, dim_1, dim_2)
!end ccm
!ccm      dcomp_input % cloud_type % d  = cld_type
                
         ! - Flags
      dcomp_input % is_land % d = cx_rebin(sfc % land_mask, dim_1, dim_2) == 1 
      dcomp_input % is_valid % d =cx_rebin( bad_pixel_mask, dim_1, dim_2) /= 1
     
       
      dcomp_input % press_sfc % d = cx_rebin( dcomp_rtm % sfc_nwp, dim_1, dim_2)
      dcomp_input % snow_class % d = sfc % snow
        
      ! - Atmospheric contents
      ! ozone column in Dobson
      dcomp_input % ozone_nwp % d = cx_rebin(dcomp_rtm % ozone_path, dim_1, dim_2)
      ! Total water Vapour above the cloud
      dcomp_input % tpw_ac % d = cx_rebin(dcomp_rtm % tpw_ac, dim_1, dim_2)
  
      ! === THE MAIN CALL of DCOMP ===  
     
        debug_mode = 0
      !   call dcomp_input % check_input (debug_mode)
        
      call dcomp_array_loop ( dcomp_input , dncomp_output , debug_mode_user = debug_mode)
     
      ! === POPULATE CLAVR-X VARIABLES FROM PIXEL_COMMON
!ccm


            ! === POPULATE CLAVR-X VARIABLES FROM PIXEL_COMMON
      
      
          
        
    
      
      call dcomp_rtm % deallocate_it()
      
      call add_to_file ( dncomp_output % cod % d, dncomp_output % cps % d , infile, path)
      

      
      

   end subroutine awg_cloud_dncomp_algorithm_iband
   
   
   subroutine add_to_file (product,prd2, file,path)
   
        use cx_hdf_write_mod, only:  &
      hdf_file_open &
      , create_sds &
      , compress_sds &
      , write_sds &
      , add_att &
      , close_sds &
      , close_file
      
      real, intent(in) :: product ( :,:)
      real, intent(in) :: prd2 ( :,:)
      logical :: first_seg = .true.
      character ( len=*), intent(in) :: file 
		character ( len=*), intent(in) ::path
      character (len = 240) :: outfile
      integer,save :: id_file
      integer, save :: sds_id, sds_id2
      
      integer :: sds_start_2d(2)
      integer :: sds_stride_2d(2) = [1,1]
      integer :: sds_edge_2d(2)
      integer :: istatus
      
		
      outfile = 'IBAND_LEVEL2_'//trim(file)//'.hdf'
      if ( first_seg ) then 
         id_file = hdf_file_open(trim(path)//trim(outfile), create=.true.)
         Sds_Id= create_sds (id_file, 'COD' , [6400,1536] , 4)
         Sds_Id2= create_sds (id_file, 'CPS' , [6400,1536] , 4)
         sds_start_2d = [0,0]
         first_seg = .false.
      end if
      
      
      sds_edge_2d = shape ( product)
       
      istatus = write_sds ( sds_id, sds_start_2d, sds_stride_2d , sds_edge_2d,  product)
      istatus = write_sds ( sds_id2, sds_start_2d, sds_stride_2d , sds_edge_2d,  prd2)
     
      
       if ( istatus /= 0 ) then
         print*,'something wrong with write sds ', istatus, sds_id
         print*,sds_id, sds_start_2d, sds_stride_2d , sds_edge_2d
         print*,'dimensions array: ',shape(product)
         stop
      end if
      sds_start_2d = [0,sds_start_2d(2)+100]
      
      if (sds_edge_2d(2) .lt. 50) then
         call close_sds (  sds_id)
         call close_sds (  sds_id2)
         call close_file (id_file)
			first_seg = .true.
     end if
   end subroutine 
   
  
  
   !---------------------------------------------------------------------------
   ! routine to set the cvs version in a global variable to write to hdf file
   !---------------------------------------------------------------------------
   subroutine set_dcomp_version()
      dcomp_version = DCOMP_RELEASE_VERSION
   end subroutine set_dcomp_version

    
end module dncomp_clavrx_bridge_mod 

