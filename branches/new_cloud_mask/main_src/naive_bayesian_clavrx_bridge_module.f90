!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: naive_bayesian_clavrx_bridge_module.f90 (src)
!       AWG_CLOUD_BAYES_BRIDGE (program)
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
!  HISTORY:
!      2014/04/06:    new interface (AW)
!      2014/04/21:  add fire mask input
!
!  GLOBAL VARIABLES:
!
!    FROM PIXEL_COMMON:
!      1. work as input
!        1.1 configuration
!             Chan_On_Flag_Default                       integer (42)
!             ancil_data_dir                             character
!             Bayesian_Cloud_Mask_Name                   character
!             num_pix                                    integer
!             num_scans_read                             integer
!        1.2 geo data:
!              lon                                       real (:,:)
!              lat                                       real (:,:)                        
!              solzen                                    real (:,:)
!              scatangle                                 real (:,:)
!              airmass                                   real (:,:)
!              glint_mask                                integer (:,:)
!        1.3 surface 
!              land                                      integer (:,:)
!              coast                                     integer (:,:)
!              snow                                      integer (:,:)
!              zsfc                                      real (:,:)
!        1.4 rtm / statistics 
!              bt_ch31_lrc                               real (:,:)
!              bt_ch31_max_3x3                           real (:,:)
!              bt_Ch31_Std_3x3                           real (:,:)
!              bt_Ch20_Std_3x3                           real (:,:)
!              Ems_Ch20_Clear_Solar_Rtm                  real (:,:)
!              Covar_Ch27_Ch31_5x5                       real (:,:)
!              ems_ch20_median_3x3                       real (:,:)
!        1.5 observations
!              ch                                        type ( observations )  
!              Ref_Ch1_Std_3x3                           real (:,:)
!              Ref_Ch1_Min_3x3                           real (:,:)
!       
!         2. work as output  
!              Posterior_Cld_Probability                 real (:,:) 
!              Bayes_Mask_Sfc_Type_Global                integer ( 3,:,:)
!              Cld_Test_Vector_Packed                    integer ( 7,:,:)
!              cld_mask                                  integer (:,:) 
!           
!     FROM CONSTANTS:
!            sym                                         type ( symbol_struct ) 
!
!     NAIVE_BAYESIAN_CLOUD_MASK_MODULE
!             cloud_mask_naive_bayes                     subroutine
!                  subroutine which retrieves cloud probability
!             cloud_mask_input_type                      type definition
!                   structure definition of cloud_mask_naive_bayes retrieval input  
!
!--------------------------------------------------------------------------------------
module naive_bayesian_clavrx_bridge_module


   ! -- MODULES USED

   use CONSTANTS , only: &
      & sym
      
   use PIXEL_COMMON, only: &
        lon &
      , lat &
      , solzen &
      , airmass &
      , scatangle &
      , glint_mask &
      , land &
      , coast &
      , snow &
      , zsfc &
      , bt_ch31_lrc &
      , bt_ch31_max_3x3 &
      , Bt_Ch31_Std_3x3 &
      , Bt_Ch20_Std_3x3 &
      , ems_Ch20_Clear_Solar_Rtm &
      , ems_ch20_median_3x3 &
      , Covar_Ch27_Ch31_5x5 &
      , ch &
      , Ref_Ch1_Std_3x3 &
      , Ref_Ch1_Min_3x3 &
      , Chan_On_Flag_Default &
      , Posterior_Cld_Probability &
      , Bayes_Mask_Sfc_Type_Global &
      , Cld_Test_Vector_Packed &
      , num_pix &
      , num_scans_read &
      , cld_mask &
      , ancil_data_dir &
      , Bayesian_Cloud_Mask_Name 
     
      
             
   use NAIVE_BAYESIAN_CLOUD_MASK_MODULE , only : &
      &   cloud_mask_naive_bayes &
      & , cloud_mask_input_type


   implicit none

   public :: AWG_CLOUD_BAYES_BRIDGE
   

contains
   !----------------------------------------------------------------------
   !
   !---------------------------------------------------------------------- 
   subroutine AWG_CLOUD_BAYES_BRIDGE()
 
      implicit none
      
     
      type ( cloud_mask_input_type ) :: mask_inp
      ! cloud mask information mask 7 bytes ( 56 bits)    
      integer :: info_flags ( 7 )
      integer :: i , j 
                 
      ! -----------    loop over pixels -----   
      line_loop: do i = 1, num_pix
         elem_loop: do  j = 1,num_scans_read
            
            if ( land (i,j) < 0 ) cycle
            
            
            mask_inp % bayesian_mask_classifier = trim(Ancil_Data_Dir)//'/naive_bayes_mask/'//trim(Bayesian_Cloud_Mask_Name) 
            
            mask_inp % geo % lat         = lat ( i,j)            
            mask_inp % geo % lon         = lon ( i,j)
            mask_inp % geo % sol_zen     = solzen ( i,j)
            mask_inp % geo % airmass     = airmass ( i,j)
            mask_inp % geo % scat_angle  = scatangle ( i,j)
            mask_inp % geo % glint       = glint_mask ( i ,j )
   
            mask_inp % sfc % land_class  = land ( i , j )
            mask_inp % sfc % coast_mask  = coast ( i , j ) 
            mask_inp % sfc % snow_class  = snow( i , j )
            mask_inp % sfc % dem         = zsfc  ( i , j )
         
            mask_inp % rtm % bt_ch31_lrc     =  Bt_Ch31_LRC ( i , j )
            mask_inp % rtm % bt_ch31_3x3_max = Bt_Ch31_Max_3x3 ( i , j )
            mask_inp % rtm % bt_ch31_3x3_std = Bt_Ch31_Std_3x3 ( i , j )
            mask_inp % rtm % bt_ch20_3x3_std = Bt_Ch20_Std_3x3( i , j )
        
            mask_inp % rtm % emis_ch31_tropo    = ch(31) % emiss_tropo (i,j)
            mask_inp % rtm % emis_ch32_tropo    = ch(32) % emiss_tropo (i,j)
            mask_inp % rtm % emis_ch20_clear    = Ems_Ch20_Clear_Solar_Rtm( i , j ) 
            mask_inp % rtm % bt_ch31_atm_sfc    = ch(31)%Bt_Toa_Clear( i , j ) 
            mask_inp % rtm % bt_ch32_atm_sfc    = ch(32)%Bt_Toa_Clear( i , j ) 
            mask_inp % rtm % ref_ch1_clear      = ch(1) % Ref_toa_clear (i,j)
				if ( chan_on_flag_default(27) == 1 ) then
            	mask_inp % rtm % bt_ch31_ch27_covar = Covar_Ch27_Ch31_5x5 ( i , j )
   				mask_inp % sat % bt_ch27             = ch(27) % bt_toa ( i , j )
				end if	
            mask_inp % sat % bt_ch20             = ch(20) % bt_toa ( i , j )
            
            mask_inp % sat % bt_ch29             = ch(29) % bt_toa ( i , j )
            mask_inp % sat % bt_ch31             = ch(31) % bt_toa ( i , j )
            mask_inp % sat % bt_ch32             = ch(32) % bt_toa ( i , j )

            mask_inp % sat % ref_ch1             = ch(1) % ref_toa ( i , j )
        
            mask_inp % sat % ref_ch2             = ch(2) % ref_toa ( i , j )
        
            mask_inp % sat % ref_ch6             = ch(6) % ref_toa ( i , j )
            mask_inp % sat % ref_ch7             = ch(7) % ref_toa ( i , j )
            mask_inp % sat % ref_ch8             = ch(8) % ref_toa ( i , j )
         
            mask_inp % sat % ref_ch26            = ch(26) % ref_toa ( i , j )
        
            mask_inp % sat % emis_ch20_3x3_mean  = ems_ch20_median_3x3 ( i , j )
        
            mask_inp % sat % ref_ch1_3x3_std     = Ref_Ch1_Std_3x3(i,j)
            mask_inp % sat % ref_ch1_3x3_min     = Ref_Ch1_Min_3x3(i,j)
                       
            mask_inp % sat % chan_on             = Chan_On_Flag_Default == 1
            
            
           
            call cloud_mask_naive_bayes ( mask_inp, Posterior_Cld_Probability ( i , j ) , info_flags )
         
              
            Bayes_Mask_Sfc_Type_Global (  i , j ) = ibits ( info_flags (3) , 0, 3 ) 
            Cld_Test_Vector_Packed ( : , i , j )  = info_flags
            
         end do elem_loop
      end do   line_loop
      
       
      
      !------------------------------------------------------------------------------------------------------------
      !--- make a cloud mask
      !------------------------------------------------------------------------------------------------------------
      Cld_Mask(:,:) = sym%CLEAR
        
      where ( Posterior_Cld_Probability >= 0.9 )
         cld_mask = sym % CLOUDY
      end where
      
      where ( Posterior_Cld_Probability >= 0.5 .and. Posterior_Cld_Probability < 0.9 )
         cld_mask = sym % PROB_CLOUDY
      end where
        
      where ( Posterior_Cld_Probability > 0.1 .and. Posterior_Cld_Probability < 0.5 )
         cld_mask = sym % PROB_CLEAR
      end where 


   end subroutine AWG_CLOUD_BAYES_BRIDGE




end module naive_bayesian_clavrx_bridge_module
