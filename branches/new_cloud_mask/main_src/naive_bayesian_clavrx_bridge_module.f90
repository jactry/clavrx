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
!
!--------------------------------------------------------------------------------------
module naive_bayesian_clavrx_bridge_module


   ! -- MODULES USED

   use CONSTANTS , only: &
      & sym
      
   use PIXEL_COMMON
       
   use CLOUD_MASK_MOD, only : &
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
      
      character (len=120):: Naive_Bayes_File_Name
      type ( cloud_mask_input_type ) :: inp
      integer :: info_flags(7)
      integer :: i , j 
      real :: erg
           
    ! -----------      
      line_loop: do i = 1, num_pix
         elem_loop: do  j = 1,num_scans_read
            
            if ( land (i,j) < 0 ) cycle
            
            
            inp % bayesian_mask_classifier = trim(Ancil_Data_Dir)//'/naive_bayes_mask/'//trim(Bayesian_Cloud_Mask_Name) 
            
            inp % geo % lat         = lat ( i,j)            
            inp % geo % lon         = lon ( i,j)
            inp % geo % sol_zen     = solzen ( i,j)
            inp % geo % airmass     = airmass ( i,j)
            inp % geo % scat_angle  = scatangle ( i,j)
            inp % geo % glint       = glint_mask ( i ,j )
   
            inp % sfc % land_class  = land ( i , j )
            inp % sfc % coast_mask  = coast ( i , j ) 
            inp % sfc % snow_class  = snow( i , j )
            inp % sfc % dem         = zsfc  ( i , j )
         
            inp % rtm % bt_ch31_lrc     =  Bt_Ch31_LRC ( i , j )
            inp % rtm % bt_ch31_3x3_max = Bt_Ch31_Max_3x3 ( i , j )
            inp % rtm % bt_ch31_3x3_std = Bt_Ch31_Std_3x3 ( i , j )
         
        
            inp % rtm % emis_ch31_tropo = ch(31) % emiss_tropo (i,j)
            inp % rtm % emis_ch32_tropo = ch(31) % emiss_tropo (i,j)
            inp % rtm % emis_ch20_clear = Ems_Ch20_Clear_Solar_Rtm( i , j ) 
            inp % rtm % bt_ch31_atm_sfc = ch(31)%Bt_Toa_Clear( i , j ) 
            inp % rtm % bt_ch32_atm_sfc = ch(32)%Bt_Toa_Clear( i , j ) 
            inp % rtm % ref_ch1_clear   = ch(1) % Ref_toa_clear (i,j)
   
            inp % sat % bt_ch20             = ch(20) % bt_toa ( i , j )
            !inp % sat % bt_ch27             = sat % chn(27) % bt ( i , j )
            inp % sat % bt_ch29             = ch(29) % bt_toa ( i , j )
            inp % sat % bt_ch31             = ch(31) % bt_toa ( i , j )
            inp % sat % bt_ch32             = ch(32) % bt_toa ( i , j )

            inp % sat % ref_ch1             = ch(1) % ref_toa ( i , j )
        
            inp % sat % ref_ch2             = ch(2) % ref_toa ( i , j )
        
            inp % sat % ref_ch6             = ch(6) % ref_toa ( i , j )
         
            inp % sat % ref_ch26            = ch(26) % ref_toa ( i , j )
        
            inp % sat % emis_ch20_3x3_mean  = ems_ch20_median_3x3 ( i , j )
        
            inp % sat % ref_ch1_3x3_std     = Ref_Ch1_Std_3x3(i,j)
            inp % sat % ref_ch1_3x3_min     = Ref_Ch1_Min_3x3(i,j)
                       
            inp % sat % chan_on ( 1:40) = Chan_On_Flag_Default(1:40) == 1
             
            call cloud_mask_naive_bayes ( inp, erg , info_flags )
         
            Posterior_Cld_Probability ( i , j ) = erg
            Bayes_Mask_Sfc_Type_Global (  i , j) = ibits ( info_flags (3) , 0, 3 ) 
            Cld_Test_Vector_Packed ( : , i , j ) = info_flags
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
