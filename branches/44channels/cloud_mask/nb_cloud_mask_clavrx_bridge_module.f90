!$Id: naive_bayesian_clavrx_bridge_module.f90 472 2014-06-30 14:54:15Z awalther $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: naive_bayesian_clavrx_bridge_module.f90 (src)
!       CLOUD_MASK_NAIVE_BAYES_BRIDGE (program)
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
!      2014/04/21:  add fire mask input (AW)
!      2014/05/07:  added diagnostic type (Denis B)
!      2014/05/09:  added version (Denis B)
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
!              lunzen                                    real (:,:)
!              scatangle                                 real (:,:)
!              airmass                                   real (:,:)
!              glint_mask                                integer (:,:)       
!              Solar_Contamination_Mask                  integer (:,:)
!              scatangle_lunar                           real (:,:)             ONLY VIIRS
!        1.3 surface 
!              land                                      integer (:,:)
!              coast                                     integer (:,:)
!              snow                                      integer (:,:)
!              zsfc                                      real (:,:)
!              city_mask                                 integer (:,:)
!        1.4 rtm / statistics 
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
!     FILE_TOOLS
!           file_test         
!
!     NAIVE_BAYESIAN_CLOUD_MASK_MODULE
!             cloud_mask_naive_bayes                     subroutine
!                  subroutine which retrieves cloud probability
!             cloud_mask_input_type                      type definition
!             cloud_mask_diagnostic                      type diagnostic
!                   structure definition of cloud_mask_naive_bayes retrieval input  
!
!--------------------------------------------------------------------------------------

module NB_CLOUD_MASK_CLAVRX_BRIDGE

   ! -- MODULES USED

   use CONSTANTS , only: &
        Sym &
      , Cloud_Mask_Version &
      , Cloud_Mask_Thresholds_Version
      
   use PIXEL_COMMON, only: &
         Nav &
      ,  Geo &
      , sfc &
      , sensor &
      , Sst_Anal_Uni &
      , Bt_Ch31_Max_3x3 &
      , Bt_Ch31_Std_3x3 &
      , Bt_Ch20_Std_3x3 &
      , Ems_Ch20_Clear_Solar_Rtm &
      , Ems_Ch20_Median_3x3 &
      , Covar_Ch27_Ch31_5x5 &
      , Ch &
      , Ref_Ch1_Std_3x3 &
      , Ref_Ch1_Min_3x3 &
      
      , Posterior_Cld_Probability &
      , Bayes_Mask_Sfc_Type_Global &
      , Cld_Test_Vector_Packed &
      , Ancil_Data_Dir &
      , Bayesian_Cloud_Mask_Name &
      , Solar_Contamination_Mask &
      , Diag_Pix_Array_1 &
      , Diag_Pix_Array_2 &
      , Diag_Pix_Array_3 &
      , Space_Mask &
      
      , Ref_ChDNB_Lunar_Std_3x3 &
      , Ref_ChDNB_Lunar_Min_3x3 &
     , cld_mask &
     
      , Ref_Min_ChI1 &
      , Ref_Max_ChI1 &
      , Ref_Mean_ChI1 &
      , Ref_Uni_ChI1 &
      , Ref_Min_ChI2 &
      , Ref_Max_ChI2 &
      , Ref_Mean_ChI2 &
      , Ref_Uni_ChI2 & 
      , Ref_Min_ChI3 &
      , Ref_Max_ChI3 &
      , Ref_Mean_ChI3 &
      , Ref_Uni_ChI3 &
      , Bt_Min_ChI4 &
      , Bt_Max_ChI4 &
      , Bt_Mean_ChI4 &
      , Bt_Uni_ChI4 &
      , Bt_Min_ChI5 &
      , Bt_Max_ChI5 &
      , Bt_Mean_ChI5 &
      , Bt_Uni_ChI5
     
   use NB_CLOUD_MASK , only : &
        NB_CLOUD_MASK_ALGORITHM &
      , Cloud_Mask_Input_Type &
      , ET_Cloudiness_Class &
      , Cloud_Mask_Diagnostic &
      , Cloud_Mask_Version_Type

   use FILE_TOOLS, only: &
        FILE_TEST

   public :: NB_CLOUD_MASK_BRIDGE
   

contains
   !----------------------------------------------------------------------
   !
   !---------------------------------------------------------------------- 
   subroutine NB_CLOUD_MASK_BRIDGE(segment_number)
 
      implicit none
     
      integer, intent(in), optional:: segment_number 
     
      type ( Cloud_Mask_Input_Type ) :: mask_inp
      type ( Cloud_Mask_Diagnostic ) :: diag
      type ( Cloud_Mask_Version_Type ) :: vers
      ! cloud mask information mask 7 bytes ( 56 bits)    
      integer :: info_flags ( 7 )
      integer :: i , j 
      ! --- set cloud mask probability thresholds based on sfc type
      real, dimension (7) :: cld_mask_probab_thresh_lo = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
      real, dimension (7) :: cld_mask_probab_thresh_mi = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
      real, dimension (7) :: cld_mask_probab_thresh_hi = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
      real :: cld_mask_probab_thresh_lo_tmp, cld_mask_probab_thresh_mi_tmp, cld_mask_probab_thresh_hi_tmp
      
      mask_inp % bayesian_mask_classifier = trim(Ancil_Data_Dir)//'/naive_bayes_mask/'//trim(Bayesian_Cloud_Mask_Name) 
      
      if ( .not. file_test ( mask_inp % bayesian_mask_classifier ) ) then
         mask_inp % bayesian_mask_classifier = trim(Ancil_Data_Dir)//'/bayes/'//trim(Bayesian_Cloud_Mask_Name)
         if ( .not. file_test ( mask_inp % bayesian_mask_classifier ) ) then
            print*,'Classifier file not there: '
            print*, 'check location: '
            print*, mask_inp % bayesian_mask_classifier
            print*,'stopping.........'
            stop
         end if
      end if
                 
                 
      ! -----------    loop over pixels -----   
      line_loop: do i = 1, ubound ( nav % lat ,1 )
         elem_loop: do  j = 1,ubound ( nav % lat ,2 )
            
            if ( Space_Mask (i,j) == 1) cycle
            
            if ( Sfc % Land (i,j) < 0 ) cycle
            
            mask_inp % geo % lat         = nav % Lat ( i , j )            
            mask_inp % geo % lon         = nav % Lon ( i , j )
            mask_inp % geo % sol_zen     = geo % Solzen ( i , j )
            mask_inp % geo % airmass     = geo % Airmass ( i , j )
            mask_inp % geo % scat_angle  = Geo % Scatangle ( i , j )
            mask_inp % geo % glint       = sfc % Glint_Mask ( i ,j ) == 1
            
            mask_inp % geo % solar_conta = Solar_Contamination_Mask ( i , j ) == 1
   
            mask_inp % sfc % land_class  = sfc % Land ( i , j )
            mask_inp % sfc % coast_mask  = sfc % Coast_mask ( i , j ) == 1 
            mask_inp % sfc % snow_class  = sfc % Snow ( i , j )
            mask_inp % sfc % sfc_type    = sfc % Sfc_Type ( i , j ) 
            mask_inp % sfc % dem         = sfc % Zsfc ( i , j )
            mask_inp % sfc % sst_anal_uni = Sst_Anal_Uni ( i, j )
            
            if ( sensor % chan_on_flag_default(1) == 1 ) then
               mask_inp % rtm % ref_ch1_clear   = ch(1) % Ref_Toa_Clear ( i , j )
               mask_inp % sat % ref_ch1         = ch(1) % Ref_Toa ( i , j )
               mask_inp % sat % ref_ch1_3x3_std = Ref_Ch1_Std_3x3 ( i , j )
               mask_inp % sat % ref_ch1_3x3_min = Ref_Ch1_Min_3x3 ( i , j )
            end if
            
            if ( sensor % chan_on_flag_default(2) == 1 ) mask_inp % sat % ref_ch2 = Ch (2) % Ref_Toa ( i , j )
            if ( sensor % chan_on_flag_default(6) == 1 ) mask_inp % sat % ref_ch6 = Ch (6) % Ref_Toa ( i , j )
            if ( sensor % chan_on_flag_default(7) == 1 ) mask_inp % sat % ref_ch7 = Ch (7) % Ref_Toa ( i , j )
            if ( sensor % chan_on_flag_default(8) == 1 ) mask_inp % sat % ref_ch8 = Ch (8) % Ref_Toa ( i , j )
            
            if ( sensor % chan_on_flag_default(20) == 1 ) then 
                mask_inp % rtm % bt_ch20_3x3_std = Bt_Ch20_Std_3x3( i , j )
                mask_inp % rtm % emis_ch20_clear = Ems_Ch20_Clear_Solar_Rtm( i , j )
                mask_inp % sat % bt_ch20         = Ch (20) % Bt_Toa ( i , j )
                mask_inp % sat % emis_ch20_3x3_mean  = Ems_Ch20_Median_3x3 ( i , j )
            end if
            
            ! -  sfc emissivity is always on 
            mask_inp % sfc % emis_ch20 =  ch(20) % sfc_emiss ( i , j )
            
            if ( sensor % chan_on_flag_default(26) == 1 ) mask_inp % sat % ref_ch26 = Ch (26) % Ref_Toa ( i , j )
            if ( sensor % chan_on_flag_default(27) == 1 ) mask_inp % sat % bt_ch27  = Ch (27) % Bt_Toa ( i , j )
            if ( sensor % chan_on_flag_default(29) == 1 ) mask_inp % sat % bt_ch29  = Ch (29) % Bt_Toa ( i , j )
            
            if ( sensor % chan_on_flag_default(31) == 1 ) then
               
               mask_inp % rtm % bt_ch31_3x3_max = Bt_Ch31_Max_3x3 ( i , j )
               mask_inp % rtm % bt_ch31_3x3_std = Bt_Ch31_Std_3x3 ( i , j )
               mask_inp % rtm % emis_ch31_tropo = Ch (31) % Emiss_Tropo ( i , j )
               mask_inp % rtm % bt_ch31_atm_sfc = Ch (31) % Bt_Toa_Clear( i , j )
               mask_inp % sat % bt_ch31         = Ch (31) % Bt_Toa ( i , j )
               if ( sensor % chan_on_flag_default(27) == 1 ) then
                  mask_inp % rtm % bt_ch31_ch27_covar = Covar_Ch27_Ch31_5x5 ( i , j )
               end if   
            end if   
               
            if (sensor %  chan_on_flag_default(32) == 1 ) then
               mask_inp % rtm % emis_ch32_tropo    = Ch (32) % Emiss_Tropo ( i , j )
               mask_inp % rtm % bt_ch32_atm_sfc    = Ch (32) % Bt_Toa_Clear( i , j ) 
               mask_inp % sat % bt_ch32            = Ch (32) % Bt_Toa ( i , j )
            end if
            
            ! - ibands of viirs
            
            if ( sensor % chan_on_flag_default(39) == 1 ) then
               mask_inp % sat % iband(1) % ref % min =   Ref_Min_ChI1 ( i , j )
               mask_inp % sat % iband(1) % ref % max =   Ref_Max_ChI1 ( i , j )
               mask_inp % sat % iband(1) % ref % mean =  Ref_Mean_ChI1 ( i , j )
               mask_inp % sat % iband(1) % ref % std =   Ref_Uni_ChI1 ( i , j )
            end if
            
            if ( sensor % chan_on_flag_default(40) == 1 ) then
               mask_inp % sat % iband(2) % ref % min =   Ref_Min_ChI2 ( i , j )
               mask_inp % sat % iband(2) % ref % max =   Ref_Max_ChI2 ( i , j )
               mask_inp % sat % iband(2) % ref % mean =  Ref_Mean_ChI2 ( i , j )
               mask_inp % sat % iband(2) % ref % std =   Ref_Uni_ChI2 ( i , j )
            end if

            if ( sensor %  chan_on_flag_default(41) == 1 ) then
               mask_inp % sat % iband(3) % ref % min =   Ref_Min_ChI3 ( i , j )
               mask_inp % sat % iband(3) % ref % max =   Ref_Max_ChI3 ( i , j )
               mask_inp % sat % iband(3) % ref % mean =  Ref_Mean_ChI3 ( i , j )
               mask_inp % sat % iband(3) % ref % std =   Ref_Uni_ChI3 ( i , j )
            end if

            if ( sensor % chan_on_flag_default(42) == 1 ) then
               mask_inp % sat % iband(4) % bt % min =   Bt_Min_ChI4 ( i , j )
               mask_inp % sat % iband(4) % bt % max =   Bt_Max_ChI4 ( i , j )
               mask_inp % sat % iband(4) % bt % mean =  Bt_Mean_ChI4 ( i , j )
               mask_inp % sat % iband(4) % bt % std =   Bt_Uni_ChI4 ( i , j )
            end if
            
            if ( sensor % chan_on_flag_default(43) == 1 ) then
               mask_inp % sat % iband(5) % bt % min =   Bt_Min_ChI5 ( i , j )
               mask_inp % sat % iband(5) % bt % max =   Bt_Max_ChI5 ( i , j )
               mask_inp % sat % iband(5) % bt % mean =  Bt_Mean_ChI5 ( i , j )
               mask_inp % sat % iband(5) % bt % std =   Bt_Uni_ChI5 ( i , j )
            end if
            
            
            
            ! - dnb cloud mask addition at night
            if ( sensor % chan_on_flag_default(44) == 1 .and. allocated( Ch (44) % Ref_Lunar_Toa ) ) then
               mask_inp % sat % ref_dnb_lunar    = Ch (44) % Ref_Lunar_Toa ( i , j )
               mask_inp % geo % lunar_zen        = geo % Lunzen ( i , j )
               mask_inp % rtm % ref_dnb_clear    = Ch (44) % Ref_Lunar_Toa_Clear( i , j )
               mask_inp % sfc % is_city          = sfc % City_Mask ( i , j ) == 1
               mask_inp % sat % ref_dnb_3x3_std  = Ref_ChDNB_Lunar_Std_3x3 ( i , j )
               mask_inp % sat % ref_dnb_3x3_min  = Ref_ChDNB_Lunar_Min_3x3 ( i , j )
               mask_inp % geo % scat_angle_lunar = geo % Scatangle_Lunar ( i , j )
               mask_inp % geo % lunar_glint_mask = sfc % Glint_Mask_Lunar ( i , j ) == 1 
            end if

            mask_inp % sat % chan_on             = sensor % Chan_On_Flag_Default == 1
             
            call NB_CLOUD_MASK_ALGORITHM ( mask_inp, Posterior_Cld_Probability ( i , j ) , info_flags &
                                         , diag , vers )
           
            Bayes_Mask_Sfc_Type_Global (  i , j ) = ibits ( info_flags (3) , 0, 3 ) 
            Cld_Test_Vector_Packed ( : , i , j )  = info_flags

            ! - save diagnostic pixels to global arrays (ONLY IF ACTIVELY DIAGNOSING)
!           Diag_Pix_Array_1 ( i , j ) = diag % diagnostic_1
!           Diag_Pix_Array_2 ( i , j ) = diag % diagnostic_2
!           Diag_Pix_Array_3 ( i , j ) = diag % diagnostic_3

            !--------------------------------------------------------------------------------------------
            !--- make a cloud mask
            !--------------------------------------------------------------------------------------------
            Cld_Mask ( i , j ) =  ET_cloudiness_class % CLEAR

            ! - based on type of srfc could be different thresholds
            if (Bayes_Mask_Sfc_Type_Global (  i , j ) > 0) then
               cld_mask_probab_thresh_lo_tmp = cld_mask_probab_thresh_lo (Bayes_Mask_Sfc_Type_Global (  i , j ))
               cld_mask_probab_thresh_mi_tmp = cld_mask_probab_thresh_mi (Bayes_Mask_Sfc_Type_Global (  i , j ))
               cld_mask_probab_thresh_hi_tmp = cld_mask_probab_thresh_hi (Bayes_Mask_Sfc_Type_Global (  i , j ))

               if ( Posterior_Cld_Probability ( i , j ) >= cld_mask_probab_thresh_hi_tmp ) then
                  Cld_Mask ( i , j ) =  ET_cloudiness_class %  CLOUDY
               end if

               if ( Posterior_Cld_Probability ( i , j ) >= cld_mask_probab_thresh_mi_tmp &
              .and. Posterior_Cld_Probability ( i , j ) < cld_mask_probab_thresh_hi_tmp ) then
                  Cld_Mask ( i , j ) =  ET_cloudiness_class %  PROB_CLOUDY
               end if

               if ( Posterior_Cld_Probability ( i , j ) > cld_mask_probab_thresh_lo_tmp &
              .and. Posterior_Cld_Probability ( i , j ) < cld_mask_probab_thresh_mi_tmp ) then
                  Cld_Mask ( i , j ) =  ET_cloudiness_class % PROB_CLEAR
               end if
            end if

            if ( Space_Mask ( i , j ) == 1) then
               Cld_Mask ( i , j ) = ET_cloudiness_class % SPACE
            end if

            if ( sfc % Land ( i , j ) < 0 .and. Space_Mask ( i , j ) /= 1) then
               Cld_Mask ( i , j ) = ET_cloudiness_class % MISSING
            end if

         end do elem_loop
      end do line_loop
      
      ! - save mask and threshold version to global variables
      Cloud_Mask_Version = vers % cloud_mask_version_id
      Cloud_Mask_Thresholds_Version = vers % cloud_mask_thresh_version_id
       
   end subroutine NB_CLOUD_MASK_BRIDGE


!------------------------------------------------------------------------------------------------------------

end module NB_CLOUD_MASK_CLAVRX_BRIDGE

