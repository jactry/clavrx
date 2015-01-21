! $Header:$
! test cloud mask
!
!
program test_cloud_mask


use nb_cloud_mask, only : &
        NB_CLOUD_MASK_ALGORITHM &
      , cloud_mask_input_type &
!      , dust_detection &
!      , fire_detection &
      , cloud_mask_diagnostic &
      , Cloud_Mask_Version_Type
      
      
type ( cloud_mask_input_type ) :: inp
integer :: info_flags(7)
real :: erg
integer :: counter
type (cloud_mask_diagnostic)  :: diag 
type (Cloud_Mask_Version_Type) :: vers

 
  


!print*,'DUST detection: ',dust_detection(-12.,12.,1.,.true.,1)

!print*,'FIRE detection: ',fire_detection ( 310.,312., 2., 1. ,77. )


!inp % bayesian_mask_classifier = &
!     '/DATA/Ancil_Data/clavrx_ancil_data/naive_bayes_mask/viirs_default_bayes_mask.nc'

inp % bayesian_mask_classifier = &
     '/data/Ancil_Data/clavrx_ancil_data/static/luts/nb_cloud_mask/viirs_default_nb_cloud_mask_lut.nc'
            
            inp % geo % lat         = 60.            
            inp % geo % lon         = 20.
            inp % geo % sol_zen     = 30.
            inp % geo % airmass     = 3.4
            inp % geo % scat_angle  = 134.
            inp % geo % glint       = 0
   
            inp % sfc % land_class  = 1
            inp % sfc % coast_mask  = 0
            inp % sfc % snow_class  = 0
            inp % sfc % dem         = 4
         
            inp % rtm % bt_ch31_lrc     =  263.4831
            inp % rtm % bt_ch31_3x3_max = 266.6003
            inp % rtm % bt_ch31_3x3_std = 1.556237
         
        
            inp % rtm % emis_ch31_tropo = 0.287442
            inp % rtm % emis_ch32_tropo = 0.287442
            inp % rtm % emis_ch20_clear = 1.073115
            inp % rtm % bt_ch31_atm_sfc = 276.2816
            inp % rtm % bt_ch32_atm_sfc = 276.1228
            inp % rtm % ref_ch1_clear   = 13.93634
   
            inp % sat % bt_ch20             = 292.84
           
            inp % sat % bt_ch29             = 290.3945
            inp % sat % bt_ch31             = 293.4831
            inp % sat % bt_ch32             = 290.9036

            inp % sat % ref_ch1             = 6.07673
        
            inp % sat % ref_ch2             = 5.67140
        
            inp % sat % ref_ch6             = 1.49225
         
            inp % sat % ref_ch26            = 1.748
        
            inp % sat % emis_ch20_3x3_mean  = 1.970
        
            inp % sat % ref_ch1_3x3_std     = 1.425
            inp % sat % ref_ch1_3x3_min     = 1.22520
                       
            inp % sat % chan_on  = .true.
                         
            call NB_CLOUD_MASK_ALGORITHM ( inp, erg , info_flags , diag , vers )      
            print*,'cloud probability: ',erg
            print*,'info flags: ', info_flags
            counter = 0
            if ( btest(info_Flags(3),4) ) counter = counter + 1
            if ( btest(info_Flags(3),5) ) counter = counter + 2
            
            print*,'test t11: ',counter

            counter = 0
            if ( btest(info_Flags(5),4) ) counter = counter + 1
            if ( btest(info_Flags(5),5) ) counter = counter + 2
            
            print*,'test emis 037: ',counter

end program test_cloud_mask
