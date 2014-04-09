! $Header:$
! test cloud mask
!
!
program test_cloud_mask


use naive_bayesian_cloud_mask_module, only : &
      &   cloud_mask_naive_bayes &
      & , cloud_mask_input_type &
      & , dust_viirs
      
      
type ( cloud_mask_input_type ) :: inp
integer :: info_flags(7)
real :: erg
integer :: counter


print*,dust_viirs(-12.,12.,1.,1,1)


inp % bayesian_mask_classifier = &
   &  '/DATA/Ancil_Data/clavrx_ancil_data/naive_bayes_mask/viirs_default_bayes_mask.txt'
            
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
                         
            call cloud_mask_naive_bayes ( inp, erg , info_flags )      
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
