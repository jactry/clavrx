!$Id:$
module ACHA_MICROPHYSICAL_MODULE

implicit none

public:: SETUP_ICE_MICROPHYSICAL_MODEL


  !  acha ice model parameters
  real, public, save:: A_BETA_11um_67um_FIT_ICE
  real, public, save:: B_BETA_11um_67um_FIT_ICE
  real, public, save:: A_BETA_11um_85um_FIT_ICE
  real, public, save:: B_BETA_11um_85um_FIT_ICE
  real, public, save:: A_BETA_11um_133um_FIT_ICE
  real, public, save:: B_BETA_11um_133um_FIT_ICE
  real, public, save:: A_Re_Beta_FIT_ICE
  real, public, save:: B_Re_Beta_FIT_ICE
  real, public, save:: C_Re_Beta_FIT_ICE
  real, public, save:: D_Re_Beta_FIT_ICE
  real, public, save:: A_Qe_065um_FIT_ICE
  real, public, save:: B_Qe_065um_FIT_ICE
  real, public, save:: C_Qe_065um_FIT_ICE
  real, public, save:: A_Qe_11um_FIT_ICE
  real, public, save:: B_Qe_11um_FIT_ICE
  real, public, save:: C_Qe_11um_FIT_ICE
  real, public, save:: A_wo_11um_FIT_ICE
  real, public, save:: B_wo_11um_FIT_ICE
  real, public, save:: C_wo_11um_FIT_ICE
  real, public, save:: A_g_11um_FIT_ICE
  real, public, save:: B_g_11um_FIT_ICE
  real, public, save:: C_g_11um_FIT_ICE

  !--- water microphysical model terms
  real, public, save:: A_BETA_11um_133um_FIT_WATER 
  real, public, save:: B_BETA_11um_133um_FIT_WATER
  real, public, save:: A_BETA_11um_85um_FIT_WATER 
  real, public, save:: B_BETA_11um_85um_FIT_WATER
  real, public, save:: A_BETA_11um_67um_FIT_WATER
  real, public, save:: B_BETA_11um_67um_FIT_WATER
  real, public, save:: A_Re_Beta_FIT_WATER
  real, public, save:: B_Re_Beta_FIT_WATER
  real, public, save:: C_Re_Beta_FIT_WATER
  real, public, save:: D_Re_Beta_FIT_WATER
  real, public, save:: A_Qe_065um_FIT_WATER
  real, public, save:: B_Qe_065um_FIT_WATER
  real, public, save:: C_Qe_065um_FIT_WATER
  real, public, save:: A_Qe_11um_FIT_WATER
  real, public, save:: B_Qe_11um_FIT_WATER
  real, public, save:: C_Qe_11um_FIT_WATER
  real, public, save:: A_wo_11um_FIT_WATER
  real, public, save:: B_wo_11um_FIT_WATER
  real, public, save:: C_wo_11um_FIT_WATER
  real, public, save:: A_g_11um_FIT_WATER
  real, public, save:: B_g_11um_FIT_WATER
  real, public, save:: C_g_11um_FIT_WATER

  contains

  subroutine SETUP_ICE_MICROPHYSICAL_MODEL(WMO_Id)

   integer, intent(in):: WMO_Id

   !--- water microphysical model derived from Mie theory by A. Heidinger
   A_BETA_11um_133um_FIT_WATER = -0.728113
   B_BETA_11um_133um_FIT_WATER = 1.743389
   A_BETA_11um_85um_FIT_WATER = 0.930569
   B_BETA_11um_85um_FIT_WATER = 0.048857
   A_BETA_11um_67um_FIT_WATER = 0.268115
   B_BETA_11um_67um_FIT_WATER = 0.702683
   A_Re_Beta_FIT_WATER = -48.36295
   B_Re_Beta_FIT_WATER = 120.40642
   C_Re_Beta_FIT_WATER = -99.39934
   D_Re_Beta_FIT_WATER =  27.70454
   A_Qe_065um_FIT_WATER =   2.48801
   B_Qe_065um_FIT_WATER =  -0.55097
   C_Qe_065um_FIT_WATER =   0.15898
   A_Qe_11um_FIT_WATER =  -0.28601
   B_Qe_11um_FIT_WATER =   1.82868
   C_Qe_11um_FIT_WATER =  -0.35691
   A_wo_11um_FIT_WATER =  -0.12614
   B_wo_11um_FIT_WATER =   0.51702
   C_wo_11um_FIT_WATER =  -0.09750
   A_g_11um_FIT_WATER =   0.11256
   B_g_11um_FIT_WATER =   0.96132
   C_g_11um_FIT_WATER =  -0.27036

   !-- ice clouds modeled as aggregate_columns (b=0.10)
   select case(WMO_Id)

     case(3:5,200:223,706:708)  !avhrr, goes-im
        A_BETA_11um_67um_FIT_ICE =    1.05090
        B_BETA_11um_67um_FIT_ICE =   -0.01637
        A_BETA_11um_85um_FIT_ICE =    1.50435
        B_BETA_11um_85um_FIT_ICE =   -0.48727
        A_BETA_11um_133um_FIT_ICE =    0.03881
        B_BETA_11um_133um_FIT_ICE =    1.01095
        A_Re_Beta_FIT_ICE =  -10.93929
        B_Re_Beta_FIT_ICE =   24.70588
        C_Re_Beta_FIT_ICE =  -17.89442
        D_Re_Beta_FIT_ICE =    4.61250
        A_Qe_065um_FIT_ICE =    2.24008
        B_Qe_065um_FIT_ICE =   -0.22154
        C_Qe_065um_FIT_ICE =    0.05269
        A_Qe_11um_FIT_ICE =   -0.25110
        B_Qe_11um_FIT_ICE =    2.30356
        C_Qe_11um_FIT_ICE =   -0.59547
        A_wo_11um_FIT_ICE =   -0.00289
        B_wo_11um_FIT_ICE =    0.58422
        C_wo_11um_FIT_ICE =   -0.16068
        A_g_11um_FIT_ICE =    0.63572
        B_g_11um_FIT_ICE =    0.39030
        C_g_11um_FIT_ICE =   -0.10809

     case(173)  !goes-np 3 chan
        A_BETA_11um_67um_FIT_ICE =    2.17382
        B_BETA_11um_67um_FIT_ICE =   -1.10936
        A_BETA_11um_85um_FIT_ICE =    2.94619
        B_BETA_11um_85um_FIT_ICE =   -1.88164
        A_BETA_11um_133um_FIT_ICE =    0.10731
        B_BETA_11um_133um_FIT_ICE =    0.93036
        A_Re_Beta_FIT_ICE = -153.69155
        B_Re_Beta_FIT_ICE =  403.09795
        C_Re_Beta_FIT_ICE = -353.06676
        D_Re_Beta_FIT_ICE =  103.95095
        A_Qe_065um_FIT_ICE =    2.24167
        B_Qe_065um_FIT_ICE =   -0.22339
        C_Qe_065um_FIT_ICE =    0.05323
        A_Qe_11um_FIT_ICE =    0.18912
        B_Qe_11um_FIT_ICE =    1.95161
        C_Qe_11um_FIT_ICE =   -0.52995
        A_wo_11um_FIT_ICE =    0.03123
        B_wo_11um_FIT_ICE =    0.52542
        C_wo_11um_FIT_ICE =   -0.13696
        A_g_11um_FIT_ICE =    0.62559
        B_g_11um_FIT_ICE =    0.40393
        C_g_11um_FIT_ICE =   -0.11503

     case(224) ! VIIRS 
        A_BETA_11um_67um_FIT_ICE =    0.95539
        B_BETA_11um_67um_FIT_ICE =    0.07902
        A_BETA_11um_85um_FIT_ICE =    1.40457
        B_BETA_11um_85um_FIT_ICE =   -0.39163
        A_BETA_11um_133um_FIT_ICE =   -0.02641
        B_BETA_11um_133um_FIT_ICE =    1.08386
        A_Re_Beta_FIT_ICE =   -8.46300
        B_Re_Beta_FIT_ICE =   18.82502
        C_Re_Beta_FIT_ICE =  -13.19113
        D_Re_Beta_FIT_ICE =    3.34004
        A_Qe_065um_FIT_ICE =    2.25389
        B_Qe_065um_FIT_ICE =   -0.23436
        C_Qe_065um_FIT_ICE =    0.05582
        A_Qe_11um_FIT_ICE =   -0.34298
        B_Qe_11um_FIT_ICE =    2.42808
        C_Qe_11um_FIT_ICE =   -0.63351
        A_wo_11um_FIT_ICE =    0.00705
        B_wo_11um_FIT_ICE =    0.59694
        C_wo_11um_FIT_ICE =   -0.16978
        A_g_11um_FIT_ICE =    0.64054
        B_g_11um_FIT_ICE =    0.37904
        C_g_11um_FIT_ICE =   -0.10336
     case default  !MODIS
        A_BETA_11um_67um_FIT_ICE =    1.54172
        B_BETA_11um_67um_FIT_ICE =   -0.49635
        A_BETA_11um_85um_FIT_ICE =    2.17033
        B_BETA_11um_85um_FIT_ICE =   -1.14189
        A_BETA_11um_133um_FIT_ICE =    0.11629
        B_BETA_11um_133um_FIT_ICE =    0.93978
        A_Re_Beta_FIT_ICE =  -39.22696
        B_Re_Beta_FIT_ICE =   97.67353
        C_Re_Beta_FIT_ICE =  -81.09035
        D_Re_Beta_FIT_ICE =   23.08110
        A_Qe_065um_FIT_ICE =    2.24408
        B_Qe_065um_FIT_ICE =   -0.23005
        C_Qe_065um_FIT_ICE =    0.05641
        A_Qe_11um_FIT_ICE =   -0.02317
        B_Qe_11um_FIT_ICE =    2.12770
        C_Qe_11um_FIT_ICE =   -0.57337
        A_wo_11um_FIT_ICE =   -0.01382
        B_wo_11um_FIT_ICE =    0.57493
        C_wo_11um_FIT_ICE =   -0.15368
        A_g_11um_FIT_ICE =    0.59489
        B_g_11um_FIT_ICE =    0.46862
        C_g_11um_FIT_ICE =   -0.14131
   end select

   end subroutine SETUP_ICE_MICROPHYSICAL_MODEL

end module ACHA_MICROPHYSICAL_MODULE
