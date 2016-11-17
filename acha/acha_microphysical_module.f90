!$Id:$
module ACHA_MICROPHYSICAL_MODULE

implicit none

integer, parameter, private:: Habit_Idx = 7

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
        include 'acha_ice_cloud_microphysical_model_noaa19.inc'

     case(173)  !goes-np 3 chan
        include 'acha_ice_cloud_microphysical_model_abi.inc'

     case(224) ! VIIRS 
        include 'acha_ice_cloud_microphysical_model_viirs.inc'

     case default  !MODIS
        include 'acha_ice_cloud_microphysical_model_modis.inc'

   end select

   end subroutine SETUP_ICE_MICROPHYSICAL_MODEL

end module ACHA_MICROPHYSICAL_MODULE
