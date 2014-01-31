!$Id: akh_cloud_type_gsip_include_2.inc,v 1.1 2013/09/12 18:23:59 wstraka Exp $
!---------------------------------------------------------------------
!
!----------------------------------------------------------------------
  Num_Elem = Num_Pix      
  Elem_End = Num_Elem
  Elem_Start = 1
  Channel_On_Flag => sat_info_gsip(1)%chanon
  Invalid_Pixel_Mask => bad_pix_mask(14,:,:)
  Vza_Idx_Rtm => ivza_rtm
  if ((Channel_On_Flag(Chan_Idx_375um)) /= 0) Ref_375um => ref7
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Rad_11um => rad14
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Rad_11um_Clear => rad_clear_ch14_rtm
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Bt_11um => bt14
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Bt_11um_Mean => Bt_Ch31_Mean_3x3
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Bt_11um_Min => Bt_Ch31_Min_3x3
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Bt_11um_Max => Bt_Ch31_Max_3x3
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Bt_11um_Clear => bt_clear_ch14_rtm 
  if ((Channel_On_Flag(Chan_Idx_12um)) /= 0) Bt_12um => bt15
  if ((Channel_On_Flag(Chan_Idx_12um)) /= 0) Bt_12um_Clear => bt_clear_ch15_rtm 
  if ((Channel_On_Flag(Chan_Idx_67um)) /= 0) Rad_H2O => rad9
  if ((Channel_On_Flag(Chan_Idx_67um)) /= 0) Rad_H2O_Clear => Rad_Atm_ch9_Rtm
  if ((Channel_On_Flag(Chan_Idx_67um)) /= 0) Bt_H2O_Clear => bt_clear_ch9_rtm
  if (((Channel_On_Flag(Chan_Idx_67um)) /= 0) .and. (Channel_On_Flag(Chan_Idx_11um) /= 0)) Covar_H2O_Window => Covar_Ch27_Ch31_5x5 
  Cirrus_Flag => gsip_pix_prod%i1_generic1
  Water_Flag => gsip_pix_prod%i1_generic2
  Ice_Probability =>  gsip_pix_prod%r4_generic8
  Tcld_H2O => gsip_pix_prod%r4_generic1
  Pcld_H2O => gsip_pix_prod%r4_generic2
  Zcld_H2O => gsip_pix_prod%r4_generic3
  Tcld_Opa => gsip_pix_prod%r4_generic4
  Pcld_Opa => gsip_pix_prod%r4_generic5
  Zcld_Opa => gsip_pix_prod%r4_generic6
  Zclr_H2O => gsip_pix_prod%r4_generic7
  Cloud_Phase => gsip_pix_prod%Cldphase
  Cloud_Type => gsip_pix_prod%Cldtype
  Cloud_Mask => gsip_pix_prod%cldmask
  Emiss_Window_Tropo => Emiss_11um_Tropo_Rtm
  Emiss_H2O_Tropo => Emiss_67um_Tropo_Rtm
  if ((Channel_On_Flag(Chan_Idx_85um)  /= 0) .and. (Channel_On_Flag(Chan_Idx_11um) /= 0)) Beta_11um_85um_Tropo => Beta_11um_85um_Tropo_Rtm
  if ((Channel_On_Flag(Chan_Idx_12um)  /= 0) .and. (Channel_On_Flag(Chan_Idx_11um) /= 0)) Beta_11um_12um_Tropo => Beta_11um_12um_Tropo_Rtm
  if ((Channel_On_Flag(Chan_Idx_133um) /= 0) .and. (Channel_On_Flag(Chan_Idx_11um) /= 0)) Beta_11um_133um_Tropo => Beta_11um_133um_Tropo_Rtm
  Number_Rtm_Levels = Nlevels_Rtm
  allocate(P_Prof(Number_Rtm_Levels))
  P_Prof = P_Std_Rtm
