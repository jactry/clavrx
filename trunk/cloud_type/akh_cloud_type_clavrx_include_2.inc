!$Id: akh_cloud_type_clavrx_include_2.inc,v 1.12 2014/01/10 01:57:03 heidinger Exp $
!---------------------------------------------------------------------
!
!----------------------------------------------------------------------
  Num_Elem = Num_Pix      
  Elem_End = Num_Elem
  Elem_Start = 1
  Channel_On_Flag => Chan_On_Flag_Default
  Invalid_Pixel_Mask => Bad_Pixel_Mask
  Vza_Idx_Rtm => Zen_Idx_Rtm
  if ((Channel_On_Flag(Chan_Idx_16um)) /= 0) Ref_16um => ch(6)%Ref_Toa
  if ((Channel_On_Flag(Chan_Idx_375um)) /= 0) Ref_375um => ch(20)%Ref_Toa
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Rad_11um => ch(31)%Rad_Toa
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Rad_11um_Clear => ch(31)%Rad_Toa_Clear
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Bt_11um => ch(31)%Bt_Toa
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Bt_11um_Mean => Bt_Ch31_Mean_3x3
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Bt_11um_Stddev => Bt_Ch31_Std_3x3
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Bt_11um_Min => Bt_Ch31_Min_3x3
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Bt_11um_Max => Bt_Ch31_Max_3x3
  if ((Channel_On_Flag(Chan_Idx_11um)) /= 0) Bt_11um_Clear => ch(31)%Bt_Toa_Clear  
  if ((Channel_On_Flag(Chan_Idx_12um)) /= 0) Bt_12um => ch(32)%Bt_Toa
  if ((Channel_On_Flag(Chan_Idx_12um)) /= 0) Bt_12um_Clear => ch(32)%Bt_Toa_Clear 
  if ((Channel_On_Flag(Chan_Idx_67um)) /= 0) Rad_H2O => ch(27)%Rad_Toa
  if ((Channel_On_Flag(Chan_Idx_67um)) /= 0) Bt_H2O => ch(27)%Bt_Toa
  if ((Channel_On_Flag(Chan_Idx_67um)) /= 0) Bt_H2O_Stddev => Bt_Ch27_Std_3x3
  if ((Channel_On_Flag(Chan_Idx_67um)) /= 0) Rad_H2O_Clear => ch(27)%Rad_Toa_Clear
  if ((Channel_On_Flag(Chan_Idx_67um)) /= 0) Bt_H2O_Clear => ch(27)%Bt_Toa_Clear
  if (((Channel_On_Flag(Chan_Idx_67um)) /= 0) .and. (Channel_On_Flag(Chan_Idx_11um) /= 0)) Covar_H2O_Window => Covar_Ch27_Ch31_5x5 
  Tcld_H2O => Tc_H2O
  Pcld_H2O => Pc_H2O
  Zcld_H2O => Zc_H2O
  Tcld_Opa => Tc_Opaque_Cloud
  Pcld_Opa => Pc_Opaque_Cloud
  Zcld_Opa => Zc_Opaque_Cloud
  Zclr_H2O => Zclr_H2O_Peak
  Cloud_Phase => Cld_Phase
  Cloud_Type => Cld_Type
  Cloud_Mask => Cld_Mask
  Emiss_Window_Tropo => ch(31)%Emiss_Tropo
  Emiss_H2O_Tropo => ch(27)%Emiss_Tropo
  if ((Channel_On_Flag(Chan_Idx_85um)  /= 0) .and. (Channel_On_Flag(Chan_Idx_11um) /= 0)) Beta_11um_85um_Tropo => Beta_11um_85um_Tropo_Rtm
  if ((Channel_On_Flag(Chan_Idx_12um)  /= 0) .and. (Channel_On_Flag(Chan_Idx_11um) /= 0)) Beta_11um_12um_Tropo => Beta_11um_12um_Tropo_Rtm
  if ((Channel_On_Flag(Chan_Idx_133um) /= 0) .and. (Channel_On_Flag(Chan_Idx_11um) /= 0)) Beta_11um_133um_Tropo => Beta_11um_133um_Tropo_Rtm
  Number_Rtm_Levels = Nlevels_Rtm
  allocate(P_Prof(Number_Rtm_Levels))
  P_Prof = P_Std_Rtm
