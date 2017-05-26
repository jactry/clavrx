  !$Id: awg_cld_hght_clavrx_include_2.inc,v 1.2 2011/08/16 18:45:07 heidinger Exp $
  !----------------------------------------------------------------------------
  ! This code contains aliasing statements between ACHA local variables
  ! and global variables in the CLAVR-x/PATMOS-x processing system
  !
  ! An analogous include file exists for the Geocat processing system
  !
  ! Begin Alias Block #2
  !----------------------------------------------------------------------------
  Number_Of_Elements = Imager@Number_Of_Elements
  Number_Of_Lines = Imager%Number_Of_Lines_Per_Segment
  Smooth_Nwp_Fields_Flag = Smooth_Nwp_Flag
  Acha_Mode_Flag = Acha%Mode
  Processing_Order => Temp_Mask
  Process_Undetected_Cloud_Flag_Local = Process_Undetected_Cloud_Flag
  Tc => Acha%Tc
  Ec => Acha%Ec
  Beta => ACHA%Beta
  Pc => ACHA%Pc
  Zc => ACHA%Zc
  Tau => ACHA%Tau
  Tc_Uncertainty => ACHA%Tc_Uncertainty
  Ec_Uncertainty => ACHA%Ec_Uncertainty
  Beta_Uncertainty => ACHA%Beta_Uncertainty
  Zc_Uncertainty => ACHA%Zc_Uncertainty
  Pc_Uncertainty => ACHA%Pc_Uncertainty
  Qf => Acha%Quality_Flag
  Elem_Idx_LRC_Clavrx => I_LRC
  Line_Idx_LRC_Clavrx => J_LRC
  Cloud_Layer => ACHA%Cld_Layer
  Cloud_Mask_Local => CLDMASK%Cld_Mask
  Cloud_Type_Local => Cld_Type
  Ch20_Surface_Emissivity => ch(20)%Sfc_Emiss
  Surface_Elevation => Sfc%Zsfc
  Latitude => Nav%Lat
  Snow_Class => Sfc%Snow
  Surface_Type => Sfc%Sfc_Type
  Elem_Idx_NWP => I_Nwp
  Line_Idx_NWP => J_Nwp
  Viewing_Zenith_Angle_Idx_RTM => Zen_Idx_Rtm
  Invalid_Data_Mask => Bad_Pixel_Mask
  Chan_On => Sensor%Chan_On_Flag_Default
  if (Chan_On(Chan_Idx_67um) == sym%YES) Bt_67um => Bt_Ch27
  if (Chan_On(Chan_Idx_85um) == sym%YES) Bt_85um => Bt_Ch29
  if (Chan_On(Chan_Idx_11um) == sym%YES) Bt_11um => Bt_Ch31
  if (Chan_On(Chan_Idx_11um) == sym%YES) Rad_11um => Rad_Ch31
  if (Chan_On(Chan_Idx_12um) == sym%YES) Bt_12um => Bt_Ch32
  if (Chan_On(Chan_Idx_133um) == sym%YES) Bt_133um => Bt_Ch33
  Cosine_Zenith_Angle => Geo%Coszen
  Lower_Cloud_Pressure => ACHA%Pc_Lower_Cloud
  Lower_Cloud_Temperature => ACHA%Tc_Lower_Cloud
  Lower_Cloud_Height => ACHA%Zc_Lower_Cloud
  Surface_Pressure => Psfc_Nwp
  Rad_Clear_67um_Local => ch(27)%Rad_Toa_Clear
  Rad_Clear_85um_Local => ch(29)%Rad_Toa_Clear
  Rad_Clear_11um_Local => ch(31)%Rad_Toa_Clear
  Rad_Clear_12um_Local => ch(32)%Rad_Toa_Clear
  Rad_Clear_133um_Local => ch(33)%Rad_Toa_Clear
  OE_Qf => Acha%OE_Quality_Flags
  Inversion_Present_Flag => Rtm(:,:)%Inversion_Level 
  !----------------------------------------------------------------------------
  ! End Alias Block #2
  !----------------------------------------------------------------------------
