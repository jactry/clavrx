!$Id: akh_cloud_type_clavrx_include_3.inc,v 1.3 2014/01/16 13:40:50 heidinger Exp $
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
     Rad_H2O_BB_Profile => Rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(27)%Rad_BB_Cloud_Profile
     Rad_11um_BB_Profile => Rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(31)%Rad_BB_Cloud_Profile
     T_Prof => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%T_prof
     Z_Prof => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Z_prof
     Tropo_Level => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level
     Sfc_Level => rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Sfc_Level
