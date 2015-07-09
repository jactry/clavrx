!$Id$
 !--- level2 indices
 integer, parameter, private:: Num_Level2_Sds = 237
 integer, private, save:: Sd_Id_Level2
 integer(kind=int4), dimension(Num_Level2_Sds), save, private:: Sds_Id_Level2
 integer(kind=int4), public, parameter:: NCDC_Attribute_Flag = 0

!--- set sds order and flags for level2 files
!
! Note, fields marked "critical for level-2b" are needed for level-2b generation
!
 integer, private, parameter::    &
     Sds_Num_Level2_Scanline = 1,        Sds_Num_Level2_scanline_Flag = 1, &   !critical for level2b
     Sds_Num_Level2_Time = 2,            Sds_Num_Level2_Time_Flag = 1, &       !critical for level2b
     Sds_Num_Level2_Bad_Scan = 3,        Sds_Num_Level2_Bad_Scan_Flag = 0, &
     Sds_Num_Level2_Bad_Pixel_Mask = 4,  Sds_Num_Level2_Bad_Pixel_Mask_Flag = 1, &  !critical for level2b
     Sds_Num_Level2_Gap_Pixel_Mask = 5,  Sds_Num_Level2_Gap_Pixel_Mask_Flag = 0, &  !critical for level2b for VIIRS
     Sds_Num_Level2_Asc_Flag = 6,        Sds_Num_Level2_Asc_Flag_Flag = 1, &        !critical for level2b (except GOES)
     Sds_Num_Level2_Meta_Data = 7,       Sds_Num_Level2_Meta_Data_Flag = 0, &
     Sds_Num_Level2_Lat = 8,             Sds_Num_Level2_Lat_Flag = 1, &         !critical for level2b
     Sds_Num_Level2_Lon = 9,             Sds_Num_Level2_Lon_Flag = 1, &         !critical for level2b
     Sds_Num_Level2_Zen = 10,            Sds_Num_Level2_Zen_Flag = 1, &         !critical for level2b
     Sds_Num_Level2_Solzen = 11,         Sds_Num_Level2_Solzen_Flag = 1, &      !critical for level2b
     Sds_Num_Level2_Relaz = 12 ,         Sds_Num_Level2_Relaz_Flag = 1, &       
     Sds_Num_Level2_Solaz = 13 ,         Sds_Num_Level2_Solaz_Flag = 1, &       
     Sds_Num_Level2_Packed_Land = 14,    Sds_Num_Level2_Packed_Land_Flag = 0, &
     Sds_Num_Level2_Ch1 = 15 ,           Sds_Num_Level2_Ch1_Flag = 1, &         !critical for level2b
     Sds_Num_Level2_Ch2 = 16 ,           Sds_Num_Level2_Ch2_Flag = 1, &
     Sds_Num_Level2_Ch3 = 17 ,           Sds_Num_Level2_Ch3_Flag = 0, &
     Sds_Num_Level2_Ch4 = 18 ,           Sds_Num_Level2_Ch4_Flag = 0, &
     Sds_Num_Level2_Ch5 = 19 ,           Sds_Num_Level2_Ch5_Flag = 0, &
     Sds_Num_Level2_Ch6 = 20 ,           Sds_Num_Level2_Ch6_Flag = 1, &
     Sds_Num_Level2_Ch7 = 21 ,           Sds_Num_Level2_Ch7_Flag = 0, &
     Sds_Num_Level2_Ch8 = 22 ,           Sds_Num_Level2_Ch8_Flag = 0, &
     Sds_Num_Level2_Ch9 = 23 ,           Sds_Num_Level2_Ch9_Flag = 0, &
     Sds_Num_Level2_Ch10 = 24 ,          Sds_Num_Level2_Ch10_Flag = 0, &
     Sds_Num_Level2_Ch11 = 25 ,          Sds_Num_Level2_Ch11_Flag = 0, &
     Sds_Num_Level2_Ch12 = 26 ,          Sds_Num_Level2_Ch12_Flag = 0, &
     Sds_Num_Level2_Ch13 = 27 ,          Sds_Num_Level2_Ch13_Flag = 0, &
     Sds_Num_Level2_Ch14 = 28 ,          Sds_Num_Level2_Ch14_Flag = 0, &
     Sds_Num_Level2_Ch15 = 29 ,          Sds_Num_Level2_Ch15_Flag = 0, &
     Sds_Num_Level2_Ch16 = 30 ,          Sds_Num_Level2_Ch16_Flag = 0, &
     Sds_Num_Level2_Ch17 = 31 ,          Sds_Num_Level2_Ch17_Flag = 0, &
     Sds_Num_Level2_Ch18 = 32 ,          Sds_Num_Level2_Ch18_Flag = 0, &
     Sds_Num_Level2_Ch19 = 33 ,          Sds_Num_Level2_Ch19_Flag = 0, &
     Sds_Num_Level2_Ch20_Ref = 34 ,      Sds_Num_Level2_Ch20_Ref_Flag = 1, &
     Sds_Num_Level2_Ch20_Bt = 35 ,       Sds_Num_Level2_Ch20_Bt_Flag = 1, &
     Sds_Num_Level2_Ch21 = 36 ,          Sds_Num_Level2_Ch21_Flag = 0, &
     Sds_Num_Level2_Ch22 = 37 ,          Sds_Num_Level2_Ch22_Flag = 0, &
     Sds_Num_Level2_Ch23 = 38 ,          Sds_Num_Level2_Ch23_Flag = 0, &
     Sds_Num_Level2_Ch24 = 39 ,          Sds_Num_Level2_Ch24_Flag = 0, &
     Sds_Num_Level2_Ch25 = 40 ,          Sds_Num_Level2_Ch25_Flag = 0, &
     Sds_Num_Level2_Ch26 = 41 ,          Sds_Num_Level2_Ch26_Flag = 0, &
     Sds_Num_Level2_Ch27 = 42 ,          Sds_Num_Level2_Ch27_Flag = 0, &
     Sds_Num_Level2_Ch28 = 43 ,          Sds_Num_Level2_Ch28_Flag = 0, &
     Sds_Num_Level2_Ch29 = 44 ,          Sds_Num_Level2_Ch29_Flag = 0, &
     Sds_Num_Level2_Ch30 = 45 ,          Sds_Num_Level2_Ch30_Flag = 0, &
     Sds_Num_Level2_Ch31 = 46 ,          Sds_Num_Level2_Ch31_Flag = 1, &       !critical for level2b
     Sds_Num_Level2_Ch32 = 47 ,          Sds_Num_Level2_Ch32_Flag = 1, &
     Sds_Num_Level2_Ch33 = 48 ,          Sds_Num_Level2_Ch33_Flag = 0, &
     Sds_Num_Level2_Ch34 = 49 ,          Sds_Num_Level2_Ch34_Flag = 0, &
     Sds_Num_Level2_Ch35 = 50 ,          Sds_Num_Level2_Ch35_Flag = 0, &
     Sds_Num_Level2_Ch36 = 51 ,          Sds_Num_Level2_Ch36_Flag = 0, &
     Sds_Num_Level2_ChDNB = 52 ,         Sds_Num_Level2_ChDNB_Flag = 0, &
     Sds_Num_Level2_Tsfc_Back = 53,      Sds_Num_Level2_Tsfc_Back_Flag = 0, &
     Sds_Num_Level2_Ch1_Std = 54 ,       Sds_Num_Level2_Ch1_Std_Flag = 1, &
     Sds_Num_Level2_Ch31_Std = 55 ,      Sds_Num_Level2_Ch31_Std_Flag = 1, &
     Sds_Num_Level2_Cldprob = 56,        Sds_Num_Level2_Cldprob_Flag = 1, &
     Sds_Num_Level2_Cld_Type = 57 ,      Sds_Num_Level2_Cld_Type_Flag = 1, &
     Sds_Num_Level2_Cld_Phase = 58 ,     Sds_Num_Level2_Cld_Phase_Flag = 1, &
     Sds_Num_Level2_Ctp = 59,            Sds_Num_Level2_Ctp_Flag = 1, &
     Sds_Num_Level2_Ctt = 60,            Sds_Num_Level2_Ctt_Flag = 1, &
     Sds_Num_Level2_Cth = 61,            Sds_Num_Level2_Cth_Flag = 1, &
     Sds_Num_Level2_Ec = 62,             Sds_Num_Level2_Ec_Flag = 1, &
     Sds_Num_Level2_Beta = 63,           Sds_Num_Level2_Beta_Flag = 1, &
     Sds_Num_Level2_Cod_Acha = 64,       Sds_Num_Level2_Cod_Acha_Flag = 1, &
     Sds_Num_Level2_Ceps_Acha = 65,      Sds_Num_Level2_Ceps_Acha_Flag = 1, &
     Sds_Num_Level2_Cod = 66,            Sds_Num_Level2_Cod_Flag = 1, &
     Sds_Num_Level2_Ceps = 67,           Sds_Num_Level2_Ceps_Flag = 1, &
     Sds_Num_Level2_Cldalb = 68,         Sds_Num_Level2_Cldalb_Flag = 1, &
     Sds_Num_Level2_Cldtrn = 69,         Sds_Num_Level2_Cldtrn_Flag = 1, &
     Sds_Num_Level2_Cldfrac = 70,        Sds_Num_Level2_Cldfrac_Flag = 1, &
     Sds_Num_Level2_Cldfrac_Uncer = 71,  Sds_Num_Level2_Cldfrac_Uncer_Flag = 1, &
     Sds_Num_Level2_Etrop = 72,          Sds_Num_Level2_Etrop_Flag = 0, &
     Sds_Num_Level2_Aot1 = 73,           Sds_Num_Level2_Aot1_Flag = 0, & !off for GOES
     Sds_Num_Level2_Olr = 74,            Sds_Num_Level2_Olr_Flag = 1, & !off for GOES
     Sds_Num_Level2_Ndvi = 75,           Sds_Num_Level2_Ndvi_Flag = 0, & !off for GOES
     Sds_Num_Level2_Tsfc = 76,           Sds_Num_Level2_Tsfc_Flag = 1, &
     Sds_Num_Level2_Tair = 77,           Sds_Num_Level2_Tair_Flag = 0, &
     Sds_Num_Level2_Trad = 78,           Sds_Num_Level2_Trad_Flag = 0, &
     Sds_Num_Level2_Rsr = 79,            Sds_Num_Level2_Rsr_Flag = 0, &
     Sds_Num_Level2_Qf1 = 80,            Sds_Num_Level2_Qf1_Flag = 0, &		!deprecated
     Sds_Num_Level2_Qf2 = 81,            Sds_Num_Level2_Qf2_Flag = 0, &		!deprecated
     Sds_Num_Level2_Cod_Dcomp_Uncer = 82,Sds_Num_Level2_Cod_Dcomp_Uncer_Flag = 1, &
     Sds_Num_Level2_Ceps_Dcomp_Uncer = 83,Sds_Num_Level2_Ceps_Dcomp_Uncer_Flag = 1, &
     Sds_Num_Level2_Cth_Acha_Uncer = 84, Sds_Num_Level2_Cth_Acha_Uncer_Flag = 1, &
     Sds_Num_Level2_Insol = 85,          Sds_Num_Level2_Insol_Flag = 1, &
     Sds_Num_Level2_Insol_Dif = 86,      Sds_Num_Level2_Insol_Dif_Flag = 1, &
     Sds_Num_Level2_Ch1_Counts = 87,     Sds_Num_Level2_Ch1_Counts_Flag = 1, &
     Sds_Num_Level2_Ch2_Counts = 88,     Sds_Num_Level2_Ch2_Counts_Flag = 1, &
     Sds_Num_Level2_Ch6_Counts = 89,     Sds_Num_Level2_Ch6_Counts_Flag = 1, &
     Sds_Num_Level2_Tpw = 90,            Sds_Num_Level2_Tpw_Flag = 0, &
     Sds_Num_Level2_Rh = 91,             Sds_Num_Level2_Rh_Flag = 0, &
     Sds_Num_Level2_Ref_Ch20_Sfc = 92,   Sds_Num_Level2_Ref_Ch20_Sfc_Flag = 0, &
     Sds_Num_Level2_Ref_Ch1_Sfc = 93,    Sds_Num_Level2_Ref_Ch1_Sfc_Flag = 0, &
     Sds_Num_Level2_Ref_Ch2_Sfc = 94,    Sds_Num_Level2_Ref_Ch2_Sfc_Flag = 0, &
     Sds_Num_Level2_Ref_Ch6_Sfc = 95,    Sds_Num_Level2_Ref_Ch6_Sfc_Flag = 0, &
     Sds_Num_Level2_Sst_Masked = 96,     Sds_Num_Level2_Sst_Masked_Flag = 0, &
     Sds_Num_Level2_Ch1_Unnorm = 97,     Sds_Num_Level2_Ch1_Unnorm_Flag = 0, &
     Sds_Num_Level2_Ch2_Unnorm = 98,     Sds_Num_Level2_Ch2_Unnorm_Flag = 0, &
     Sds_Num_Level2_Ch6_Unnorm = 99,     Sds_Num_Level2_Ch6_Unnorm_Flag = 0, &
     Sds_Num_Level2_Ch1_Clear = 100,     Sds_Num_Level2_Ch1_Clear_Flag = 0, &
     Sds_Num_Level2_Ch1_Mean = 101,      Sds_Num_Level2_Ch1_Mean_Flag = 0, &
     Sds_Num_Level2_Sst_Unmasked = 102,  Sds_Num_Level2_Sst_Unmasked_Flag = 0, &
     Sds_Num_Level2_Wnd_Spd = 103,       Sds_Num_Level2_Wnd_Spd_Flag = 0, &
     Sds_Num_Level2_Wnd_Dir = 104,       Sds_Num_Level2_Wnd_Dir_Flag = 0, &
     Sds_Num_Level2_Ch1_Dark = 105,      Sds_Num_Level2_Ch1_Dark_Flag = 0, &
     Sds_Num_Level2_Cwp = 106,           Sds_Num_Level2_Cwp_Flag = 1, &
     Sds_Num_Level2_Wnd_Spd_Cld_Top= 107, Sds_Num_Level2_Wnd_Spd_Cld_Top_Flag = 0, &
     Sds_Num_Level2_Wnd_Dir_Cld_Top= 108, Sds_Num_Level2_Wnd_Dir_Cld_Top_Flag = 0, &
     Sds_Num_Level2_Temp_Tropo = 109,    Sds_Num_Level2_Temp_Tropo_Flag = 0, &
     Sds_Num_Level2_Cld_Mask =  110,     Sds_Num_Level2_Cld_Mask_Flag = 1, &
     Sds_Num_Level2_Sfc_Type = 111,      Sds_Num_Level2_Sfc_Type_Flag = 1, &
     Sds_Num_Level2_Coast_Mask = 112,    Sds_Num_Level2_Coast_Mask_Flag = 0, &
     Sds_Num_Level2_Land_Mask = 113,     Sds_Num_Level2_Land_Mask_Flag = 1, &
     Sds_Num_Level2_Snow_Mask = 114,     Sds_Num_Level2_Snow_Mask_Flag = 1, &
     Sds_Num_Level2_Cod_Dcomp_Qf = 115,  Sds_Num_Level2_Cod_Dcomp_Qf_Flag = 0, &
     Sds_Num_Level2_Ceps_Dcomp_Qf = 116, Sds_Num_Level2_Ceps_Dcomp_Qf_Flag = 0, &
     Sds_Num_Level2_Cth_Acha_Qf = 117,   Sds_Num_Level2_Cth_Acha_Qf_Flag = 0, &
     Sds_Num_Level2_Ec_Acha_Qf = 118,    Sds_Num_Level2_Ec_Acha_Qf_Flag = 0, &
     Sds_Num_Level2_Beta_Acha_Qf = 119,  Sds_Num_Level2_Beta_Acha_Qf_Flag = 0, &
     Sds_Num_Level2_Aot2 = 120,          Sds_Num_Level2_Aot2_Flag = 1, &
     Sds_Num_Level2_Aot6 = 121,          Sds_Num_Level2_Aot6_Flag = 1, &
     Sds_Num_Level2_Aot_QF = 122,        Sds_Num_Level2_Aot_QF_Flag = 1, &
     Sds_Num_Level2_Zsfc = 123,          Sds_Num_Level2_Zsfc_Flag = 0, &
     Sds_Num_Level2_Cth_H2O = 124,       Sds_Num_Level2_Cth_H2O_Flag = 0, &
     Sds_Num_Level2_Cth_Opa = 125,       Sds_Num_Level2_Cth_Opa_Flag = 0, &
     Sds_Num_Level2_Cld_Mask_Aux = 126,  Sds_Num_Level2_Cld_Mask_Aux_Flag = 0, & !off for GOES
     Sds_Num_Level2_Cld_Type_Aux = 127,  Sds_Num_Level2_Cld_Type_Aux_Flag = 0, & !off for GOES
     Sds_Num_Level2_Ch20_Clear = 128,    Sds_Num_Level2_Ch20_Clear_Flag = 0, &
     Sds_Num_Level2_Ch27_Clear = 129,    Sds_Num_Level2_Ch27_Clear_Flag = 0, &
     Sds_Num_Level2_Ch28_Clear = 130,    Sds_Num_Level2_Ch28_Clear_Flag = 0, &
     Sds_Num_Level2_Ch29_Clear = 131,    Sds_Num_Level2_Ch29_Clear_Flag = 0, &
     Sds_Num_Level2_Ch30_Clear = 132,    Sds_Num_Level2_Ch30_Clear_Flag = 0, &
     Sds_Num_Level2_Ch31_Clear = 133,    Sds_Num_Level2_Ch31_Clear_Flag = 1, &
     Sds_Num_Level2_Ch32_Clear = 134,    Sds_Num_Level2_Ch32_Clear_Flag = 0, &
     Sds_Num_Level2_Ch33_Clear = 135,    Sds_Num_Level2_Ch33_Clear_Flag = 0, &
     Sds_Num_Level2_Glintzen = 136,      Sds_Num_Level2_Glintzen_Flag = 0, &
     Sds_Num_Level2_Scatzen = 137,       Sds_Num_Level2_Scatzen_Flag = 0, &
     Sds_Num_Level2_Bayes_Sfc_Type = 138,Sds_Num_Level2_Bayes_Sfc_Type_Flag = 0, &
     Sds_Num_Level2_Acha_Order = 139,    Sds_Num_Level2_Acha_Order_Flag = 0, &
     Sds_Num_Level2_Acha_Inver = 140,    Sds_Num_Level2_Acha_Inver_Flag = 0, &
     Sds_Num_Level2_Cth_Top = 141,       Sds_Num_Level2_Cth_Top_Flag = 0, &
     Sds_Num_Level2_Cth_Base = 142,      Sds_Num_Level2_Cth_Base_Flag = 0, &
     Sds_Num_Level2_Cld_Tests = 143,     Sds_Num_Level2_Cld_Tests_Flag = 0, &
     Sds_Num_Level2_Ref_Max_ChI1 = 144,  Sds_Num_Level2_Ref_Max_ChI1_Flag = 0, &
     Sds_Num_Level2_Ref_Min_ChI1 = 145,  Sds_Num_Level2_Ref_Min_ChI1_Flag = 0, &
     Sds_Num_Level2_Ref_Mean_ChI1 = 146, Sds_Num_Level2_Ref_Mean_ChI1_Flag = 0, &
     Sds_Num_Level2_Ref_Uni_ChI1 = 147,  Sds_Num_Level2_Ref_Uni_ChI1_Flag = 0, &
     Sds_Num_Level2_Ref_Max_ChI2 = 148,  Sds_Num_Level2_Ref_Max_ChI2_Flag = 0, &
     Sds_Num_Level2_Ref_Min_ChI2 = 149,  Sds_Num_Level2_Ref_Min_ChI2_Flag = 0, &
     Sds_Num_Level2_Ref_Mean_ChI2 = 150, Sds_Num_Level2_Ref_Mean_ChI2_Flag = 0, &
     Sds_Num_Level2_Ref_Uni_ChI2 = 151,  Sds_Num_Level2_Ref_Uni_ChI2_Flag = 0, &
     Sds_Num_Level2_Ref_Max_ChI3 = 152,  Sds_Num_Level2_Ref_Max_ChI3_Flag = 0, &                                                                       
     Sds_Num_Level2_Ref_Min_ChI3 = 153,  Sds_Num_Level2_Ref_Min_ChI3_Flag = 0, &                                                                       
     Sds_Num_Level2_Ref_Mean_ChI3 = 154, Sds_Num_Level2_Ref_Mean_ChI3_Flag = 0, &                                                                      
     Sds_Num_Level2_Ref_Uni_ChI3 = 155,  Sds_Num_Level2_Ref_Uni_ChI3_Flag = 0, &                                                                       
     Sds_Num_Level2_Bt_Max_ChI4 = 156,   Sds_Num_Level2_Bt_Max_ChI4_Flag = 0, &                                                                        
     Sds_Num_Level2_Bt_Min_ChI4 = 157,   Sds_Num_Level2_Bt_Min_ChI4_Flag = 0, &                                                                        
     Sds_Num_Level2_Bt_Mean_ChI4 = 158,  Sds_Num_Level2_Bt_Mean_ChI4_Flag = 0, &                                                                       
     Sds_Num_Level2_Bt_Uni_ChI4 = 159,   Sds_Num_Level2_Bt_Uni_ChI4_Flag = 0, &
     Sds_Num_Level2_Bt_Max_ChI5 = 160,   Sds_Num_Level2_Bt_Max_ChI5_Flag = 0, &
     Sds_Num_Level2_Bt_Min_ChI5 = 161,   Sds_Num_Level2_Bt_Min_ChI5_Flag = 0, &
     Sds_Num_Level2_Bt_Mean_ChI5 = 162,  Sds_Num_Level2_Bt_Mean_ChI5_Flag = 0, &
     Sds_Num_Level2_Bt_Uni_ChI5 = 163,   Sds_Num_Level2_Bt_Uni_ChI5_Flag = 0, &
     Sds_Num_Level2_Ndvi_White = 164,    Sds_Num_Level2_Ndvi_White_Flag = 0, &
     Sds_Num_Level2_Glint_Mask = 165,    Sds_Num_Level2_Glint_Mask_Flag = 1, &
     Sds_Num_Level2_Psfc_Back = 166,     Sds_Num_Level2_Psfc_Back_Flag = 0, &
     Sds_Num_Level2_Pmsl_Back = 167,     Sds_Num_Level2_Pmsl_Back_Flag = 0, &
     Sds_Num_Level2_Kindex = 168,        Sds_Num_Level2_Kindex_Flag = 0, &
     Sds_Num_Level2_Cwp_Nwp = 169,       Sds_Num_Level2_Cwp_Nwp_Flag = 1, &
     Sds_Num_Level2_Cfrac_Nwp = 170,     Sds_Num_Level2_Cfrac_Nwp_Flag = 1, &
     Sds_Num_Level2_Pc_Nwp = 171,        Sds_Num_Level2_Pc_Nwp_Flag = 0, &
     Sds_Num_Level2_Ncld_Nwp = 172,      Sds_Num_Level2_Ncld_Nwp_Flag = 0, &
     Sds_Num_Level2_Cld_Type_Nwp = 173,  Sds_Num_Level2_Cld_Type_Nwp_Flag = 0, &
     Sds_Num_Level2_Rain_Rate = 174,     Sds_Num_Level2_Rain_Rate_Flag = 0, &
     Sds_Num_Level2_ChDNB_Lunar = 175,   Sds_Num_Level2_ChDNB_Lunar_Flag = 0, &
     Sds_Num_Level2_Adj_Pix_Cld_Mask = 176, Sds_Num_Level2_Adj_Pix_Cld_Mask_Flag = 0, &
     Sds_Num_Level2_Dcomp_Quality = 177, Sds_Num_Level2_Dcomp_Quality_Flag = 1, &
     Sds_Num_Level2_Dcomp_Info = 178,    Sds_Num_Level2_Dcomp_Info_Flag = 1, &
     Sds_Num_Level2_Acha_Quality = 179,  Sds_Num_Level2_Acha_Quality_Flag = 1, &
     Sds_Num_Level2_Acha_Info = 180,     Sds_Num_Level2_Acha_Info_Flag = 1, &
     Sds_Num_Level2_Diag1 = 181,         Sds_Num_Level2_Diag1_Flag = 0, &
     Sds_Num_Level2_Diag2 = 182,         Sds_Num_Level2_Diag2_Flag = 0, &
     Sds_Num_Level2_Diag3 = 183,         Sds_Num_Level2_Diag3_Flag = 0, &
     Sds_Num_Level2_Cod_Nlcomp = 184,    Sds_Num_Level2_Cod_Nlcomp_Flag = 0, &
     Sds_Num_Level2_Ceps_Nlcomp = 185,   Sds_Num_Level2_Ceps_Nlcomp_Flag = 0, &
     Sds_Num_Level2_Cod_Nlcomp_Uncer = 186,  Sds_Num_Level2_Cod_Nlcomp_Uncer_Flag = 0, &
     Sds_Num_Level2_Ceps_Nlcomp_Uncer = 187, Sds_Num_Level2_Ceps_Nlcomp_Uncer_Flag = 0, &
     Sds_Num_Level2_Nlcomp_Quality = 188,  Sds_Num_Level2_Nlcomp_Quality_Flag = 0, &
     Sds_Num_Level2_Nlcomp_Info = 189,    Sds_Num_Level2_Nlcomp_Info_Flag = 0, &
     Sds_Num_Level2_Ctt_Acha_Uncer = 190, Sds_Num_Level2_Ctt_Acha_Uncer_Flag = 0, &
     Sds_Num_Level2_Lunzen = 191,         Sds_Num_Level2_Lunzen_Flag = 0, &
     Sds_Num_Level2_LunRelaz = 192,       Sds_Num_Level2_LunRelaz_Flag = 0,&
     Sds_Num_Level2_Lunaz = 193,          Sds_Num_Level2_Lunaz_Flag = 0, &
     Sds_Num_Level2_Latpc = 194,          Sds_Num_Level2_Latpc_Flag = 0, &
     Sds_Num_Level2_Lonpc = 195,          Sds_Num_Level2_Lonpc_Flag = 0, &
     Sds_Num_Level2_Sataz = 196,          Sds_Num_Level2_Sataz_Flag = 0, &
     Sds_Num_Level2_CldInsol = 197,       Sds_Num_Level2_CldInsol_Flag = 1, &
     Sds_Num_Level2_CldInsol_Dif = 198,   Sds_Num_Level2_CldInsol_Dif_Flag = 1, &
     Sds_Num_Level2_LCL_Nwp = 199,        Sds_Num_Level2_LCL_Nwp_Flag = 0,   &
     Sds_Num_Level2_CCL_Nwp = 200,        Sds_Num_Level2_CCL_Nwp_Flag = 0, &
     Sds_Num_Level2_Alt = 201,            Sds_Num_Level2_Alt_Flag = 0, &
     Sds_Num_Level2_Ch31_Rad_Atm = 202,   Sds_Num_Level2_Ch31_Rad_Atm_Flag = 0, &
     Sds_Num_Level2_Ch31_Trans_Atm = 203, Sds_Num_Level2_Ch31_Trans_Atm_Flag = 0, &
     Sds_Num_Level2_Ch31_Sfc_Emiss = 204, Sds_Num_Level2_Ch31_Sfc_Emiss_Flag = 0, &
     Sds_Num_Level2_Ch20_Sfc_Emiss = 205, Sds_Num_Level2_Ch20_Sfc_Emiss_Flag = 0, &
     Sds_Num_Level2_Ch31_Rad_Atm_Dwn = 206, Sds_Num_Level2_Ch31_Rad_Atm_Dwn_Flag = 0, &
     Sds_Num_Level2_Ozone = 207,          Sds_Num_Level2_Ozone_Flag = 0, &
     Sds_Num_Level2_Acha_Cost = 208,      Sds_Num_Level2_Acha_Cost_Flag = 0, &
     Sds_Num_Level2_Hcld_Dcomp = 209,     Sds_Num_Level2_Hcld_Dcomp_Flag = 0, &
     Sds_Num_Level2_Cdnc_Dcomp = 210,     Sds_Num_Level2_Cdnc_Dcomp_Flag = 0, &
     Sds_Num_Level2_Bt375_Snd = 211,      Sds_Num_Level2_Bt375_Snd_Flag = 0, &
     Sds_Num_Level2_Bt11_Snd = 212,       Sds_Num_Level2_Bt11_Snd_Flag = 0, &
     Sds_Num_Level2_Bt12_Snd = 213,       Sds_Num_Level2_Bt12_Snd_Flag = 0, &
     Sds_Num_Level2_High_Cld = 214,       Sds_Num_Level2_High_Cld_Flag = 0, &
     Sds_Num_Level2_Mid_Cld = 215,        Sds_Num_Level2_Mid_Cld_Flag = 0, &
     Sds_Num_Level2_Low_Cld = 216,        Sds_Num_Level2_Low_Cld_Flag = 0, &
     Sds_Num_Level2_Dust = 217,           Sds_Num_Level2_Dust_Flag = 0, &
     Sds_Num_Level2_Smoke = 218,          Sds_Num_Level2_Smoke_Flag = 0, &
     Sds_Num_Level2_Shadow = 219,         Sds_Num_Level2_Shadow_Flag = 0, &
     Sds_Num_Level2_Fire = 220,           Sds_Num_Level2_Fire_Flag = 0, & 
     Sds_Num_Level2_Zc_Lower = 221,       Sds_Num_Level2_Zc_Lower_Flag = 0, &
     Sds_Num_Level2_Pc_Lower = 222,       Sds_Num_Level2_Pc_Lower_Flag = 0, &
     Sds_Num_Level2_Cod_Mask = 223,       Sds_Num_Level2_Cod_Mask_Flag = 0, &
     Sds_Num_Level2_Ch37 = 224,           Sds_Num_Level2_Ch37_Flag = 0, &
     Sds_Num_Level2_Ch38 = 225,           Sds_Num_Level2_Ch38_Flag = 0, &
     Sds_Num_Level2_Ch37_Clear = 226,     Sds_Num_Level2_Ch37_Clear_Flag = 0, &
     Sds_Num_Level2_Ch38_Clear = 227,     Sds_Num_Level2_Ch38_Clear_Flag = 0, &
     Sds_Num_Level2_Cth_Aux = 228,        Sds_Num_Level2_Cth_Aux_Flag = 0, &
     Sds_Num_Level2_Cth_Sndr = 229,       Sds_Num_Level2_Cth_Sndr_Flag = 0, &
     Sds_Num_Level2_Cld_Layer = 230,      Sds_Num_Level2_Cld_Layer_Flag = 0, &
     Sds_Num_Level2_Conv_Mask = 231,      Sds_Num_Level2_Conv_Mask_Flag = 0, &
     Sds_Num_Level2_Beta_11_67 = 232,     Sds_Num_Level2_Beta_11_67_Flag = 0, &
     Sds_Num_Level2_Beta_11_85 = 233,     Sds_Num_Level2_Beta_11_85_Flag = 0, &
     Sds_Num_Level2_Beta_11_12 = 234,     Sds_Num_Level2_Beta_11_12_Flag = 0, &
     Sds_Num_Level2_Beta_11_13 = 235,     Sds_Num_Level2_Beta_11_13_Flag = 0, &
     Sds_Num_Level2_Ctp_Top = 236,        Sds_Num_Level2_Ctp_Top_Flag = 0, &
     Sds_Num_Level2_Ctp_Base = 237,       Sds_Num_Level2_Ctp_Base_Flag = 0
