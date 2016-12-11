!$Id$
 !--- level2 indices
 integer, parameter, private:: Num_Level2_Sds = 290
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
     Sds_Num_Level2_Bad_Scan = 3,        Sds_Num_Level2_Bad_Scan_Flag = 1, &
     Sds_Num_Level2_Bad_Pixel_Mask = 4,  Sds_Num_Level2_Bad_Pixel_Mask_Flag = 1, &  !critical for level2b
     Sds_Num_Level2_Gap_Pixel_Mask = 5,  Sds_Num_Level2_Gap_Pixel_Mask_Flag = 1, &  !critical for level2b for VIIRS
     Sds_Num_Level2_Asc_Flag = 6,        Sds_Num_Level2_Asc_Flag_Flag = 1, &        !critical for level2b (except GOES)
     Sds_Num_Level2_Meta_Data = 7,       Sds_Num_Level2_Meta_Data_Flag = 1, &
     Sds_Num_Level2_Lat = 8,             Sds_Num_Level2_Lat_Flag = 1, &         !critical for level2b
     Sds_Num_Level2_Lon = 9,             Sds_Num_Level2_Lon_Flag = 1, &         !critical for level2b
     Sds_Num_Level2_Zen = 10,            Sds_Num_Level2_Zen_Flag = 1, &         !critical for level2b
     Sds_Num_Level2_Solzen = 11,         Sds_Num_Level2_Solzen_Flag = 1, &      !critical for level2b
     Sds_Num_Level2_Relaz = 12 ,         Sds_Num_Level2_Relaz_Flag = 1, &       
     Sds_Num_Level2_Solaz = 13 ,         Sds_Num_Level2_Solaz_Flag = 1, &       
     Sds_Num_Level2_Packed_Land = 14,    Sds_Num_Level2_Packed_Land_Flag = 1, &
     Sds_Num_Level2_Ch1 = 15 ,           Sds_Num_Level2_Ch1_Flag = 1, &         !critical for level2b
     Sds_Num_Level2_Ch2 = 16 ,           Sds_Num_Level2_Ch2_Flag = 0, &
     Sds_Num_Level2_Ch3 = 17 ,           Sds_Num_Level2_Ch3_Flag = 0, &
     Sds_Num_Level2_Ch4 = 18 ,           Sds_Num_Level2_Ch4_Flag = 0, &
     Sds_Num_Level2_Ch5 = 19 ,           Sds_Num_Level2_Ch5_Flag = 0, &
     Sds_Num_Level2_Ch6 = 20 ,           Sds_Num_Level2_Ch6_Flag = 0, &
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
     Sds_Num_Level2_Ch26 = 41 ,          Sds_Num_Level2_Ch26_Flag = 1, &
     Sds_Num_Level2_Ch27 = 42 ,          Sds_Num_Level2_Ch27_Flag = 1, &
     Sds_Num_Level2_Ch28 = 43 ,          Sds_Num_Level2_Ch28_Flag = 1, &
     Sds_Num_Level2_Ch29 = 44 ,          Sds_Num_Level2_Ch29_Flag = 1, &
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
     Sds_Num_Level2_Cldprob = 56,        Sds_Num_Level2_Cldprob_Flag = 0, &
     Sds_Num_Level2_Cld_Type = 57 ,      Sds_Num_Level2_Cld_Type_Flag = 0, &
     Sds_Num_Level2_Cld_Phase = 58 ,     Sds_Num_Level2_Cld_Phase_Flag = 0, &
     Sds_Num_Level2_Ctp = 59,            Sds_Num_Level2_Ctp_Flag = 0, &
     Sds_Num_Level2_Ctt = 60,            Sds_Num_Level2_Ctt_Flag = 0, &
     Sds_Num_Level2_Cth = 61,            Sds_Num_Level2_Cth_Flag = 0, &
     Sds_Num_Level2_Ec = 62,             Sds_Num_Level2_Ec_Flag = 0, &
     Sds_Num_Level2_Beta = 63,           Sds_Num_Level2_Beta_Flag = 0, &
     Sds_Num_Level2_Cod_Acha = 64,       Sds_Num_Level2_Cod_Acha_Flag = 0, &
     Sds_Num_Level2_Ceps_Acha = 65,      Sds_Num_Level2_Ceps_Acha_Flag = 0, &
     Sds_Num_Level2_Cod = 66,            Sds_Num_Level2_Cod_Flag = 0, &
     Sds_Num_Level2_Ceps = 67,           Sds_Num_Level2_Ceps_Flag = 0, &
     Sds_Num_Level2_Cldalb = 68,         Sds_Num_Level2_Cldalb_Flag = 0, &
     Sds_Num_Level2_Cldtrn = 69,         Sds_Num_Level2_Cldtrn_Flag = 0, &
     Sds_Num_Level2_Cldfrac = 70,        Sds_Num_Level2_Cldfrac_Flag = 0, &
     Sds_Num_Level2_Cldfrac_Uncer = 71,  Sds_Num_Level2_Cldfrac_Uncer_Flag = 0, &
     Sds_Num_Level2_Etrop = 72,          Sds_Num_Level2_Etrop_Flag = 1, &
     Sds_Num_Level2_Aot1 = 73,           Sds_Num_Level2_Aot1_Flag = 0, & !off for GOES
     Sds_Num_Level2_Olr = 74,            Sds_Num_Level2_Olr_Flag = 0, & !off for GOES
     Sds_Num_Level2_Ndvi = 75,           Sds_Num_Level2_Ndvi_Flag = 0, & !off for GOES
     Sds_Num_Level2_Tsfc = 76,           Sds_Num_Level2_Tsfc_Flag = 1, &
     Sds_Num_Level2_Tair = 77,           Sds_Num_Level2_Tair_Flag = 0, &
     Sds_Num_Level2_Trad = 78,           Sds_Num_Level2_Trad_Flag = 0, &
     Sds_Num_Level2_Rsr = 79,            Sds_Num_Level2_Rsr_Flag = 0, &
     Sds_Num_Level2_Qf1 = 80,            Sds_Num_Level2_Qf1_Flag = 0, &		!deprecated
     Sds_Num_Level2_Qf2 = 81,            Sds_Num_Level2_Qf2_Flag = 0, &		!deprecated
     Sds_Num_Level2_Cod_Dcomp_Uncer = 82,Sds_Num_Level2_Cod_Dcomp_Uncer_Flag = 1, &
     Sds_Num_Level2_Ceps_Dcomp_Uncer = 83,Sds_Num_Level2_Ceps_Dcomp_Uncer_Flag = 0, &
     Sds_Num_Level2_Cth_Acha_Uncer = 84, Sds_Num_Level2_Cth_Acha_Uncer_Flag = 0, &
     Sds_Num_Level2_Insol = 85,          Sds_Num_Level2_Insol_Flag = 0, &
     Sds_Num_Level2_Insol_Dif = 86,      Sds_Num_Level2_Insol_Dif_Flag = 0, &
     Sds_Num_Level2_Ch1_Counts = 87,     Sds_Num_Level2_Ch1_Counts_Flag = 0, &
     Sds_Num_Level2_Ch2_Counts = 88,     Sds_Num_Level2_Ch2_Counts_Flag = 0, &
     Sds_Num_Level2_Ch6_Counts = 89,     Sds_Num_Level2_Ch6_Counts_Flag = 0, &
     Sds_Num_Level2_Tpw = 90,            Sds_Num_Level2_Tpw_Flag = 1, &
     Sds_Num_Level2_Rh = 91,             Sds_Num_Level2_Rh_Flag = 0, &
     Sds_Num_Level2_Ref_Ch20_Sfc = 92,   Sds_Num_Level2_Ref_Ch20_Sfc_Flag = 0, &
     Sds_Num_Level2_Ref_Ch1_Sfc = 93,    Sds_Num_Level2_Ref_Ch1_Sfc_Flag = 0, &
     Sds_Num_Level2_Ref_Ch2_Sfc = 94,    Sds_Num_Level2_Ref_Ch2_Sfc_Flag = 0, &
     Sds_Num_Level2_Ref_Ch6_Sfc = 95,    Sds_Num_Level2_Ref_Ch6_Sfc_Flag = 0, &
     Sds_Num_Level2_Sst_Masked = 96,     Sds_Num_Level2_Sst_Masked_Flag = 0, &
     Sds_Num_Level2_Ch1_Unnorm = 97,     Sds_Num_Level2_Ch1_Unnorm_Flag = 0, &
     Sds_Num_Level2_Ch2_Unnorm = 98,     Sds_Num_Level2_Ch2_Unnorm_Flag = 0, &
     Sds_Num_Level2_Ch6_Unnorm = 99,     Sds_Num_Level2_Ch6_Unnorm_Flag = 0, &
     Sds_Num_Level2_Ch1_Clear = 100,     Sds_Num_Level2_Ch1_Clear_Flag = 1, &
     Sds_Num_Level2_Ch1_Mean = 101,      Sds_Num_Level2_Ch1_Mean_Flag = 1, &
     Sds_Num_Level2_Sst_Unmasked = 102,  Sds_Num_Level2_Sst_Unmasked_Flag = 1, &
     Sds_Num_Level2_Wnd_Spd = 103,       Sds_Num_Level2_Wnd_Spd_Flag = 0, &
     Sds_Num_Level2_Wnd_Dir = 104,       Sds_Num_Level2_Wnd_Dir_Flag = 0, &
     Sds_Num_Level2_Ch1_Dark = 105,      Sds_Num_Level2_Ch1_Dark_Flag = 0, &
     Sds_Num_Level2_Cwp = 106,           Sds_Num_Level2_Cwp_Flag = 0, &
     Sds_Num_Level2_Wnd_Spd_Cld_Top= 107, Sds_Num_Level2_Wnd_Spd_Cld_Top_Flag = 0, &
     Sds_Num_Level2_Wnd_Dir_Cld_Top= 108, Sds_Num_Level2_Wnd_Dir_Cld_Top_Flag = 0, &
     Sds_Num_Level2_Temp_Tropo = 109,    Sds_Num_Level2_Temp_Tropo_Flag = 1, &
     Sds_Num_Level2_Cld_Mask =  110,     Sds_Num_Level2_Cld_Mask_Flag = 1, &
     Sds_Num_Level2_Sfc_Type = 111,      Sds_Num_Level2_Sfc_Type_Flag = 1, &
     Sds_Num_Level2_Coast_Mask = 112,    Sds_Num_Level2_Coast_Mask_Flag = 1, &
     Sds_Num_Level2_Land_Mask = 113,     Sds_Num_Level2_Land_Mask_Flag = 1, &
     Sds_Num_Level2_Snow_Mask = 114,     Sds_Num_Level2_Snow_Mask_Flag = 1, &
     Sds_Num_Level2_Cod_Dcomp_Qf = 115,  Sds_Num_Level2_Cod_Dcomp_Qf_Flag = 0, &
     Sds_Num_Level2_Ceps_Dcomp_Qf = 116, Sds_Num_Level2_Ceps_Dcomp_Qf_Flag = 0, &
     Sds_Num_Level2_Cth_Acha_Qf = 117,   Sds_Num_Level2_Cth_Acha_Qf_Flag = 0, &
     Sds_Num_Level2_Ec_Acha_Qf = 118,    Sds_Num_Level2_Ec_Acha_Qf_Flag = 0, &
     Sds_Num_Level2_Beta_Acha_Qf = 119,  Sds_Num_Level2_Beta_Acha_Qf_Flag = 0, &
     Sds_Num_Level2_Aot2 = 120,          Sds_Num_Level2_Aot2_Flag = 0, &
     Sds_Num_Level2_Aot6 = 121,          Sds_Num_Level2_Aot6_Flag = 0, &
     Sds_Num_Level2_Aot_QF = 122,        Sds_Num_Level2_Aot_QF_Flag = 0, &
     Sds_Num_Level2_Zsfc = 123,          Sds_Num_Level2_Zsfc_Flag = 0, &
     Sds_Num_Level2_Cth_H2O = 124,       Sds_Num_Level2_Cth_H2O_Flag = 0, &
     Sds_Num_Level2_Cth_Opa = 125,       Sds_Num_Level2_Cth_Opa_Flag = 0, &
     Sds_Num_Level2_Cld_Mask_Aux = 126,  Sds_Num_Level2_Cld_Mask_Aux_Flag = 0, & !off for GOES
     Sds_Num_Level2_Cld_Type_Aux = 127,  Sds_Num_Level2_Cld_Type_Aux_Flag = 0, & !off for GOES
     Sds_Num_Level2_Cld_Phase_Aux = 128, Sds_Num_Level2_Cld_Phase_Aux_Flag = 0, & !off for GOES
     Sds_Num_Level2_Ch20_Clear = 129,    Sds_Num_Level2_Ch20_Clear_Flag = 1, &
     Sds_Num_Level2_Ch27_Clear = 130,    Sds_Num_Level2_Ch27_Clear_Flag = 1, &
     Sds_Num_Level2_Ch28_Clear = 131,    Sds_Num_Level2_Ch28_Clear_Flag = 0, &
     Sds_Num_Level2_Ch29_Clear = 132,    Sds_Num_Level2_Ch29_Clear_Flag = 0, &
     Sds_Num_Level2_Ch30_Clear = 133,    Sds_Num_Level2_Ch30_Clear_Flag = 0, &
     Sds_Num_Level2_Ch31_Clear = 134,    Sds_Num_Level2_Ch31_Clear_Flag = 1, &
     Sds_Num_Level2_Ch32_Clear = 135,    Sds_Num_Level2_Ch32_Clear_Flag = 0, &
     Sds_Num_Level2_Ch33_Clear = 136,    Sds_Num_Level2_Ch33_Clear_Flag = 0, &
     Sds_Num_Level2_Glintzen = 137,      Sds_Num_Level2_Glintzen_Flag = 1, &
     Sds_Num_Level2_Scatzen = 138,       Sds_Num_Level2_Scatzen_Flag = 1, &
     Sds_Num_Level2_Bayes_Sfc_Type = 139,Sds_Num_Level2_Bayes_Sfc_Type_Flag = 0, &
     Sds_Num_Level2_Acha_Order = 140,    Sds_Num_Level2_Acha_Order_Flag = 0, &
     Sds_Num_Level2_Acha_Inver = 141,    Sds_Num_Level2_Acha_Inver_Flag = 0, &
     Sds_Num_Level2_Cth_Top = 142,       Sds_Num_Level2_Cth_Top_Flag = 0, &
     Sds_Num_Level2_Cth_Base = 143,      Sds_Num_Level2_Cth_Base_Flag = 0, &
     Sds_Num_Level2_Cld_Tests = 144,     Sds_Num_Level2_Cld_Tests_Flag = 1, &
     Sds_Num_Level2_Ref_Max_ChI1 = 145,  Sds_Num_Level2_Ref_Max_ChI1_Flag = 0, &
     Sds_Num_Level2_Ref_Min_ChI1 = 146,  Sds_Num_Level2_Ref_Min_ChI1_Flag = 0, &
     Sds_Num_Level2_Ref_Mean_ChI1 = 147, Sds_Num_Level2_Ref_Mean_ChI1_Flag = 0, &
     Sds_Num_Level2_Ref_Uni_ChI1 = 148,  Sds_Num_Level2_Ref_Uni_ChI1_Flag = 0, &
     Sds_Num_Level2_Ref_Max_ChI2 = 149,  Sds_Num_Level2_Ref_Max_ChI2_Flag = 0, &
     Sds_Num_Level2_Ref_Min_ChI2 = 150,  Sds_Num_Level2_Ref_Min_ChI2_Flag = 0, &
     Sds_Num_Level2_Ref_Mean_ChI2 = 151, Sds_Num_Level2_Ref_Mean_ChI2_Flag = 0, &
     Sds_Num_Level2_Ref_Uni_ChI2 = 152,  Sds_Num_Level2_Ref_Uni_ChI2_Flag = 0, &
     Sds_Num_Level2_Ref_Max_ChI3 = 153,  Sds_Num_Level2_Ref_Max_ChI3_Flag = 0, &
     Sds_Num_Level2_Ref_Min_ChI3 = 154,  Sds_Num_Level2_Ref_Min_ChI3_Flag = 0, &
     Sds_Num_Level2_Ref_Mean_ChI3 = 155, Sds_Num_Level2_Ref_Mean_ChI3_Flag = 0, &
     Sds_Num_Level2_Ref_Uni_ChI3 = 156,  Sds_Num_Level2_Ref_Uni_ChI3_Flag = 0, &
     Sds_Num_Level2_Bt_Max_ChI4 = 157,   Sds_Num_Level2_Bt_Max_ChI4_Flag = 0, &
     Sds_Num_Level2_Bt_Min_ChI4 = 158,   Sds_Num_Level2_Bt_Min_ChI4_Flag = 0, &
     Sds_Num_Level2_Bt_Mean_ChI4 = 159,  Sds_Num_Level2_Bt_Mean_ChI4_Flag = 0, &
     Sds_Num_Level2_Bt_Uni_ChI4 = 160,   Sds_Num_Level2_Bt_Uni_ChI4_Flag = 0, &
     Sds_Num_Level2_Bt_Max_ChI5 = 161,   Sds_Num_Level2_Bt_Max_ChI5_Flag = 0, &
     Sds_Num_Level2_Bt_Min_ChI5 = 162,   Sds_Num_Level2_Bt_Min_ChI5_Flag = 0, &
     Sds_Num_Level2_Bt_Mean_ChI5 = 163,  Sds_Num_Level2_Bt_Mean_ChI5_Flag = 0, &
     Sds_Num_Level2_Bt_Uni_ChI5 = 164,   Sds_Num_Level2_Bt_Uni_ChI5_Flag = 0, &
     Sds_Num_Level2_Ndvi_White = 165,    Sds_Num_Level2_Ndvi_White_Flag = 0, &
     Sds_Num_Level2_Glint_Mask = 166,    Sds_Num_Level2_Glint_Mask_Flag = 1, &
     Sds_Num_Level2_Psfc_Back = 167,     Sds_Num_Level2_Psfc_Back_Flag = 0, &
     Sds_Num_Level2_Pmsl_Back = 168,     Sds_Num_Level2_Pmsl_Back_Flag = 0, &
     Sds_Num_Level2_Kindex = 169,        Sds_Num_Level2_Kindex_Flag = 0, &
     Sds_Num_Level2_Cwp_Nwp = 170,       Sds_Num_Level2_Cwp_Nwp_Flag = 0, &
     Sds_Num_Level2_Cfrac_Nwp = 171,     Sds_Num_Level2_Cfrac_Nwp_Flag = 0, &
     Sds_Num_Level2_Pc_Nwp = 172,        Sds_Num_Level2_Pc_Nwp_Flag = 0, &
     Sds_Num_Level2_Ncld_Nwp = 173,      Sds_Num_Level2_Ncld_Nwp_Flag = 0, &
     Sds_Num_Level2_Cld_Type_Nwp = 174,  Sds_Num_Level2_Cld_Type_Nwp_Flag = 0, &
     Sds_Num_Level2_Rain_Rate = 175,     Sds_Num_Level2_Rain_Rate_Flag = 0, &
     Sds_Num_Level2_ChDNB_Lunar = 176,   Sds_Num_Level2_ChDNB_Lunar_Flag = 0, &
     Sds_Num_Level2_Adj_Pix_Cld_Mask = 177, Sds_Num_Level2_Adj_Pix_Cld_Mask_Flag = 0, &
     Sds_Num_Level2_Dcomp_Quality = 178,  Sds_Num_Level2_Dcomp_Quality_Flag = 0, &
     Sds_Num_Level2_Dcomp_Info = 179,    Sds_Num_Level2_Dcomp_Info_Flag = 0, &
     Sds_Num_Level2_Acha_Quality = 180,  Sds_Num_Level2_Acha_Quality_Flag = 0, &
     Sds_Num_Level2_Acha_Info = 181,     Sds_Num_Level2_Acha_Info_Flag = 0, &
     Sds_Num_Level2_Diag1 = 182,         Sds_Num_Level2_Diag1_Flag = 0, &
     Sds_Num_Level2_Diag2 = 183,         Sds_Num_Level2_Diag2_Flag = 0, &
     Sds_Num_Level2_Diag3 = 184,         Sds_Num_Level2_Diag3_Flag = 0, &
     Sds_Num_Level2_Cod_Nlcomp = 185,    Sds_Num_Level2_Cod_Nlcomp_Flag = 0, &
     Sds_Num_Level2_Ceps_Nlcomp = 186,   Sds_Num_Level2_Ceps_Nlcomp_Flag = 0, &
     Sds_Num_Level2_Cod_Nlcomp_Uncer = 187,  Sds_Num_Level2_Cod_Nlcomp_Uncer_Flag = 0, &
     Sds_Num_Level2_Ceps_Nlcomp_Uncer = 188, Sds_Num_Level2_Ceps_Nlcomp_Uncer_Flag = 0, &
     Sds_Num_Level2_Nlcomp_Quality = 189, Sds_Num_Level2_Nlcomp_Quality_Flag = 0, &
     Sds_Num_Level2_Nlcomp_Info = 190,    Sds_Num_Level2_Nlcomp_Info_Flag = 0, &
     Sds_Num_Level2_Ctt_Acha_Uncer = 191, Sds_Num_Level2_Ctt_Acha_Uncer_Flag = 1, &
     Sds_Num_Level2_Lunzen = 192,         Sds_Num_Level2_Lunzen_Flag = 0, &   
     Sds_Num_Level2_LunRelaz = 193,       Sds_Num_Level2_LunRelaz_Flag = 0,&       
     Sds_Num_Level2_Lunaz = 194,          Sds_Num_Level2_Lunaz_Flag = 0, &
     Sds_Num_Level2_Latpc = 195,          Sds_Num_Level2_Latpc_Flag = 0, &
     Sds_Num_Level2_Lonpc = 196,          Sds_Num_Level2_Lonpc_Flag = 0, &
     Sds_Num_Level2_Sataz = 197,          Sds_Num_Level2_Sataz_Flag = 1, &
     Sds_Num_Level2_CldInsol = 198,       Sds_Num_Level2_CldInsol_Flag = 0, &
     Sds_Num_Level2_CldInsol_Dif = 199,   Sds_Num_Level2_CldInsol_Dif_Flag = 0, &
     Sds_Num_Level2_LCL_Nwp = 200,        Sds_Num_Level2_LCL_Nwp_Flag = 0,   &
     Sds_Num_Level2_CCL_Nwp = 201,        Sds_Num_Level2_CCL_Nwp_Flag = 0, &
     Sds_Num_Level2_Alt = 202,            Sds_Num_Level2_Alt_Flag = 0, &
     Sds_Num_Level2_Ch31_Rad_Atm = 203,   Sds_Num_Level2_Ch31_Rad_Atm_Flag = 0, &
     Sds_Num_Level2_Ch31_Trans_Atm = 204, Sds_Num_Level2_Ch31_Trans_Atm_Flag = 0, &
     Sds_Num_Level2_Ch31_Sfc_Emiss = 205, Sds_Num_Level2_Ch31_Sfc_Emiss_Flag = 0, &
     Sds_Num_Level2_Ch20_Sfc_Emiss = 206, Sds_Num_Level2_Ch20_Sfc_Emiss_Flag = 0, &
     Sds_Num_Level2_Ch31_Rad_Atm_Dwn = 207, Sds_Num_Level2_Ch31_Rad_Atm_Dwn_Flag = 0, &
     Sds_Num_Level2_Ozone = 208,          Sds_Num_Level2_Ozone_Flag = 0, &
     Sds_Num_Level2_Acha_Cost = 209,      Sds_Num_Level2_Acha_Cost_Flag = 0, &
     Sds_Num_Level2_Aux_Cost = 210,       Sds_Num_Level2_Aux_Cost_Flag = 0, &
     Sds_Num_Level2_Hcld_Dcomp = 211,     Sds_Num_Level2_Hcld_Dcomp_Flag = 0, &
     Sds_Num_Level2_Cdnc_Dcomp = 212,     Sds_Num_Level2_Cdnc_Dcomp_Flag = 0, &
     Sds_Num_Level2_Bt375_Snd = 213,      Sds_Num_Level2_Bt375_Snd_Flag = 0, &
     Sds_Num_Level2_Bt11_Snd = 214,       Sds_Num_Level2_Bt11_Snd_Flag = 0, &
     Sds_Num_Level2_Bt12_Snd = 215,       Sds_Num_Level2_Bt12_Snd_Flag = 0, &
     Sds_Num_Level2_High_Cld = 216,       Sds_Num_Level2_High_Cld_Flag = 0, &
     Sds_Num_Level2_Mid_Cld = 217,        Sds_Num_Level2_Mid_Cld_Flag = 0, &
     Sds_Num_Level2_Low_Cld = 218,        Sds_Num_Level2_Low_Cld_Flag = 0, &
     Sds_Num_Level2_Dust = 219,           Sds_Num_Level2_Dust_Flag = 0, &
     Sds_Num_Level2_Smoke = 220,          Sds_Num_Level2_Smoke_Flag = 0, &
     Sds_Num_Level2_Shadow = 221,         Sds_Num_Level2_Shadow_Flag = 0, &
     Sds_Num_Level2_Fire = 222,           Sds_Num_Level2_Fire_Flag = 0, &
     Sds_Num_Level2_Zc_Lower = 223,       Sds_Num_Level2_Zc_Lower_Flag = 0, &
     Sds_Num_Level2_Pc_Lower = 224,       Sds_Num_Level2_Pc_Lower_Flag = 0, &
     Sds_Num_Level2_Cod_Mask = 225,       Sds_Num_Level2_Cod_Mask_Flag = 0, &
     Sds_Num_Level2_Cod_Aux = 226,        Sds_Num_Level2_Cod_Aux_Flag = 0, &
     Sds_Num_Level2_Ch37 = 227,           Sds_Num_Level2_Ch37_Flag = 0, &
     Sds_Num_Level2_Ch38 = 228,           Sds_Num_Level2_Ch38_Flag = 0, &
     Sds_Num_Level2_Ch45 = 229,           Sds_Num_Level2_Ch45_Flag = 0, &
     Sds_Num_Level2_Ch37_Clear = 230,     Sds_Num_Level2_Ch37_Clear_Flag = 0, &
     Sds_Num_Level2_Ch38_Clear = 231,     Sds_Num_Level2_Ch38_Clear_Flag = 0, &
     Sds_Num_Level2_Ch45_Clear = 232,     Sds_Num_Level2_Ch45_Clear_Flag = 0, &
     Sds_Num_Level2_Cth_Aux = 233,        Sds_Num_Level2_Cth_Aux_Flag = 0, &
     Sds_Num_Level2_Cth_Sndr = 234,       Sds_Num_Level2_Cth_Sndr_Flag = 0, &
     Sds_Num_Level2_Cld_Layer = 235,      Sds_Num_Level2_Cld_Layer_Flag = 0, &
     Sds_Num_Level2_Conv_Prob = 236,      Sds_Num_Level2_Conv_Prob_Flag = 0, &
     Sds_Num_Level2_Beta_11_67 = 237,     Sds_Num_Level2_Beta_11_67_Flag = 0, &
     Sds_Num_Level2_Beta_11_85 = 238,     Sds_Num_Level2_Beta_11_85_Flag = 0, &
     Sds_Num_Level2_Beta_11_12 = 239,     Sds_Num_Level2_Beta_11_12_Flag = 0, &
     Sds_Num_Level2_Beta_11_13 = 240,     Sds_Num_Level2_Beta_11_13_Flag = 0, &
     Sds_Num_Level2_Ctp_Top = 241,        Sds_Num_Level2_Ctp_Top_Flag = 0, &
     Sds_Num_Level2_Ctp_Top1_Aux = 242,   Sds_Num_Level2_Ctp_Top1_Aux_Flag = 0, &
     Sds_Num_Level2_Ctp_Top2_Aux = 243,   Sds_Num_Level2_Ctp_Top2_Aux_Flag = 0, &
     Sds_Num_Level2_Ctp_Base = 244,       Sds_Num_Level2_Ctp_Base_Flag = 0, &
     Sds_Num_Level2_Supercool_Prob = 245, Sds_Num_Level2_Supercool_Prob_Flag = 0, &
     Sds_Num_Level2_Ctp_Acha_Uncer = 246, Sds_Num_Level2_Ctp_Acha_Uncer_Flag = 0, &
     Sds_Num_Level2_Ctp_Aux_Uncer1 = 247, Sds_Num_Level2_Ctp_Aux_Uncer1_Flag = 0, &
     Sds_Num_Level2_Ctp_Aux_Uncer2 = 248, Sds_Num_Level2_Ctp_Aux_Uncer2_Flag = 0, &
     Sds_Num_Level2_Ch34_Clear = 249,     Sds_Num_Level2_Ch34_Clear_Flag = 0, &
     Sds_Num_Level2_Ch35_Clear = 250,     Sds_Num_Level2_Ch35_Clear_Flag = 0, &
     Sds_Num_Level2_Ch36_Clear = 251,     Sds_Num_Level2_Ch36_Clear_Flag = 0, &
     Sds_Num_Level2_Ch27_CSBT = 252,      Sds_Num_Level2_Ch27_CSBT_Flag = 0, &
     Sds_Num_Level2_Ch28_CSBT = 253,      Sds_Num_Level2_Ch28_CSBT_Flag = 0, &
     Sds_Num_Level2_Ch37_CSBT = 254,      Sds_Num_Level2_Ch37_CSBT_Flag = 0, &
     Sds_Num_Level2_Ch33_CSBT = 255,      Sds_Num_Level2_Ch33_CSBT_Flag = 0, &
     Sds_Num_Level2_Base_Alt = 256,       Sds_Num_Level2_Base_Alt_Flag = 0, &
     Sds_Num_Level2_Ch31_Max_3x3 = 257,   Sds_Num_Level2_Ch31_Max_3x3_Flag = 0, &
     Sds_Num_Level2_Ch20_Emiss = 258,     Sds_Num_Level2_Ch20_Emiss_Flag = 0, &
     Sds_Num_Level2_Ch20_Emiss_Clear=259, Sds_Num_Level2_Ch20_Emiss_Clear_Flag = 0, &
     Sds_Num_Level2_Ch20_Emiss_Median=260,Sds_Num_Level2_Ch20_Emiss_Median_Flag = 0, &
     Sds_Num_Level2_Ch1_Min_3x3 = 261,    Sds_Num_Level2_Ch1_Min_3x3_Flag = 0, &
     Sds_Num_Level2_Bt_11_67_Covar = 262, Sds_Num_Level2_Bt_11_67_Covar_Flag = 0, &
     Sds_Num_Level2_Btd_Ch31_Ch32 = 263,  Sds_Num_Level2_Btd_Ch31_Ch32_Flag = 0, &
     Sds_Num_Level2_Ctt_Opa = 264,        Sds_Num_Level2_Ctt_Opa_Flag = 0, &
     Sds_Num_Level2_Ch20_Bt_Median = 265, Sds_Num_Level2_Ch20_Bt_Median_Flag = 0, &
     Sds_Num_Level2_Tsfc_Uni_Back = 266,  Sds_Num_Level2_Tsfc_Uni_Back_Flag = 0, &
     Sds_Num_Level2_Ceps_Aux = 267,       Sds_Num_Level2_Ceps_Aux_Flag = 0, &
     Sds_Num_Level2_Sndr_Fov = 268,       Sds_Num_Level2_Sndr_Fov_Flag = 0, &
     Sds_Num_Level2_Sndr_X = 269,         Sds_Num_Level2_Sndr_X_Flag = 0, &
     Sds_Num_Level2_Sndr_Y = 270,         Sds_Num_Level2_Sndr_Y_Flag = 0, &
     Sds_Num_Level2_Sndr_Mask = 271,      Sds_Num_Level2_Sndr_Mask_Flag = 0, &
     Sds_Num_Level2_Inver_Strength = 272, Sds_Num_Level2_Inver_Strength_Flag = 0, &
     Sds_Num_Level2_Inver_Top = 273,      Sds_Num_Level2_Inver_Top_Flag = 0, &
     Sds_Num_Level2_Inver_Base = 274,     Sds_Num_Level2_Inver_Base_Flag = 0, &
     Sds_Num_Level2_Ch31_Rad = 275 ,      Sds_Num_Level2_Ch31_Rad_Flag = 0, &
     Sds_Num_Level2_Cbh_Qf = 276,         Sds_Num_Level2_Cbh_Qf_Flag = 0, &
     Sds_Num_Level2_Cld_Phase_IR = 277,   Sds_Num_Level2_Cld_Phase_IR_Flag = 0, &
     Sds_Num_Level2_Cld_Type_IR = 278,    Sds_Num_Level2_Cld_Type_IR_Flag = 0, &
     Sds_Num_Level2_Tc_Lower = 279,       Sds_Num_Level2_Tc_Lower_Flag = 0, &
     Sds_Num_Level2_Tpw_Ac = 280,         Sds_Num_Level2_Tpw_Ac_Flag = 0, &
     Sds_Num_Level2_Ch1_Clr_BCM = 281,    Sds_Num_Level2_Ch1_Clr_BCM_Flag = 1, &
     Sds_Num_Level2_Ch1_Clr_Min_BCM = 282,  Sds_Num_Level2_Ch1_Clr_Min_BCM_Flag = 1, &
     Sds_Num_Level2_Ch1_Clr_Max_BCM = 283,  Sds_Num_Level2_Ch1_Clr_Max_BCM_Flag = 1, &
     Sds_Num_Level2_Ch1_Clr_Std_BCM = 284,  Sds_Num_Level2_Ch1_Clr_Std_BCM_Flag = 1, &
     Sds_Num_Level2_BTD_11_12_NWC = 285,    Sds_Num_Level2_BTD_11_12_NWC_Flag = 1, &
     Sds_Num_Level2_EMISS_39_NWC = 286,     Sds_Num_Level2_EMISS_39_NWC_Flag = 1, &
     Sds_Num_Level2_EMISS_39_Clr_BCM = 287, Sds_Num_Level2_EMISS_39_Clr_BCM_Flag = 1, &
     Sds_Num_Level2_etrop_11um_LRC = 288,   Sds_Num_Level2_etrop_11um_LRC_Flag = 1, &
     Sds_Num_Level2_sfc_Ch1_Clr_fill = 289,   Sds_Num_Level2_sfc_Ch1_Clr_fill_flag = 1
