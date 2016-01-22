!----------------------------------------------------------------------
!--- viewing geometry thresholds
!----------------------------------------------------------------------
 REAL (KIND=REAL4), PRIVATE, PARAMETER:: Day_Sol_Zen_Thresh = 87.0
 REAL (KIND=REAL4), PRIVATE, PARAMETER:: Refl_Test_Sol_Zen_Thresh = 75.0
 REAL (KIND=REAL4), PRIVATE, PARAMETER:: Refl_Test_Air_Mass_Thresh = 6.0
 REAL (KIND=REAL4), PRIVATE, PARAMETER:: RVCT_Sol_Zen_Thresh = 70.0
 REAL (KIND=REAL4), PRIVATE, PARAMETER:: RVCT_Scat_Zen_Thresh = 90.0
 REAL (KIND=REAL4), PRIVATE, PARAMETER:: Glint_Zen_Thresh = 40.0


 REAL (KIND=REAL4), PRIVATE, PARAMETER:: Water_Freezing_Point = 273.15
 REAL (KIND=REAL4), PRIVATE, PARAMETER:: Max_Glint_Clr_Obs_Chn14_BTDIFF = 5.0
 REAL (KIND=REAL4), PRIVATE, PARAMETER:: Max_Glint_Clr_Rel_Refl2_Stddev_Thresh = 0.10
 REAL (KIND=REAL4), PARAMETER, PRIVATE:: Rayleigh_Optical_Depth_Chn2 = 0.05
 REAL (KIND=REAL4), PARAMETER, PRIVATE:: Rayleigh_Single_Scatter_Albedo = 1.00

!--- Aerosol
 Real(kind=real4), parameter, private:: Aerosol_Single_Scatter_Albedo_Chn2 = 1.0
 Real(kind=real4), parameter, private:: Aerosol_Optical_Depth_Chn2 = 0.1
 Real(kind=real4), parameter, private:: Aerosol_Asymmetry_Parameter = 0.8

!----------------------------------------------------------------------
! cloud detection thresholds
!----------------------------------------------------------------------
 REAL(KIND=REAL4), PRIVATE, PARAMETER:: BT_Chn14_Max_FMFT_Thresh = 310.0, &
                                        BTDIFF_Chn14_Chn15_Min_FMFT_Thresh = 0.0, &
                                        PFMFT_Land_Thresh = 2.50, &
                                        PFMFT_Ocean_Thresh = 0.80, &
                                        PFMFT_Snow_Thresh = 1.0, &
                                        PFMFT_Cold_Surface_Thresh = 1.0, &
                                        NFMFT_Land_Thresh = 0.60, &
                                        NFMFT_Ocean_Thresh = 0.60, &
                                        NFMFT_Snow_Thresh = 2.00, &
					NFMFT_BTD_Chn14_Chn15_Max_Thresh = 1.50, &
                                        Beta_Tropo_Chn14_Chn15_Min_SWCIRRUS_Thresh = 0.95, &
                                        Beta_Tropo_Chn14_Chn15_Max_SWCIRRUS_Thresh = 1.30, &
                                        BT_Chn14_Max_SWCIRRUS_Thresh = 310.0, &
                                        Emiss_Tropo_Chn14_Max_SWCIRRUS_Thresh = 0.95, &
                                        Emiss_Tropo_Chn14_Min_SWCIRRUS_Thresh = 0.05, &
                                        BT_Chn14_Max_TMFT_Thresh = 310.0, &
                                        BTDIFF_Chn7_Chn14_Min_TMFT_Thresh = -10.0, &
                                        TMFT_Land_Thresh =  5.0, & !10.0, &
                                        TMFT_Ocean_Thresh = 3.0, & !5.0, &
                                        BT_Chn14_Clr_Uni_Thresh_Land = 1.1, &
                                        BT_Chn14_Clr_Uni_Thresh_Ocn = 0.60, &
                                        Refl_Chn2_Clr_Uni_Thresh_Land = 0.20, &
                                        Refl_Chn2_Clr_Uni_Thresh_Ocean = 1.0, &
                                        Sfc_Temp_Uni_NWP_Thresh = 4.0, &
                                        Emiss_Chn14_Tropo_Ocn_Thresh = 0.05,  &   
                                        Emiss_Chn14_Tropo_Land_Thresh = 0.30, &     
                                        Emiss_Chn14_Tropo_Snow_Thresh = 0.30, &    
                                        Emiss_Chn14_Tropo_Cold_Surface_Thresh = 0.50, &    
                                        Emiss_Chn14_Tropo_Desert_Thresh = 0.40, &   
                                        Emiss_Chn14_Tropo_LRC_Ocn_Thresh = 0.28,  &   
                                        Emiss_Chn14_Tropo_LRC_Land_Thresh = 0.30, &    
                                        Emiss_Chn14_Tropo_LRC_Desert_Thresh = 0.40, &    
                                        Emiss_Chn14_Tropo_LRC_Snow_Thresh = 0.50, &    
                                        Emiss_Chn14_Tropo_LRC_Cold_Surface_Thresh = 0.50, &    
                                        EMISS4_Emiss_Chn7_High_Ocn_Thresh = 0.22, &
                                        EMISS4_Emiss_Chn7_High_Land_Thresh = 0.40, & !0.20,  &
                                        EMISS4_Emiss_Chn7_High_Snow_Thresh = 0.40, & !0.10,  &
                                        EMISS4_BT_Chn14_Max_Thresh = 310.0,  &
                                        NIR_Chn5_Refl_Thresh_Snow = 15.0, &         !NIRREF threshold on band 5 
                                        NIR_Chn7_Refl_Thresh_Snow = 15.0, &         !NIRREF threshold on band 7 
                                        CIRREF_Thresh = 5.0, &                      !CIRREF threshold on band 4
                                        CIRREF_Sfc_Hgt_Limit = 2000.0, &            !maximum sfc alt (m) to apply CIRREF
                                        NIRREF_Sfc_Hgt_Limit = 1000.0, &            !maximum sfc alt (m) to apply NIRREF
                                        ULST_Emiss_Chn7_Clr_Min_Thresh = 0.85, &
                                        ULST_Emiss_Chn7_Clr_Max_Thresh = 1.25, &
                                        ULST_Emiss_Chn7_Max_Thresh = 1.05, &
                                       !ULST_Emiss_Chn7_Diff_Thresh_Land_GOES = 0.10, &
                                       !ULST_Emiss_Chn7_Diff_Thresh_Ocean_GOES = 0.05, &
                                       !ULST_Emiss_Chn7_Diff_Thresh_Land_SEVIRI = 0.20, &
                                       !ULST_Emiss_Chn7_Diff_Thresh_Ocean_SEVIRI = 0.20, &
                                        ULST_Emiss_Chn7_Diff_Thresh_Land = 0.12, &
                                        ULST_Emiss_Chn7_Diff_Thresh_Ocean = 0.08, &
                                        ULST_Emiss_Chn7_Diff_Thresh_Snow = 0.12, &
                                        ULST_Emiss_Tropo_Chn14_Min_Thresh = -1.00,        &
                                        ULST_BT_Chn14_Max_Thresh = 310.00,        &
                                        RFMFT_Hi_Land_Thresh = 0.80, &
                                        RFMFT_Hi_Ocean_Thresh = 0.30, &
                                        RFMFT_BT_Chn14_Max_Thresh=300.0, &
                                        RFMFT_BTDIFF_Chn14_Chn15_Min_Thresh=1.0, &
                                        MECIDA_Chn10_Low_Pass_Thresh=0.5, &
                                        MECIDA_Chn16_Thresh1=253.0, &
                                        MECIDA_Chn16_Thresh2=233.0

 REAL(KIND=REAL4), PRIVATE, PARAMETER:: RTCT_Ocn_Thresh = 3.2   
 REAL(KIND=REAL4), PRIVATE, PARAMETER:: RTCT_Land_Thresh = 4.1

!----------------------------------------------------------------------
! other thresholds
!----------------------------------------------------------------------
 REAL (KIND=REAL4), PRIVATE, PARAMETER:: Refl_Chn2_Snow_Thresh = 10.0    !minimum albedo to allow snow to exist
 REAL (KIND=REAL4), PRIVATE, PARAMETER:: BT_Chn14_Snow_Thresh = 277.0   !maximum temperature to allow snow to exist
 REAL (KIND=REAL4), PRIVATE, PARAMETER:: BT_Chn14_15Min_Temporal_Offset = 2.0  !maximum temperature to allow snow to exist
