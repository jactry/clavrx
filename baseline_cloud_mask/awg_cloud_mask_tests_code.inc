!$Id: awg_cloud_mask_tests_code.inc,v 1.29 2011/04/15 16:38:46 wstraka Exp $
!====================================================================
! Subroutine Name: Set_Cmask_Thresholds
!
! Function:
!   Set Thresholds for tests that involve the 3.9 micron channel
!
! Description:
!   This Subroutine sets the thresholds for the Emiss_4 and ULST for a 
!   given satellite. This is needed because the 3.9 micron channel have
!   dIFferent spectral responce functions for various instruments.
!
! Calling Sequence:
!   CALL Set_Cmask_Thresholds(satellite_name)
!
! Inputs: Sat_name
!
! Outputs: None (thresholds in include are set)
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. In
!        this case, the thresholds for each satellite.
!
! Restrictions:  None
!
!====================================================================

 SUBROUTINE Set_Cmask_Thresholds(Sat_name, Algo_Name)
    CHARACTER(*), INTENT(IN) :: Sat_name
    CHARACTER(*), INTENT(inout) :: Algo_Name
    CHARACTER(len=256) :: Algo_Name_Tmpy

    Algo_Name_Tmpy = TRIM(trim(Algo_Name)//"abi")
    
   ! Set thresholds for AVHRR  
   IF ((index(TRIM(Sat_name),'AVHRR')) /= 0) THEN
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_AVHRR
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_AVHRR
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_AVHRR
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_AVHRR
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_AVHRR
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_AVHRR
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_AVHRR

   ! Set thresholds for EOS TERRA MODIS
   ELSE IF ((index(TRIM(Sat_name),'Terra')) /= 0) THEN
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_MODIS
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_MODIS
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_MODIS
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_MODIS
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_MODIS
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_MODIS
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_MODIS

   ! Set thresholds for EOS AQUA MODIS
   ELSE IF ((index(TRIM(Sat_name),'Aqua')) /= 0) THEN
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_MODIS
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_MODIS
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_MODIS
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_MODIS
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_MODIS
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_MODIS
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_MODIS
         
   ! Set thresholds for MSG 
   ELSE IF ((index(TRIM(Sat_name),'Meteosat')) /= 0) THEN
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_SEVIRI
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_SEVIRI
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_SEVIRI
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_SEVIRI
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_SEVIRI
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_SEVIRI
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_SEVIRI
      Algo_Name_Tmpy = TRIM(trim(Algo_Name)//"seviri")

   ! Set thresholds for MTSAT
   ELSE IF ((index(TRIM(Sat_name),'MTSAT')) /= 0) THEN
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_MTSAT
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_MTSAT
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_MTSAT
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_MTSAT
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_MTSAT
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_MTSAT
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_MTSAT

   ! Set thresholds for GOES (I-P, ABI)
   ELSE
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_GOES
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_GOES
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_GOES
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_GOES
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_GOES
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_GOES
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_GOES
   ENDIF 
   
   
    Algo_Name = trim(Algo_Name_Tmpy)
   
 END SUBROUTINE Set_Cmask_Thresholds


!====================================================================
! Function Name: Atmospheric_Scattering_Chn2
!
! Function:
!   Computes the Single Scatteing Reflectance in Chn2 due to assumed Rayleigh 
!   and Aerosol.  Background values are used.  Should be replaced by CRTM.
!
! Description:
!   Based on the cosine of the satellite zenith angle, cosine of the solar
!   zenith angle, the cosine of the scattering angle and the air mass factor,
!   the single scattering reflectance for the 0.64 micron channel is calcuated
!   and returned.
!   
!
! Calling Sequence:
!   Refl_Sing_Scat = Atmospheric_Scattering_Chn2(Cos_Sat_Zen,  &
!                                                Cos_Sol_Zen, &
!                                                Cos_Scat_Zen, &
!                                                Air_Mass_Factor)
!
! Inputs:
!   cosine of the satellite zenith angle
!   cosine of the solar zenith angle
!   cosine of the scattering angle
!   air mass factor
!
! Outputs: Function returns the atmosphericly corrected 0.64 micron reflectance (Refl_Atm_Chn2)
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

FUNCTION Atmospheric_Scattering_Chn2(Cos_Sat_Zen, &
                                     Cos_Sol_Zen, &
                                     Cos_Scat_Zen, &
                                     Aerosol_Optical_Depth_Chn2, &
                                     Air_Mass_Factor) &
                                     RESULT(Refl_Atm_Chn2)

  REAL(KIND=REAL4), INTENT(IN):: Cos_Sat_Zen
  REAL(KIND=REAL4), INTENT(IN):: Cos_Sol_Zen
  REAL(KIND=REAL4), INTENT(IN):: Cos_Scat_Zen
  REAL(KIND=REAL4), INTENT(IN):: Air_Mass_Factor
  REAL(KIND=REAL4), INTENT(IN):: Aerosol_Optical_Depth_Chn2
  REAL(KIND=REAL4):: Refl_Atm_Chn2
  REAL(KIND=REAL4):: Transmission
  REAL(KIND=REAL4):: Phase_Function_Ray
  REAL(KIND=REAL4):: Phase_Function_Aer
  REAL(KIND=REAL4):: Total_Optical_Depth
  REAL(KIND=REAL4):: Single_Scatter_Albedo
  REAL(KIND=REAL4):: Phase_Function
                              
         !--- initialize
         Refl_Atm_Chn2 = 0.0

         !--- check for day conditions 
         IF (Cos_Sol_Zen > 0.0) THEN

            !--- compute phase functions for rayleigh and aerosol
            Phase_Function_Ray = 0.75*(1.0 + Cos_Scat_Zen**2)

            Phase_Function_Aer = (1.0 - Aerosol_Asymmetry_Parameter**2) / &
                                 ( (1.0 + Aerosol_Asymmetry_Parameter**2 -  &
                                   2.0*Aerosol_Asymmetry_Parameter*Cos_Scat_Zen)**(1.5) )

            Total_Optical_Depth = Rayleigh_Optical_Depth_Chn2 + Aerosol_Optical_Depth_Chn2

            !--- single scatter albedo
            Single_Scatter_Albedo =  &
                  (Aerosol_Single_Scatter_Albedo_Chn2*Aerosol_Optical_Depth_Chn2 +  &
                   Rayleigh_Optical_Depth_Chn2)/ ( Total_Optical_Depth )

            !--- compute effective phase function
            Phase_Function = (Aerosol_Single_Scatter_Albedo_Chn2*Aerosol_Optical_Depth_Chn2*Phase_Function_Aer +  &
                              Rayleigh_Optical_Depth_Chn2*Phase_Function_Ray)/ &
                             (Single_Scatter_Albedo * Total_Optical_Depth)


            Transmission = exp(-Air_Mass_Factor * Total_Optical_Depth)

            Refl_Atm_Chn2 = 100.0*Single_Scatter_Albedo * Phase_Function / &
                             (4.0*Air_Mass_Factor*Cos_Sat_Zen*Cos_Sol_Zen) * (1.0 - Transmission )

         ENDIF

         RETURN

END FUNCTION Atmospheric_Scattering_Chn2


!====================================================================
! Function Name: Compute_Clear_Sky_Scatter
!
! Function:
!   Computes the single scater and aerosol reflectance assuming that the
!   gas is mixed in with scattering
!
! Description:
   
!
! Calling Sequence:
!   Refl_Sing_Scat = Compute_Clear_Sky_Scatter(Aerosol_Optical_Depth_Chn2, &
!                                     Aerosol_Single_Scatter_Albedo_Chn2, &
!                                     Aerosol_Asymmetry_Parameter, &
!                                     Rayleigh_Optical_Depth_Chn2, &
!                                     Sat_Name, & 
!                                     TPW, &
!                                     Total_Ozone_Path_NWP, &
!                                     Scat_Zen, &
!                                     Cos_Sat_Zen, &
!                                     Cos_Sol_Zen, &
!                                     Sfc_Alb_View, &
!                                     Sfc_Alb_Sun, &
!                                     Transmission_Sing_Scat, &
!                                     Refl_Sing_Scat)
!
! Inputs:
!   Aerosol Optical depth at 0.64 micron
!   Aerosol Single Scatter Albedeo at 0.64 micron
!   Aerosol Asymmetry Parameter
!   Rayleigh Optical Depth at 0.64 micron
!   Satellite name
!   Total Precipitable Water
!   Total_Ozone_Path
!   scattering angle
!   cosine of the satellite zenith angle
!   cosine of the solar zenith angle
!   Clear sky Surface Albedo 
!   Clear Sky Surface Albedo
!   Single Scattering Transmission Coefficent
!
! Outputs: Single Scatter reflectance
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================


 SUBROUTINE Compute_Clear_Sky_Scatter(Aerosol_Optical_Depth_Chn2, &
                                     Aerosol_Single_Scatter_Albedo_Chn2, &
                                     Aerosol_Asymmetry_Parameter, &
                                     Rayleigh_Optical_Depth_Chn2, &
                                     Sat_Name, & 
                                     TPW, &
                                     Total_Ozone_Path_NWP, &
                                     Scat_Zen, &
                                     Cos_Sat_Zen, &
                                     Cos_Sol_Zen, &
                                     Sfc_Alb_View, &
                                     Sfc_Alb_Sun, &
                                     Transmission_Sing_Scat, &
                                     Refl_Sing_Scat)

   REAL, INTENT(IN):: Aerosol_Optical_Depth_Chn2
   REAL, INTENT(IN):: Aerosol_Single_Scatter_Albedo_Chn2
   REAL, INTENT(IN):: Aerosol_Asymmetry_Parameter
   REAL, INTENT(IN):: Rayleigh_Optical_Depth_Chn2
   CHARACTER (len=*), INTENT(IN):: Sat_Name 
   REAL, INTENT(IN):: TPW
   REAL, INTENT(IN):: Total_Ozone_Path_NWP
   REAL, INTENT(IN):: Scat_Zen
   REAL, INTENT(IN):: Cos_Sat_Zen
   REAL, INTENT(IN):: Cos_Sol_Zen
   REAL, INTENT(IN):: Sfc_Alb_View
   REAL, INTENT(IN):: Sfc_Alb_Sun
   REAL, INTENT(out):: Transmission_Sing_Scat
   REAL, INTENT(out):: Refl_Sing_Scat
   REAL:: Air_Mass_Factor
   REAL:: Aero_Phase_Funct
   REAL:: Ray_Phase_Funct
   REAL:: OD_Gas
   REAL:: OD_Total
   REAL:: OD_Scat_Total
   REAL:: OD_Iso_Total
   REAL:: Trans_Iso_Total_View
   REAL:: Trans_Iso_Total_Sun
   REAL:: OD_Iso_Scat_Total
   REAL:: Cos_Scat_Zen
   REAL:: Eff_Phase_Funct
   REAL:: Single_Scat_Alb
   REAL:: Refl_Sing_Scat_a
   REAL:: Refl_Sing_Scat_b
   REAL:: Refl_Sing_Scat_c
   REAL:: OD_ozone
   REAL:: OD_h2o
   REAL, DIMENSION(3):: OD_ozone_coef
   REAL, DIMENSION(3):: OD_h2o_coef

   !--- Set Default Gas Transmission Coefficients
   OD_ozone_coef = (/0.000566454,8.25224e-05,1.94007e-08/)
   OD_h2o_coef = (/  0.000044758, 0.00264790,-0.0000713698/)


   !--- Set Satellite SpecIFic Gas Transmission Coefficients
   IF ((index(TRIM(Sat_Name),'Meteosat')) /= 0) THEN
       OD_ozone_coef = (/0.000566454,8.25224e-05,1.94007e-08/)
       OD_h2o_coef = (/  0.000044758, 0.00264790,-0.0000713698/)
   ENDIF

   !--- compute cosine of scattering angle
   Cos_Scat_Zen = cos(Scat_Zen * dtor)

   !--- compute gaseous optical depth and transmission
   OD_h2o = OD_h2o_coef(1) + OD_h2o_coef(2)*TPW + OD_h2o_coef(3) * TPW**2
   
   OD_ozone = OD_ozone_coef(1) + OD_ozone_coef(2) * Total_Ozone_Path_NWP + &
              OD_ozone_coef(3) * Total_Ozone_Path_NWP**2

   OD_Gas = OD_h2o + OD_ozone

   !-- compute Rayleigh phase function
   Air_Mass_Factor = 1.0 / Cos_Sat_Zen + 1.0 / Cos_Sol_Zen
   
   Ray_Phase_Funct = 0.75 * (1.0 + Cos_Scat_Zen**2)

   !--- compute total transmission
   OD_Total = Aerosol_Optical_Depth_Chn2 + Rayleigh_Optical_Depth_Chn2 + OD_Gas
   
   Transmission_Sing_Scat = exp(-OD_Total * Air_Mass_Factor)

   OD_Iso_Total = (1.0 - Aerosol_Asymmetry_Parameter) * Aerosol_Optical_Depth_Chn2 + &
                   Rayleigh_Optical_Depth_Chn2 + OD_Gas
   
   Trans_Iso_Total_View = exp(-OD_Iso_Total / Cos_Sat_Zen)
   
   Trans_Iso_Total_Sun = exp(-OD_Iso_Total / Cos_Sol_Zen)

   !--- compute total scattering optical depth
   OD_Scat_Total = Aerosol_Single_Scatter_Albedo_Chn2 *&
                   Aerosol_Optical_Depth_Chn2 + Rayleigh_Optical_Depth_Chn2
   
   OD_Iso_Scat_Total = Aerosol_Single_Scatter_Albedo_Chn2 * (1.0 - Aerosol_Asymmetry_Parameter) * &
                    Aerosol_Optical_Depth_Chn2 + Rayleigh_Optical_Depth_Chn2

   !--- single scatter albedo
   Single_Scat_Alb = (Aerosol_Single_Scatter_Albedo_Chn2 * Aerosol_Optical_Depth_Chn2 + &
                     Rayleigh_Optical_Depth_Chn2) / ( OD_Total )

   !aerosol phase function (Henyey-Greenstein)
   Aero_Phase_Funct = (1.0 - Aerosol_Asymmetry_Parameter**2) / &
                    ( (1.0 + Aerosol_Asymmetry_Parameter**2 -  &
                       2.0 * Aerosol_Asymmetry_Parameter*Cos_Scat_Zen)**(1.5) )

   !--- compute effective phase function
   Eff_Phase_Funct = (Aerosol_Single_Scatter_Albedo_Chn2 * Aerosol_Optical_Depth_Chn2 * &
   Aero_Phase_Funct + Rayleigh_Optical_Depth_Chn2 * Ray_Phase_Funct) / (OD_Scat_Total)

   !--- compute single scatter reflectance (0-100%)
   Refl_Sing_Scat_a = Single_Scat_Alb * Eff_Phase_Funct / (4.0 * Air_Mass_Factor * &
                        Cos_Sat_Zen * Cos_Sol_Zen) * (1.0 - &
                        Transmission_Sing_Scat )

   Refl_Sing_Scat_b = (OD_Iso_Scat_Total / (2.0*Cos_Sol_Zen)) * &
                        Trans_Iso_Total_View * Sfc_Alb_View

   Refl_Sing_Scat_c = (OD_Iso_Scat_Total / (2.0 * Cos_Sat_Zen)) * &
                       Trans_Iso_Total_Sun * Sfc_Alb_Sun

   Refl_Sing_Scat = 100.0 * (Refl_Sing_Scat_a + Refl_Sing_Scat_b + Refl_Sing_Scat_c)

 END SUBROUTINE Compute_Clear_Sky_Scatter

!====================================================================
! Function Name: RUT_Routine
!
! Function:
!   Reflectance Uniformity Routine (RUT) Test
!
! Description:
!   Computes the Reflectance Uniformity Routine (RUT) test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RUT_Routine (&
!                              Is_Land, &
!                              Is_Snow, &
!                              Is_Coast, &
!                              Refl_Chn2_Clear(Elem_Idx,Line_Idx), &
!                              Refl_Chn2_Stddev_3x3(Elem_Idx,Line_Idx), &
!                              Sol_Zen)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel coast (YES/NO)
!   0.64 miron clear sky reflectance
!   Standard Deviation of the 0.64 micron reflectance over a 3x3 box
!   Solar zenith angle
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test (Test_Result)
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION RUT_Routine(Is_Land, &
                     Is_Snow, &
                     Is_Coast, &
                     Refl_Chn2_Clear, &
                     Refl_Chn2_Stddev_3x3, &
                     Sol_Zen) &
                     RESULT(Test_Result)

       INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
       INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
       INTEGER(KIND=INT4), INTENT(IN) :: Is_Coast
       REAL(KIND=REAL4), INTENT(IN) :: Refl_Chn2_Clear
       REAL(KIND=REAL4), INTENT(IN) :: Refl_Chn2_Stddev_3x3
       REAL(KIND=REAL4), INTENT(IN) :: Sol_Zen
       INTEGER(KIND=int1) :: Test_Result

       REAL(KIND=REAL4) :: Test_Threshold

       Test_Result = sym%NO

       IF (Is_Snow == sym%NO .AND. Is_Coast == sym%NO) THEN

          IF (Sol_Zen < RUT_SOL_ZEN_THRESH) THEN

              !--- compute threshold
              IF (Is_Land == sym%YES) THEN
                 Test_Threshold = MAX(0.5,Refl_Chn2_Clear * REFL_CHN2_CLR_UNI_THRESH_LAND)
              ELSE
                 Test_Threshold = REFL_CHN2_CLR_UNI_THRESH_OCEAN
              ENDIF

              !--- apply test
              IF (Refl_Chn2_Stddev_3x3 > Test_Threshold) THEN
                     Test_Result = sym%YES
              ENDIF

           ENDIF

       ENDIF

       RETURN

 END FUNCTION RUT_Routine

!====================================================================
! Function Name: TUT_Routine
!
! Function:
!   Thermal Uniformity Routine (TUT) Test
!
! Description:
!   Computes the Thermal Uniformity Routine (TUT) test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = TUT_Routine (&
!                              Is_Land, &
!                              Is_Coast, &
!                              BT_Chn14_Stddev_3x3, &
!                              Sfc_Hgt_Stddev_3x3)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel coast (YES/NO)
!   Standard Deviation of the 11 micron BT over a 3x3 box
!   Standard Deviation of the surface height over a 3x3 box
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION TUT_Routine(Is_Land, &
                      Is_Coast, &
                      BT_Chn14_Stddev_3x3, &
                      Sfc_Hgt_Stddev_3x3)  &
                      RESULT (Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land 
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Coast
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Stddev_3x3
   REAL(KIND=REAL4), INTENT(IN) :: Sfc_Hgt_Stddev_3x3
   INTEGER(KIND=int1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold

   Test_Result = sym%NO

   IF (Is_Coast == sym%NO) THEN

         !
         !7K/km is the adiabatic lapse rate
         !
         Test_Threshold = 3.0 * 7.0*Sfc_Hgt_Stddev_3X3/1000.0  

         IF (Is_Land == sym%YES) THEN

      IF (BT_Chn14_Stddev_3x3 > BT_CHN14_CLR_UNI_THRESH_LAND  &
          +Test_Threshold) THEN

         Test_Result = sym%YES

      ENDIF

         ELSE

      IF (BT_Chn14_Stddev_3x3 > BT_CHN14_CLR_UNI_THRESH_OCN  &
           +Test_Threshold) THEN

         Test_Result = sym%YES

      ENDIF

         ENDIF

   ENDIF

   RETURN

 END FUNCTION TUT_Routine

!====================================================================
! Function Name: RTCT_Routine
!
! Function:
!   Relative Thermal Contrast Test
!
! Description:
!   Computes the Relative Thermal Contrast Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RTCT_Routine(&
!                              Is_Land,  &
!                              Is_Coast,  &
!                              Is_Snow, &
!                              Is_Cold_Surface, &
!                              BT_Chn14, &
!                              BT_Chn14_Min_3x3(Elem_Idx,Line_Idx), &
!                              BT_Chn14_Max_3x3(Elem_Idx,Line_Idx), &
!                              Sfc_Hgt_Stddev_3x3(Elem_Idx,Line_Idx))
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel coast (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel cold surface (YES/NO)
!   11 micron BT 
!   Minimum of the 11 micron BT in a 3x3 box
!   Maximum of the 11 micron BT in a 3x3 box
!   Standard Deviation of the surface height over a 3x3 box
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION RTCT_Routine(Is_Land, &
                      Is_Coast,  &
                      Is_Snow, &
                      Is_Cold_Surface, &
                      BT_Chn14, &
                      BT_Chn14_Min_3x3, &
                      BT_Chn14_Max_3x3, &
                      Sfc_Hgt_Stddev_3x3) &
                      RESULT (Test_Result)

   INTEGER (KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER (KIND=INT4), INTENT(IN) :: Is_Coast
   INTEGER (KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER (KIND=INT4), INTENT(IN) :: Is_Cold_Surface
   REAL (KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL (KIND=REAL4), INTENT(IN) :: BT_Chn14_Min_3x3
   REAL (KIND=REAL4), INTENT(IN) :: BT_Chn14_Max_3x3
   REAL (KIND=REAL4), INTENT(IN) :: Sfc_Hgt_Stddev_3x3
   INTEGER (KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold

   Test_Result = sym%NO

   IF (Is_Cold_Surface == sym%NO) THEN

       IF (Is_Land == sym%NO) THEN

            Test_Threshold = RTCT_OCN_THRESH

       ELSE

            Test_Threshold = RTCT_LAND_THRESH

       ENDIF

       !
       !7K/km is the adiabatic lapse rate
       !
       Test_Threshold = Test_Threshold + 3.0 * 7.0*Sfc_Hgt_Stddev_3X3/1000.0

       IF ((Is_Coast == sym%NO) .AND.  &
           (Is_Snow == sym%NO) .AND.  &
           (BT_Chn14_Min_3x3 <= 300.0) .AND.  &
           (BT_Chn14_Max_3x3 - BT_Chn14 > Test_Threshold)) THEN

           Test_Result = sym%YES

       ENDIF

   ENDIF

 END FUNCTION RTCT_Routine


!====================================================================
! Function Name: ETROP_Routine
!
! Function:
!   11 micron Tropospheric Emissivity Test
!
! Description:
!   Computes the 11 micron Tropospheric Emissivity Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  ETROP_Routine( &
!                               Is_Snow, &
!                               Is_Land, &
!                               Is_Coast, &
!                               Is_Desert, &
!                               Is_Cold_Surface, &
!                               Land_Type, &
!                               BT_Chn14, &
!                               BT_Chn14_Clr, &
!                               Emiss_Tropo_Chn14(Elem_Idx,Line_Idx), &
!                               Emiss_Tropo_Chn14_LRC, &
!                               BT_Chn14_Stddev_3x3(Elem_Idx,Line_Idx))
!
! Inputs:
!   Is pixel snow (YES/NO)
!   Is pixel land (YES/NO)
!   Is pixel coast (YES/NO)
!   Is pixel desert (YES/NO)
!   Is pixel cold surface (YES/NO)
!   Is pixel land type
!   11 micron BT 
!   Clear sky 11 micron BT 
!   11 micron tropopause emissivity
!   11 micron tropopause emissivity at the Local radiative center
!   Standard Deviation of the 11 micron BT over a 3x3 box
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION ETROP_Routine(Is_Snow, &
                        Is_Land, &
                        Is_Coast, &
                        Is_Desert, &
                        Is_Cold_Surface, &
                        Land_Type, &
                        BT_Chn14, &
                        BT_Chn14_Clr, &
                        Emiss_Tropo_Chn14, &
                        Emiss_Tropo_Chn14_LRC, &
                        BT_Chn14_Stddev_3x3) &
                        RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Coast
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Desert
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Cold_Surface
   INTEGER(KIND=INT4), INTENT(IN) :: Land_Type
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Clr
   REAL(KIND=REAL4), INTENT(IN) :: Emiss_Tropo_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: Emiss_Tropo_Chn14_LRC
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Stddev_3x3
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Value
   REAL(KIND=REAL4) :: Test_Threshold
   REAL(KIND=REAL4) :: Test_LRC_Threshold

   Test_Result = sym%NO

   IF ((BT_Chn14 > 170.0) .AND.  &
       (BT_Chn14 < 315.0) .AND.  &
       (BT_Chn14_Clr > 240.0)) THEN

       !--- assign threshold
       Test_Threshold = EMISS_CHN14_TROPO_LAND_THRESH
       Test_LRC_Threshold = EMISS_CHN14_TROPO_LRC_LAND_THRESH
       
       IF (Is_Land == sym%NO .AND. Is_Coast == sym%NO) THEN 
         Test_Threshold = EMISS_CHN14_TROPO_OCN_THRESH
         Test_LRC_Threshold = EMISS_CHN14_TROPO_LRC_OCN_THRESH
       ENDIF
       IF (Is_Snow == sym%YES) THEN 
         Test_Threshold = EMISS_CHN14_TROPO_SNOW_THRESH
         Test_LRC_Threshold = EMISS_CHN14_TROPO_LRC_SNOW_THRESH
       ENDIF
       IF (Is_Desert == sym%YES) THEN 
         Test_Threshold = EMISS_CHN14_TROPO_DESERT_THRESH
         Test_LRC_Threshold = EMISS_CHN14_TROPO_LRC_DESERT_THRESH
       ENDIF
       IF (Is_Cold_Surface == sym%YES) THEN 
         Test_Threshold = EMISS_CHN14_TROPO_COLD_SURFACE_THRESH
         Test_LRC_Threshold = EMISS_CHN14_TROPO_LRC_COLD_SURFACE_THRESH
       ENDIF

       !--- select value to test
       Test_Value = Emiss_Tropo_Chn14

       !--- apply test
       IF (Test_Value > Test_Threshold) THEN
         Test_Result = sym%YES
       ENDIF

       !--- apply LRC portion of test
       IF (Emiss_Tropo_Chn14_LRC /= Missing_Value_REAL4) THEN

         !--- select value to test
         Test_Value = Emiss_Tropo_Chn14_LRC
        
         !--- apply test
         IF (Test_Value > Test_LRC_Threshold) THEN
          Test_Result = sym%YES
         ENDIF

       ENDIF



      !--- apply a tigher threshold for highly nonunIForm pixels
      Test_Threshold = 0.20
      IF ((Is_Coast == sym%NO) .AND. (BT_Chn14_Stddev_3x3 > 0.5)) THEN
         IF (Test_Value > Test_Threshold) THEN
             Test_Result = sym%YES
         ENDIF
      ENDIF 


      !--- perform a restoral near land where sst field is often errorneous
      !--- select value to test
      Test_Value = Emiss_Tropo_Chn14

      IF (Test_Result == sym%YES) THEN
         IF ((Test_Value < 0.20) .AND.  &
            (BT_Chn14_Stddev_3x3 < 1.0) .AND. &
            (Land_Type /= sym%LAND) .AND.  &
            (Land_Type /= sym%DEEP_OCEAN)) THEN

            Test_Result = sym%NO

         END IF 
      END IF 

   ENDIF

   RETURN

 END FUNCTION ETROP_Routine
 
 !====================================================================
! Function Name: PFMFT_Routine
!
! Function:
!   Positive Split Window DIFferences Test (PFMFT) 
!
! Description:
!   Computes the Positive Split Window DIFferences Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = PFMFT_Routine( &
!                      Is_Land, &
!                      Is_Snow, &
!                      Is_Cold_Surface, &
!                      BT_Chn14, &
!                      BT_Chn15, &
!                      BT_Chn14_Clr, &
!                      BT_Chn15_Clr, &
!                      BT_Chn14_Stddev_3x3)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel cold surface (YES/NO)
!   11 micron BT 
!   12 micron BT 
!   Clear sky 11 micron BT 
!   Clear sky 12 micron BT 
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION PFMFT_Routine( &
                      Is_Land, &
                      Is_Snow, &
                      Is_Cold_Surface, &
                      BT_Chn14, &
                      BT_Chn15, &
                      BT_Chn14_Clr, &
                      BT_Chn15_Clr, &
                      BT_Chn14_Stddev_3x3) &
                      RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Cold_Surface
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn15
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Clr
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn15_Clr
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Stddev_3x3
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold
   REAL(KIND=REAL4) :: Test_Value

   Test_Result = sym%NO

   IF((BT_Chn14 < BT_CHN14_MAX_FMFT_THRESH) .AND. &
      (BT_Chn14_Clr - BT_Chn15_Clr > BTDIFF_CHN14_CHN15_MIN_FMFT_THRESH) .AND. &
      (BT_Chn14_Stddev_3x3 > PFMFT_BT_CHN14_STDDEV_3x3_THRESH)) THEN

      IF ((BT_Chn14 > 270.0) .AND. (BT_Chn14_Clr > 270.0)) then
        Test_Value = (BT_Chn14 - BT_Chn15) -  &
                     (BT_Chn14_Clr - BT_Chn15_Clr) *(BT_Chn14 - 260.0) / &
                     (BT_Chn14_Clr - 260.0)
      ELSE
        Test_Value = (BT_Chn14 - BT_Chn15)
      ENDIF
 
 
      !-- set appropriate threshold
      IF (Is_Land == sym%NO) THEN
            Test_Threshold = PFMFT_OCEAN_THRESH
      ELSE
            Test_Threshold = PFMFT_LAND_THRESH
      ENDIF
      IF (Is_Snow == sym%YES) THEN
            Test_Threshold = PFMFT_SNOW_THRESH
      ENDIF
      IF (Is_Cold_Surface == sym%YES) THEN
            Test_Threshold = PFMFT_COLD_SURFACE_THRESH
      ENDIF

      !--- Apply Test
      IF (Test_Value > Test_Threshold) THEN
          Test_Result = sym%YES
      ENDIF
 
   ENDIF
 
   RETURN

 END FUNCTION PFMFT_Routine
 
!====================================================================
! Function Name: NFMFT_Routine
!
! Function:
!   Test for Negative Split Window DIFferences (NFMFT)
!
! Description:
!   Computes the Negative Split Window DIFferences Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = NFMFT_Routine( &
!                      Is_Land, &
!                      Is_Snow, &
!                      BT_Chn14, &
!                      BT_Chn15, &
!                      BT_Chn14_Clr, &
!                      BT_Chn15_Clr)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel desert (YES/NO)
!   11 micron BT 
!   12 micron BT 
!   Clear sky 11 micron BT 
!   Clear sky 12 micron BT 
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION NFMFT_Routine( &
                      Is_Land, &
                      Is_Snow, &
                      Is_Desert, &
                      BT_Chn14, &
                      BT_Chn15, &
                      BT_Chn14_Clr, &
                      BT_Chn15_Clr) &
                      RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Desert
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn15
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Clr
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn15_Clr
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold
   REAL(KIND=REAL4) :: Test_Value

   Test_Result = sym%NO

   !--- skip this test for elevated values of BT_Chn14-BT_Chn15
   IF (BT_Chn14 - BT_Chn15 > NFMFT_BTD_CHN14_CHN15_MAX_THRESH) THEN
      RETURN
   ENDIF

   !-- set appropriate threshold
   Test_Threshold = NFMFT_LAND_THRESH
   IF (Is_Land == sym%NO) THEN
         Test_Threshold = NFMFT_OCEAN_THRESH
   ENDIF
   IF (Is_Desert == sym%YES) THEN
        Test_Threshold = NFMFT_DESERT_THRESH
   ENDIF
   IF (Is_Snow == sym%YES) THEN
        Test_Threshold = NFMFT_SNOW_THRESH
   ENDIF

   !--- Compute Value to Test
   Test_Value = (BT_Chn14 - BT_Chn15) -  &
                (BT_Chn14_Clr - BT_Chn15_Clr)

   !--- Apply Test
   IF (Test_Value < Test_Threshold) THEN
       Test_Result = sym%YES
   ENDIF

 END FUNCTION NFMFT_Routine

!====================================================================
! Function Name: RFMFT_Routine
!
! Function:
!   Relative Split Window DIFferences Test (RFMFT) 
!
! Description:
!   Computes the Relative Split Window DIFferences Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RFMFT_Routine( &
!                      Is_Land, &
!                      Is_Coast, &
!                      BT_Chn14, &
!                      BT_Chn15, &
!                      BTD_Chn14_Chn15_NWC)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel coast (YES/NO)
!   11 micron BT 
!   12 micron BT 
!   11 - 12 micron BT DIFf at the nearest warm center
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION RFMFT_Routine( &
                      Is_Land, &
                      Is_Coast, &
                      BT_Chn14, &
                      BT_Chn15, &
                      BTD_Chn14_Chn15_NWC) &
                      RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Coast
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn15
   REAL(KIND=REAL4), INTENT(IN) :: BTD_Chn14_Chn15_NWC
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Value
   REAL(KIND=REAL4) :: Test_Threshold

   !--- initialize
   Test_Result = sym%NO

   !--- do not apply over hot land due to NWP surface temperature biases
   IF ((Is_Land == sym%YES) .AND. (BT_Chn14 > RFMFT_BT_CHN14_MAX_THRESH)) THEN
        Test_Result = sym%NO
        RETURN
   ENDIF

   IF (BT_Chn14 - BT_Chn15 < RFMFT_BTDIFF_CHN14_CHN15_MIN_THRESH) THEN
        Test_Result = sym%NO
        RETURN
   ENDIF

   IF (Is_Coast == sym%YES) THEN
        Test_Result = sym%NO
        RETURN
   ENDIF

   !--- test for high departures - cirrus

   !--- pick thresholds
   Test_Threshold = RFMFT_HI_LAND_THRESH
   IF (Is_Land == sym%NO) THEN
        Test_Threshold = RFMFT_HI_OCEAN_THRESH
   ENDIF

   !--- apply test
   Test_Value = abs((BT_Chn14 - BT_Chn15) - BTD_Chn14_Chn15_NWC)

   IF (Test_Value > Test_Threshold) THEN
        Test_Result = sym%YES
   ENDIF

 END FUNCTION RFMFT_Routine


!====================================================================
! Function Name: TEMPIR_Routine
!
! Function:
!   Temporal IR Test (TEMPIR) 
!
! Description:
!   Computes the Temporal IR Test (YES/NO), as described in the ABI cloud
!   mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = TEMPIR_Routine( &
!                          BT_Chn14, &
!                          BT_Chn14_Clr, &
!                          BT_Chn14_15min, &
!                          BT_Chn14_15min_Clr)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel coast (YES/NO)
!   11 micron BT 
!   Clear sky 11 micron BT 
!   11 micron BT from 15 minutes ago for same pixel
!   Clear sky 11 micron BT from 15 minutes ago for same pixel
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION  TEMPIR_Routine( &
                          BT_Chn14, &
                          BT_Chn14_Clr, &
                          BT_Chn14_15min, &
                          BT_Chn14_15min_Clr) &
                          RESULT(Test_Result)

   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Clr
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_15min
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_15min_Clr
   REAL(KIND=REAL4) :: Test_Threshold
   INTEGER(KIND=INT1) :: Test_Result

   !--- initialize output
   Test_Result = sym%NO

   IF ((BT_Chn14_15min /= Missing_Value_REAL4) .AND. &
       (BT_Chn14_15min_Clr /= Missing_Value_REAL4) .AND. &
       (BT_Chn14_15min < 330.0) .AND. &
       (BT_Chn14_15min_Clr < 330.0)) THEN

      !--- make threshold
      Test_Threshold = (BT_Chn14_15min_Clr - BT_Chn14_Clr) + &
                        TEMPIR_BT_CHN14_15MIN_TEMPORAL_OFFSET

      !--- apply test
      IF ((BT_Chn14_15min - BT_Chn14) > Test_Threshold)  THEN

        Test_Result = sym%YES
 
      ENDIF

   ENDIF

   RETURN

 END FUNCTION TEMPIR_Routine

!====================================================================
! Function Name: RGCT_Routine
!
! Function:
!   Reflectance Gross Constrast Test (RGCT)
!
! Description:
!   Computes the Reflectance Gross Constrast Test (RGCT) (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RGCT_Routine( &
!                       Is_Snow, &
!                       Is_Glint, &
!                       Sol_Zen, &
!                       Refl_Chn2, &
!                       RGCT_Threshold)
!
! Inputs:
!   Is pixel snow (YES/NO)
!   Is pixel glint (YES/NO)
!   Solar zenith angle
!   Air mass factor
!   0.64 miron reflectance
!   RGCT threshold for that pixel
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION RGCT_Routine( &
                       Is_Snow, &
                       Is_Glint, &
                       Sol_Zen, &
                       Refl_Chn2, &
                       RGCT_Threshold) &
                       RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Glint
   REAL(KIND=REAL4), INTENT(IN) :: Sol_Zen
   REAL(KIND=REAL4), INTENT(IN) :: Refl_Chn2
   REAL(KIND=REAL4), INTENT(IN) :: RGCT_Threshold
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold


   !--- initialize
   Test_Result = sym%NO


   !--- apply test
   IF (Is_Snow == sym%NO .AND. Is_Glint == sym%NO) THEN

        IF (Sol_Zen < RGCT_Sol_Zen_Thresh) THEN 

            Test_Threshold = RGCT_Threshold

            IF (Refl_Chn2 > Test_Threshold) THEN

                Test_Result =  sym%YES
  
            ENDIF

        ENDIF

   ENDIF

   RETURN

 END FUNCTION RGCT_Routine

!====================================================================
! Function Name: RVCT_Routine
!
! Function:
!   Relative Visible Contrast Test (RVCT)
!
! Description:
!   Computes the Relative Visible Contrast Test (RGCT) (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RVCT_Routine( &
!                        Is_Coast, &
!                        Is_Snow, &
!                        Is_Land, &
!                        Refl_Chn2, &
!                        Scat_Zen, &
!                        Sol_Zen, &
!                        Refl_Chn2_Clear_Stddev_3x3, &
!                        Refl_Chn2_Min_3x3)
!
! Inputs:
!   Is pixel coast (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel land(YES/NO)
!   Is pixel glint (YES/NO)
!   0.64 miron reflectance
!   Scattering zenith angle
!   Solar zenith angle
!   Standard Deviation of the clear sky 0.64 micron reflectance over a 3x3 box
!   Minima of the 0.64 micron reflectance over a 3x3 box
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION RVCT_Routine( &
                        Is_Coast, &
                        Is_Snow, &
                        Is_Land, &
                        Refl_Chn2, &
                        Scat_Zen, &
                        Sol_Zen, &
                        Refl_Chn2_Clear_Stddev_3x3, &
                        Refl_Chn2_Min_3x3) &
                        RESULT(Test_Result)

   INTEGER(KIND=INT4),INTENT(IN) :: Is_Coast
   INTEGER(KIND=INT4),INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4),INTENT(IN) :: Is_Land
   REAL(KIND=REAL4),INTENT(IN) :: Refl_Chn2
   REAL(KIND=REAL4),INTENT(IN) :: Scat_Zen
   REAL(KIND=REAL4),INTENT(IN) :: Sol_Zen
   REAL(KIND=REAL4),INTENT(IN) :: Refl_Chn2_Clear_Stddev_3x3
   REAL(KIND=REAL4),INTENT(IN) :: Refl_Chn2_Min_3x3
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold

   !--- initialize
   Test_Result = sym%NO

   !--- derive threshold
   Test_Threshold = 999.0
   IF ( (Is_Coast == sym%NO) .AND. (Is_Snow == sym%NO) .AND.  &
        (Scat_Zen > RVCT_SCAT_ZEN_THRESH) .AND. &
        (Sol_Zen < RVCT_SOL_ZEN_THRESH)) THEN

      Test_Threshold = 5.0 

      IF (Is_Land == sym%YES) THEN
           Test_Threshold = Test_Threshold + 5.0 + 4.0*Refl_Chn2_Clear_Stddev_3x3 
      ENDIF

      IF (Is_Land == sym%NO .OR. Refl_Chn2_Clear_Stddev_3x3 <= 0.0) THEN
           Test_Threshold = 10.0
      ENDIF

      !--- apply test
      IF (Refl_Chn2 - Refl_Chn2_Min_3x3 > Test_Threshold) THEN
           Test_Result = sym%YES
      ENDIF

   ENDIF

   RETURN

 END FUNCTION RVCT_Routine
 



!====================================================================
! Function Name: NIRREF_Chn5_Routine
!
! Function:
!   1.6 micron (Ch4) Near-IR Reflectance Test
!
! Description:
!   Computes the 1.6 micron (Ch4) Near-IR Reflectance Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  NIRREF_Chn5_Routine( &
!                         Is_Coast, &
!                         Is_Snow, &
!                         Sol_Zen, &
!                         Sfc_Hgt, &
!                         NDSI, &
!                         Refl_Chn5)
!
! Inputs:
!   Is pixel coast (YES/NO)
!   Is pixel snow (YES/NO)
!   Solar zenith angle
!   Surface height
!   Normalized Difference Snow Index
!   1.5 micron reflection
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION NIRREF_Chn5_Routine(&
                         Is_Coast, &
                         Is_Snow, &
                         Sol_Zen, &
                         Sfc_Hgt, &
                         NDSI, &
                         Refl_Chn5) &
                         RESULT(Test_Result)

   INTEGER(KIND=INT4),INTENT(IN) :: Is_Coast
   INTEGER(KIND=INT4),INTENT(IN) :: Is_Snow
   REAL(KIND=REAL4),INTENT(IN) :: Sol_Zen
   REAL(KIND=REAL4),INTENT(IN) :: Sfc_Hgt
   REAL(KIND=REAL4),INTENT(IN) :: NDSI
   REAL(KIND=REAL4),INTENT(IN) :: Refl_Chn5
   INTEGER(KIND=INT1) :: Test_Result

    !--- initialize
    Test_Result = sym%NO

    !--- apply test
    IF ( (Is_Coast == sym%NO) .AND.  &
         (Is_Snow == sym%YES) .AND.  &
         (Sol_Zen < NIRREF_CHN5_SOL_ZEN_THRESH) .AND.  &
         (NDSI < NIRREF_NDSI_THRESH_SNOW) .AND.  &
         (Sfc_Hgt < NIRREF_SFC_HGT_LIMIT)) THEN

         IF (Refl_Chn5 > NIRREF_CHN5_REFL_THRESH_SNOW) THEN

            Test_Result = sym%YES

         ENDIF

   ENDIF

   RETURN

 END FUNCTION NIRREF_Chn5_Routine
 
!====================================================================
! Function Name: NIRREF_Chn7_Routine
!
! Function:
!   3.9 micron (Ch7) Near-IR Reflectance Test
!
! Description:
!   Computes the 3.9 micron (Ch7) Near-IR Reflectance Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  NIRREF_Chn7_Routine( &
!                         Is_Snow, &
!                         Sol_Zen, &
!                         Sfc_Hgt, &
!                         Refl_Chn7)
!
! Inputs:
!   Is pixel snow (YES/NO)
!   Solar zenith angle
!   Surface height
!   3.9 micron reflectance
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION NIRREF_Chn7_Routine(&
                         Is_Snow, &
                         Sol_Zen, &
                         Sfc_Hgt, &
                         Refl_Chn7) &
                         RESULT(Test_Result)

   INTEGER(KIND=INT4),INTENT(IN) :: Is_Snow
   REAL(KIND=REAL4),INTENT(IN) :: Sol_Zen
   REAL(KIND=REAL4),INTENT(IN) :: Sfc_Hgt
   REAL(KIND=REAL4),INTENT(IN) :: Refl_Chn7
   INTEGER(KIND=INT1) :: Test_Result

    !--- initialize
    Test_Result = sym%NO

    !--- apply test
    IF ( (Is_Snow == sym%YES) .AND.  &
         (Sol_Zen < NIRREF_CHN7_SOL_ZEN_THRESH) .AND.  &
         (Sfc_Hgt < NIRREF_SFC_HGT_LIMIT)) THEN

         IF (Refl_Chn7 > NIRREF_CHN7_REFL_THRESH_SNOW) THEN

            Test_Result = sym%YES

         ENDIF

   ENDIF

   RETURN

 END FUNCTION NIRREF_Chn7_Routine

!====================================================================
! Function Name: CIRREF_Routine
!
! Function:
!   Cirrus Reflectance Test (1.38 micron), CIRREF
!
! Description:
!   Computes the Cirrus Reflectance Test (1.38 micron) (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  CIRREF_Routine( &
!                         Is_Snow, &
!                         Refl_Chn4, &
!                         Sfc_Hgt_Max_3x3, &
!                         CIRREF_Sfc_Hgt_Limit, &
!                         Sol_Zen, &
!                         Refl_Test_Sol_Zen_Thresh)!
! Inputs:
!   Is pixel snow (YES/NO)
!   1.38 micron reflectance
!   Minima of the surface height over a 3x3 box
!   CIRREF surface height limit threshold
!   Solar zenith angle
!   reflectance test solar zenith threshold
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION  CIRREF_Routine( &
                         Is_Snow, &
                         Refl_Chn4, &
                         Sfc_Hgt_Max_3x3, &
                         Sol_Zen) &
                         RESULT(Test_Result)

   INTEGER(KIND=INT4) :: Is_Snow
   REAL(KIND=REAL4), INTENT(IN) :: Refl_Chn4
   REAL(KIND=REAL4), INTENT(IN) :: Sfc_Hgt_Max_3x3
   REAL(KIND=REAL4), INTENT(IN) :: Sol_Zen
   INTEGER(KIND=INT1) :: Test_Result

   !--- initialize
   Test_Result = sym%NO

   !--- apply test
   IF ( (Is_Snow == sym%NO) .AND.             &
        (Sfc_Hgt_Max_3X3 < CIRREF_SFC_HGT_LIMIT) .AND. &
        (Sol_Zen < CIRREF_SOL_ZEN_THRESH)) THEN

      IF (Refl_Chn4 > CIRREF_THRESH) THEN

         Test_Result = sym%YES

      ENDIF

   ENDIF

   RETURN

 END FUNCTION CIRREF_Routine
 

!====================================================================
! Function Name: EMISS4_Routine
!
! Function:
!   4 micron (Ch7) Emissivity (EMISS4) Cloud Test 
!
! Description:
!   Computes the 4 micron (Ch7) Emissivity (EMISS4) Cloud Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  EMISS4_Routine(Is_Glint,  &
!                         Is_Coast, &
!                         Is_Land, &
!                         Is_Snow, &
!                         Is_Desert, &
!                         BT_Chn14, &
!                         Emiss_Chn7,  &
!                         Emiss_Chn7_Clr
! Inputs:
!   Is pixel glint (YES/NO)
!   Is pixel coast(YES/NO)
!   Is pixel land (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel Desert (YES/NO)
!   11 micron BT
!   3.9 micron emissvity
!   Clear sky 3.9 micron emissivity
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================
 FUNCTION EMISS4_Routine(Is_Glint,  &
                         Is_Land, &
                         Is_Snow, &
                         Is_Desert, &
                         BT_Chn14, &
                         Emiss_Chn7,  &
                         Emiss_Chn7_Clr, &
                         Sfc_Emiss_Chn7_RTM) &
                         RESULT (Test_Result)
  INTEGER(KIND=INT4), INTENT(IN):: Is_Glint
  INTEGER(KIND=INT4), INTENT(IN):: Is_Land
  INTEGER(KIND=INT4), INTENT(IN):: Is_Snow
  INTEGER(KIND=INT4), INTENT(IN):: Is_Desert
  REAL(KIND=REAL4), INTENT(IN):: BT_Chn14
  REAL(KIND=REAL4), INTENT(IN):: Emiss_Chn7
  REAL(KIND=REAL4), INTENT(IN):: Emiss_Chn7_Clr
  REAL(KIND=REAL4), INTENT(IN):: Sfc_Emiss_Chn7_RTM
  INTEGER(KIND=INT1):: Test_Result
  REAL(KIND=REAL4):: Test_Value
  REAL(KIND=REAL4):: Test_Threshold

  Test_Result = sym%NO

   !
   !---avoid glint
   !
   IF (Is_Glint == sym%YES) THEN   
      RETURN
   ENDIF

   !
   !---avoid warm pixels
   !
   IF (BT_Chn14 > EMISS4_BT_CHN14_MAX_THRESH) THEN   
      RETURN
   ENDIF

   !--- Determine threshold based on surface condition
   Test_Threshold = EMISS4_Emiss_Chn7_Ocn_Thresh
   IF (Is_Land == sym%Yes) THEN
       Test_Threshold = EMISS4_Emiss_Chn7_Land_Thresh
   ENDIF
   IF (Is_Desert == sym%Yes) THEN
       Test_Threshold = EMISS4_Emiss_Chn7_Desert_Thresh
   ENDIF
   IF (Is_Snow == sym%Yes) THEN
       Test_Threshold = EMISS4_Emiss_Chn7_Snow_Thresh
   ENDIF

   !--- Augment in presence of very low emissive surfaces
   IF (Sfc_Emiss_Chn7_RTM < EMISS4_SFC_EMISS_CHN7_THRESH) then
       Test_Threshold = Test_Threshold + 0.5
   ENDIF

   !--- Compute value to test
   Test_Value = (Emiss_Chn7 - Emiss_Chn7_Clr) / Emiss_Chn7_Clr

   !--- apply test
   IF (Test_Value > Test_Threshold) THEN

       Test_Result = sym%YES

   ENDIF


   RETURN

END FUNCTION EMISS4_Routine


!====================================================================
! Function Name: ULST_Routine
!
! Function:
!   UnIForm Low Stratus Test - ULST
!
! Description:
!   Computes the UnIForm Low Stratus Test (ULST) (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  ULST_Routine(&
!                       Is_Land, &
!                       Is_Day, &
!                       Is_Snow, &
!                       BT_Chn14,  &
!                       BT_Chn14_Stddev_3x3,  &
!                       Emiss_Chn7,  &
!                       Emiss_Chn7_Clr, &
!                       Emiss_Tropo_Chn14)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel day (YES/NO)
!   Is pixel snow (YES/NO)
!   11 micron BT
!   Standard Deviation of the 11 micron BT over a 3x3 box
!   3.9 micron emissvity
!   Clear sky 3.9 micron emissivity
!   11 micron tropopause emissvity
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION ULST_Routine(&
                       Is_Land, &
                       Is_Day, &
                       Is_Snow, &
                       Is_Cold_Surface, &
                       BT_Chn14,  &
                       Emiss_Chn7,  &
                       Emiss_Chn7_Clr, &
                       Emiss_Chn7_NWC, &
                       Sfc_Emiss_Chn7) &
                       RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Day
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Cold_Surface
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: Emiss_Chn7
   REAL(KIND=REAL4), INTENT(IN) :: Emiss_Chn7_Clr
   REAL(KIND=REAL4), INTENT(IN) :: Emiss_Chn7_NWC
   REAL(KIND=REAL4), INTENT(IN) :: Sfc_Emiss_Chn7
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold 

   !--- initialize
   Test_Result = sym%NO

   !--- exclude day time data
   IF (Is_Day == sym%YES) THEN
        RETURN
   ENDIF

   !--- exclude day time data
   IF (Is_Cold_Surface == sym%YES) THEN
        RETURN
   ENDIF

   !--- exclude very warm pixels
   IF (BT_Chn14 > ULST_BT_CHN14_MAX_THRESH) THEN
        RETURN
   ENDIF

   !--- exclude pixels with large Emiss_Chn7
   IF (Emiss_Chn7 > ULST_EMISS_CHN7_MAX_THRESH) THEN
        RETURN
   ENDIF

   !--- exclude deserts
   IF (Sfc_Emiss_Chn7 < ULST_EMISS_CHN7_SFC_THRESH) THEN
        RETURN
   ENDIF

   !--- exclude values whose clear values are suspect
   IF ((Emiss_Chn7_Clr < ULST_EMISS_CHN7_CLR_MIN_THRESH) .OR. &
       (Emiss_Chn7_Clr > ULST_EMISS_CHN7_CLR_MAX_THRESH)) THEN
        RETURN
   ENDIF

!  !--- Apply Test for Difference with NWC
   Test_Threshold = ULST_EMISS_CH7_NWC_THRESH
   IF (( Emiss_Chn7_NWC > 0) .AND. &
       ( Emiss_Chn7_NWC - Emiss_Chn7 > Test_Threshold)) THEN
       Test_Result = sym%YES
   ENDIF

   !--- Set Threshold
   Test_Threshold = ULST_EMISS_CHN7_DIFF_THRESH_OCEAN
   IF (Is_Land == sym%YES) THEN
        Test_Threshold = ULST_EMISS_CHN7_DIFF_THRESH_LAND
   ENDIF
   IF (Is_Snow == sym%YES) THEN
        Test_Threshold = ULST_EMISS_CHN7_DIFF_THRESH_SNOW
   ENDIF

   !--- Apply Test for DIFferences Relative to Clear Sky
   IF (Emiss_Chn7_Clr - Emiss_Chn7 > Test_Threshold) THEN
        Test_Result = sym%YES
   ENDIF


   !--- Bulk Threshold
   IF (Emiss_Chn7 < 0.80) THEN
        Test_Result = sym%YES
   ENDIF

   RETURN

 END FUNCTION ULST_Routine

!====================================================================
! Function Name: TERM_THERM_STAB_Routine
!
! Function:
!    Maia terminator stability test
!
! Description:
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = Term_Therm_Stab_Routine(&
!            Is_Land, &
!            Sol_Zen, &
!            BT_Chn11, &
!            BT_Chn11_1Hr, &
!            BT_Chn14, &
!            BT_Chn14_1Hr, &
!            BT_Chn15, &
!            BT_Chn15_1Hr, &
!            Cmask_1Hr)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Solar Zenith Angle
!   Current and 1 hour previous 8.5 micron BT 
!   Current and 1 hour previous 11 micron BT 
!   Current and 1 hour previous 12 micron BT 
!   1 hour previous Cloud Mask
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
! Reference: MAIA paper (need to insert)
!
!====================================================================
 FUNCTION Term_Therm_Stab_Routine(&
             Is_Land, &
             Sol_Zen, &
             BT_Chn11, &
             BT_Chn11_1Hr, &
             BT_Chn14, &
             BT_Chn14_1Hr, &
             BT_Chn15, &
             BT_Chn15_1Hr, &
             Cmask_1Hr) &
      RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN):: Is_Land
   REAL(KIND=REAL4), INTENT(IN):: Sol_Zen
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn11
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn11_1Hr
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn14_1Hr
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn15
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn15_1Hr
   INTEGER(KIND=INT4), INTENT(IN):: Cmask_1Hr
   INTEGER(KIND=INT1):: Test_Result
   REAL(KIND=REAL4):: BT14_DIFf 
   REAL(KIND=REAL4):: Test_Value 
   REAL(KIND=REAL4):: Threshold


   !--- initialize
   Test_Result = sym%NO

   !--- exclude day time data
   IF (Sol_Zen < TERM_THERM_STAB_SOLZEN_MIN_THRESH .OR.  &
       Sol_Zen > TERM_THERM_STAB_SOLZEN_MAX_THRESH) THEN
        RETURN
   ENDIF
   
   !--- exclude large B14_DIFF
   BT14_DIFf = ABS(BT_Chn14 - BT_Chn14_1Hr)
   IF (BT14_DIFf > TERM_THERM_STAB_BT14_DIFF_THRESH) THEN
        RETURN
   ENDIF
  
   IF (Is_Land == sym%YES) THEN
        Test_Value = ABS((BT_Chn14_1Hr - BT_Chn11_1Hr) - &
            (BT_Chn14 - BT_Chn11))
        Threshold = TERM_THERM_STAB_BTD_14_11_THRESH
   
   ELSE
        Test_Value = ABS((BT_Chn14_1Hr - BT_Chn15_1Hr) - &
            (BT_Chn14 - BT_Chn15))
        Threshold = TERM_THERM_STAB_BTD_14_15_THRESH
   ENDIF
   
   !--- Apply Test 
   IF (Test_Value < Threshold .and. Cmask_1Hr == sym%CLOUDY) THEN
        Test_Result = sym%YES
   ENDIF


   RETURN

 END FUNCTION Term_Therm_Stab_Routine


!====================================================================
! Function Name: CIRH2O_Routine
!
! Function:
!    Maia terminator stability test
!
! Description:
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = CIRH2O_Routine(&
!
! Inputs:
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
! Reference: None
!
!====================================================================


 FUNCTION CIRH2O_Routine( &
                         Sfc_Hgt, &
                         WaterVapor_Window_Correlation, &
                         WaterVapor_Stddev, &
                         Window_Stddev, &
                         Total_Precipitable_Water, &
                         Cos_Sen_Zen)&
                         RESULT(Test_Result) 

   REAL(KIND=REAL4), INTENT(IN):: Sfc_Hgt
   REAL(KIND=REAL4), INTENT(IN):: WaterVapor_Window_Correlation
   REAL(KIND=REAL4), INTENT(IN):: WaterVapor_Stddev
   REAL(KIND=REAL4), INTENT(IN):: Window_Stddev
   REAL(KIND=REAL4), INTENT(IN):: Total_Precipitable_Water
   REAL(KIND=REAL4), INTENT(IN):: Cos_Sen_Zen
   INTEGER(KIND=INT1):: Test_Result

   Test_Result = sym%NO

   IF (Total_Precipitable_Water / Cos_Sen_Zen < CIRH2O_TPW_THRESH) THEN
      RETURN
   ENDIF

   IF (WaterVapor_Window_Correlation > CIRH2O_CORRELATION_THRESH .AND.  &
       WaterVapor_Stddev > CIRH2O_BT_CHN10_STDDEV_THRESH .AND. &
       Window_Stddev > CIRH2O_BT_CHN14_STDDEV_THRESH .AND. &
       Sfc_Hgt < CIRH2O_SFC_HGT_LIMIT) THEN
        Test_Result = sym%YES
   ENDIF


 END FUNCTION CIRH2O_Routine 
 
!====================================================================
! Function Name: MECIDA_Routine
!
! Function:
!   MeCiDA Water Vapor Test
!
! Description:
!   Computes the Meteosat Second Generation Cirrus Detection Algorithm (MeCiDA) test #5 (YES/NO), 
!   as described in the ABI cloud mask ATBD and 
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  MECIDA_Routine( &
!                         BT_Chn10_Mean_15x15, &
!                         BT_Chn10, &
!                         BT_Chn16, &
!                         MECIDA_Chn10_Low_Pass_Thresh, &
!                         MECIDA_Chn16_Thresh1,   &
!                         MECIDA_Chn16_Thresh2)
!
! Inputs:
!   Mean of the 7.4 micron BT over a 3x3 box
!   7.4 micron BT
!   13 micron BT
!   MeCiDA 7.4 micron low pass threshold
!   MeCiDA 13 micron, 1st threshold
!   MeCiDA 13 micron, 2nd threshold
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
! Reference: Krebs, W., Mannstein, H., Bugliaro, L., and Mayer, B.: Technical note: 
!            A new day- and night-time Meteosat Second Generation Cirrus Detection
!            Algorithm MeCiDA, Atmos. Chem. Phys.,7, 6145-6159, 2007
!
!====================================================================
 FUNCTION MECIDA_Routine( &
                         BT_Chn10_Mean_15x15, &
                         BT_Chn10, &
                         BT_Chn16, &
                         MECIDA_Chn10_Low_Pass_Thresh, &
                         MECIDA_Chn16_Thresh1,   &
                         MECIDA_Chn16_Thresh2)   &
                         RESULT(Test_Result)

   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn10
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn16
   REAL(KIND=REAL4), INTENT(IN) :: MECIDA_Chn10_Low_Pass_Thresh
   REAL(KIND=REAL4), INTENT(IN) :: MECIDA_Chn16_Thresh1
   REAL(KIND=REAL4), INTENT(IN) :: MECIDA_Chn16_Thresh2
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn10_Mean_15x15
   INTEGER(KIND=INT1) :: Test_Result

   !--- initialize output
   Test_Result = sym%NO

   IF ((BT_Chn10 /= Missing_Value_REAL4) .AND. &
       (BT_Chn16 /= Missing_Value_REAL4)) THEN

      !--- apply test
      IF ((BT_Chn10_Mean_15x15 - BT_Chn10) > MECIDA_Chn10_Low_Pass_Thresh .AND. &
          (BT_CHn16 < MECIDA_Chn16_Thresh1)) THEN

        Test_Result = sym%YES

      ENDIF

      IF (BT_CHn16 < MECIDA_Chn16_Thresh2) THEN

        Test_Result = sym%YES

      ENDIF

   ENDIF

   RETURN

 END FUNCTION MECIDA_Routine


!====================================================================
! Function Name: MODIS_Land_BTD_Routine
!
! Function:
!   MODIS night time BTD difference using the 11 micron and 3.9 micron BTD
!
! Description:
!   Computes the MODIS night time BTD difference test #5
!   as described in the ABI cloud mask ATBD and the MOD35 ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  MODIS_Land_BTD_Routine( &
!                         Is_Land, &
!                         Is_DAY, &
!                         BT_Chn7, &
!                         BT_Chn14, &
!                         TPW)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel day (YES/NO)
!   3.9 micron BT
!   11 micron BT
!   TPW from NWP
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
! Reference: MODIS MOD35 Collection 6 ATBD
!
!====================================================================

 FUNCTION MODIS_Land_BTD_Routine(&
                       Is_Land, &
                       Is_Day, &
                       BT_Chn7,  &
                       BT_Chn14,  &
                       TPW) &
                       RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Day
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn7
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: TPW
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Value 
   REAL(KIND=REAL4) :: Test_Threshold 

   Test_Result = sym%NO
   
   Test_Value = BT_Chn14 - BT_Chn7
   
   Test_Threshold = (-0.0077) + (1.1234 *TPW) - (0.3403 * (TPW**2)) - 2

   !
   !---Only a night time test
   !
   IF (Is_Day == sym%YES) THEN   
      RETURN
   ENDIF

   !
   !---Only a land test
   !
   IF (Is_Land == sym%NO) THEN   
      RETURN
   ENDIF
   
   IF (Test_Value < Test_Threshold) THEN
   
        Test_Result = sym%YES
   
   ENDIF

   RETURN

 END FUNCTION MODIS_Land_BTD_Routine

!====================================================================
! Function Name: Pearson_Corr
!
! Function:
!    Compute the Pearson Correlation Coefficient for two mxn arrays
!
! Description: Pearson's product-moment coefficient
!   
! Calling Sequence: BT_WV_BT_Window_Corr(Elem_Idx,Line_Idx) = Pearson_Corr( &
!                       sat%bt10(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                       sat%bt14(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                       sat%Bad_Pixel_Mask(10,Array_Right:Array_Left,Array_Top:Array_Bottom), &
!                       sat%Bad_Pixel_Mask(14,Array_Right:Array_Left,Array_Top:Array_Bottom), &
!                      Array_Width, Array_Hgt)
!   
!
! Inputs:
!   Array 1
!   Array 2
!   Elem_size
!   Line_size
!
! Outputs: 
!   Pearson Correlation coefficent
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: Standard definition for Pearson correlation
!
!====================================================================

FUNCTION Pearson_Corr(Array_One,Array_Two,Bad_Pixel_One, &
                      Bad_Pixel_Two, &
                      Array_Width,Array_Hght) RESULT(Pearson_Corr_Coeff)
   REAL(KIND=REAL4), INTENT(IN), DIMENSION(:,:):: Array_One
   REAL(KIND=REAL4), INTENT(IN), DIMENSION(:,:):: Array_Two
   INTEGER(KIND=INT1), INTENT(IN), DIMENSION(:,:):: Bad_Pixel_One
   INTEGER(KIND=INT1), INTENT(IN), DIMENSION(:,:):: Bad_Pixel_Two
   INTEGER(KIND=INT4), INTENT(IN):: Array_Width
   INTEGER(KIND=INT4), INTENT(IN):: Array_Hght
   REAL(KIND=REAL4), DIMENSION(Array_Width,Array_Hght):: Pearson_Corr_Term_1
   REAL(KIND=REAL4), DIMENSION(Array_Width,Array_Hght):: Pearson_Corr_Term_2
   REAL(KIND=REAL8):: Pearson_Corr_Top_Term_1
   REAL(KIND=REAL8):: Pearson_Corr_Top_Term_2
   REAL(KIND=REAL8):: Pearson_Corr_Bottom_Term_1
   REAL(KIND=REAL8):: Pearson_Corr_Bottom_Term_2
   REAL(KIND=REAL4):: Pearson_Corr_Coeff
   REAL(KIND=REAL8):: Mean_Array_One
   REAL(KIND=REAL8):: Mean_Array_Two
   REAL(KIND=REAL8):: Sum_Array_One
   REAL(KIND=REAL8):: Sum_Array_Two

   !--- skip computation for pixel arrays with any missing data
   IF (sum(Bad_Pixel_One) > 0 .OR. sum(Bad_Pixel_Two) > 0) THEN
      Pearson_Corr_Coeff = Missing_Value_Real4
      RETURN
   ENDIF


   Sum_Array_One = sum(Array_One)
   Sum_Array_Two = sum(Array_Two)

   Mean_Array_One = Sum_Array_One / (Array_Width*Array_Hght)
   Mean_Array_Two = Sum_Array_Two / (Array_Width*Array_Hght)

   Pearson_Corr_Term_1 = Array_One - Mean_Array_One
   Pearson_Corr_Term_2 = Array_Two - Mean_Array_Two

   Sum_Array_One = sum(Pearson_Corr_Term_1)
   Sum_Array_Two = sum(Pearson_Corr_Term_2)

   Mean_Array_One = 0.0
   Mean_Array_Two = 0.0

   Pearson_Corr_Top_Term_1 = sum(Pearson_Corr_Term_1*Pearson_Corr_Term_2)
   
   Pearson_Corr_Top_Term_2 = (Sum_Array_One*Sum_Array_Two) / (Array_Width*Array_Hght)
   
   Pearson_Corr_Bottom_Term_1 = sum(Pearson_Corr_Term_1**2) - &
                                ((Sum_Array_One)**2) / (Array_Width*Array_Hght)
                                 
   Pearson_Corr_Bottom_Term_2 = sum(Pearson_Corr_Term_2**2) - &
                                ((Sum_Array_Two)**2) / (Array_Width*Array_Hght)

   Pearson_Corr_Coeff = (Pearson_Corr_Top_Term_1 - Pearson_Corr_Top_Term_2) / &
                         sqrt(Pearson_Corr_Bottom_Term_1 * &
                              Pearson_Corr_Bottom_Term_2)
   
 END FUNCTION Pearson_Corr
 
!====================================================================
! Function Name: Term_Refl_Norm
!
! Function:
!    Renormalize reflectances to improve performance near the terminator 
! using the parameteization given by Li et. al. 2006
!
! Description: Renormalizes reflectances in the terminator region
!   
! Calling Sequence: Refl_Chn2 = Term_Refl_Norm(Cos_Sol_Zen,Refl_Chn2)
!   
!
! Inputs:
!   Cosine of the Solar Zenith Angle
!   Channel 2 reflectance
!
! Outputs: 
!   Renormalized reflectance
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: Li et. al. 2006
!
!====================================================================
 FUNCTION Term_Refl_Norm(Cos_Sol_Zen,Reflectance)  &
          RESULT(Reflectance_Normalized)

   REAL(KIND=REAL4), INTENT(IN):: Cos_Sol_Zen
   REAL(KIND=REAL4), INTENT(IN):: Reflectance
   REAL(KIND=REAL4):: Reflectance_Normalized
   REAL(KIND=REAL4):: Norm_Param

   Reflectance_Normalized = Reflectance * Cos_Sol_Zen

   Norm_Param = 24.35 / (2*Cos_Sol_Zen + sqrt(498.5225*(Cos_Sol_Zen**2) + 1) )

   Reflectance_Normalized = Reflectance_Normalized*Norm_Param 

 END FUNCTION Term_Refl_Norm
