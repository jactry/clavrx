! Created on Wed Mar  8 20:52:42 2017
! by awalther
! on machine: kuu.local
select case(trim(name))
case("scan_line_number")
   if (allocated ( scan_number)) then
      data_dim1_dtype3 = scan_number(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("scan_line_time")
   if (allocated ( utc_scan_time_hours)) then
      data_dim1_dtype4 = utc_scan_time_hours(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("bad_scan_line_flag")
   if (allocated ( bad_scan_flag)) then
      data_dim1_dtype1 = bad_scan_flag(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("asc_des_flag")
   if (allocated ( nav%ascend)) then
      data_dim1_dtype1 = nav%ascend(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("bad_pixel_mask")
   if (allocated ( bad_pixel_mask)) then
      data_dim2_dtype1 = bad_pixel_mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("gap_pixel_mask")
   if (allocated ( gap_pixel_mask)) then
      data_dim2_dtype1 = gap_pixel_mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("diagnostic_1")
   if (allocated ( diag_pix_array_1)) then
      data_dim2_dtype4 = diag_pix_array_1(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("diagnostic_2")
   if (allocated ( diag_pix_array_2)) then
      data_dim2_dtype4 = diag_pix_array_2(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("diagnostic_3")
   if (allocated ( diag_pix_array_3)) then
      data_dim2_dtype4 = diag_pix_array_3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("packed_pixel_meta_data")
case("latitude")
   if (allocated ( nav%lat)) then
      data_dim2_dtype_r4 = nav%lat(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("longitude")
   if (allocated ( nav%lon)) then
      data_dim2_dtype_r4 = nav%lon(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("latitude_pc")
   if (allocated ( nav%lat_pc)) then
      data_dim2_dtype_r4 = nav%lat_pc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("longitude_pc")
   if (allocated ( nav%lon_pc)) then
      data_dim2_dtype_r4 = nav%lon_pc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("sensor_zenith_angle")
   if (allocated ( geo%satzen)) then
      data_dim2_dtype_r4 = geo%satzen(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("solar_zenith_angle")
   if (allocated ( geo%solzen)) then
      data_dim2_dtype_r4 = geo%solzen(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("relative_azimuth_angle")
   if (allocated ( geo%relaz)) then
      data_dim2_dtype_r4 = geo%relaz(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("solar_azimuth_angle")
   if (allocated ( geo%solaz)) then
      data_dim2_dtype_r4 = geo%solaz(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("sensor_azimuth_angle")
   if (allocated ( geo%sataz)) then
      data_dim2_dtype_r4 = geo%sataz(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("glint_zenith_angle")
   if (allocated ( geo%glintzen)) then
      data_dim2_dtype_r4 = geo%glintzen(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("scattering_angle")
   if (allocated ( geo%scatangle)) then
      data_dim2_dtype_r4 = geo%scatangle(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("lunar_zenith_angle")
   if (allocated ( geo%lunzen)) then
      data_dim2_dtype_r4 = geo%lunzen(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("lunar_relative_azimuth_angle")
   if (allocated ( geo%lunrelaz)) then
      data_dim2_dtype_r4 = geo%lunrelaz(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("lunar_azimuth_angle")
   if (allocated ( geo%lunaz)) then
      data_dim2_dtype_r4 = geo%lunaz(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("glint_mask")
   if (allocated ( sfc%glint_mask)) then
      data_dim2_dtype1 = sfc%glint_mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("coast_mask")
   if (allocated ( sfc%coast_mask)) then
      data_dim2_dtype1 = sfc%coast_mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("surface_type")
   if (allocated ( sfc%sfc_type)) then
      data_dim2_dtype1 = sfc%sfc_type(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("land_class")
   if (allocated ( sfc%land)) then
      data_dim2_dtype1 = sfc%land(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("snow_class")
   if (allocated ( sfc%snow)) then
      data_dim2_dtype1 = sfc%snow(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("surface_elevation")
   if (allocated ( sfc%zsfc)) then
      data_dim2_dtype_r4 = sfc%zsfc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom")
   if (allocated ( ch(1)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(1) ) cycle 
      data_dim2_dtype_r4 = ch(1)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_86um_nom")
   if (allocated ( ch(2)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(2) ) cycle 
      data_dim2_dtype_r4 = ch(2)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_47um_nom")
   if (allocated ( ch(3)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(3) ) cycle 
      data_dim2_dtype_r4 = ch(3)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_55um_nom")
   if (allocated ( ch(4)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(4) ) cycle 
      data_dim2_dtype_r4 = ch(4)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_1_2um_nom")
   if (allocated ( ch(5)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(5) ) cycle 
      data_dim2_dtype_r4 = ch(5)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_1_60um_nom")
   if (allocated ( ch(6)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(6) ) cycle 
      data_dim2_dtype_r4 = ch(6)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_2_1um_nom")
   if (allocated ( ch(7)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(7) ) cycle 
      data_dim2_dtype_r4 = ch(7)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_41um_nom")
   if (allocated ( ch(8)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(8) ) cycle 
      data_dim2_dtype_r4 = ch(8)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_44um_nom")
   if (allocated ( ch(9)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(9) ) cycle 
      data_dim2_dtype_r4 = ch(9)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_1_38um_nom")
   if (allocated ( ch(26)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(26) ) cycle 
      data_dim2_dtype_r4 = ch(26)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_sol_dnb_nom")
   if (allocated ( ch(44)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(44) ) cycle 
      data_dim2_dtype_r4 = ch(44)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_lunar_dnb_nom")
   if (allocated ( ch(44)%ref_lunar_toa)) then
      if (.NOT. sensor%chan_on_flag_default(44) ) cycle 
      data_dim2_dtype_r4 = ch(44)%ref_lunar_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_3_75um_nom")
   if (allocated ( ch(20)%ref_toa)) then
      if (.NOT. sensor%chan_on_flag_default(20) ) cycle 
      data_dim2_dtype_r4 = ch(20)%ref_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_3_75um_nom")
   if (allocated ( ch(20)%bt_toa)) then
      if (.NOT. sensor%chan_on_flag_default(20) ) cycle 
      data_dim2_dtype_r4 = ch(20)%bt_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_3_9um_nom")
   if (allocated ( ch(22)%Bt_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(22) ) cycle 
      data_dim2_dtype_r4 = ch(22)%Bt_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_8_5um_nom")
   if (allocated ( ch(27)%bt_toa)) then
      if (.NOT. sensor%chan_on_flag_default(27) ) cycle 
      data_dim2_dtype_r4 = ch(27)%bt_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_11_0um_nom")
   if (allocated ( ch(31)%bt_toa)) then
      if (.NOT. sensor%chan_on_flag_default(31) ) cycle 
      data_dim2_dtype_r4 = ch(31)%bt_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_12_0um_nom")
   if (allocated ( ch(32)%bt_toa)) then
      if (.NOT. sensor%chan_on_flag_default(32) ) cycle 
      data_dim2_dtype_r4 = ch(32)%bt_toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom_min_3x3")
   if (allocated ( ref_Ch1_Min_3x3)) then
      if ( .NOT. sensor%chan_on_flag_default(1) ) cycle 
      data_dim2_dtype_r4 = ref_Ch1_Min_3x3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom_stddev_3x3")
   if (allocated ( Ref_Ch1_Std_3x3)) then
      if ( .NOT. sensor%chan_on_flag_default(1) ) cycle 
      data_dim2_dtype_r4 = Ref_Ch1_Std_3x3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_11_0um_nom_stddev_3x3")
   if (allocated ( Ref_Ch1_Std_3x3)) then
      if ( .NOT. sensor%chan_on_flag_default(1) ) cycle 
      data_dim2_dtype_r4 = Ref_Ch1_Std_3x3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_11_0um_nom_max_3x3")
   if (allocated ( Bt_Ch31_Max_3x3)) then
      if ( .NOT. sensor%chan_on_flag_default(31) ) cycle 
      data_dim2_dtype_r4 = Bt_Ch31_Max_3x3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_probability")
   if (allocated ( CLDMASK%posterior_cld_probability)) then
      data_dim2_dtype_r4 = CLDMASK%posterior_cld_probability(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_mask")
   if (allocated ( CLDMASK%cld_mask)) then
      data_dim2_dtype1 = CLDMASK%cld_mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("adj_pix_cloud_mask")
   if (allocated ( CLDMASK%Adj_Pix_Cld_Mask)) then
      data_dim2_dtype1 = CLDMASK%Adj_Pix_Cld_Mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_mask_test_packed_results")
case("bayes_mask_sfc_type")
   if (allocated ( CLDMASK%Bayes_Mask_Sfc_Type)) then
      data_dim2_dtype1 = CLDMASK%Bayes_Mask_Sfc_Type(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("dust_mask")
   if (allocated ( dust_mask)) then
      data_dim2_dtype1 = dust_mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("smoke_mask")
   if (allocated ( smoke_mask)) then
      data_dim2_dtype1 = smoke_mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("shadow_mask")
   if (allocated ( shadow_mask)) then
      data_dim2_dtype1 = shadow_mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("fire_mask")
   if (allocated ( fire_mask)) then
      data_dim2_dtype1 = fire_mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_type")
   if (allocated ( cld_type)) then
      data_dim2_dtype1 = cld_type(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_phase")
   if (allocated ( cld_phase)) then
      data_dim2_dtype1 = cld_phase(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_phase_aux")
   if (allocated ( Cld_Phase_Aux)) then
      data_dim2_dtype1 = Cld_Phase_Aux(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_mask_aux")
   if (allocated ( CLDMASK%Cld_Mask_Aux)) then
      data_dim2_dtype1 = CLDMASK%Cld_Mask_Aux(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_type_aux")
   if (allocated ( Cld_Type_Aux)) then
      data_dim2_dtype1 = Cld_Type_Aux(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_height_aux")
   if (allocated ( zc_aux)) then
      data_dim2_dtype_r4 = zc_aux(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_height_sndr")
   if (allocated ( zc_co2)) then
      data_dim2_dtype_r4 = zc_co2(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_press_acha")
   if (allocated ( acha%pc)) then
      data_dim2_dtype_r4 = acha%pc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_temp_acha")
   if (allocated ( acha%tc)) then
      data_dim2_dtype_r4 = acha%tc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_height_acha")
   if (allocated ( acha%zc)) then
      data_dim2_dtype_r4 = acha%zc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_height_top_acha")
   if (allocated ( acha%zc_top)) then
      data_dim2_dtype_r4 = acha%zc_top(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_press_top_acha")
   if (allocated ( ACHA%Pc_Top)) then
      data_dim2_dtype_r4 = ACHA%Pc_Top(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_press_top_l1_aux")
   if (allocated ( Pc_Top1_Aux)) then
      data_dim2_dtype_r4 = Pc_Top1_Aux(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_press_top_l2_aux")
   if (allocated ( Pc_Top2_Aux)) then
      data_dim2_dtype_r4 = Pc_Top2_Aux(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_height_base_acha")
   if (allocated ( acha%zc_base)) then
      data_dim2_dtype_r4 = acha%zc_base(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_press_base_acha")
   if (allocated ( ACHA%Pc_Base)) then
      data_dim2_dtype_r4 = ACHA%Pc_Base(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_height_lower_acha")
   if (allocated ( ACHA%Lower_zc)) then
      data_dim2_dtype_r4 = ACHA%Lower_zc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_pressure_lower_acha")
   if (allocated ( ACHA%Lower_pc)) then
      data_dim2_dtype_r4 = ACHA%Lower_pc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_altitude_acha")
   if (allocated ( acha%alt)) then
      data_dim2_dtype_r4 = acha%alt(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_base_altitude_acha")
   if (allocated ( ACHA%Base_Alt)) then
      data_dim2_dtype_r4 = ACHA%Base_Alt(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("acha_processing_order")
   if (allocated ( acha%processing_order)) then
      data_dim2_dtype1 = acha%processing_order(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("acha_inversion_flag")
   if (allocated ( ACHA%Inversion_Flag)) then
      data_dim2_dtype1 = ACHA%Inversion_Flag(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_layer")
   if (allocated ( CCL%Cld_Layer)) then
      data_dim2_dtype1 = CCL%Cld_Layer(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("conv_cloud_probability")
   if (allocated ( ACHA%Conv_Cld_Prob)) then
      data_dim2_dtype_r4 = ACHA%Conv_Cld_Prob(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("supercooled_cloud_probability")
   if (allocated ( ACHA%Supercooled_Cld_Prob)) then
      data_dim2_dtype_r4 = ACHA%Supercooled_Cld_Prob(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_height_h2o")
   if (allocated ( Zc_H2O)) then
      data_dim2_dtype_r4 = Zc_H2O(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_height_opaque")
   if (allocated ( Zc_opaque_cloud)) then
      data_dim2_dtype_r4 = Zc_opaque_cloud(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_temp_opaque")
   if (allocated ( Tc_Opaque_Cloud)) then
      data_dim2_dtype_r4 = Tc_Opaque_Cloud(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_emiss_acha")
   if (allocated ( acha%ec)) then
      data_dim2_dtype_r4 = acha%ec(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_beta_acha")
   if (allocated ( acha%beta)) then
      data_dim2_dtype_r4 = acha%beta(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_height_uncer_acha")
   if (allocated ( acha%zc_uncertainty)) then
      data_dim2_dtype_r4 = acha%zc_uncertainty(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_temp_uncer_acha")
   if (allocated ( acha%tc_uncertainty)) then
      data_dim2_dtype_r4 = acha%tc_uncertainty(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_press_uncer_acha")
   if (allocated ( ACHA%Pc_Uncertainty)) then
      data_dim2_dtype_r4 = ACHA%Pc_Uncertainty(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_press_uncer1_aux")
   if (allocated ( Pc_Uncertainty1_Aux)) then
      data_dim2_dtype_r4 = Pc_Uncertainty1_Aux(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_press_uncer2_aux")
   if (allocated ( Pc_Uncertainty2_Aux)) then
      data_dim2_dtype_r4 = Pc_Uncertainty2_Aux(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_temp_acha_qf")
   if (allocated ( acha%oe_quality_flags)) then
      data_dim2_dtype1 = acha%oe_quality_flags(1,:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_emiss_acha_qf")
   if (allocated ( acha%oe_quality_flags)) then
      data_dim2_dtype1 = acha%oe_quality_flags(2,:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_beta_acha_qf")
   if (allocated ( acha%oe_quality_flags)) then
      data_dim2_dtype1 = acha%oe_quality_flags(3,:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_opd_acha")
   if (allocated ( acha%tau)) then
      data_dim2_dtype_r4 = acha%tau(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_reff_acha")
   if (allocated ( acha%reff)) then
      data_dim2_dtype_r4 = acha%reff(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("acha_quality")
   if (allocated ( ACHA%Packed_Quality_Flags)) then
      data_dim2_dtype1 = ACHA%Packed_Quality_Flags(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("acha_info")
   if (allocated ( ACHA%Packed_Meta_Data_Flags)) then
      data_dim2_dtype1 = ACHA%Packed_Meta_Data_Flags(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cost_acha")
   if (allocated ( ACHA%Cost)) then
      data_dim2_dtype_r4 = ACHA%Cost(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cost_aux")
   if (allocated ( cost_aux)) then
      data_dim2_dtype_r4 = cost_aux(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_opd_dcomp")
   if (allocated ( tau_dcomp)) then
      data_dim2_dtype_r4 = tau_dcomp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_opd_aux")
   if (allocated ( Tau_Aux)) then
      data_dim2_dtype_r4 = Tau_Aux(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_reff_dcomp")
   if (allocated ( reff_dcomp)) then
      data_dim2_dtype_r4 = reff_dcomp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_reff_aux")
   if (allocated ( Reff_Aux)) then
      data_dim2_dtype_r4 = Reff_Aux(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_opd_dcomp_unc")
   if (allocated ( tau_dcomp_cost)) then
      data_dim2_dtype_r4 = tau_dcomp_cost(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_reff_dcomp_unc")
   if (allocated ( reff_dcomp_cost)) then
      data_dim2_dtype_r4 = reff_dcomp_cost(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_opd_dcomp_qf")
   if (allocated ( tau_dcomp_qf)) then
      data_dim2_dtype1 = tau_dcomp_qf(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_reff_dcomp_qf")
   if (allocated ( reff_dcomp_qf)) then
      data_dim2_dtype1 = reff_dcomp_qf(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("dcomp_quality")
   if (allocated ( Dcomp_quality_flag)) then
      data_dim2_dtype1 = Dcomp_quality_flag(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("dcomp_info")
   if (allocated ( dcomp_info_flag)) then
      data_dim2_dtype2 = dcomp_info_flag(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("insolation_dcomp")
   if (allocated ( Insolation_Dcomp)) then
      data_dim2_dtype_r4 = Insolation_Dcomp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("insolation_diffuse_dcomp")
   if (allocated ( Insolation_Dcomp_Diffuse)) then
      data_dim2_dtype_r4 = Insolation_Dcomp_Diffuse(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cdnc_dcomp")
   if (allocated ( cdnc_dcomp)) then
      data_dim2_dtype_r4 = cdnc_dcomp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("hcld_dcomp")
   if (allocated ( hcld_dcomp)) then
      data_dim2_dtype_r4 = hcld_dcomp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_opd_nlcomp")
   if (allocated ( Tau_Nlcomp)) then
      data_dim2_dtype_r4 = Tau_Nlcomp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_reff_nlcomp")
   if (allocated ( Reff_Nlcomp)) then
      data_dim2_dtype_r4 = Reff_Nlcomp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_opd_nlcomp_unc")
   if (allocated ( Tau_Nlcomp_Cost)) then
      data_dim2_dtype_r4 = Tau_Nlcomp_Cost(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_reff_nlcomp_unc")
   if (allocated ( Reff_Nlcomp_Cost)) then
      data_dim2_dtype_r4 = Reff_Nlcomp_Cost(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("nlcomp_quality")
   if (allocated ( Nlcomp_quality_flag)) then
      data_dim2_dtype1 = Nlcomp_quality_flag(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("nlcomp_info")
   if (allocated ( Nlcomp_Info_Flag)) then
      data_dim2_dtype2 = Nlcomp_Info_Flag(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_albedo_0_65um_nom")
   if (allocated ( cloud_063um_albedo)) then
      data_dim2_dtype_r4 = cloud_063um_albedo(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_transmission_0_65um_nom")
   if (allocated ( cloud_063um_transmission_solar)) then
      data_dim2_dtype_r4 = cloud_063um_transmission_solar(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_fraction")
   if (allocated ( CCL%Cloud_Fraction)) then
      data_dim2_dtype_r4 = CCL%Cloud_Fraction(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("high_cloud_fraction")
   if (allocated ( CCL%high_cloud_fraction)) then
      data_dim2_dtype_r4 = CCL%high_cloud_fraction(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("mid_cloud_fraction")
   if (allocated ( CCL%mid_cloud_fraction)) then
      data_dim2_dtype_r4 = CCL%mid_cloud_fraction(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("low_cloud_fraction")
   if (allocated ( CCL%low_cloud_fraction)) then
      data_dim2_dtype_r4 = CCL%low_cloud_fraction(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_fraction_uncertainty")
   if (allocated ( CCL%Cloud_Fraction_Uncer)) then
      data_dim2_dtype_r4 = CCL%Cloud_Fraction_Uncer(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("emiss_3_75um_nom")
   if (allocated ( Ems_Ch20)) then
      if ( .NOT. sensor%chan_on_flag_default(2) ) cycle 
      data_dim2_dtype_r4 = Ems_Ch20(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("emiss_3_75um_nom_clear")
   if (allocated ( Ems_Ch20_Clear_Rtm)) then
      if ( .NOT. sensor%chan_on_flag_default(2) ) cycle 
      data_dim2_dtype_r4 = Ems_Ch20_Clear_Rtm(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("emiss_3_75um_nom_median_3x3")
   if (allocated ( Ems_Ch20_Median_3x3)) then
      if ( .NOT. sensor%chan_on_flag_default(2) ) cycle 
      data_dim2_dtype_r4 = Ems_Ch20_Median_3x3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_3_75um_nom_median_3x3")
   if (allocated ( Bt_Ch20_median_3x3)) then
      if ( .NOT. sensor%chan_on_flag_default(2) ) cycle 
      data_dim2_dtype_r4 = Bt_Ch20_median_3x3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("diff_ch31_ch32_bt_ch31_max_3x3")
   if (allocated ( Btd_Ch31_Ch32_Bt_Ch31_Max_3x3)) then
      if ( .NOT. sensor%chan_on_flag_default(31) ) cycle 
      data_dim2_dtype_r4 = Btd_Ch31_Ch32_Bt_Ch31_Max_3x3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("emissivity_11um_tropopause")
   if (allocated ( ch(31)%Emiss_Tropo)) then
      if (.NOT. sensor%chan_on_flag_default(31) ) cycle 
      data_dim2_dtype_r4 = ch(31)%Emiss_Tropo(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("beta_11um_85um_tropopause")
   if (allocated ( Beta_11um_85um_Tropo_Rtm)) then
      data_dim2_dtype_r4 = Beta_11um_85um_Tropo_Rtm(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("beta_11um_12um_tropopause")
   if (allocated ( Beta_11um_12um_Tropo_Rtm)) then
      data_dim2_dtype_r4 = Beta_11um_12um_Tropo_Rtm(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("aot_0_65um_nom")
   if (allocated ( aot1)) then
      data_dim2_dtype_r4 = aot1(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("aot_qf")
   if (allocated ( aot_qf)) then
      data_dim2_dtype1 = aot_qf(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("olr")
   if (allocated ( olr)) then
      data_dim2_dtype_r4 = olr(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("ndvi_sfc")
   if (allocated ( ndvi_sfc)) then
      data_dim2_dtype_r4 = ndvi_sfc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("ndvi_sfc_white_sky_nwp")
   if (allocated ( Ndvi_Sfc_White_Sky)) then
      data_dim2_dtype_r4 = Ndvi_Sfc_White_Sky(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("surface_temperature_retrieved")
   if (allocated ( Tsfc_Retrieved)) then
      data_dim2_dtype_r4 = Tsfc_Retrieved(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("atmos_rad_11_0um_nom_nwp")
   if (allocated ( ch(31)%Rad_Atm)) then
      if (.NOT. sensor%chan_on_flag_default(31) ) cycle 
      data_dim2_dtype_r4 = ch(31)%Rad_Atm(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("atmos_trans_11_0um_nom_nwp")
   if (allocated ( ch(31)%Trans_Atm)) then
      if (.NOT. sensor%chan_on_flag_default(31) ) cycle 
      data_dim2_dtype_r4 = ch(31)%Trans_Atm(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("sfc_emiss_3_75um_nom_nwp")
   if (allocated ( ch(20)%Sfc_Emiss)) then
      if (.NOT. sensor%chan_on_flag_default(20) ) cycle 
      data_dim2_dtype_r4 = ch(20)%Sfc_Emiss(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("sfc_emiss_11_0um_nom_nwp")
   if (allocated ( ch(31)%Sfc_Emiss)) then
      if (.NOT. sensor%chan_on_flag_default(31) ) cycle 
      data_dim2_dtype_r4 = ch(31)%Sfc_Emiss(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("surface_air_temperature_nwp")
   if (allocated ( Tair_Nwp_Pix)) then
      data_dim2_dtype_r4 = Tair_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("surface_temperature_background")
   if (allocated ( Tsfc_Nwp_Pix)) then
      data_dim2_dtype_r4 = Tsfc_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("sst_background_uni_3x3")
   if (allocated ( Sst_Anal_Uni)) then
      data_dim2_dtype_r4 = Sst_Anal_Uni(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("surface_relative_humidity_nwp")
   if (allocated ( Rh_Nwp_Pix)) then
      data_dim2_dtype_r4 = Rh_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("surface_pressure_background")
   if (allocated ( Psfc_Nwp_Pix)) then
      data_dim2_dtype_r4 = Psfc_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("mean_sealevel_pressure_background")
   if (allocated ( Pmsl_Nwp_Pix)) then
      data_dim2_dtype_r4 = Pmsl_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("k_index_nwp")
   if (allocated ( K_Index_Nwp_Pix)) then
      data_dim2_dtype_r4 = K_Index_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_water_path_nwp")
   if (allocated ( Cwp_Nwp_Pix)) then
      data_dim2_dtype_r4 = Cwp_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_fraction_nwp")
   if (allocated ( Cfrac_Nwp_Pix)) then
      data_dim2_dtype_r4 = Cfrac_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_press_nwp")
   if (allocated ( Pc_Nwp_Pix)) then
      data_dim2_dtype_r4 = Pc_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("number_cloud_layers_nwp")
   if (allocated ( Ncld_Layers_Nwp_Pix)) then
      data_dim2_dtype1 = Ncld_Layers_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_type_nwp")
   if (allocated ( Cld_Type_Nwp_Pix)) then
      data_dim2_dtype1 = Cld_Type_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("tropopause_temperature_nwp")
   if (allocated ( Ttropo_Nwp_Pix)) then
      data_dim2_dtype_r4 = Ttropo_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("lcl_nwp")
   if (allocated ( LCL_Height_Nwp_Pix)) then
      data_dim2_dtype_r4 = LCL_Height_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("ccl_nwp")
   if (allocated ( CCL_Height_Nwp_Pix)) then
      data_dim2_dtype_r4 = CCL_Height_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom_counts")
   if (allocated ( Ch1_Counts)) then
      if ( .NOT. sensor%chan_on_flag_default(1) ) cycle 
      data_dim2_dtype2 = Ch1_Counts(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_86um_nom_counts")
   if (allocated ( Ch2_Counts)) then
      if ( .NOT. sensor%chan_on_flag_default(2) ) cycle 
      data_dim2_dtype2 = Ch2_Counts(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_1_60um_nom_counts")
   if (allocated ( Ch6_Counts)) then
      if ( .NOT. sensor%chan_on_flag_default(6) ) cycle 
      data_dim2_dtype2 = Ch6_Counts(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("total_precipitable_water_nwp")
   if (allocated ( Tpw_Nwp_Pix)) then
      data_dim2_dtype_r4 = Tpw_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("total_column_ozone_nwp")
   if (allocated ( ozone_nwp_pix)) then
      data_dim2_dtype_r4 = ozone_nwp_pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_3_75um_nom_atmos_corr")
   if (allocated ( ch(20)%Ref_Sfc)) then
      if (.NOT. sensor%chan_on_flag_default(20) ) cycle 
      data_dim2_dtype_r4 = ch(20)%Ref_Sfc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom_atmos_corr")
   if (allocated ( ch(1)%Ref_Sfc)) then
      if (.NOT. sensor%chan_on_flag_default(1) ) cycle 
      data_dim2_dtype_r4 = ch(1)%Ref_Sfc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_86um_nom_atmos_corr")
   if (allocated ( ch(2)%Ref_Sfc)) then
      if (.NOT. sensor%chan_on_flag_default(2) ) cycle 
      data_dim2_dtype_r4 = ch(2)%Ref_Sfc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_1_60um_nom_atmos_corr")
   if (allocated ( ch(6)%Ref_Sfc)) then
      if (.NOT. sensor%chan_on_flag_default(6) ) cycle 
      data_dim2_dtype_r4 = ch(6)%Ref_Sfc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_2_1um_nom_atmos_corr")
   if (allocated ( ch(7)%Ref_Sfc)) then
      if (.NOT. sensor%chan_on_flag_default(7) ) cycle 
      data_dim2_dtype_r4 = ch(7)%Ref_Sfc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom_clear_sky")
   if (allocated ( ch(1)%Ref_Toa_Clear)) then
      if (.NOT. sensor%chan_on_flag_default(1) ) cycle 
      data_dim2_dtype_r4 = ch(1)%Ref_Toa_Clear(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_3_75um_nom_clear_sky")
   if (allocated ( ch(20)%Bt_Toa_Clear)) then
      if (.NOT. sensor%chan_on_flag_default(20) ) cycle 
      data_dim2_dtype_r4 = ch(20)%Bt_Toa_Clear(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_8_5um_nom_clear_sky")
   if (allocated ( ch(27)%Bt_Toa_Clear)) then
      if (.NOT. sensor%chan_on_flag_default(27) ) cycle 
      data_dim2_dtype_r4 = ch(27)%Bt_Toa_Clear(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_11_0um_nom_clear_sky")
   if (allocated ( ch(31)%Bt_Toa_Clear)) then
      if (.NOT. sensor%chan_on_flag_default(31) ) cycle 
      data_dim2_dtype_r4 = ch(31)%Bt_Toa_Clear(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_12_0um_nom_clear_sky")
   if (allocated ( ch(32)%Bt_Toa_Clear)) then
      if (.NOT. sensor%chan_on_flag_default(32) ) cycle 
      data_dim2_dtype_r4 = ch(32)%Bt_Toa_Clear(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom_mean_3x3")
   if (allocated ( Ref_Ch1_Mean_3x3)) then
      if ( .NOT. sensor%chan_on_flag_default(1) ) cycle 
      data_dim2_dtype_r4 = Ref_Ch1_Mean_3x3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("pixel_sst_unmasked")
   if (allocated ( Sst_Unmasked)) then
      data_dim2_dtype_r4 = Sst_Unmasked(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("wind_speed_10m_nwp")
   if (allocated ( Wnd_Spd_10m_Nwp_Pix)) then
      data_dim2_dtype_r4 = Wnd_Spd_10m_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("wind_direction_10m_nwp")
   if (allocated ( Wnd_Dir_10m_Nwp_Pix)) then
      data_dim2_dtype_r4 = Wnd_Dir_10m_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom_dark")
   if (allocated ( Ref_Ch1_Dark_Composite)) then
      if ( .NOT. sensor%chan_on_flag_default(1) ) cycle 
      data_dim2_dtype_r4 = Ref_Ch1_Dark_Composite(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_water_path")
   if (allocated ( Cwp_Dcomp)) then
      data_dim2_dtype_r4 = Cwp_Dcomp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("rain_rate")
   if (allocated ( Rain_Rate_Dcomp)) then
      data_dim2_dtype_r4 = Rain_Rate_Dcomp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("wind_speed_cloud_top_nwp")
   if (allocated ( Wnd_Spd_Cld_Top_Nwp_Pix)) then
      data_dim2_dtype_r4 = Wnd_Spd_Cld_Top_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("wind_direction_cloud_top_nwp")
   if (allocated ( Wnd_Dir_Cld_Top_Nwp_Pix)) then
      data_dim2_dtype_r4 = Wnd_Dir_Cld_Top_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_opd_dcomp_1")
   if (allocated ( tau_dcomp_1)) then
      data_dim2_dtype_r4 = tau_dcomp_1(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_opd_dcomp_2")
   if (allocated ( tau_dcomp_2)) then
      data_dim2_dtype_r4 = tau_dcomp_2(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_opd_dcomp_3")
   if (allocated ( tau_dcomp_3)) then
      data_dim2_dtype_r4 = tau_dcomp_3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_reff_dcomp_1")
   if (allocated ( reff_dcomp_1)) then
      data_dim2_dtype_r4 = reff_dcomp_1(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_reff_dcomp_2")
   if (allocated ( reff_dcomp_2)) then
      data_dim2_dtype_r4 = reff_dcomp_2(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_reff_dcomp_3")
   if (allocated ( reff_dcomp_3)) then
      data_dim2_dtype_r4 = reff_dcomp_3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("albedo_sfc_white_sky_0_65um")
   if (allocated ( ch(1)%sfc_ref_white_sky)) then
      if (.NOT. sensor%chan_on_flag_default(1) ) cycle 
      data_dim2_dtype_r4 = ch(1)%sfc_ref_white_sky(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("albedo_sfc_white_sky_1_60um")
   if (allocated ( ch(6)%sfc_ref_white_sky)) then
      if (.NOT. sensor%chan_on_flag_default(6) ) cycle 
      data_dim2_dtype_r4 = ch(6)%sfc_ref_white_sky(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("albedo_sfc_white_sky_2_20um")
   if (allocated ( ch(7)%sfc_ref_white_sky)) then
      if (.NOT. sensor%chan_on_flag_default(7) ) cycle 
      data_dim2_dtype_r4 = ch(7)%sfc_ref_white_sky(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_95um_nom")
   if (allocated ( ch(17)%Ref_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(17) ) cycle 
      data_dim2_dtype_r4 = ch(17)%Ref_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_93um_nom")
   if (allocated ( ch(18)%Ref_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(18) ) cycle 
      data_dim2_dtype_r4 = ch(18)%Ref_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_94um_nom")
   if (allocated ( ch(19)%Ref_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(19) ) cycle 
      data_dim2_dtype_r4 = ch(19)%Ref_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_6_7um_nom")
   if (allocated ( ch(27)%Bt_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(27) ) cycle 
      data_dim2_dtype_r4 = ch(27)%Bt_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_7_3um_nom")
   if (allocated ( ch(28)%Bt_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(28) ) cycle 
      data_dim2_dtype_r4 = ch(28)%Bt_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_9_7um_nom")
   if (allocated ( ch(30)%Bt_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(30) ) cycle 
      data_dim2_dtype_r4 = ch(30)%Bt_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_13_3um_nom")
   if (allocated ( ch(33)%Bt_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(33) ) cycle 
      data_dim2_dtype_r4 = ch(33)%Bt_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_13_6um_nom")
   if (allocated ( ch(34)%Bt_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(34) ) cycle 
      data_dim2_dtype_r4 = ch(34)%Bt_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_13_9um_nom")
   if (allocated ( ch(35)%Bt_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(35) ) cycle 
      data_dim2_dtype_r4 = ch(35)%Bt_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_14_2um_nom")
   if (allocated ( ch(36)%Bt_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(36) ) cycle 
      data_dim2_dtype_r4 = ch(36)%Bt_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("csbt_mask_6_7um_nom")
   if (allocated ( ch(27)%CSBT_Mask)) then
      if (.NOT. sensor%chan_on_flag_default(27) ) cycle 
      data_dim2_dtype1 = ch(27)%CSBT_Mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("csbt_mask_7_3um_nom")
   if (allocated ( ch(28)%CSBT_Mask)) then
      if (.NOT. sensor%chan_on_flag_default(28) ) cycle 
      data_dim2_dtype1 = ch(28)%CSBT_Mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("csbt_mask_13_3um_nom")
   if (allocated ( ch(33)%CSBT_Mask)) then
      if (.NOT. sensor%chan_on_flag_default(33) ) cycle 
      data_dim2_dtype1 = ch(33)%CSBT_Mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_11um_vs_67um_covar_5x5")
   if (allocated ( Covar_Ch27_Ch31_5x5)) then
      if ( .NOT. sensor%chan_on_flag_default(27) ) cycle 
      data_dim2_dtype_r4 = Covar_Ch27_Ch31_5x5(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom_iband_max")
   if (allocated ( Ref_Max_ChI1)) then
      data_dim2_dtype_r4 = Ref_Max_ChI1(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom_iband_min")
   if (allocated ( Ref_Min_ChI1)) then
      data_dim2_dtype_r4 = Ref_Min_ChI1(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom_iband_mean")
   if (allocated ( Ref_Mean_ChI1)) then
      data_dim2_dtype_r4 = Ref_Mean_ChI1(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_65um_nom_iband_uni")
   if (allocated ( Ref_Uni_ChI1)) then
      data_dim2_dtype_r4 = Ref_Uni_ChI1(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_86um_nom_iband_max")
   if (allocated ( Ref_Max_ChI2)) then
      data_dim2_dtype_r4 = Ref_Max_ChI2(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_86um_nom_iband_min")
   if (allocated ( Ref_Min_ChI2)) then
      data_dim2_dtype_r4 = Ref_Min_ChI2(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_86um_nom_iband_mean")
   if (allocated ( Ref_Mean_ChI2)) then
      data_dim2_dtype_r4 = Ref_Mean_ChI2(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_0_86um_nom_iband_uni")
   if (allocated ( Ref_Uni_ChI2)) then
      data_dim2_dtype_r4 = Ref_Uni_ChI2(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_1_61um_nom_iband_max")
   if (allocated ( Ref_Max_ChI3)) then
      data_dim2_dtype_r4 = Ref_Max_ChI3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_1_61um_nom_iband_min")
   if (allocated ( Ref_Min_ChI3)) then
      data_dim2_dtype_r4 = Ref_Min_ChI3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_1_61um_nom_iband_mean")
   if (allocated ( Ref_Mean_ChI3)) then
      data_dim2_dtype_r4 = Ref_Mean_ChI3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("refl_1_61um_nom_iband_uni")
   if (allocated ( Ref_Uni_ChI3)) then
      data_dim2_dtype_r4 = Ref_Uni_ChI3(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_3_74um_nom_iband_max")
   if (allocated ( Bt_Max_ChI4)) then
      data_dim2_dtype_r4 = Bt_Max_ChI4(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_3_74um_nom_iband_min")
   if (allocated ( Bt_Min_ChI4)) then
      data_dim2_dtype_r4 = Bt_Min_ChI4(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_3_74um_nom_iband_mean")
   if (allocated ( Bt_Mean_ChI4)) then
      data_dim2_dtype_r4 = Bt_Mean_ChI4(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_3_74um_nom_iband_uni")
   if (allocated ( Bt_Uni_ChI4)) then
      data_dim2_dtype_r4 = Bt_Uni_ChI4(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_11_0um_nom_iband_max")
   if (allocated ( Bt_Max_ChI5)) then
      data_dim2_dtype_r4 = Bt_Max_ChI5(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_11_0um_nom_iband_min")
   if (allocated ( Bt_Min_ChI5)) then
      data_dim2_dtype_r4 = Bt_Min_ChI5(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_11_0um_nom_iband_mean")
   if (allocated ( Bt_Mean_ChI5)) then
      data_dim2_dtype_r4 = Bt_Mean_ChI5(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("temp_11_0um_nom_iband_uni")
   if (allocated ( Bt_Uni_ChI5)) then
      data_dim2_dtype_r4 = Bt_Uni_ChI5(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("radiance_11_0um_nom")
   if (allocated ( ch(31)%Rad_Toa)) then
      if (.NOT. sensor%chan_on_flag_default(31) ) cycle 
      data_dim2_dtype_r4 = ch(31)%Rad_Toa(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_type_ir")
   if (allocated ( Cld_Type_Ir)) then
      data_dim2_dtype1 = Cld_Type_Ir(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cloud_phase_ir")
   if (allocated ( Cld_Phase_Ir)) then
      data_dim2_dtype1 = Cld_Phase_Ir(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_height_sounder")
   if (allocated ( Zc_Co2)) then
      data_dim2_dtype_r4 = Zc_Co2(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_temp_lower_acha")
   if (allocated ( ACHA%Lower_Tc)) then
      data_dim2_dtype_r4 = ACHA%Lower_Tc(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("cld_base_hgt_qf")
   if (allocated ( ACHA%base_Quality_Flag)) then
      data_dim2_dtype1 = ACHA%base_Quality_Flag(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("total_precipitable_water_nwp_above_cloud")
   if (allocated ( Tpw_Above_Cloud_Nwp_Pix)) then
      data_dim2_dtype_r4 = Tpw_Above_Cloud_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("inversion_strength_nwp")
   if (allocated ( Inversion_Strength_Nwp_Pix)) then
      data_dim2_dtype_r4 = Inversion_Strength_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("inversion_base_height_nwp")
   if (allocated ( Inversion_Base_Nwp_Pix)) then
      data_dim2_dtype_r4 = Inversion_Base_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
case("inversion_top_height_nwp")
   if (allocated ( Inversion_Top_Nwp_Pix)) then
      data_dim2_dtype_r4 = Inversion_Top_Nwp_Pix(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)
   end if
end select
