module cr_data_pool_mod
   
   use cx_sfc_mod, only: sfc_main_type
   use cx_geo_mod, only:lon_lat_index , geo_type
   use cx_nwp_mod, only:nwp_main_type
   use cx_sat_mod,only: sat_main_type
   use cx_imp_mod, only: imp_main_type
   use cx_rtm_mod, only: rtm_main_type
   
   type ( sfc_main_type ) :: sfc_obj_g
   type ( geo_type) :: geo_obj_g
   type ( nwp_main_type ) :: nwp_obj_g
   type ( sat_main_type ) :: sat_obj_g
   type ( imp_main_type ) :: imp_obj_g
   type ( rtm_main_type ) :: rtm_obj_g
   !type ( rtm_main_type ) :: rtm_obj_g
   !type ( sen_main_type ) :: sen_obj_g


end module cr_data_pool_mod
