! $Id:$
!
!   this tool reads HDF LUTs
!    
!
!    3 subroutines
!        1. read_hdf_dcomp_dims
!        2. read_hdf_dcomp_data_rfl
!        3. read_hdf_dcomp_data_ems

module dcomp_lut_hdf_mod
   
   use cx_sds_io_mod,only: &
      cx_sds_finfo &
      , cx_sds_read &
      , cx_sds_read_raw
   
   use cx_sds_type_definitions_mod, only: &
      cx_sds_type &
      , cx_att_type &
      , cx_sds_data_type &
      , MAXNCNAM
   

   implicit none

contains
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------   
   subroutine read_hdf_dcomp_dims ( file &
      , sat_zen &
      , sol_zen &
      , rel_azi &
      , cod &
      , cps )
      
      character ( len =*), intent(in) :: file
      
      real , allocatable, intent(out) :: sat_zen ( :)
      real , allocatable,  intent(out) :: sol_zen ( :)
      real , allocatable, intent(out) :: rel_azi ( :)
      real , allocatable,  intent(out) :: cod ( :)
      real , allocatable, intent(out) :: cps ( :) 
      
      if ( cx_sds_read ( file, 'sensor_zenith_angle', sat_zen) < 0 ) goto 9999
      if ( cx_sds_read ( file, 'solar_zenith_angle', sol_zen) < 0 ) goto 9999
      if ( cx_sds_read ( file, 'relative_azimuth_angle', rel_azi) < 0 ) goto 9999
      if ( cx_sds_read ( file, 'log10_optical_depth', cod) < 0 ) goto 9999
      if ( cx_sds_read ( file, 'log10_eff_radius', cps) < 0 ) goto 9999
      

9999 continue
   end    subroutine read_hdf_dcomp_dims
   
   ! ----------------------------------------------------------------
   !
   !  Reads reflectance tables
   !    input : hdf_file
   !    output: cloud_albedo
   !            cloud_transmission
   !            cloud_spherical_albedo
   !            cloud_reflectance
   !
   ! ----------------------------------------------------------------
   subroutine read_hdf_dcomp_data_rfl (file, cld_alb, cld_trn, cld_sph_alb, cld_refl )
      
      character ( len =*) :: file
      
      real , allocatable,  intent(out) :: cld_alb ( :,:,:)
      real , allocatable,  intent(out) :: cld_trn ( :,:,:)
      real , allocatable, intent(out) :: cld_sph_alb ( :,:)
      real ,  allocatable, intent(out) :: cld_refl ( :,:,:,:,:)
      
    
      if ( cx_sds_read ( file, 'albedo', cld_alb) < 0 ) goto 9999
      if ( cx_sds_read ( file, 'transmission', cld_trn) < 0 ) goto 9999
      if ( cx_sds_read ( file, 'spherical_albedo', cld_sph_alb) < 0 ) goto 9999
      if ( cx_sds_read ( file, 'reflectance', cld_refl) < 0 ) goto 9999
     
9999 continue   
   
   end subroutine read_hdf_dcomp_data_rfl

   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------      
   subroutine read_hdf_dcomp_data_ems (hdf_file, cld_ems, cld_trn_ems)
      character ( len =*) , intent(in) :: hdf_file
      
      real , allocatable, intent(out) :: cld_ems ( :,:,:)
      real , allocatable, intent(out) :: cld_trn_ems ( :,:,:)
      
      integer :: nsds
      
      
      
      integer , parameter :: N_PARAMS_EMS = 2
      
      
      if ( cx_sds_read ( hdf_file, 'cloud_emissivity', cld_ems) < 0 ) goto 9999
      if ( cx_sds_read ( hdf_file, 'cloud_transmission', cld_trn_ems) < 0 ) goto 9999
      
      
      
     
      ! cld_ems = reshape(psd%r4values,[9,29,45])
         
9999 continue     
      ! cld_trn_ems = reshape(psd%r4values,[9,29,45])
   end subroutine read_hdf_dcomp_data_ems

end module dcomp_lut_hdf_mod
