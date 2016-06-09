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

   use cx_hdf_read_mod , only : &
         & hdf_sds &
         , hdf_data &
         , hdf_get_file_sds &
         , MAXNCNAM
   implicit none

contains
   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------   
   subroutine read_hdf_dcomp_dims ( hdf_file &
      , sat_zen &
      , sol_zen &
      , rel_azi &
      , cod &
      , cps )
      
      character ( len =*), intent(in) :: hdf_file
      
      real , intent(out) :: sat_zen ( 45)
      real , intent(out) :: sol_zen ( 45)
      real , intent(out) :: rel_azi ( 45)
      real , intent(out) :: cod ( 29)
      real , intent(out) :: cps ( 9) 
      
      integer :: nsds
      character(len=MAXNCNAM), dimension(:), allocatable :: &
       sds_name

      type(hdf_sds), dimension(:), allocatable, target :: &
       sds                          ! Tableau des structures des SDS extraits

      integer :: isds
      integer  :: i                                      
       
      type(hdf_sds), pointer                           :: &
       ps                           ! Pointeur sur la structure du SDS courant 
      
      type(hdf_data), pointer                          :: &
       psd                        ! Pointeur sur les données du SDS courant
   
   
      allocate ( sds_name ( 5)) 
      sds_name =(/ character(len=30) :: 'sensor_zenith_angle'  &
            , 'solar_zenith_angle' &
            , 'relative_azimuth_angle' &
            , 'log10_optical_depth' &
            , 'log10_eff_radius'/)
      
       
      if (hdf_get_file_sds(hdf_file, nsds, sds, nsdsn = 5, sds_name = sds_name) < 0) then
         print*,'hdf file not readable ', trim(hdf_file)
         stop  
      end if 
        
      ps => sds(1); psd=> ps%data
      sat_zen = psd%r4values 
       
      ps => sds(2); psd=> ps%data
      sol_zen = psd%r4values 
      
      ps => sds(3); psd=> ps%data
       rel_azi = psd%r4values 
      
      ps => sds(4); psd=> ps%data
      cod = psd%r4values 
      
      ps => sds(5); psd=> ps%data      
      cps = psd%r4values 
      
      psd => null()
      
      deallocate ( sds_name)     

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
   subroutine read_hdf_dcomp_data_rfl (hdf_file, cld_alb, cld_trn, cld_sph_alb, cld_refl )
      
      character ( len =*) :: hdf_file
      
      real , intent(out) :: cld_alb ( 9,29,45)
      real , intent(out) :: cld_trn ( 9,29,45)
      real , intent(out) :: cld_sph_alb ( 9,29)
      real , intent(out) :: cld_refl ( 9,29,45,45,45)
      
      integer :: nsds
      character(len =MAXNCNAM), allocatable :: sds_name(:)
      type(hdf_sds) , allocatable, target :: sds(:)
      
      type(hdf_sds) , pointer :: ps => null()
      type(hdf_data), pointer :: psd => null()
      integer , parameter :: N_PARAMS = 4
      
      integer :: i , last , first 
   
      allocate ( sds_name ( N_PARAMS) )
      sds_name = (/ character (len =20) :: 'albedo' , 'transmission' , 'spherical_albedo', 'reflectance'  /)
      
      if ( hdf_get_file_sds ( hdf_file, nsds , sds , nsdsn = N_PARAMS, sds_name = sds_name ) < 0 ) then
         print*,'hdf file not readable'
         stop
      end if   
      deallocate ( sds_name )
   
      ps => sds(1); psd=> ps%data
      cld_alb = reshape(psd%r4values,[9,29,45])
      
      ps => sds(2); psd=>ps%data
      cld_trn =reshape (psd%r4values,[9,29,45])
      
      ps => sds(3); psd=> ps%data
      cld_sph_alb = reshape ( psd%r4values, [9,29] )
      
      ps => sds(4); psd=>ps%data
      
      ! reconstruct 5d array
      do i = 1 , 45 
         first = (i -1 ) *  9 * 29 * 45 * 45 + 1
         last = first + (9 * 29 * 45 * 45) - 1
         cld_refl(:,:,:,:,i) = reshape( psd%r4values(first : last  ),[9,29,45,45])
      end do 
       
   
   
   end subroutine read_hdf_dcomp_data_rfl

   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------      
   subroutine read_hdf_dcomp_data_ems (hdf_file, cld_ems, cld_trn_ems)
      character ( len =*) , intent(in) :: hdf_file
      
      real , intent(out) :: cld_ems ( 9,29,45)
      real , intent(out) :: cld_trn_ems ( 9,29,45)
      
      integer :: nsds
      
      type(hdf_sds) , allocatable, target :: sds(:)
      character(len =MAXNCNAM), allocatable :: sds_name_ems(:)
      type(hdf_sds) , pointer :: ps => null()
      type(hdf_data), pointer :: psd => null()
      
      integer , parameter :: N_PARAMS_EMS = 2
      
      
      allocate ( sds_name_ems ( N_PARAMS_EMS) )
      sds_name_ems = (/ character(len=20) :: 'cloud_emissivity' , 'cloud_transmission' /)
      
      if ( hdf_get_file_sds ( hdf_file, nsds , sds , nsdsn = N_PARAMS_EMS, sds_name = sds_name_ems ) < 0 ) then
         print*,'hdf file not readable'
         stop
      end if   
      deallocate ( sds_name_ems )
      
      ps => sds(1); psd=> ps%data
      cld_ems = reshape(psd%r4values,[9,29,45])
         
      ps => sds(2); psd=> ps%data
      cld_trn_ems = reshape(psd%r4values,[9,29,45])
   end subroutine read_hdf_dcomp_data_ems

end module dcomp_lut_hdf_mod
