! $Id$
module muri_mod
   use cx_sds_io_mod, only: &
      cx_sds_finfo &
      , cx_sds_read
   
   use cx_sds_type_definitions_mod
   
   type muri_lut_type
      real, allocatable :: sol(:,:)
      real, allocatable :: sat(:,:)
      real, allocatable :: azi(:,:)
      real, allocatable :: app_refl(:,:,:,:,:,:,:)
      real, allocatable :: aot_aer(:,:,:,:,:,:,:) 
      real :: aot_aer_fine (8,6,4)
      real :: aot_aer_coarse (8,6,5)
   end type muri_lut_type
   
   type(muri_lut_type) :: lut
   
   
   type muri_input_type
      real :: rfl(6)
   
   end type muri_input_type
   
   type ( muri_input_type) :: inp
   
   !integer, parameter :: MAXNCDIM = 32
   !integer, parameter :: MAXNCNAM = 128

contains
   subroutine start
      real :: ref(6)
      real :: lat,lon,sol,sat,azi
      character (len =200) :: lut
      call read_lut
   
      print*,'klopp'
   end subroutine start
   
   
   subroutine read_lut
      character (len =400)::lut_file
      integer :: istatus
      integer :: ftype
      integer :: nsds
      integer :: natt
      character ( len = MAXNCNAM), allocatable :: sds_name(:)
      character ( len = MAXNCNAM), allocatable :: att_name(:)
      real,allocatable :: sol_zen_ang(:,:)
      real,allocatable :: sat_zen_ang(:,:)
      
      
      lut_file = trim('/DATA/AHI_AEROSOL/AHI_Aerosol_LUT/AHI_Ocean_Aerosol_LUT_v1.hdf')
      istatus = cx_sds_finfo ( trim(lut_file), ftype, nsds,sds_name, natt,att_name)
      print*,sds_name
      istatus = cx_sds_read ( trim(lut_file),'Solar_Zenith_Angles',lut % sol)
      istatus = cx_sds_read ( trim(lut_file),'View_Zenith_Angles', lut % sat)
      istatus = cx_sds_read ( trim(lut_file),'Relative_Azimuth_Angles', lut % azi)
      
      istatus = cx_sds_read( trim(lut_file), 'Apparent_Reflectance', lut%app_refl)
      istatus = cx_sds_read( trim(lut_file), 'Apparent_Reflectance', lut%app_refl)
      istatus = cx_sds_read ( trim(lut_file),'Aer_AOT_total', lut%aot_aer )
      
      lut % aot_aer_fine = lut%aot_aer(1,1,1,:,:,1,1:4)
      lut % aot_aer_coarse = lut%aot_aer(1,1,1,:,:,1,5:9)
      
      print*,lut% sol, shape(lut%sol)
      !print*,aer_tot
      
   end subroutine read_lut


end module muri_mod


