! $Id$

module muri_forward_mod
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
      
      contains
      procedure :: read_lut
   end type muri_lut_type
   
   type(muri_lut_type) :: lut
   type muri_fwd_type
      real :: rfl(6)
      real :: jacobians
   end type muri_fwd_type
   
   
contains
   subroutine muri_forward ( state, fwd)
      real, intent(in) :: state
      type ( muri_fwd_type ), intent(out) :: fwd
   
      
      call lut%read_lut
      
      
      
      fwd % rfl(1) = 2.1
   
   end subroutine
   
   subroutine read_lut(this)
      class(muri_lut_type ) :: this
      
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
      !print*,sds_name
      istatus = cx_sds_read ( trim(lut_file),'Solar_Zenith_Angles',this % sol)
      istatus = cx_sds_read ( trim(lut_file),'View_Zenith_Angles', this % sat)
      istatus = cx_sds_read ( trim(lut_file),'Relative_Azimuth_Angles', this % azi)
      
      istatus = cx_sds_read( trim(lut_file), 'Apparent_Reflectance', this%app_refl)
      istatus = cx_sds_read( trim(lut_file), 'Apparent_Reflectance', this%app_refl)
      istatus = cx_sds_read ( trim(lut_file),'Aer_AOT_total', this % aot_aer )
      
      lut % aot_aer_fine =   this % aot_aer(1,1,1,:,:,1,1:4)
      lut % aot_aer_coarse = this % aot_aer(1,1,1,:,:,1,5:9)
      
      print*,lut% sol
      print*,'shape sol: ', shape(lut%sol)
      print*,' '
      print*,this % aot_aer(:,3,4,2,2,1,2)
      print*,'done'
      !print*,aer_tot
      
   end subroutine read_lut
   

end module  muri_forward_mod
