! $Id$


module muri_lut_mod
   use cx_sds_type_definitions_mod
   
   use cx_sds_io_mod, only: &
      cx_sds_finfo &
      , cx_sds_read
   
   implicit none 
   type muri_lut_type
      logical :: is_read = .false.
      real, allocatable :: sol(:,:)
      real, allocatable :: sat(:,:)
      real, allocatable :: azi(:,:)
      real, allocatable :: ws(:,:)
      real  :: aot_550nm (8)
      real, allocatable :: app_refl(:,:,:,:,:,:,:)
      real, allocatable :: aot_aer(:,:,:,:,:,:,:) 
      
      real :: aot_aer_fine (8,6,4)
      real :: aot_aer_coarse (8,6,5)
      real :: refl_fine (8,6,4)
      real :: refl_coarse (8,6,5)
      
      contains
      procedure ::read_lut => muri_lut_type__read_lut 
 
   end type muri_lut_type
   
   type ( muri_lut_type)  :: lut
   
   
contains
  !
   !
   !
   subroutine muri_lut_type__read_lut(this)
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
      logical :: file_exists
      real, allocatable :: temp_2d_real(:,:)
      
      
      if ( this % is_read) return
      
      lut_file = trim('/DATA/AHI_AEROSOL/AHI_Aerosol_LUT/AHI_Ocean_Aerosol_LUT_v1.hdf')
      INQUIRE(file = lut_file,EXIST=file_exists)
      if ( .not. file_exists) then 
         print*,'LUT file not there stopping'
         stop
      
      end if
      istatus = cx_sds_finfo ( trim(lut_file), ftype, nsds,sds_name, natt,att_name)
     
      istatus = cx_sds_read ( trim(lut_file),'Solar_Zenith_Angles',this % sol)
      istatus = cx_sds_read ( trim(lut_file),'View_Zenith_Angles', this % sat)
      istatus = cx_sds_read ( trim(lut_file),'Relative_Azimuth_Angles', this % azi)
      istatus = cx_sds_read ( trim(lut_file),'Wind_Speed', this % ws)
      istatus = cx_sds_read ( trim(lut_file),'AOT_at_550nm',temp_2d_real)
      this % aot_550nm(:) = temp_2d_real (1,:)
      istatus = cx_sds_read( trim(lut_file), 'Apparent_Reflectance', this%app_refl)
     
      istatus = cx_sds_read ( trim(lut_file),'Aer_AOT_total', this % aot_aer )
      
      
      ! - this is okay LUTs replicate results for some dimensions
      this % aot_aer_fine =   this % aot_aer(1,1,1,:,:,1,1:4)
      this % aot_aer_coarse = this % aot_aer(1,1,1,:,:,1,5:9)
      
      !TODO this has to be interpolated
      this % refl_fine =   this % app_refl(1,1,1,:,:,1,1:4)
      this % refl_coarse = this % app_refl(1,1,1,:,:,1,5:9)
      this % is_read = .true.
      
   end subroutine muri_lut_type__read_lut
   !
   !
   !
  ! subroutine muri_lut_type__sub_table ( this, sol,sat,azi,ws )
      
   
   
 !  end subroutine muri_lut_type__sub_table
   
  
   



end module muri_lut_mod





