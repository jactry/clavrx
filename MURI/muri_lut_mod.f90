! $Id$


module muri_lut_mod
   use cx_sds_type_definitions_mod
   
   use cx_sds_io_mod, only: &
      cx_sds_finfo &
      , cx_sds_read
   
   implicit none 
   type muri_lut_type
      logical :: is_read = .false.
      real, allocatable :: sol(:)
      real, allocatable :: sat(:)
      real, allocatable :: azi(:)
      real, allocatable :: ws(:)
      real  :: aot_550nm (8)
      real, allocatable :: app_refl(:,:,:,:,:,:,:)
      real, allocatable :: aot_aer(:,:,:,:,:,:,:) 
      
      real :: aot_aer_fine (8,6,4)
      real :: aot_aer_coarse (8,6,5)
      real :: refl_fine (8,6,4)
      real :: refl_coarse (8,6,5)
      
      contains
      procedure ::read_lut => muri_lut_type__read_lut 
      procedure ::sub_table => muri_lut_type__sub_table 
 
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
     
      istatus = cx_sds_read ( trim(lut_file),'Solar_Zenith_Angles', temp_2d_real)
      allocate ( this %sol, source = temp_2d_real(1,:))
      
      istatus = cx_sds_read ( trim(lut_file),'View_Zenith_Angles',temp_2d_real )
      allocate ( this %sat, source = temp_2d_real(1,:))
      
      istatus = cx_sds_read ( trim(lut_file),'Relative_Azimuth_Angles', temp_2d_real)
      allocate ( this %azi, source = temp_2d_real(1,:))
      
      istatus = cx_sds_read ( trim(lut_file),'Wind_Speed', temp_2d_real)
      allocate ( this %ws, source = temp_2d_real(1,:))
      
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
   subroutine muri_lut_type__sub_table ( this, sol,sat,azi,ws )
      use aw_lib_array
      class (  muri_lut_type ) :: this
      real, intent(in) :: sol
      real, intent(in) :: sat
      real, intent(in) :: azi
      real, intent(in) :: ws 
      
      integer,parameter :: idp = selected_int_kind(13)
      integer,parameter :: sp = selected_real_kind(p=6,r=37)
      integer,parameter :: dp = selected_real_kind(p=15,r=307)
      real ,allocatable :: temp_4d(:,:,:,:)
      
      integer :: i,j,k
      
      
      allocate (temp_4d(size(this % sol), size(this % sat),size(this % azi),size(this % ws) ))
      do i = 1, 8
         do j = 1, 6
            do k = 1, 4
               temp_4d = this % app_refl(:,:,:,i,j,:,k)
               this % refl_fine (i,j,k) = interp4d (this % sol, this % sat, this % azi, this%ws,temp_4d,sol,sat,azi,ws )
                temp_4d = this % aot_aer(:,:,:,i,j,:,k)
               this % aot_aer_fine (i,j,k) = interp4d (this % sol, this % sat, this % azi, this%ws,temp_4d,sol,sat,azi,ws )
               
            end do
         end do
      end do 
      
      
       do i = 1, 8
         do j = 1, 6
            do k = 1, 5
               temp_4d = this % app_refl(:,:,:,i,j,:,k+4)               
               this % refl_coarse (i,j,k) = interp4d (this % sol, this % sat, this % azi, this%ws,temp_4d,sol,sat,azi,ws )
               temp_4d = this % aot_aer(:,:,:,i,j,:,k+4)
               this % aot_aer_coarse (i,j,k) = interp4d (this % sol, this % sat, this % azi, this%ws,temp_4d,sol,sat,azi,ws )
            end do
         end do
      end do   
        
      
   
   end subroutine muri_lut_type__sub_table
   
  
   



end module muri_lut_mod





