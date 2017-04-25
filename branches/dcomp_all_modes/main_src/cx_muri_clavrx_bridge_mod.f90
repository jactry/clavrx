!
!
!            mod    ahi
! 0.47       1      3   
! 0.51       2      4
! 0.64       3      1  
! 0.86       4      2  
! 1.6        5      6
! 2.2        6      7
!
!
module cx_muri_clavrx_bridge_mod

    use muri_interface_mod , only: &
         muri_in_array_type &
       , muri_out_array_type
   use pixel_common, only: &
      geo,ch,cldmask,sfc
   
   type cx_muri_structure
      logical :: is_set = .false.
      real,allocatable :: aod(:,:)
      integer :: dim(2)
        contains
        procedure :: allocate =>  cx_muri_structure__allocate
        procedure :: deallocate =>  cx_muri_structure__deallocate
   end type cx_muri_structure
   
   type(cx_muri_structure),save :: muri
   
   
contains

   subroutine cx_muri_algorithm(dim1,dim2) 
      implicit none
      integer, intent(in) :: dim1,dim2
      type(muri_in_array_type),save :: input
      type(muri_out_array_type),save :: output
      integer :: i 
      integer, parameter :: ahi_map_modis(6) = [3,4,1,2,6,7]
     
     
      
      ! it is more efficient to allocate only for first segment
      
     
      ! - this checks if allocation is needed 
      call input % allocate ( dim1,dim2)
    
      input % sol = geo % solzen
      input % sat = geo % satzen
      input % azi = geo % relaz
      input % windspeed = 2.0
     
      do i=1,6 
         
         input % ref(i,:,:) = ch(ahi_map_modis(i))%ref_toa/100.
        
      end do
      
      input % do_it = CLDMASK%cld_mask(1:dim1,1:dim2) == 0 .AND. sfc%land(1:dim1,1:dim2) .NE. 1
      
      ! - this checks if allocation is needed 
      call output % allocate(dim1,dim2)
      call  muri_array_loop (input, output )
      
      
      !call muri % allocate (dim1,dim2)
      muri % aod = output % aot
      !call output % deallocate
      ! call input % deallocate
      
      
   end subroutine cx_muri_algorithm
   
   
   
   subroutine cx_muri_structure__allocate(this, dim1, dim2)
      class(cx_muri_structure) :: this
      integer, intent(in) :: dim1,dim2
      
      if ( this % is_set ) then 
         if (dim1 .eq. this % dim(1) .and. dim2 .eq. this % dim(2) ) return 
      end if
      
      call this % deallocate
      
      allocate ( this % aod ( dim1,dim2))
      this % is_set = .true.
   end subroutine cx_muri_structure__allocate
   
   subroutine cx_muri_structure__deallocate(this)
      class(cx_muri_structure) :: this
      if ( allocated (this % aod) ) deallocate ( this % aod)
      this % is_set = .false.
   end subroutine cx_muri_structure__deallocate

end module cx_muri_clavrx_bridge_mod
