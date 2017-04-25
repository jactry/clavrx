! $Id$
!
!
module muri_interface_mod

   type muri_in_array_type
      logical :: is_allocated = .false.
      integer :: dim(2)
      real,allocatable :: sol(:,:)
      real,allocatable :: sat(:,:)
      real,allocatable :: azi(:,:)
      real,allocatable :: ref(:,:,:)
      real,allocatable :: windspeed(:,:)
      logical, allocatable :: do_it(:,:)
      
      contains
      procedure :: allocate=>muri_in_array_type__allocate
      procedure :: deallocate=>muri_in_array_type__deallocate
   
   end type muri_in_array_type
   
   type muri_out_array_type
      logical :: is_allocated = .false.
      integer :: dim(2)
      real, allocatable :: aot(:,:)
      real, allocatable :: aot_channel(:,:,:)
      real, allocatable :: fmf(:,:)
      integer, allocatable :: fm_mode(:,:)
      integer, allocatable :: cm_mode(:,:)
      
      
      contains
      procedure :: allocate=>muri_out_array_type__allocate
      procedure :: deallocate =>muri_out_array_type__deallocate
   
   end type muri_out_array_type

contains

   subroutine muri_in_array_type__allocate(this,dim1,dim2)
      class(muri_in_array_type) :: this
      integer :: dim1,dim2
     
      
      ! - first check if already allocated with correct dims
      
      if ( this % is_allocated ) then 
         if ( dim1 .eq. this % dim(1) .and.  dim2 .eq. this % dim(2) ) return
      end if
      call this % deallocate()
      this % dim(1) = dim1
      this % dim(2) = dim2
      allocate ( this % sol( dim1,dim2))
      allocate ( this % sat( dim1,dim2))
      allocate ( this % azi( dim1,dim2))
      allocate ( this % windspeed( dim1,dim2))
      allocate ( this % do_it(dim1,dim2))
      
     
      allocate ( this % ref( 6,dim1,dim2))
     this % is_allocated = .true.
      
   
   end subroutine  muri_in_array_type__allocate  
   
   
   subroutine muri_in_array_type__deallocate(this)
      class(muri_in_array_type) :: this
      this % dim = [0,0]
     
      if (allocated (this % sol) ) deallocate ( this % sol)
      if (allocated (this % sat) ) deallocate ( this % sat)
      if (allocated (this % azi) ) deallocate ( this % azi)
      if (allocated (this % windspeed) ) deallocate ( this % windspeed)
      if (allocated (this % do_it) ) deallocate ( this % do_it)
      if (allocated (this % ref) ) deallocate ( this % ref)
       
       this % is_allocated = .false.
   
   end subroutine  muri_in_array_type__deallocate 
   
   
   
   !
   !
   !
   subroutine muri_out_array_type__allocate(this,dim1,dim2)
      class(muri_out_array_type) :: this
      integer, intent(in) :: dim1,dim2
      print*,'aou is allocated : ',this % is_allocated
      if ( this % is_allocated ) then
          if ( dim1 .eq. this % dim(1) .and.  dim2 .eq. this % dim(2) ) return      
      end if 
      print*,'out allocate again'
      call this % deallocate()
      
      this % dim(1) = dim1
      this % dim(2) = dim2
      allocate ( this % aot( dim1,dim2))
      allocate ( this % aot_channel(6, dim1,dim2))
      allocate ( this % fmf( dim1,dim2))
      allocate ( this % fm_mode( dim1,dim2))
      allocate ( this % cm_mode( dim1,dim2))
      this % aot = -999.
      this % is_allocated = .true.
   
   end  subroutine muri_out_array_type__allocate
   
   !
   !
   !
      subroutine muri_out_array_type__deallocate(this)
      class(muri_out_array_type) :: this
      integer :: test
      
      this % dim = [0,0]
      
      print*,allocated (this% aot)
      
      if (allocated (this% aot) ) deallocate ( this % aot, STAT =  test)
      print*,test
      print*,allocated(this% aot_channel)
      if (allocated (this% aot_channel) ) deallocate ( this % aot_channel)
      print*,allocated (this% fmf)
      if (allocated (this% fmf) ) deallocate ( this % fmf)
      print*,allocated (this% fm_mode)
      if (allocated (this% fm_mode) ) deallocate ( this % fm_mode)
      print*,allocated (this% cm_mode)
      if (allocated (this% cm_mode) ) deallocate ( this % cm_mode)
      
      this % is_allocated = .false.

   end subroutine  muri_out_array_type__deallocate 

end module muri_interface_mod
