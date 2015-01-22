! collection of array tools
!
! 

module cx_tools_array_mod
   
   
   type array_stats_type
      logical :: is_allocated
      real, allocatable :: mean (:,:)
      real, allocatable :: min (:,:)
      real, allocatable :: max (:,:)
      real, allocatable :: std (:,:)  
   contains
      procedure :: allocate_it 
      procedure :: deallocate_it  
      procedure :: initialize 
   end type
   
   public:: array_statistics
   

contains
   subroutine allocate_it ( self, dim_1, dim_2)
      class ( array_stats_type) :: self
      integer, intent(in) :: dim_1, dim_2
      
      if ( self % is_allocated .eqv. .false.) then
         allocate ( self % mean ( dim_1, dim_2))
         allocate ( self % min ( dim_1, dim_2))
         allocate ( self % max ( dim_1, dim_2))
         allocate ( self % std ( dim_1, dim_2))
      end if
      self % is_allocated = .true.
   
   end subroutine allocate_it
   
   subroutine deallocate_it ( self)
      class ( array_stats_type) :: self
      
      
      if ( self % is_allocated ) then
         deallocate ( self % mean)
         deallocate ( self % min )
         deallocate ( self % max )
         deallocate ( self % std )
      end if
      self % is_allocated = .false.
   
   end subroutine deallocate_it
   
   subroutine initialize ( self )
      class ( array_stats_type ) :: self
      self % mean = 0.0
      self % std = 0.0
      self % min = 0.0
      self % max = 0.0
   
   end subroutine initialize
   
   ! ------------------------------------------------
   !
   ! -------------------------------------------------
   subroutine array_statistics ( &
       data_array & 
      , n_box &
      , out  &
      , mask_in  ) 
      ! this routine should provide mean , minvalue, maxvalue and standard deviation in a 
      !  nxn environment
   
      real , intent(in) :: data_array(:,:)
      integer , intent(in) :: n_box
      type(array_stats_type) :: out
      logical , intent(in) , optional :: mask_in(:,:)
           
      logical, allocatable :: mask ( :,:)
      real, dimension(n_box,n_box) :: box
      logical, dimension(n_box,n_box) :: mask_box
      integer :: dim_1 , dim_2
      integer  :: n_box_half 
      real :: sum_box
      integer :: count_box
      
      
      ! - executable
           
      n_box_half = (n_box -1)/ 2
      dim_1 = size(data_array,1)
      dim_2 = size(data_array,2)
      
      call out % allocate_it ( dim_1, dim_2)
      call out % initialize()
      allocate ( mask ( dim_1, dim_2) )
      mask = .true.
          
      if ( present ( mask_in) ) mask = mask_in
      
      line_loop: do j = 1 , dim_1 
         elem_loop: do i = 1, dim_2 
            if ( mask ( j , i) .eqv. .false. ) cycle elem_loop
            
            ! edge requires this
            i0 = max ( i - n_box_half , 1)
            i1 = min ( i + n_box_half , dim_2)
            j0 = max ( j - n_box_half , 1)
            j1 = min ( j + n_box_half , dim_1)
            
            ! -init local box
            mask_box = .false.
            box = -999.
            
            mask_box( 1 : 1+j1-j0, 1: 1+i1-i0) = mask( j0 : j1 , i0 : i1 )
            box( 1 : 1+j1-j0, 1: 1+i1-i0) = data_array ( j0 : j1 , i0 : i1 )
            sum_box = sum (  box , mask_box )
            count_box = max ( 1, count ( mask_box ) )  
            
            out % mean ( j , i ) = sum_box /  count_box
            out % max ( j , i ) = maxval (  box , mask_box ) 
            out % min ( j , i ) = minval (  box , mask_box ) 
            out % std ( j , i ) = sqrt ( sum ( box  ** 2 ,mask ) &
               & / count_box  -  out % mean ( j , i ) ** 2 )         
            
         end do elem_loop
      end do line_loop   
      
      deallocate ( mask)
   
   end subroutine array_statistics
   
   !  ------------------------------------------------
   !
   ! -------------------------------------------------
   subroutine local_radiative_center ( &
               data_array &
               , min_val &
               , max_val &   
               , dum &
               , mask_in  ) 

      real, intent(in) :: data_array(:,:)
      real, intent(in) :: min_val , max_val
      logical, optional , intent (in) :: mask_in (:,:)
      real,  intent(out)  :: dum(:,:)
      
      integer, parameter :: steps_max = 10
   
      integer :: dim_1 , dim_2
      real, dimension(3,3):: local_box , local_box_grad
      integer :: idx_1 , idx_2
      integer :: idx_1_step, idx_2_step
      integer :: start_1 , start_2
      integer :: end_1 , end_2
      integer :: high_grad(2)
      
      ! -executable
      
    
      dim_1 = size(data_array,1)
      dim_2 = size(data_array,2)
      
      dum = -999. 
      start_1 = 1
      start_2 = 1
      end_1 = dim_1
      end_2 = dim_2
   
      
      loop_1: do idx_1 = start_1 + 1 , end_1 - 1
         loop_2: do idx_2 = start_1 + 1 , end_2 - 1
      
            idx_1_step = idx_1
            idx_2_step = idx_2
            
            ! - loop over steps along gradient
            loop_step: do idx_step = 1 ,  steps_max  
               
               local_box = data_array ( idx_1_step - 1 : idx_1_step + 1 &
                                   & ,  idx_2_step - 1 : idx_2_step + 1  )   
               
               
               ! - if lowest value cycle
               if (minval(local_box) == data_array (idx_1_step,idx_2_step ) ) then
                  dum(idx_1,idx_2) = data_array ( idx_1_step, idx_2_step)
                  exit loop_step
               end if
               
               
               ! - find the next  pixel
               local_box_grad =  local_box - data_array (idx_1_step,idx_2_step )               
               high_grad =  minloc ( local_box_grad )
               
               idx_1_step = idx_1_step + high_grad(1) - 2
               idx_2_step = idx_2_step + high_grad(2) - 2
               
               ! - is this pixela has a lower value tham the min value
               if ( data_array ( idx_1_step, idx_2_step ) < min_val ) then
                  dum(idx_1,idx_2) = data_array ( idx_1_step, idx_2_step)
                  exit loop_step
               end if 
              
               ! - set to this value if edge value
               if (idx_1_step == 1 .or. idx_2_step == 1 &
                     .or. idx_1_step == end_1 .or. idx_2_step == end_2) then
                  dum(idx_1,idx_2) = data_array ( idx_1_step, idx_2_step)
                  exit loop_step
               end if
               
            end do loop_step
         end do loop_2
      end do loop_1   
   
  
   end subroutine local_radiative_center



end module cx_tools_array_mod


