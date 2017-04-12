!
!
!
!
!

module aw_lib_array
   use lib_array,only: &
      locate &
      , interp2d
      
   implicit none
   save
   private
   
   integer,parameter :: idp = selected_int_kind(13)
   integer,parameter :: sp = selected_real_kind(p=6,r=37)
   integer,parameter :: dp = selected_real_kind(p=15,r=307)


   public :: interp4d
   interface interp4d
      module procedure interp4d_sp
      module procedure interp4d_dp
     
   end interface interp4d

contains
 
  
    function interp4d_dp(w,x,y,z,array,w0,x0,y0,z0,bounds_error,fill_value) result(value)
    ! Bilinar interpolation of array = f(w,x,y,z) at (w0,x0,y0,z0)

    implicit none

    real(dp),intent(in) :: w(:),x(:),y(:),z(:),array(:,:,:,:),w0,x0,y0,z0

    logical,intent(in),optional :: bounds_error
    ! whether to raise an out of bounds error

    real(dp),intent(in),optional :: fill_value
    ! value for out of bounds if bounds_error is .false.

    real(dp) :: value,norm
    integer :: i1,i2,j1,j2
    integer :: k1,k2,m1,m2
    
    logical :: bounds_error_tmp
    real(dp) :: fill_value_tmp
    
    integer :: ii,jj
    real(dp),allocatable :: temp_array_2d_interp(:,:)
    real(dp),allocatable :: temp_2d(:,:)
    
    if(present(bounds_error)) then
       bounds_error_tmp = bounds_error
    else
       bounds_error_tmp = .true.
    end if

    if(.not.bounds_error_tmp) then
       if(present(fill_value)) then
          fill_value_tmp = fill_value
       else
          fill_value_tmp = 0._dp
       end if
    end if

    if(size(w).ne.size(array,1)) stop "w does not match array"
    if(size(x).ne.size(array,2)) stop "x does not match array"
    if(size(y).ne.size(array,3)) stop "y does not match array"
    if(size(z).ne.size(array,4)) stop "z does not match array"  
   
    i1 = locate(w,w0) ; i2 = i1 + 1
    j1 = locate(x,x0) ; j2 = j1 + 1
    k1 = locate(y,y0) ; k2 = k1 + 1
    m1 = locate(z,z0) ; m2 = m1 + 1

    if(i1==-1) then
       if(bounds_error_tmp) then
          write(0,'("ERROR: Interpolation out of bounds : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') w0,w(1),w(size(w))
          stop
       else
          value = fill_value_tmp
          return
       end if
    end if

 if(j1==-1) then
       if(bounds_error_tmp) then
          write(0,'("ERROR: Interpolation out of bounds : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') x0,x(1),x(size(x))
          stop
       else
          value = fill_value_tmp
          return
       end if
    end if
    
    
       if(k1==-1) then
       if(bounds_error_tmp) then
          write(0,'("ERROR: Interpolation out of bounds : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') y0,y(1),y(size(y))
          stop
       else
          value = fill_value_tmp
          return
       end if
    end if

 if(m1==-1) then
       if(bounds_error_tmp) then
          write(0,'("ERROR: Interpolation out of bounds : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') z0,z(1),z(size(z))
          stop
       else
          value = fill_value_tmp
          return
       end if
    end if
    
   
    ! first reduce from 4d to 2d
    allocate ( temp_array_2d_interp(size(w),size(x))  )
    allocate ( temp_2d(size(y),size(z))  )
    do ii = 1, size(w)
      do jj = 1, size(x)
      
         temp_2d = array(ii,jj,:,:)
        
         temp_array_2d_interp(ii,jj) = interp2d (y,z,temp_2d , y0,z0)!,bounds_error=.true.,fill_value = -999.)
        
      end do
    end do   
    
    value = interp2d (w,x,temp_array_2d_interp,w0,x0)!,bounds_error=.true.,fill_value = -999.)
   

  end function interp4d_dp
  
    function interp4d_sp(w,x,y,z,array,w0,x0,y0,z0,bounds_error,fill_value) result(value)
    ! Bilinar interpolation of array = f(w,x,y,z) at (w0,x0,y0,z0)

    implicit none

    real(sp),intent(in) :: w(:),x(:),y(:),z(:),array(:,:,:,:),w0,x0,y0,z0

    logical,intent(in),optional :: bounds_error
    ! whether to raise an out of bounds error

    real(dp),intent(in),optional :: fill_value
    ! value for out of bounds if bounds_error is .false.

    real(sp) :: value,norm
    integer :: i1,i2,j1,j2
    integer :: k1,k2,m1,m2
    
    logical :: bounds_error_tmp
    real(dp) :: fill_value_tmp
    
    integer :: ii,jj
    real(sp),allocatable :: temp_array_2d_interp(:,:)

    if(present(bounds_error)) then
       bounds_error_tmp = bounds_error
    else
       bounds_error_tmp = .true.
    end if

    if(.not.bounds_error_tmp) then
       if(present(fill_value)) then
          fill_value_tmp = fill_value
       else
          fill_value_tmp = 0._dp
       end if
    end if
   
    if(size(w).ne.size(array,1)) stop "w does not match array in 4d_interp"
    if(size(x).ne.size(array,2)) stop "x does not match array in 4d_interp"
    if(size(y).ne.size(array,3)) stop "y does not match array in 4d_interp"
    if(size(z).ne.size(array,4)) stop "z does not match array in 4d_interp"  
   
    i1 = locate(w,w0) ; i2 = i1 + 1
    j1 = locate(x,x0) ; j2 = j1 + 1
    k1 = locate(y,y0) ; k2 = k1 + 1
    m1 = locate(z,z0) ; m2 = m1 + 1

    if(i1==-1) then
       if(bounds_error_tmp) then
          write(0,'("ERROR: Interpolation out of bounds 1.DIM : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') w0,w(1),w(size(w))
          stop
       else
          value = fill_value_tmp
          return
       end if
    end if

 if(j1==-1) then
       if(bounds_error_tmp) then
          write(0,'("ERROR: Interpolation out of bounds 2.DIM: ",ES11.4," in [",ES11.4,":",ES11.4,"]")') x0,x(1),x(size(x))
          stop
       else
          value = fill_value_tmp
          return
       end if
    end if
    
    
       if(k1==-1) then
       if(bounds_error_tmp) then
          write(0,'("ERROR: Interpolation out of bounds 3.DIM: ",ES11.4," in [",ES11.4,":",ES11.4,"]")') y0,y(1),y(size(y))
          stop
       else
          value = fill_value_tmp
          return
       end if
    end if

 if(m1==-1) then
       if(bounds_error_tmp) then
          write(0,'("ERROR: Interpolation out of bounds 4.DIM: ",ES11.4," in [",ES11.4,":",ES11.4,"]")') z0,z(1),z(size(z))
          stop
       else
          value = fill_value_tmp
          return
       end if
    end if
    
   
    ! first reduce from 4d to 2d
    allocate ( temp_array_2d_interp(size(w),size(x))  )
    do ii = 1, size(w)
      do jj = 1, size(x)
      
         
         temp_array_2d_interp(ii,jj) = interp2d (y,z, array(ii,jj,:,:) , y0,z0,bounds_error=.true.,fill_value = -999.)
      
      end do
    end do   
    
    value = interp2d(w,x,temp_array_2d_interp,w0,x0) !,bounds_error=.true.,fill_value = -999.)
   
   
   
         
            

  end function interp4d_sp
  
 
   

end module aw_lib_array
