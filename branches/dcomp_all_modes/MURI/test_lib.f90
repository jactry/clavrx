program test_lib

use aw_lib_array


   integer,parameter :: idp = selected_int_kind(13)
   integer,parameter :: sp = selected_real_kind(p=6,r=37)
   integer,parameter :: dp = selected_real_kind(p=15,r=307)
   
  real(dp) :: w(3),x(3),y(3),z(3),array(3,3,3,3),w0,x0,y0,z0 ,fv
   real :: val
   integer :: ii
 
   
   
   w0=0.3
   x0=1.4
   y0 =2.3
   z0 = 7.9
   
   w=(/ 0.,1.,2./)
   x=(/ 0.,1.,2./)
   y=(/ 0.,1.,4./)
   z=(/ 6.,7.,8./)
   
   
   array(2,2,2,2) = 1.3
   array(1,2,2,2) = 3.3
   array(2,2,2,2) = 13.3
   
   do ii = 0,100 
   w0 = ii/10.
   fv=-999.
   val = interp4d(w,x,y,z,array,w0,x0,y0,z0,bounds_error=.false.,fill_value=fv)

   print*,w0,val
   end do
end program test_lib
