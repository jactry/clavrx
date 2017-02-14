! $Id$
!
!

!program test_it

!interface
!   subroutine view2d(x)
!   real,dimension(:,:) :: x
!   end subroutine view2d
!end interface

!real::lon(100,100)=6.
!real::lat(100,100)=9.
!real , dimension(100,100) ::dat 
!dat(50,:)=88
!dat(60,:)=2228.



!call view2d(dat)



!end program test_it

subroutine aw_image_idl ( dat , lat,lon, minv ,maxv,titlev)
   real, intent(in), dimension( : ,:)  :: dat
   real, intent(in), dimension( : ,:)  :: lat
   real, intent(in), dimension( : ,:)  :: lon
   real, intent(in), optional :: minv ,maxv
   character ( len =*) , optional :: titlev
  integer :: x1,x2
  character(5) :: s1,s2
  character(len=200) :: idl_call
  character(len=20) :: title
  real :: minvalue , maxvalue
  
  !
  title='Clavrx Diagnostics'
  minvalue = minval(dat)
  maxvalue = maxval(dat)
  
  if ( present(minv)) minvalue = minv
  if ( present(maxv)) maxvalue = maxv
  if ( present ( titlev)) title = titlev
   call system ( 'rm -f aw_image_temp.dat')
   open ( 13, file ='aw_image_temp.dat', action='write')
  
         
            write ( 13 , * ) dat
     
   close (13)
      call system ( 'rm -f aw_image_temp_lat.dat')
   open ( 13, file ='aw_image_temp_lat.dat', action='write')
  
         
            write ( 13 , * ) lat
     
   close (13)
      call system ( 'rm -f aw_image_temp_lon.dat')
   open ( 13, file ='aw_image_temp_lon.dat', action='write')
  
         
            write ( 13 , * ) lon
     
   close (13)
   
   
   open(14,  file='aw_image_control.dat',action='write')
      write(14,*) minvalue
      write(14,*) maxvalue
      write(14,*) title
   close(14)
  
   x1=ubound (dat,1)
   x2=ubound (dat,2)

   write(s1, '(I5)' ) x1
    write(s2, '(I5)' ) x2
print*,s1,s2   
   idl_call = 'idl -quiet -e "aw_image_fortran" -args '//trim(s1)//' '//trim(s2)//'&' 
   call system(trim(idl_call)  ) 
   

end subroutine aw_image_idl



