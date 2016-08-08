! $Header$

module cx_array_tools_mod

   implicit none
   
   interface cx_rebin
      module procedure cx_rebin_real, cx_rebin_i1
   end interface


contains

 ! thjos function redo the IDL rebin function
   ! this should be moved to a module asasp.
   function cx_rebin_real ( array, d1, d2 )
      real, intent(in) :: array(:,:)
      integer, intent(in) :: d1, d2
      real, allocatable :: cx_rebin_real(:,:)
      
      integer :: dim_array (2)
      integer :: i, j
      
      dim_array = shape ( array)
      allocate (cx_rebin_real(d1,d2))
      
      do i = 1 , dim_array(1)
         do j = 1, d2/2
            cx_rebin_real(2*i-1 : 2*i, 2*j-1 : 2*j) = array(i,j)
         
         end do
      end do   
   
   
   end function cx_rebin_real

 ! thjos function redo the IDL rebin function
   ! this should be moved to a module asasp.
   function cx_rebin_i1 ( array, d1, d2 )
      integer(kind=1), intent(in) :: array(:,:)
      integer, intent(in) :: d1, d2
      integer(kind=1), allocatable :: cx_rebin_i1(:,:)
      
      integer :: dim_array (2)
      integer :: i, j
      
      dim_array = shape ( array)
      allocate (cx_rebin_i1(d1,d2))
      
      do i = 1 , dim_array(1)
         do j = 1, d2/2
            cx_rebin_i1(2*i-1 : 2*i, 2*j-1 : 2*j) = array(i,j)
         
         end do
      end do   
   
   
   end function cx_rebin_i1
end module cx_array_tools_mod
