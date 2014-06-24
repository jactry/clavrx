! $Id: dcomp_math_tools_mod.f90 385 2014-06-05 20:54:18Z awalther $
!
!  HISTORY: 06/05/2014: change of filename 
!
module dcomp_math_tools_mod

   integer , public :: debug_mode
contains

   !Subroutine to find the inverse of a square matrix
   !Author : Louisda16th a.k.a Ashwith J. Rego
   !Reference : Algorithm has been well explained in:
   !http://math.uww.edu/~mcfarlat/inverse.htm           
   !http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
   SUBROUTINE findinv(matrix, inverse, n, errorflag)
	   IMPLICIT NONE
	   !Declarations
	   INTEGER, INTENT(IN) :: n
	   INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
	   REAL, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
	   REAL, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
	
	   LOGICAL :: FLAG = .TRUE.
	   INTEGER :: i, j, k
	   REAL :: m
	   REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
	   !Augment input matrix with an identity matrix
	   DO i = 1, n
		   DO j = 1, 2*n
			   IF (j <= n ) THEN
				   augmatrix(i,j) = matrix(i,j)
			   ELSE IF ((i+n) == j) THEN
				   augmatrix(i,j) = 1
			   Else
				   augmatrix(i,j) = 0
			   END IF
		   END DO
	   END DO
	
	   !Reduce augmented matrix to upper traingular form
	   DO k =1, n-1
		   IF (augmatrix(k,k) == 0) THEN
			   FLAG = .FALSE.
			   DO i = k+1, n
				   IF (augmatrix(i,k) /= 0) THEN
					   DO j = 1,2*n
						   augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					   END DO
					   FLAG = .TRUE.
					   EXIT
				   END IF
				   IF (FLAG .EQV. .FALSE.) THEN
					   !PRINT*, "Matrix is non - invertible"
					   inverse = 0
					   errorflag = -1
					   return
				   END IF
			   END DO
		   END IF
		   DO j = k+1, n			
			   m = augmatrix(j,k)/augmatrix(k,k)
			   DO i = k, 2*n
				   augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			   END DO
		   END DO
	   END DO
	
	   !Test for invertibility
	   DO i = 1, n
		   IF (augmatrix(i,i) == 0) THEN
			   !PRINT*, "Matrix is non - invertible"
			   inverse = 0
			   errorflag = -1
			   return
		   END IF
	   END DO
	
	   !Make diagonal elements as 1
	   DO i = 1 , n
		   m = augmatrix(i,i)
		   DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		   END DO
	   END DO
	
	   !Reduced right side half of augmented matrix to identity matrix
	   DO k = n-1, 1, -1
		   DO i =1, k
		      m = augmatrix(i,k+1)
			   DO j = k, (2*n)
				   augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			   END DO
		   END DO
	   END DO				
	
	   !store answer
	   DO i =1, n
		   DO j = 1, n
			   inverse(i,j) = augmatrix(i,j+n)
		   END DO
	   END DO
	   errorflag = 0
   END SUBROUTINE findinv 

   ! ---------------------------

   subroutine dcomp_interpolation_weight ( count, value, data_in, weight_out, index_out , near_index)
  
      integer, intent ( in ) :: count
      real , intent ( in ) :: value
      real , dimension(:), intent( in ) :: data_in
	  
      real , intent( out ) , optional :: weight_out
      integer, intent ( out ) , optional :: index_out 
	   integer, intent ( out ) , optional :: near_index 
      integer :: index2
      real :: weight
      integer :: index
  
      call locate (data_in, count, value, index)
       
      index  = max ( 1 , min ( count - 1 , index ) )
      index2 = max ( 1 , min ( count , index + 1 ) )
      
      weight = (value - data_in( index )) / ( data_in( index2 ) - data_in(index))
	 
      if ( weight < 0. ) then
         index = 1 
		 weight = 0.
      end if
	  
	   if ( present (near_index ) ) near_index = nint ( index + weight )
      if ( present ( weight_out)) weight_out = weight
      if ( present ( index_out)) index_out = index
  
   end subroutine dcomp_interpolation_weight

   !-------------------------------------------------------------------------
   ! subroutine LOCATE(xx, n, x, j)
   ! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
   !--------------------------------------------------------------------------
   subroutine LOCATE(xx, n, x, j)

!   Arguments
      integer,                        intent(in)  :: n
      integer,                        intent(out) :: j
      real ,               intent(in)  :: x
      real , dimension(:), intent(in)  :: xx

!   Local variables
      integer :: i, jl, jm, ju

      jl = 0
      ju = n + 1
      do i = 1, 2*n
         if (ju-jl <= 1) then
            exit
         end if
         jm = (ju + jl) / 2
         if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
            jl = jm
         else
            ju = jm
         end if
      end do
      if (x == xx(1)) then
         j=1
      else if (x == xx(n)) then
         j = n - 1
      else
         j = jl
      endif

   end subroutine LOCATE
   
   
   ! ---------------------------------------------------------------------------------
   !
   !  INTERPOLATE_2D
   !
   !  Linear Interpolation for a 2 x 2 
   ! 
   !  Returns interpolated value of a 2D array with 2 elements for each dimension
   !
   !  INPUT: 
   !     table:      3D array
   !     Wgt_Dim1,Wgt_Dim2 : weights for each dimension
   !
   ! ----------------------------------------------------------------------------------
   subroutine interpolate_2d(table,Wgt_Dim1, Wgt_Dim2 , delta_1 , delta_2, value , dval_d1, dval_d2 ) 
      real, dimension(2,2), intent(in) :: table
      real, intent(in) :: Wgt_Dim1 , Wgt_Dim2
	   real, intent(in) :: delta_1 , delta_2
      
	   real , intent ( out )  :: value
	   real, intent ( out ) :: dval_d1, dval_d2
     

      !- locals
     
      real :: r , s
	  
	   r = wgt_dim1
	   s = wgt_dim2 


      value  =     ( 1.0 - r ) * ( 1.0 - s ) * table( 1 , 1 ) &
               & + r * ( 1.0 - s )          * table( 2 , 1 ) &
               & + (1.0 - r ) * s           * table( 1 , 2 ) &
				   & + r * s                    * table( 2 , 2 )

      dval_d2 =  ( ( 1.0 - r ) * ( table( 1 , 2 ) - table( 1 , 1 ) )   &
               & +  r * ( table( 2 , 2 ) - table( 2 , 1 ) ) ) /  &
               &  delta_2

      dval_d1 =   ( ( 1.0 - s ) * ( table( 2 , 1 ) - table( 1 , 1 ) )  &
                 &    + s * ( table( 2 , 2 ) - table( 1 , 2 ) ) ) /  &
                 &       delta_1 


                  
   end subroutine interpolate_2d  


end module dcomp_math_tools_mod
