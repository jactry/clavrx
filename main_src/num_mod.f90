! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: num_mod.f90 (src)
!       NUMERICAL_ROUTINES (program)
!
! PURPOSE: library of useful numerical functions
!
! Description: 
!              Note, routines that only appear in the volcanic
!              ash functions are in ash_num_mod.f90.  This was
!              done to minimize code delivery to NCDC for PATMOS-x.
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! public:: 
!   LOCATE
!   JULIAN
!   COMPUTE_MONTH
!   COMPUTE_DAY
!   VAPOR
!   VAPOR_ICE
!   INVERT_2x2
!   INVERT_3x3
!   INVERT_4x4
!   FIND_BOUNDS
!   PACK_BYTES
!   COMPUTE_TIME_HOURS
!   COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES
!   GRADIENT_MEANDER
!   COMPUTE_MEDIAN
!   COMPUTE_MEDIAN_SEGMENT
!
!--------------------------------------------------------------------------------------
 module NUMERICAL_ROUTINES
  use CONSTANTS
  
  implicit none
  private
  public:: LOCATE, &
          
           INVERT_MATRIX,  &
           INVERT_2x2,  &
           INVERT_3x3,  &
           INVERT_4x4,  &
           INVERT_DIAGONAL,  &
           FIND_BOUNDS,  &
           
           COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES,  &
           GRADIENT_MEANDER, &
           COMPUTE_MEDIAN, &
           COMPUTE_MEDIAN_SEGMENT, &
           
           COVARIANCE
           
  contains

!----------------------------------------------------------------------
! subroutine COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES
!
! Routine to compute spatial uniformity for every pixel in an array
! based on its nxn neighborhood
!
! z - the input array
! n - the size of the box (1=3x3, 2=5x5, ...)
! bad_mask - mask array (only values with sym%YES will contribute)
! imin - starting x-index of array
! imax - ending x-index of array
!
! z_mean - mean of z over the (2n+1 x 2n+1) box
! z_min - minimum of z over the (2n+1 x 2n+1) box
! z_max - maximum of z over the (2n+1 x 2n+1) box
! z_std - standard deviation of z over the (2n+1 x 2n+1) box
! i_loc_of_max - i index of the maximum value
!----------------------------------------------------------------------
subroutine COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                          z, &
                                          n, &
                                          uni_land_mask_flag, &
                                          bad_mask, &
                                          land_mask,  &
                                          imin, &
                                          imax, &
                                          jmin, &
                                          jmax, &
                                          z_mean, &
                                          z_min,  &
                                          z_max,  &
                                          z_std,  &
                                          i_loc_of_max,  &
                                          j_loc_of_max,  &
                                          i_loc_of_min,  &
                                          j_loc_of_min)

  real(kind=real4), dimension(:,:), intent(in):: z
  integer(kind=int1), dimension(:,:), intent(in):: bad_mask
  integer(kind=int1), dimension(:,:), intent(in):: land_mask
  integer, intent(in):: uni_land_mask_flag
  real(kind=real4), dimension(:,:), intent(out):: z_mean
  real(kind=real4), dimension(:,:), intent(out):: z_min
  real(kind=real4), dimension(:,:), intent(out):: z_max
  real(kind=real4), dimension(:,:), intent(out):: z_std
  integer, intent(in):: n
  integer, intent(in):: imin
  integer, intent(in):: imax
  integer, intent(in):: jmin
  integer, intent(in):: jmax
  integer, dimension(:,:), intent(out):: i_loc_of_min
  integer, dimension(:,:), intent(out):: j_loc_of_min
  integer, dimension(:,:), intent(out):: i_loc_of_max
  integer, dimension(:,:), intent(out):: j_loc_of_max
  integer:: i
  integer:: j
  integer:: i1
  integer:: i2
  integer:: j1
  integer:: j2
  integer:: ii
  integer:: jj
  integer:: n_good

  !--- initialize to missing
  z_mean = missing_value_real4
  z_min = missing_value_real4
  z_max = missing_value_real4
  z_std = missing_value_real4
  i_loc_of_max = missing_value_int1
  j_loc_of_max = missing_value_int1
  i_loc_of_min = missing_value_int1
  j_loc_of_min = missing_value_int1

  line_loop: do j = jmin, jmax - jmin + 1

  !--- set limits of NxN array in the j-direction
  j1 = max(jmin,j-n)   !top index of local array
  j2 = min(jmax,j+n)   !bottom index of local array

  element_loop: do i = imin, imax - imin + 1


     !--- initial checks for this pixel
     if (bad_mask(i,j) == sym%YES .or. z(i,j) == Missing_Value_Real4) then
         z_mean(i,j) = Missing_Value_Real4
         z_std(i,j) = Missing_Value_Real4
         z_min(i,j) = Missing_Value_Real4
         z_max(i,j) = Missing_Value_Real4
         cycle
      endif


     !--- set limits of NxN array in the i-direction
     i1 = max(imin,i-n)   !left index of local array
     i2 = min(imax,i+n)   !right index of local array

     !--- initialize
     z_mean(i,j) = 0.0
     z_std(i,j) = 0.0
     z_min(i,j) = huge(z_min(i,j))
     z_max(i,j) = -1.0*huge(z_min(i,j))
     n_good = 0

     !--- go through each element in NxN array
     sub_line_loop: do jj = j1,j2
      sub_element_loop: do ii = i1,i2

        if (bad_mask(ii,jj) == sym%YES) then
          cycle
        endif

        if (z(ii,jj) == missing_value_real4) then
          cycle
        endif
       
        if ((uni_land_mask_flag == sym%YES) .and. &
             (land_mask(i,j) /= land_mask(ii,jj))) then    
          cycle
        endif

        n_good = n_good + 1
        z_mean(i,j) = z_mean(i,j) + z(ii,jj)
        z_std(i,j) = z_std(i,j) + z(ii,jj)**2
 
        if (z(ii,jj) < z_min(i,j)) then
           z_min(i,j) = z(ii,jj)
           i_loc_of_min(i,j) = ii
           j_loc_of_min(i,j) = jj
        endif

        if (z(ii,jj) > z_max(i,j)) then
           z_max(i,j) = z(ii,jj)
           i_loc_of_max(i,j) = ii
           j_loc_of_max(i,j) = jj
        endif

       end do sub_element_loop
     end do sub_line_loop

     !--- if any good pixel found, compute mean and standard deviation
     if (n_good > 0) then 
       z_mean(i,j) = z_mean(i,j) / n_good
       z_std(i,j) = sqrt(max(0.0,(z_std(i,j)/n_good - z_mean(i,j)**2)))
     endif
      
 end do element_loop

end do line_loop

end subroutine COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES

!-------------------------------------------------------------------------
! subroutine LOCATE(xx, n, x, j)
! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
!--------------------------------------------------------------------------
  subroutine LOCATE(xx, n, x, j)

!   Arguments
    integer,                        intent(in)  :: n
    integer,                        intent(out) :: j
    real (kind=ipre),               intent(in)  :: x
    real (kind=ipre), dimension(:), intent(in)  :: xx

!   Local variables
    integer :: i, jl, jm, ju

    jl = 0
    ju = n + 1
    do i = 1, 2*n
       if (ju-jl <= 1) then
          exit
       endif
       jm = (ju + jl) / 2
       if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if (x == xx(1)) then
       j=1
    else if (x == xx(n)) then
       j = n - 1
    else
       j = jl
    endif

  end subroutine LOCATE



!--------------------------------------------------------------------------
! Invert a square matrix
!
! follows:
! http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
!--------------------------------------------------------------------------
function INVERT_MATRIX(Matrix, Matrix_Inv, Matrix_Size) RESULT(Status)

  real(kind=real4), dimension(:,:), intent(in) :: Matrix
  real(kind=real4), dimension(:,:), intent(OUT) :: Matrix_Inv
  integer(kind=INT4), intent(in) :: Matrix_Size

  integer(kind=INT4) :: Status
  integer(kind=INT4) :: Singular_Flag
  integer:: i,j,ni, nj
  real(kind=real4) :: Zero
  logical:: Diag_Flag

  Status = Sym%SUCCESS

  Zero = epsilon(Matrix(1,1))

  !---- check for a square matrix
  ni = size(Matrix,1)
  nj = size(Matrix,2)
  if (ni /= nj) then
     print *, "size mismatch in invert matrix"
     Status = sym%FAILURE
     return
  endif
  !---- check for a diagonal matrix
  Diag_Flag = .true.
  do i = 1, ni
    do j = 1, nj
       if (i == j) cycle
       if (abs(Matrix(i,j)) > Zero) then
         Diag_Flag = .false.
         exit
       endif
    enddo
  enddo
 
  if (Diag_Flag) then
      call INVERT_DIAGONAL(Matrix, Matrix_Inv, Singular_Flag)
      if (Singular_Flag == 1) THEN
        print *, "diagonal failure"
        Status = Sym%FAILURE
      endif

  else

  !---- if not diagonal matrix, do full inversion
  select case(Matrix_Size)

    case(2)

      call INVERT_2x2(Matrix, Matrix_Inv, Singular_Flag)
      if (Singular_Flag == 1) THEN
        print *, "non-diagonal 2x2 failure"
        Status = Sym%FAILURE
      endif

    case(3)

      call INVERT_3x3(Matrix, Matrix_Inv, Singular_Flag)
      if (Singular_Flag == 1) THEN
        print *, "non-diagonal 3x3 failure"
        Status = Sym%FAILURE
      endif

    case(4)

      call INVERT_4x4(Matrix, Matrix_Inv, Singular_Flag)
      if (Singular_Flag == 1) THEN
        print *, "non-diagonal 4x4 failure"
        Status = Sym%FAILURE
      endif

    case default

      Status = Sym%FAILURE

      print*,"Invalid matrix size."

  end select

  endif

  return

end function Invert_Matrix

!--------------------------------------------------------------------------
! subroutine INVERT_2x2(A,A_inv,ierr)
!
! Matrix Inversion for a 2x2 matrix
!
! A - input matrix
! A_inv - the inverse of A (output)
! ierr - output flag that report a singular matrix and a failure
!--------------------------------------------------------------------------
subroutine INVERT_2x2(A,A_inv,ierr)
  real, dimension(:,:), intent(in):: A
  real, dimension(:,:), intent(out):: A_inv
  real:: determinant
  integer, intent(out):: ierr

!--- compute determinant
  ierr = 0
  determinant = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  if (determinant == 0.0) then
!       print *, "Singular Matrix in Invert 2x2"
        ierr = 1
  endif

!--- compute inverse
  A_inv(1,1) = A(2,2)
  A_inv(1,2) = -A(1,2)
  A_inv(2,1) = -A(2,1)
  A_inv(2,2) = A(1,1)
  A_inv = A_inv / determinant

end subroutine INVERT_2x2
!--------------------------------------------------------------------------
! subroutine INVERT_3x3(A,A_inv,ierr)
!
! Matrix Inversion for a 3x3 matrix
!--------------------------------------------------------------------------
subroutine INVERT_3x3(AA,AA_inv,ierr)
  real, dimension(:,:), intent(in):: AA
  real, dimension(:,:), intent(out):: AA_inv
  integer, intent(out):: ierr
  real(kind=real8) :: determinant
  real(kind=real8), dimension(3,3):: A, A_inv

  ierr = 0

  A = real(AA,kind=real8)

  !--- compute determinant
  determinant = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3)) - &
                A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3)) + &
                A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
  if (determinant == 0.0) then
        ierr = 1
  endif

  !--- compute inverse
  A_inv(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
  A_inv(1,2) = A(1,3)*A(3,2) - A(3,3)*A(1,2)
  A_inv(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)
  A_inv(2,1) = A(2,3)*A(3,1) - A(3,3)*A(2,1)
  A_inv(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
  A_inv(2,3) = A(1,3)*A(2,1) - A(2,3)*A(1,1)
  A_inv(3,1) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
  A_inv(3,2) = A(1,2)*A(3,1) - A(3,2)*A(1,1)
  A_inv(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

  A_inv = A_inv / determinant

  AA_inv = real(A_inv, kind=real4)

end subroutine INVERT_3x3
!--------------------------------------------------------------------------
! subroutine INVERT_4x4(A,A_inv,ierr)
!
! Matrix Inversion for a 4x4 matrix
!--------------------------------------------------------------------------
subroutine INVERT_4x4(AA,AA_inv,ierr)
  real, dimension(:,:), intent(in):: AA
  real, dimension(:,:), intent(out):: AA_inv
  integer, intent(out):: ierr
  real(kind=real8) :: determinant
  real(kind=real8), dimension(4,4):: A, A_inv

  ierr = 0

  A = real(AA,kind=real8)

  !--- compute determinant
  determinant = A(1,1)*(A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3)) + &
                A(1,2)*(A(2,1)*A(3,4)*A(4,3) + A(2,3)*A(3,1)*A(4,4) + A(2,4)*A(3,3)*A(4,1)) + &
                A(1,3)*(A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2)) + &
                A(1,4)*(A(2,1)*A(3,3)*A(4,2) + A(2,2)*A(3,1)*A(4,3) + A(2,3)*A(3,2)*A(4,1)) - &
                A(1,1)*(A(2,2)*A(3,4)*A(4,3) + A(2,3)*A(3,2)*A(4,4) + A(2,4)*A(3,3)*A(4,2)) - &
                A(1,2)*(A(2,1)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,3)) - &
                A(1,3)*(A(2,1)*A(3,4)*A(4,2) + A(2,2)*A(3,1)*A(4,4) + A(2,4)*A(3,2)*A(4,1)) - &
                A(1,4)*(A(2,1)*A(3,2)*A(4,3) + A(2,2)*A(3,3)*A(4,1) + A(2,3)*A(3,1)*A(4,2))

  if (determinant == 0.0) then
        !print *, "Singular Matrix in Invert 4x4 ", determinant
        ierr = 1
  endif

  !--- compute inverse
  A_inv(1,1) = A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) - A(2,2)*A(3,4)*A(4,3) - A(2,3)*A(3,2)*A(4,4) - A(2,4)*A(3,3)*A(4,2)
  A_inv(1,2) = A(1,2)*A(3,4)*A(4,3) + A(1,3)*A(3,2)*A(4,4) + A(1,4)*A(3,3)*A(4,2) - A(1,2)*A(3,3)*A(4,4) - A(1,3)*A(3,4)*A(4,2) - A(1,4)*A(3,2)*A(4,3)
  A_inv(1,3) = A(1,2)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,2) + A(1,4)*A(2,2)*A(4,3) - A(1,2)*A(2,4)*A(4,3) - A(1,3)*A(2,2)*A(4,4) - A(1,4)*A(2,3)*A(4,2)
  A_inv(1,4) = A(1,2)*A(2,4)*A(3,3) + A(1,3)*A(2,3)*A(3,4) + A(1,4)*A(2,3)*A(3,2) - A(1,2)*A(2,3)*A(3,4) - A(1,3)*A(2,4)*A(3,2) - A(1,4)*A(2,2)*A(3,3)

  A_inv(2,1) = A(2,1)*A(3,4)*A(4,3) + A(2,3)*A(3,1)*A(4,4) + A(2,4)*A(3,3)*A(4,1) - A(2,1)*A(3,3)*A(4,4) - A(2,3)*A(3,4)*A(4,1) - A(2,4)*A(3,1)*A(4,3)
  A_inv(2,2) = A(1,1)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,1) + A(1,4)*A(3,1)*A(4,3) - A(1,1)*A(3,4)*A(4,3) - A(1,3)*A(3,1)*A(4,4) - A(1,4)*A(3,3)*A(4,1)
  A_inv(2,3) = A(1,1)*A(2,4)*A(4,3) + A(1,3)*A(2,1)*A(4,4) + A(1,4)*A(2,3)*A(4,1) - A(1,1)*A(2,3)*A(4,4) - A(1,3)*A(2,4)*A(4,1) - A(1,4)*A(2,1)*A(4,3)
  A_inv(2,4) = A(1,1)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,1) + A(1,4)*A(2,1)*A(3,3) - A(1,1)*A(2,4)*A(3,3) - A(1,3)*A(2,1)*A(3,4) - A(1,4)*A(2,3)*A(3,1)

  A_inv(3,1) = A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2) - A(2,1)*A(3,4)*A(4,2) - A(2,2)*A(3,1)*A(4,4) - A(2,4)*A(3,2)*A(4,1)
  A_inv(3,2) = A(1,1)*A(3,4)*A(4,2) + A(1,2)*A(3,1)*A(4,4) + A(1,4)*A(3,2)*A(4,1) - A(1,1)*A(3,2)*A(4,4) - A(1,2)*A(3,4)*A(4,1) - A(1,4)*A(3,1)*A(4,2)
  A_inv(3,3) = A(1,1)*A(2,2)*A(4,4) + A(1,2)*A(2,4)*A(4,1) + A(1,4)*A(2,1)*A(4,2) - A(1,1)*A(2,4)*A(4,2) - A(1,2)*A(2,1)*A(4,4) - A(1,4)*A(2,2)*A(4,1)
  A_inv(3,4) = A(1,1)*A(2,4)*A(3,2) + A(1,2)*A(2,1)*A(3,4) + A(1,4)*A(2,2)*A(3,1) - A(1,1)*A(2,2)*A(3,4) - A(1,2)*A(2,4)*A(3,1) - A(1,4)*A(2,1)*A(3,2)

  A_inv(4,1) = A(2,1)*A(3,3)*A(4,2) + A(2,2)*A(3,1)*A(4,3) + A(2,3)*A(3,2)*A(4,1) - A(2,1)*A(3,2)*A(4,3) - A(2,2)*A(3,3)*A(4,1) - A(2,3)*A(3,1)*A(4,2)
  A_inv(4,2) = A(1,1)*A(3,2)*A(4,3) + A(1,2)*A(3,3)*A(4,1) + A(1,3)*A(3,1)*A(4,2) - A(1,1)*A(3,3)*A(4,2) - A(1,2)*A(3,1)*A(4,3) - A(1,3)*A(3,2)*A(4,1)
  A_inv(4,3) = A(1,1)*A(2,3)*A(4,2) + A(1,2)*A(2,1)*A(4,3) + A(1,3)*A(2,2)*A(4,1) - A(1,1)*A(2,2)*A(4,3) - A(1,2)*A(2,3)*A(4,1) - A(1,3)*A(2,1)*A(4,2)
  A_inv(4,4) = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) - A(1,3)*A(2,2)*A(3,1)

  A_inv = A_inv / determinant

  AA_inv = real(A_inv, kind=real4)

end subroutine INVERT_4x4
!--------------------------------------------------------------------------
! subroutine INVERT_DIAGONAL(A,A_inv,ierr)
!
! Matrix Inversion for a nxn diagnonal matrix
!--------------------------------------------------------------------------
subroutine INVERT_DIAGONAL(A,A_inv,ierr)
  real, dimension(:,:), intent(in):: A
  real, dimension(:,:), intent(out):: A_inv
  integer, intent(out):: ierr
  integer:: i, j, ni, nj
  real:: determinant
  real:: zero

  zero = epsilon(zero)
  ierr = 0
  A_inv = 0.0

  ni = size(A,1)
  nj = size(A,2)
 
  if(ni /= nj) then 
    ierr = 1
    return
  endif

  do i = 1, ni
     j = i
     if (abs(A(i,j)) > zero) then
        A_inv(i,i) = 1.0 / A(i,i)
     else
       ierr = 1
       return
     endif
  enddo

end subroutine INVERT_DIAGONAL

!--------------------------------------------------------------------
! subroutine FIND_BOUNDS(lat,lon,wlon,elon,slat,nlat,dateline_flg)
!
! This subroutine reads in a satellite data set containing lat/lon
! values and finds the western/eastern-most longitude and
! southern/northern-most latitude values.
!--------------------------------------------------------------------

subroutine FIND_BOUNDS(lat,lon,wlon,elon,slat,nlat,dateline_flg)

!Inputs and outputs are as follows:

!Inputs:
! lat    - array of satellite latitudes
! lon    - array of satellite longitudes

!Outputs
! wlon   - western-most longitude value in satellite data
! elon   - eastern-most longitude value in satellite data
! slat   - southern-most latitude value in satellite data
! nlat   - northern-most latitude value in satellite data

!Declare calling parameters

  real(kind=real4),dimension(:,:),intent(in) :: lat,lon
  real(kind=real4),intent(out) :: wlon,elon,slat,nlat
  integer(kind=int1), intent(out) :: dateline_flg

  integer(kind=int4) :: astatus, nx, ny
  real(kind=real4), dimension(:,:), allocatable :: dum

  nx = size(lon,1)
  ny = size(lon,2)

  allocate(dum(nx,ny),stat=astatus)
  if (astatus /= 0) then
    print "(a,'Not enough memory to allocate dummy longitude array.')",EXE_PROMPT
    stop
  endif

  dum = lon

  nlat = missing_value_real4
  slat = missing_value_real4
  elon = missing_value_real4
  wlon = missing_value_real4

  nlat = maxval(lat, mask = lat >= -90.0 .and. lat <= 90.0)
  slat = minval(lat, mask = lat >= -90.0 .and. lat <= 90.0)

  elon = maxval(lon, mask = lon >= -180.0 .and. lon <= 180.0)
  wlon = minval(lon, mask = lon >= -180.0 .and. lon <= 180.0)

  !Check if intl. date line is crossed
  dateline_flg = 0
  if (elon > 160.0 .and. wlon < -160) then
    where(lon < 0.0 .and. lon > missing_value_real4)
      dum = lon + 360.0
    endwhere
    !write(*,*) "international date line is assumed to be crossed"
    elon = maxval(dum, mask = lon >= -180.0 .and. lon <= 180.0)
    wlon = minval(dum, mask = lon >= -180.0 .and. lon <= 180.0)
    dateline_flg = 1
  endif

  deallocate(dum, stat=astatus)
  if (astatus /= 0) then
    print "(a,'Error deallocating dummy longitude array.')",EXE_PROMPT
    stop
  endif

end subroutine FIND_BOUNDS

!------------------------------------------------------------------------
! subroutine PACK_BYTES_I1(input_bytes,bit_depth,output_bytes)
! subroutine PACK_BYTES_I2(input_bytes,bit_depth,output_bytes)
!
! Routines to pack individual bytes into a single byte
!
! input:
! input_bytes - vector of bytes to be packed into output_byte
! bit_start - vector of bit starting positions (1-7) for each input byte
! bit_depth - vector of bit depths (1-7) for each input byte (total can not exceed 8)
!
! output: 
! output_byte - byte variable that holds the bit values of the input_bytes - can be i1 or i2
!
! local
!  n_in - number of elements of input vectors
!  i_in - index of input_bytes (1-n_in)
!  n_out - number of elements of output vectors
!  i_out - index of output_bytes (1-n_out)
!
! Note:
! 1.  if the input byte has information in bits greater then bit depth - they are removed 
!
!
! Example, pack an input_byte wth  bit_start = 2 and bit depth 3 
!
! input byte
!           x x x
! _ _ _ _ _ _ _ _
!
! result of first ishft
! x x x
! _ _ _ _ _ _ _ _
!
! result of second ishft
!       x x x
! _ _ _ _ _ _ _ _
!
! Author: Andrew Heidinger
!
! Version History:  
! February 2006 - Created
!-----------------------------------------------------------------------------------

!--- This Version packs into one byte words
   subroutine PACK_BYTES_I1(input_bytes,bit_depth,output_bytes)
    integer(kind=int1), dimension(:), intent(in):: input_bytes
    integer(kind=int4), dimension(:), intent(in):: bit_depth
    integer(kind=int1), dimension(:), intent(out):: output_bytes
    integer(kind=int1):: bit_start, bit_end, bit_offset
    integer(kind=int1):: temp_byte
    integer:: n_in,i_in,n_out,i_out
    integer, parameter:: word_bit_depth = 8
  
!--- determine size of vectors
    n_in = size(input_bytes)
    n_out = size(output_bytes)

!--- reset output byte
   output_bytes = 0

!--- initialize
   bit_offset = 0
   bit_start = 0
   bit_end = 0
   i_out = 1

!--- loop through input bytes
   do i_in = 1, n_in

!--- determine starting and ending bit locations
     bit_start = bit_offset + 1
     bit_end = bit_start + bit_depth(i_in) - 1

!--- determine if this input byte will fit on current output byte, if not go to next
     if (bit_end > word_bit_depth) then
      i_out = i_out + 1
      bit_offset = 0
      bit_start = bit_offset + 1
      bit_end = bit_start + bit_depth(i_in) - 1
     endif

!--- check for exceeding the space allowed for the packed bytes
     if (i_out > n_out) then
       print *, "ERROR: Insufficient space for bit packing" 
       return
     endif

!--- place input byte into correct position
     temp_byte =0
     temp_byte = ishft(input_bytes(i_in),word_bit_depth-bit_depth(i_in))   !first ishft
     temp_byte = ishft(temp_byte,bit_end - word_bit_depth)                 !second ishft

!--- modify output byte
     output_bytes(i_out) = output_bytes(i_out) + temp_byte

!--- update bit offset
     bit_offset = bit_offset + bit_depth(i_in)

   enddo

  end subroutine  PACK_BYTES_I1

!--- This Version packs into two byte words
   subroutine PACK_BYTES_I2(input_bytes,bit_depth,output_bytes)
    integer(kind=int1), dimension(:), intent(in):: input_bytes
    integer(kind=int4), dimension(:), intent(in):: bit_depth
    integer(kind=int2), dimension(:), intent(out):: output_bytes
    integer(kind=int1):: bit_start, bit_end, bit_offset
    integer(kind=int2):: temp_byte                 
    integer:: n_in,i_in,n_out,i_out
    integer, parameter:: word_bit_depth = 16

!--- determine size of vectors
    n_in = size(input_bytes)
    n_out = size(output_bytes)

!--- reset output byte
   output_bytes = 0

!--- initialize
   bit_offset = 0
   bit_start = 0
   bit_end = 0
   i_out = 1

!--- loop through input bytes
   do i_in = 1, n_in

!--- determine starting and ending bit locations
     bit_start = bit_offset + 1
     bit_end = bit_start + bit_depth(i_in) - 1

!--- determine if this input byte will fit on current output byte, if not go to next
     if (bit_end > word_bit_depth) then
      i_out = i_out + 1
      bit_offset = 0
      bit_start = bit_offset + 1
      bit_end = bit_start + bit_depth(i_in) - 1
     endif

!--- check for exceeding the space allowed for the packed bytes
     if (i_out > n_out) then
       print *, "ERROR: Insufficient space for bit packing"
       return
     endif

!--- place input byte into correct position
     temp_byte =0
     temp_byte = ishft(INT(input_bytes(i_in),kind=INT2),word_bit_depth-bit_depth(i_in))   !first ishft
     temp_byte = ishft(temp_byte,bit_end - word_bit_depth)                 !second ishft

!--- modify output byte
     output_bytes(i_out) = output_bytes(i_out) + temp_byte

!--- update bit offset
     bit_offset = bit_offset + bit_depth(i_in)

   enddo

  end subroutine  PACK_BYTES_I2
 !----------------------------------------------------------------------
 !  Local Linear Radiative Center
 !----------------------------------------------------------------------
 subroutine GRADIENT_MEANDER(Meander_Flag, &
                             Grid_Data, &
                             Element_Start, Number_Of_Elements, & 
                             Line_Start, Number_Of_Lines, & 
                             Max_Grad_Distance, &
                             Grad_Flag,  &
                             Missing_LRC_Value, &
                             Skip_LRC_Mask, &
                             Min_Grid_Data_Valid, Max_Grid_Data_Valid, &
                             ielem_LRC, iline_LRC)

  integer, intent(in):: Meander_Flag
  real (kind=real4), intent(in), dimension(:,:) :: Grid_Data
  integer (kind=int4), intent(in):: Element_Start
  integer (kind=int4), intent(in):: Number_of_Elements
  integer (kind=int4), intent(in):: Line_Start
  integer (kind=int4), intent(in):: Number_of_Lines
  integer (kind=int4), intent(in):: Max_Grad_Distance
  integer (kind=int4), intent(in):: Grad_Flag
  integer (kind=int4), intent(in):: Missing_LRC_Value
  integer (kind=int1), intent(in), dimension(:,:):: Skip_LRC_Mask
  real (kind=real4), intent(in):: Min_Grid_Data_Valid
  real (kind=real4), intent(in):: Max_Grid_Data_Valid
  integer (kind=int4), intent(out), dimension(:,:):: ielem_LRC
  integer (kind=int4), intent(out), dimension(:,:):: iline_LRC
  real, dimension(3,3):: Grad_Array
  integer, dimension(2):: Grad_Indices
  integer:: ielem
  integer:: iline
  integer:: ielem_Previous
  integer:: iline_Previous
  integer:: ielem_Next
  integer:: iline_Next
  real:: Grad_Temp
  integer:: Element_End
  integer:: Line_End
  integer:: ipoint
  integer:: ielem_dir
  integer:: iline_dir
 
  Element_End = Number_of_Elements + Element_Start - 1
  Line_End = Number_of_Lines + Line_Start - 1

  !--- initialize
  ielem_LRC = Missing_LRC_Value
  iline_LRC = Missing_LRC_Value

!----------------------------------------------------------------------
! loop through pixels in segment
!----------------------------------------------------------------------
Element_Loop:  do ielem = Element_Start+1, Element_End-1
Line_Loop:    do iline = Line_Start+1, Line_End-1

      !--- skip data due to mask
      if (Skip_LRC_Mask(ielem,iline) == sym%YES) cycle

      !-- check for out of bounds data
      if (Grad_Flag ==  1 .and. Grid_Data(ielem,iline) < Min_Grid_Data_Valid) cycle
      if (Grad_Flag ==  -1 .and. Grid_Data(ielem,iline) > Max_Grid_Data_Valid) cycle

      !-- check for data that already meets LRC criteria
      if ((Grad_Flag ==  1 .and. Grid_Data(ielem,iline) > Max_Grid_Data_Valid) .or. &
          (Grad_Flag ==  -1 .and. Grid_Data(ielem,iline) < Min_Grid_Data_Valid)) then
              ielem_LRC(ielem,iline) = ielem
              iline_LRC(ielem,iline) = iline
      endif

      !--- initialize previous variables
      ielem_Previous = ielem
      iline_Previous = iline

      !---- go long gradient and check for a reversal or saturation
      do ipoint = 1,Max_Grad_Distance

        !--- compute local gradient, find strongest gradient in 3x3 array and compute direction
        if (ipoint == 1 .or. Meander_Flag == sym%YES) then

         !--- construct 3x3 array for analysis
         Grad_Array =  &
           Grid_Data(ielem_Previous-1:ielem_Previous+1,iline_Previous-1:iline_Previous+1) -  &
           Grid_Data(ielem_Previous,iline_Previous)

         !--- look for bad data
         if (minval(Grad_Array) == Missing_Value_Real4) exit 

         !--- compute local gradients, find strongest gradient
         if (Grad_Flag == 1) then
          Grad_Indices = maxloc(Grad_Array)
         else
          Grad_Indices = minloc(Grad_Array)
         endif 

         !--- compute direction
         ielem_Dir = Grad_Indices(1)  - 2
         iline_Dir = Grad_Indices(2)  - 2

         !--- check for pixels that are located at  minima/maxima
         if (ielem_Dir == 0 .and. iline_Dir == 0) then
           ielem_LRC(ielem,iline) = ielem_Previous
           iline_LRC(ielem,iline) = iline_Previous
           exit
         endif

        endif

        !-- select next point on the path
        ielem_Next = ielem_Previous + ielem_Dir
        iline_Next = iline_Previous + iline_Dir

        !--- check for hitting segment boundaries
        if (ielem_Next == Element_Start .or. ielem_Next == Element_End .or. &
             iline_Next == Line_Start .or. iline_Next == Line_End) then
              ielem_LRC(ielem,iline) = ielem_Previous
              iline_LRC(ielem,iline) = iline_Previous
              exit
         endif

         !--- check for hitting bad data
         if (Skip_LRC_Mask(ielem_Next,iline_Next) == sym%YES) then
              ielem_LRC(ielem,iline) = ielem_Previous
              iline_LRC(ielem,iline) = iline_Previous
              exit
         endif

         !--- check for sign reversal
         if (Meander_Flag == sym%NO) then

          Grad_Temp = Grid_Data(ielem_Next,iline_Next) -  &
                      Grid_Data(ielem_Previous,iline_Previous)

          if (Grad_Flag * Grad_Temp < 0) then
              ielem_LRC(ielem,iline) = ielem_Previous
              iline_LRC(ielem,iline) = iline_Previous
              exit
          endif
         endif

         !--- check for saturation
         if (Grad_Flag == 1 .and. Grid_Data(ielem_Next,iline_Next) > Max_Grid_Data_Valid) then
              ielem_LRC(ielem,iline) = ielem_Next
              iline_LRC(ielem,iline) = iline_Next
              exit
         endif
         if (Grad_Flag == -1 .and. Grid_Data(ielem_Next,iline_Next) < Min_Grid_Data_Valid) then
              ielem_LRC(ielem,iline) = ielem_Next
              iline_LRC(ielem,iline) = iline_Next
              exit
         endif

         !--- store position
         ielem_Previous = ielem_Next
         iline_Previous = iline_Next

      enddo

    end do Line_Loop
  end do Element_Loop

end subroutine GRADIENT_MEANDER

!==============================================================
! subroutine COMPUTE_MEDIAN(z,mask,z_median,z_mean,z_std_median)
!
! Median filter
!==============================================================
subroutine COMPUTE_MEDIAN(z,mask,z_median,z_mean,z_std_median)

! The purpose of this function is to find 
! median (emed), minimum (emin) and maximum (emax)
! for the array elem with nelem elements. 

 real, dimension(:,:), intent(in):: z
 real, intent(out):: z_median
 real, intent(out):: z_mean
 real, intent(out):: z_std_median
 integer(kind=int1), dimension(:,:), intent(in):: mask 
 integer:: i,j,k,nx,ny,nelem
 real, dimension(:), allocatable::x
 real(kind=real4):: u

 z_median = missing_value_real4
 z_std_median = missing_value_real4
 z_mean = missing_value_real4

 nx = size(z,1)
 ny = size(z,2)

 nelem = nx * ny

 allocate(x(nelem))
 x = 0.0

 k = 0
 do i = 1, nx
   do j = 1, ny
      if (mask(i,j) == sym%NO .and. z(i,j) /= missing_value_real4) then
           k = k + 1   
           x(k) = z(i,j)
      endif
  enddo
 enddo

 nelem = k
   
 if (nelem < 1) then
     if (allocated(x)) deallocate(x)
     return
 endif 
!--- sort the array into ascending order
  do i=1,nelem-1
   do j=i+1,nelem
    if(x(j)<x(i))then
     u=x(j)
     x(j)=x(i)
     x(i)=u
    end if   
   end do
  end do

!---- pick the median
  if(mod(nelem,2)==1)then
   i=nelem/2+1
   z_median=x(i)
  else  
   i=nelem/2
   z_median=(x(i)+x(i+1))/2
   end if

!--- compute standard deviation wrt median
  z_mean = sum(x(1:nelem))/nelem
  z_std_median = sqrt(sum((x(1:nelem) - z_median)**2) / nelem)


! if (z_std_median > 60.0) then 
!         print *, "big std median ", z_std_median, nelem, x(1:nelem)
!         print *, "z_nxn = ", z
! endif

  if (allocated(x)) deallocate(x)

end subroutine COMPUTE_MEDIAN

!----------------------------------------------------------------------
! subroutine COMPUTE_MEDIAN_SEGMENT(z,mask,n,imin,imax,jmin,jmax,
!                                   z_median,z_std_median)
!
! Compute standard deviaion of an array wrt to the median
!----------------------------------------------------------------------
subroutine COMPUTE_MEDIAN_SEGMENT(z,mask,n,imin,imax,jmin,jmax, &
                                  z_median, &
                                  z_std_median)
  real(kind=real4), dimension(:,:), intent(in):: z
  integer(kind=int1), dimension(:,:), intent(in):: mask
  real(kind=real4), dimension(:,:), intent(out):: z_std_median
  real(kind=real4), dimension(:,:), intent(out):: z_median
! real(kind=real4), dimension(:,:), intent(out):: z_mean
  integer, intent(in):: n
  integer, intent(in):: imin
  integer, intent(in):: imax
  integer, intent(in):: jmin
  integer, intent(in):: jmax
  integer:: i
  integer:: j
  integer:: i1
  integer:: i2
  integer:: j1
  integer:: j2
  real(kind=real4) :: z_mean

  do i = imin, imax
    do j = jmin, jmax

     j1 = max(jmin,j-n)   !top index of local array
     j2 = min(jmax,j+n)   !bottom index of local array
     i1 = max(imin,i-n)   !left index of local array
     i2 = min(imax,i+n)   !right index of local array

     !--- compute median
     call COMPUTE_MEDIAN(z(i1:i2,j1:j2),mask(i1:i2,j1:j2),z_median(i,j), &
                         z_mean,z_std_median(i,j))

     enddo
  enddo

end subroutine COMPUTE_MEDIAN_SEGMENT
!====================================================================
! Function Name: Pearson_Corr
!
! Function:
!    Compute the Pearson Correlation Coefficient for two mxn arrays
!
! Description: Pearson's product-moment coefficient
!   
! Calling Sequence: BT_WV_BT_Window_Corr(Elem_Idx,Line_Idx) = Pearson_Corr( &
!                       sat%bt10(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                       sat%bt14(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                      Array_Width, Array_Hgt)
!   
!
! Inputs:
!   Array 1
!   Array 2
!   Elem_size
!   Line_size
!
! Outputs: 
!   Pearson Correlation coefficent
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: Standard definition for Pearson correlation
!
!====================================================================
function Pearson_Corr(Array_One,Array_Two,Array_Width,Array_Hght,Invalid_Data_Mask)  &
         RESULT(Pearson_Corr_Coeff)
   real(kind=real4), intent(in), dimension(:,:):: Array_One
   real(kind=real4), intent(in), dimension(:,:):: Array_Two
   integer(kind=INT4), intent(in):: Array_Width
   integer(kind=INT4), intent(in):: Array_Hght
   integer(kind=INT1), intent(in), dimension(:,:):: Invalid_Data_Mask
   real(kind=real4), dimension(Array_Width,Array_Hght):: Pearson_Corr_Term_1
   real(kind=real4), dimension(Array_Width,Array_Hght):: Pearson_Corr_Term_2
   real(kind=real8):: Pearson_Corr_Top_Term_1
   real(kind=real8):: Pearson_Corr_Top_Term_2
   real(kind=real8):: Pearson_Corr_Bottom_Term_1
   real(kind=real8):: Pearson_Corr_Bottom_Term_2
   real(kind=real4):: Pearson_Corr_Coeff
   real(kind=real8):: Mean_Array_One
   real(kind=real8):: Mean_Array_Two
   real(kind=real8):: Sum_Array_One
   real(kind=real8):: Sum_Array_Two

   !--- skip computation for pixel arrays with any missing data
   if (sum(Invalid_Data_Mask) > 0) then
      Pearson_Corr_Coeff = Missing_Value_Real4
      return
   endif

   Sum_Array_One = sum(Array_One)
   Sum_Array_Two = sum(Array_Two)

   Mean_Array_One = Sum_Array_One / (Array_Width*Array_Hght)
   Mean_Array_Two = Sum_Array_Two / (Array_Width*Array_Hght)

   Pearson_Corr_Term_1 = Array_One - Mean_Array_One
   Pearson_Corr_Term_2 = Array_Two - Mean_Array_Two

   Sum_Array_One = sum(Pearson_Corr_Term_1)
   Sum_Array_Two = sum(Pearson_Corr_Term_2)

   Mean_Array_One = 0.0
   Mean_Array_Two = 0.0

   Pearson_Corr_Top_Term_1 = sum(Pearson_Corr_Term_1*Pearson_Corr_Term_2)
   
   Pearson_Corr_Top_Term_2 = (Sum_Array_One*Sum_Array_Two) / (Array_Width*Array_Hght)
   
   Pearson_Corr_Bottom_Term_1 = sum(Pearson_Corr_Term_1**2) - &
                                ((Sum_Array_One)**2) / (Array_Width*Array_Hght)
                                 
   Pearson_Corr_Bottom_Term_2 = sum(Pearson_Corr_Term_2**2) - &
                                ((Sum_Array_Two)**2) / (Array_Width*Array_Hght)

   Pearson_Corr_Coeff = (Pearson_Corr_Top_Term_1 - Pearson_Corr_Top_Term_2) / &
                         sqrt(Pearson_Corr_Bottom_Term_1 * &
                              Pearson_Corr_Bottom_Term_2)
   
 end function Pearson_Corr

!====================================================================
! Function Name: Covariance
!
! Function:
!    Compute the Covariance for two mxn arrays
!
! Description: Covariance = E(XY) - E(X)*E(Y)
!   
! Calling Sequence: BT_WV_BT_Window_Covar(Elem_Idx,Line_Idx) = Covariance( &
!                       sat%bt10(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                       sat%bt14(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                      Array_Width, Array_Hgt)
!   
!
! Inputs:
!   Array 1 - the first array (X)
!   Array 2 - the second array (Y)
!   Elem_size
!   Line_size
!
! Outputs: 
!   Covariance of X and Y
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: Standard definition for the Covariance Computation
!
!====================================================================
function Covariance(Array_One,Array_Two,Array_Width,Array_Hght,Invalid_Data_Mask) &
         RESULT(Covar_Array_One_Array_Two)
   real(kind=real4), intent(in), dimension(:,:):: Array_One
   real(kind=real4), intent(in), dimension(:,:):: Array_Two
   integer(kind=INT4), intent(in):: Array_Width
   integer(kind=INT4), intent(in):: Array_Hght
   integer(kind=INT1), intent(in), dimension(:,:):: Invalid_Data_Mask

   real(kind=real8):: Mean_Array_One
   real(kind=real8):: Mean_Array_Two
   real(kind=real8):: Mean_Array_One_x_Array_Two
   real(kind=real8):: Sum_Array_One
   real(kind=real8):: Sum_Array_Two
   real(kind=real8):: Sum_Array_One_x_Array_Two
   real(kind=real4):: Covar_Array_One_Array_Two

   !--- skip computation for pixel arrays with any missing data
   if (sum(Invalid_Data_Mask) > 0) then
      Covar_Array_One_Array_Two = Missing_Value_Real4
      return
   endif

   Sum_Array_One = sum(Array_One)
   Sum_Array_Two = sum(Array_Two)

   Mean_Array_One = Sum_Array_One / (Array_Width*Array_Hght)
   Mean_Array_Two = Sum_Array_Two / (Array_Width*Array_Hght)

   Sum_Array_One_x_Array_Two = sum(Array_One*Array_Two)
   Mean_Array_One_x_Array_Two = Sum_Array_One_x_Array_Two / (Array_Width*Array_Hght)
   
   Covar_Array_One_Array_Two  = Mean_Array_One_x_Array_Two - &
                                Mean_Array_One * Mean_Array_Two 
   
 end function Covariance
 





!------------------------------------------------------------------------------------- 

end module NUMERICAL_ROUTINES
