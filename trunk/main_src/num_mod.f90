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
  
  public:: LOCATE, &
           JULIAN, &
           COMPUTE_MONTH, &
           COMPUTE_DAY, &
           VAPOR, &
           VAPOR_ICE, &
           INVERT_MATRIX,  &
           INVERT_2x2,  &
           INVERT_3x3,  &
           INVERT_4x4,  &
           INVERT_DIAGONAL,  &
           FIND_BOUNDS,  &
           COMPUTE_TIME_HOURS, &
           COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES,  &
           GRADIENT_MEANDER, &
           COMPUTE_MEDIAN, &
           COMPUTE_MEDIAN_SEGMENT, &
           LEAP_YEAR_FCT, &
           WIND_SPEED, &
           WIND_DIRECTION, &
           COUNTSUBSTRING, &
           SPLIT_STRING, &
           REPLACE_CHAR_IN_STRG, &
           svdcmp, pythag, &
           GEN_INVERSE           
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

!-------------------------------------------------
! subroutine JULIAN(iday,imonth,iyear,jday)
! compute julian day
! input:
!         iday - integer day
!         imonth - integer month
!         iyear - integer year (2 or four digits)
!         
! output : jday - julian day
!--------------------------------------------------
 subroutine JULIAN(iday,imonth,iyear,jday)

!-- Computes julian day (1-365/366)
        integer, intent(in)::  iday,imonth,iyear
        integer, intent(out):: jday
        integer::  j
        integer, dimension(12)::  jmonth

        jmonth = reshape ((/31,28,31,30,31,30,31,31,30,31,30,31/),(/12/))

        jday = iday
        if (modulo(iyear,4) == 0) then
            jmonth(2)=29
        endif

        do j = 1,imonth-1
           jday = jday + jmonth(j)
        end do

   end subroutine JULIAN

!--------------------------------------------
! compute the month
!---------------------------------------------
 function COMPUTE_MONTH(jday,ileap) result(month)
   integer, intent(in):: ileap
   integer, intent(in):: jday
   integer:: month

   month = 0
   if (jday < 32) then
     month = 1
   elseif (jday < 60+ileap) then
     month = 2
   elseif (jday < 91+ileap) then
     month = 3
   elseif (jday < 121+ileap) then
     month = 4
   elseif (jday < 152+ileap) then
     month = 5
   elseif (jday < 182+ileap) then
     month = 6
   elseif (jday < 213+ileap) then
     month = 7
   elseif (jday < 244+ileap) then
     month = 8
   elseif (jday < 274+ileap) then
     month = 9
   elseif (jday < 305+ileap) then
     month = 10
   elseif (jday < 335+ileap) then
     month = 11
   else
     month = 12
   endif

 end function COMPUTE_MONTH
!--------------------------------------------
! compute the day
!---------------------------------------------
 function COMPUTE_DAY(jday,ileap) result(day)
   integer, intent(in):: ileap
   integer, intent(in):: jday
   integer:: day

   if (jday < 32) then
     day = jday
   elseif (jday < 60) then
     day = jday - 31
   elseif ((jday == 60).and.(ileap == 1)) then
     day = jday - 31
   elseif ((jday == 60).and.(ileap == 0)) then
     day = jday - 59
   elseif (jday < 91+ileap) then
     day = jday - (59 + ileap)
   elseif (jday < 121+ileap) then
     day = jday - (90 + ileap)
   elseif (jday < 152+ileap) then
     day = jday - (120 + ileap)
   elseif (jday < 182+ileap) then
     day = jday - (151 + ileap)
   elseif (jday < 213+ileap) then
     day = jday - (181 + ileap)
   elseif (jday < 244+ileap) then
     day = jday - (212 + ileap)
   elseif (jday < 274+ileap) then
     day = jday - (243 + ileap)
   elseif (jday < 305+ileap) then
     day = jday - (273 + ileap)
   elseif (jday < 335+ileap) then
     day = jday - (304 + ileap)
   else
     day = jday - (334 + ileap)
   endif

   return

 end function COMPUTE_DAY

 !---- function to return the system in fractional hours
 function COMPUTE_TIME_HOURS()   result(time_hours)
   character(len=8):: system_date
   character(len=10):: system_time
   character(len=5):: system_time_zone
   integer, dimension(8):: system_time_value

   real:: time_hours

   call DATE_AND_TIME(system_date,system_time,system_time_zone, system_time_value)

   time_hours = real(system_time_value(5)) +  &
                     (real(system_time_value(6)) + &
                      real(system_time_value(7) + &
                      real(system_time_value(8))/1000.0)/60.0)/60.0
   return

 end function COMPUTE_TIME_HOURS

!----------------------------------------------------------------
! functions to compute some needed water vapor parameters
!----------------------------------------------------------------
 function VAPOR(T) result(es)
                                                                     
!  T in Kelvin                                                          
!  es in mbar

  implicit none
  real, intent (in) :: T
  real :: es

   es = 6.112 * exp(17.67 * (T-273.16) / (T - 29.66))

  return 
end function VAPOR

!---- saturation vapor pressure for ice
function VAPOR_ICE(T) result(es)
   implicit none
   real, intent(in):: T
   real:: es
     es = 6.1078 * exp(21.8745584 * (T-273.16) / (T - 7.66))
  return 
end function VAPOR_ICE

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

      call GEN_INVERSE(Matrix,Matrix_Inv,Matrix_Size,Matrix_Size,Singular_Flag)
      if (Singular_Flag == 1) THEN
        print *, "non-diagonal failure, matrix size greater than 4"
        Status = Sym%FAILURE
        !print*,"Invalid matrix size."
      endif

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
 
!---------------------------------------------------------------------
! function leap_year_fct(yyyy) result(leap_flg)
!
! Function to determine if an input year is a leap year.
!---------------------------------------------------------------------

function leap_year_fct(yyyy) result(leap_flg)
  integer, intent(in) :: yyyy
  integer :: leap_flg

  leap_flg = 0

  if ((modulo(yyyy,4) == 0 .and. modulo(yyyy,100) /= 0) .or. &
       modulo(yyyy,400) == 0) leap_flg = 1

  return

end function leap_year_fct

!---------------------------------------------------------------------
! SUBPROGRAM:  W3FC05        EARTH U,V WIND COMPONENTS TO DIR AND SPD
!   PRGMMR: CHASE            ORG: NMC421      DATE:88-10-26
!
! ABSTRACT: GIVEN THE TRUE (EARTH ORIENTED) WIND COMPONENTS
!   COMPUTE THE WIND DIRECTION AND SPEED.
!   INPUT WINDS AT THE POLE ARE ASSUMED TO FOLLOW THE WMO
!   CONVENTIONS, WITH THE OUTPUT DIRECTION COMPUTED IN ACCORDANCE
!   WITH WMO STANDARDS FOR REPORTING WINDS AT THE POLE.
!   (SEE OFFICE NOTE 241 FOR WMO DEFINITION.)
!
! PROGRAM HISTORY LOG:
!   81-12-30  STACKPOLE, JOHN
!   88-10-19  CHASE, P.   ALLOW OUTPUT VALUES TO OVERLAY INPUT
!   89-01-21  R.E.JONES   CONVERT TO MICROSOFT FORTRAN 4.10
!   90-06-11  R.E.JONES   CONVERT TO SUN FORTRAN 1.3
!   91-03-30  R.E.JONES   SiliconGraphics FORTRAN
!
! USAGE:    call W3FC05 (U, V, DIR, SPD)
!
!   INPUT ARGUMENT LIST:
!     U        - real*4 EARTH-ORIENTED U-COMPONENT
!     V        - real*4 EARTH-ORIENTED V-COMPONENT
!
!   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
!     DIR      - real*4 WIND DIRECTION, DEGREES.  VALUES WILL
!                BE FROM 0 TO 360 INCLUSIVE.
!     SPD      - real*4 WIND SPEED IN SAME UNITS AS INPUT
!---------------------------------------------------------------------

subroutine WIND_SPEED_AND_DIRECTION(u,v,dir,spd)
  real, intent(in) :: u
  real, intent(in) :: v
  real, intent(out) :: spd
  real, intent(out) :: dir
 

  real, parameter:: SPDTST = 1.0e-10
  real, parameter:: RTOD = 57.2957795
  real, parameter:: dchalf = 180.0

  spd = Missing_Value_Real4
  dir = Missing_Value_Real4
  if (u == Missing_Value_Real4 .or. v == Missing_Value_Real4) then
     return
  endif 
  spd = sqrt(u * u + v * v)
  if (spd < SPDTST) THEN
        dir = 0.0
  else 
        dir = atan2(u,v) * RTOD + DCHALF
  endif

  return

end subroutine WIND_SPEED_AND_DIRECTION
!-------------------------------------------------------------------------------------
! elemental funcions wind_speed and wind_direction taken from W3FC03 from UCAR
!
! W3FC05 header documentation follows:
!
! ABSTRACT: GIVEN THE TRUE (EARTH ORIENTED) WIND COMPONENTS
!   COMPUTE THE WIND DIRECTION AND SPEED.
!   INPUT WINDS AT THE POLE ARE ASSUMED TO FOLLOW THE WMO
!   CONVENTIONS, WITH THE OUTPUT DIRECTION COMPUTED IN ACCORDANCE
!   WITH WMO STANDARDS FOR REPORTING WINDS AT THE POLE.
!   (SEE OFFICE NOTE 241 FOR WMO DEFINITION.)
!
! PROGRAM HISTORY LOG:
!   81-12-30  STACKPOLE, JOHN
!   88-10-19  CHASE, P.   ALLOW OUTPUT VALUES TO OVERLAY INPUT
!   89-01-21  R.E.JONES   CONVERT TO MICROSOFT FORTRAN 4.10
!   90-06-11  R.E.JONES   CONVERT TO SUN FORTRAN 1.3
!   91-03-30  R.E.JONES   SiliconGraphics FORTRAN
!
! USAGE:    CALL W3FC05 (U, V, DIR, SPD)
!
!   INPUT ARGUMENT LIST:
!     U        - real*4 EARTH-ORIENTED U-COMPONENT
!     V        - real*4 EARTH-ORIENTED V-COMPONENT
!
!   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
!     DIR      - real*4 WIND DIRECTION, DEGREES.  VALUES WILL
!                BE FROM 0 TO 360 INCLUSIVE.
!     SPD      - real*4 WIND SPEED IN SAME UNITS AS INPUT
!-------------------------------------------------------------------------------------
real elemental function WIND_SPEED ( u ,v )

  real, intent(in)::  u
  real, intent(in):: v

  if (u == Missing_Value_Real4 .or. v == Missing_Value_Real4) then
     wind_speed = Missing_Value_Real4
  else
     wind_speed = sqrt(u * u + v * v)
  endif

end function wind_speed
!-------------------------------------------------------------------------------------
! taken from W3FC03 from UCAR
!-------------------------------------------------------------------------------------
real elemental function WIND_DIRECTION ( u ,v )

  real, intent(in)::  u
  real, intent(in):: v
  real, parameter:: rtod = 57.2957795
  real, parameter:: dchalf = 180.0
  real, parameter:: spdtst = 1.0e-10

  if (u == Missing_Value_Real4 .or. v == Missing_Value_Real4) then
     wind_direction = Missing_Value_Real4
  else
     if (abs(u) < spdtst .and. abs(v) < spdtst) then
      wind_direction =  0.0
     else
      wind_direction = atan2(u,v) * rtod + dchalf
     endif
  endif

end function wind_direction

!------------------------------------------------------------------------------------- 

!-------------------------------------------------------------------------------------                                                   
! count how many characters in the string 
!-------------------------------------------------------------------------------------                                                                 
function COUNTSUBSTRING (s1, s2) result(c)
  character(*), intent(in) :: s1, s2
  integer :: c, p, posn

  c = 0
  if(len(s2) == 0) return
  p = 1
  do
    posn = index(s1(p:), s2)
    if(posn == 0) return
    c = c + 1
    p = p + posn + len(s2)
  end do

end function

!-------------------------------------------------------------------------------------                                                   
! splits a string to substrings, returns array
!-------------------------------------------------------------------------------------

function SPLIT_STRING (str, separ, dims, word) result(error_status)
  character(*), intent(in) :: str
  character(*), intent(in) :: separ
  integer, intent(in) :: dims
  integer :: error_status
  character(100), dimension(:), allocatable, intent (out) :: word
  integer :: pos1, pos2, i

  error_status = 0
  pos1 = 1
  if (allocated (word) ) deallocate (word)
  allocate (word (dims))
  do i = 1, dims
    pos2 = INDEX(str(pos1:), separ)
    if (pos2 == 0) THEN
       word(i) = str(pos1:)
       EXIT
    endif
    word(i) = str(pos1:pos1+pos2-2)
    pos1 = pos2+pos1
 enddo

end function

!------------------------------------------------------------------------------------- 

function REPLACE_CHAR_IN_STRG (string_inout,target_char,substring_char, &
                    what_del) result(error_status)
 character(*), intent(inout) :: string_inout
 character(*), intent(in) :: target_char, substring_char, what_del
 integer :: indx, error_status

 error_status = 0
 indx = index(string_inout,target_char)
 if (indx > 0) then
    if (trim(what_del) == 'before') then
       string_inout = trim(substring_char)//trim(string_inout(indx+1:))
    elseif (trim(what_del) == 'after') then
       string_inout = trim(string_inout(:indx-1))//trim(substring_char)
    else
       error_status = 1
    endif
 else
    error_status = 1
    return
 endif

endfunction 

!------------------------------------------------------------------------------------- 

!***************************************************************
!** compute general inverse of a matrix using svdcmp from 
!** Numerical Recipes
!** input:
!**        A - a matrix whose general inverse is sought
!**        m,p dimensions of matrix A   
!** output:
!**        A_inv - general inerse of A
!**        ierr - error code from svdcmp
!***************************************************************
subroutine GEN_INVERSE(A,A_inv,m,p,ierr)
 real (kind=ipre), intent(in), dimension(:,:):: A
 integer, intent(in):: m,p
 integer, intent(out):: ierr
 real (kind=ipre), intent(out), dimension(:,:):: A_inv
 integer:: i
 real (kind=ipre), dimension(p,p):: W,V
 real (kind=ipre), dimension(m,p):: U
 real (kind=ipre), dimension(p):: ww
 U = A
 if ( p == 0 ) then
   print *,"p = 0 in gen_inv",m,p,shape(W)
   stop
 endif
 call svdcmp(U,m,p,m,p,ww,V,ierr)
 W=0.0
  do i = 1,p
    W(i,i)=1.0/ww(i)
  enddo 
 A_inv = matmul(V,matmul(W,transpose(U)))
end subroutine GEN_INVERSE
!----------------------------------------------------------------

 subroutine svdcmp(a,m,n,mp,np,w,v,ierr)
   integer, intent(in):: m,mp,n,np
   integer, intent(out):: ierr
!  real (kind=ipre), intent(inout), dimension(:,:):: a(mp,np),v(np,np),w(np)
   real (kind=ipre), intent(inout), dimension(:,:):: a
   real (kind=ipre), intent(out), dimension(:,:):: v
   real (kind=ipre), intent(out), dimension(:):: w
   integer, parameter:: NMAX=500
!CU    USES pythag
   integer:: i,its,j,jj,k,l,nm
   real (kind=ipre):: anorm,c,f,g,h,s,sscale,x,y,z
   real (kind=ipre), dimension(NMAX):: rv1
      ierr = 0
      g=0.0
      sscale=0.0
      anorm=0.0
      twentyfive: do i=1,n
        l=i+1
        rv1(i)=sscale*g
        g=0.0
        s=0.0
        sscale=0.0
        if(i<=m)then
          do  k=i,m
            sscale=sscale+abs(a(k,i))
          enddo
          if(sscale/=0.0)then
            do k=i,m
              a(k,i)=a(k,i)/sscale
              s=s+a(k,i)*a(k,i)
            end do
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do j=l,n
              s=0.0
              do k=i,m
                s=s+a(k,i)*a(k,j)
              end do
              f=s/h
              do k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
              end do
            end do
            do k=i,m
              a(k,i)=sscale*a(k,i)
            end do
          endif
        endif
        w(i)=sscale *g
        g=0.0
        s=0.0
        sscale=0.0
        if((i<=m).and.(i/=n))then
          do k=l,n
            sscale=sscale+abs(a(i,k))
          enddo
          if(sscale/=0.0)then
            do k=l,n
              a(i,k)=a(i,k)/sscale
              s=s+a(i,k)*a(i,k)
            enddo 
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do k=l,n
              rv1(k)=a(i,k)/h
            enddo 
            do j=l,m
              s=0.0
              do k=l,n
                s=s+a(j,k)*a(i,k)
              enddo
              do k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
              enddo
            enddo
            do k=l,n
              a(i,k)=sscale*a(i,k)
            enddo
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
       end do twentyfive
      do i=n,1,-1
        if(i < n)then
          if(g /= 0.0)then
            do j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
            enddo
            do j=l,n
              s=0.0
              do k=l,n
                s=s+a(i,k)*v(k,j)
              enddo
              do k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
              enddo
            enddo
          endif
          do j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
          enddo
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
      enddo
      do i=min(m,n),1,-1
        if (i == 0) then
          print *,"error in svdcmp",m,n,i
        endif
        l=i+1
        g=w(i)
        do j=l,n
          a(i,j)=0.0
        enddo
        if(g /= 0.0)then
          g=1.0/g
          do j=l,n
            s=0.0
            do k=l,m
              s=s+a(k,i)*a(k,j)
            enddo
            f=(s/a(i,i))*g
            do k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
            enddo
          enddo
          do j=i,m
            a(j,i)=a(j,i)*g
          enddo
        else
          do j= i,m
            a(j,i)=0.0
          end do
        endif
        a(i,i)=a(i,i)+1.0
      enddo
      fortynine: do k=n,1,-1
        fortyeight: do its=1,30
          do l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm) == anorm)  then
!              goto 2
               exit
            endif
            if((abs(w(nm))+anorm) == anorm) then
!                goto 1
                 exit
            endif
          enddo
  if((abs(rv1(l))+anorm) /= anorm)  then
!1      c=0.0
        c=0.0
        s=1.0
        do i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm) == anorm) then
!              goto 2
               exit
            endif
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
            enddo
         enddo
       endif
!2         z=w(k)
          z=w(k)
          if (l == k)then
            if(z < 0.0)then
              w(k)=-z
              do j=1,n
                v(j,k)=-v(j,k)
              enddo
            endif
!           goto 3
            exit
          endif
!         if(its.eq.30) pause 'no convergence in svdcmp'
          if(its == 30) then
             ierr = 1
             return
          endif
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0_ipre)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
            enddo
            z=pythag(f,h)
            w(j)=z
            if(z/=0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
            end do
          end do
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
        end do fortyeight
!3       continue
 end do fortynine
 end subroutine svdcmp
!***********************************************************************

 function pythag(a,b) result(pyth)
   real (kind=ipre), intent(in):: a,b
   real (kind=ipre):: pyth
   real (kind=ipre):: absa,absb
     absa=abs(a)
     absb=abs(b)
       if (absa > absb) then
           pyth = absa * sqrt(1.0+(absb/absa)**2)
       else
           if (absb == 0) then
               pyth = 0.0
           else
               pyth = absb*sqrt(1.0 + (absa/absb)**2)
           endif
       endif
 end function pythag

!------------------------------------------------------------------------------------- 
end module NUMERICAL_ROUTINES

