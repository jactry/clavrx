!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: sort_module.f90 (src)
!       Sort_Module (program)
!
! PURPOSE: Shell sort with Fortran 90 Array Sections
!
! DESCRIPTION: 
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
!--------------------------------------------------------------------------------------
  module Sort_Module
    implicit none

  contains

    subroutine Sort (Array)
      ! (C) Copyright 1995, Loren P. Meissner
      ! Permission is given to copy or use for any purpose, provided that
      !   these 3 comment lines are included with the source program.
      integer, intent(in out), dimension(:) :: Array
      integer :: N, H, I
 ! start subroutine Sort
      N = Size (Array)
      H = 1
      do
        H = 3 * H + 1
        if (H > N) exit
      end do
      do
        H = H / 3
        do I = 1, H
          call InSort(Array(I: N : H)) ! Array section
        end do
        if (H <= 1) exit
      end do
      return
    end subroutine Sort

    subroutine InSort (Array) ! Straight insertion
      integer, intent(in out), dimension(:) :: Array
      integer :: N, I, J
      real :: Next
  ! start subroutine InSort
      N = Size (Array)
      do I = 1, N - 1
        Next = Array(I + 1)
        do J = I, 1, -1
          if (Array(J) <= Next) exit
          Array(J + 1) = Array(J)
        end do
        Array(J + 1) = Next
      end do
      return
    end subroutine InSort

  end module Sort_Module
