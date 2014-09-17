!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: Compare_Float_Numbers.f90 (src)
!       Compare_Float_Numbers (program)
!
! PURPOSE: Module containing routines to perform equality and relational
!          comparisons on floating point numbers.
!
! DESCRIPTION:
!
! AUTHORS: Paul van Delst, CIMSS/SSEC 01-Apr-2003
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

MODULE Compare_Float_Numbers


  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module usage
  USE CONSTANTS
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: Compare_Float
  PUBLIC :: OPERATOR (.EqualTo.)
  PUBLIC :: OPERATOR (.GreaterThan.)
  PUBLIC :: OPERATOR (.LessThan.)


  ! ---------------------
  ! Procedure overloading
  ! ---------------------

  INTERFACE Compare_Float
    MODULE PROCEDURE Compare_Real_Single
    MODULE PROCEDURE Compare_Real_Double
    MODULE PROCEDURE Compare_Complex_Single
    MODULE PROCEDURE Compare_Complex_Double
  END INTERFACE Compare_Float

  INTERFACE OPERATOR (.EqualTo.)
    MODULE PROCEDURE Is_Equal_To_Single
    MODULE PROCEDURE Is_Equal_To_Double
  END INTERFACE OPERATOR (.EqualTo.)

  INTERFACE OPERATOR (.GreaterThan.)
    MODULE PROCEDURE Is_Greater_Than_Single
    MODULE PROCEDURE Is_Greater_Than_Double
  END INTERFACE OPERATOR (.GreaterThan.)

  INTERFACE OPERATOR (.LessThan.)
    MODULE PROCEDURE Is_Less_Than_Single
    MODULE PROCEDURE Is_Less_Than_Double
  END INTERFACE OPERATOR (.LessThan.)


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Module RCS Id string
  CHARACTER( * ), PRIVATE, PARAMETER :: MODULE_RCS_ID = &
    '$Id$'


CONTAINS


!----------------------------------------------------------------------------------
! NAME:
!       .EqualTo.
!
! PURPOSE:
!       Relational operator to test the equality of REAL operands.
!
! CALLING SEQUENCE:
!       IF ( x .EqualTo. y ) THEN
!         .....
!       END IF
!
! OPERANDS:
!       x, y:        Two congruent floating point data objects to compare.
!                    UNITS:      N/A
!                    TYPE:       REAL(kind=real4)   [ == default real]
!                                  OR
!                                REAL(kind=real8)
!                    DIMENSION:  Scalar, or any allowed rank array.
!
! OPERATOR RESULT:
!       (x .EqualTo. y)    The result is a logical value indicating whether
!                          the operands are equal to within numerical precision
!                          UNITS:      N/A
!                          TYPE:       LOGICAL
!                          DIMENSION:  Same as operands.
!
! PROCEDURE:
!       The test performed is
!
!         ABS( x - y ) < SPACING( MAX(ABS(x),ABS(y)) )
!
!       If the result is .TRUE., the numbers are considered equal.
!
!----------------------------------------------------------------------------------

  ELEMENTAL FUNCTION Is_Equal_To_Single( x, y ) RESULT( Equal_To )
    REAL(kind=real4), INTENT(IN)  :: x, y
    LOGICAL :: Equal_To
    Equal_To = ABS(x-y) < SPACING( MAX(ABS(x),ABS(y)) )
  END FUNCTION Is_Equal_To_Single

  ELEMENTAL FUNCTION Is_Equal_To_Double( x, y ) RESULT( Equal_To )
    REAL(kind=real8), INTENT(IN)  :: x, y
    LOGICAL :: Equal_To
    Equal_To = ABS(x-y) < SPACING( MAX(ABS(x),ABS(y)) )
  END FUNCTION Is_Equal_To_Double


!----------------------------------------------------------------------------------
! NAME:
!       .GreaterThan.
!
! PURPOSE:
!       Relational operator to test if one REAL operand is greater than another.
!
! CALLING SEQUENCE:
!       IF ( x .GreaterThan. y ) THEN
!         .....
!       END IF
!
! OPERANDS:
!       x, y:        Two congruent floating point data objects to compare.
!                    UNITS:      N/A
!                    TYPE:       REAL(kind=real4)   [ == default real]
!                                  OR
!                                REAL(kind=real8)
!                    DIMENSION:  Scalar, or any allowed rank array.
!
! OPERATOR RESULT:
!       (x .GreaterThan. y)    The result is a logical value indicating whether
!                              the operand x is greater than y by more than
!                              the spacing between representable floating point
!                              numbers.
!                              UNITS:      N/A
!                              TYPE:       LOGICAL
!                              DIMENSION:  Same as operands.
!
! PROCEDURE:
!       The test performed is
!
!         ( x - y ) >= SPACING( MAX(ABS(x),ABS(y)) )
!
!       If the result is .TRUE., x is considered greater than y.
!
!----------------------------------------------------------------------------------

  ELEMENTAL FUNCTION Is_Greater_Than_Single( x, y ) RESULT ( Greater_Than )
    REAL(kind=real4), INTENT(IN) :: x, y
    LOGICAL :: Greater_Than
    IF ( (x-y) >= SPACING( MAX( ABS(x), ABS(y) ) ) ) THEN
      Greater_Than = .TRUE.
    ELSE
      Greater_Than = .FALSE.
    END IF
  END FUNCTION Is_Greater_Than_Single


  ELEMENTAL FUNCTION Is_Greater_Than_Double( x, y ) RESULT ( Greater_Than )
    REAL(kind=real8), INTENT(IN) :: x, y
    LOGICAL :: Greater_Than
    IF ( (x-y) >= SPACING( MAX( ABS(x), ABS(y) ) ) ) THEN
      Greater_Than = .TRUE.
    ELSE
      Greater_Than = .FALSE.
    END IF
  END FUNCTION Is_Greater_Than_Double


!----------------------------------------------------------------------------------
! NAME:
!       .LessThan.
!
! PURPOSE:
!       Relational operator to test if one REAL operand is less than another.
!
! CALLING SEQUENCE:
!       IF ( x .LessThan. y ) THEN
!         .....
!       END IF
!
! OPERANDS:
!       x, y:        Two congruent floating point data objects to compare.
!                    UNITS:      N/A
!                    TYPE:       REAL(kind=real4)   [ == default real]
!                                  OR
!                                REAL(kind=real8)
!                    DIMENSION:  Scalar, or any allowed rank array.
!
! OPERATOR RESULT:
!       (x .LessThan. y)    The result is a logical value indicating whether
!                           the operand x is less than y by more than the
!                           spacing between representable floating point
!                           numbers.
!                           UNITS:      N/A
!                           TYPE:       LOGICAL
!                           DIMENSION:  Same as operands.
!
! PROCEDURE:
!       The test performed is
!
!         ( y - x ) >= SPACING( MAX(ABS(x),ABS(y)) )
!
!       If the result is .TRUE., x is considered less than y.
!
!----------------------------------------------------------------------------------

  ELEMENTAL FUNCTION Is_Less_Than_Single( x, y ) RESULT ( Less_Than )
    REAL(kind=real4), INTENT(IN) :: x, y
    LOGICAL :: Less_Than
    IF ( (y-x) >= SPACING( MAX( ABS(x), ABS(y) ) ) ) THEN
      Less_Than = .TRUE.
    ELSE
      Less_Than = .FALSE.
    END IF
  END FUNCTION Is_Less_Than_Single


  ELEMENTAL FUNCTION Is_Less_Than_Double( x, y ) RESULT ( Less_Than )
    REAL(kind=real8), INTENT(IN) :: x, y
    LOGICAL :: Less_Than
    IF ( (y-x) >= SPACING( MAX( ABS(x), ABS(y) ) ) ) THEN
      Less_Than = .TRUE.
    ELSE
      Less_Than = .FALSE.
    END IF
  END FUNCTION Is_Less_Than_Double


!----------------------------------------------------------------------------------
! NAME:
!       Compare_Float
!
! PURPOSE:
!       Function to compare floating point scalars and arrays with adjustible
!       precision tolerance.
!
! CALLING SEQUENCE:
!       Result = Compare_Float( x, y,     &  ! Input
!                               ULP = ULP )  ! Optional input
!
! INPUT ARGUMENTS:
!       x, y:        Two congruent floating point data objects to compare.
!                    UNITS:      N/A
!                    TYPE:       REAL(kind=real4)   [ == default real]
!                                  OR
!                                REAL(kind=real8)
!                                  OR
!                                COMPLEX(kind=real4)
!                                  OR
!                                COMPLEX(kind=real8)
!                    DIMENSION:  Scalar, or any allowed rank array.
!                    ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       ULP:         Unit of data precision. The acronym stands for "unit in
!                    the last place," the smallest possible increment or decrement
!                    that can be made using a machine's floating point arithmetic.
!                    A 0.5 ulp maximum error is the best you could hope for, since
!                    this corresponds to always rounding to the nearest representable
!                    floating-point number. Value must be positive - if a negative
!                    value is supplied, the absolute value is used.
!                    If not specified, the default value is 1.
!                    UNITS:      N/A
!                    TYPE:       INTEGER
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: OPTIONAL, INTENT(IN)
!                  
! FUNCTION RESULT:
!       Result:      The return value is a logical value indicating whether
!                    the inputs are equal (to within the required precision)
!                    .TRUE.  - if the floating point numbers are equal to
!                              within the specified tolerance. 
!                    .FALSE. - if the floating point numbers are different.
!                    UNITS:      N/A
!                    TYPE:       LOGICAL
!                    DIMENSION:  Scalar
!
! PROCEDURE:
!       The test performed is
!
!         ABS( x - y ) < ( ULP * SPACING( MAX(ABS(x),ABS(y)) ) )
!
!       If the result is .TRUE., the numbers are considered equal.
!
!       The intrinsic function SPACING(x) returns the absolute spacing of numbers
!       near the value of x,
!
!                      {     EXPONENT(x)-DIGITS(x)
!                      {  2.0                        for x /= 0
!         SPACING(x) = {
!                      {  
!                      {  TINY(x)                    for x == 0
!
!       The ULP optional argument scales the comparison.
!
!       James Van Buskirk and James Giles suggested this method for floating
!       point comparisons in the comp.lang.fortran newsgroup.
!
!       For complex numbers, the same test is applied to both the real and
!       imaginary parts and each result is ANDed.
!
!----------------------------------------------------------------------------------

  ELEMENTAL FUNCTION Compare_Real_Single( x, y, ulp ) RESULT( Compare )
    REAL(kind=real4),      INTENT(IN)  :: x
    REAL(kind=real4),      INTENT(IN)  :: y
    INTEGER, OPTIONAL, INTENT(IN)  :: ulp
    LOGICAL :: Compare
    REAL(kind=real4) :: Rel
    Rel = 1.0_real4
    IF ( PRESENT( ulp ) ) THEN
      Rel = REAL( ABS(ulp), kind=real4 )
    END IF
    Compare = ABS(x-y) < ( Rel * SPACING( MAX(ABS(x),ABS(y)) ) )
  END FUNCTION Compare_Real_Single


  ELEMENTAL FUNCTION Compare_Real_Double( x, y, ulp ) RESULT( Compare )
    REAL(kind=real8),      INTENT(IN)  :: x
    REAL(kind=real8),      INTENT(IN)  :: y
    INTEGER, OPTIONAL, INTENT(IN)  :: ulp
    LOGICAL :: Compare
    REAL(kind=real8) :: Rel
    Rel = 1.0_real8
    IF ( PRESENT( ulp ) ) THEN
      Rel = REAL( ABS(ulp), kind=real8 )
    END IF
    Compare = ABS( x-y ) < ( Rel * SPACING( MAX(ABS(x),ABS(y)) ) )
  END FUNCTION Compare_Real_Double


  ELEMENTAL FUNCTION Compare_Complex_Single( x, y, ulp ) RESULT( Compare )
    COMPLEX(kind=real4),   INTENT(IN)  :: x
    COMPLEX(kind=real4),   INTENT(IN)  :: y
    INTEGER, OPTIONAL, INTENT(IN)  :: ulp
    LOGICAL :: Compare
    REAL(kind=real4) :: xr, xi
    REAL(kind=real4) :: yr, yi
    xr=REAL(x,kind=real4); xi=AIMAG(x)
    yr=REAL(y,kind=real4); yi=AIMAG(y)
    Compare = Compare_Real_Single(xr,yr,ulp=ulp) .AND. &
              Compare_Real_Single(xi,yi,ulp=ulp)
  END FUNCTION Compare_Complex_Single


  ELEMENTAL FUNCTION Compare_Complex_Double( x, y, ulp ) RESULT( Compare )
    COMPLEX(kind=real8),   INTENT(IN)  :: x
    COMPLEX(kind=real8),   INTENT(IN)  :: y
    INTEGER, OPTIONAL, INTENT(IN)  :: ulp
    LOGICAL :: Compare
    REAL(kind=real8) :: xr, xi
    REAL(kind=real8) :: yr, yi
    xr=REAL(x,kind=real8); xi=AIMAG(x)
    yr=REAL(y,kind=real8); yi=AIMAG(y)
    Compare = Compare_Real_Double(xr,yr,ulp=ulp) .AND. &
              Compare_Real_Double(xi,yi,ulp=ulp)
  END FUNCTION Compare_Complex_Double

END MODULE Compare_Float_Numbers
