module GOES_IMAGER_ROUTINES
 use CONSTANTS
 use GOES_AREA_IO_MODULE
 implicit none
 public:: LMODEL,LPOINT,COMP_ES,COMP_LP,INST2E,GPOINT,TIME50!,GATT
 private

 contains
!========================================================================
!     I N S T 2 E
!========================================================================
      SUBROUTINE INST2E( R, P, Y, A, AT )
!
!     AUTHOR:         KELLY DEAN
!
!     CREATED:        October 1994
!   
!     DEVELOPED FOR:  CIRA/COLORADO STATE UNIVERSITY
!
!     PURPOSE:        
!	Procedure INST2E accepts the single precision roll, pitch and yaw
!    angles of an instrument and returns the double precision instrument
!    to earth coordinates transformation matrix.
!
!     REVISION:       0.0
!
!     REFERENCES:
!                     OTHER DOCUMENTS
!
!     COMMENTS:
!               Adapted from Igor Levine program for Integral Systems, Inc.
!
!     ARGUMENTS:
!       NAME:   TYPE:       PURPOSE:                        IN/OUT:
!	 R      REAL*8	    Roll angle (rad)		     IN
!	 P      REAL*8      Pitch angle (rad)		     IN
!	 Y      REAL*8	    Yaw angle (rad)		     IN
!	 A      REAL*8	    Spacecraft to ECEF coordinates  OUT
!				 transformation matrix
!	 AT     REAL*8	    Instrument to ECEF coordinates  OUT
!				 transformation matrix
!    
!     VARIABLES:
!       NAME:    PURPOSE:
!       ****************  INTEGER    *****************
!	I    Indices
!	J    Indices
!       ****************  REAL*8     *****************
!	RPY  Instrument to body coordinates transformation matrix
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     VARIABLE DECLARATION SECTION:
      REAL*8	  A(3,3),AT(3,3),R,RPY(3,3),P,Y
      INTEGER*4   I,J
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
!     Compute instrument to body coordinates transformation matrix
!     by using a small angle approximation of  trigonometric function
!     of the roll, pitch and yaw.
!
      RPY(1,1) = 1.0D0 - 0.5D0 * ( P * P + Y * Y )
      RPY(1,2) = -Y
      RPY(1,3) = P
      RPY(2,1) = Y + P * R
      RPY(2,2) = 1.0D0 - 0.5D0 * ( Y * Y + R * R )
      RPY(2,3) = -R
      RPY(3,1) = -P + R * Y
      RPY(3,2) = R + P * Y
      RPY(3,3) = 1.0D0 - 0.5D0 * ( P * P + R * R )
!
!     Multiplication of matrices A and RPY
!
      DO I = 1,3
       DO J = 1,3
        AT(I,J) = A(I,1)*RPY(1,J)+A(I,2)*RPY(2,J)+A(I,3)*RPY(3,J)
       ENDDO
      ENDDO
!
!
      RETURN
      END SUBROUTINE INST2E
!======================================================================
!     T I M E 5 0
!======================================================================
      REAL*8 FUNCTION TIME50(btim)
!
!     AUTHOR:         Garrett Campbell and Kelly Dean
!
!     CREATED:        October 1994
!   
!     DEVELOPED FOR:  CIRA/COLORADO STATE UNIVERSITY
!
!     PURPOSE:        
!	Function TIME50 will take the epoch time from the GVAR NAVstr%and
!      convert it to minutes from January 1, 1950.  NOTE - Epoch time in
!      the NAVstr%is not in the same format as other BCD times.C
!
!     REVISION:       0.0
!
!     ARGUMENTS:
!       NAME:   TYPE:       PURPOSE:                          IN/OUT:
!	btim     BYTE        Binary coded data (BCD) time      IN
!    
!     FUNCTIONS:
!       NAME:   TYPE:       PURPOSE:                        LIBRARY:
!        MOD     INTEGER     Returns a remainder              Intrinsic
!
!       NAME:    PURPOSE:
!       ****************  INTEGER    *****************
!       day_100   Part of day extracted from BCD
!       day_10    Part of day extracted from BCD
!       day_1     Part of day extracted from BCD
!       hour_10   Part of hour extracted from BCD
!       hour_1    Part of hour extracted from BCD
!       min_10    Part of minute extracted from BCD
!       min_1     Part of minute extracted from BCD
!       NY        YEAR
!       ND        DAY OF YEAR
!       NH        HOUR
!       NM        MINUTE
!       ibt
!       J         Loop control variable
!       year_1000 Part of year extracted from BCD
!       year_100  Part of year extracted from BCD
!       year_10   Part of year extracted from BCD
!       year_1    Part of year extracted from BCD
!       ****************  REAL*8     *****************
!       S        SECONDS - Double precision
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     VARIABLE DECLARATION SECTION:
!
      byte btim(8),bt
      INTEGER NY,ND,NH,NM,J
      INTEGER year_1000, year_100, year_10, year_1
      INTEGER day_100, day_10, day_1, hour_10, hour_1, min_10, min_1
      INTEGER sec_10, sec_1, msec_100, msec_10, msec_1
      integer ibt
      REAL flywheel
      REAL*8 S
!
!     Equivalence DECLARATION SECTION:
!
      equivalence (bt,ibt)
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
!    Extract the Binary Coded Time into separate year, Julian day, hour
!    minutes, seconds.
!
!     bt = btim(1)
      day_1   = ibt/16
      hour_10 = mod(ibt,16)
      bt = btim(2)
      day_100  = mod(ibt/16,8)
      day_10   = mod(ibt,16)
      flywheel = mod(ibt,2)
      bt = btim(3)
      year_10   = ibt/16
      year_1   = mod(ibt,16)
      bt = btim(4)
      year_1000 = ibt/16
      year_100  = mod(ibt,16)
      bt = btim(5)
      msec_10 = ibt/16
      msec_1  = mod(ibt,16)
      bt = btim(6)
      sec_1    = ibt/16
      msec_100 = mod(ibt,16)
      bt = btim(7)
      min_1  = ibt/16
      sec_10 = mod(ibt,16)
      bt = btim(8)
      hour_1 = ibt/16
      min_10 = mod(ibt,16)
!
!     Make the year, Julian day, hour, minute, and seconds.
!
      ny = year_1000 * 1000 + year_100 * 100 + year_10 * 10 + year_1
      nd = day_100 * 100 + day_10 * 10 + day_1 
      nh = hour_10 * 10 + hour_1
      nm = min_10 * 10 + min_1
      s  = sec_10 * 10.0D0 + sec_1 +                               &
          msec_100 * 0.1D0 + msec_10 * 0.01D0 + msec_1 * 0.001D0

!
!     HERE WE CONVERT INTEGER YEAR AND DAY OF YEAR TO NUMBER OF                 
!     DAYS FROM 0 HOUR UT, 1950 JAN. 1.0                                        
!     THIS CONVERTION IS BASED ON AN ALGORITHM BY FLIEGEL AND VAN               
!     FLANDERN, COMM. OF ACM, VOL.11, NO. 10, OCT. 1968 (P.657)                 
!
      j = nd + 1461 * (ny + 4799) / 4 - 3 *     &
         ( ( ny + 4899 ) / 100 ) / 4 - 2465022
!
!    Compute time in minutes from January 1.0, 1950 as double precision.
!
      TIME50 = j * 1440.0D0 + nh * 60.0D0 + nm + s / 60.0D0
!
!
      RETURN
      END FUNCTION TIME50
!=======================================================================
!    L M O D E L
!=======================================================================
      SUBROUTINE LMODEL( T, TU, NAVstr, RLAT, RLON )
!
!     AUTHOR:         KELLY DEAN
!
!     CREATED:        October 1994
!   
!     DEVELOPED FOR:  CIRA/COLORADO STATE UNIVERSITY
!
!     PURPOSE:        
!	Procedure LModel accepts an input time and a set of O&A parameters
!	and computes position of the satellite, the attitude angles and
!	attitudes misalignment and the instrument to earth fixed coordinates
!	transformation matrix.
!
!	This procedure computes the position of the satellite and the
!	attitude of the imager or sounder. The calculations are based
!	on the Oats orbit and attitude model represented by the O&A
!	parameter set in NAVstr%
!
!     REVISION:       0.0
!
!     REFERENCES:
!                 Part of this code was adapted from Igor Levine work
!            for Integal System, Inc.
!
!     ARGUMENTS:
!       NAME:   TYPE:       PURPOSE:                                IN/OUT:
!	T        REAL*8      Input time from Jan 1, 1950 (Minutes)   IN
!	TU       REAL*8      Epoch time from Jan 1, 1950 (Minutes)   IN
!	RLAT     REAL*8	     Subsatellite Geodetic latitude (rad)    OUT
!	RLON     REAL*8	     Subsatellite Geodetic Longitude (rad)   OUT
!    
!     SUBROUTINES:
!       NAME:   PURPOSE:                                    LIBRARY:
!        INST2E   Computes instrument to earth coordinates      GVARnav
!
!     FUNCTIONS:
!       NAME:   TYPE:       PURPOSE:                        LIBRARY:
!        DATAN    REAL*8      Arc tangent (double precision)        Intrinsic
!        DATAN2   REAL*8      Arc Tangent (double precision)        Intrinsic
!        DCOS     REAL*8      Cosine (double precision)             Intrinsic
!        DSIN     REAL*8      Sine (double precision)               Intrinsic
!        DTAN     REAL*8      tagent (double precision)             Intrinsic
!        GATT     REAL*8      Compute attitude and misalignment angle  GVARnav
!
!     VARIABLES:
!       NAME:    PURPOSE:
!       ****************  REAL*8       *****************
!       XS     NORMALIZED S/C POSITION IN ECEF COORDINATES              
!       BT    ECEF TO INSTRUMENT COORDINATES TRANSFORMATION            
!       Q3    USED IN SUBROUTINE LPOINT                                
!       PITCH PITCH ANGLES OF INSTRUMENT (RAD)            
!       ROLL  ROLL ANGLES OF INSTRUMENT (RAD)            
!       YAW   YAW ANGLES OF INSTRUMENT (RAD)             
!       PMA  PITCH MISALIGNMENTS OF INSTRUMENT (RAD)         
!       RMA  ROLL MISALIGNMENTS OF INSTRUMENT (RAD)         
!	R    Normalized satellite distance (km)
!	TS   Time from EPOCH (minutes)
!	B    Spacecraft to earth fixed coordinates transmation matrix
!	TE   Exponential time delay from EPOCH (minutes)
!	PHI  Subsatellite geocentric latitude (rad)
!	DR   Radial distance from the nominal (km)
!	PSI  Orbital yaw (rad)
!	LAM  IMC longitude (rad)
!	U    Argument of latitude (rad)
!	SU   DSIN(U)
!	CU   DCOS(U)                                             
!	SI   Sine of the orbit inclination
!	CI   Cosine of the orbit inclination
!	SLAT Sine of geocentric latitude
!	ASC  Longitude of the ascending node (rad)
!	SA   Sine of ASC
!	CA   Cosine of ASC
!	SYAW Sine of the orbit yaw
!	WA   Solar orbit angle (rad)
!	W    Orbit angle (rad)
!	SW   DSIN(W)
!	CW   DCOS(W)
!	S2W  DSIN(2*W)
!	C2W  DCOS(2*W)
!	SW1  DSIN(0.927*W)
!	CW1  DCOS(0.927*W)
!	SW3  Sine of 1.9268*W
!	CW3  Cosine of 1.9268*W
!	DLAT Change in sine of geocentric latitude
!	DYAW Change in sine of orbit yaw
!	A1   Work area
!	A2   Work area
!	XS   S/C position in ECEF coordinates
!
!     COMMON BLOCKS:
!       NAME:     CONTENTS:
!       ELCOMM     Instrument position and attitude variables and 
!                   transformation matrix
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     CONSTANT DECLARATION SECTION:
!
!     include 'kdf.inc'
!
!     VARIABLE DECLARATION SECTION:
!
      INTEGER IMCstatus
      REAL*8 T, TU
      REAL*8 REC(336)
      REAL*8 PI, DEG, RAD, NOMORB, AE, FER, AEBE2, AEBE3, AEBE4
      REAL*8 RLAT, RLON, R, TS, TE, PHI, DR, PSI, LAM, U, SU, CU
      REAL*8 SI, CI, SLAT, ASC, SA, CA, SYAW, WA, W, SW, CW, S2W, C2W
      REAL*8 SW1, CW1, SW3, CW3, DLAT, DYAW, A1, A2
      REAL*8 B(3,3), BT(3,3), XS(3)
      REAL*8 Q3, PITCH, ROLL, YAW
      REAL*8 PMA, RMA
      type(GVAR_NAV):: NAVstr
!
!     FUNCTION DECLARATION SECTION:
!
      REAL*8 DATAN,DATAN2,DCOS,DSIN,DTAN!,GATT
!
!     COMMON BLOCKS:
!
      COMMON /ELCOMM/ xs, bt, q3, pitch, roll, yaw, pma, rma 
!
!     INITIALIZATIONS: (Description mathematical and earth-related constants)
!
!      PI     = 3.141592653589793D0
      PI     = 4.0*DATAN(1.0D0)
      DEG    = 180D0 / PI
      RAD    = PI / 180D0  	! Degrees to radians conversion (PI/180)
      NOMORB = 42164.365D0 	!  Nominal radial distance of satellite (km)
      AE     = 6378.137D0  	!  Earth equatorial radius (km)
      FER    = 1.0D0 - ( 6356.7533D0 / AE )  ! Earth flattening coefficient 
      AEBE2  = 1.0D0 / (1.0D0 - FER )**2 
      AEBE3  = AEBE2 - 1. 
      AEBE4  = ( 1.0D0 - FER )**4-1.
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
!     Determine the IMC status
! 
!      IMCstatus = IBITS(NAVstr%stat,6,1)  
!
!     Assign referenec values to the subsatellite longitude and
!     latitude, the radial distance and the orbit yaw.
!                                                                               
      LAM = NAVstr%ref_long * 1.0D-7
      DR  = NAVstr%ref_rad_dist
      PHI = NAVstr%ref_lat
      PSI = NAVstr%ref_orb_yaw
!
!     Assign reference values to the attitudes and misalignments
!
      ROLL  = NAVstr%ref_att_roll
      PITCH = NAVstr%ref_att_pitch
      YAW   = NAVstr%ref_att_yaw
      RMA   = 0.0
      PMA   = 0.0
!
!     IF IMC_active is OFF, compute changes in the satellite orbit
!
      IF ( IMCstatus .NE. 0 ) THEN
!	PRINT *, ' IMC turned off............'
!
!       Compute time since EPOCH (minutes)
!
        TS = T - TU
!
!       Compute orbite angle and the related trigonometric functions.
!       earth rotational rate (.729115E-4 rad/sec).
!
        W   = 0.729115e-4 * 60.0D0 * TS
        SW  = DSIN(W)
        CW  = DCOS(W)
        SW1 = DSIN(0.927D0*W)
        CW1 = DCOS(0.927D0*W)
        S2W = DSIN(2.0D0*W)
        C2W = DCOS(2.0D0*W)
        SW3 = DSIN(1.9268D0*W)
        CW3 = DCOS(1.9268D0*W)
!
!     Computes change in the IMC_active longitude from the reference.
!
        LAM = LAM + ( NAVstr%ref_long_change(1) * 1.0D-7 ) +             &
      	    ( ( NAVstr%ref_long_change(2) * 1.0D-7 ) +             &
      	      ( NAVstr%ref_long_change(3) * 1.0D-7 ) * W ) * W +             &
              ( NAVstr%ref_long_change(10) * 1.0D-7 ) * SW1 +             &
      	      ( NAVstr%ref_long_change(11) * 1.0D-7 ) * CW1 +             &
      	    ( ( NAVstr%ref_long_change(4) * 1.0D-7 )  * SW  +             &
      	      ( NAVstr%ref_long_change(5) * 1.0D-7 ) * CW +             &
              ( NAVstr%ref_long_change(6) * 1.0D-7 ) * S2W +             &
      	      ( NAVstr%ref_long_change(7) * 1.0D-7 ) * C2W +             &
      	      ( NAVstr%ref_long_change(8) * 1.0D-7 ) * SW3 +             &
      	      ( NAVstr%ref_long_change(9) * 1.0D-7 ) * CW3 +             &
              W * ( ( NAVstr%ref_long_change(12) * 1.0D-7 ) * SW +             &
      	      ( NAVstr%ref_long_change(13) * 1.0D-7 ) * CW ) ) * 2.0D0
!
!       Computes change in radial distance from the reference (km)
!
        DR = DR + ( NAVstr%ref_rad_dist_change(1) * 1.0D-7 ) +             &
      		  ( NAVstr%ref_rad_dist_change(2) * 1.0D-7 ) * CW  +             &
      		  ( NAVstr%ref_rad_dist_change(3) * 1.0D-7 ) * SW  +             &
           	  ( NAVstr%ref_rad_dist_change(4) * 1.0D-7 ) * C2W +             &
      		  ( NAVstr%ref_rad_dist_change(5) * 1.0D-7 ) * S2W +             &
      		  ( NAVstr%ref_rad_dist_change(6) * 1.0D-7 ) * CW3 +             &
      		  ( NAVstr%ref_rad_dist_change(7) * 1.0D-7 ) * SW3 +             &
      		  ( NAVstr%ref_rad_dist_change(8) * 1.0D-7 ) * CW1 +             &
      		  ( NAVstr%ref_rad_dist_change(9) * 1.0D-7 ) * SW1 +             &
      		  W * ( ( NAVstr%ref_rad_dist_change(10) * 1.0D-7 ) * CW +             &
      		  ( NAVstr%ref_rad_dist_change(11) * 1.0D-7 ) * SW )
!
!       Computes the sine of the change in the geocentric latitude.
!
        DLAT = ( NAVstr%sine_lat(1) * 1.0D-7 ) +             &
      	       ( NAVstr%sine_lat(2) * 1.0D-7 ) * CW  +             &
      	       ( NAVstr%sine_lat(3) * 1.0D-7 ) * SW  +             &
      	       ( NAVstr%sine_lat(4) * 1.0D-7 ) * C2W +             &
      	       ( NAVstr%sine_lat(5) * 1.0D-7 ) * S2W +             &
      	       W * ( ( NAVstr%sine_lat(6) * 1.0D-7 ) * CW +             &
      	       ( NAVstr%sine_lat(7) * 1.0D-7 ) * SW ) +             &
      	       ( NAVstr%sine_lat(8) * 1.0D-7 ) * CW1 +             &
      	       ( NAVstr%sine_lat(9) * 1.0D-7 ) * SW1
!
!	Computes geocentric latitude by using an expansion for arcsine.
!
        PHI = PHI + DLAT * ( 1.0D0 + DLAT * DLAT / 6.0D0 )
!
!	Computes sine of the change in the orbit yaw.
!
        DYAW = ( NAVstr%sine_orb_yaw(1) * 1.0D-7 ) +             &
      	       ( NAVstr%sine_orb_yaw(2) * 1.0D-7 ) * SW  +             &
      	       ( NAVstr%sine_orb_yaw(3) * 1.0D-7 ) * CW  +             &
      	       ( NAVstr%sine_orb_yaw(4) * 1.0D-7 ) * S2W +             &
      	       ( NAVstr%sine_orb_yaw(5) * 1.0D-7 ) * C2W +             &
      	       W * ( ( NAVstr%sine_orb_yaw(6) * 1.0D-7 ) * SW +             &
      	       ( NAVstr%sine_orb_yaw(7) * 1.0D-7 ) * CW ) +             &
      	       ( NAVstr%sine_orb_yaw(8) * 1.0D-7 ) * SW1 +             &
      	       ( NAVstr%sine_orb_yaw(9) * 1.0D-7 ) * CW1
!
!	Computes the orbit yaw by using an expansion for arcsine.
!
        PSI = PSI + DYAW * ( 1.0D0 + DYAW * DYAW / 6.0D0 )
!      ELSE
!C         WRITE(6,*) ' IMC is turned on .......... >',IMCstatus
      ENDIF
!
!     Conversion of the IMC_active longitude and orbit yaw to the subsatellite
!     longitude and the orbit inclination (REF: GOES-PCC-TM-2473). Inputs
!     required for earth location and gridding 
!
      SLAT = DSIN(PHI)
      SYAW = DSIN(PSI)

      SI = SLAT**2 + SYAW**2
      CI = DSQRT(1.0D0 - SI )
      SI = DSQRT(SI)
      IF ( SYAW .NE. 0.0D0 ) THEN
        U = DATAN2(SLAT,SYAW)
      ELSE IF (SLAT .GT. 0.0D0 ) THEN
        U = 1.570796D0
      ELSE IF (SLAT .LT. 0.0D0 ) THEN
        U = 4.712389D0
      ELSE 
        U = LAM
      ENDIF
!
      SU = DSIN(U)
      CU = DCOS(U)
!
!     Computes longitude of the ascending node.
!
      ASC = LAM - U
      SA = DSIN(ASC)
      CA = DCOS(ASC)
!
!     Computes the subsatellite geographic latitude (rad)
!
      RLAT = DATAN(AEBE2*DTAN(PHI))
!
!     Computes the subsatellite geographic longitude (rad)
!
      RLON = ASC + DATAN2(CI*SU,CU)
!
!     Computes the spacecraft to earth fixed coordinates transformation matrix.
!
!         (VECTOR IN ECEF COORDINATES) = B * (VECTOR IN S/C COORDINATES)
!
      B(1,2) = -SA * SI
      B(2,2) = CA * SI
      B(3,2) = -CI
      B(1,3) = -CA * CU + SA * SU * CI
      B(2,3) = -SA * CU - CA * SU * CI
      B(3,3) = -SLAT
      B(1,1) = -CA * SU - SA * CU * CI
      B(2,1) = -SA * SU + CA * CU * CI
      B(3,1) = CU * SI
!
!     Computes the normalized spacecraft position vector in earth fixed
!     coordinates - XS.
!
      R = (NOMORB + DR) / AE
      XS(1) = -B(1,3) * R
      XS(2) = -B(2,3) * R
      XS(3) = -B(3,3) * R
!
!     Precomputes Q3 ( Used in LPoint ).
!
      Q3 = XS(1)**2 + XS(2)**2 + AEBE2 * XS(3)**2 - 1.0D0
!
!     Computes the attitudes and misalignments IF IMC_active is OFF
!
      IF ( IMCstatus .NE. 0 ) THEN
!	PRINT *, ' IMC turned off............'
!
!     Computes the solar orbit angle
!                                                                               
         WA = ( NAVstr%solar_rate * 1.0D-7 ) * TS
!
!     Computes the difference between current time, TS, and the
!     exponential time. Note that both times are since EPOCH.
!
         TE = TS - ( NAVstr%exp_start_time * 1.0D-7 )
!
!     Computes ROLL + ROLL Misalignment
!
         ROLL = ROLL + GATT(NAVstr%roll_att,WA,TE)
!
!     Computes Pitch + Pitch Misalignment
!                                                                               
         PITCH = PITCH + GATT(NAVstr%pitch_att,WA,TE)
!
!     Computes YAW
!
         YAW = YAW + GATT(NAVstr%yaw_att,WA,TE)
!
!     Computes roll misalignment
!                                                                               
         RMA = GATT(NAVstr%roll_misalgn,WA,TE)
!
!     Computes pitch misalignment
!
         PMA = GATT(NAVstr%pitch_misalgn,WA,TE)
!
!     Apply the Earth Sensor compensation IF needed.
!
         IF ( TS .GE. ( NAVstr%start_time * 1.0D-2 ) ) THEN
                ROLL = ROLL + NAVstr%IMC_corr_roll * 1.0D-7
                PITCH = PITCH + NAVstr%IMC_corr_pitch * 1.0D-7
                YAW = YAW + NAVstr%IMC_corr_yaw * 1.0D-7
         ENDIF
!      ELSE
!C         WRITE(6,*) ' IMC is turned on .......... >',IMCstatus
      ENDIF
!
!     Computes the instrument to earth fixed coordinates transformation
!     matrix - BT
!
      CALL INST2E( ROLL, PITCH, YAW, B, BT)
!
      RETURN

      CONTAINS
 
!=======================================================================
!     G A T T
!=======================================================================
      REAL*8 FUNCTION GATT( IMGR_REP, WA, TE)
!
!     AUTHOR:         KELLY DEAN
!
!     CREATED:        October 1994
!   
!     DEVELOPED FOR:  CIRA/COLORADO STATE UNIVERSITY
!
!     PURPOSE:        
!	   This function computes an attitude/misalignment angle from
!	 a given subset of the O&A parameters.
!
!     REVISION:       0.0
!
!     ARGUMENTS:
!       NAME:   TYPE:       PURPOSE:                        IN/OUT:
!        IMGR_REP  STRUCTURE
!        TE        REAL*8     Input exponential time          IN 
!                             delay from epoch (minutes)
!	 WA	   REAL*8     Input solar orbit angle (rad)   IN
!    
!     FUNCTIONS:
!       NAME:   TYPE:       PURPOSE:                        LIBRARY:
!        DCOS    REAL*8      Cosine ( Double precision )      INTRINSIC
!
!     VARIABLES:
!       NAME:    PURPOSE:
!       ****************  INTEGER    *****************
!       l      Loop control variable
!       m      Temporary variable for order of monomial sinusoids
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     INCLUDE DECLARATION SECTION:
!
!     include 'kdf.inc'
!
!     VARIABLE DECLARATION SECTION:
!
      INTEGER l,m
      REAL*8  TE, WA
      type (IMGR_RP) :: IMGR_REP
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
      GATT = IMGR_REP%mean_att_ang_const * 1.0D-7
!
!	Computes the exponential term.
!
      IF ( TE .GE. 0.0D0 ) THEN 
       GATT = GATT + (IMGR_REP%exp_mag * 1.0D-7 ) *     &
            EXP(-te / ( IMGR_REP%exp_time_const * 1.0D-2 ))
      ENDIF
!
!	Calculation of sinusoids.
!
      DO l = 1, IMGR_REP%num_sinu_per_angle
         GATT = GATT + ( IMGR_REP%sinusoid(l)%mag_sinu * 1.0D-7 ) *    &
               DCOS(wa * l +                                           &
               ( IMGR_REP%sinusoid(l)%phase_ang_sinu * 1.0D-7 ) )
      ENDDO
!
!	Computes monomial sinusoids.
!
      DO l = 1, IMGR_REP%num_mono_sinu
          m = IMGR_REP%monomial(l)%order_mono_sinu
          GATT = GATT + (IMGR_REP%monomial(l)%mag_mono_sinu * 1.0D-7) *    &
            (wa-(IMGR_REP%monomial(l)%ang_from_epoch * 1.0D-7) )**m *    &
            DCOS( IMGR_REP%monomial(l)%order_appl_sinu * wa +     &
            ( IMGR_REP%monomial(l)%phase_ang_sinu * 1.0D-7 ) )
       ENDDO
!
!
       RETURN
       END FUNCTION GATT

      END SUBROUTINE LMODEL
!***********************************************************************        
!***********************************************************************        
!**                                                                             
!**   INTEGRAL SYSTEMS, INC.                                                    
!**                                                                             
!***********************************************************************        
!**                                                                             
!**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
!**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
!**   ROUTINE   : LPOINT                                                        
!**   SOURCE    : F.LPOINT                                                      
!**   LOAD NAME : ANY                                                           
!**   PROGRAMMER: IGOR LEVINE                                                   
!**                                                                             
!**   VER.    DATA    BY   COMMENT                                              
!**   ----  --------  ---  ---------------------------------------------        
!**   A     01/09/89  IL   INITIAL CREATION                                     
!**   A     06/02/89  IL   COORDINATE AXES CHANGED ACCORDING TO                 
!**                        FORD'S DEFINITION IN SDAIP, DRL504-01                
!**                                                                             
!***********************************************************************        
!**                                                                             
!**   THIS SUBROUTINE CONVERTS THE INSTRUMENT ELEVATION AND SCAN                
!**   ANGLES TO THE RELATED GEOGRAPHIC LATITUDE AND LONGITUDE.                  
!**                                                                             
!***********************************************************************        
!**                                                                             
!**   CALLED BY       : ANY                                                     
!**   COMMONS MODIFIED: NONE                                                    
!**   INPUTS          : NONE                                                    
!**   OUTPUTS         : NONE                                                    
!**   ROUTINES CALLED : NONE                                                    
!**                                                                             
!***********************************************************************        
!***********************************************************************        
      SUBROUTINE LPOINT(ALPHA,ZETA,RLAT,RLON,IERR)                              
!                                                                               
!     CALLING PARAMETERS                                                        
!                                                                               
      REAL*8   ALPHA                                                            
!                             ELEVATION ANGLE (RAD)                             
      REAL*8   ZETA                                                             
!                             SCAN ANGLE (RAD)                                  
      REAL*8   RLAT                                                             
!                             LATITUDE IN RADIANS (OUTPUT)                      
      REAL*8   RLON                                                             
!                             LONGITUDE IN RADIANS (OUTPUT)                     
      INTEGER IERR                                                              
!                             OUTPUT STATUS; 0 - POINT ON THE EARTH             
!                             FOUND, 1 - INSTRUMENT POINTS OFF EARTH            
!                                                                               
!     LOCAL VARIABLES                                                           
!                                                                               
      REAL*8 G1(3)                                                              
!                          POINTING VECTOR IN EARTH CENTERED COORDINATES        
      REAL*8 H                                                                  
!                          SLANT DISTANCE TO THE EARTH POINT (KM)               
      REAL*8 Q1,Q2,D                                                            
!                          WORK SPACE                                           
      REAL*8 G(3)                                                               
!                          POINTING VECTOR IN INSTRUMENT COORDINATES            
      REAL*8 U(3)                                                               
!                          COORDINATES OF THE EARTH POINT (KM)                  
      REAL*8 SA,CA,DA,DZ,D1,CZ                                                  
!                                     WORK SPACE                                
!                                                                               
!     INCLUDE FILES                                                             
!                                                                               
      REAL*8 PI                                                                 
           PARAMETER (PI=3.141592653589793D0)                                   
      REAL*8 DEG                                                                
           PARAMETER (DEG=180.D0/PI)                                            
      REAL*8 RAD                                                                
           PARAMETER (RAD=PI/180.D0)                                            
!                    DEGREES TO RADIANS CONVERSION PI/180                       
      REAL*8 NOMORB                                                             
           PARAMETER (NOMORB=42164.365D0)                                       
!                    NOMINAL RADIAL DISTANCE OF SATELLITE (km)                  
      REAL*8 AE                                                                 
           PARAMETER (AE=6378.137D0)                                            
!                    EARTH EQUATORIAL RADIUS (km)                               
      REAL*8 FER                                                                
           PARAMETER (FER=1.D0-(6356.7533D0/AE))                                
!                    EARTH FLATTENING COEFFICIENT = 1-(BE/AE)                   
      REAL*4 AEBE2                                                              
           PARAMETER (AEBE2=1.D0/(1.D0-FER)**2)                                 
      REAL*4 AEBE3                                                              
           PARAMETER (AEBE3=AEBE2-1.)                                           
      REAL*4 AEBE4                                                              
           PARAMETER (AEBE4=(1.D0-FER)**4-1.)
      REAL*8 XS(3)                                                              
!                      NORMALIZED S/C POSITION IN ECEF COORDINATES              
      REAL*8 BT(3,3)                                                            
!                      ECEF TO INSTRUMENT COORDINATES TRANSFORMATION            
      REAL*8  Q3                                                                
!                      USED IN SUBROUTINE LPOINT                                
      REAL*8 PITCH,ROLL,YAW                                                     
!                          PITCH,ROLL,YAW ANGLES OF INSTRUMENT (RAD)            
      REAL*8 PMA,RMA                                                            
!                          PITCH,ROLL MISALIGNMENTS OF INSTRUMENT (RAD)         
         COMMON /ELCOMM/ XS,BT,Q3,PITCH,ROLL,YAW,PMA,RMA
!***********************************************************************        
      IERR=1                                                                    
!                                                                               
!     COMPUTES TRIGONOMETRIC FUNCTIONS OF THE SCAN AND ELEVATION                
!     ANGLES CORRECTED FOR THE ROLL AND PITCH MISALIGNMENTS                     
!                                                                               
      CA=DCOS(ALPHA)
      SA=DSIN(ALPHA)
      DA=ALPHA-PMA*SA*(1.0D0+DTAN(ZETA))-RMA*(1.0D0-CA)
      DZ=ZETA+RMA*SA                                                            
!                              CORRECTED SCAN ANGLE                             
      CZ=DCOS(DZ)
!                                                                               
!     COMPUTES POINTING VECTOR IN INSTRUMENT COORDINATES                        
!                                                                               
      G(1)=DSIN(DZ)
      G(2)=-CZ*DSIN(DA)
      G(3)=CZ*DCOS(DA)
!                                                                               
!     TRANSFORMS THE POINTING VECTOR TO EARTH FIXED COORDINATES                 
!                                                                               
      G1(1)=BT(1,1)*G(1)+BT(1,2)*G(2)+BT(1,3)*G(3)                              
      G1(2)=BT(2,1)*G(1)+BT(2,2)*G(2)+BT(2,3)*G(3)                              
      G1(3)=BT(3,1)*G(1)+BT(3,2)*G(2)+BT(3,3)*G(3)                              
!                                                                               
!     COMPUTES COEFFICIENTS AND SOLVES A QUADRATIC EQUATION TO                  
!     FIND THE INTERSECT OF THE POINTING VECTOR WITH THE EARTH                  
!     SURFACE                                                                   
!                                                                               
      Q1=G1(1)**2+G1(2)**2+AEBE2*G1(3)**2                                       
      Q2=XS(1)*G1(1)+XS(2)*G1(2)+AEBE2*XS(3)*G1(3)                              
      D=Q2*Q2-Q1*Q3                                                             
      IF (DABS(D).LT.1.D-9) D=0.0D0
!                                                                               
!     IF THE DISCIMINANTE OF THE EQUATION, D, IS NEGATIVE, THE                  
!     INSTRUMENT POINTS OFF THE EARTH                                           
!                                                                               
      IF (D.LT.0.0D0) THEN
         RLAT=999999.0D0
         RLON=999999.0D0
         RETURN                                                                 
      END IF                                                                    
      D=DSQRT(D)                                                                
!                                                                               
!     SLANT DISTANCE FROM THE SATELLITE TO THE EARTH POINT                      
!                                                                               
      H=-(Q2+D)/Q1                                                              
!                                                                               
!     CARTESIAN COORDINATES OF THE EARTH POINT                                  
!                                                                               
      U(1)=XS(1)+H*G1(1)                                                        
      U(2)=XS(2)+H*G1(2)                                                        
      U(3)=XS(3)+H*G1(3)                                                        
!                                                                               
!     SINUS OF GEOCENTRIC LATITUDE                                              
!                                                                               
      D1=U(3)/DSQRT(U(1)**2+U(2)**2+U(3)**2)
!                                                                               
!     GEOGRAPHIC (GEODETIC) COORDINATES OF THE POINT                            
!                                                                               
      RLAT=DATAN(AEBE2*D1/DSQRT(1.0D0-D1*D1))
      RLON=DATAN2(U(2),U(1))
      IERR=0                                                                    
      RETURN                                                                    
      END  SUBROUTINE LPOINT                                                                   
!***********************************************************************        
!***********************************************************************        
!**                                                                             
!**   INTEGRAL SYSTEMS, INC.                                                    
!**                                                                             
!***********************************************************************        
!**                                                                             
!**   PROJECT   : OPERATIONS GROUND EQUIPMENT FOR GOES-NEXT                     
!**   SYSTEM    : EARTH LOCATION USERS GUIDE                                    
!**   ROUTINE   : GPOINT                                                        
!**   SOURCE    : F.GPOINT                                                      
!**   LOAD NAME : ANY                                                           
!**   PROGRAMMER: IGOR LEVINE                                                   
!**                                                                             
!**   VER.    DATA    BY   COMMENT                                              
!**   ----  --------  ---  ---------------------------------------------        
!**   A     12/10/87  IL   INITIAL CREATION                                     
!**   A     06/10/88  IL   REPLACED ASIN WITH ATAN TO SAVE TIME                 
!**   A     06/02/89  IL   COORDINATE AXES CHANGED ACCORDING TO                 
!**                        FORD'S DEFINITION IN SDAIP, DRL 504-01               
!**                                                                             
!***********************************************************************        
!**                                                                             
!**   THIS SUBROUTINE CONVERTS GEOGRAPHIC LATITUDE AND LONGITUDE                
!**   TO THE RELATED ELEVATION AND SCAN ANGLES.                                 
!**                                                                             
!***********************************************************************        
!**                                                                             
!**   CALLED BY       : ANY                                                     
!**   COMMONS MODIFIED: NONE                                                    
!**   INPUTS          : NONE                                                    
!**   OUTPUTS         : NONE                                                    
!**   ROUTINES CALLED : NONE                                                    
!**                                                                             
!***********************************************************************        
!***********************************************************************        
      SUBROUTINE GPOINT(RLAT,RLON,ALF,GAM,IERR)                                 
!                                                                               
!     CALLING PARAMETERS                                                        
!                                                                               
      REAL*8   RLAT	! GEOGRAPHIC LATITUDE IN RADIANS (INPUT)            
      REAL*8   RLON	! GEOGRAPHIC LONGITUDE IN RADIANS (INPUT)           
      REAL*8   ALF	! ELEVATION ANGLE IN RADIANS (OUTPUT)               
      REAL*8   GAM	! SCAN ANGLE IN RADIANS (OUTPUT)                    
      INTEGER IERR	! OUTPUT STATUS
!                             0 - SUCCESSFUL COMPLETION,         
!                             1 - POINT WITH GIVEN LAT/LON IS INVISIBLE         
!
!     LOCAL VARIABLES                                                           
!
      REAL*8 F(3)	! POINTING VECTOR IN EARTH CENTERED COORDINATES
      REAL*8 FT(3)	! POINTING VECTOR IN INSTRUMENT COORDINATES
      REAL*8 U(3)	! COORDINATES OF THE EARTH POINT (KM)
      REAL*8 SING,SLAT,W1,W2  ! WORK SPACE
!                                                                               
!     INCLUDE FILES                                                             
!                                                                               
      REAL*8 PI                                                                 
           PARAMETER (PI=3.141592653589793D0)                                   
      REAL*8 DEG                                                                
           PARAMETER (DEG=180.D0/PI)                                            
      REAL*8 RAD                                                                
           PARAMETER (RAD=PI/180.D0)                                            
!                    DEGREES TO RADIANS CONVERSION PI/180                       
      REAL*8 NOMORB                                                             
           PARAMETER (NOMORB=42164.365D0)                                       
!                    NOMINAL RADIAL DISTANCE OF SATELLITE (km)                  
      REAL*8 AE                                                                 
           PARAMETER (AE=6378.137D0)                                            
!                    EARTH EQUATORIAL RADIUS (km)                               
      REAL*8 FER                                                                
           PARAMETER (FER=1.D0-(6356.7533D0/AE))                                
!                    EARTH FLATTENING COEFFICIENT = 1-(BE/AE)                   
      REAL*4 AEBE2                                                              
           PARAMETER (AEBE2=1.D0/(1.D0-FER)**2)                                 
      REAL*4 AEBE3                                                              
           PARAMETER (AEBE3=AEBE2-1.)                                           
      REAL*4 AEBE4                                                              
           PARAMETER (AEBE4=(1.D0-FER)**4-1.)
      REAL*8 XS(3)                                                              
!                      NORMALIZED S/C POSITION IN ECEF COORDINATES              
      REAL*8 BT(3,3)                                                            
!                      ECEF TO INSTRUMENT COORDINATES TRANSFORMATION            
      REAL*8  Q3                                                                
!                      USED IN SUBROUTINE LPOINT                                
      REAL*8 PITCH,ROLL,YAW                                                     
!                          PITCH,ROLL,YAW ANGLES OF INSTRUMENT (RAD)            
      REAL*8 PMA,RMA                                                            
!                          PITCH,ROLL MISALIGNMENTS OF INSTRUMENT (RAD)         
         COMMON /ELCOMM/ XS,BT,Q3,PITCH,ROLL,YAW,PMA,RMA
!***********************************************************************        
!                                                                               
!     COMPUTES SINUS OF GEOGRAPHIC (GEODETIC) LATITUDE                          
!                                                                               
      SING=DSIN(RLAT)
      W1=AEBE4*SING*SING                                                        
!                                                                               
!     SINUS OF THE GEOCENTRIC LATITUDE                                          
!                                                                               
      SLAT=((0.375D0*W1-0.5D0)*W1+1.0D0)*SING/AEBE2
!                                                                               
!     COMPUTES LOCAL EARTH RADIUS AT SPECIFIED POINT                            
!                                                                               
      W2=SLAT*SLAT                                                              
      W1=AEBE3*W2                                                               
      W1=(0.375D0*W1-0.5D0)*W1+1.D0
!                                                                               
!     COMPUTES CARTESIAN COORDINATES OF THE POINT                               
!                                                                               
      U(3)=SLAT*W1                                                              
      W2=W1*DSQRT(1.0D0-W2)
      U(1)=W2*DCOS(RLON)
      U(2)=W2*DSIN(RLON)
!                                                                               
!     POINTING VECTOR FROM SATELLITE TO THE EARTH POINT                         
!                                                                               
      F(1)=U(1)-XS(1)                                                           
      F(2)=U(2)-XS(2)                                                           
      F(3)=U(3)-XS(3)                                                           
      W2=U(1)*SNGL(F(1))+U(2)*SNGL(F(2))+           & 
         U(3)*SNGL(F(3))*AEBE2                                                  
!                                                                               
!     VERIFIES VISIBILITY OF THE POINT                                          
!                                                                               
      IF (W2.GT.0.0D0) THEN
!                               INVISIBLE POINT ON THE EARTH                    
                   IERR=1                                                       
                   ALF=99999.0D0
                   GAM=99999.0D0
                   RETURN                                                       
       END IF                                                                   
!                                                                               
!     CONVERTS POINTING VECTOR TO INSTRUMENT COORDINATES                        
!                                                                               
      FT(1)=BT(1,1)*F(1)+BT(2,1)*F(2)+BT(3,1)*F(3)                              
      FT(2)=BT(1,2)*F(1)+BT(2,2)*F(2)+BT(3,2)*F(3)                              
      FT(3)=BT(1,3)*F(1)+BT(2,3)*F(2)+BT(3,3)*F(3)                              
!                                                                               
!     CONVERTS POINTING VECTOR TO SCAN AND ELEVATION ANGLES AND                 
!     CORRECTS FOR THE ROLL AND PITCH MISALIGNMENTS                             
!                                                                               
      GAM=ATAN(FT(1)/SQRT(FT(2)**2+FT(3)**2))                                   
      ALF=-DATAN(FT(2)/FT(3))
      W1=DSIN(ALF)
      W2=DCOS(ALF)
      ALF=ALF+RMA*(1.0D0-W2)+PMA*W1*(1.0D0+DTAN(GAM))
      GAM=GAM-RMA*W1                                                            
      IERR=0                                                                    
      RETURN                                                                    
      END  SUBROUTINE GPOINT                                                                      
!=======================================================================
!     C O M P _ E S
!=======================================================================
      SUBROUTINE COMP_ES( NAVstr, line, pixel, elev, scan )
!
!     AUTHOR:         KELLY DEAN
!
!     CREATED:        January 1995
!   
!     DEVELOPED FOR:  CIRA/COLORADO STATE UNIVERSITY
!
!     PURPOSE:        
!		Compute the elevation and scan angles related to the
!	satellite line and pixel numbers.
!
!     REVISION:       1.0
!
!     ARGUMENTS:
!       NAME:   TYPE:       PURPOSE:                        IN/OUT:
!	 NAVstr   Structure   Navigation information         IN
! 	 line     REAL*8      Satellite line number          IN
! 	 pixel    REAL*8      Satellite pixel number         IN
! 	 elev     REAL*8      Elevation angle (rad)          OUT
! 	 scan     REAL*8      Scan angle (rad)               OUT
!    
!     CONSTANTS:
!       NAME:    PURPOSE:
!       ****************  REAL*8     *****************
!	elvln   Elevation angle per detector line (rad)
!	elvmax  Bounds in elevation
!	scnmax  Bounds in scan angle
!	scnpx   Scan angle per pixel (rad)
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     INCLUDE DECLARATION SECTION:
!
!     include 'kdf.inc'
!
!     CONSTANT DECLARATION SECTION:
!
      INTEGER incmax(2) /6136,2805/
      REAL*8  elvmax(2) / 0.2208960D0, 0.22089375D0/
      REAL*8  elvln(2)  /28.0D-6,    280.0D-6/
      REAL*8  elvinc(2) / 8.0D-6,     17.5D-6/
      REAL*8  scnmax(2) / 0.245440D0,  0.2454375D0 /
      REAL*8  scnpx(2)  /16.0D-6,    280.0D-6/
      REAL*8  scninc(2) /16.0D-6,     35.0D-6/
!
!     VARIABLE DECLARATION SECTION:
!
      type(GVAR_NAV):: NAVstr
      REAL*8 line, pixel, elev, scan
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
      IF ( NAVstr%instr .EQ. 1 ) THEN 
!	Recompute elevation and scan biases based on user inputs of
!	cycles and increments obtained from GVAR.
       elvmax(NAVstr%instr) = ( NAVstr%ns_cyl *      &
                               incmax(NAVstr%instr) +      &
                               NAVstr%ns_inc ) *      &
                               elvinc(NAVSTR%instr)
       scnmax(NAVstr%instr) = ( NAVstr%ew_cyl *     &
                                incmax(NAVstr%instr) +     & 
                                NAVstr%ew_inc ) *      &
                                scninc(NAVstr%instr)
!      Compute elevation angle (rad)
       elev = elvmax(NAVstr%instr) + (4.50 - line) * elvln(NAVstr%instr) 
!      Compute scan angle (rad)
       scan = (pixel - 1.0) * scnpx(NAVstr%instr) - scnmax(NAVstr%instr) 
      ELSE IF( NAVstr%instr .EQ. 2 ) THEN
!	Recompute elevation and scan biases based on user inputs of
!	cycles and increments obtained from GVAR.
       elvmax(NAVstr%instr)  = ( (9 - NAVstr%ns_cyl) *    &
                                incmax(NAVstr%instr) -     &
                                NAVstr%ns_inc ) *     &
                                elvinc(NAVstr%instr)
       scnmax(NAVstr%instr) = ( NAVstr%ew_cyl *     &
                                incmax(NAVstr%instr) + & 
                                NAVstr%ew_inc ) *  &
                                scninc(NAVstr%instr)
!      Compute elevation angle (rad)
       elev = elvmax(NAVstr%instr) + (2.50 - line) * elvln(NAVstr%instr) 
!      Compute scan angle (rad)
       scan = (pixel - 1.0)*scnpx(NAVstr%instr)-scnmax(NAVstr%instr)
      ELSE
!      Unknown instrument.....
       elev = 0.0D0
       scan = 0.0D0
      ENDIF
!
!
      RETURN
      END SUBROUTINE COMP_ES
!======================================================================
!     C O M P _ L P
!======================================================================
      SUBROUTINE COMP_LP( NAVstr, ELEV, SCAN, RL, RP )
!
!     AUTHOR:         KELLY DEAN
!
!     CREATED:        January 1995
!   
!     DEVELOPED FOR:  CIRA/COLORADO STATE UNIVERSITY
!
!     PURPOSE:        
!	  Subroutine COMP_LP converts elevation and scan angles to the
!       fractional line and pixel numbers.
!
!     REVISION:       0.0
!
!     ARGUMENTS:
!       NAME:   TYPE:       PURPOSE:                        IN/OUT:
!	 NAVstr   Structure   Navigation information         IN
!	 ELEV	  REAL*8     Elevation angle (rad)           IN
!	 SCAN	  REAL*8:    Scan angle (rad)                IN
!	 RL	  REAL*8     Line Number                     OUT
!	 RP	  REAL*8     Pixel Number                    OUT
!
!     CONSTANTS:
!       NAME:    PURPOSE:
!       ****************  REAL*8       *****************
!  	elvln    Elevation angle per detector line (rad)
! 	elvmax   Bounds in elevation
!	scnmax   Bounds in scan angle
!	scnpx    Scan angle per pixel (rad)
!
!     -------------------------------------------------------
!     --------------  ACTUAL CODE STARTS HERE  --------------
!     -------------------------------------------------------
!
!     INCLUDE DECLARATION SECTION:
!
!     include 'kdf.inc'
!
!     CONSTANT DECLARATION SECTION:
!
      INTEGER incmax(2) /6136,2805/
      REAL*8  elvmax(2) / 0.2208960D0, 0.22089375D0/
      REAL*8  elvln(2)  /28.0D-6,    280.0D-6/
      REAL*8  elvinc(2) / 8.0D-6,     17.5D-6/
      REAL*8  scnmax(2) / 0.245440D0,  0.2454375D0 /
      REAL*8  scnpx(2)  /16.0D-6,    280.0D-6/
      REAL*8  scninc(2) /16.0D-6,     35.0D-6/
!
!     VARIABLE DECLARATION SECTION:
!
      type(GVAR_NAV):: NAVstr
      REAL*8 ELEV, SCAN, RL, RP
!
!     ******--------------------------------------------******
!     ******--------  MAIN BODY STARTS HERE  -----------******
!     ******--------------------------------------------******
!
      IF ( NAVstr%instr .EQ. 1 ) THEN
!	Recompute elevation and scan biases based on user inputs of
!	cycles and increments obtained from GVAR.
       elvmax(NAVstr%instr) = ( NAVstr%ns_cyl *      &
                               incmax(NAVstr%instr) +      &
                               NAVstr%ns_inc ) *      &
                               elvinc(NAVSTR%instr)
       scnmax(NAVstr%instr) = ( NAVstr%ew_cyl *     &
                                incmax(NAVstr%instr) +      &
                                NAVstr%ew_inc ) *           &
                                scninc(NAVstr%instr)
!       Compute fractional line number.
      RL = ( ELVMAX(NAVstr%instr) - ELEV ) / ELVLN(NAVstr%instr) 
      RL = RL + 4.5D0
!       Compute fractional pixel number.
       RP = ( SCNMAX(NAVstr%instr) + SCAN ) / SCNPX(NAVstr%instr)+1.0D0 
      ELSE IF ( NAVstr%instr .EQ. 2 ) THEN
!       Recompute elevation and scan biases based on user inputs of
!       cycles and increments obtained from GVAR.
        elvmax(NAVstr%instr)  = ( (9 - NAVstr%ns_cyl) *    &
                                incmax(NAVstr%instr) -     &
                                NAVstr%ns_inc ) *     &
                                elvinc(NAVstr%instr)   
        scnmax(NAVstr%instr) = ( NAVstr%ew_cyl *     &
                                incmax(NAVstr%instr) +    & 
                                NAVstr%ew_inc ) *     &
                                scninc(NAVstr%instr)
!       Compute fractional line number.
        RL = ( ELVMAX(NAVstr%instr) - ELEV ) / ELVLN(NAVstr%instr) 
        RL = RL + 2.5D0
!       Compute fractional pixel number.
        RP = ( SCNMAX(NAVstr%instr)+SCAN ) / SCNPX(NAVstr%instr)+1.0D0 
      ELSE
!      Unknown instrument.....
       RL = 0.0D0
       RP = 0.0D0
      ENDIF
!
      RETURN
      END SUBROUTINE COMP_LP
end module GOES_IMAGER_ROUTINES
