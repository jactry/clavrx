!$Id: conver.f90 9 2014-01-31 08:19:35Z awalther $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: conver.f90 (src)
!       Conver (program)
!
! PURPOSE:
!     Converts narrowband reflectance to broadband albedo using spectral
!     conversions derived from simulated measurements of satellites
!     and anisotropic correction factors obtained from NIMBUS-7 ERB
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
! REVISION HISTORY:
! CVS info:
!
!  $Id: conver.f90 9 2014-01-31 08:19:35Z awalther $
!
!  $Log: conver.f90,v $
!  Revision 1.4.2.2  2014/01/26 04:48:32  heidinger
!  updated
!
!  Revision 1.2.2.2  2014/01/24 21:47:59  mhiley
!  Format headers according to CDR General Software Coding Standards
!
!  Revision 1.2.2.1  2013/04/18 16:17:41  dbotambekov
!  change eq to eqv for gfort
!
!  Revision 1.2  2011/08/16 18:45:27  heidinger
!  resolved conflicts
!
!  Revision 1.1.2.2  2011/02/22 23:20:43  heidinger
!  updated
!
!  Revision 1.1.2.1  2011/02/22 23:10:04  heidinger
!  added
!
!  Revision 1.2  2004/12/20 18:17:45  mikew
!  Added CVS tags.
!
! REFERENCES:  Laszlo, I., H. Jacobowitz and A. Gruber, 1988:
!              The relative merits of narrowband channels for
!              estimating broadband albedos,
!              J. Atmos. Oceanic Tech. 5, pp.757-773.   (LJG)
!
!              Taylor, V.R, and L.L. Stowe, 1986: Revised
!              reflectance and emission models from NIMBUS-7 ERB
!              data, Proceedings of the Sixth Conference on
!              Atmospheric Radiation, AMS, Williamsburg, Va.,
!              May 13-16.                                (TS)
!
!              Wydick, J.E, P.A. Davis and A. Gruber, 1987:
!              Estimation of Broadband Planetary Albedo from
!              Operational Narrowband Satellite Measurements,
!              NOAA Tech. Rep. NESDIS 27.               (WDG)
!
!     NOTE1: Since it is assumed that the narrowband reflectance
!            -Reflnb- corresponds to either a completely clear or
!            a completely cloudy pixel, only the clear and the
!            overcast subgroups of the NIMBUS-7 ERB reflectance
!            models are used.
!
!     NOTE2: Narrow- to broadband albedo conversion factors for
!            the NOAA-7/CH1 over snow is from WDG.
!
!     ROUTINE CALLED:  LOCATE
!
!
!   I N P U T :
!                  Clear  : true, if clear case, false if cloudy
!                  Flname : path and name containing bdrf data
!                  Isatid : satellite id.
!                            1 - noaa-7/ch1
!                            2 - goes-6/ch1
!                            3 - goes-5/ch1
!                            4 - meteosat-2/ch1
!                            5 - gms-2/ch1
!                            6 - goes-8/ch1
!                  Isrtyp : surface type index (1-water, 2-vegetation,
!                            3-desert, 4-snow/ice)
!                  Nobdrc : if true, no bidirectional correction is done
!                  Nontob : if true, assumes that broadband albedo is
!                              the same as the narrowband reflectance
!                  Reflnb : narrowband reflectance ( 0.0 - 1.0 )
!                  Relaz  : relative azimuth angle (degrees)
!                              measured from the half-plane containing
!                              the sun to the half-plane containing the
!                              satellite (Range: 0-180)
!                  Satmu  : cosine of satellite zenith angle
!                  Sunmu  : cosine of solar zenith angle
!
!   O U T P U T :  Bbalb  : broadband albedo ( 0.0 - 1.0 )
!                  Noretr : true, if no narrow- to broadband conversion
!                           exists for satellite -isatid-
!
!   I N T E R N A L   V A R I A B L E S:
!
!   Aniscr(Iaz,Ist,Isn,Ity) : Iaz=1 to Nazima; Ist=1 to Nsatan; Isn=1 to
!                   Nsunan; Ity=1 to Ntypes; bidirectional
!                   correction factors
!   Azimut(Iaz)       : Iaz=1 to Nazima; azimuth angles in bdr
!                   tables
!   Dimerr            : true if dimension error is detected
!   First             : true on first entry, false otherwise
!   Head(Ity)         : Ity=1 to Ntypes; header for describing
!                   table
!   Iscene            : scene id (1-clear ocean; 2-clear veget;
!                   3-clear desert; 4-clear snow, 5-overcast)
!   Nazima            : no. of relative azimuth angles in tables
!   Nsatan            : no. of satellite zenith angles in tables
!   Nsatel            : no. of satellites for which narrow- to
!                   broadband conversion is implemented
!   Nscene            : no. of scene types for which narrow- to
!                   broadband conversion is implemented
!   Nsunan            : no. of solar zenith angles in tables
!   Ntypes            : no. of scene types in tables
!   Offset(Iscene,Isatid)   : Iscene=1 to Nscene, Isatid=1 to Nsatel;
!                             offset in the narrow- to broadband
!                             conversion equation
!   Slope(Iscene,Isatid)    : as above, but slope
!
!   LOCAL SYMBOLIC DIMENSIONS :
!
!   Maxazm : max no. of azimuth angles in bdr tables
!   Maxsat : max no. of satellite zenith angles in bdr tables
!   Maxsun : max no. of solar zenith angle values in bdr tables
!   Maxtyp : max no. of scene types in bdr tables
!--------------------------------------------------------------------------------------

    SUBROUTINE Conver( Unit16,    &
                       Nontob, Nobdrc, Reflnb, Isrtyp, Clear, Isatid, &
                       Relaz, Satmu, Sunmu, Bbalb, Noretr, Flname, &
                       Glint )

! ..
! .. Parameters ..
      INTEGER, PARAMETER :: Maxazm=8, Maxsat=7, Maxsun=10, &
        Maxtyp=12, Nsatel=6, Nscene=5, Maxdim=3
! ..
! .. Scalar Arguments ..
      INTEGER :: Unit16
      REAL :: Bbalb, Glint, Reflnb, Relaz, Satmu, Sunmu
      INTEGER :: Isatid, Isrtyp
      LOGICAL :: Clear, Nobdrc, Nontob, Noretr
      CHARACTER (*) :: Flname
! ..
! .. Local Scalars ..
      INTEGER :: I, Iazm(2), Iscene, Isnm(2), Istm(2), Ity, Ntypes
      INTEGER, SAVE :: Nazima, Nsatan, Nsunan
      LOGICAL :: Dimerr, Adhoc
      LOGICAL, SAVE :: First
      LOGICAL:: Wrtdim
      INTEGER K, J, Mm, Ndim, Iaz, Ist, Isn
      REAL Fraz,Frst,Frsn,Anifac
! ..
! .. Local Arrays ..
      REAL, SAVE :: Aniscr(Maxazm,Maxsat,Maxsun,Maxtyp), &
        Azimut(Maxazm+1), Satzna(Maxsat+1), Sunzna(Maxsun+1), &
        Azim(Maxazm), Satz(Maxsat), Sunz(Maxsun)
      REAL :: Offset(Nscene,Nsatel), Slope(Nscene,Nsatel)
      CHARACTER (80) :: Head(Maxtyp)
      REAL :: Array(2**Maxdim), Fr(Maxdim)
! ..
! .. Data Statements ..
      DATA First/ .TRUE./
      DATA (Slope(I,1),I=1,Nscene)/0.902, 0.779, 0.804, 0.760, 0.780/
      DATA (Slope(I,2),I=1,Nscene)/0.883, 0.816, 0.800, -0.155, 0.786/
      DATA (Slope(I,3),I=1,Nscene)/0.883, 0.816, 0.795, -0.156, 0.786/
      DATA (Slope(I,4),I=1,Nscene)/0.945, 0.945, 0.909, 0.803, 0.787/
      DATA (Slope(I,5),I=1,Nscene)/0.861, 0.826, 0.776, -0.151, 0.783/
      DATA (Slope(I,6),I=1,Nscene)/0.841, 0.794, 0.838, 0.715, 0.930/
      DATA (Offset(I,1),I=1,Nscene)/0.01426, 0.06831, 0.02819, 0.0083, &
        0.05004/
      DATA (Offset(I,2),I=1,Nscene)/0.01100, 0.05192, 0.02611, 0.7856, &
        0.04646/
      DATA (Offset(I,3),I=1,Nscene)/0.01051, 0.05170, 0.02708, 0.7857, &
        0.04583/
      DATA (Offset(I,4),I=1,Nscene)/0.01591, 0.00524, 0.06137, &
        -0.0098, 0.04421/
      DATA (Offset(I,5),I=1,Nscene)/0.00941, 0.04584, 0.02983, 0.7831, &
        0.04397/
      DATA (Offset(I,6),I=1,Nscene)/0.01422, 0.01659, 0.01715, &
        0.13650, -0.05430/
! ..

      Noretr = .FALSE.
      Bbalb = Reflnb
!                       ** ASSUME BROADBAND ALBEDO IS THE SAME AS
!                       ** THE NARROWBAND REFLECTANCE

      IF (Nontob .AND. Nobdrc) RETURN

!                       ** DO NARROW- TO BROADBAND CONVERSION

      IF ( .NOT. Nontob .AND. (Isatid<1 .OR. Isatid>6)) THEN

!                       ** DO NOT RETRIEVE IF UNKNOWN SATELLITE

        Noretr = .TRUE.
        CALL Errmsg(Unit16,'CONVER--UNKNOWN SATELLITE',.FALSE.)
        RETURN

      END IF

      IF (Reflnb<1E-6) RETURN

      IF (Clear) THEN

        Iscene = Isrtyp

      ELSE

        Iscene = 5

      END IF

      IF ( .NOT. Nontob) Bbalb = Slope(Iscene,Isatid)*Reflnb + &
        Offset(Iscene,Isatid)

      IF ( Bbalb.LE.0.0 ) Bbalb = Reflnb

      IF (Nobdrc) RETURN

!                       ** APPLY BIDIRECTIONAL CORRECTIONS

      IF (First) THEN

        First = .FALSE.

!                       **  READ ANISOTROPIC CORR FACTORS FROM TABLE

        OPEN (Unit=9,File=Flname,Status='OLD')

        READ (9,Fmt=90000,End=10,Err=20) Ntypes, Nazima, Nsatan, &
          Nsunan

        Dimerr = .FALSE.

        IF (Nazima>Maxazm) Dimerr = Wrtdim('MAXAZM',Nazima)
        IF (Nsatan>Maxsat) Dimerr = Wrtdim('MAXSAT',Nsatan)
        IF (Nsunan>Maxsun) Dimerr = Wrtdim('MAXSUN',Nsunan)
        IF (Ntypes>Maxtyp) Dimerr = Wrtdim('MAXTYP',Ntypes)
        IF (Dimerr) CALL Errmsg(Unit16,'BDRFAC--DIMENSION ERROR(S)',.TRUE.)

        READ (9,Fmt=90010,End=10,Err=20) (Azimut(Iaz),Iaz=1,Nazima+1)
        READ (9,Fmt=90010,End=10,Err=20) (Satzna(Ist),Ist=1,Nsatan+1)
        READ (9,Fmt=90010,End=10,Err=20) (Sunzna(Isn),Isn=1,Nsunan+1)

        DO Ity = 1, Ntypes

          READ (9,Fmt='(A)',End=10,Err=20) Head(Ity)

          DO Isn = 1, Nsunan

            DO Ist = 1, Nsatan

              READ (9,Fmt=90020,End=10,Err=20) &
                (Aniscr(Iaz,Ist,Isn,Ity),Iaz=1,Nazima)

            END DO

          END DO

        END DO

        CLOSE (Unit=9)

      END IF
      
      DO Iaz = 1, Nazima
        Azim( Iaz ) = ( Azimut(Iaz) + Azimut(Iaz+1) ) / 2.
      END DO
      DO Ist = 1, Nsatan
        Satz( Ist ) = ( Satzna(Ist) + Satzna(Ist+1) ) / 2
      END DO
      DO Isn = 1, Nsunan
        Sunz( Isn ) = ( Sunzna(Isn) + Sunzna(Isn+1) ) / 2
      END DO

!                       ** Assume symmetry with respect to the plane
!                       ** containing the sun and the zenith

      IF (Relaz>180.0) Relaz = 360.0 - Relaz

!                       ** LOCATE ANGLES IN THEIR TABLE VALUES

      CALL Locate( Relaz, Azim, Nazima, Iazm(1), Iazm(2), Fraz )
      CALL Locate( Satmu, Satz, Nsatan, Istm(1), Istm(2), Frst )
      CALL Locate( Sunmu, Sunz, Nsunan, Isnm(1), Isnm(2), Frsn )

!                       ** ASSOCIATE SCENE TYPES

      IF (Isrtyp==1) THEN

        Ity = 1

      ELSE IF (Isrtyp==2) THEN

        Ity = 2

      ELSE IF (Isrtyp==3) THEN

        Ity = 4

      ELSE IF (Isrtyp==4) THEN

        Ity = 3

      END IF

      IF ( .NOT. Clear) Ity = 12

      Ndim  = 3
      
      Fr( 1 ) = Fraz
      Fr( 2 ) = Frst
      Fr( 3 ) = Frsn
      
      Mm = 0
      
      DO K = 1, 2
         DO J = 1, 2
            DO I = 1, 2
              
               Mm = Mm + 1
               Array( Mm ) = Aniscr( Iazm(I), Istm(J), Isnm(K), Ity )
               
            END DO
         END DO
      END DO

      
      Anifac = Rintnd( Fr, Ndim, Array )

      Bbalb = Bbalb / anifac

!             ** Do additional angular corrections; this is only an
!             ** ad hoc empirical correction based on adjusting the
!             ** constants until the glint more or less disappears 
!     
      Adhoc = .false.
      IF ( Adhoc ) THEN
        
         IF ( Clear .AND. Iscene == 1 .AND. & 
              Glint.LE.30.0 .AND. Glint.GT.20.0 ) Bbalb = Bbalb / 1.2
         IF ( Clear .AND. Iscene == 1 .AND. & 
            Glint.LE.20.0 .AND. Glint.GT.10.0 ) Bbalb = Bbalb / 1.35
         IF ( Clear .AND. Iscene == 1 .AND. &
            Glint.LE.10.0 ) Bbalb = Bbalb / 1.45
         IF ( Clear .AND. Iscene == 1 .AND. &
            Glint.LE.20.0 .AND. Satmu.GE.0.95 .AND. Relaz.GE.150.0 ) &
            Bbalb = Bbalb / 1.5

         IF ( Clear .AND. Iscene == 1 .AND. &
            Glint.LE.35.0 .AND. Glint.GT.20.0 .AND.  &
            Relaz.LE.160.0 .AND. Relaz.GT.130.0 ) &
            Bbalb = Bbalb / 1.35
      
         IF ( Bbalb.LE.0.0 .OR. Bbalb.GE.1.0 ) then
           Noretr = .TRUE.
          !print *, Iscene, Anifac, Reflnb, Bbalb
         end if

      END IF

      RETURN

10    WRITE (6,Fmt='(//,1X,A)') 'CONVER--PREMAT END OF BDR FILE'
20    WRITE (6,Fmt='(//,1X,A)') 'CONVER--ERROR READING BDR TABLE'

      STOP

90000 FORMAT (4I3)
90010 FORMAT (11F6.3)
90020 FORMAT (8F10.5)

      CONTAINS

         FUNCTION Rintnd( Fr, Ndim, Array )
!
!        Linearly interpolates in multi-dimensions
!
!   INPUT VARIABLES
!     Array(I) : I=1 to 2**Ndim, a one-dimensional array containing
!                  the data to be interpolated on
!     Fr(I)    : I=1 to Ndim, fractional weight in dimension -I-
!     Ndim     : dimension of interpolation
!
!   OUTPUT VARIABLE
!     Rintnd   : interpolated value

! ..
! .. Scalar Arguments ..
      INTEGER :: Ndim
! ..
! .. Array Arguments ..
      REAL :: Array(2**Ndim), Fr(Ndim)
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: Maxdim = 32
! ..
! .. Local Scalars ..
      INTEGER :: I, J, K
      LOGICAL :: Dimerr
      LOGICAL, SAVE :: First
! ..
! .. Local Arrays ..
      REAL :: Ytmp(Maxdim)
! ..
! .. Data Statements ..
      DATA First/ .TRUE./
! ..
! .. Function Return Value ..
      REAL :: Rintnd
! ..
      IF (First) THEN
!                   ** Check local dimension

        Dimerr = .FALSE.
        First = .FALSE.

        IF (Maxdim<2**Ndim) Dimerr = Wrtdim('Maxdim',2**Ndim)

        IF (Dimerr) CALL Errmsg(Unit16,'Rintnd--DIMENSION ERROR',.TRUE.)

      END IF

      DO I = 1, 2**Ndim

        Ytmp(I) = Array(I)

      END DO

      DO I = 1, Ndim

        K = 1

        DO J = 1, 2**(Ndim-I)

          Ytmp(J) = Fr(I)*Ytmp(K) + (1.-Fr(I))*Ytmp(K+1)
          K = K + 2

        END DO

      END DO

      Rintnd = Ytmp(1)

      RETURN

    END FUNCTION Rintnd

    END SUBROUTINE Conver

