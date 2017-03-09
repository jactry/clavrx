C$Id: sasrab.f 1433 2015-10-30 16:43:33Z cmolling $
c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: sasrab.f (src)
c
c PURPOSE: Satellite Algorithm for Shortwave RAdiation Budget
c
c DESCRIPTION: 
c
c AUTHORS:
c  Andrew Heidinger, Andrew.Heidinger@noaa.gov
c  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
c  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
c  William Straka, CIMSS, wstraka@ssec.wisc.edu
c
c COPYRIGHT
c THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
c DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
c FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
c INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
c EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
c DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
c THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
c SUPPORT TO USERS.
c
c--------------------------------------------------------------------------

C**********************************************************************C
      SUBROUTINE  Sasrab( Auxpath,
     &                    Glat, Glon, Year, Month, Jday, Gmtime,
     &                    Snowfr, Cdfrac, Cldref, Clrref, Cmpref, Ozone,
     &                    Pwater, Solmu, Satmu, Relaz, Glint, Isatid,
     &                    Misval, Mxfpr, Mxint, Flxall, Flxclr, Aertau,
     &                    Cldtau, Nfpar, Nintv )

C**********************************************************************C
C
C            Satellite Algorithm for Shortwave RAdiation Budget
C
C     VERSION 1.0
C     Last revision: Oct 05 2002
C
C======================================================================C
C
C INPUT
C
C   Cdfrac   : cloud cover fraction (0-1)
C   Cldref   : mean narrowband reflectance from cloudy pixels (0-1)
C   Clrref   : mean narrowband reflectance from clear pixels (0-1)
C   Cmpref   : narrowband reflectance from clear(dark) sky composite (0-1)
C   Glat     : latitude of center of grid (degrees)
C                 range: -90.0 (south) - + 90.0 (north)
C   Glint    : glint angle (degrees)
C   Glon     : longitude of center of grid (degrees)
C                 range: -180.0 (west) - + 180.0 (east)
C   Gmtime   : Greenwich Mean Time (hours)
C   Isatid   : satellite ID; see subroutine Conver for valid values
C   Jday     : day of year (1-365/366)
C   Misval   : value representing missing data
C   Month    : month (1-12)
C   Mxfpr    : max number of flux outputs; first dimension of Flxall and
C                  Flxclr
C   Mxint    : max number of output spectral intervals; second dimension
C                  of Flxall and Flxclr
C   Ozone    : column amount of ozone in atm-cm
C   Pwater   : column amount of precipitable water (cm)
C   Relaz    : relative azimuth angle (degrees)
C   Snowfr   : snow cover fraction (0-1)
C   Satmu    : cosine of satellite zenith angle
C   Solmu    : cosine of the solar zenith angle
C   Year     : year (e.g., 2001)
C
C OUTPUT
C
C   Aertau   : aerosol optical depth at 0.55 microns
C   Cldtau   : cloud optical depth at 0.55 microns
C   Flxall(Ip,Is) : Ip=1 to Npfar, Is=1 to Nintv; all-sky flux (W/m^2)
C   Flxclr(Ip,Is) : Ip=1 to Npfar, Is=1 to Nintv; clear-sky flux (W/m^2)
C   Nfpar    : number of flux-type parameters
C   Nintv    : number of spectral intervals
C
C+---------------------------------------------------------------------+
C   ROUTINES CALLED (IN ORDER): Chkset, Prtinp, Reftra, Sundat, Getmod,
C                               Auxdat, Insest
C+---------------------------------------------------------------------+
C
C+---------------------------------------------------------------------+
C     L O C A L   V A R I A B L E S:
C
C   Iaer          : aerosol model id
C   Ilat          : latitude index of gridpoint in Matthes' maps
C   Ilon          : longitude index of gridpoint in Matthes' maps
C   Isrmod        : surface albedo model id
C   Isrtyp        : surface type id (1-water, 2-vegetation, 3-desert,
C                                    4-Snow/ice)
C   Itb           : pointer indicating the refl/transt able to be used
C                       for the current grid
C   Ntable        : no. of refl/trans (RT) tables
C   Ntaua(Ib)     : Ib=1 to Ntable; no. of aerosol optical depths in
C                        refl/trans (RT) tables
C   Ntauc(Ib)     : Ib=1 to Ntable; no. of cloud optical depths in
C                        refl/trans (RT) tables
C   Nwlint(Ib)    : Ib=1 to Ntable; no. of spectral interval in RT table
C   Aertau        : aerosol optical depth at 0.55 microns
C   Cldtau        : cloud optical depth at 0.55 microns
C   Rsuntb(Ib)    : Ib=1 to Ntable; sun-earth distance factor in
C                        refl/trans (RT) tables
C   Sclaer        : scale factor for aerosol optical depth; when writing
C                       out in table, the optical depth is multiplied by
C                       -Sclaer- (set in subroutine *Chkset*).
C   Sclcld        : as -sclaer-, but for cloud optical depth
C   Solctb(Ib)    : ib=1 to ntable; solar constant (W/m sq) in RT table
C   Solirr(Il,Ib) : il=1 to nwlint, ib=1 to ntable; spectral solar
C                        Irradiance (W/m sq) at toa in refl/trans table
C   Srfddw        : diffuse surface down-flux (W/m sq)
C   Srfdnw        : surface down-flux (W/m sq) (direct+diffuse)
C   Srfupw        : surface up-flux (W/m sq)
C   Sundis        : sun-earth distance factor for current day
C   Aertab(It,Ib) : It=1 to Ntaua(Ib), Ib=1 to Ntable; aerosol optical
C                      depth values at o.55 microns in refl/trans table
C   Cldtab(It,Ib) : It=1 to Ntauc(Ib), Ib=1 to Ntable; cloud optical
C                      depth values at 0.55 microns in refl/trans table
C   Tauclm        : climatological value of aerosol optical depth at
C                      0.55 microns for current grid
C   Toadnw        : top of the atmosphere down-flx (W/m sq)
C   Toaupw        : top of the atmosphere up-flx (W/m sq)
C   Utau          : if true, read in user supplied monthly-mean visible
C                     aerosol optical depth data
C   Waveln(Il,Ib) : Il=1 to Nwlint(Ib)+1, Ib=1 to Ntable; boundaries of
C                   spectral intervals (microns)
C+---------------------------------------------------------------------+
C
C     LOCAL SYMBOLIC DIMENSIONS 1:
C
C     Maxtab  = max no. of refl/trans tables
C     Maxtau  = max no. of optical depth values in refl/trans tables
C     Maxwvl  = max no. of spectral intervals in refl/trans tables
C+---------------------------------------------------------------------+
C

C     .. Parameters ..
      INTEGER  Maxtab, Maxtau, Maxwvl, Maxflx, Maxotr, Maxint, Maxfpr
      PARAMETER  ( Maxtab=2, Maxtau=15, Maxwvl=5, Maxflx=27, Maxotr=5,
     &           Maxint=3, Maxfpr=5 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER * 120 Auxpath
      INTEGER Get_Lun
      INTEGER Unit16
      INTEGER Unit1
      REAL  Aertau, Cdfrac, Cldbba, Cldref, Cldtau, Clrbba, Clrref,
     &      Cmpbba, Cmpref, Glat, Glint, Glon, Gmtime, Misval, Ozone,
     &      Pwater, Relaz, Satmu, Snowfr, Solmu
      INTEGER  Isatid, Jday, Month, Mxfpr, Mxint, Nfpar, Nintv,
     &         Year
C     ..
C     .. Array Arguments ..
      REAL  Flxall( Mxfpr, Mxint ), Flxclr( Mxfpr, Mxint )
C     ..
C     .. Local Scalars ..
      REAL  Declin, Smumax, Smumin, Solcon, Sundis, Tauclm, Vernum,
     &      Vertab
      INTEGER   
     &         Ilat, Ilon, Isrmod, Isrtyp, Itb, Mxnh2o, Mxno3, Mxnsun,
     &         Mxntab, Mxntau, Mxnwvl, Nflux, Nother, Ntable, Prvday
      LOGICAL  Climat, Extrap, First, Matthe, Newday, Nobdrc, Nontob,
     &         Prnt, Satalb, Utau
      CHARACTER  * 8 Cdate
      CHARACTER  * 10 Ctime
      CHARACTER  * 80 Header
      CHARACTER  * 256 Fnaer, Fnbdr, Fnrt, Fnvt
C     ..
C     .. Local Arrays ..
      REAL  Aertab( Maxtau, Maxtab ), Cldtab( Maxtau, Maxtab ),
     &      Rsuntb( Maxtab ), Solctb( Maxtab ),
     &      Solirr( Maxwvl, Maxtab ), Waveln( Maxwvl + 1, Maxtab )
      INTEGER  Iaer( Maxtab ), Iw1( Maxint ), Iw2( Maxint ),
     &         Ntaua( Maxtab ), Ntauc( Maxtab ), Nwlint( Maxtab )
C     ..
C     .. External Subroutines ..
      EXTERNAL  Auxdat, Chkset, Getmod, Insest, Prtinf, Reftra,
     &          Sundat
C     ..
C     .. Save statement ..
      SAVE  Itb, Isrtyp, Tauclm, Smumax, Smumin, Extrap, Solcon, Vernum,
     &      Nobdrc, Nontob, Climat, Satalb, Matthe, Prnt, Utau, Nflux,
     &      Nother, Fnrt, Fnvt, Fnbdr, Fnaer, First, Prvday, Ntable,
     &      Unit16, Unit1, Iaer
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC  Aint, Int, Mod
C     ..
C     .. Intrinsic Subroutines ..
      INTRINSIC  Date_and_time
C     ..
C     .. Data statements ..
      DATA  Prvday / - 1 / , First / .true. /
C     ..


      IF ( First ) THEN

         Vernum = 1.0

         Solcon = 1372.6
         Smumin = 0.2
         Smumax = 1.0

         Extrap = .true.
         Nobdrc = .false.
         Nontob = .false.
         Climat = .false.
         Satalb = .true.
         Matthe = .false.
         Prnt   = .false.
         Utau   = .false.

         Header = 'SASRAB adopted for new GSIP'

C                       ** Set up path and names of auxiliary input data
C
C***** AKH modified PATH was Auxpath
         Fnrt   = trim(Auxpath)//'reftra.dat'
         Fnvt   = trim(Auxpath)//'land_cover_umd_8km.dat'
         Fnbdr  = trim(Auxpath)//'bdrcor.dat'
         Fnaer  = trim(Auxpath)//'aertab.dat'

C                       ** Set up path and names of output files
C
C----- heidinger moved this to the main code - before call in Insol
C        Fnlst  = 'Results/yyyymmdddhhmm.lst'
C        WRITE ( Fnlst(9:21), FMT='(I4.4,I2.2,I3.3,2I2.2)' ) Year,
C    &     Month, Jday, Int( Gmtime ), Int( Mod(Gmtime,Aint(Gmtime)) *
C    &     60 )
C        OPEN ( Unit=16, File=Fnlst, Form='FORMATTED', Status='UNKNOWN' )
C---------------------------------------------------------------------------
       CALL Date_and_time( Date=Cdate, Time=Ctime )
       Unit16 = Get_Lun()
       OPEN(Unit=Unit16,file = 'SASRAB_LogFile', status='REPLACE')
       WRITE (Unit16, FMT= * ) 'Starting date and time: ', Cdate( 1:4 ),
     &     '/', Cdate( 5:6 ), '/', Cdate( 7:8 ), ' ', Ctime( 1:2 ), ':',
     &     Ctime( 3:4 ), ':', Ctime( 5:10 )

C
C                       ** Check dimensions and variables, and
C                       ** do various setup tasks
C
         CALL Chkset( Unit16,
     &                Smumin, Smumax, Solcon, Matthe, Satalb, Maxflx,
     &                Maxotr, Maxfpr, Maxint, Mxfpr, Mxint, Nflux,
     &                Nother, Nfpar, Nintv, Iw1, Iw2 )
C
C                       ** Print information on options and parameters
C

         CALL Prtinf( Unit16, 
     &                Smumin, Smumax, Solcon, Climat, Nontob, Nobdrc,
     &                Satalb, Matthe, Extrap, Prnt, Header, Vernum,
     &                Misval, Utau )
C
C                       ** Read refl-trans table(s) and print info
C                       ** about them
C
         CALL Reftra( Unit16, 
     &                Maxtab, Maxwvl, Maxtau, Iaer, Ntable, Ntaua,
     &                Ntauc, Nwlint, Rsuntb, Solirr, Solctb, Aertab,
     &                Cldtab, Waveln, Mxntau, Mxnsun, Mxnwvl, Mxnh2o,
     &                Mxno3, Mxntab, Solcon, Vertab, Fnrt )

         First  = .false.

      END IF

      Newday = Jday .NE. Prvday
      Prvday = Jday

C                       ** Calculate sun-earth distance factor and
C                       ** declination for current day
C
      IF ( Newday ) CALL Sundat( Jday, Sundis, Declin )
C
C
C                       ** Select aerosol and surface albedo model
C                       ** representative for current grid, and get
C                       ** climatological aerosol optical depth
C

!--------------------------------------------------------------------
! Heidinger Mod
!--------------------------------------------------------------------
!     Isrtyp = cell_sfc_type

      Unit1 = Get_Lun()

      CALL Getmod( Unit16,Unit1,  
     &             Fnvt, Fnaer, Utau, Glat, Glon, Month, Ntable, Iaer,
     &             Itb, Isrmod, Isrtyp, Tauclm, Ilat, Ilon )
C
C                       ** Obtain auxiliary atmospheric data
C
      CALL Auxdat( Climat, Jday, Glat, Pwater, Ozone )
C
C                       ** Estimate instantaneous values of toa/surface
C                       ** fluxes and aerosol/cloud optical depths
C
      CALL Insest( Unit16, 
     &             Cdfrac, Clrref, Cldref, Cmpref, Pwater, Ozone,
     &             Isatid, Snowfr, Solmu, Satmu, Relaz, Glint,
     &             Solirr(1,Itb), Aertab(1,Itb), Cldtab(1,Itb),
     &             Waveln(1,Itb), Jday, Glat, Glon, Gmtime, 
     &             Ntaua(Itb), Ntauc(Itb), Nwlint(Itb), Itb, Isrmod,
     &             Isrtyp, Tauclm, Sundis, Declin, Nontob, Matthe,
     &             Satalb, Nobdrc, Prnt, Extrap, Smumin, Smumax, Ilat,
     &             Ilon, Mxntau, Mxnsun, Mxnwvl, Mxnh2o, Mxno3, Mxntab,
     &             Maxfpr, Maxint, Misval, Nfpar, Nintv, Iw1, Iw2,
     &             Fnbdr, Flxall, Flxclr, Aertau, Cldtau, Cmpbba,
     &             Clrbba, Cldbba )


      RETURN

 9000 FORMAT ( 1X, '>>> PROCESSING DATA FOR DAY:', 5X, I3 )

      END

C**********************************************************************C
      SUBROUTINE  Aertab( Unit16, 
     $                    Flname, Glat, Glon, Month, 
     $                    Mistau, Aertau )
C**********************************************************************C
C
C        Reads equal-angle monthly mean visible aerosol optical depth
C        data, and finds the aerosol optical depth for the current
C        location. Returns missing value if aerosol optical depth for
C        current location is not available.
C
C  INPUT:
C
c    Unit16 : Unit number for diagnostic output
C    Flname : name of file containing aerosol optical depth data
C    Glat   : latitude (deg)
C    Glon   : longitude (deg)
C    Month  : month (1-12)
C
C  OUTPUT:
C
C    Aertau : monthly-mean aerosol optical depth
C    Mistau : value representing missing aerosol data in table of
C               monthly-mean aerosol optical thicknesses
C
C  LOCAL VARIABLES:
C
C    Clat1  : center latitude (deg) of first cell in aerosol table
C    Clatn  : center latitude (deg) of last cell in aerosol table
C    Clon1  : center longitude (deg) of first cell in aerosol table
C    Clonn  : center longitude (deg) of last cell in aerosol table
C    Dlat
C      Dlon : latitude and longitude resolution (deg)
C    Ilat
C      Ilon : latitude and longitude indices of current location in
C                aerosol table
C    Maxlat : max number of latitudes in aerosol table
C    Maxlon : max number of longitudes in aerosol table
C    Month1 : first month in table
C    Monthn : last month in table
C    Nlat   : number of latitudes in aerosol table
C    Nlon   : number of longitudes in aerosol table
C    Nmonth : number of months in aerosol table
C+---------------------------------------------------------------------+
C
C     .. Parameters ..
      INTEGER  Maxlat, Maxlon
      PARAMETER  ( Maxlat=51, Maxlon=111 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER  Flname * ( * )
      INTEGER  Month
      REAL  Aertau, Glat, Glon, Mistau
C     ..
C     .. Local Scalars ..
      CHARACTER  Hdline * 80
      LOGICAL  Newmon
      INTEGER  Ilat, Ilon, Irec, Lat, Lon, Month1, Monthn, Nlat, Nlon,
     &         Nmonth, Nskip, Pmonth
      REAL  Clat1, Clatn, Clon1, Clonn, Dlat, Dlon, Length
      INTEGER Get_Lun
      INTEGER Unit16, Unit2
C     ..
C     .. Local Arrays ..
      REAL  Tau( Maxlon, Maxlat )
C     ..
C     .. Save statement ..
      SAVE  Pmonth, Tau
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC  Int, Mod, Trim
C     ..
C     .. Data statements ..

      DATA  Pmonth / -1 /
C     ..
      Newmon = Month .NE. Pmonth

      IF ( Newmon ) THEN
C
C                   ** Read data for all cells for the current month
C
         Pmonth = Month

         Unit2 = Get_Lun()
         OPEN (Unit=Unit2,File=Flname, Form='UNFORMATTED', Status='OLD',
     &        Access='DIRECT', Recl=111*4, ERR=20 )

         READ ( Unit=Unit2, Rec=1 ) Hdline
         READ ( Unit=Unit2, Rec=2 ) Nlat, Nlon, Nmonth, 
     &     Clat1, Clatn, Clon1,
     &     Clonn, Month1, Monthn, Dlat, Dlon, Mistau

         Nskip  = 2 + ( Month-1 ) * Nlat
         Irec   = Nskip

         DO 10 Lat = 1, Nlat
            Irec   = Irec + 1

            READ (Unit=Unit2, Rec=Irec ) ( Tau(Lon,Lat), Lon=1, Nlon )
10       CONTINUE

         CLOSE ( Unit=Unit2 )

      END IF

      IF ( Glat.LT.Clat1-Dlat/2. .OR. Glat.GT.Clatn+Dlat/2. .OR.
     &     Glon.LT.Clon1-Dlon/2. .OR. Glon.GT.Clonn+Dlon/2. ) THEN
C
C                    ** No data for current location in aerosol table
C
         Aertau = Mistau
         RETURN

      END IF

C                    ** Find latitude and longitude indices
C
      Length = Glat - ( Clat1-Dlat/2. )
      Ilat   = Int( Length/Dlat )
      IF ( Mod(Length,Dlat).NE.0 ) Ilat   = Ilat + 1
      IF ( Ilat == 0 ) Ilat   = 1

      Length = Glon - ( Clon1-Dlon/2. )
      Ilon   = Int( Length/Dlon )
      IF ( Mod(Length,Dlon).NE.0 ) Ilon   = Ilon + 1
      IF ( Ilon == 0 ) Ilon   = 1

      Aertau = Tau( Ilon, Ilat )

      RETURN

20    WRITE (*, FMT='(/,1X, 2A)' ) 'Error opening file: ',
     &  Trim( Flname )
      STOP

      END

C**********************************************************************C
      SUBROUTINE  Auxdat( Climat, Jday, Glat, Pwater, Ozone )
C**********************************************************************C
C
C     Gets climatological values of precip water and/or ozone amounts if
C     variable -Climat- is 'true'. If -Climat- is 'false' and both ozone
C     and water vapor amounts are provided (Pwater>=0, Ozone>=0) simply
C     returns to calling program, otherwise gets missing value from
C     climatology.
C
C     REFERENCE:  KNEIZYS, ET. AL., ATMOSPHERIC TRANSMITTANCE/RADIANCE--
C                 COMPUTER CODE LOWTRAN5, AFGL-TR-80-0067 (FEB. 1980)
C
C
C     INPUT:  CLIMAT : IF TRUE, USE CLIMATOLOGICAL OZONE AND WATER VAPOR
C             GLAT   : LATITUDE IN DEGREES
C             JDAY   : DAY OF YEAR
C             PWATER : COLUMN AMOUNT OF WATER VAPOR IN G/SQ CM
C                      (IF CLIMAT=FALSE)
C             OZONE  : COLUMN AMOUNT OF OZONE IN ATM-CM
C                      (IF CLIMAT=FALSE)
C
C     OUTPUT: PWATER : COLUMN AMOUNT OF WATER VAPOR IN G/SQ CM
C                      (IF CLIMAT=TRUE, OR PWATER<0 ON INPUT)
C             OZONE  : COLUMN AMOUNT OF OZONE IN ATM-CM
C                      (IF CLIMAT=TRUE, OR OZONE<0 ON INPUT)
C
C     Called by- Sasrab
C     Calls- none
C+---------------------------------------------------------------------+
C
C     .. Scalar Arguments ..
      LOGICAL  Climat
      INTEGER  Jday
      REAL  Glat, Ozone, Pwater
C     ..
C     .. Local Scalars ..
      LOGICAL  Amjjas, Jfmond
      REAL  H2oamt, O3amt
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC  Abs
C     ..

C                       ** NO CLIMATOLOGY REQUIRED AND BOTH WATER VAPOR
C                       ** AND OZONE ARE PROVIDED -> RETURN
C
      IF ( .NOT.Climat .AND. Pwater.GE.0.0 .AND. Ozone.GE.0.0 ) RETURN
C
C                       ** TROPICAL MODEL
C
      IF ( Glat.GE.-25.0 .AND. Glat.LE.+25.0 ) THEN
         O3amt  = 0.246
         H2oamt = 4.12
         GO TO 10

      END IF
C
      Amjjas = Jday .GE. 91 .AND. Jday .LE. 273
      Jfmond = Jday .LT. 91 .OR. Jday .GT. 273
C
C                       ** MIDLATITUDE MODEL
C
      IF ( Abs(Glat).GT.25.0 .AND. Abs(Glat).LE.55.0 ) THEN
         IF ( Glat.GT.0.0 .AND. Amjjas .OR.
     &        Glat.LT.0.0 .AND. Jfmond ) THEN
C
C                                         ** SUMMER
            O3amt  = 0.318
            H2oamt = 2.92

         ELSE
C                                         ** WINTER
            O3amt  = 0.396
            H2oamt = 0.85

         END IF

         GO TO 10

      END IF
C
C                       ** SUBARCTIC MODEL
C
      IF ( Abs(Glat).GT.55.0 ) THEN
         IF ( Glat.GT.0.0 .AND. Amjjas .OR.
     &        Glat.LT.0.0 .AND. Jfmond ) THEN
C
C                                         ** SUMMER
            O3amt  = 0.344
            H2oamt = 2.09

         ELSE
C                                         ** WINTER
            O3amt  = 0.478
            H2oamt = 0.42

         END IF

      END IF
C
C                       ** USE CLIMATOLOGY
C
10    IF ( Climat ) THEN
         Pwater = H2oamt
         Ozone  = O3amt
         RETURN

      END IF
C
C                       ** REPLACE MISSING DATA FROM CLIMATOLOGY
C
      IF ( Pwater.LT.0.0 ) Pwater = H2oamt
      IF ( Ozone.LT.0.0 ) Ozone  = O3amt
C
      RETURN

      END

C**********************************************************************C
      SUBROUTINE  Chkset( Unit16,
     &                    Smumin, Smumax, Solcon, Matthe, Satalb,
     &                    Maxflx, Maxotr, Maxfpr, Maxint, Mxfpr, Mxint,
     &                    Nflux, Nother, Nfpar, Nintv, Iw1, Iw2 )
C**********************************************************************C
C
C       Checks dimensions of internal arrays and sets certain variables.
C
C  Nflux       : no. of fluxes
C  Nfpar       : no. of flux-type parameters
C  Nintv       : no. of spectral intervals for output
C  Nother      : no. of non-flux type parameters
C
C  Called by- Sasrab
C  Calls- Wrtbad, Wrtdim, Errmsg
C+---------------------------------------------------------------------+
C
C     .. Scalar Arguments ..
      INTEGER Unit16
      LOGICAL  Matthe, Satalb
      INTEGER  Maxflx, Maxfpr, Maxint, Maxotr, Mxint, Mxfpr, Nflux,
     &         Nfpar, Nintv,   Nother
      REAL  Smumax, Smumin, Solcon
C     ..
C     .. Array Arguments ..
      INTEGER  Iw1( * ), Iw2( * )
C     ..
C     .. Local Scalars ..
      LOGICAL  Dimerr, Inperr
C     ..
C     .. External Functions ..
      LOGICAL  Wrtbad, Wrtdim
      EXTERNAL  Wrtbad, Wrtdim
C     ..
C     .. External Subroutines ..
      EXTERNAL  Errmsg
C     ..
C
      Inperr = .FALSE.
      Dimerr = .FALSE.
C
      IF ( Smumin.LT.0.0 ) Inperr = Wrtbad( Unit16, 'SMUMIN' )
      IF ( Smumax.LT.Smumin .OR. Smumax.GT. 1.0 )
     &  Inperr = Wrtbad( Unit16,'SMUMAX' )
      IF ( Solcon.LE.0.0 ) Inperr = Wrtbad(Unit16, 'SOLCON' )
      IF ( Matthe .AND. Satalb ) Inperr=Wrtbad( Unit16,'MATTHE/SATALB')

      IF ( Mxfpr.LT.Nfpar ) Dimerr = Wrtdim( Unit16,'Mxfpr', Nfpar )
      IF ( Mxint.LT.Nintv ) Dimerr = Wrtdim( Unit16,'Mxint', Nintv )
      IF ( Maxfpr.LT.Nfpar ) Dimerr = Wrtdim(Unit16, 'Maxfpr', Nfpar )
      IF ( Maxint.LT.Nintv ) Dimerr = Wrtdim( Unit16,'Maxint', Nintv )
      IF ( Maxflx.LT.Nflux ) Dimerr = Wrtdim( Unit16,'Maxflx', Nflux )
      IF ( Maxotr.LT.Nother ) Dimerr = Wrtdim(Unit16, 'Maxotr', Nother )

      IF ( Inperr .OR .Dimerr ) CALL Errmsg
     &   ( Unit16,'SASRAB--INPUT AND/OR DIMENSION ERRORS', .TRUE. )


      Nflux  = 2 * Nfpar * Nintv - 3
      Nother = 3

C                       ** Set up starting and ending indexes for
C                       ** output spectral intervals (1=sw, 2=visible,
C                       ** 3=nir)

      Iw1( 1 ) = 1
      Iw2( 1 ) = 5
      Iw1( 2 ) = 2
      Iw2( 2 ) = 4
      Iw1( 3 ) = 5
      Iw2( 3 ) = 5

      RETURN

      END


C**********************************************************************C
      SUBROUTINE  Fluxes( Unit16,
     &                    Trndir, Trndif, Reflec, Sphtrn, Sphref,
     &                    Solcon, Rmu, Rsun, Nwlint, Albdir, Albdif,
     &                    Nintv, Maxpar, Iw1, Iw2, Flux )
C**********************************************************************C
C
C     Computes up- and down-fluxes at the top and bottom of an
C     atmosphere above a surface with non-zero albedo.
C
C
C     I N P U T   V A R I A B L E S:
C
C     Albdif(Iw)  : Iw=1 to Nwlint; diffuse surface albedo
C     Albdir(Iw)  : Iw=1 to Nwlint; direct surface albedo
C     Iw1(In)     : In=1 to Nintv; starting index for spectral sum
C     Iw2(In)     : In=1 to Nintv; ending index for spectral sum
C     Maxpar      : first dimension of array -Flux-
C     Nintv       : no. of spectral intervals for output
C     Nwlint      : no. of spectral intervals in refl/trans tables
C     Reflec(Iw)  : Iw=1 to Nwlint; diffuse toa atmospheric reflectance
C                      (does not include surface effect)
C     Rmu         : cosine of solar zenith angle
C     Rsun        : earth-sun distance correction Fr
C     Solcon(Iw)  : Iw=1 to Nwlint; toa solar irradiance at mean sun-
C                      earth distance
C     Sphref(Iw)  : Iw=1 to Nwlint; spherical atmospheric reflectance
C                      (does not include surface effect)
C     Sphtrn(Iw)  : Iw=1 to Nwlint; spherical atmospheric transmittance
C                      (does not include surface effect)
C     Trndif(Iw)  : Iw=1 to Nwlint; diffuse atmospheric transmittance
C                      (does not include surface effect)
C     Trndir(Iw)  : Iw=1 to Nwlint; direct atmospheric transmittance
C                      (does not include surface effect)
C
C     O U T P U T   V A R I A B L E S:
C
C     Flux(Ip,In) : Ip=1 to 5, In=1 to Nintv, flux (W/m sq) integrated
C                      to output spectral intervals
C                      Ip=1, top of atmosphere down-flux
C                      Ip=2, top of atmosphere up-flux
C                      Ip=3, surface down-flux
C                      Ip=4, surface up-flux
C                      Ip=5, diffuse surface down-flux
C                      In=1, shortwave
C                      In=2, visible (PAR)
C                      In=3, near infrared
C
C
C     I N T E R N A L   V A R I A B L E S:
C
C     Dimerr : logical flag, true if dimension error is detected
C     Dirpar : direct par (W/m sq)
C     First  : true on first entry, false thereafter
C     Rp(Iw) : Iw=1 to Nwlint; diffuse reflectance including the effect
C                of surface
C     Sddf(Iw): Iw=1 to Nwlint; spectral surface diffuse downward flux
C     Sdf(Iw) : Iw=1 to Nwlint; spectral surface downward flux
C     Suf(Iw) : Iw=1 to Nwlint; spectral surface upward flux
C     Tdf(Iw) : Iw=1 to Nwlint; spectral TOA downward flux
C     Tp(Iw)  : Iw=1 to Nwlint; diffuse transmittance including the
C                 effect of surface
C     Tuf(Iw) : Iw=1 to Nwlint; spectral TOA upward flux
C
C     LOCAL SYMBOLIC DIMENSION:
C
C     Maxwvl : max no. of spectral intervals
C
C     Called by- Insest
C     Calls- Wrtdim, Errmsg, Sums
C+---------------------------------------------------------------------+
C
C     .. Parameters ..
      INTEGER Unit16
      INTEGER  Maxwvl
      PARAMETER  ( Maxwvl=5 )
C     ..
C     .. Scalar Arguments ..
      INTEGER  Maxpar, Nintv, Nwlint
      REAL  Rmu, Rsun
C     ..
C     .. Array Arguments ..
      INTEGER  Iw1( * ), Iw2( * )
      REAL  Albdif( * ), Albdir( * ), Flux( Maxpar, * ), Reflec( * ),
     &      Solcon( * ), Sphref( * ), Sphtrn( * ), Trndif( * ),
     &      Trndir( * )
C     ..
C     .. Local Scalars ..
      LOGICAL  Dimerr, First
      INTEGER  In, Iw
C     ..
C     .. Local Arrays ..
      REAL  Fr( Maxwvl ), Rp( Maxwvl ), Sddf( Maxwvl ), Sdf( Maxwvl ),
     &      Suf( Maxwvl ), Tdf( Maxwvl ), Tp( Maxwvl ), Tuf( Maxwvl )
C     ..
C     .. External Functions ..
      LOGICAL  Wrtdim
      REAL  Sums
      EXTERNAL  Wrtdim, Sums
C     ..
C     .. External Subroutines ..
      EXTERNAL  Errmsg
C     ..
C     .. Save statement ..
      SAVE  First
C     ..
C     .. Data statements ..
      DATA  First / .TRUE. /
C     ..
C
C
      IF ( First ) THEN
C                       ** Check local dimension
         First  = .FALSE.
         Dimerr = .FALSE.
         IF ( Nwlint.GT.Maxwvl ) Dimerr = Wrtdim( Unit16,'Maxwvl'
     &      , Nwlint )
         IF ( Dimerr ) CALL Errmsg(Unit16, 'FLUXES--DIM ERROR', .TRUE. )

      END IF

C                       ** Calculate refl/trans functions for the
C                       ** 'planetary problem' (add surface reflection)
C
      DO 10 Iw = 1, Nwlint
         Fr( Iw ) = 1. / ( 1.-Albdif(Iw) *Sphref(Iw) ) *
     &              ( Albdir(Iw) *Trndir(Iw)+Albdif(Iw) *Trndif(Iw) )
         Rp( Iw ) = Reflec( Iw ) + Fr( Iw ) * Sphtrn( Iw )
         Tp( Iw ) = Trndif( Iw ) + Fr( Iw ) * Sphref( Iw )
10    CONTINUE

C                       ** Calculate spectral fluxes
C
      DO 20 Iw = 1, Nwlint
         Tdf( Iw ) = Rsun * Rmu * Solcon( Iw )
         Tuf( Iw ) = Rsun * Rmu * Solcon( Iw ) * Rp( Iw )
         Sdf( Iw ) = Rsun * Rmu * Solcon( Iw ) * ( Trndir(Iw)+Tp(Iw) )
         Suf( Iw ) = Rsun * Rmu * Solcon( Iw ) *
     &               ( Trndir(Iw) *Albdir(Iw)+Tp(Iw) *Albdif(Iw) )
         Sddf( Iw ) = Rsun * Rmu * Solcon( Iw ) * Tp( Iw )
20    CONTINUE

C                      ** Integrate spectrally to get fluxes in output
C                      ** spectral intervals, and save the different
C                      ** type of fluxes in a singel array
C
      DO 30 In = 1, Nintv
         Flux( 1, In ) = Sums(Unit16,Tdf, Iw1(In), Iw2(In), Nwlint)
         Flux( 2, In ) = Sums(Unit16,Tuf, Iw1(In), Iw2(In), Nwlint)
         Flux( 3, In ) = Sums(Unit16,Sdf, Iw1(In), Iw2(In), Nwlint)
         Flux( 4, In ) = Sums(Unit16,Suf, Iw1(In), Iw2(In), Nwlint)
         Flux( 5, In ) = Sums(Unit16,Sddf, Iw1(In), Iw2(In), Nwlint)
30    CONTINUE

      RETURN

      END

C**********************************************************************C
      SUBROUTINE  Getalb( Unit16,
     &                    Trndir, Trndif, Reflec, Sphtrn, Sphref,
     &                    Tauaer, Tauclm, Solirr, Toaalb, Diralb,
     &                    Difalb, Nwlint, Ntauar, Mxwavl, Noretr )
C**********************************************************************C
C
C     Version 3. (May 1995)
C
C     Estimates spectral surface albedos from the broadband toa albedo
C     derived from the clear sky composite radiance, assuming aerosol
C     optical depth corresponds to climatology, and the surface albedo
C     has the same spectral dependence as that of the reference model.
C     The surface albedo is estimated by solving the equation below for
C     Scldir and Scldif.
C
C     TOAALB = SUM (FROM I=1 TO NWLINT) OF { REFLE(I) * Solirr(I) +
C                       1 / ( 1 - SCLDIF * DIFALB(I) * SPHRE(I) ) *
C                                 ( SCLDIR * DIRALB(I) * DIRTR(I) +
C                                 SCLDIF * DIFALB(I) * DIFTR(I) ) *
C                                          SPHTR(I) * Solirr(I) } /
C              SUM (FROM I=1 TO NWLINT) OF Solirr(I)
C
C     See the list below for the variables in the equation.
C 
C     It is assumed that Scldif = Scldir.
C
C     NOTE: If albedo cannot be retrieved for some reason (bright
C           model atmosphere, bright satellite albedo, etc.) the
C           reference albedos are returned.
C           Both direct and diffuse components of the surface
C           albedo are scaled with the same amount, even for ocean
C           (in the previous versions diffuse albedo of water was
C           not scaled).
C
C
C     ROUTINE CALLED:  LOCATE
C
C
C     I N P U T   V A R I A B L E S:
C
C     DIFALB(IL)    : IL=1 TO NWLINT; REFERENCE DIFFUSE SURFACE ALBEDO
C     DIRALB(IL)    : IL=1 TO NWLINT; REFERENCE DIRECT SURFACE ALBEDO
C     MXWAVL        : MAX NO. OF SPECTRAL INTERVALS IN CALLING ROUTINE
C     NTAUAR        : NO. OF AEROSOL OPTICAL DEPTH VALUES
C     NWLINT        : NO. OF SPECTRAL INTERVALS
C     REFLEC(IL,IT) : IL=1 TO NWLINT, IT=1 TO NTAUAR; CLEAR SKY DIFFUSE
C                     REFLECTANCE OF THE ATMOSPHERE ALONE (NO EFFECT OF
C                     THE SURFACE IS INCLUDED)
C     Solirr(IL)    : IL=1 TO NWLINT; TOA SPECTRAL SOLAR IRRAD (W/M SQ)
C     SPHREF(IL,IT) : IL=1 TO NWLINT, IT=1 TO NTAUAR; CLEAR SKY SPHERIC
C                     REFLECTANCE OF THE ATMOSPHERE ALONE
C     SPHTRN(IL,IT) : IL=1 TO NWLINT, IT=1 TO NTAUAR; CLEAR SKY SPHERIC
C                     TRANSMITTANCE OF THE ATMOSPHERE ALONE
C     TAUAER(IT)    : IT=1 TO NTAUA; TABLE VALUES OF AEROSOL OPTICAL
C                     DEPTH
C     TAUCLM        : CLIMATOLOGICAL AEROSOL OPTICAL DEPTH AT 0.55 MICR
C     TOAALB        : BROADBAND ALBEDO FROM CLEAR SKY COMPOSITE RADIANCE
C     TRNDIF(IL,IT) : IL=1 TO NWLINT, IT=1 TO NTAUAR; CLEAR SKY DIFFUSE
C                     TRANSMITTANCE OF THE ATMOSPHERE ALONE
C     TRNDIR(IL,IT) : IL=1 TO NWLINT, IT=1 TO NTAUAR; CLEAR SKY DIRECT
C                     TRANSMITTANCE OF THE ATMOSPHERE ALONE
C
C
C     O U T P U T  V A R I A B L E S:
C
C     DIFALB(IL)  : IL=1 TO NWLINT; REFERENCE DIFFUSE SURFACE ALBEDO
C                   SCALED BY A FACTOR DERIVED FROM SATELLITE
C     DIRALB(IL)  : IL=1 TO NWLINT; REFERENCE DIRECT SURFACE ALBEDO
C                   SCALED BY A FACTOR DERIVED FROM SATELLITE
C     NORETR      : TRUE IF SURFACE ALBEDO IS NOT RETRIEVED
C
C
C     I N T E R N A L   V A R I A B L E S:
C
C     DIFTR(IL) : IL=1 TO NWLINT, DIFFUSE ATMOSPHERIC TRANSMITTANCE FOR
C                 AEROSOL OPTICAL DEPTH -TAUCLM-. (WITHOUT THE EFFECT OF
C                 THE SURFACE)
C     DIRTR(IL) : IL=1 TO NWLINT, AS -DIFTR-, BUT DIRECT TRANSMITTANCE
C     REFLE(IL) : IL=1 TO NWLINT, DIFFUSE ATMOSPHERIC REFLECTANCE FOR
C                 AEROSOL OPTICAL DEPTH -TAUCLM-. (WITHOUT THE EFFECT OF
C                 THE SURFACE)
C     Refsw     : shortwave TOA albedo of the atmosphere
C     SCLDIF    : SCALING FACTOR DERIVED FROM THE SATELLITE MEASUREMENT.
C                 SPECTRAL DIFFUSE ALBEDOS OF THE REFERENCE MODEL ARE
C                 MULTIPLIED BY -SCLDIF-
C     SCLDIR    : AS -SCLDIF-, BUT FOR DIRECT ALBEDOS
C     Solsw     : Shortwave solar irradiance
C     SPHRE(IL) : IL=1 TO NWLINT, SPHERICAL REFLECTANCE FOR AEROSOL
C                 OPTICAL DEPTH -TAUCLM-. (SURFACE IS NOT INCLUDED)
C     SPHTR(IL) : IL=1 TO NWLINT, AS -SPHRE-, BUT TRANSMITTANCE
C
C
C     LOCAL SYMBOLIC DIMENSION:
C
C     MAXWVL : MAX NO. OF SPECTRAL INTERVALS
C
C     Called by- Insest
C     Calls- Wrtdim, Errmsg, Locate
C+---------------------------------------------------------------------+
C
C     .. Parameters ..
      INTEGER Unit16
      INTEGER  Maxwvl
      PARAMETER  ( Maxwvl=5 )
C     ..
C     .. Scalar Arguments ..
      LOGICAL  Noretr
      INTEGER  Mxwavl, Ntauar, Nwlint
      REAL  Tauclm, Toaalb
C     ..
C     .. Array Arguments ..
      REAL  Difalb( * ), Diralb( * ), Reflec( Mxwavl, * ), Solirr( * ),
     &      Sphref( Mxwavl, * ), Sphtrn( Mxwavl, * ), Tauaer( * ),
     &      Trndif( Mxwavl, * ), Trndir( Mxwavl, * )
C     ..
C     .. Local Scalars ..
      LOGICAL  Dimerr, First
      INTEGER  I, Il, It1, It2, L
      REAL  Ftau, Refsw, Scale, Solsw, Term1, Term2
C     ..
C     .. Local Arrays ..
      REAL  Diftr( Maxwvl ), Dirtr( Maxwvl ), Refle( Maxwvl ),
     &      Sphre( Maxwvl ), Sphtr( Maxwvl )
C     ..
C     .. External Functions ..
      LOGICAL  Wrtdim
      EXTERNAL  Wrtdim
C     ..
C     .. External Subroutines ..
      EXTERNAL  Errmsg, Locate
C     ..
C     .. Save statement ..
      SAVE  First
C     ..
C     .. Data statements ..
      DATA  First / .TRUE. /
C     ..

      Noretr = .FALSE.

      IF ( First ) THEN
C
C                       ** Check local dimension
C
         First  = .FALSE.
         Dimerr = .FALSE.
         IF ( Nwlint.GT.Maxwvl) Dimerr=Wrtdim(Unit16, 'MAXWVL',Nwlint)
         IF ( Dimerr ) CALL Errmsg(Unit16, 'GETALB -- DIMENSION ERROR',
     &                              .TRUE. )

      END IF
C
C                       ** Locate-TAUCLM- in array -TAUAER-
C
      CALL Locate( Tauclm, Tauaer, Ntauar, It1, It2, Ftau )
C
C                       ** Interpolate reflectances and transmittances
C                       ** for -TAUCLM-
C
      DO 10 L = 1, Nwlint
         Refle( L ) = Ftau * Reflec( L, It1 ) +
     &                ( 1.-Ftau ) * Reflec( L, It2 )
         Dirtr( L ) = Ftau * Trndir( L, It1 ) +
     &                ( 1.-Ftau ) * Trndir( L, It2 )
         Diftr( L ) = Ftau * Trndif( L, It1 ) +
     &                ( 1.-Ftau ) * Trndif( L, It2 )
         Sphre( L ) = Ftau * Sphref( L, It1 ) +
     &                ( 1.-Ftau ) * Sphref( L, It2 )
         Sphtr( L ) = Ftau * Sphtrn( L, It1 ) +
     &                ( 1.-Ftau ) * Sphtrn( L, It2 )
10    CONTINUE
C
C                      ** Calculate shortwave TOA albedo of atmosphere
C                      ** without surface
      Solsw  = 0.0
      Refsw  = 0.0
      DO 20 Il = 1, Nwlint
         Solsw  = Solsw + Solirr( Il )
         Refsw  = Refsw + Refle( Il ) * Solirr( Il )
20    CONTINUE
      Refsw  = Refsw / Solsw
C
C                      ** Check if model atmosphere is brighter than the
C                      ** satellite observed TOA albedo. In this case
C                      ** surface albedo cannot be retrieved, so use
C                      ** the reference values.
C
      IF ( Refsw.GE.Toaalb ) RETURN
C
C                      ** Calculate scaling factor
C
      Term1  = 0.0
      Term2  = 0.0
      DO 30 Il = 1, Nwlint
         Term1  = Term1 + Solirr( Il ) * ( Toaalb-Refle(Il) ) *
     &            Difalb( Il ) * Sphre( Il )
         Term2  = Term2 + Solirr( Il ) *
     &            ( Diralb(Il) *Dirtr(Il)+Difalb(Il) *Diftr(Il) ) *
     &            Sphtr( Il )
30    CONTINUE
      Scale  = ( Toaalb-Refsw ) * Solsw / ( Term1+Term2 )

C
C                      ** Make sure scaling of spectral albedos does
C                      ** not lead to non-physical (less than or equal
C                      ** to zero, greater than or equal to one) values.
C                      ** In case it does, use unscaled reference
C                      ** values instead.
C
      DO 40 I = 1, Nwlint

         IF ( Scale*Diralb(I).LE.0.0 .OR. Scale*Difalb(I).LE.0.0 .OR.
     &        Scale*Diralb(I).GE.1.0 .OR. Scale*Difalb(I).GE.1.0 ) THEN
            RETURN

         END IF

40    CONTINUE
C
C                      ** Scale reference spectral albedos
C
      DO 50 I = 1, Nwlint
         Diralb( I ) = Scale * Diralb( I )
         Difalb( I ) = Scale * Difalb( I )
50    CONTINUE

      RETURN

      END

C**********************************************************************C
      SUBROUTINE  Insest( Unit16, 
     &                    Cdfrac, Clrref, Cldref, Cmpref, Pwater,
     &                    Ozone, Isatid, Snowfr, Smu, Satmu, Relaz,
     &                    Glint, Solirr, Aertab, Cldtab, Waveln, Jday,
     &                    Glat, Glon, Time, Ntaua, Ntauc, Nwlint,
     &                    Itable, Isrmod, Isrtyp, Tauclm, Sundis,
     &                    Declin, Nontob, Matthe, Satalb, Nobdrc, Prnt,
     &                    Extrap, Smumin, Smumax, Ilat, Ilon, Mxntau,
     &                    Mxnsun, Mxnwvl, Mxnh2o, Mxno3, Mxntab, Mxnfpr,
     &                    Mxnint, Misval, Nflux, Nintv, Iw1, Iw2,
     &                    Fnbdr, Flxall, Flxclr, Tauclr, Taucld,
     &                    Cmpbba, Clrbba, Cldbba)
C**********************************************************************C
C
C     1) Estimates instantaneous values of:
C         a) clear sky toa and surface fluxes (both up and down), as
C            well as aerosol optical depth, from toa albedo derived from
C            the average of the clear pixel radiances;
C         b) cloudy sky toa and surface fluxes (both up and down), as
C            well as cloud optical depth, from toa albedo derived from
C            the average of the cloudy pixel radiances.
C     2) Computes an average of the fluxes by weighting their clear and
C        cloudy values by the number of the clear and cloudy pixels in
C        the cell.
C     3) If -Prnt=true-, prints relevent input parameters, toa
C        narrowband reflectance and broadband albedo from clear and
C        cloudy sky satellite measurements, cloud cover, estimated
C        values of the upward, downward and net fluxes and albedos at
C        the top of the atmosphere and at the surface, separately for
C        clear, cloudy and mixed skies. In addition, estimated values
C        of the aerosol optical depth for clear sky and cloud optical
C        depth for cloudy sky are also reported.
C
C     NOTE: For a better representation of clear-sky fluxes, clear-sky
C           flux is obtained from the clear-sky composite radiance
C           whenever the clear-sky radiance is missing.

C     I N P U T  V A R I A B L E S :
C
C     Aertab(It)  : It=1 to Ntaua; table values of aerosol optical depth
C     Cdfrac      : cloud cover fraction
C     Cldref      : mean narrowband reflectance from cloudy pixels
C     Cldtab(It)  : It=1 to Ntauc; table values of cloud optical depths
C     Clrref      : mean narrowband reflectance from clear pixels
C     Cmpref      : narrowband reflectance from clear sky composite
C     Declin      : solar declination in degrees
C     Extrap      : if true, fluxes are extrapolated when no match of
C                   the toa albedos is found. If false, fluxes are not
C                   retrieved
C     Fnbdr       : path and name containing TOA bdrf data
C     Glat        : latitude of the gridpoint (-90 to +90 deg)
C     Glon        : longitude of the gridpoint (-180 to +180 deg)
C     Ilat        : latitude index of gridpoint in maps of matthes
C     Ilon        : same as above, but longitude index
C     Isatid      : satellite id code
C     Isrmod      : surface model id
C     Isrtyp      : integer number indicating the surface type
C                   ( 1-water, 2-vegetation, 3-desert, 4-snow/ice )
C     Itable      : id of refl-trans table to be used
C     Iw1(In)     : In=1 to Nintv; starting index for spectral sum
C     Iw2(In)     : In=1 to Nintv; ending index for spectral sum
C     Jday        : day of the year
C     Matthe      : if true, reference surface albedo is scaled to the
C                   seasonal surface albedo data of matthews
C     Misval      : value representing missing data
C     Mxnfpr      : max no. of flux parameters
C     Mxnh2o      : max no. of water vapor amounts in common block
C                   *tables* in subroutine *reftra*
C     Mxnint      : max no. of output spectral intervals
C     Mxno3       : max no. of ozone amounts in common block
C                   *tables* in subroutine *reftra*
C     Mxnsun      : max no. of solar zenith angles in common block
C                   *tables* in subroutine *reftra*
C     Mxntab      : max no. of refl/trans tables  in common block
C                   *tables* in subroutine *reftra*
C     Mxntau      : max no. of optical depths in common block
C                   *tables* in subroutine *reftra*
C     Mxnwvl      : max no. of spectral intervals in common block
C                   *tables* in subroutine *reftra*
C     Nflux       : no. of flux parameters
C     Nintv       : no. of spectral intervals for output
C     Nobdrc      : if true, no bdr correction is done
C     Nontob      : if true, no narrow-to-broad conversion of the toa
C                   albedos is performed
C     Ntaua       : no. of aerosol optical depth values in table
C     Ntauc       : no. of cloud optical depth values in table
C     Nwlint      : no. of spectral intervals
C     Ozone       : column amount of ozone in atm-cm
C     Prnt        : if true, prints instantaneous values of fluxes and
C                   optical depths
C     Pwater      : amount of precipitable water (cm)
C     Relaz       : relative azimuth angle
C     Satalb      : if true, surface albedo is retrieved from satellite
C     Satmu       : cosine of satellite zenith angle
C     Smu         : cosine of the solar zenith angle
C     Smumax      : max value of solar zenith cosine to be considered
C     Smumin      : min value of solar zenith cosine to be considered
C     Snowfr      : snow/ice cover fraction (range: 0-1)
C     Solirr(Il)  : Il=1 to Nwlint(Ib); spectral solar constant from
C                   ref/trans table
C     Sundis      : sun-earth distance factor for the day
C     Tauclm      : climatological aerosol optical depth
C     Time        : local mean time in hours
C     Waveln(Il)  : Il=1 to Nwlint(Ib)+1; boundaries of spectral
C                   intervals in ref/tran tables (microns)
C
C     O U T P U T  V A R I A B L E S :
C
C     Flxall(Ip,Iv) : Ip=1 to Nflux, Iv=1 to Nintv; all-sky fluxes
C     Flxclr(Ip,Iv) : Ip=1 to Nflux, Iv=1 to Nintv; clear-sky fluxes
C     Taucld        : estimated cloud optical depth at 0.55 microns
C     Tauclr        : estimated aerosol optical depth at 0.55 microns
C
C        Note: Flxall and Flxclr are in W/m2, the indexes Ip and Iv are:
C                      Ip=1, top of atmosphere down-flux
C                      Ip=2, top of atmosphere up-flux
C                      Ip=3, surface down-flux
C                      Ip=4, surface up-flux
C                      Ip=5, diffuse surface down-flux
C                         Iv=1, shortwave
C                         Iv=2, visible (PAR)
C                         Iv=3, near infrared
C
C
C   I N T E R N A L  V A R I A B L E S:
C
C   Ccover         : cloud cover (in percent)
C   Cldbba         : broadband toa albedo derived from -Cldref-
C   Cldref         : narrowband toa reflectance from cloudy pixels
C   Clmtau         : climatological aerosol optical depth
C   Clrbba         : broadband toa albedo derived from -Clrref/Cmpref-
C   Clrref         : narrowband toa reflectance from clear pixels
C   Cmpbba         : broadband toa albedo derived from -Cmpref-
C   Cmpref         : narrowband toa reflectance from clear sky composite
C   Difalb(Il)     : Il=1 to Nwlint; spectral diffuse surface refl
C   Diralb(Il)     : Il=1 to Nwlint; spectral direct surface refl
C   Dimerr         : true, if dimension error is detected
C   Diftra(Il,It)  : Il=1 to Nwlint, It=1 to Ntaua; clear sky diffuse
C                    transmittance interpolated from refl/trans table
C                    to -Smu-, -Pwater- and -Ozone-.
C   Diftrc(Il,It)  : as above, but cloudy sky diffuse transmittance
C   Dirtra(Il,It)  : as above, but clear sky direct transmittance
C   Dirtrc(Il,It)  : as above, but cloudy sky direct transmittance
C   Isurft         : surface type id
C   Jdprev         : value of -Jday- in the previous entry
C   Newday         : true, if new day
C   Nomcld         : true if match between the toa measured and modeled
C                    albedos was not found for the cloudy scene
C   Nomclr         : true if match between the toa measured and modeled
C                    albedos was not found for the clear scene
C   Noretr         : true if retrieval of fluxes is not possible
C   Polnit         : true if polar night, false otherwise
C   Reflea(Il,It)  : as above, but clear sky diffuse reflectance
C   Reflec(Il,It)  : as above, but cloudy sky diffuse reflectance
C   Sphrea(Il,It)  : as above, but clear sky spherical reflectance
C   Sphrec(Il,It)  : as above, but cloudy sky spherical reflectance
C   Sphtra(Il,It)  : as above, but clear sky spherical transmittance
C   Sphtrc(Il,It)  : as above, but cloudy sky spherical transmittance
C
C     LOCAL SYMBOLIC DIMENSIONS:
C
C     Maxflx  =   max no. of flux parameters
C     Maxint  =   max no. of spectral intervals for output
C     Maxtau  =   max no. of optical depths
C     Maxwvl  =   max no. of spectral intervals
C
C   Called by- Sasrab
C   Calls- Wrtdim, Errmsg, Intpol, Refsur, Conver, Getalb, Matalb,
C          Matcha, Fluxes
C+---------------------------------------------------------------------+
C
C     .. Parameters ..
      INTEGER Unit16
      INTEGER Maxtau, Maxwvl, Maxflx, Maxint
      PARAMETER  ( Maxtau=15, Maxwvl=5, Maxflx=5, Maxint=3 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER  Fnbdr * ( * )
      LOGICAL  Extrap, Matthe, Nobdrc, Nontob, Prnt, Satalb
      INTEGER  Ilat, Ilon, Isatid, Isrmod, Isrtyp, Itable,
     &         Jday, Mxnfpr, Mxnh2o, Mxnint, Mxno3, Mxnsun, Mxntab,
     &         Mxntau, Mxnwvl, Nflux, Nintv, Ntaua, Ntauc, Nwlint
      REAL Glint
      REAL  Cdfrac, Cldref, Clrref, Cmpref, Declin, Glat, Glon,
     &      Misval, Ozone, Pwater, Relaz, Satmu, Smu, Smumax, Smumin,
     &      Snowfr, Sundis, Taucld, Tauclm, Tauclr, Time
C     ..
C     .. Array Arguments ..
      INTEGER  Iw1( * ), Iw2( * )
      REAL  Cldtab( * ), Aertab( * ), Solirr( * ), Waveln( * ),
     &      Flxall( Mxnfpr, Mxnint ), Flxclr( Mxnfpr, Mxnint )
C     ..
C     .. Local Scalars ..
      LOGICAL  Badinp, Dimerr, Docld, Doclr, First, Newday, Nomcld,
     &         Nomclr, Noretr, North, Pole, Polnit, Winter
      INTEGER  In, Ip, Isurft, Jdprev
      REAL  Ccover, Cldbba,  Clmtau, Clrbba, Cmpbba
C     ..
C     .. Local Arrays ..
      REAL  Difalb( Maxwvl ), Difref( Maxwvl ),
     &      Diftra( Maxwvl, Maxtau ), Diftrc( Maxwvl, Maxtau ),
     &      Diftrn( Maxwvl ), Diralb( Maxwvl ),
     &      Dirtra( Maxwvl, Maxtau ), Dirtrc( Maxwvl, Maxtau ),
     &      Dirtrn( Maxwvl ), Flxcld( Maxflx, Maxint ),
     &      Reflea( Maxwvl, Maxtau ), Reflec( Maxwvl, Maxtau ),
     &      Sphrea( Maxwvl, Maxtau ), Sphrec( Maxwvl, Maxtau ),
     &      Sphref( Maxwvl ), Sphtra( Maxwvl, Maxtau ),
     &      Sphtrc( Maxwvl, Maxtau ), Sphtrn( Maxwvl )
C     ..
      INTEGER Get_Lun
      INTEGER Unit18

C     .. External Functions ..
      LOGICAL  Wrtdim
      EXTERNAL  Wrtdim
C     ..
C     .. External Subroutines ..
      EXTERNAL  Conver, Errmsg, Fluxes, Getalb, Intpol, Matalb, Matcha,
     &          Refsur
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC  Max
C     ..
C     .. Save statement ..
      SAVE  First, Jdprev, Diftra, Dirtra, Reflea, Sphrea, Sphtra,
     &      Diftrc, Dirtrc, Reflec, Sphrec, Sphtrc
C     ..
C     .. Data statements ..
      DATA  First / .TRUE. / , Jdprev / -1 /
C     ..

C                       ** CHECK LOCAL DIMENSIONS
      IF ( First ) THEN
         First  = .FALSE.
         Dimerr = .FALSE.
         IF (Maxtau.LT.Max(Ntaua,Ntauc)) Dimerr=Wrtdim(unit16,'Maxtau',
     &        Max(Ntaua,Ntauc) )
         IF ( Maxwvl.LT.Nwlint )Dimerr=Wrtdim(Unit16,'Maxwvl', Nwlint )
         IF ( Maxflx.LT.Nflux ) Dimerr = Wrtdim(unit16,'Maxflx', Nflux )
         IF ( Maxint.LT.Nintv ) Dimerr = Wrtdim(unit16,'Maxint', Nintv )
         IF ( Dimerr ) CALL Errmsg(Unit16, 'INSEST--DIMENSION ERROR(S)',
     &                              .TRUE. )
         IF ( Prnt ) THEN
            WRITE (Unit16, FMT=9000 )
            WRITE (Unit16, FMT=9010 )

         END IF

      END IF
C
C                       ** INITIALIZE RELEVANT VARIABLES
      Noretr = .FALSE.
      Nomclr = .FALSE.
      Nomcld = .FALSE.
      Tauclr = Misval
      Taucld = Misval
      DO 20 In = 1, Nintv
         DO 10 Ip = 1, Nflux
            Flxall( Ip, In ) = Misval
            Flxclr( Ip, In ) = Misval
10       CONTINUE
20    CONTINUE
      Cmpbba = Misval
      Clrbba = Misval
      Cldbba = Misval
C
C                       ** SEE IF POLAR NIGHT
C
      Polnit = .FALSE.
      IF ( Declin.GT.0.0 ) Polnit = Glat .LT. ( Declin-90.0 )
      IF ( Declin.LT.0.0 ) Polnit = Glat .GT. ( Declin+90.0 )

      IF ( Polnit ) THEN
C                       ** Set fluxes to zero
         DO 40 In = 1, Nintv
            DO 30 Ip = 1, Nflux
               Flxall( Ip, In ) = 0.0
               Flxclr( Ip, In ) = 0.0
30          CONTINUE
40       CONTINUE

         RETURN

      END IF
C
C                       ** DO NOT RETRIEVE FLUXES/OPTICAL DEPTHS IF:
C                       ** 1) THERE ARE NO OBSERVATIONS (BOTH -Clrref-
C                       **    AND -Cldref- ARE MISSING);
C                       ** 2) SUN IS TOO CLOSE TO HORIZON OR NIGHTTIME;
C                       ** 3) THERE IS NO CLEAR COMPOSITE
C
      IF ( Clrref == Misval .AND. Cldref == Misval ) RETURN
      IF ( Smu.LT.Smumin .OR. Smu.GT.Smumax ) RETURN
      IF ( Cmpref == Misval ) RETURN

C            ** Make sure all radiances and number of pixels are valid

      Badinp = .FALSE.
      IF ( Clrref.LE.0.0 .AND. Cdfrac.NE.1. ) Badinp = .TRUE.
      IF ( Cldref.LE.0.0 .AND. Cdfrac.NE.0. ) Badinp = .TRUE.
      IF ( Satalb .AND. Cmpref.LE.0.0 ) Badinp = .TRUE.
      IF ( Badinp ) THEN
         WRITE (Unit16, FMT= * ) Cdfrac, Clrref
         WRITE (Unit16, FMT= * ) Cdfrac, Cldref
         WRITE (Unit16, FMT= * ) Satalb, Cmpref, Time
         CALL Errmsg( Unit16, 'Insest-Bad input data', .FALSE. )
         RETURN

      END IF
C
C                       ** INTERPOLATE TABLE VALUES OF TRANSMISSION
C                       ** AND REFLECTANCE FOR THE ACTUAL VALUES OF
C                       ** SOLAR ZENITH ANGLE, WATER VAPOR AND OZONE
C
      Doclr  = Clrref .NE. Misval .OR. Satalb
      Docld  = Cldref .NE. Misval
      CALL Intpol( Unit16,
     &             Smu, Pwater, Ozone, Itable, Doclr, Docld, Nwlint,
     &             Ntaua, Ntauc, Maxwvl, Mxntau, Mxnsun, Mxnwvl, Mxnh2o,
     &             Mxno3, Mxntab, Dirtra, Diftra, Reflea, Sphtra,
     &             Sphrea, Dirtrc, Diftrc, Reflec, Sphtrc, Sphrec )
C
C                       ** GET SPECTRAL REFLECTANCES OF REFERENCE
C                       ** SURFACE (BRIEGLEB'S MODELS+)
C
      CALL Refsur( Smu, Isrmod, Snowfr, Nwlint, Waveln, Itable, Diralb,
     &             Difalb, Noretr )
      IF ( Noretr ) RETURN

      Isurft = Isrtyp
      Clmtau = Tauclm
C
C                       ** SET SURFACE TYPE -ISURFT- (AND -CLMTAU-)
C                       ** TO SNOW IF SNOW IS PRESENT
C
      IF ( Snowfr.GT.0.0 ) THEN
         Isurft = 4
         Pole   = Glat .GE. 66.0 .OR. Glat .LE. -66.0
         North  = Glat .GT. 0.0
         Winter = Jday .LT. 90 .OR. Jday .GT. 304
         IF ( Pole ) Clmtau = 0.05
         IF ( North .AND. Pole .AND. Winter ) Clmtau = 0.2

      END IF


      IF ( Satalb ) THEN
C                       ** ESTIMATE SPECTRAL SURFACE ALBEDO FROM
C                       ** SATELLITE MEASUREMENT
C
C                                 * GET TOA BROADBAND ALBEDO FROM CLEAR
C                                 * SKY COMPOSITE NARROWBAND REFLECTANCE
C
         CALL Conver( Unit16,
     &                Nontob, Nobdrc, Cmpref, Isurft, .TRUE., Isatid,
     &                Relaz, Satmu, Smu, Cmpbba, Noretr, Fnbdr, Glint )

         IF ( Noretr ) RETURN
C
C                                 * GET SURFACE ALBEDO BY CORRECTING FOR
C                                 * ATMOSPHERIC EFFECTS
C
         CALL Getalb( Unit16, 
     &                Dirtra, Diftra, Reflea, Sphtra, Sphrea, Aertab,
     &                Clmtau, Solirr, Cmpbba, Diralb, Difalb, Nwlint,
     &                Ntaua, Maxwvl, Noretr )
         IF ( Noretr ) RETURN

      ELSE IF ( Matthe ) THEN
C
C                       ** SCALE ALBEDO OF REFERENCE SURFACE-MODEL TO
C                       ** THAT OF MATTHEWS
C
         CALL Matalb( Unit16,
     &                Jday, Ilat, Ilon, Glat, Nwlint, Solirr, Snowfr,
     &                Isrmod, Waveln, Diralb, Difalb, Noretr )
         IF ( Noretr ) RETURN

      END IF

C                       ** Initialize clear fluxes
C
      DO 60 In = 1, Nintv
         DO 50 Ip = 1, Nflux
            Flxclr( Ip, In ) = Misval
50       CONTINUE
60    CONTINUE

      IF ( Clrref == Misval ) THEN
         Tauclr = Misval
         Clrref = Misval
         Clrbba = Misval

      ELSE

C                       ** C L E A R   S C E N E
C
C                       ** CONVERT CLEAR-SKY NARROWBAND REFLECTANCE
C                       ** INTO BROADBAND  PLANETARY ALBEDO -CLRBBA-
C
         CALL Conver( Unit16,
     &                Nontob, Nobdrc, Clrref, Isurft, .TRUE., Isatid,
     &                Relaz, Satmu, Smu, Clrbba, Noretr, Fnbdr, Glint )
         IF ( Noretr ) RETURN

C                       ** MATCH MEASURED AND MODELED ALBEDOS TO GET
C                       ** A E R O S O L  OPTICAL DEPTH AND  C L E A R
C                       ** S K Y  REFL/TRANS FUNCTIONS
C
         CALL Matcha( Unit16,
     &                Ntaua, Aertab, Clrbba, Dirtra, Diftra, Reflea,
     &                Sphtra, Sphrea, Diralb, Difalb, Solirr, Smu,
     &                Sundis, Nwlint, Maxwvl, Clmtau, Extrap, Tauclr,
     &                Dirtrn, Diftrn, Difref, Sphref, Sphtrn, Nomclr )

C                       ** ADD SURFACE REFLECTANCE AND INTEGRATE
C                       ** SPECTRALLY TO OBTAIN CLEAR SKY FLUXES AT
C                       ** THE BOUNDARIES
C
         IF ( Extrap .OR. .NOT.Nomclr ) CALL Fluxes( Unit16,
     &        Dirtrn, Diftrn,
     &        Difref, Sphtrn, Sphref, Solirr, Smu, Sundis, Nwlint,
     &        Diralb, Difalb, Nintv, Mxnfpr, Iw1, Iw2, Flxclr )

      END IF

C                       ** Initialize cloudy fluxes
C
      DO 80 In = 1, Nintv
         DO 70 Ip = 1, Nflux
            Flxcld( Ip, In ) = Misval
70       CONTINUE
80    CONTINUE

      IF ( Cldref == Misval ) THEN
         Taucld = Misval
         Cldref = Misval
         Cldbba = Misval

      ELSE

C                       ** C L O U D Y   S C E N E
C
C                       ** CONVERT CLOUDY-SKY NARROWBAND REFLECTANCE
C                       ** INTO BROADBAND PLANETARY ALBEDO -CLDBBA-
C
         CALL Conver( Unit16,
     &                Nontob, Nobdrc, Cldref, Isurft, .FALSE., Isatid,
     &                Relaz, Satmu, Smu, Cldbba, Noretr, Fnbdr, Glint )
         IF ( Noretr ) RETURN

C                       ** MATCH MEASURED AND MODELED ALBEDOS TO GET
C                       ** C L O U D   OPTICAL DEPTH AND  C L O U D Y
C                       ** S K Y  REFL/TRANS FUNCTIONS
C
         CALL Matcha( Unit16,
     &                Ntauc, Cldtab, Cldbba, Dirtrc, Diftrc, Reflec,
     &                Sphtrc, Sphrec, Diralb, Difalb, Solirr, Smu,
     &                Sundis, Nwlint, Maxwvl, 15.0, Extrap, Taucld,
     &                Dirtrn, Diftrn, Difref, Sphref, Sphtrn, Nomcld )

C                       ** ADD SURFACE REFLECTANCE AND INTEGRATE
C                       ** SPECTRALLY TO OBTAIN CLOUDY SKY FLUXES AT
C                       ** THE BOUNDARIES
C
         IF ( Extrap .OR. .NOT.Nomcld ) CALL Fluxes( Unit16,
     &        Dirtrn, Diftrn,
     &        Difref, Sphtrn, Sphref, Solirr, Smu, Sundis, Nwlint,
     &        Diralb, Difalb, Nintv, Maxflx, Iw1, Iw2, Flxcld )

      END IF

      IF ( Extrap .OR. (.NOT.Nomclr.AND..NOT.Nomcld) ) THEN

C                       ** COMPUTE WEIGHTED AVERAGE OF CLEAR AND CLOUDY
C                       ** RETRIEVALS
C
         DO 100 In = 1, Nintv
            DO 90 Ip = 1, Nflux
               Flxall( Ip, In ) = (1.-Cdfrac) * Flxclr(Ip,In)+
     &                                 Cdfrac * Flxcld(Ip,In) 
90          CONTINUE
100      CONTINUE

      END IF
      
C                       ** PRINT INSTANTANEOUS FLUX AND OPTICAL DEPTH
C
      IF ( Prnt ) THEN

         Newday = Jday .NE. Jdprev
         IF ( Newday ) THEN
            Jdprev = Jday
            WRITE (Unit16, FMT=9020 ) Jday

         END IF

         IF ( Cdfrac.NE.Misval ) Ccover = Cdfrac * 100.

         WRITE (Unit16, FMT=9030 )
         IF ( Clrref.NE.Misval ) THEN
            WRITE (Unit16, FMT=9040 ) Glat, Glon, Time, Smu, Ozone, Pwater,
     &        Clrref, Clrbba, Tauclr, ( (Flxclr(Ip,In),Ip=1,Nflux),
     &        In=1, Nintv )

         ELSE IF ( Clrref == Misval ) THEN
            WRITE (Unit16, FMT=9050 ) Glat, Glon, Time, Smu, Ozone, Pwater,
     &        Clrref, Clrbba, Tauclr, ( (Flxclr(Ip,In),Ip=1,Nflux),
     &        In=1, Nintv )

         END IF

         IF ( Cldref.NE.Misval ) THEN
            WRITE (Unit16, FMT=9060 ) Cldref, Cldbba, Taucld,
     &        ( (Flxcld(Ip,In),Ip=1,Nflux), In=1, Nintv )

         ELSE IF ( Cldref == Misval ) THEN
            WRITE (Unit16, FMT=9070 ) Cldref, Cldbba, Taucld,
     &        ( (Flxcld(Ip,In),Ip=1,Nflux), In=1, Nintv )

         END IF

        
         WRITE (Unit16, FMT=9080 ) ( (Flxall(Ip,In),Ip=1,Nflux), In=1,
     &     Nintv ), Ccover

      END IF

      RETURN



9000  FORMAT ( /, /, '1', 176 ('*'), /, 55X,
     &       'FLUX AND OPTICAL DEPTH ESTMATES', /, 1X, 176 ('*'), / )
9010  FORMAT ( /, T2, 'BBAL--TOA BROADBAND ALBEDO', T46,
     &       'LON---LONGITUDE (DEGREES)', T99,
     &       'SRFDDN-DIFFUSE SURFACE DOWNWARD FLUX (W/SQ M)', /, T2,
     &       'CCOVR-CLOUD COVER (%)', T46,
     &       'NBRE--TOA NARROWBAND REFLECTANCE', T99,
     &       'SRFUP-SURFACE UPWARD FLUX (W/SQ M)', /, T2,
     &       'CSZA--COS( SOLAR ZENITH ANGLE )', T46,
     &       'O3----OZONE AMOUNT (ATM-CM)', T99, 'TAU---OPTICAL DEPTH',
     &       /, T2, 'GMT---NOMINAL GMT (HOURS)', T46,
     &       'SCENE-SCENE TYPE', T99, 'TOADN-TOA DOWNWARD FLUX (W/SQ M)'
     &       , /, T2, 'H2O---WATER VAPOR AMOUNT (CM)', T46,
     &       'SRFDN-SURFACE DOWNWARD FLUX (W/SQ M)', T99,
     &       'TOAUP-TOA UPWARD FLUX (W/SQ M)', /, T2,
     &       'LAT---LATITUDE (DEGREES)' )
9020  FORMAT ( /, /, 1X, '|', 75 ('-'), '>', 3X, 'DAY OF YEAR=', I4, 3X,
     &       '<', 75 ('-'), '|' )
9030  FORMAT ( /, 58X, 2X, '<---------  Shortwave  --------->', 2X,
     &       '<----------  Visible  ---------->', 2X,
     &       '<-------  Near Infrared -------->', /, 1X,
     &       '   LAT     LON UTIME  CSZA    O3   H2O',
     &       '  NBRE  BBAL    TAU', 3 (
     &       '  TOADN  TOAUP  SRFDN  SRFUP SRFDDN'), '   SCENE CCOVR' )
9040  FORMAT ( 1X, F6.2, F8.2, F6.2, 3F6.3, 2F6.3, F7.2, 15F7.1,
     &       '   CLEAR' )
9050  FORMAT ( 1X, F6.2, F8.2, F6.2, 3F6.3, 2F6.0, F7.0, 15F7.0,
     &       '   CLEAR' )
9060  FORMAT ( 39X, 2F6.3, F7.2, 15F7.1, '  CLOUDY' )
9070  FORMAT ( 39X, 2F6.0, F7.0, 15F7.0, '  CLOUDY' )
9080  FORMAT ( 58X, 15F7.1, '   MIXED(Flxall)', F6.1 )
      END

C**********************************************************************C
      SUBROUTINE  Intpol( Unit16,
     &                    Smu, Pwater, Ozone, Ib, Doclr, Docld, Nwlint,
     &                    Ntaua, Ntauc, Mxwavl, Mxntau, Mxnsun, Mxnwvl,
     &                    Mxnh2o, Mxno3, Mxntab, Dirtra, Diftra, Reflea,
     &                    Sphtra, Sphrea, Dirtrc, Diftrc, Reflec,
     &                    Sphtrc, Sphrec )
C**********************************************************************C
C
C     Interpolates table values of transmission and reflectance for
C     actual value of solar zenith angle, water vapor and ozone amount.
C     The results are only functions of the optical depth.
C
C
C   I N P U T   V A R I A B L E S :
C
C   Doclr      : if true, interpolate clear table
C   Docld      : if true, interpolate cloudy table
C   Ib         : table index for selecting 'right' refl/trans
C   Mxnh2o     : max no. of water vapor amounts in common block
C            *tables* in subroutine *reftra*
C   Mxno3      : max no. of ozone amounts in common block
C            *tables* in subroutine *reftra*
C   Mxnsun     : max no. of solar zenith angles in common block
C            *tables* in subroutine *reftra*
C   Mxntab     : max no. of refl/trans tables  in common block
C            *tables* in subroutine *reftra*
C   Mxntau     : max no. of optical depths in common block
C            *tables* in subroutine *reftra*
C   Mxnwvl     : max no. of spectral intervals in common block
C            *tables* in subroutine *reftra*
C   Mxwavl     : max no. of spectral intervals in calling routine
C   Ntaua      : no. of aerosol optical depth values
C   Ntauc      : no. of cloud optical depth values
C   Nwlint     : no. of spectral intervals
C   Ozone      : actual amount of ozone (atm-cm)
C   Pwater     : actual amount of water vapor (cm)
C   Smu        : actual value of cosine of solar zenith angle
C
C   INPUT VARIABLES TRANSFERRED IN COMMON BLOCK -TABLES-:
C
C   Musun,  Wavap,  O3amt,  Scafac, Nmusun, Nwavap,
C   No3amt, Iatm,   Cldftr, Cldrtr, Cldfre, Clspre,
C   Clsptr, Crdftr, Crdrtr, Crdfre, Crspre, Crsptr
C
C   FOR A DESCRIPTION OF THESE VARIABLES SEE -SUBROUTINE Reftra-
C
C
C   O U T P U T   V A R I A B L E S :
C
C   Diftra(L,N)  : L=1 to Nwlint, N=1 to Ntaua, diffuse transmissivity,
C   Dirtra(L,N)  :                direct  transmissivity,
C   Reflea(L,N)  :                diffuse reflectivity,
C   Sphrea(L,N)  :                spherical reflectivity
C   Sphtra(L,N)  :                spherical transmissivity
C            interpolated from table values of clear cases for
C            the actual values of -Smu, Pwater and Ozone-
C
C   Diftrc(L,N)  : L=1 to Nwlint, N=1 to Ntauc,
C   Dirtrc(L,N)  :
C   Reflec(L,N)  :     the same as above, but from cloudy cases of
C   Sphrec(L,M)  :     table
C   Sphtrc(L,N)  :
C
C
C     LOCAL SYMBOLIC DIMENSIONS:
C
C     Maxh2o : max no. of water vapor amounts
C     Maxo3  : max no. of ozone amounts
C     Maxsun : max no. of solar zenith angles
C     Maxtab : max no. of refl/trans tables
C     Maxtau : max no. of optical depths
C     Maxwvl : max no. of spectral intervals
C
C   Called by- Insest
C   Calls- Wrtdm2, Errmsg, Locate, Rint3d
C+---------------------------------------------------------------------+
C
C     .. Parameters ..
      INTEGER Unit16
      INTEGER  Maxtau, Maxsun, Maxwvl, Maxh2o, Maxo3, Maxtab
      PARAMETER  ( Maxtau=15, Maxsun=9, Maxwvl=5, Maxh2o=5, Maxo3=4,
     &           Maxtab=2 )
C     ..
C     .. Scalar Arguments ..
      LOGICAL  Docld, Doclr
      INTEGER  Ib, Mxnh2o, Mxno3, Mxnsun, Mxntab, Mxntau, Mxnwvl,
     &         Mxwavl, Ntaua, Ntauc, Nwlint
      REAL  Ozone, Pwater, Smu
C     ..
C     .. Array Arguments ..
      REAL  Diftra( Mxwavl, * ), Diftrc( Mxwavl, * ),
     &      Dirtra( Mxwavl, * ), Dirtrc( Mxwavl, * ),
     &      Reflea( Mxwavl, * ), Reflec( Mxwavl, * ),
     &      Sphrea( Mxwavl, * ), Sphrec( Mxwavl, * ),
     &      Sphtra( Mxwavl, * ), Sphtrc( Mxwavl, * )
C     ..
C     .. Arrays in Common ..
      INTEGER  Iatm( Maxtab ), Nmusun( Maxtab ), No3amt( Maxtab ),
     &         Nwavap( Maxtab )
      INTEGER  * 2 Cldfre( Maxwvl, Maxsun, Maxtau, Maxh2o, Maxo3,
     &             Maxtab ), Cldftr( Maxwvl, Maxsun, Maxtau, Maxh2o,
     &             Maxo3, Maxtab ), Cldrtr( Maxwvl, Maxsun, Maxtau,
     &             Maxh2o, Maxo3, Maxtab ), Clspre( Maxwvl, Maxsun,
     &             Maxtau, Maxh2o, Maxo3, Maxtab ),
     &             Clsptr( Maxwvl, Maxsun, Maxtau, Maxh2o, Maxo3,
     &             Maxtab ), Crdfre( Maxwvl, Maxsun, Maxtau, Maxh2o,
     &             Maxo3, Maxtab ), Crdftr( Maxwvl, Maxsun, Maxtau,
     &             Maxh2o, Maxo3, Maxtab ), Crdrtr( Maxwvl, Maxsun,
     &             Maxtau, Maxh2o, Maxo3, Maxtab ),
     &             Crspre( Maxwvl, Maxsun, Maxtau, Maxh2o, Maxo3,
     &             Maxtab ), Crsptr( Maxwvl, Maxsun, Maxtau, Maxh2o,
     &             Maxo3, Maxtab )
      REAL  Musun( Maxsun, Maxtab ), O3amt( Maxo3, Maxtab ),
     &      Scafac( Maxtab ), Wavap( Maxh2o, Maxtab )
C     ..
C     .. Local Scalars ..
      LOGICAL  Dimerr, First, Newh2o, Newmu, Newo3, Newtab
      INTEGER  I1, I2, Iprb, J1, J2, L, M1, M2, N
      REAL  Diftr, Dirtr, Fmu, Fo3, Fpw, Prh2o, Prmu, Pro3, Refle,
     &      Sphre, Sphtr
C     ..
C     .. External Functions ..
      LOGICAL  Wrtdm2
      REAL  Rint3d
      EXTERNAL  Wrtdm2, Rint3d
C     ..
C     .. External Subroutines ..
      EXTERNAL  Errmsg, Locate
C     ..
C     .. Common blocks ..
      COMMON / Tables / Musun, Wavap, O3amt, Scafac, Nmusun, Nwavap,
     &       No3amt, Iatm, Cldftr, Cldrtr, Cldfre, Clspre, Clsptr,
     &       Crdftr, Crdrtr, Crdfre, Crspre, Crsptr
C     ..
C     .. Save statement ..
      SAVE  First, Prmu, Prh2o, Pro3, M1, M2, I1, I2, J1, J2, Fmu, Fpw,
     &      Fo3, Iprb
C     ..
C     .. Data statements ..
      DATA  First / .TRUE. / , Prmu / -1.0 / , Prh2o / -1.0 / ,
     &      Pro3 / -1.0 / , Iprb / -1 /
C     ..
C
      IF ( First ) THEN
C
C                       ** CHECK DIMENSIONS IN COMMON BLOCK *TABLES*
         First  = .FALSE.
         Dimerr = .FALSE.
         IF ( Maxtau.NE.Mxntau ) Dimerr = Wrtdm2( Unit16,'MAXTAU'
     &      , Mxntau )
         IF ( Maxsun.NE.Mxnsun ) Dimerr = Wrtdm2( Unit16,'MAXSUN'
     &      , Mxnsun )
         IF ( Maxwvl.NE.Mxnwvl ) Dimerr = Wrtdm2( Unit16,'MAXWVL'
     &      , Mxnwvl )
         IF ( Maxh2o.NE.Mxnh2o ) Dimerr = Wrtdm2( Unit16,'MAXH2O'
     &      , Mxnh2o )
         IF ( Maxo3.NE.Mxno3 ) Dimerr = Wrtdm2( Unit16,'MAXO3 ', Mxno3 )
         IF ( Maxtab.NE.Mxntab ) Dimerr = Wrtdm2( Unit16,'MAXTAB'
     &      , Mxntab )
         IF ( Dimerr ) 
     &      CALL Errmsg(Unit16, 'INTPOL--DIM ERROR(S)',.TRUE.)

      END IF

      Newtab = Ib .NE. Iprb
      IF ( Newtab ) Iprb   = Ib
C
C                       ** LOCATE -SMU-
C
      Newmu  = Smu .NE. Prmu
      IF ( Newmu ) THEN
         CALL Locate( Smu, Musun(1,Ib), Nmusun(Ib), M1, M2, Fmu )
         Prmu   = Smu

      END IF
C
C                       ** LOCATE -PWATER-
C
      Newh2o = Pwater .NE. Prh2o
      IF ( Newh2o ) THEN
         CALL Locate( Pwater, Wavap(1,Ib), Nwavap(Ib), I1, I2, Fpw )
         Prh2o  = Pwater

      END IF
C
C                       ** LOCATE -OZONE-
C
      Newo3  = Ozone .NE. Pro3
      IF ( Newo3 ) THEN
         CALL Locate( Ozone, O3amt(1,Ib), No3amt(Ib), J1, J2, Fo3 )
         Pro3   = Ozone

      END IF

      IF ( .NOT.Newmu .AND. .NOT.Newh2o .AND. .NOT.Newo3 .AND.
     &     .NOT.Newtab ) RETURN


      IF ( Doclr ) THEN
C
C                       ** INTERPOLATE LINEARLY ON CLEAR REFL-TRANS
C                       ** VALUES
C
         DO 20 N = 1, Ntaua
            DO 10 L = 1, Nwlint
               Dirtr  = Rint3d( Fmu, Fpw, Fo3, Crdrtr(L,M1,N,I1,J1,Ib),
     &                  Crdrtr(L,M2,N,I1,J1,Ib),
     &                  Crdrtr(L,M1,N,I2,J1,Ib),
     &                  Crdrtr(L,M2,N,I2,J1,Ib),
     &                  Crdrtr(L,M1,N,I1,J2,Ib),
     &                  Crdrtr(L,M2,N,I1,J2,Ib),
     &                  Crdrtr(L,M1,N,I2,J2,Ib),
     &                  Crdrtr(L,M2,N,I2,J2,Ib) )
               Diftr  = Rint3d( Fmu, Fpw, Fo3, Crdftr(L,M1,N,I1,J1,Ib),
     &                  Crdftr(L,M2,N,I1,J1,Ib),
     &                  Crdftr(L,M1,N,I2,J1,Ib),
     &                  Crdftr(L,M2,N,I2,J1,Ib),
     &                  Crdftr(L,M1,N,I1,J2,Ib),
     &                  Crdftr(L,M2,N,I1,J2,Ib),
     &                  Crdftr(L,M1,N,I2,J2,Ib),
     &                  Crdftr(L,M2,N,I2,J2,Ib) )
               Refle  = Rint3d( Fmu, Fpw, Fo3, Crdfre(L,M1,N,I1,J1,Ib),
     &                  Crdfre(L,M2,N,I1,J1,Ib),
     &                  Crdfre(L,M1,N,I2,J1,Ib),
     &                  Crdfre(L,M2,N,I2,J1,Ib),
     &                  Crdfre(L,M1,N,I1,J2,Ib),
     &                  Crdfre(L,M2,N,I1,J2,Ib),
     &                  Crdfre(L,M1,N,I2,J2,Ib),
     &                  Crdfre(L,M2,N,I2,J2,Ib) )
               Sphre  = Rint3d( Fmu, Fpw, Fo3, Crspre(L,M1,N,I1,J1,Ib),
     &                  Crspre(L,M2,N,I1,J1,Ib),
     &                  Crspre(L,M1,N,I2,J1,Ib),
     &                  Crspre(L,M2,N,I2,J1,Ib),
     &                  Crspre(L,M1,N,I1,J2,Ib),
     &                  Crspre(L,M2,N,I1,J2,Ib),
     &                  Crspre(L,M1,N,I2,J2,Ib),
     &                  Crspre(L,M2,N,I2,J2,Ib) )
               Sphtr  = Rint3d( Fmu, Fpw, Fo3, Crsptr(L,M1,N,I1,J1,Ib),
     &                  Crsptr(L,M2,N,I1,J1,Ib),
     &                  Crsptr(L,M1,N,I2,J1,Ib),
     &                  Crsptr(L,M2,N,I2,J1,Ib),
     &                  Crsptr(L,M1,N,I1,J2,Ib),
     &                  Crsptr(L,M2,N,I1,J2,Ib),
     &                  Crsptr(L,M1,N,I2,J2,Ib),
     &                  Crsptr(L,M2,N,I2,J2,Ib) )
C
C                       ** SCALE BACK REFL/TRANS FUNCTIONS
C
               Dirtra( L, N ) = Dirtr / Scafac( Ib )
               Diftra( L, N ) = Diftr / Scafac( Ib )
               Reflea( L, N ) = Refle / Scafac( Ib )
               Sphrea( L, N ) = Sphre / Scafac( Ib )
               Sphtra( L, N ) = Sphtr / Scafac( Ib )
C
10          CONTINUE
20       CONTINUE

      END IF

      IF ( Docld ) THEN
C
C                       ** INTERPOLATE LINEARLY ON CLOUDY REFL-TRANS
C                       ** VALUES
C
         DO 40 N = 1, Ntauc
            DO 30 L = 1, Nwlint
               Dirtr  = Rint3d( Fmu, Fpw, Fo3, Cldrtr(L,M1,N,I1,J1,Ib),
     &                  Cldrtr(L,M2,N,I1,J1,Ib),
     &                  Cldrtr(L,M1,N,I2,J1,Ib),
     &                  Cldrtr(L,M2,N,I2,J1,Ib),
     &                  Cldrtr(L,M1,N,I1,J2,Ib),
     &                  Cldrtr(L,M2,N,I1,J2,Ib),
     &                  Cldrtr(L,M1,N,I2,J2,Ib),
     &                  Cldrtr(L,M2,N,I2,J2,Ib) )
               Diftr  = Rint3d( Fmu, Fpw, Fo3, Cldftr(L,M1,N,I1,J1,Ib),
     &                  Cldftr(L,M2,N,I1,J1,Ib),
     &                  Cldftr(L,M1,N,I2,J1,Ib),
     &                  Cldftr(L,M2,N,I2,J1,Ib),
     &                  Cldftr(L,M1,N,I1,J2,Ib),
     &                  Cldftr(L,M2,N,I1,J2,Ib),
     &                  Cldftr(L,M1,N,I2,J2,Ib),
     &                  Cldftr(L,M2,N,I2,J2,Ib) )
               Refle  = Rint3d( Fmu, Fpw, Fo3, Cldfre(L,M1,N,I1,J1,Ib),
     &                  Cldfre(L,M2,N,I1,J1,Ib),
     &                  Cldfre(L,M1,N,I2,J1,Ib),
     &                  Cldfre(L,M2,N,I2,J1,Ib),
     &                  Cldfre(L,M1,N,I1,J2,Ib),
     &                  Cldfre(L,M2,N,I1,J2,Ib),
     &                  Cldfre(L,M1,N,I2,J2,Ib),
     &                  Cldfre(L,M2,N,I2,J2,Ib) )
               Sphre  = Rint3d( Fmu, Fpw, Fo3, Clspre(L,M1,N,I1,J1,Ib),
     &                  Clspre(L,M2,N,I1,J1,Ib),
     &                  Clspre(L,M1,N,I2,J1,Ib),
     &                  Clspre(L,M2,N,I2,J1,Ib),
     &                  Clspre(L,M1,N,I1,J2,Ib),
     &                  Clspre(L,M2,N,I1,J2,Ib),
     &                  Clspre(L,M1,N,I2,J2,Ib),
     &                  Clspre(L,M2,N,I2,J2,Ib) )
               Sphtr  = Rint3d( Fmu, Fpw, Fo3, Clsptr(L,M1,N,I1,J1,Ib),
     &                  Clsptr(L,M2,N,I1,J1,Ib),
     &                  Clsptr(L,M1,N,I2,J1,Ib),
     &                  Clsptr(L,M2,N,I2,J1,Ib),
     &                  Clsptr(L,M1,N,I1,J2,Ib),
     &                  Clsptr(L,M2,N,I1,J2,Ib),
     &                  Clsptr(L,M1,N,I2,J2,Ib),
     &                  Clsptr(L,M2,N,I2,J2,Ib) )
C
C                       ** SCALE BACK REFL/TRANS FUNCTIONS
C
               Dirtrc( L, N ) = Dirtr / Scafac( Ib )
               Diftrc( L, N ) = Diftr / Scafac( Ib )
               Reflec( L, N ) = Refle / Scafac( Ib )
               Sphrec( L, N ) = Sphre / Scafac( Ib )
               Sphtrc( L, N ) = Sphtr / Scafac( Ib )

30          CONTINUE
40       CONTINUE

      END IF

      RETURN

      END

C**********************************************************************C
      SUBROUTINE  Locate( X, Xt, Nxt, I1, I2, Fr )
C**********************************************************************C
C
C     LOCATES -X- IN ARRAY -XT-, SO THAT -X- IS INSIDE OF AN INTERVAL OF
C     TWO ADJACENT -XT- VALUES: XT(I1) AND XT(I2). RETURNS I1, I2 AND
C     FRACTION FR (SEE OUTPUT).
C
C     INPUT:    NXT    : NUMBER OF VALUES IN ARRAY -XT-
C               X      : VARIABLE TO BE LOCATED
C               XT(I)  : I=1 TO NXT, TABLE VALUES OF VARIABLE -X-
C
C     OUTPUT:   FR     : ( X - XT(I2) ) / ( XT(I1) - XT(I2) )
C               I1     : INDEX SO THAT XT(I1) .LE.(.GE.) X
C               I2     : INDEX SO THAT XT(I2) .GT.(.LT.) X
C
C      Called by- Locate, Conver, Getalb
C      Calls- none
C+---------------------------------------------------------------------+

C     .. Scalar Arguments ..
      INTEGER  I1, I2, Nxt
      REAL  Fr, X
C     ..
C     .. Array Arguments ..
      REAL  Xt( * )
C     ..
C     .. Local Scalars ..
      INTEGER  I, Imax, Imin
      REAL  Xmax, Xmaxp, Xmin, Xminp, Zero
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC  Abs, Amax1, Amin1
C     ..
C     .. Data statements ..
C
      DATA  Zero / 1E-3 /
C     ..
C
C                       ** FIND MAX AND MIN OF -XT-
      Xminp  = 1.E36
      Xmaxp  = -1.E36
      DO 10 I = 1, Nxt
         Xmax   = Amax1( Xt(I), Xmaxp )
         Xmin   = Amin1( Xt(I), Xminp )
         IF ( Xmax.GT.Xmaxp ) Imax   = I
         IF ( Xmin.LT.Xminp ) Imin   = I
         Xmaxp  = Xmax
         Xminp  = Xmin
10    CONTINUE
C
      I1     = 1
      I2     = Nxt
C                       ** EXTRAPOLATE FOR THE SPECIAL CASES OF:
C                            * A) ONLY ONE -X-
C
      IF ( Nxt == 1 ) THEN
         Fr     = 1.0
         RETURN

      END IF
C                            * B) -X- IS SMALLER THAN ANY OF -XT-S
C
      IF ( X.LT.Xmin ) THEN
         I1     = 1
         I2     = Imin
         Fr     = 0.0
         RETURN
C                            * C) -X- IS GREATER THAN ANY OF -XT-S
C
      ELSE IF ( X.GT.Xmax ) THEN
         I1     = Imax
         I2     = Nxt
         Fr     = 1.0
         RETURN

      END IF
C
C                       ** NORMAL CASE: FIND I1 AND I2, SO THAT -X-
C                       ** IS IN THE INTERVAL OF ( XT(I1), XT(I2) )
C
      Fr     = 1.0
      DO 20 I = 2, Nxt
         IF ( (Xt(I-1).LE.X.AND.X.LE.Xt(I)) .OR.
     &        (Xt(I-1).GE.X.AND.X.GE.Xt(I)) ) THEN
            I1     = I - 1
            I2     = I
            IF ( Abs(Xt(I1)-Xt(I2)).GT.Zero ) Fr     = ( X-Xt(I2) ) /
     &           ( Xt(I1)-Xt(I2) )
            RETURN

         END IF

20    CONTINUE

      END

C**********************************************************************C
      SUBROUTINE  Matalb( Unit16, 
     &                    Jday, Ilat, Ilon, Lat, Nwlint, Solcon, Snowfr,
     &                    Isrmod, Waveln, Diralb, Difalb, Noretr )
C**********************************************************************C
C
C     Reads seasonal snow-free land-albedo data of e. matthews.
C     Scales the spectral reflectances of the reference model to yield a
C     broadband albedo reported in the albedo maps.
C
C
C     REFERENCE:   MATTHEWS, E., 1985: ATLAS OF ARCHIVED VEGETATION,
C                  LAND-USE AND SEASONAL ALBEDO DATA SETS. NASA
C                  TECHNICAL MEMORANDUM 86199, FEBRUARY 1985.
C
C
C     I N P U T :
C                  Unit16     : logical unit number for output
C                  Ilat      : LATITUDE INDEX OF GRIDPOINT IN MATTHES'
C                              MAPS
C                  Ilon      : SAME AS ABOVE, BUT LONGITUDE INDEX
C                  Jday      : DAY OF YEAR
C                  Diralb(I) : I=1 TO NWLINT; SURFACE ALBEDO OF THE
C                              REFERENCE MODEL FOR DIRECT RADIATION
C                  Difalb(I) : I=1 TO NWLINT; SURFACE ALBEDO OF THE
C                              REFERENCE MODEL FOR DIFFUSE RADIATION
C                  Lat       : latitude (deg)
C                  Nwlint    : NO. OF SPECTRAL INTVS IN CALLING ROUTINE
C                  Snowfr    : SNOW COVER FRACTION (range: 0-1)
C                  Solcon(I) : I=1 TO NWLINT; SPECTRAL SOLAR IRRADIANCE
C                              (W/M SQ) AT THE TOP OF THE ATMOSPHERE
C                  Waveln(I) : I=1 TO NWLINT+1,; boundaries of spectral
C                              intervals in ref/tran tables (microns) 
C
C     O U T P U T:
C                  Diralb(I) : I=1 TO NWLINT; SCALED SURFACE ALBEDO FOR
C                              DIRECT RADIATION
C                  Difalb(I) : I=1 TO NWLINT; SCALED SURFACE ALBEDO FOR
C                              DIFFUSE  RADIATION
C                  Noretr    : TRUE, IF ALBEDO RETRIEVAL IS NOT POSSIBLE
C
C     I N T E R N A L   V A R I A B L E S:
C
C     ALBMAT(I,J) : I=1 TO 360, J=1 TO 180; SEASONAL SNOW-FREE ALBEDO
C                   DATA OF MATTHEWS ON 1X1 DEGREE GRID
C     NEWMON      : TRUE IF CURRENT SEASON DIFFERS FROM THE PREVIOUS
C                   ONE, FALSE OTHERWISE
C     PRVMON      : VALUE OF -MONTH- IN PREVIOUS ENTRY
C     MONTH(LU)   : LU=1 TO 4; NAME OF MONTH REPRESENTING A SEASON
C
C     Called by- Insest
C     Calls- Smumon, Refsur
C+---------------------------------------------------------------------+
C
C     .. Scalar Arguments ..
      LOGICAL  Noretr
      INTEGER  Ilat, Ilon, Jday, Nwlint, Unit16
      REAL     lat, Snowfr
      INTEGER Get_Lun
      INTEGER Unit3
C     ..
C     .. Array Arguments ..
      REAL  Difalb( * ), Diralb( * ), Solcon( * ), Waveln( * )
C     ..
C     .. Local Scalars ..
      CHARACTER  Prvmon * 6, Flname * 10
      LOGICAL  Newmon
      INTEGER  I, J, Lu
      INTEGER Imonth, Isrmod, Idummy
      REAL Smuszn, Albref
      REAL  Albedo, Dnwflx, Sclfac, Upwflx
C     ..
C     .. Local Arrays ..
      CHARACTER  Month( 4 ) * 6
      REAL  Albmat( 360, 180 ), Dirasn( 5 ), Difasn( 5 )
C     ..
C     .. External Functions ..
      REAL  Smumon
      EXTERNAL  Smumon
C     ..
C     .. External Subroutines ..
      EXTERNAL  Refsur
C     ..
C     .. Save statement ..
      SAVE  Albmat, Prvmon
C     ..
C     .. Data statements ..
      DATA  Prvmon / '      ' /
      DATA  Month / 'janalb', 'apralb', 'julalb', 'octalb' /
C     ..
C
C                       ** DETERMINE SEASON FOR READING PROPER
C                       ** ALBEDO DATA
C
      IF ( (1.LE.Jday.AND.Jday.LE.59) .OR. Jday.GE.335 ) THEN
         Lu     = 1
         Imonth = 1

      ELSE IF ( 60.LE.Jday .AND. Jday.LE.151 ) THEN
         Lu     = 2
         Imonth = 4

      ELSE IF ( 152.LE.Jday .AND. Jday.LE.243 ) THEN
         Lu     = 3
         Imonth = 7

      ELSE IF ( 244.LE.Jday .AND. Jday.LE.334 ) THEN
         Lu     = 4
         Imonth = 10

      END IF

      Newmon = Month( Lu ) .NE. Prvmon
      IF ( Newmon ) THEN
         Prvmon = Month( Lu )
C
C                    ** READ SEASONAL ALBEDO DATA OF MATTHEWS
C
         Flname = Month( Lu )//'.dat'
         Unit3 = Get_Lun()
         OPEN (Unit=Unit3, File='Auxdata/'//Flname, Form='FORMATTED',
     &        Status='OLD', Err=1000 )
         DO 10 J = 1, 180
            READ ( Unit=Unit3, FMT='(250F6.2,110F6.2)', END=40, ERR=50 )
     &           ( Albmat(I,J), I=1, 360 )
10       CONTINUE
         CLOSE ( Unit=Unit3 )

      END IF

      Albedo = Albmat( Ilon, Ilat ) / 100.0

C                     ** Only snow-free land albedos are provided in
C                     ** the maps of Matthews; return
C
      IF ( Albedo.LE.0.0 .OR. Snowfr.GT.0.0 ) RETURN
      
C                     ** Calculate seasonal average solar zenith cosine
       
      Smuszn = Smumon( Imonth, Lat )
      
C                     ** Get reference spectral surface albedos for
C                     ** seasonal average solar zenith cosine
                   
      CALL Refsur( Smuszn, Isrmod, Snowfr, Nwlint, Waveln, Idummy,
     &             Dirasn, Difasn, Noretr )
      IF ( Noretr ) RETURN
      
C                     ** Estimate broadband reference albedos by
C                     ** averaging spectral data weighted by -Solcon-.
C                     ** (Assume weights for direct and diffuse
C                     ** components as 0.95 and 0.05, respectively.)
C
      Dnwflx = 0.0
      Upwflx = 0.0
      DO 20 I = 1, Nwlint
         Dnwflx = Dnwflx + Solcon( I )
         Upwflx = Upwflx + Solcon( I ) * 
     &            ( 0.95 * Dirasn(I) + 0.05 * Difasn(I) )
20    CONTINUE
      Albref = Upwflx / Dnwflx

      IF ( Albref.GT.0.0 ) THEN
         Sclfac = Albedo / Albref
   
      ELSE
         Sclfac = 1.0
         
      END IF
      
      DO 30 I = 1, Nwlint
         Diralb( I ) = Sclfac * Diralb( I )
         Difalb( I ) = Sclfac * Difalb( I )
30    CONTINUE

      RETURN

40    WRITE (Unit16,FMT='(//,1X,A)' )'MATALB--PREMAT END OF ALBEDO FILE'
50    WRITE (Unit16, FMT='(//,1X,A)')'MATALB--ERROR READING SRF ALBEDOS'
      STOP

1000  WRITE (*, FMT='(/,1X, 2A)' ) 'Error opening file: ', 
     $                                  Trim(Flname)
      STOP

      END

C**********************************************************************C
      REAL FUNCTION  Smumon( Month, Glat )
C**********************************************************************C
C
C        Calls function Avermu to calculate the monthly average solar
C        zenith angles for several latitudes. Daily averages are
C        calculated for each day of the month first, then a simple
C        arithmetic average of the daily values is obtained.
C
C        NOTE: Each year is assumed to have 365 days (no leap year)

C     INPUT VARIABLES:
C     Glat         : latitude (degrees)
C     Month        : month
C
C     INTERNAL VARIABLES:
C     Iday         : day of year (Jan 1 = 1, Dec 31 = 365)
C     Ndays(Month) : Month=1 to 12; number of days in month
C     Solmu(Month) : Month=1 to 12; average solar zenith angle for the
C                    month
C     Fulday       : if true, daily average refers to the time interval
C                    of 24 hours; if false, the time interval is from
C                    sunrise to sunset
C
C     Called by- Matalb
C     Calls- Avermu
C+---------------------------------------------------------------------+
C
C     .. Scalar Arguments ..
      INTEGER  Month
      REAL  Glat
C     ..
C     .. Local Scalars ..
      LOGICAL  Fulday
      INTEGER  Iday, Iday1, Im
C     ..
C     .. Local Arrays ..
      INTEGER  Ndays( 12 )
C     ..
C     .. External Functions ..
      REAL  Avermu
      EXTERNAL  Avermu
C     ..
C     .. Data statements ..
      DATA  Ndays / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
C     ..

      Fulday = .FALSE.

      Smumon = 0.0

      Iday1 = 0
      DO 10 Im = 1, Month-1
         Iday1 = Iday1 + Ndays( Im )
10    CONTINUE
      Iday1 = Iday1 + 1
      
      DO 20 Iday = Iday1, Iday1 - 1 + Ndays( Month )
         Smumon = Smumon + Avermu( Iday, Glat, Fulday )
20    CONTINUE

      Smumon = Smumon / Ndays( Month )


      RETURN
      END

C**********************************************************************C
      REAL  FUNCTION  Avermu( Iday, Glat, Fulday )
C**********************************************************************C
C
C       Calculates the daily-average value of the cosine of the
C       solar zenith angle for a given day and latitude. The
C       average is obtained by integrating the cosine of solar
C       zenith angle from sunrise to sunset, and by the dividing
C       the integral by the legth of the day. Depending on variable
C       Fulday the length of the day is either 24 hours or from sunrise
C       to sunset.
C
C       Note: Since the integration is from sunrise to sunset,
C             average cosines are set to zero for regions experiencing
C             polar night.

C       Input:
C
C         Iday   : day of year. valid range: 1-365. (jan 1 = 1)
C         Glat   : latitude (in degrees). valid range: -90 to +90.
C                  north is positive, south is negative.
C         Fulday : if true, daily average refers to the time interval
C                  of 24 hours; if false, the time interval is from
C                  sunrise to sunset

C       Output:
C
C         Avermu : daily-average of cosine of solar zenith angle

C       Local:
C
C         Daylen : length of day in radians
C         Decl   : declination of sun in degrees
C         Sunris : time of sunrise in hours
C         Sunset : time of sunset in hours
C
C      Called by- Smumon
C      Calls- Sunda2
C+---------------------------------------------------------------------+

C     .. Scalar Arguments ..
      LOGICAL  Fulday
      INTEGER  Iday
      REAL  Glat
C     ..
C     .. Local Scalars ..
      LOGICAL  First
      REAL  Daylen, Decl, Declr, Glatr, Glon, Gmt, Halfd, Halfr, Pi,
     &      Rsun, Solaza, Solmu, Solzen, Sunris, Sunset, Tstime
C     ..
C     .. External Subroutines ..
      EXTERNAL  Sunda2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC  Asin, Cos, Sin
C     ..
C     .. Save statement ..
      SAVE  First
C     ..
C     .. Data statements ..
      DATA  First / .TRUE. /
C     ..

      IF ( First ) THEN
         Pi     = 2.0 * Asin( 1.0 )
         First  = .FALSE.

      END IF

C                  ** Call subroutine Sundat to get time of sunrise and
C                  ** sunset, and declination. The variables Glon and
C                  ** Gmt are used in subroutine Sundat to calculate
C                  ** the true solar time, but they are not needed for
C                  ** calculating the daily average solar zenith angle.
C                  ** They are set to zero here. Subroutine Sundat also
C                  ** returns several other variables which are not used
C                  ** here (Rsun, Solmu, Solzen, Solaza).

      Glon   = 0.0
      Gmt    = 0.0

      CALL Sunda2( Iday, Glat, Glon, Gmt, Rsun, Decl, Sunris, Sunset,
     &             Tstime, Solmu, Solzen, Solaza )
     
C                  ** Calculate length of half day in degrees
C
      Halfd  = ( Sunset-Sunris ) / 2.0 * 360.0 / 24.0

C                  ** Convert degrees to radians
C
      Declr  = Decl * Pi / 180.
      Glatr  = Glat * Pi / 180.
      Halfr  = Halfd * Pi / 180.

C                  ** Calculate the integral of solar zenith angle
C                  ** cosines from sunrise to sunset. This is done
C                  ** analitically.

      Avermu = ( 2.0 * (Halfr*Sin(Declr) *Sin(Glatr)+
     &         Cos(Declr) *Cos(Glatr) *Sin(Halfr)) )

C                  ** Set length of day
C
      IF ( Fulday ) THEN
         Daylen = 2. * Pi

      ELSE

         Daylen = 2. * Halfr

      END IF

C                  ** Calculate daily average value
C
      IF ( Daylen.NE.0.0 ) THEN
         Avermu = Avermu / Daylen

      ELSE

         Avermu = 0.0

      END IF

      RETURN

      END

C**********************************************************************C
      SUBROUTINE  Sunda2( Iday, Glat, Glon, Gmt, Rsun, Decl, Sunris,
     &                    Sunset, Tstime, Solmu, Solzen, Solaza )
C**********************************************************************C
C
C       For a given location, day of the year, time of the day, computes
C       the earth-sun distance factor, the declination, zenith and
C       azimuth angles of the sun, the time of sunrise and sunset, as
C       well as the local time (true solar time).
C
C     >> Note: refraction is not included.
C
C        Reference:
C        Spencer, J. W., 1971: Fourier series representation of the
C                              position of the sun, Search 2 (5), 172
C
C+---------------------------------------------------------------------+
C                            I N P U T
C   Iday    =  day of the year (1 on 1 january)
C   Glat    =  latitude of location in degrees
C   Glon    =  longitude of location in degrees
C   Gmt     =  grennwich meridian time in hours
C
C   >> Note: longitude is negative (positive) west (east) of
C            Greenwich. Latitude is positive (negative) on
C            northern (southern) hemisphere.
C+---------------------------------------------------------------------+
C                           O U T P U T
C   Decl    =  declination of sun in degrees
C   Tstime  =  local time (true solar time) in hours
C   Rsun    =  earth-sun distance factor ( square of the ratio of the
C              mean distance to the actual one)
C   Solaza  =  azimuth angle of sun (in degrees) measured from local
C              north eastward along the horizon (0.le.Solaza.le.+360).
C              (local north is definied by the half-plane of the local
C              meridian containing the north pole between the zenith and
C              nadir
C   Solzen  =  solar zenith angle in degrees
C   Solmu   =  cosine of solar zenith angle
C   Sunris  =  time of sunrise (true solar time) in hours
C   Sunset  =  time of sunset in hours
C+---------------------------------------------------------------------+
C                    L O C A L   V A R I A B L E S
C   Iprday  =  previous day
C   Dangle  =  day angle (in radians)
C   Declr   =  declination of the sun in radians
C   Hdleng  =  half day length
C   Eqtime  =  equation of time in hours
C   Pi      =  3.1415927...
C
C   Called by- Avermu
C   Calls- none
C+---------------------------------------------------------------------+

C     .. Scalar Arguments ..
      INTEGER  Iday
      REAL  Decl, Glat, Glon, Gmt, Rsun, Solaza, Solmu, Solzen, Sunris,
     &      Sunset, Tstime
C     ..
C     .. Local Scalars ..
      LOGICAL  First
      INTEGER  Iprday
      REAL  Cosaza, Dangle, Declr, Dpr, Eqtime, Glatr, Hdleng, Pi, Rpd,
     &      Zero
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC  Abs, Acos, Asin, Cos, Sign, Sin, Tan
C     ..
C     .. Save statement ..
      SAVE  First, Iprday, Declr, Eqtime, Hdleng, Pi, Rpd, Dpr
C     ..
C     .. Data statements ..
      DATA  Iprday / -1 /, Zero / 1.E-5 /, First / .TRUE. /
C     ..
C
      IF ( First ) THEN
         Pi     = 2.0 * Asin( 1.0 )
         Rpd    = Pi / 180.
         Dpr    = 180. / Pi
         First = .FALSE.

      END IF
C
C
C                       ** Compute quantities kept constant throughout
C                       ** the day
C
      IF ( Iday.NE.Iprday ) THEN
         Iprday = Iday
         Dangle = 2. * Pi * ( Iday-1 ) / 365.
         Rsun   = 1.00011 + 0.034221 * Cos( Dangle ) +
     &            0.00128 * Sin( Dangle ) +
     &            0.000719 * Cos( 2. *Dangle ) +
     &            0.000077 * Sin( 2. *Dangle )
         Declr  = ( 0.006918-0.399912 *Cos(Dangle)+
     &            0.070257 *Sin(Dangle)-0.006758 *Cos(2. *Dangle)+
     &            0.000907 *Sin(2. *Dangle)-0.002697 *Cos(3. *Dangle)+
     &            0.00148 *Sin(3. *Dangle) )

         Eqtime = ( 0.000075+0.001868 *Cos(Dangle)-
     &            0.032077 *Sin(Dangle)-0.014615 *Cos(2. *Dangle)-
     &            0.04089 *Sin(2. *Dangle) ) * 12. / Pi

      END IF
            
      Glatr  = Glat * Rpd

C                       ** Calculate length of half-day
C
      IF ( Pi/2.-Abs(Glatr).LE.Abs(Declr)+Zero ) THEN

C                       ** Sun is above horizon all day long
         IF ( Glatr*Declr.GE.0.0 ) Hdleng = 12.0
C                       ** Sun is below horizon all day long
         IF ( Glatr*Declr.LT.0.0 ) Hdleng = 0.0

      ELSE

         Hdleng = Acos( -Tan(Glatr) *Tan(Declr) ) * Dpr / 15.

      END IF

      Sunris = 12. - Hdleng
      Sunset = 12. + Hdleng

C                       ** Handle special case: latitude = +/- 90 deg.
C                       ** azimuth is measured from eastward from
C                       ** date-line.
C
      IF ( Abs(90.-Abs(Glat)).LE.Zero ) THEN
         Tstime = Eqtime + Gmt
         IF ( Tstime.LT.0.0 ) Tstime = Tstime + 24.0
         IF ( Tstime.GT.24.0 ) Tstime = Tstime - 24.0
         Solmu  = Sin( Declr ) * Sin( Glatr )
         IF ( Solmu.GE.0 ) THEN
            Sunris = 0.0
            Sunset = 24.0

         ELSE

            Sunris = 12.0
            Sunset = 12.0

         END IF

         Solaza = Tstime * 15.
         Solzen = Acos( Solmu ) * Dpr
         Decl   = Declr * Dpr
         RETURN

      END IF
C
C                      ** Handle regular case: (latitude .lt. +90 deg
C                      ** and latitude .gt. -90 deg)
C
      Tstime = 24. / 360. * Glon + Eqtime + Gmt
      IF ( Tstime.LT.0.0 ) Tstime = Tstime + 24.0
      IF ( Tstime.GT.24.0 ) Tstime = Tstime - 24.0

      Solmu  = Sin( Glatr ) * Sin( Declr ) +
     &         Cos( Glatr ) * Cos( Declr ) *
     &         Cos( (12.-Tstime) *15. *Rpd )
     
      IF ( Abs(1.-Abs(Solmu)).LE.Zero ) THEN
C
C                       ** Solar zenith angle equals zero or 180 deg.
C
         Solmu  = Sign( 1.0, Solmu )
         Solzen = Acos( Solmu )
         Solaza = Pi / 2.
         
      ELSE

         Solzen = Acos( Solmu )
         Cosaza = ( Sin(Declr)-Sin(Glatr) *Solmu ) /
     &            ( Cos(Glatr) *Sin(Solzen) )
         
         IF ( Cosaza.GT.0.0 ) THEN
            Cosaza = Min( Cosaza, 1.0 )
         ELSE IF ( Cosaza.LT.0.0 ) THEN
            Cosaza = Max( Cosaza, -1.0 )
         ELSE
            Cosaza = 0.0
         END IF
         
         Solaza = Acos( Cosaza )

      END IF

      IF ( Tstime.GT.12.0 ) Solaza = 2. * Pi - Solaza
C
C                       ** Convert radians to degrees
C
      Solaza = Solaza * Dpr
      Solzen = Solzen * Dpr
      Decl   = Declr * Dpr

      RETURN

      END

C**********************************************************************C
      SUBROUTINE  Matcha( Unit16,
     &                    Ntau, Tau, Bbalb, Drttau, Dfttau, Reftau,
     &                    Spttau, Sprtau, Albdir, Albdif, Solcon, Rmu,
     &                    Rsun, Nwlint, Mxwavl, Tauclm, Extrap, Opdept,
     &                    Dirtrn, Diftrn, Reflec, Sphref, Sphtrn,
     &                    Nomatc )
C**********************************************************************C
C
C     1) For each optical depth value in table computes the broadband
C        toa albedo by adding the effect of the surface reflection to
C        the atmospheric reflection;
C     2) Matches the toa measured albedo -Bbalb- with modeled ones:
C        finds the interval of the optical depth array -Tau-  for which
C        the corresponding toa albedos contain the measured -Bbalb-;
C     3) Interpolates optical depth and reflection/transmission from
C        their array values.
C
C   I N P U T :
C
C   Albdif(Il)    : Il=1 to Nwlint, diffuse surface albedo
C   Albdir(Il)    : Il=1 to Nwlint, direct surface albedo
C   Bbalb       : broadband planetary albedo from satellite
C   Dfttau(Il,It) : Il=1 to Nwlint, It=1 to Ntau, diffuse transmissivity
C                 of the standard problem (no surface reflection)
C   Drttau(Il,It) : as above, but direct transmissivity
C   Extrap        : logical flag to indicate whether to extrapolate
C                 refl-trans functions when match of the toa albedos
C                 is not found.
C                 = true, extrapolate refl-trans functions;
C                 = false, do not retrieve refl-trans functions.
C   Mxwavl        : max no. of spectral intervals
C   Ntau          : no. of optical depth values
C   Nwlint        : no. of spectral intervals
C   Reftau(Il,It) : as -Dfttau-, but diffuse reflectivity
C   Rmu           : cosine of solar zenith angle
C   Rsun          : sun-earth distance correction factor
C   Solcon(Il)    : Il=1 to Nwlint; toa solar irradiance at mean sun-
C                 earth distance
C   Sprtau(Il,It) : as -Dfttau-, but spherical reflectivity
C   Spttau(Il,It) : as -Dfttau-, but spherical transmissivity
C   Tau(i)      : optical depth values in the refl-trans table
C   Tauclm      : climatological optical depth
C
C   O U T P U T :
C
C   Diftrn(Il) : Il=1 to Nwlint; interpolated/extrapolated diffuse
C              atmospheric transmittance derived from matching of
C              the toa albedo (does not include surface reflection)
C   Dirtrn(Il) : Il=1 to Nwlint; as above, but direct transmittance
C   Nomatc     : false, if match of albedos is found and retrieval of
C                   refl-trans functions and optical depth are
C                   successful; true otherwise.
C   Opdept     : estimated visible optical depth
C   Reflec(Il) : Il=1 to Nwlint; as -diftrn-, but diffuse reflectivity
C   Sphref(Il) : Il=1 to Nwlint; as -diftrn-, but spherical reflectivity
C   Sphtrn(Il) : Il=1 to Nwlint; as -diftrn-, but spherical transmissiv
C
C
C   I N T E R N A L   V A R I A B L E S:
C
C   Albmax     : max value of -Albtau-
C   Albmin     : min value of -Albtau-
C   Albtau(It) : It=1 to Ntau; toa broadband albedo computed from the
C              refl-trans functions. Includes the effect of the
C              surface reflection.
C   Dimerr     : true if dimension error is detected
C   First      : true on first entry, false thereafter
C   Islope(It) : It=1 to Ntau-1; integer variable signaling the
C              constancy (0), increase (1) or decrease (-1) of
C              -Albtau-, the albedo with tau.
C   I1(it)     : It=1 to Ntau-1; left-boundary of interval containing
C              -Bbalb-
C   I2(it)     : It=1 to Ntau-1; right-boundary of interval containing
C              -Bbalb-
C   Kmax       : sequence number of interval containing -Albmax-
C   Kmin       : sequence number of interval containing -Albmin-
C   Slope      : slope for interpolation of refl-trans functions and
C              optical depth
C   Tmptau(It) : It=1 to Ntau-1; array holding candidate values of
C              -Opdept-
C   Tmpslp(It) : It=1 to Ntau-1; array holding candidate values of
C                -Slope-
C
C     LOCAL SYMBOLIC DIMENSION:
C
C     Maxtau : max no. of optical depth values
C
C   Called by- Insest
C   Calls- Wrtdim, Errmsg
C+---------------------------------------------------------------------+
C
C     .. Parameters ..
      INTEGER  Unit16
      INTEGER  Maxtau
      PARAMETER  ( Maxtau=15 )
C     ..
C     .. Scalar Arguments ..
      LOGICAL  Extrap, Nomatc
      INTEGER  Mxwavl, Ntau, Nwlint
      REAL  Bbalb, Opdept, Rmu, Rsun, Tauclm
C     ..
C     .. Array Arguments ..
      REAL  Albdif( * ), Albdir( * ), Dfttau( Mxwavl, * ), Diftrn( * ),
     &      Dirtrn( * ), Drttau( Mxwavl, * ), Reflec( * ),
     &      Reftau( Mxwavl, * ), Solcon( * ), Sphref( * ), Sphtrn( * ),
     &      Sprtau( Mxwavl, * ), Spttau( Mxwavl, * ), Tau( * )
C     ..
C     .. Local Scalars ..
      LOGICAL  Dimerr, First
      INTEGER  I, Il, In, Ip, It, K1, K2, Kk, Kmax, Kmin, Ncross, Npt
      REAL  Albmax, Albmin, Amaxpr, Aminpr, Dif, Distn, Distp, Dwnflx,
     &      Fc, Replan, Slope, Upwflx, Zero
C     ..
C     .. Local Arrays ..
      INTEGER  I1( Maxtau-1 ), I2( Maxtau-1 ), Islope( Maxtau-1 )
      REAL  Albtau( Maxtau ), Tmpslp( Maxtau-1 ), Tmptau( Maxtau-1 )
C     ..
C     .. External Functions ..
      LOGICAL  Wrtdim
      EXTERNAL  Wrtdim
C     ..
C     .. External Subroutines ..
      EXTERNAL  Errmsg
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC  Abs, Amax1, Amin1
C     ..
C     .. Save statement ..
      SAVE  First
C     ..
C     .. Data statements ..
      DATA  Zero / 1E-4 / , First / .TRUE. /
C     ..
C
C
      IF ( First ) THEN
         First  = .FALSE.
         Dimerr = .FALSE.
         IF ( Maxtau.LT.Ntau ) Dimerr = Wrtdim( Unit16,'MAXTAU', Ntau )
         IF ( Dimerr ) 
     &       CALL Errmsg( Unit16,'MATCHA--DIMENSION ERROR', .TRUE. )

      END IF
C                       ** FOR EACH OPTICAL DEPTH, COMPUTE THE BROADBAND
C                       ** TOA ALBEDO FOR THE 'PLANETARY PROBLEM' (ADD
C                       ** THE EFFECT OF SURFACE REFLECTION)
      Dwnflx = 0.0
      DO 10 Il = 1, Nwlint
         Dwnflx = Dwnflx + Rsun * Rmu * Solcon( Il )
10    CONTINUE
      DO 30 It = 1, Ntau
         Upwflx = 0.0
         DO 20 Il = 1, Nwlint
            Fc     = ( Albdir(Il) *Drttau(Il,It)+
     &               Albdif(Il) *Dfttau(Il,It) ) /
     &               ( 1.-Albdif(Il) *Sprtau(Il,It) )
            Replan = Reftau( Il, It ) + Fc * Spttau( Il, It )
            Upwflx = Upwflx + Replan * Rsun * Rmu * Solcon( Il )
20       CONTINUE
         Albtau( It ) = Upwflx / Dwnflx
30    CONTINUE
C
C                       ** OBTAIN MAX AND MIN VALUES OF -ALBTAU-
C
      Kmin   = 1
      Kmax   = 1
      Albmin = Albtau( 1 )
      Albmax = Albtau( 1 )
      DO 40 It = 2, Ntau
         Aminpr = Albmin
         Amaxpr = Albmax
         Albmin = Amin1( Albtau(It), Albmin )
         Albmax = Amax1( Albtau(It), Albmax )
         IF ( Albmin.NE.Aminpr ) Kmin   = It
         IF ( Albmax.NE.Amaxpr ) Kmax   = It
40    CONTINUE
C
C                       ** CHECK IF -BBALB- IS OUT OF RANGE OF -ALBMIN-
C                       ** -ALBMAX-. IF YES, AND -EXTRAP=FALSE- DO NOT
C                       ** RETRIEVE ANYTHING.
C
      Nomatc = Bbalb .LT. Albmin .OR. Bbalb .GT. Albmax
      IF ( Nomatc .AND. .NOT.Extrap ) RETURN
C
C                       ** NO MATCH FOUND, BUT USE EXTREME VALUES ANYWAY
C
      IF ( Bbalb.LT.Albmin ) THEN
         K1     = Kmin
         K2     = 1
         Slope  = 0.0
         Opdept = Tau( K1 )
         GO TO 120

      ELSE IF ( Bbalb.GT.Albmax ) THEN
         K1     = Kmax
         K2     = 1
         Slope  = 0.0
         Opdept = Tau( K1 )
         GO TO 120

      END IF
C                       ** CHECK IF ONLY ONE -TAU- OR IF CURVE IS
C                       ** HORIZONTAL. IF YES, SINCE IN THIS CASE -TAU-
C                       ** IS NOT WELL DEFINED, USE CLIMATOLOGICAL -TAU-
C                       ** TO GET RETRIEVALS OF REFL/TRANS FUNCTIONS.
C
      IF ( Ntau == 1 .OR. Albmax-Albmin.LE.Zero ) THEN
C
         CALL Errmsg(Unit16,'MATCHA--TAU IS FROM CLIMATOLOGY',.FALSE.)
C
C                       ** FIND INTERVAL SO THAT TAU(I1).LE.TAUCLM.GT.
C                       ** TAU(I2)
         Npt    = 1
         DO 50 It = 2, Ntau
            IF ( Tauclm.LT.Tau(It) ) THEN
               I1( 1 ) = It - 1
               I2( 1 ) = It
               GO TO 90

            END IF

50       CONTINUE

      END IF
C                       ** IN EACH INTERVAL, SEE IF -ALBTAU- IS:
C                       ** 1) CONSTANT  (SLOPE=0);
C                       ** 2) INCREASES (SLOPE.GT.0);
C                       ** 3) DECREASES (SLOPE.LT.0).
C
      DO 60 It = 2, Ntau
         In     = It - 1
         Dif    = Albtau( It ) - Albtau( It-1 )
         IF ( Dif.GT.0.0 ) THEN
            Islope( In ) = 1

         ELSE IF ( Dif.LT.0.0 ) THEN
            Islope( In ) = -1

         END IF

         IF ( Abs(Dif).LE.Zero ) Islope( In ) = 0
60    CONTINUE

C                       ** COUNT THE NUMBER OF POSSIBLE CROSSING POINTS
      Ncross = 1
      DO 70 In = 2, Ntau - 1
         IF ( Islope(In).NE.Islope(In-1) ) Ncross = Ncross + 1
         IF ( Islope(In) == 0 .AND. Islope(In-1) == 
     &        0 ) Ncross = Ncross + 1
70    CONTINUE
C
C                       ** FIND INTERVAL FOR WHICH -ALBTAU(I1).LE.
C                       ** BBALB.LE.ALBTAU(I2)-
      Npt    = 0
      DO 80 It = 2, Ntau
         IF ( (Albtau(It-1).LE.Bbalb.AND.Bbalb.LE.Albtau(It)) .OR.
     &        (Albtau(It-1).GE.Bbalb.AND.Bbalb.GE.Albtau(It)) ) THEN
            Npt    = Npt + 1
            I1( Npt ) = It - 1
            I2( Npt ) = It
C                       ** ALL CROSS-POINTS ARE FOUND, LEAVE LOOP
            IF ( Npt == Ncross ) GO TO 90

         END IF

80    CONTINUE
C                       ** OBTAIN OPTICAL DEPTH ESTIMATES BY LINEARLY
C                       ** INTERPOLATING INTERVAL VALUES
C
90    DO 100 Ip = 1, Npt
         IF ( Abs(Albtau(I1(Ip))-Albtau(I2(Ip))).LE.Zero ) THEN
            Tmpslp( Ip ) = 0.5

         ELSE

            Tmpslp( Ip ) = ( Bbalb-Albtau(I1(Ip)) ) /
     &                     ( Albtau(I2(Ip))-Albtau(I1(Ip)) )

         END IF

         Tmptau( Ip ) = Tmpslp( Ip ) * Tau( I2(Ip) ) +
     &                  ( 1.-Tmpslp(Ip) ) * Tau( I1(Ip) )
100   CONTINUE

      Kk     = 1

      IF ( Npt.GT.0 ) THEN
C
C                       ** MORE THAN ONE CROSS-POINT
C                       ** SELECT -TAU- INTERVAL CLOSEST TO -TAUCLM-
C
         Distp  = Abs( Tauclm-Tmptau(1) )
         DO 110 Ip = 2, Npt
            Distn  = Abs( Tauclm-Tmptau(Ip) )
            IF ( Distn.LT.Distp ) Kk     = Ip
            Distp  = Distn
110      CONTINUE

      END IF

      K1     = I1( Kk )
      K2     = I2( Kk )
      Slope  = Tmpslp( Kk )
      Opdept = Tmptau( Kk )

      !--- begin error checking added by Andy Heidinger
      if ((K1 .lt. 1) .or. (K2 .lt. 1)) then
        RETURN
      endif
      !--- end error checking added by Andy Heidinger

120   DO 130 I = 1, Nwlint
         Dirtrn( I ) = Slope * Drttau( I, K2 ) +
     &                 ( 1.-Slope ) * Drttau( I, K1 )
         Diftrn( I ) = Slope * Dfttau( I, K2 ) +
     &                 ( 1.-Slope ) * Dfttau( I, K1 )
         Reflec( I ) = Slope * Reftau( I, K2 ) +
     &                 ( 1.-Slope ) * Reftau( I, K1 )
         Sphtrn( I ) = Slope * Spttau( I, K2 ) +
     &                 ( 1.-Slope ) * Spttau( I, K1 )
         Sphref( I ) = Slope * Sprtau( I, K2 ) +
     &                 ( 1.-Slope ) * Sprtau( I, K1 )
130   CONTINUE

      RETURN

      END

C**********************************************************************C
      SUBROUTINE  Prtinf( Unit16,
     &                    Smumin, Smumax, Solcon, Climat, Nontob,
     &                    Nobdrc, Satalb, Matthe, Extrap, Prnt, Header, 
     &                    Vernum, Misval, Utau )
C**********************************************************************C
C
C       Prints information on options used and parameters set
C
C     Called by- Sasrab
C     Calls- none

C     .. Scalar Arguments ..
      INTEGER Unit16
      CHARACTER  Header * ( * )
      LOGICAL  Climat, Extrap, Matthe, Nobdrc, Nontob, Prnt, Satalb,
     &         Utau
      REAL  Misval, Smumax, Smumin, Solcon, Vernum
C     ..

      WRITE (Unit16, FMT=9000 ) Vernum
      WRITE (Unit16, FMT='(/,2A)' ) ' (1) HEADER= ', Header
      WRITE (Unit16, FMT=9010 ) Climat, Nontob, Nobdrc, Satalb, Matthe,
     &  Extrap, Prnt, Misval
      WRITE (Unit16, FMT=9050 ) Smumin, Smumax, Solcon

      IF ( Utau ) THEN
         WRITE (Unit16, FMT='(/,A)' )
     &         ' (4) AEROSOL CLIMATOLOGY: User supplied aerosol '//
     &         'optical depth is used'
      ELSE
         WRITE (Unit16, FMT='(/,A)' ) ' (4) AEROSOL CLIMATOLOGY: ' //
     &    'SRA aerosol optical depth is used'
      END IF

      RETURN

9000  FORMAT ( 1X, 125 ('*'), /, 34X,
     &       'SHORTWAVE RADIATION BUDGET SATELLITE ALGORITHM, VERSION ',
     &       F5.2, /, 1X, 125 ('*') )
9010  FORMAT ( /, ' (2) OPTIONS:', /, 5X, 'Climat=', L2, ',  Nontob=',
     &       L2, ',  Nobdrc=', L2, ',  Satalb=', L2, ',  Matthe=', L2,
     &       ',  Extrap=', L2, ',  Prnt=', L2, ',  Misval=', F10.3 )
9050  FORMAT ( /, ' (3) PARAMETERS:', /, 5X, 'Smumin=', F5.2,
     &       ',  Smumax=', F5.2, ',  Solcon=', F6.1 )
      END

C**********************************************************************C
      SUBROUTINE  Refsur( Sunmu, Isrmod, Snowfr, Nwlint, Waveln, Itable,
     &                    Diralb, Difalb, Noretr )
C**********************************************************************C
C
C     Set band surface albedos of reference-models
C
C     REFERENCES:
C              BMRH: Briegleb, B. P., P. Minnis, V. Ramanathan and
C                    E. Harrison, 1986: Comparison of regional clear sky
C                    albedos inferred from satellite observations and
C                    model calculations. J. Climate Appl. Meteor.,
C                    25, 214-226.
C
C              WW:   Wiscombe, W.J., and S. Warren, 1980: A model for
C                    the spectral albedo of snow. I: Pure snow, J.
C                    Atmos. Sci., 37, 2712-2733.
C
C     I N P U T :
C
C     Isrmod    :  surface model selector
C                  =-1, snow/ice model.
C                  = 1, model 1,
C                  = 2, model 2a,
C                  = 3, model 2b,
C                  = 4, model 3,
C                  = 5, model 4,
C                  = 6, model 5,
C                  = 7, model 6,
C                  = 8, model 7,
C                  = 9, model 8,
C                  =10, model 9  of Briegleb et al.
C                  =11, ocean albedo model of Briegleb et al.
C     Itable    :  ref/tran table id
C     Nwlint    :  no. of spectral intervals in calling routine
C     Sunmu     :  cosine of solar zenith angle
C     Snowfr    :  snow/ice cover fraction (range:0-1)
C     Waveln(Il)  :  Il=1 to Nwlint+1; boundaries of spectral
C                    intervals in ref/tran tables (microns)
C
C
C     O U T P U T:
C
C     Difalb(I) : I=1 to Nint, surface albedo for diffuse radiation
C     Diralb(I) : I=1 to Nint, surface albedo for direct radiation
C     Noretr  : true, if surface albedo model does not exist
C
C     I N T E R N A L   V A R I A B L E S:
C
C     Asympa(I)     : I=1 to Nint; asymmetry parameter
C     D1(Isrmod)    : Isrmod=1 to Nmod; empirical constant for the first
C                     component in the models of bmrh.
C     D2(Isrmod)    : same as above, but for the second component
C     Frac(Isrmod)  : Isrmod=1 to Nmod; fraction of the first component
C                     in the two component albedo model of BMRH
C     Onemom(I)     : I=1 to Nint; single scattering coalbedo (1-omega)
C     Snowdf(I)     : I=1 to Nint; diffuse snow albedo
C     Snowdr(I)     : I=1 to Nint; direct snow albedo
C     R01(Isrmod,I) : Isrmod=1 to Nmod, I=1 to Nint; reflectivity of the
C                   first component for 60 degree solar zenith angle
C     R02(Isrmod,I) : same as above, but for the second component
C     Wvl(I)        : wavelengths of boundaries of spectral intervals in
C                     surface albedo models (microns)
C
C     LOCAL SYMBOLIC DIMENSIONS:
C
C     Nmod  : no. of land albedo models
C     Nint  : no. of spectral intervals in albedo models
C
C     Called by- Insest, Matalb
C     Calls- Errmsg, Wrtbad
C 
C+---------------------------------------------------------------------+
C
C     .. Parameters ..
      INTEGER  Nmod, Nint
      PARAMETER  ( Nmod=10, Nint=5 )
C     ..
C     .. Scalar Arguments ..
      LOGICAL  Noretr
      INTEGER  Isrmod, Itable, Nwlint
      REAL  Snowfr, Sunmu
C     ..
C     .. Array Arguments ..
      REAL  Difalb( * ), Diralb( * ), Waveln( * )
C     ..
C     .. Local Scalars ..
      LOGICAL  Inperr, Newtab
      INTEGER  I
      INTEGER Unit16
      REAL  Astar, Asystr, Bstar, Ome, Omestr, P, Prvtab, Rcomp1,
     &      Rcomp2, Refdif, Refdir, Zai
C     ..
C     .. Local Arrays ..
      REAL  Asympa( Nint ), D1( Nmod ), D2( Nmod ), Frac( Nmod ),
     &      Onemom( Nint ), R01( Nmod, Nint ), R02( Nmod, Nint ),
     &      Snowdf( Nint ), Snowdr( Nint ), Wvl( Nint+1 )
C     ..
C     .. External Functions ..
      LOGICAL  Wrtbad
      EXTERNAL  Wrtbad
C     ..
C     .. External Subroutines ..
      EXTERNAL  Errmsg
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC  Alog, Sqrt
C     ..
C     .. Save statement ..
      SAVE  Prvtab
C     ..
C     .. Data statements ..
      DATA  Wvl / 0.2, 0.4, 0.5, 0.6, 0.7, 4.0 /
C               ** 'MIXED FARMING & TALL GRASSLAND        (1)'
      DATA  ( R01(1,I), I=1, Nint ) / 2 * 0.04, 2 * 0.08, 0.24 /
      DATA  ( R02(1,I), I=1, Nint ) / 2 * 0.05, 2 * 0.10, 0.30 /
      DATA  Frac( 1 ) / 0.8 / , D1( 1 ) / 0.4 / , D2( 1 ) / 0.4 /
C               ** 'TALL/MED GRASS & EVERGREEN SHRUBLAND (2A)'
      DATA  ( R01(2,I), I=1, Nint ) / 2 * 0.05, 2 * 0.10, 0.30 /
      DATA  ( R02(2,I), I=1, Nint ) / 2 * 0.08, 2 * 0.17, 0.35 /
      DATA  Frac( 2 ) / 0.8 / , D1( 2 ) / 0.4 / , D2( 2 ) / 0.1 /
C               ** 'SHORT GRASS, MEADOW AND SHRUBLAND    (2B)'
      DATA  ( R01(3,I), I=1, Nint ) / 2 * 0.08, 2 * 0.17, 0.35 /
      DATA  ( R02(3,I), I=1, Nint ) / 2 * 0.05, 2 * 0.10, 0.30 /
      DATA  Frac( 3 ) / 0.8 / , D1( 3 ) / 0.1 / , D2( 3 ) / 0.4 /
C               ** 'EVERGREEN FOREST                      (3)'
      DATA  ( R01(4,I), I=1, Nint ) / 2 * 0.03, 2 * 0.04, 0.20 /
      DATA  ( R02(4,I), I=1, Nint ) / 2 * 0.05, 2 * 0.10, 0.22 /
      DATA  Frac( 4 ) / 0.9 / , D1( 4 ) / 0.1 / , D2( 4 ) / 0.1 /
C               ** 'MIXED DECIDUOUS EVERGREEN FOREST      (4)'
      DATA  ( R01(5,I), I=1, Nint ) / 2 * 0.03, 2 * 0.06, 0.30 /
      DATA  ( R02(5,I), I=1, Nint ) / 2 * 0.03, 2 * 0.04, 0.20 /
      DATA  Frac( 5 ) / 0.5 / , D1( 5 ) / 0.1 / , D2( 5 ) / 0.1 /
C               ** 'DECIDUOUS FOREST                      (5)'
      DATA  ( R01(6,I), I=1, Nint ) / 2 * 0.03, 2 * 0.06, 0.30 /
      DATA  ( R02(6,I), I=1, Nint ) / 2 * 0.07, 2 * 0.13, 0.28 /
      DATA  Frac( 6 ) / 0.9 / , D1( 6 ) / 0.1 / , D2( 6 ) / 0.1 /
C               ** 'TROPICAL EVERGREEN BROADLEAVED FOREST (6)'
      DATA  ( R01(7,I), I=1, Nint ) / 2 * 0.03, 2 * 0.04, 0.20 /
      DATA  ( R02(7,I), I=1, Nint ) / 2 * 0.05, 2 * 0.10, 0.22 /
      DATA  Frac( 7 ) / 0.9 / , D1( 7 ) / 0.1 / , D2( 7 ) / 0.1 /
C               ** 'MEDIUM/TALL GRASSLAND AND WOODLAND    (7)'
      DATA  ( R01(8,I), I=1, Nint ) / 2 * 0.03, 2 * 0.05, 0.25 /
      DATA  ( R02(8,I), I=1, Nint ) / 2 * 0.05, 2 * 0.10, 0.30 /
      DATA  Frac( 8 ) / 0.8 / , D1( 8 ) / 0.1 / , D2( 8 ) / 0.4 /
C               ** 'DESERT                                (8)'
      DATA  ( R01(9,I), I=1, Nint ) / 2 * 0.28, 2 * 0.42, 0.50 /
      DATA  ( R02(9,I), I=1, Nint ) / 2 * 0.15, 2 * 0.25, 0.40 /
      DATA  Frac( 9 ) / 0.5 / , D1( 9 ) / 0.4 / , D2( 9 ) / 0.4 /
C               ** 'TUNDRA                                (9)'
      DATA  ( R01(10,I), I=1, Nint ) / 2 * 0.04, 2 * 0.10, 0.25 /
      DATA  ( R02(10,I), I=1, Nint ) / 2 * 0.07, 2 * 0.13, 0.28 /
      DATA  Frac( 10 ) / 0.5 / , D1( 10 ) / 0.1 / , D2( 10 ) / 0.1 /
C                       ** OPTICAL PARAMETERS FOR SNOW-GRAIN SIZE OF
C                       ** 50 MICRONS (PURE, FRESH SNOW)
      DATA  Onemom / 2.22E-6, 2.6E-6, 6.4E-6, 1.6E-5, 2.7E-3 /
      DATA  Asympa / 0.885, 3 * 0.888, 0.894 /
      DATA  Prvtab / -1 /
C     ..
C
C                       ** SEE IF ALBEDO MODEL EXISTS
C
      IF ( Isrmod.LT.-1 .OR. Isrmod.GT.11 ) THEN
       CALL Errmsg(Unit16,  
     &      'REFSUR--SURFACE MODEL NOT IMPLEMENTED',.FALSE.)
       Noretr = .TRUE.
       RETURN
      END IF
C
C                       ** CHECK SPECTRAL INTERVALS
C
      Newtab = Itable .NE. Prvtab
      IF ( Newtab ) THEN
         Prvtab = Itable
         Inperr = .FALSE.
         IF ( Nwlint.NE.Nint ) Inperr = Wrtbad( Unit16,'NWLINT' )
         DO 10 I = 1, Nint + 1
            IF ( Wvl(I).NE.Waveln(I) ) Inperr = Wrtbad( Unit16
     &           ,'WAVELN' )
10       CONTINUE
         IF ( Inperr ) CALL Errmsg( Unit16, 
     &                     'REFSUR--INCOMPATIBLE SPEC INTVLS',
     &                              .TRUE. )

      END IF
C
      IF ( Isrmod == -1 .OR. Snowfr.GT.0.0 ) THEN
C
C                       ** COMPUTE ALBEDOS OF SEMI-INFINITE SNOW.
C                       ** EQV. 4 AND 7 OF WW.
C
         DO 20 I = 1, Nint
            Ome    = 1. - Onemom( I )
            Omestr = ( 1.-Asympa(I) *Asympa(I) ) * Ome /
     &               ( 1.-Asympa(I) *Asympa(I) *Ome )
            Asystr = Asympa( I ) / ( 1.+Asympa(I) )
            Astar  = 1. - Omestr * Asystr
            Bstar  = Asystr / Astar
            Zai    = Sqrt( 3. *Astar* (1.-Omestr) )
            P      = 2. * Zai / ( 3. *Astar )
C
            Snowdr( I ) = ( Omestr* (1.-Bstar*Zai*Sunmu) ) /
     &                    ( (1.+P) * (1.+Zai*Sunmu) )
            Snowdf( I ) = ( (1.+Bstar)/ (Zai*Zai) * (Zai-Alog(1.+Zai))-
     &                    Bstar/2. ) * ( 2. *Omestr ) / ( 1.+P )
            Diralb( I ) = Snowdr( I )
            Difalb( I ) = Snowdf( I )
20       CONTINUE

      END IF
C
      IF ( Isrmod.GE.1 .AND. Isrmod.LE.10 ) THEN
C
C                       ** LAND ALBEDO MODELS OF BRIEGLEB ET AL.
C                       ** EQV. 7 OF BMRH
C
         DO 30 I = 1, Nint
            Rcomp1 = R01( Isrmod, I ) * ( 1.+D1(Isrmod) ) /
     &               ( 1.+2. *D1(Isrmod) *Sunmu )
            Rcomp2 = R02( Isrmod, I ) * ( 1.+D2(Isrmod) ) /
     &               ( 1.+2. *D2(Isrmod) *Sunmu )
            Diralb( I ) = Frac( Isrmod ) * Rcomp1 +
     &                    ( 1.-Frac(Isrmod) ) * Rcomp2
            Difalb( I ) = Diralb( I )
30       CONTINUE

      ELSE IF ( Isrmod == 11 ) THEN
C
C                       ** OCEAN ALBEDO MODEL OF BRIEGLEB ET AL.
C                       ** EQV. 8 OF BMRH
C
         Refdir = ( 2.6/ (Sunmu**1.7+0.065)+
     &            15. * (Sunmu-0.1) * (Sunmu-0.5) * (Sunmu-1.0) ) / 100.
         Refdif = 0.06
         DO 40 I = 1, Nint
            Diralb( I ) = Refdir
            Difalb( I ) = Refdif
40       CONTINUE

      END IF

      IF ( Snowfr.GT.0.0 ) THEN
C
C                       ** USE SNOW COVERAGE DATA TO COMPUTE A
C                       ** WEIGHTED AVERAGE OF THE ALBEDOS
C
         DO 50 I = 1, Nint
            Diralb( I ) = Snowfr * Snowdr( I ) +
     &                    ( 1.-Snowfr ) * Diralb( I )
            Difalb( I ) = Snowfr * Snowdf( I ) +
     &                    ( 1.-Snowfr ) * Difalb( I )
50       CONTINUE

      END IF

      RETURN

      END

C**********************************************************************C
      SUBROUTINE  Reftra( Unit16,
     &                    Mxtabs, Mxwvls, Mxtaus, Iaer, Ntable, Ntaua,
     &                    Ntauc, Nwlint, Rsuntb, Solirr, Solctb, Tauaer,
     &                    Taucld, Waveln, Mxntau, Mxnsun, Mxnwvl,
     &                    Mxnh2o, Mxno3, Mxntab, Solcon, Vertab,
     &                    Flname )
C**********************************************************************C
C
C     READS REFLECTION-TRANSMISSION TABLES FROM LOGICAL UNIT 8.
C     PRINTS INFORMATION ABOUT ATMOSPHERIC MODELS USED IN TABLES.
C     EACH TABLE CONTAINES THE SPECTRAL VALUES OF REFLECTION,
C     TRANSMISSION, SPHERICAL ALBEDO AND TRANSMISSION IN A SCATTERING-
C     ABSORBING ATMOSPHERE OVER A NON-REFLECTING SURFACE. IN EACH TABLE
C     THERE ARE TWO DATA BLOCKS. THE FIRST BLOCK CONTAINS RESULTS FOR
C     CLOUDY ATMOSPHERES, WHILE THE SECOND ONE IS FOR CLEAR ATMOSPHERES.
C     IN BOTH DATA BLOCKS, RESULTS ARE STORED FOR THE SAME SET OF SOLAR
C     ZENITH ANGLES, WATER VAPOR AND OZONE AMOUNTS. IN ADDITION, FOR
C     EACH COMBINATION OF THE SOLAR ZENITH ANGLE, WATER VAPOR AND OZONE
C     AMOUNT, IN BLOCK 1 THE CLOUD OPTICAL DEPTH, IN BLOCK 2 THE AEROSOL
C     OPTICAL DEPTH IS VARIED.
C
C   I N P U T   V A R I A B L E S:
C
C   Flname :  path and name of file containing ref/tran data
C   MXTABS :  MAX NO. OF TABLES IN CALLING PROGRAM
C   MXTAUS :  MAX NO. OF OPTICAL DEPTHS IN CALLING PROGRAM
C   MXWVLS :  MAX NO. OF SPECTRAL INTERVALS IN CALLING PROGRAM
C
C   O U T P U T   V A R I A B L E S :
C
C   IAER(IB)     : IB=1 TO NTABLE; AEROSOL MODEL ID
C                       1 - MAR-I
C                       2 - MAR-II
C                       3 - CONT-I
C                       4 - CONT-II
C   IATM(IB)     : IB=1 TO NTABLE; ATMOSPHERIC MODEL ID
C                       1 - TROPICAL
C                       2 - MIDLATITUDE SUMMER
C                       3 - MIDLATITUDE WINTER
C                       4 - SUBARCTIC SUMMER
C                       5 - SUBARCTIC WINTER
C   MXNH2O       : MAX NO. OF PRECIP. WATER VALUES IN TABLES
C   MXNO3        : MAX NO. OF OZONE AMOUNT VALUES IN TABLES
C   MXNSUN       : MAX NO. OF SUN ANGLES IN TABLES
C   MXNTAB       : MAX NO. OF REFL/TRANS TABLES
C   MXNTAU       : MAX NO. OF OPTICAL DEPTHS IN TABLES
C   MXNWVL       : MAX NO. OF SPECTRAL INTERVALS IN TABLES
C   MUSUN(IM,IB) : IM=1 TO NMUSUN(IB), IB=1 TO NTABLE; SOLAR ZENITH
C                    ANGLES
C   NMUSUN(IB)   : IB=1 TO NTABLE; NO. OF SOLAR ZENITH ANGLES
C   NO3AMT(IB)   : IB=1 TO NTABLE; NO. OF OZONE AMOUNT VALUES
C   NTABLE       : NO. OF TABLES
C   NTAUA(IB)    : IB=1 TO NTABLE; NO. OF AEROSOL OPTICAL DEPTH VALUES
C   NTAUC(IB)    : IB=1 TO NTABLE; NO. OF CLOUD OPTICAL DEPTH VALUES
C   NTAUMX(IB)   : IB=1 TO NTABLE; MAX( NTAUA(IB), NTAUC(IB) )
C   NWAVAP(IB)   : IB=1 TO NTABLE; NO. OF PRECIP. WATER VALUES
C   NWLINT(IB)   : IB=1 TO NTABLE; NO. OF SPECTRAL INTERVALS
C   O3AMT(IO,IB) : IO=1 TO NO3AMT(IB), IB=1 TO NTABLE; OZONE AMOUNTS
C   RSUNTB(IB)   : IB=1 TO NTABLE; SUN-EARTH DISTANCE FACTOR
C   SCAFAC(IB)   : IB=1 TO NTABLE; SCALE FACTOR . REFL/TRANS VALUES IN
C                  TABLES ARE MULTIPLIED BY -SCAFAC-
C   SOLCTB(IB)   : IB=1 TO NTABLE; SOLAR CONSTANT (W/SQ M)
C   SOLIRR(IL,IB): IL=1 TO NWLINT(IB), IB=1 TO NTABLE; TOA
C                    SPECTRAL SOLAR IRRADIANCE (W/M2) AT SUN-EARTH
C                    DISTANCE RSUNTB(IB)
C   TAUAER(IT,IB): IT=1 TO NTAUA(IB), IB=1 TO NTABLE; AEROSOL OPTICAL
C                    DEPTH VALUES
C   TAUCLD(IT,IB): IT=1 TO NTAUC(IB), IB=1 TO NTABLE; CLOUD OPTICAL
C                    DEPTH VALUES
C   WAVAP(IW,IB) : IW=1 TO NWAVAP(IB), IB=1 TO NTABLE; PRECIPITABLE
C                    WATER VALUES
C   WAVELN(IL,IB): IL=1 TO NWLINT(IB)+1, IB=1 TO NTABLE; BOUNDARIES OF
C                  SPECTRAL INTERVALS
C
C            REFL/TRANS FUNCTIONS FOR   C L O U D Y   CASES
C
C   CLDRTR(IL,IM,IT,IW,IO,IB) : IL=1 TO NWLINT(IB), IT=1 TO NTAUC(IB),
C                               IW=1 TO NWAVAP(IB), IO=1 TO NO3AMT(IB),
C                               IB=1 TO NTABLE; DIRECT TRANSMISSIVITY
C   CLDFTR(L,IM,IT,IW,IO,IB)  : DIFFUSE TRANSMISSIVITY
C   CLDFRE(L,IM,IT,IW,IO,IB)  : BUT DIFFUSE REFLECTIVITY
C   CLSPTR(L,IM,IT,IW,IO,IB)  : BUT SPHERICAL TRANSMISSIVITY
C   CLSPRE(L,IM,IT,IW,IO,IB)  : BUT SPHERICAL REFLECTIVITY
C
C            REFL/TRANS FUNCTIONS FOR   C L E A R   CASES
C
C   CRDRTR(L,IM,IT,IW,IO,IB)  : IL=1 TO NWLINT(IB), IT=1 TO NTAUC(IB),
C                               IW=1 TO NWAVAP(IB), IO=1 TO NO3AMT(IB),
C                               IB=1 TO NTABLE; DIRECT TRANSMISSIVITY
C   CRDFTR(L,IM,IT,IW,IO,IB)  : DIFFUSE TRANSMISSIVITY
C   CRDFRE(L,IM,IT,IW,IO,IB)  : BUT DIFFUSE REFLECTIVITY
C   CRSPTR(L,IM,IT,IW,IO,IB)  : BUT SPHERICAL TRANSMISSIVITY
C   CRSPRE(L,IM,IT,IW,IO,IB)  : BUT SPHERICAL REFLECTIVITY
C
C
C   I N T E R N A L  V A R I A B L E S :
C
C   CLBOT      : CLOUD BOTTOM HEIGHT (KM)
C   CLOUDY     : TRUE, IF RESULTS IN TABLE ARE FOR CLOUDY CASE
C   CLTOP      : CLOUD TOP HEIGHT (KM)
C   CUMAER     : CUMULATIVE AEROSOL OPTICAL DEPTH AT 0.55 MICRONS
C   CUMH2O     : CUMULATIVE WATER VAPOR AMOUNT (G/SQ CM)
C   CUMRAY     : CUMULATIVE RAYLEIGH OPTICAL DEPTH AT 0.55 MICRONS
C   HEADER     : CHARACTER VARIABLE DESCRIBING THE CONTENT OF THE TABLE
C   INPERR     : TRUE, IF DIMENSIONS OF ARRAYS ARE IN ERROR
C   NOO3       : IF TRUE, REFL/TRANS VALUES IN TABLES DOES   N O T
C                         INCLUDE THE EFFECT OF OZONE ABSORPTION
C   NORAY      : IF TRUE, REFL/TRANS VALUES IN TABLES DOES   N O T
C                         INCLUDE MOLECULAR SCATTERING
C   NOSCAT     : IF TRUE, O N L Y  GASEOUS ABSORPTION IS INCLUDED IN
C                         REFL/TRANS VALUES IN TABLES
C   NZ         : NO. OF VERTICAL LEVELS IN MODEL
C   NTAU       : NO. OF OPTICAL DEPTH VALUES IN CURRENT BLOCK
C   O3AMT      : TOTAL OZONE AMOUNT (ATM-CM)
C   P          : PRESSURE AT ATMOSPHERIC LEVELS (MB)
C   T          : TEMPERATURE AT ATMOSPHERIC LEVELS (K)
C   Z          : ALTITUDE OF ATMOSPHERIC LEVELS (KM)
C
C+---------------------------------------------------------------------+
C     LOCAL SYMBOLIC DIMENSIONS:
C
C     MAXH2O  =   MAX NO. OF PRECIP. WATER VALUES
C     MAXO3   =   MAX NO. OF OZONE AMOUNT VALUES
C     MAXSUN  =   MAX NO. OF SUN ANGLES
C     MAXTAB  =   MAX NO. OF REFL/TRANS TABLES
C     MAXTAU  =   MAX NO. OF OPTICAL DEPTHS
C     MAXWVL  =   MAX NO. OF SPECTRAL INTERVALS
C
C   Called by- Sasrab
C   Calls- Wrtdim, Errmsg
C+---------------------------------------------------------------------+
C
C     .. Parameters ..
      INTEGER  Maxtau, Maxsun, Maxwvl, Maxh2o, Maxo3, Maxtab
      PARAMETER  ( Maxtau=15, Maxsun=9, Maxwvl=5, Maxh2o=5, Maxo3=4,
     &           Maxtab=2 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER  Flname * ( * )
      INTEGER  Mxnh2o, Mxno3, Mxnsun, Mxntab, Mxntau, Mxnwvl, Mxtabs,
     &         Mxtaus, Mxwvls, Ntable
      REAL  Solcon, Vertab
C     ..
C     .. Array Arguments ..
      INTEGER  Iaer( * ), Ntaua( * ), Ntauc( * ), Nwlint( * )
      REAL  Rsuntb( * ), Solctb( * ), Solirr( Mxwvls, * ),
     &      Tauaer( Mxtaus, * ), Taucld( Mxtaus, * ),
     &      Waveln( Mxwvls+1, * )
C     ..
C     .. Arrays in Common ..
      INTEGER  Iatm( Maxtab ), Nmusun( Maxtab ), No3amt( Maxtab ),
     &         Nwavap( Maxtab )
      INTEGER  * 2 Cldfre( Maxwvl, Maxsun, Maxtau, Maxh2o, Maxo3,
     &             Maxtab ), Cldftr( Maxwvl, Maxsun, Maxtau, Maxh2o,
     &             Maxo3, Maxtab ), Cldrtr( Maxwvl, Maxsun, Maxtau,
     &             Maxh2o, Maxo3, Maxtab ), Clspre( Maxwvl, Maxsun,
     &             Maxtau, Maxh2o, Maxo3, Maxtab ),
     &             Clsptr( Maxwvl, Maxsun, Maxtau, Maxh2o, Maxo3,
     &             Maxtab ), Crdfre( Maxwvl, Maxsun, Maxtau, Maxh2o,
     &             Maxo3, Maxtab ), Crdftr( Maxwvl, Maxsun, Maxtau,
     &             Maxh2o, Maxo3, Maxtab ), Crdrtr( Maxwvl, Maxsun,
     &             Maxtau, Maxh2o, Maxo3, Maxtab ),
     &             Crspre( Maxwvl, Maxsun, Maxtau, Maxh2o, Maxo3,
     &             Maxtab ), Crsptr( Maxwvl, Maxsun, Maxtau, Maxh2o,
     &             Maxo3, Maxtab )
      REAL  Musun( Maxsun, Maxtab ), O3amt( Maxo3, Maxtab ),
     &      Scafac( Maxtab ), Wavap( Maxh2o, Maxtab )
C     ..
C     .. Local Scalars ..
      CHARACTER  Header * 80, Subhad * 80
      LOGICAL  Cloudy, Dimer1, Dimer2, Noo3, Noray, Noscat
      INTEGER  Ib, Il, Im, Io, It, Iw, K, Ll, Lv, Ncldtb, Nclrtb, Ntau,
     &         Nwl, Nz
      REAL  Clbot, Cltop, Cumaer, Cumh2o, Cumo3, Cumray, P, T, Z
C     ..
C     .. Local Arrays ..
      INTEGER  Ntaumx( Maxtab )
      INTEGER Get_Lun
      INTEGER Unit16
      INTEGER Unit18
C     ..
C     .. External Functions ..
      LOGICAL  Wrtdim
      EXTERNAL  Wrtdim
C     ..
C     .. External Subroutines ..
      EXTERNAL  Errmsg
C     ..
C     .. Common blocks ..
      COMMON / Tables / Musun, Wavap, O3amt, Scafac, Nmusun, Nwavap,
     &       No3amt, Iatm, Cldftr, Cldrtr, Cldfre, Clspre, Clsptr,
     &       Crdftr, Crdrtr, Crdfre, Crspre, Crsptr
C     ..

      Mxntab = Maxtab
      Mxntau = Maxtau
      Mxnsun = Maxsun
      Mxnh2o = Maxh2o
      Mxnwvl = Maxwvl
      Mxno3  = Maxo3

      Unit18 = Get_Lun()
      OPEN ( Unit=Unit18, File=Flname, Form='FORMATTED', Status='OLD',
     &       ERR=1000 )
      
      Dimer1 = .FALSE.
      Dimer2 = .FALSE.

      WRITE (Unit16, FMT='(/,1X,A)' ) '>>> READING REFL-TRANS TABLE(S)'
      WRITE (Unit16, FMT='(/,A)' ) ' (5) REFL-TRANS TABLES:'
      READ (Unit18, FMT= *, END=100, ERR=110 ) Ntable, Vertab
      WRITE (Unit16, FMT=9010 ) Ntable, Vertab

      DO 70 Ib = 1, Ntable

         READ (Unit18, FMT='(A)' ) Header
         READ (Unit18, FMT= *, END=100, ERR=110 ) Noo3, Noray, Noscat
         READ (Unit18, FMT= *, END=100, ERR=110 ) Nwlint( Ib ),
     &     No3amt( Ib ), Nwavap( Ib ), Nmusun( Ib ), Ntaumx( Ib ),
     &     Scafac( Ib )
C
C                            * CHECK DIMENSIONS
C
         IF ( Ntable.GT.Mxtabs ) Dimer1 = Wrtdim( Unit16,'MXTABS', 
     &        Ntable )
         IF ( Nwlint(Ib).GT.Mxwvls ) Dimer1 = Wrtdim( Unit16,'MXWVLS',
     &        Nwlint(Ib) )
         IF ( Ntaumx(Ib).GT.Mxtaus ) Dimer1 = Wrtdim( Unit16,'MXTAUS',
     &        Ntaumx(Ib) )
         IF ( Dimer1 ) CALL Errmsg(Unit16, 'Reftra--DIMENSION ERROR(S)',
     &                              .TRUE. )

         IF ( Ntable.GT.Maxtab ) Dimer2 = Wrtdim( Unit16, 'MAXTAB', 
     &        Ntable )
         IF ( Nwlint(Ib).GT.Maxwvl ) Dimer2 = Wrtdim( Unit16,'MAXWVL',
     &        Nwlint(Ib) )
         IF ( Nwavap(Ib).GT.Maxh2o ) Dimer2 = Wrtdim( Unit16,'MAXH2O',
     &        Nwavap(Ib) )
         IF ( No3amt(Ib).GT.Maxo3 ) Dimer2 = Wrtdim( Unit16,'MAXO3 ',
     &        No3amt(Ib) )
         IF ( Nmusun(Ib).GT.Maxsun ) Dimer2 = Wrtdim( Unit16,'MAXSUN',
     &        Nmusun(Ib) )
         IF ( Ntaumx(Ib).GT.Maxtau ) Dimer2 = Wrtdim( Unit16,'MAXTAU',
     &        Ntaumx(Ib) )
         IF ( Dimer2 ) CALL Errmsg(Unit16,'Reftra--DIMENSION ERROR(S)',
     &                              .TRUE. )

         READ (Unit18, FMT= *, END=100, ERR=110 ) ( Waveln(Il,Ib), Il=1,
     &     Nwlint(Ib)+1 )
         READ (Unit18, FMT= *, END=100, ERR=110 ) ( Solirr(Il,Ib), Il=1,
     &     Nwlint(Ib) ), Rsuntb( Ib ), Solctb( Ib )
         READ (Unit18, FMT= *, END=100, ERR=110 ) ( Wavap(Iw,Ib), Iw=1,
     &     Nwavap(Ib) )
         READ (Unit18, FMT= *, END=100, ERR=110 ) ( O3amt(Io,Ib), Io=1,
     &     No3amt(Ib) )
         READ (Unit18, FMT= *, END=100, ERR=110 ) ( Musun(Im,Ib), Im=1,
     &     Nmusun(Ib) )
         READ (Unit18, FMT= *, END=100, ERR=110 ) Iatm( Ib ),
     &     Iaer( Ib ), Cumo3

         WRITE (Unit16, FMT=9020 ) Ib, Header
         WRITE (Unit16, FMT=9030 ) Nwlint( Ib ), No3amt( Ib ),
     &     Nwavap( Ib ), Nmusun( Ib )
         WRITE (Unit16, FMT=9040 ) ( Waveln(Il-1,Ib), 
     &     Waveln(Il,Ib), Il=2, Nwlint(Ib)+1 )
         WRITE (Unit16, FMT=9050 ) ( Solirr(Il,Ib), Il=1, Nwlint(Ib) )
         WRITE (Unit16, FMT=9060 ) Rsuntb( Ib ), Solctb( Ib )
         WRITE (Unit16, FMT=9070 ) ( O3amt(Io,Ib), Io=1, No3amt(Ib) )
         WRITE (Unit16, FMT=9080 ) ( Wavap(Iw,Ib), Iw=1, Nwavap(Ib) )
         WRITE (Unit16, FMT=9090 ) ( Musun(Im,Ib), Im=1, Nmusun(Ib) )
         WRITE (Unit16, FMT=9120 ) Iatm( Ib ), Iaer( Ib ), Cumo3
         IF ( Noo3 ) WRITE (Unit16,FMT='(A)' )' NO OZONE ABSORPTION'
         IF ( Noray ) WRITE(Unit16,FMT='(A)' )' NO MOLECULAR SCATTERING'
         IF ( Noscat ) WRITE (Unit16, FMT='(A)' ) ' NO SCATTERING'
C
C                            * READ CLEAR AND CLOUDY BLOCKS OF TABLE
C
         Ncldtb = 0
         Nclrtb = 0
         DO 60 K = 1, 2
            READ (Unit18, FMT='(A)', END=100, ERR=110 ) Subhad
            READ (Unit18, FMT= *, END=100, ERR=110 ) Cloudy, Ntau, Nz

            WRITE (Unit16, FMT=9130 ) K, Subhad, Nz, Ntau
            IF ( Cloudy ) THEN
               Ncldtb = Ncldtb + 1
               Ntauc( Ib ) = Ntau
               READ (Unit18, FMT= *, END=100, ERR=110 ) ( Taucld(It,Ib),
     &           It=1, Ntau )
               READ (Unit18, FMT= *, END=100, ERR=110 ) Cltop, Clbot
               WRITE (Unit16, FMT=9100 ) ( Taucld(It,Ib), It=1, Ntau )
               WRITE (Unit16, FMT=9110 ) Cltop, Clbot

            ELSE

               Nclrtb = Nclrtb + 1
               Ntaua( Ib ) = Ntau
               READ (Unit18, FMT= *, END=100, ERR=110 ) ( Tauaer(It,Ib),
     &           It=1, Ntau )
               WRITE (Unit16, FMT=9100 ) ( Tauaer(It,Ib), It=1, Ntau )

            END IF

            WRITE (Unit16, FMT=9140 )
            DO 10 Lv = 1, Nz
               READ (Unit18, FMT= *, END=100, ERR=110 ) Ll, Z, P, T,
     &           Cumh2o, Cumray, Cumaer
               WRITE (Unit16, FMT=9150 ) Ll, Z, P, T, Cumh2o, Cumray,
     &           Cumaer
10          CONTINUE

            Nwl    = Nwlint( Ib )
            DO 50 Io = 1, No3amt( Ib )
               DO 40 Iw = 1, Nwavap( Ib )
                  DO 30 It = 1, Ntau
                     DO 20 Im = 1, Nmusun( Ib )
                        IF ( Cloudy ) THEN
                           READ (Unit18, FMT=9000, END=100,
     &                       ERR=110 ) ( Cldrtr(Il,Im,It,Iw,Io,Ib),
     &                       Il=1, Nwl ), ( Cldftr(Il,Im,It,Iw,Io,Ib),
     &                       Il=1, Nwl ), ( Cldfre(Il,Im,It,Iw,Io,Ib),
     &                       Il=1, Nwl ), ( Clsptr(Il,Im,It,Iw,Io,Ib),
     &                       Il=1, Nwl ), ( Clspre(Il,Im,It,Iw,Io,Ib),
     &                       Il=1, Nwl )

                        ELSE

                           READ (Unit18, FMT=9000, END=100,
     &                       ERR=110 ) ( Crdrtr(Il,Im,It,Iw,Io,Ib),
     &                       Il=1, Nwl ), ( Crdftr(Il,Im,It,Iw,Io,Ib),
     &                       Il=1, Nwl ), ( Crdfre(Il,Im,It,Iw,Io,Ib),
     &                       Il=1, Nwl ), ( Crsptr(Il,Im,It,Iw,Io,Ib),
     &                       Il=1, Nwl ), ( Crspre(Il,Im,It,Iw,Io,Ib),
     &                       Il=1, Nwl )

                        END IF

20                   CONTINUE
30                CONTINUE
40             CONTINUE
50          CONTINUE

60       CONTINUE

C                            * CHECK IF BOTH CLEAR AND CLOUDY TABLES
C                            * ARE PRESENT
         IF ( Nclrtb.NE.1 .AND. Ncldtb.NE.
     &        1 ) CALL Errmsg( Unit16, 
     &        'Reftra--NO CLEAR/CLOUDY TABLE FOUND',
     &        .TRUE. )

70    CONTINUE

      CLOSE ( Unit=Unit18 )
C
C                       ** OBTAIN TOA SPECTRAL IRRADIANCES AT MEAN
C                       ** EARTH-SUN DISTANCE (1 AU)
C
      DO 90 Ib = 1, Ntable
         DO 80 Il = 1, Nwlint( Ib )
            Solirr( Il, Ib ) = Solirr( Il, Ib ) * Solcon /
     &                         Solctb( Ib ) / Rsuntb( Ib )
80       CONTINUE
90    CONTINUE

      RETURN


100   WRITE (Unit16, FMT='(//,1X,A)' ) 'Reftra--PREMATURE END OF FILE'
110   WRITE (Unit16, FMT='(//,1X,A)' )
     &  'Reftra--ERROR READING REF/TRAN TABLE(S)'
      STOP

1000  WRITE (*, FMT='(/,1X, 2A)' ) 'Error opening file: ', 
     &                                   Trim(Flname)
      STOP

9000  FORMAT ( 25I3 )
9010  FORMAT ( /, ' NO. OF REFL-TRANS TABLE(S) TO BE USED :', I4, 5X,
     &       ' VERSION: ', F5.2 )
9020  FORMAT ( /, 1X, 80 ('-'), /, ' REFL-TRANS TABLE NO.', I3, ':', 1X,
     &       A, /, 1X, 80 ('-') )
9030  FORMAT ( /, ' CONTAINES:', I5, ' SPECTRAL INTERVAL(S),', /, 11X,
     &       I5, ' OZONE AMOUNT(S),', /, 11X, I5,
     &       ' WATER VAPOR VALUE(S),', /, 11X, I5,
     &       ' SOLAR ZENITH ANGLE VALUE(S)' )
9040  FORMAT ( /, ' SPEC INTVL(S) (MICR):', 7X,
     &       5 ('(',F5.2,'-',F5.2,')',2X) )
9050  FORMAT ( ' SPEC SOLAR IRRAD (W/SQ M): ', 5 (F11.3,4X) )
9060  FORMAT ( ' SUN-EARTH DIST FACTOR: ', F6.3, ',  SOLAR CONSTANT ',
     &       '(W/SQ M): ', F10.2 )
9070  FORMAT ( ' OZONE AMOUNT(S) (ATM-CM):        ', 5F8.2 )
9080  FORMAT ( ' WATER VAPOR AMOUNT(S) (G/SQ CM): ', 5F8.2 )
9090  FORMAT ( ' COS (SOLAR ZENITH ANGLE):        ', 9F8.2 )
9100  FORMAT ( ' OPTICAL DEPTH(S):', 5 (10F7.2,/,18X) )
9110  FORMAT ( ' CLOUD TOP HEIGHT:', F6.2, ' KM,  CLOUD BOTTOM HEIGHT:',
     &       F6.2, ' KM' )
9120  FORMAT ( /, ' ATMOSPHERIC ID:', I3, ',   AEROSOL ID:', I3,
     &       ',   REF OZONE AMOUNT (ATM-CM):', F7.3 )
9130  FORMAT ( /, /, ' >>> SUB-TABLE', I2, ':', 1X, A, /, ' HAS', I4,
     &       ' LEVELS AND ', I3, ' OPTICAL DEPTH VALUES.', / )
9140  FORMAT ( /, 1X, 15 ('-'), ' REF ATM VERTICAL STRUCTURE ',
     &       15 ('-'), /,
     &     '    HEIGHT    PRESS   TEMP    CUM H2O    CUM RAY    CUM AER'
     &       , /,
     &     '      (KM)     (MB)    (K)  (G/SQ CM)        TAU        TAU'
     &       , / )
9150  FORMAT ( I4, 0P, F6.2, F9.2, F7.1, 1P, 3E11.3 )

      END

C**********************************************************************C
      REAL  FUNCTION  Rint3d( Frx, Fry, Frz, X1y1z1, X2y1z1, X1y2z1,
     &                        X2y2z1, X1y1z2, X2y1z2, X1y2z2, X2y2z2 )
C**********************************************************************C
C
C        INTERPOLATES LINEARLY IN THREE DIMENSIONS
C
C     Called by-  Intpol
C     Calls- none
C
C     .. Scalar Arguments ..
      INTEGER  * 2 X1y1z1, X1y1z2, X1y2z1, X1y2z2, X2y1z1, X2y1z2,
     &             X2y2z1, X2y2z2
      REAL  Frx, Fry, Frz
C     ..
C     .. Local Scalars ..
      REAL  Y1z1, Y1z2, Y2z1, Y2z2, Z1, Z2
C     ..
      Y1z1   = Frx * X1y1z1 + ( 1.-Frx ) * X2y1z1
      Y2z1   = Frx * X1y2z1 + ( 1.-Frx ) * X2y2z1
      Y1z2   = Frx * X1y1z2 + ( 1.-Frx ) * X2y1z2
      Y2z2   = Frx * X1y2z2 + ( 1.-Frx ) * X2y2z2
C
      Z1     = Fry * Y1z1 + ( 1.-Fry ) * Y2z1
      Z2     = Fry * Y1z2 + ( 1.-Fry ) * Y2z2
C
      Rint3d = Frz * Z1 + ( 1.-Frz ) * Z2
C
      RETURN

      END

C**********************************************************************C
      REAL  FUNCTION  Sums( Unit16, Flux, Iw1, Iw2, Nwltot )
C**********************************************************************C

C        Sums up elements of the NwlTot-element array Flux between
C        elements Iw1 and Iw2
C
C     Called by- Fluxes
C     Calls- Errmsg

C     .. Scalar Arguments ..
      INTEGER Unit16
      INTEGER  Iw1, Iw2, Nwltot
C     ..
C     .. Array Arguments ..
      REAL  Flux( * )
C     ..
C     .. Local Scalars ..
      LOGICAL  Inperr
      INTEGER  Iw
C     ..
C     .. External Subroutines ..
      EXTERNAL  Errmsg
C     ..
      Inperr = ( Iw1.LT.1 ) .OR. ( Iw2.GT.Nwltot )
      IF ( Inperr ) CALL Errmsg(Unit16, 'Sums-Input error(s)', .TRUE. )

      Sums = 0.0
      DO 10 Iw = Iw1, Iw2
         Sums = Sums + Flux( Iw )
10    CONTINUE

      RETURN

      END

C**********************************************************************C
      SUBROUTINE  Sundat( Iday, Rsun, Declin )
C**********************************************************************C
C
C     FOR A GIVEN DAY OF THE YEAR COMPUTES THE EARTH-SUN DISTANCE FACTOR
C     AND THE DECLINATION
C
C     REFERENCE:
C        SPENCER, J. W., 1971: FOURIER SERIES REPRESENTATION OF THE
C                              POSITION OF THE SUN. SEARCH 2 (5), 172
C
C   I N P U T   V A R I A B L E:
C
C   IDAY    =  DAY OF THE YEAR (1 ON 1 JANUARY)
C
C   O U T P U T   V A R I A B L E S :
C
C   RSUN    =  EARTH-SUN DISTANCE FACTOR ( SQUARE OF THE RATIO OF THE
C              MEAN DISTANCE TO THE ACTUAL ONE)
C   DECLIN  =  DECLINATION OF SUN (DEGREES)
C
C   L O C A L   V A R I A B L E S:
C
C   DANGLE  =  DAY ANGLE (IN RADIANS)
C   PI      = 3.1418....
C
C   Called by- Sasrab
C   Calls- none
C+---------------------------------------------------------------------+
C     .. Scalar Arguments ..
      INTEGER  Iday
      REAL  Declin, Rsun
C     ..
C     .. Local Scalars ..
      LOGICAL  First
      REAL  Dangle, Declr, Pi
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC  Asin, Cos, Sin
C     ..
C     .. Save statement ..
      SAVE  First, Pi
C     ..
C     .. Data statements ..
      DATA  First / .TRUE. /
C     ..
C
      IF ( First ) Pi     = 2. * Asin( 1.0 )
      First  = .FALSE.
C
      Dangle = 2. * Pi * ( Iday-1 ) / 365.
      Rsun   = 1.00011 + 0.034221 * Cos( Dangle ) +
     &         0.00128 * Sin( Dangle ) + 0.000719 * Cos( 2. *Dangle ) +
     &         0.000077 * Sin( 2. *Dangle )
      Declr  = ( 0.006918-0.399912 *Cos(Dangle)+0.070257 *Sin(Dangle)-
     &         0.006758 *Cos(2. *Dangle)+0.000907 *Sin(2. *Dangle)-
     &         0.002697 *Cos(3. *Dangle)+0.00148 *Sin(3. *Dangle) )
      Declin = Declr * 180. / Pi
C
      RETURN

      END

c **********************************************************************
c ********** Warning and Error-handling routines ***********************
c **********************************************************************

      SUBROUTINE  ERRMSG( Unit16, MESSAG, FATAL )
C
C        PRINT OUT A WARNING OR ERROR MESSAGE;  ABORT IF ERROR
C
      INTEGER UNIT16
      LOGICAL       FATAL, ONCE
      CHARACTER*(*) MESSAG
      INTEGER       MAXMSG, NUMMSG
      SAVE          MAXMSG, NUMMSG, ONCE
      DATA NUMMSG / 0 /,  MAXMSG / 100 /,  ONCE / .FALSE. /
C
C
      IF ( FATAL )  THEN
         WRITE (Unit16, '(2A)' )  ' ******* ERROR >>>>>> ', MESSAG
         STOP
      END IF
C
      NUMMSG = NUMMSG + 1
      IF ( NUMMSG.GT.MAXMSG )  THEN
         IF ( .NOT.ONCE )  WRITE (Unit16,99 )
         ONCE = .TRUE.
      ELSE
         WRITE (Unit16, '(2A)' )  ' ******* WARNING >>>>>> ', MESSAG
      END IF
C
      RETURN
C
   99 FORMAT( ///,' >>>>>>  TOO MANY WARNING MESSAGES --  ',
     $   'THEY WILL NO LONGER BE PRINTED  <<<<<<<', /// )
      END

      LOGICAL FUNCTION  WRTBAD ( Unit16, VARNAM )
C
C          WRITE NAMES OF ERRONEOUS VARIABLES AND RETURN 'TRUE'
C
C      INPUT :   VARNAM = NAME OF ERRONEOUS VARIABLE TO BE WRITTEN
C                         ( CHARACTER, ANY LENGTH )
C ----------------------------------------------------------------------
      CHARACTER*(*)  VARNAM
      INTEGER        MAXMSG, NUMMSG, Unit16
      SAVE  NUMMSG, MAXMSG
      DATA  NUMMSG / 0 /,  MAXMSG / 50 /
C
      WRTBAD = .TRUE.
      NUMMSG = NUMMSG + 1
      WRITE (Unit16, '(3A)' )  ' ****  INPUT VARIABLE  ', VARNAM,
     $                     '  IN ERROR  ****'
      IF ( NUMMSG == MAXMSG )
     $   CALL  ERRMSG (Unit16, 'TOO MANY INPUT ERRORS.  ABORTING...$', 
     $                 .TRUE. )
      RETURN
      END

      LOGICAL FUNCTION  WRTDIM ( Unit16, DIMNAM, MINVAL )
C
C          WRITE NAME OF TOO-SMALL SYMBOLIC DIMENSION AND
C          THE VALUE IT SHOULD BE INCREASED TO;  RETURN 'TRUE'
C
C      INPUT :  DIMNAM = NAME OF SYMBOLIC DIMENSION WHICH IS TOO SMALL
C                        ( CHARACTER, ANY LENGTH )
C               MINVAL = VALUE TO WHICH THAT DIMENSION SHOULD BE
C                        INCREASED (AT LEAST)
C ----------------------------------------------------------------------
      CHARACTER*(*)  DIMNAM
      INTEGER        MINVAL
      INTEGER        Unit16
C
C
      WRITE (Unit16, '(3A,I7)' )  ' ****  SYMBOLIC DIMENSION  ', DIMNAM,
     $                     '  SHOULD BE INCREASED TO AT LEAST ', MINVAL
      WRTDIM = .TRUE.
      RETURN
      END

      LOGICAL FUNCTION  WRTDM2 ( Unit16, DIMNAM, NVAL )
C
C          WRITE NAME OF BAD SYMBOLIC DIMENSION AND
C          THE VALUE IT SHOULD BE EQUAL TO;  RETURN 'TRUE'
C
C      INPUT :  DIMNAM = NAME OF BAD SYMBOLIC DIMENSION
C                        ( CHARACTER, ANY LENGTH )
C               NVAL   = VALUE TO WHICH THAT DIMENSION SHOULD BE EQUAL
C
C ----------------------------------------------------------------------
      CHARACTER*(*)  DIMNAM
      INTEGER        NVAL
      INTEGER        Unit16
C
C
      WRITE (Unit16, '(3A,I7)' )  ' ****  SYMBOLIC DIMENSION  ', DIMNAM,
     $                        '  SHOULD BE EQUAL TO ', NVAL
      WRTDM2 = .TRUE.
      RETURN
      END

c **********************************************************************
c ********** End of Warning and Error-handling routines ****************
c **********************************************************************
