!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: getmod.f90 (src)
!       Getmod (program)
!
! PURPOSE:
!     Reads in vegetation type data of matthews; associates surface type
!     with albedo models of Briegleb et al.; selects appropriate aerosol
!     model from tables, and assigns climatological monthly-mean visible
!     aerosol optical depth for current location.
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
!  CVS info:
! 
!   $Id$
! 
!   $Log: getmod.f90,v $
!   Revision 1.5.2.2  2014/01/26 04:48:33  heidinger
!   updated
!
!   Revision 1.3.2.2  2014/01/24 21:48:00  mhiley
!   Format headers according to CDR General Software Coding Standards
!
!   Revision 1.3.2.1  2013/04/18 16:17:42  dbotambekov
!   change eq to eqv for gfort
! 
!   Revision 1.3  2011/10/20 20:40:59  heidinger
!   print statements
! 
!   Revision 1.2  2011/08/16 18:45:27  heidinger
!   resolved conflicts
! 
!   Revision 1.1.2.3  2011/08/01 03:24:41  heidinger
!   added unit1getmod.f90
! 
!   Revision 1.1.2.2  2011/07/11 19:03:26  heidinger
!   uni16 bugs
! 
!   Revision 1.1.2.1  2011/02/22 23:09:42  heidinger
!   added
! 
!   Revision 1.2  2004/12/20 18:26:35  mikew
!   Added CVS tags.
!
! REFERENCES: 1) Briegleb, B. P., P. Minnis, V. Ramanathan and
!                E. Harrison, 1986: Comparison of regional clear sky
!                albedos inferred from satellite observations and
!                model calculations. J. Climate Appl. Meteor.,
!                25, 214-226.
!
!             2) Matthews, E., 1985: Atlas of archived vegetation,
!                land-use and seasonal albedo data sets. NASA
!                Technical Memorandum 86199, February 1985.
!
!
!   I N P U T:
!
!       Fnaer    : path and name of file containing aerosol optical
!                     depth data
!       Fnvt     : path and name of file containing surface type data
!       Glat     : latitude  of gridpoint (-90 to +90)
!       Glon     : longitude of gridpoint (-180 to +180)
!       Iaer(Ib) : Ib=1 to Ntable, aerosol model id in the refl/trans
!                  table (1=MAR-I; 2=MAR-II; 3=CONT-I; 4=CONT-II,
!                  NOTE: aerosol models 2 and 4 are not implemented in
!                  this version)
!       Ntable   : no. of refl-tran tables
!       Month    : month
!       Utau     : if true, read in user supplied monthly-mean visible
!                     aerosol optical depth data
!
!   O U T P U T :
!
!       Ilat     : latitude index of gridpoint in Matthes' maps
!       Ilon     : longitude index of gridpoint in Matthes' maps
!       Isrmod   : surface albedo model id.
!                         -1  - snow/ice
!                          1  - mixed farming, tall grassland
!                          2  - tall/medium grassland, shrubland
!                          3  - short grassland, meadow and shrubland
!                          4  - evergreen forest
!                          5  - mixed deciduous forest
!                          6  - deciduous forest
!                          7  - tropical evergreen broadleaved forest
!                          8  - medium/tall grassland, woodlands
!                          9  - desert
!                          10 - tundra
!                          11 - water
!       Isrtyp    : surface type; 1=ocean; 2=land; 3=desert, 4=snow/ice
!       Itable    : sequence number of table containing the 'right'
!                   atmosphere and aerosol model
!       Tauclm    : climatological value of aerosol optical depth at
!                   0.55 microns
!
!   I N T E R N A L  V A R I A B L E S :
!
!       First       : true on first entry, false thereafter
!       Frslat      : latitude of center of 1st cell of Matthews' grid
!       Frslon      : longitude of center of 1st cell of Matthews' grid
!       Ialbmo(Iv)  : Iv=1 to 32, array containing the albedo model id
!                     of briegleb et al.
!       Kaer        : aerosol id associated with the surface type
!                     (1=MAR-I; 2=MAR-II; 3=CONT-I; 4=CONT-II)
!       Kaerp       : value of -Kaer- in the previous pass
!       Ksrmp       : value of -Isrmod- in the previous pass
!       Mistau      : value representing missing aerosol optical depth
!       Newaer      : true if current aerosol model differs from the
!                     previous one, false otherwise
!       Newatm      : true if current atmospheric model differs from the
!                     previous one, false otherwise
!       Newmod      : true if current surface model differs from the
!                     previous one, false otherwise
!       Vegtyp(J,I) : J=1 to 180, I=1 to 360, vegetation type of an
!                     1 degree by 1 degree (lon/lat) area
!
!
!--------------------------------------------------------------------------------------
subroutine Getmod(Unit16, Unit1, &
                  Fnvt, Fnaer, Utau, Glat, Glon, Month, Ntable, &
                  Iaer, Itable, Isrmod, Isrtyp, Tauclm, Ilat, Ilon)

   ! Modules
   use CONSTANTS

   ! Parameters
   integer, parameter :: Numlat = 2880
   integer, parameter :: Numlon = 5760

   ! Input Arguments
   integer,           intent(in)    :: Unit16
   integer,           intent(in)    :: Unit1
   character (len=*), intent(in)    :: Fnaer, Fnvt
   logical,           intent(in)    :: Utau
   integer,           intent(in)    :: Month, Ntable, Iaer(*)
   real,              intent(in)    :: Glat
   real,              intent(inout) :: Glon

   ! Output Arguments
   integer,          intent(out) :: Ilat, Ilon, Isrmod, Isrtyp, Itable
   real,             intent(out) :: Tauclm

   ! Local Scalars
   logical, save :: First = .true.
   logical :: Newaer
   integer :: Ib, Iv
   integer :: status = 0
   integer, save :: Kaer, Kaerp, Ksrmp
   real :: Dlat = 0.0625
   real :: Dlon = 0.0625
   real :: Frslat = 89.96875
   real :: Frslon = -179.96875
   real :: Mistau

   ! Local Arrays
   integer :: Ialbmo(0:13) = (/ 11, 4, 7, 6, 8, 5, 8, 3, 3, 3, 2, 1, 9, 10 /)
   integer (kind=int1), save :: Vegtyp(Numlon, Numlat)


   if (First) then

      First = .false.

      !    Initialize aerosol and surface-albedo model ID's
      Kaerp = -9999
      Ksrmp = -9999

      !    Read UMD landcover data
      open(Unit=Unit1, File=trim(Fnvt), Access='DIRECT', Status='OLD', &
          Action='READ', Recl=Numlon, Iostat=status)
      if (status .ne. 0) then
        write (6, FMT='(/,1X, 2A)') 'Error opening file: ', trim(Fnvt)
        stop
      end if
      !    read(Unit=Unit1, Rec=1)
      do Ilat = 1, Numlat
         read(Unit=Unit1, Rec=Ilat, iostat=status) Vegtyp(1:Numlon, Ilat)

         if (status .ne. 0) stop 'GETMOD--ERROR READING VEG TYPE DATA'
      end do
      close(Unit=Unit1)

   end if

   if (Glon .gt. 180.) Glon = Glon - 360.

   ! Get lat/lon indeces for location
   Ilat = int( ((Frslat+Dlat/2.)-Glat)/Dlat ) + 1
   Ilat = min( Ilat, Numlat )
   Ilon = int( (Glon - (Frslon-Dlon/2.))/Dlon ) + 1
   Ilon = min( Ilon, Numlon )

   ! Find vegetation type and albedo of location
   Iv = Vegtyp(Ilon, Ilat)

   ! Associate surface type with one of the land albedo models of 
   ! Briegleb et al.
   Isrmod = Ialbmo(Iv)

   ! Find surface type for current albedo model
   Ksrmp  = Isrmod
   Isrtyp = 2
   if (Isrmod == 11) then
      Isrtyp = 1
   else if (Isrmod == 9) then
       Isrtyp = 3
   else if (Isrmod == -1) then
      Isrtyp = 4
   end if

   ! Get monthly mean aerosol optical depth and aerosol model
   Mistau = -999.0
   Tauclm = Mistau

   if (Utau) call Aertab(Unit16, Fnaer, Glat, Glon, Month, Mistau, Tauclm)

   if (Isrtyp == 1) then
      Kaer = 1
      if (Tauclm == Mistau) Tauclm = 0.05
   else if (Isrtyp == 2) then
      Kaer = 3
      if (Tauclm == Mistau) Tauclm = 0.05
   else if (Isrtyp == 3) then
      Kaer = 3
      if (Tauclm == Mistau) Tauclm = 0.05
   else if (Isrtyp == 4) then
      Kaer = 1
      if (Tauclm == Mistau) Tauclm = 0.05
   end if

   Newaer = Kaer .ne. Kaerp

   ! Use previous aerosol model for this location
   if (.not. Newaer) return

   ! New model for this location. Search thru refl/trans tables until 
   ! a match for -Iaer- is found.
   Kaerp = Kaer

   do Ib = 1, Ntable
      if (Kaer == Iaer(Ib)) then
         Itable = Ib
         return
      end if
   enddo

   call Errmsg('GETMOD--ATMOSPHERIC MODEL NOT FOUND', .true.)

  return

end subroutine Getmod
