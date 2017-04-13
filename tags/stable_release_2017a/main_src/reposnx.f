c--------------------------------------------------------------------------
c Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
c
c NAME: reposnx.f (src)
c
c PURPOSE: 
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

      SUBROUTINE REPOSN (DTS,SLATS,SLONGS,TIMERR,IMAX,JMAX)
c
c     To reposition AVHRR or other data arrays to correct timing error:
c     ---------------------------------------------------------
*     Inputs:
*     dts(imax,jmax) - an array of FOV times in MJDN, possibly erroneous
*     slats(imax,jmax) - an array of FOV geodetic latitudes to be adjusted
*     slongs(imax,jmax) - an array of FOV longitudes to be adjusted
*     timerr - Timing error in seconds; a positive value means the spacecraft
*       clock is too fast.
*     imax,jmax - the cross-scan and along-scan dimensions of the above arrays
c
*     The times, latitudes, and longitudes are adjusted and returned by
c     this routine.
c     -----------------------------------------------------
c     This routine was originally written in Meteorological Fortran
c     (MeteFor), an extension of Fortran-77, in order to utilize the
c     vector and matrix notation available in MeteFor.  This
c     routine is also maintained in MeteFor.  Some of the original
c     MeteFor code may appear in statements which have been commented
c     out.  The original MeteFor source (suffix .hlf) is more
c     readable and self-documenting.  See the document MeteFor.doc.
c     --------------------------------------------------------
      implicit double precision(d)
      parameter(d2=2.d+0)
*
      double precision dtmid(2),dts(imax,jmax)
      real slats(imax,jmax),slongs(imax,jmax)
*     vector (f)vcoord,(f)vec4,(f)vconic,vmid(2)
      real vmid(3,2)
      double precision rq81,rq82
      real vq41(3),vq42(3)
      real  v1(3),v2(3),vam(3),vc(3),vdl(3)
cxr   vcoord,vec4,angbtw,erad84,vconic
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      k3 = imax/3
      k13 = k3 + 1
      k23 = imax - k3
      re = 6371.
*
      do   2 j=1,2
         jj = j
         jg = 1 + (jmax - 1)*(jj - 1)
         dtmid(j) = (dts(k13,jg) + dts(k23,jg))/d2
*
*      v1 = vcoord(dts(k13,jg),vec4(slats(k13,jg),slongs(k13,jg),re),
*    x  'DLC')
      call vec4(slats(k13,jg),slongs(k13,jg),re,vq41)
      call vcoord(dts(k13,jg),vq41,'DLC',v1)
*
*      v2 = vcoord(dts(k23,jg),vec4(slats(k23,jg),slongs(k23,jg),re),
*    x  'DLC')
      call vec4(slats(k23,jg),slongs(k23,jg),re,vq41)
      call vcoord(dts(k23,jg),vq41,'DLC',v2)
*
*     vmid(j) = .5*(v1 + v2)
      vq41(1) = v1(1)+v2(1)
      vq41(2) = v1(2)+v2(2)
      vq41(3) = v1(3)+v2(3)
      vq42(1) = vq41(1)*.5
      vq42(2) = vq41(2)*.5
      vq42(3) = vq41(3)*.5
      vmid(1,j) = vq42(1)
      vmid(2,j) = vq42(2)
      vmid(3,j) = vq42(3)
    2 continue
*
*     angvel = angbtw(vmid(1),vmid(2))/(dtmid(2) - dtmid(1))
      rq81=dtmid(2)
      rq82=dtmid(1)
      rq41=angbtw(vmid(1,1),vmid(1,2))
      angvel=rq41/(rq81-rq82)
*
*     vam = vmid(1) <x> vmid(2)
      vam(1) = vmid(2,1)*vmid(3,2) - vmid(2,2)*vmid(3,1)
      vam(2) = vmid(3,1)*vmid(1,2) - vmid(3,2)*vmid(1,1)
      vam(3) = vmid(1,1)*vmid(2,2) - vmid(1,2)*vmid(2,1)
      terr = timerr/86400.
      ang = angvel*terr
*
      do   4 j=1,jmax
        do   6 i=1,imax
          slat = slats(i,j)
          re = erad84(slat)
*
*         vc = vcoord(dts(i,j), vec4(slat,slongs(i,j),re), 'DLC')
          call vec4(slat,slongs(i,j),re,vq41)
          call vcoord(dts(i,j),vq41,'DLC',vc)
*
          dts(i,j)=dts(i,j)-terr
*
*      vc = vconic(&, vam, ang)
      call vconic(vc,vam,ang,vc)
*      vdl = vcoord(dts(i,j), vc, 'CDL')
      call vcoord(dts(i,j),vc,'CDL',vdl)
       slats(i,j) = vdl(1)
       slongs(i,j) = vdl(2)
    6 continue
    4 continue
*
      return
      end
c     rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
      Subroutine VCOORD(dt,vecin,cmode,  vecout)
*
cxr   sine,cosine,geodet,geocen,arctan,verneq
! -----------------------------------------------------------------
! DOCUMENTATION:
*
! Cooperative Institute for Meteorological Satellite Studies
! University of Wisconsin
! 1225 W Dayton Street
! Madison, Wisconsin  53706      USA
*
! Name: VCOORD; a subroutine-type subprogram
!     This is a vector-valued function under MeteFor.
*
! Inputs: Time, real(double), in Julian Day Number
!         or Modified Julian Day Number (MJDN)
!     a real(single) 3-dimensional vector in some coordinate system
!     a character*3 variable stating the desired transformation
!     12 Jan 00: Modified Julian Day Number (from epoch of
!     1 Jan 70) is also accepted.
*
! Output: a 3-dimensional real(single) transformed vector
*
! Usage:
!     To inter-convert among various coordinate systems:
!      cmode = LLC         Geocentric Lat/Long to Celestial
!      cmode = LLT         Geocentric Lat/Long to Terrestrial
!      cmode = CT          Celestial to Terrestrial
!      cmode = DLC         Geodetic Lat/Long to Celestial
!      cmode = DLT         Geodetic Lat/Long to Terrestrial
*
!     and the inverses of these operations:  (CLL, TLL, TC, DLT, TDL)
*
!     Note regarding time:  An input value for time (dt) is required
!     only if the celestial coordinate system is involved, either as
!     input or output coordinates.  If the transformation is
!     strictly between Lat/Long/CD on the one hand, and terrestrial
!     on the other, a dummy value of time (e.g. zero) can be used.
*
! Modules:  verneq, sine, cosine, arctan
*
! Comment:  The three coordinate systems treated by this routine are
!   defined as follows:
*
!   Celestial: The x-coordinate is directed from the center of
!   the earth toward the vernal equinox.  The y-coordinate is
!   directed 90 deg east of the x-coordinate, and the z-coordinate
!   is directed from the earth's center toward the north pole.
!   The three coordinates form a right-hand orthonormal system.
*
!   Terrestrial: The x-coordinate is directed from the earth's center
!   toward the Greenwich meridian; the y-coordinate is directed
!   90 deg east of the x-coordinate; the z coordinate is directed
!   from the earth's center toward the north pole.  The terrestrial
!   coordinates are a right-hand orthonormal system.
*
!   Latitude/Longitude:  The three components are latitude in
!   degrees, positive north; longitude in degrees, positive east;
!   and radius vector in arbitrary units from the earth's center
!   to a given point.  These coordinates are not orthogonal
!   vectors.  The coordinates are vectors only in the sense that
!   they are three ordered real numbers.  The latitude may be either
!   geocentric or geodetic, depending on the conversion mode.
*
! Programmer: F W Nagle
!  ---------------------------------------------------------------------
      double precision dt,deqr2,deqrad
      character*3 cmodes(10),cmode
      dimension vecin(3),vecout(3),vecinp(3)
      save
      data cmodes/'LLC', 'LLT', 'CT ', 'CLL', 'TLL', 'TC ',
     x  'DLC', 'DLT', 'CDL', 'TDL'/
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      kdetic = 0
      lm = 0
*
      do   2 j=1,10
      lm = j
      if(cmodes(j)(1:3) == cmode(1:3)) go to 8
    2 continue
*
      write(*, '('' VCOORD: Unrecognized Transformation: '',
     1  a3)') CMODE
      return
8     jtest = (lm-1)*(lm-3)*(lm-4)*(lm-6)*(lm-7)*(lm-9)
      velong= 0.
*
      if(jtest == 0)then
      if(dt.gt.2.6d+6.or.dt.lt.0.d+0) go to 920
      if(dt.gt.100000.d+0.and.dt.lt.2.4d+6) go to 920
        velong= verneq(dt)
      end if
*
      do   4 j=1,3
        vecinp(j) = vecin(j)
    4 continue
*
      if(lm == 7.or.lm == 8)then
        lm=lm-6
        vecinp(1)=geocen(vecinp(1))
      end if
*
      if(lm.ge.9)then
*       lm = & - 5
        lm=lm-5
        kdetic = 1
      end if
*
      if(lm.ge.4)lm=3-lm
      if(iabs(lm).gt.1) go to 50
      if(lm.lt.0) go to 20
*
10    if(abs(vecinp(1)).gt.90.)then
      write(*, '('' VCOORD: Bad input latitude: '', f12.2)')
     1  vecinp(1)
      go to 922
      end if
*
      deqrad = vecinp(3) * cosine(vecinp(1))
      clong = vecinp(2) - velong
      vecout(1) = deqrad * cosine(clong)
      vecout(2) = deqrad * sine(clong)
      vecout(3) = vecinp(3) * sine(vecinp(1))
      return
*
20    clong = arctan(vecinp(1), vecinp(2))
      deqr2 = dble(vecinp(1))**2 + dble(vecinp(2))**2
      eqr = dsqrt(deqr2)
      vecout(1) = arctan(eqr, vecinp(3))
      vecout(2) = clong + velong
      if(vecout(2).lt.-180.)vecout(2) = vecout(2)+360.
      if(vecout(2).gt.180.)vecout(2) = vecout(2)-360.
      vecout(3) = dsqrt(deqr2 + dble(vecinp(3))**2)
      if(kdetic.gt.0)vecout(1)=geodet(vecout(1))
      return
*
50    if(iabs(lm).gt.2) go to 100
      if(lm.lt.0) go to 20
      go to 10
*
100   if(lm.lt.0) go to 110
*
102   vecout(3) = vecinp(3)
	clong = arctan(vecinp(1), vecinp(2))
	glong = clong + velong
	deqrad = dsqrt(dble(vecinp(1))**2 + dble(vecinp(2))**2)
	vecout(1) = deqrad * cosine(glong)
	vecout(2) = deqrad * sine(glong)
	return
*
110   velong = -velong
	go to 102
*
920   write(*, '('' VCOORD: The input time is wrong: '', 1pd20.10)')
     1  dt
*
922   vecout(1) = 0.
	vecout(2) = 0.
	vecout(3) = 0.
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
c     Vector*4 Function VEC4(x,y,z)
	Subroutine VEC4(x,y,z,  vec)
c     To return a vector whose components are the three given
c     values x,y,z
*
	real vec(3)
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	vec(1) = x
	vec(2) = y
	vec(3) = z
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Function ANGBTW(V1,V2)
*
*     To compute the Real*4 angle in degrees between two  Vector*4:
*
	Implicit None
	Double Precision dot,rq81,da1,da2,darg,d0,dtest,dsqrt,dabs
	Double Precision d1,dsign,dc,dacos
	Real angbtw
*
	parameter(dc = 180.d+0/3.141592653589793d+0, d0=0.d+0,
     1  d1=1.d+0)
	double precision vv1(3),vv2(3)
	real  v1(3),v2(3)
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	vv1(1) = v1(1)
	vv1(2) = v1(2)
	vv1(3) = v1(3)
*
	vv2(1) = v2(1)
	vv2(2) = v2(2)
	vv2(3) = v2(3)
*
	dot=vv1(1)*vv2(1)+vv1(2)*vv2(2)+vv1(3)*vv2(3)
	da1=vv1(1)*vv1(1)+vv1(2)*vv1(2)+vv1(3)*vv1(3)
	da2=vv2(1)*vv2(1)+vv2(2)*vv2(2)+vv2(3)*vv2(3)
	darg = da1*da2
*
	if(darg.le.d0)then
	write(*, '('' ANGBTW: an argument is zero '')')
	angbtw = 360.
*
	else
	dtest = dot/dsqrt(darg)
	if(dabs(dtest).gt.d1)dtest=dsign(d1,dtest)
	angbtw = dc*dacos(dtest)
	end if
*
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Function ERAD84(gdlat)
*
*     28 Aug 04; added comment
*     9 Jan 03; uses the spheroid WGS 84
*     5 Jan 03; Liam Gumley reports that MODIS uses the WGS
c     spheroid of 1984, in which the Equatorial radius is
c     6378.137 km, and eccentricity .081819190843, which implies
c     a polar radius of 6356.75231424498 km, i.e.
*     b2 = a2(1 - e2)
*
c     To compute the local radius of the earth in kilometers as a
c     function of geodetic latitude:
*     See EARCOMP.hlf
*
	Implicit Double Precision (d)
*
	Parameter (dmajor=6378.137d+0, dminor=6356.752314245d+0,
     1  b = dminor/dmajor, dmaj2=dmajor*dmajor,
     2  dratio=dminor/dmajor, degrad=3.1415926535898d+0/180.d+0)
*     This is the WGS84 Spheroid.
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	if(abs(gdlat) .gt. 90.) then
	write(*, '('' ERAD84: Given Latitude Is Out of Range: '',
     1  1pe16.4)') gdlat
	erad84 = -1.
	end if
*
	gclat = geocen(gdlat)
*     converts to geocentric
	dxe = dmajor * dcos(degrad * gclat)
	dxe2 = dxe*dxe
	dye = dratio*dsqrt(dmaj2 - dxe2)
	erad84 = dsqrt(dxe2 + dye*dye)
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Subroutine VCONIC(vg,vplane,rot,vqonic)
*
cxr   sine,cosine,unit4
*
*   This routine was originally written in Meteorological Fortran
c   (Metefor).  Some of the original code still appears in statements
c   which have been commented out.  The original Metefor routine
c   (*.hlf) is better annotated.
*
*     21 Aug 98; shortened slightly
*     Imagine a cone whose axis is the vector VCONE.  Let the
c     vector VG be an element of the cone directed away from the
c     apex of the cone.  This vector-valued function rotates the
c     vector VG rightward (clockwise) about the axis VCONE, as
c     VCONE is pointed toward the viewer. If the cone is degenerate,
c     i.e. if it is a plane, then this function is the same as
c     VROT4.
*
c       More crudely, imagine looking into the interior of an empty
c     ice cream cone.  The vector VCONE is the axis of the cone,
c     directed toward the viewer.  Let the vector VG
c     be any element of the cone extending from the apex of
c     the cone toward the lip.  If the viewer now rotates the cone
c     'rot' degrees clockwise from his perspective, this function
c     computes the new vector value of VG after the rotation.
*     If rot<0, then the rotation is counter-clockwise.
*     The cone may be merely a plane, in which case VCONE and VG are
c     normal.  A practical application would be determining the direction
c     of view of a conically-scanning satellite.
*
	Implicit None
	real vq41(3),vq42(3),vq43(3),vq44(3),vq45(3)
	real  uv(3),vplane(3),ur(3),vg(3),uh(3),sine,cosine,rot
	real rq41,rq42,rq43,rq44,vqonic(3)
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
*     uv = unit4(vplane)
	call unit4(vplane,uv)
*
*     ur = unit4(vg ** uv)
	vq41(1) = vg(2)*uv(3) - uv(2)*vg(3)
	vq41(2) = vg(3)*uv(1) - uv(3)*vg(1)
	vq41(3) = vg(1)*uv(2) - uv(1)*vg(2)
	call unit4(vq41,ur)
*
*     uh = unit4(uv ** ur)
	vq41(1) = uv(2)*ur(3) - ur(2)*uv(3)
	vq41(2) = uv(3)*ur(1) - ur(3)*uv(1)
	vq41(3) = uv(1)*ur(2) - ur(1)*uv(2)
	call unit4(vq41,uh)
*
*     vconic = (uh*vg) * (ur*sine(rot) + uh*cosine(rot)) +
*    1  (vg*uv) * uv
	rq41=sine(rot)
	rq42=cosine(rot)
	rq43=uh(1)*vg(1)+uh(2)*vg(2)+uh(3)*vg(3)
	vq41(1) = ur(1)*rq41
	vq41(2) = ur(2)*rq41
	vq41(3) = ur(3)*rq41
	vq42(1) = uh(1)*rq42
	vq42(2) = uh(2)*rq42
	vq42(3) = uh(3)*rq42
	vq43(1) = vq41(1)+vq42(1)
	vq43(2) = vq41(2)+vq42(2)
	vq43(3) = vq41(3)+vq42(3)
	rq44 = vg(1)*uv(1)+vg(2)*uv(2)+vg(3)*uv(3)
	vq44(1) = vq43(1)*rq43
	vq44(2) = vq43(2)*rq43
	vq44(3) = vq43(3)*rq43
	vq45(1) = uv(1)*rq44
	vq45(2) = uv(2)*rq44
	vq45(3) = uv(3)*rq44
	vqonic(1) = vq44(1)+vq45(1)
	vqonic(2) = vq44(2)+vq45(2)
	vqonic(3) = vq44(3)+vq45(3)
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Function SINE(xx)
*     Argument in degrees
*
	Implicit Double Precision(d)
	Parameter( dc = 3.141592653589793d+0/180.d+0)
c     xxxxxxxxxxxxxxxxxxxxxxxxxx
	dxx = xx
	ds = dsin(dc*dxx)
	sine = ds
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Function COSINE(xx)
*
	Implicit Double Precision(d)
	Parameter( dc = 3.141592653589793d+0/180.d+0)
c     xxxxxxxxxxxxxxxxxxxxxxxxx
	dxx = xx
	ds = dcos(dc*dxx)
	cosine = ds
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Function GEODET(gcl)
*
*     To convert geocentric latitude to geodetic:
*
*  Given an elliptical earth with semi-major axis = 1, and
*  semi-minor axis = b, then from the equation of an ellipse,
*  we have
*
*  b2 x2 + y2 = b2                                  (1)
*  y = b sqrt(1 - x2)                               (2)
*  sqrt(1 - x2) = y/b                               (3)
*  dy/dx = -bx/sqrt(1 - x2)        from (2)         (4)
*
*  But the slope dy/dx of a point on the ellipse is the negative
*  reciprocal of the tangent of the geodetic latitude, i.e.
*
*  tan D = sqrt(1 - x2)/bx          where D = geodetic latitude
*  bx tan D = sqrt(1 - x2) = y/b    using (3)
*  b2 tan D = y/x
*
*  On the other hand, y/x is the tangent of the geocentric
c  latitude.  Hence
*
c   b2 tan D = tan C,   or    tan D = (tan C)/b2    (5)
*
c   which can be used to interconvert latitudes in either direction.
*
	Parameter (smajx=6378388., sminx=6356911.946,
     1  b=sminx/smajx, b2 = b*b)
c     These values of semi-major and semi-minor axis were taken
c     from  the Encyclopedia Britannica (see Geodesy).
*
*     These values are also shown as the
*     "International Ellipsoid" as defined in
*     "American Practical Navigator" (Bowditch), Pub. No. 9,
*     Volume I, p. 1119, Defense Mapping Agency.
*
	save
cxr   sine,cosine,arctan
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	gc = abs(gcl)
*
	if(gc.ge.90.)then
	 gd = 90.
	  else
	 tanc = sine(gc)/cosine(gc)
	 gd = arctan(b2, tanc)
	end if
*
	geodet = gd*sign(1., gcl)
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Function GEOCEN(gdl)
*
*     See GEODET for explanation of this algorithm.
*
	Implicit None
	Real gdl,gd,abs,gc,tand,sine,cosine,arctan,b2,geocen,sign,
     1  smajx,sminx,b
*
	parameter (smajx=6378388., sminx=6356911.946,
     1  b=sminx/smajx, b2 = b*b)
cxr   sine,arctan,cosine
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	gd = abs(gdl)
*
	if(gd.ge.90.)then
	gc = 90.
	 else
	tand = sine(gd)/cosine(gd)
	gc = arctan(1., b2*tand)
	end if
*
	geocen = gc*sign(1., gdl)
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Function ARCTAN(xx,yy)
*
c     This routine was originally written in Meteorological Fortran
c     (MeteFor).  Some of the original code may appear in statements
c     which have been commented out.  The original MeteFor source
c     code may be more self-documenting and readable.
*
c     19 May 06; on laptop
*     SFS; output is in degrees
*     29 Aug 01; see Approximations for Digital Computers, Hastings.
*
	Implicit None
	Double Precision dpi,d1,dpio4,d0,dr2d,dc1,dc3,dc5,dc7,dc9
	Double Precision dc11,dc13,dc15,dx,dy,dabs,dd,dq,da,da2,dr
	Real xx,yy,arctan
*
	Parameter(dpi=3.14159265358979d+0, d1=1.d+0,
     1  dpio4=dpi/4.d+0, d0=0.d+0, dr2d=180.d+0/dpi,
     1  dc1=.9999993329d+0, dc3=-.3332985605d+0,
     2  dc5=.1994653599d+0, dc7=-.1390853351d+0, dc9=.964200441d-1,
     3  dc11=-.559098861d-1, dc13=.218612288d-1, dc15=-.4054058d-2)
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	dx = xx
	dy = yy
*
	if(dabs(dx).le.d0) then
	if(dabs(dy).le.d0) then
	write(*, '('' ARCTAN: both arguments are zero.'')')
	dd = 720.d+0
	end if
*
	dd = 90.d+0
	if(dy.lt.d0)dd = -90.d+0
*
	else
	dq  = dabs(dy/dx)
	da = (dq - d1)/(dq + d1)
	da2 = da*da
	dr = dpio4 + da*(dc1 + da2*(dc3 + da2*(dc5 + da2*(dc7 + da2*(dc9
     1  + da2*(dc11 + da2*(dc13 + da2*dc15)))))))
	dd = dr2d*dr
*
	if(dx.gt.d0)then
	if(dy.lt.d0)dd = -dd
	 else
	dd = 180.d+0 - dd
	if(dy.lt.d0) dd = -dd
	end if
	end if
*
	arctan = dd
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Function VERNEQ(dt)
*
cxr   dfrac
*     To return the longitude (positive east) of the Vernal Equinox
c     at the given time dt, expressed either as JDN or MJDN:
*     See page 37  Montenbruck & Pfleger
*
*     Compared with old VERNEQ.for  9 Dec 91;  constant error of
c     .003 deg;  looks OK.
*
	Implicit None
	Double Precision dt,mjd,mjd0,dint,ut,t,gst,gstc,dfrac,ve
	Real verneq
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	mjd = dt - 2400000.5d0
	if(dt .lt. 1721060.d+0) mjd = mjd + 2440588
	mjd0 = dint(mjd)
	ut = 24.d0 * (mjd - mjd0)
	t = (mjd0 - 51544.5d0)/36525.d0
	gst = 6.697374558d0 + 1.0027379093d0 * ut
     1 + (8640184.812866d0+(.093104d0 - 6.2d-6*t)*t)*t/3600.
	gstc = dfrac(gst/24.d0)
	ve = 360. * (1.d0 - gstc)
	if(ve .gt. 180.d0) ve = ve - 360.
	verneq = ve
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
c     Vector Function UNIT4(VEC)
	Subroutine Unit4(vec,uvec)
*
*     2 Apr 03; to compute the norm of an arbitrary vector:
c
	real vec(3),uvec(3)
	equivalence(jr2,r2),(sum2,jsum2)
*     data j$01/z'20000000'/             dec 536870912
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	sum2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
*
	if(sum2 .gt. 0.)then
	jr2 = jsum2/2 + 536870912
*
	do 2 n = 1,6
	r1 = r2
2     r2 = .5*(r1 + sum2/r1)
*
	do 4 j = 1,3
4     uvec(j) = vec(j)/r2
*
	    else
	do 6 j = 1,3
6     uvec(j) = 0.
*
	write(*, '('' UNIT4: Warning - the given vector argument has '',
     1  ''zero length.'')')
	end if
*
	return
	end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Double Precision Function DFRAC(ditem)
*
*     The function returns the fraction from the algebraically
*     lesser integer, i.e.
*
c     1.7 ===> .7   and   -1.7 ===> .3
*
	Implicit Double Precision (a-h,o-z)
cxr   dflt
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	df = ditem
*
	if(dabs(df) .gt. 2147483647.9999d0) then
	write(*, '('' DFRAC: argument is too large: '', d20.10)') df
	dfrac = 999999999.
*
	  else
	item = df
	df = df - dflt(item)
	if(df .lt. 0.d+0) df = df + 1.d+0
	dfrac = df
	end if
*
      return
      end
*      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
      Double Precision Function DFLT(narg)
c     DFLOAT does not work on some Unix platforms
      Double Precision dd
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      dd = narg
      dflt = dd
      return
      end
c     rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
      Double Precision Function DTMJDN(myrx,mon,mday,jhour,jmin,jsec)
c
c     To compute the Modified Julian Day Number (MJDN), based on
c     the epoch of 12z 1 January 1970:
c
*     The Modified Julian Day Number is the Julian Day Number
c     of a given date, with the JDN of 12z 1 Jan 1970 subtracted.
c
*      The last 3 arguments may be either REAL*4 or INTEGER*4.
c
*     This routine should not be used for dates earlier than
c     1 January 1970.
c
      Implicit None
      Double Precision djul,dhour,dmin,dsec
      Real hour,fmin,sec
      Integer khour,kmin,ksec,imax,myr,myrx,jhour,jmin,jsec,kalday
      Integer mday,jday,month,mon,jcode,m100,mod,m400,m4,m,ndm
      Integer julday
*
      equivalence(khour,hour),(kmin,fmin),(ksec,sec)
      save
      data imax/16777216/
cxr   julday
c     xxxxxxxxxxxxxxxxxxxxxxxxxxx
      myr = myrx
      khour = jhour
      kmin = jmin
      ksec = jsec
      if(khour .lt. imax)hour = khour
      if(kmin .lt. imax)fmin = kmin
      if(ksec .lt. imax)sec = ksec
      kalday = mday
      jday = kalday
      month = mon
      if(mon .ne. 0) go to 10
*
      jcode = 1115212
      m100 = mod(myr,100)
      m400 = mod(myr,400)
      m4 = mod(myr,4)
      if(m400 == 0 .or. (m4 == 0 .and. m100.ne.0))jcode = 1115208
*
      do 12 m = 1,12
      month = m
      ndm = 31 - mod(jcode,4)
      jday = jday - ndm
      if(jday.le.0) go to 20
*
12    jcode = jcode/4
*
20    kalday = jday + ndm
*
      if(month.gt.12.or.kalday.gt.31)then
      write(*, '('' DTMJDN: month & day:'', 2i8)') month,kalday
      stop 99
      end if
*
10    djul = julday(myr, month, kalday) - 2440588
      dhour = hour
      dmin = fmin
      dsec = sec
      djul = djul + (dhour + dmin/60.d+0 + dsec/3600.d+0)/24.d+0
      dtmjdn = djul - .5d+0
      if(djul .le. 0.d+0)dtmjdn = 0.d+0
      return
      end
c     rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
      Integer Function JULDAY(jyear,month,jday)
c
c     1 Jun 07; logic somewhat modernized
*     27 Jun 98; 2-digit years less than 70 are construed as
c     21st century.   jyear may optionally be 4 digits to avoid
c     ambiguity.
*
c     To get the Julian day number of the day beginning
c     at 12 Z of the given calendar date.  If jyear is given with
c     only two digits, the 20th century is assumed if it is greater
c     than 69, and the 21st century otherwise.  The routine is
c     valid for any date of the Gregorian calendar later than
c     1 Jan 100 AD, and earlier than 1 Jan 4000 AD, since
c     years divisible by 4000 are common years in the Gregorian
c     calendar. If month = 0, then jday is interpreted as the numeric
c     day of the year, 1-366.
cxr   julday
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      if(jyear.lt.0 .or. month.lt.0 .or. month.gt.12) go to 160
      if(jday.le.0 .or. jday.gt.366) go to 160
      if(month.ne.0 .and. jday.gt.31) go to 160
c
      jyr = jyear
      if(jyr .lt. 70) jyr = jyr + 2000
      if(jyr .lt. 100) jyr = jyr + 1900
      kxtra = 1
      kwads = jyr/400
      jyrup = 400*kwads
      nquads = (jyr - jyrup)/4
      julian = 1721060 + 146097*kwads + 1461*nquads - nquads/25
c      1721060 is the J.D.N. of 1 Jan 0000 AD.
      if(mod(nquads,25) /= 0 .or. nquads == 0) go to 110
c     in a 'short' quadrennium
*
      julian = julian + 1
      kxtra = 0
c     'JULIAN' now points to the first day of the quadrennium con-
c     taining the given date.
c
110   jyrup = jyrup + 4*nquads
      nyears = jyr - jyrup
*
*     changed here 1 Jun 07
      if(nyears.ge.1 .or. (mod(jyr,100) == 0 .and. mod(jyr,400) /= 0))
     x  then
c      Common year
	 kode = 1115212
	 julian = julian + 365*nyears + kxtra
c
	   else
c      Leap year
	 kode = 1115208
	end if
c
	if(month .le. 1) go to 140
c
	do 130 j = 2,month
	julian = julian + 31 - mod(kode,4)
130   kode = kode/4
c
140   julday = julian + jday - 1
	return
c
160   write(6,161) jyear,month,jday
161   format(' Error in arguments to JULDAY ',3i12)
	julday = 0
	return
	end
c     rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Subroutine TCIVIL(dtime,  nyr,nmo,nda,nhr,min,sec,  nthday)
c
*     SFS
*     TCIVIL - to convert Double Precision time in JDN or MJDN to time
c     expressed in civil units.  Nthday varies from 1 to 366.
*     This differs from TINVER in that seconds are returned as
c     REAL*4, and hence may be fractional.
c     A four-digit year nyr is returned.
*
	Implicit None
	Double Precision dtime,dt,dfrac,dflt
	Real secs,fksecs,frac,fnsec,sec
	Integer jdt,jdtx,nyr,nmo,nda,nthday,ksecs,nhr,min,nsec
*
	character*4 cmon
cxr   julinv,dflt
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	if(dtime .le. 0.d+0) go to 800
*
	dt = dtime + .5d+0
	jdt = dt
	jdtx = jdt
	if(dt.lt.1721060.d+0)jdtx = jdtx + 2440588
*     An input smaller than 1721060 (beginning of Christian era)
c     is presumed to be a Modified Julian Day Number.
	call julinv(jdtx,  nyr,nmo,nda,  cmon,nthday)
	dfrac = dt - dflt(jdt)
	secs = 86400.d+0 * dfrac + 1.d-6
	ksecs = secs
	nhr = ksecs/3600
	min = (ksecs - 3600*nhr)/60
	nsec = ksecs - 3600*nhr - 60*min
	fksecs = ksecs
	frac = secs - fksecs
	fnsec = nsec
	sec = fnsec + frac
	if(sec.ge.60.)sec=59.99999
	return
*
800   nyr = 0
	nmo = 0
	nda = 0
	nhr = 0
	min = 0
	sec = 0
	write(*,'('' WARNING: Input to TCIVIL is bad: '',1pd20.10)')dtime
	return
	end
c     rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
	Subroutine JULINV(julin, ky,km,kd, cmon,nth)
*
*     To convert the Julian Day Number (JDN) to civil date:
*     7 Aug 98; tested using JULTEST
*     This routine works for years 4 AD to 3999 AD.
*     A 4-digit year ky is returned.  The returned ASCII month is of
*     the form ' Jan', i.e. character*4 with a lead blank.  nth
*     is the day of the year (1-366).
*
	Implicit None
	Integer jdn,julin,ndce,kwads,nrem,long,ncents,kq0001,nquads
	Integer kq0002,kyear,kq0003,kym1,leap,mod,kode,nth,mm,m
	Integer ndm,kd,km,ky
	character*4 cmon,cmons(12)
	integer jqmax/z'7fffffff'/
*
	data cmons/' Jan',' Feb',' Mar',' Apr',' May',' Jun',' Jul',
     1  ' Aug',' Sep',' Oct',' Nov',' Dec'/
c     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	jdn = julin
	if(jdn.lt.1721060 .or. jdn.gt.3182410) go to 900
*
	ndce = jdn - 1721060
	kwads = ndce/146097
	nrem = ndce - 146097*kwads
	long = 1
	ncents = 0
*
	do   2 kq0001 = 1,jqmax
*      nrem = & - (36524 + long)
	 nrem = nrem-(36524+long)
	 if(nrem.lt.0) go to   4
*
	 long = 0
2     ncents = ncents+1
*
4     nrem = nrem + (36524 + long)
	nquads = 0
*
	do   6 kq0002 = 1,jqmax
	nrem = nrem-(1460+long)
	if(nrem.lt.0) go to   8
*
	 long = 1
6     nquads = nquads+1
*
8     nrem = nrem + (1460 + long)
	kyear = 400*kwads + 100*ncents + 4*nquads
	long = 1
*
	do  10 kq0003 = 1,jqmax
	nrem = nrem-(365+long)
	if(nrem.lt.0) go to  12
*
	 long = 0
10    kyear = kyear+1
*
12    nrem = nrem+(365+long)
	kym1 = kyear - 1
	leap = (mod(kym1,4)+1)/4 - (mod(kym1,100)+1)/100
     1  + (mod(kym1,400)+1)/400
	kode = 1115212 - 4*leap
	nth = nrem + 1
*
	do  14 m = 1,12
	 mm = m
	 ndm = 31 - mod(kode,4)
	 nrem = nrem-ndm
	 if(nrem.lt.0) go to  16
*
14    kode = kode/4
*
16    nrem = nrem+ndm
	kd = nrem + 1
	km = mm
	cmon = cmons(mm)
	ky = kyear
	return
*
900   write(*, '('' JULINV: JDN out of range:'', i14)') jdn
	return
	end
