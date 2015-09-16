module cx_pfaast_tools_mod

   use cx_pfaast_constants_mod

   contains
   
!ccccccccc
   subroutine calpir(t_avg_ref,  amt_wet_ref, amt_ozo_ref, &
                             t_avg,      amt_wet,     amt_ozo, &
                             p_avg,      sec_theta,   n_layers, &
                             n_dry_pred, n_wet_pred,  n_ozo_pred, &
                             n_con_pred, &
                             pred_dry,   pred_wet,    pred_ozo, &
                             pred_con)
!c ... version of 18.03.03

!c  PURPOSE:

!c    Routine to calculate the predictors for the dry (temperature), 
!c      wet and ozone components of a fast transmittance model for a
!c      scanning satellite based instrument.

!c  REFERENCES:

!c    AIRS FTC package science notes and software, S. Hannon and L. Strow,
!c      Uni. of Maryland, Baltimore County (UMBC)

!c  CREATED:

!    19-Sep-1996 HMW

!  ARGUMENTS:

!      Input
!    -----------
!     t_avg_ref  - REAL*4 reference layer average temperature array (K)

!    amt_wet_ref - REAL*4 reference water vapour amount array (k.mol)/cm^2

!    amt_ozo_ref - REAL*4 reference ozone amount array (k.mol)/cm^2

!      t_avg     - REAL*4 layer average temperature array (K)

!     amt_wet    - REAL*4 water vapour amount array (k.mol)/cm^2

!     amt_ozo    - REAL*4 ozone amount array (k.mol)/cm^2

!      p_avg     - REAL*4 layer average pressure array (mb)

!    sec_theta   - REAL*4 secant of the zenith angle array

!     n_layers   - INT*4 Number of atmospheric layers

!    n_dry_pred  - INT*4 number of dry (temperature) predictors

!    n_wet_pred  - INT*4 number of water vapour predictors

!    n_ozo_pred  - INT*4 number of ozone predictors

!    n_con_pred  - INT*4 number of water vapour continuum predictors

!      Output
!    -----------
!     pred_dry   - REAL*4 dry gas (temperature) predictor matrix

!     pred_wet   - REAL*4 water vapour predictor matrix

!     pred_ozo   - REAL*4 ozone predictor matrix

!     pred_con   - REAL*4 water vapour continuum predictor matrix

!  COMMENTS:

!    Levels or Layers?
!    -----------------
!      Profile data is input at a number of *LAYERS*.

!    Layer Numbering pt. A
!    ---------------------
!      Layer 1   => Atmosphere between LEVELs 1 & 2
!      Layer 2   => Atmosphere between LEVELs 2 & 3
!                        .
!                        .
!                        .
!      Layer L-1 => Atmosphere between LEVELs L-1 & L

!    Layer Numbering pt. B
!    ---------------------
!      For the HIS instrument, Layer 1 is at the top of the atmosphere
!        and Layer L-1 is at the surface.    

!    Layer Numbering pt. C
!    ---------------------
!      In this routine the number of *LAYERS* is passed in the argument
!        list, _not_ the number of LEVELS.  This was done to improve
!        the readability of this code, i.e. loop from 1->L(ayers) 
!        rather than from 1->L(evels)-1.

!=======================================================================


      implicit none

!------------------------------------------------------------------------
!                             Arguments
!------------------------------------------------------------------------

! -- Input

      integer*4 n_layers, &
               n_dry_pred, n_wet_pred, n_ozo_pred, n_con_pred

      real*4    t_avg_ref(*), amt_wet_ref(*), amt_ozo_ref(*), &
               t_avg(*),     amt_wet(*),     amt_ozo(*), &
               p_avg(*),     sec_theta(*)

! -- Output

      real*4    pred_dry(n_dry_pred, *), &
               pred_wet(n_wet_pred, *), &
               pred_ozo(n_ozo_pred, *), &
               pred_con(n_con_pred, *)

!------------------------------------------------------------------------
!                           Local variables
!------------------------------------------------------------------------

! -- Parameters

      integer*4 max_layers
      parameter ( max_layers = 100 )

      integer*4 max_dry_pred, max_wet_pred, max_ozo_pred, max_con_pred
      parameter ( max_dry_pred = 8, &
                 max_wet_pred = 13, &
                 max_ozo_pred = 9, &
                 max_con_pred = 4 )

! -- Scalars

      integer*4 l

! -- Arrays

!     ....Pressure
      real*4    p_dp(max_layers), &
              p_norm(max_layers)

!     ....Temperature
      real*4    delta_t(max_layers), &
               t_ratio(max_layers), &
               pw_t_ratio(max_layers)      ! Pressure weighted

!     ....Water vapour
      real*4    wet_ratio(max_layers), &
               pw_wet(max_layers),     &    ! Pressure weighted
               pw_wet_ref(max_layers), &   ! Pressure weighted
               pw_wet_ratio(max_layers)    ! Pressure weighted

!     ....Ozone
      real*4    ozo_ratio(max_layers),    &
               pw_ozo_ratio(max_layers), &  ! Pressure weighted
               pow_t_ratio(max_layers)     ! Pressure/ozone weighted

!************************************************************************
!                         ** Executable code **
!************************************************************************

!------------------------------------------------------------------------
!                   -- Check that n_layers is o.k. --
!------------------------------------------------------------------------

      if( n_layers .gt. max_layers )then
        write(*,'(/10x,''*** calpir : n_layers > max_layers'')')
        stop
      end if 

!------------------------------------------------------------------------
!         -- Check that numbers of predictors is consistent --
!------------------------------------------------------------------------

!     ---------------------------------
!     # of dry (temperature) predictors
!     ---------------------------------

      if( n_dry_pred .ne. max_dry_pred )then
        write(*,'(/10x,''*** calpir : invalid n_dry_pred'')')
        stop
      end if 

!     ----------------------------
!     # of water vapour predictors
!     ----------------------------

      if( n_wet_pred .ne. max_wet_pred )then
        write(*,'(/10x,''*** calpir : invalid n_wet_pred'')')
        stop
      end if 

!     ---------------------
!     # of ozone predictors
!     ---------------------

      if( n_ozo_pred .ne. max_ozo_pred )then
        write(*,'(/10x,''*** calpir : invalid n_ozo_pred'')')
        stop
      end if 

!     --------------------------------------
!     # of water vapour continuum predictors
!     --------------------------------------

      if( n_con_pred .ne. max_con_pred )then
        write(*,'(/10x,''*** calpir : invalid n_con_pred'')')
        stop
      end if 

!------------------------------------------------------------------------
!         -- Calculate ratios, offsets, etc, for top layer --
!------------------------------------------------------------------------

!     ------------------
!     Pressure variables
!     ------------------

      p_dp(1)   = p_avg(1) * ( p_avg(2) - p_avg(1) )
      p_norm(1) = 0.0

!     ---------------------
!     Temperature variables
!     ---------------------

      delta_t(1)    = t_avg(1) - t_avg_ref(1)
      t_ratio(1)    = t_avg(1) / t_avg_ref(1)
      pw_t_ratio(1) = 0.0

!     ----------------
!     Amount variables
!     ----------------

!     ....Water vapour
 
      wet_ratio(1)    = amt_wet(1) / amt_wet_ref(1)
      pw_wet(1)       = p_dp(1) * amt_wet(1)
      pw_wet_ref(1)   = p_dp(1) * amt_wet_ref(1)
      pw_wet_ratio(1) = wet_ratio(1)

!     ....Ozone

      ozo_ratio(1)    = amt_ozo(1) / amt_ozo_ref(1)
      pw_ozo_ratio(1) = 0.0
      pow_t_ratio(1)  = 0.0

!------------------------------------------------------------------------
!         -- Calculate ratios, offsets, etc, for all layers --
!------------------------------------------------------------------------

      do l = 2, n_layers

!       ------------------
!       Pressure variables
!       ------------------

        p_dp(l) = p_avg(l) * ( p_avg(l) - p_avg(l-1) )
        p_norm(l) = p_norm(l-1) + p_dp(l)

!       ---------------------
!       Temperature variables
!       ---------------------

        delta_t(l)    = t_avg(l) - t_avg_ref(l)
        t_ratio(l)    = t_avg(l) / t_avg_ref(l)
        pw_t_ratio(l) = pw_t_ratio(l-1) + ( p_dp(l) * t_ratio(l-1) )

!       ----------------
!       Amount variables
!       ----------------

!       ..Water vapour

        wet_ratio(l)  = amt_wet(l) / amt_wet_ref(l)
        pw_wet(l)     = pw_wet(l-1) + ( p_dp(l) * amt_wet(l) )
        pw_wet_ref(l) = pw_wet_ref(l-1) + ( p_dp(l) * amt_wet_ref(l) )
        
!       ..Ozone

        ozo_ratio(l)    = amt_ozo(l) / amt_ozo_ref(l)
        pw_ozo_ratio(l) = pw_ozo_ratio(l-1) + &
                           ( p_dp(l) * ozo_ratio(l-1) )
        pow_t_ratio(l)  = pow_t_ratio(l-1) + &
                           ( p_dp(l) * ozo_ratio(l-1) * delta_t(l-1) )

      end do

!------------------------------------------------------------------------
!              -- Scale the pressure dependent variables --
!------------------------------------------------------------------------

      do l = 2, n_layers

        pw_t_ratio(l)   = pw_t_ratio(l) / p_norm(l)
        pw_wet_ratio(l) = pw_wet(l) / pw_wet_ref(l)
        pw_ozo_ratio(l) = pw_ozo_ratio(l) / p_norm(l)
        pow_t_ratio(l)  = pow_t_ratio(l) / p_norm(l)
 
      end do

!------------------------------------------------------------------------
!                     -- Load up predictor arrays --
!------------------------------------------------------------------------

      do l = 1, n_layers

!       ----------------------
!       Temperature predictors
!       ----------------------
         
        pred_dry(1,l) = sec_theta(l)
        pred_dry(2,l) = sec_theta(l) * sec_theta(l)
        pred_dry(3,l) = sec_theta(l) * t_ratio(l)
        pred_dry(4,l) = pred_dry(3,l) * t_ratio(l)
        pred_dry(5,l) = t_ratio(l)
        pred_dry(6,l) = t_ratio(l) * t_ratio(l)
        pred_dry(7,l) = sec_theta(l) * pw_t_ratio(l)
        pred_dry(8,l) = pred_dry(7,l) / t_ratio(l) 

!       -----------------------
!       Water vapour predictors
!       -----------------------

        pred_wet(1,l)  = sec_theta(l) * wet_ratio(l)
        pred_wet(2,l)  = sqrt( pred_wet(1,l) )
        pred_wet(3,l)  = pred_wet(1,l) * delta_t(l)
        pred_wet(4,l)  = pred_wet(1,l) * pred_wet(1,l)
        pred_wet(5,l)  = abs( delta_t(l) ) * delta_t(l) * pred_wet(1,l)
        pred_wet(6,l)  = pred_wet(1,l) * pred_wet(4,l)
        pred_wet(7,l)  = sec_theta(l) * pw_wet_ratio(l)
        pred_wet(8,l)  = pred_wet(2,l) * delta_t(l)
        pred_wet(9,l)  = sqrt( pred_wet(2,l) )
        pred_wet(10,l) = pred_wet(7,l) * pred_wet(7,l)
        pred_wet(11,l) = sqrt( pred_wet(7,l) )
        pred_wet(12,l) = pred_wet(1,l)
! +++ old
!       pred_wet(13,l) = pred_wet(2,l)
! +++ new
        pred_wet(13,l) = pred_wet(1,l) / pred_wet(10,l)

!       ----------------
!       Ozone predictors
!       ----------------

        pred_ozo(1,l) = sec_theta(l) * ozo_ratio(l)
        pred_ozo(2,l) = sqrt( pred_ozo(1,l) )
        pred_ozo(3,l) = pred_ozo(1,l) * delta_t(l)
        pred_ozo(4,l) = pred_ozo(1,l) * pred_ozo(1,l)
        pred_ozo(5,l) = pred_ozo(2,l) * delta_t(l)
        pred_ozo(6,l) = sec_theta(l) * pw_ozo_ratio(l)
        pred_ozo(7,l) = sqrt( pred_ozo(6,l) ) * pred_ozo(1,l)
        pred_ozo(8,l) = pred_ozo(1,l) * pred_wet(1,l)
        pred_ozo(9,l) = sec_theta(l) * pow_t_ratio(l) * pred_ozo(1,l)

!       ---------------------------------
!       Water vapour continuum predictors
!       ---------------------------------

        pred_con(1,l) = sec_theta(l) * wet_ratio(l) / &
                       ( t_ratio(l) * t_ratio(l) )  
        pred_con(2,l) = pred_con(1,l) * pred_con(1,l) / sec_theta(l)
        pred_con(3,l) = sec_theta(l) * wet_ratio(l) / t_ratio(l) 
        pred_con(4,l) = pred_con(3,l) * wet_ratio(l)

      end do
         
      return
      end subroutine calpir
!!ccccccc
      subroutine conpir( p, t, w, o, n_levels, i_dir, &
                              p_avg, t_avg, w_amt, o_amt)
! ... version of 19.09.96

!  PURPOSE:

!    Function to convert atmospheric water vapour (g/kg) and ozone (ppmv)
!      profiles specified at n_levels layer BOUNDARIES to n_levels-1
!      integrated layer amounts of units (k.moles)/cm^2.  The average
!      LAYER pressure and temperature are also returned.

!  REFERENCES:

!    AIRS LAYERS package science notes, S. Hannon and L. Strow, Uni. of
!      Maryland, Baltimore County (UMBC)

!  CREATED:

!    19-Sep-1996 HMW

!  ARGUMENTS:

!     Input
!    --------
!       p     - REAL*4 pressure array (mb)

!       t     - REAL*4 temperature profile array (K)

!       w     - REAL*4 water vapour profile array (g/kg)

!       o     - REAL*4 ozone profile array (ppmv)

!    n_levels - INT*4 number of elements used in passed arrays

!     i_dir   - INT*4 direction of increasing layer number

!                 i_dir = +1, Level(1) == p(top)         } satellite/AC
!                             Level(n_levels) == p(sfc)  }    case

!                 i_dir = -1, Level(1) == p(sfc)         } ground-based
!                             Level(n_levels) == p(top)  }    case

!     Output
!    --------
!     p_avg   - REAL*4 average LAYER pressure array (mb)

!     t_avg   - REAL*4 average LAYER temperature (K)

!     w_amt   - REAL*4 integrated LAYER water vapour amount array (k.moles)/cm^2

!     o_amt   - REAL*4 integrated LAYER ozone amount array (k.moles)/cm^2

!  ROUTINES:

!    Subroutines:
!    ------------
!      gphite      - calculates geopotential height given profile data.

!    Functions:
!    ----------
!      NONE

!  COMMENTS:

!    Levels or Layers?
!    -----------------
!      Profile data is input at a number of *LEVELS*.  Number densitites
!        are calculated for *LAYERS* that are bounded by these levels.
!        So, for L levels there are L-1 layers.

!    Layer Numbering
!    ---------------
!      Layer 1   => Atmosphere between LEVELs 1 & 2
!      Layer 2   => Atmosphere between LEVELs 2 & 3
!                        .
!                        .
!                        .
!      Layer L-1 => Atmosphere between LEVELs L-1 & L

!=======================================================================

!-----------------------------------------------------------------------
!              -- Prevent implicit typing of variables --
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!                           -- Arguments --
!-----------------------------------------------------------------------

! -- Arrays

      real*4    p(*), t(*), w(*), o(*), &
               p_avg(*), t_avg(*), w_amt(*), o_amt(*)

! -- Scalars

      integer*4 n_levels, i_dir

!-----------------------------------------------------------------------
!                         -- Local variables --
!-----------------------------------------------------------------------

! -- Parameters

      integer*4 max_levels
      parameter ( max_levels = 101 )         ! Maximum number of layers

      real*4    r_equator, r_polar, r_avg
      parameter ( r_equator = 6.378388e+06,  & ! Earth radius at equator
                 r_polar   = 6.356911e+06,  & ! Earth radius at pole
                 r_avg     = 0.5*(r_equator+r_polar) )

      real*4    g_sfc
      parameter ( g_sfc = 9.80665 )         ! Gravity at surface

      real*4    rho_ref
      parameter ( rho_ref = 1.2027e-12 )    ! Reference air "density"

      real*4    mw_dryair, mw_h2o, mw_o3
      parameter ( mw_dryair = 28.97,    &    ! Molec. wgt. of dry air (g/mol)
                 mw_h2o    = 18.0152,   &    ! Molec. wgt. of water
                 mw_o3     = 47.9982 )     ! Molec. wgt. of ozone

      real*4    R_gas, R_air
      parameter ( R_gas = 8.3143,        &   ! Ideal gas constant (J/mole/K)
                 R_air = 0.9975*R_gas )    ! Gas constant for air (worst case) 

! -- Scalars

      integer*4 l, l_start, l_end, l_indx

      real*4    rho1, rho2, p1, p2, w1, w2, o1, o2, z1, z2, &
               c_avg, g_avg, z_avg, w_avg, o_avg, &
               dz, dp, r_hgt, wg, og, A, B

! -- Arrays

      real*4    z(max_levels),        &      ! Pressure heights (m)
               g(max_levels),         &     ! Acc. due to gravity (m/s/s)
               mw_air(max_levels),     &    ! Molec. wgt. of air (g/mol)
               rho_air(max_levels),   &     ! air mass density (kg.mol)/m^3
               c(max_levels),          &    ! (kg.mol.K)/(N.m)
               w_ppmv(max_levels)          ! h2o LEVEL amount (ppmv)

!***********************************************************************
!                         ** Executable code **
!***********************************************************************

!-----------------------------------------------------------------------
!           -- Calculate initial values of pressure heights --
!-----------------------------------------------------------------------

      call gphite( p, t, w, 0.0, n_levels, i_dir, z)

!-----------------------------------------------------------------------
!      -- Set loop bounds for direction sensitive calculations --
!      -- so loop iterates from surface to the top             --
!-----------------------------------------------------------------------

      if( i_dir .gt. 0 )then

!       --------------------
!       Data stored top down
!       --------------------

        l_start = n_levels
        l_end   = 1

      else

!       ---------------------
!       Data stored bottom up
!       ---------------------

        l_start = 1
        l_end   = n_levels

      end if

!-----------------------------------------------------------------------
!          -- Air molecular mass and density, and gravity --
!          -- as a function of LEVEL                      --
!-----------------------------------------------------------------------

!     -----------------------
!     Loop from bottom to top
!     -----------------------

      do l = l_start, l_end, -1*i_dir

!       ---------------------------------
!       Convert water vapour g/kg -> ppmv
!       ---------------------------------

        w_ppmv(l) = 1.0e+03 * w(l) * mw_dryair / mw_h2o

!       -----------------------------------------
!       Calculate molecular weight of air (g/mol)
!       ----------------------------------------

        mw_air(l) = ( ( 1.0 - (w_ppmv(l)/1.0e+6) ) * mw_dryair ) + &
                  ( ( w_ppmv(l)/1.0e+06 ) * mw_h2o )

!       ----------------
!       Air mass density
!       ----------------

        c(l) = 0.001 * mw_air(l) / R_air    ! 0.001 factor for g -> kg
        rho_air(l) = c(l) * p(l) / t(l)

!       -------
!       Gravity
!       -------

        r_hgt = r_avg + z(l)                !  m
        g(l) = g_sfc -               &       !  m/s^2
              g_sfc*( 1.0 - ( (r_avg*r_avg)/(r_hgt*r_hgt) ) )

      end do
 
!-----------------------------------------------------------------------
!                        -- LAYER quantities --
!-----------------------------------------------------------------------

!     -----------------------
!     Loop from bottom to top
!     -----------------------

      do l = l_start, l_end+i_dir, -1*i_dir

!       -------------------------------------------------------
!       Determine output array index.  This is done so that the
!       output data is always ordered from 1 -> L-1 regardless
!       of the orientation of the input data.  This is true by
!       default only for the bottom-up case.  For the top down
!       case no correction would give output layers from 2 -> L
!       -------------------------------------------------------

        if( i_dir .gt. 0 )then

          l_indx = l - 1

        else

          l_indx = l

        end if

!       ---------------------------------------
!       Assign current layer boundary densities
!       ---------------------------------------
 
        rho1 = rho_air(l)
        rho2 = rho_air(l-i_dir)
 
!       ---------
!       Average c
!       ---------

        c_avg = ( (rho1*c(l)) + (rho2*c(l-i_dir)) ) / ( rho1 + rho2 )

!       ---------
!       Average t
!       ---------

        t_avg(l_indx) =  &
               ( (rho1*t(l)) + (rho2*t(l-i_dir)) ) / ( rho1 + rho2 )

!       ---------
!       Average p
!       ---------

        p1 = p(l)
        p2 = p(l-i_dir)

        z1 = z(l)
        z2 = z(l-i_dir)

        dp = p2 - p1

        A = log(p2/p1) / (z2-z1)
        B = p1 / exp(A*z1)

        p_avg(l_indx) = dp / log(p2/p1)

!       ------------------------------------------------
!       LAYER thickness (rather long-winded as it is not
!       assumed the layers are thin) in m. Includes
!       correction for altitude/gravity.
!       ------------------------------------------------

!       ...Initial values
        g_avg = g(l)
        dz = -1.0 * dp * t_avg(l_indx) / ( g_avg*c_avg*p_avg(l_indx) )

!       ...Calculate z_avg
        z_avg = z(l) + ( 0.5*dz )

!       ...Calculate new g_avg
        r_hgt = r_avg + z_avg 
        g_avg = g_sfc - g_sfc*( 1.0 - ( (r_avg*r_avg)/(r_hgt*r_hgt) ) )

!       ...Calculate new dz
        dz = -1.0 * dp * t_avg(l_indx) / ( g_avg*c_avg*p_avg(l_indx) )

!       ----------------------------------------
!       Calculate LAYER amounts for water vapour
!       ----------------------------------------

        w1 = w_ppmv(l)
        w2 = w_ppmv(l-i_dir)

        w_avg =  ( (rho1*w1) + (rho2*w2) ) / ( rho1+rho2 )

        w_amt(l_indx) =  &
           rho_ref * w_avg * dz * p_avg(l_indx) / t_avg(l_indx)

!       ---------------------------------
!       Calculate LAYER amounts for ozone
!       ---------------------------------

        o1 = o(l)
        o2 = o(l-i_dir)

        o_avg =  ( (rho1*o1) + (rho2*o2) ) / ( rho1+rho2 )

        o_amt(l_indx) =  &
           rho_ref * o_avg * dz * p_avg(l_indx) / t_avg(l_indx)

      end do

      return
      end subroutine conpir
!!ccccccc
      subroutine gphite( p, t, w, z_sfc, n_levels, i_dir, z)
! .... version of 18.05.00

! PURPOSE:

!  Routine to compute geopotential height given the atmospheric state.
!    Includes virtual temperature adjustment.

! CREATED:

!  19-Sep-1996 Received from Hal Woolf, recoded by Paul van Delst
!  18-May-2000 Logic error related to z_sfc corrected by Hal Woolf

!  ARGUMENTS:

!     Input
!    --------
!       p     - REAL*4 pressure array (mb)

!       t     - REAL*4 temperature profile array (K)

!       w     - REAL*4 water vapour profile array (g/kg)

!     z_sfc   - REAL*4 surface height (m).  0.0 if not known.

!    n_levels - INT*4 number of elements used in passed arrays

!     i_dir   - INT*4 direction of increasing layer number

!                 i_dir = +1, Level(1) == p(top)         } satellite/AC
!                             Level(n_levels) == p(sfc)  }    case

!                 i_dir = -1, Level(1) == p(sfc)         } ground-based
!                             Level(n_levels) == p(top)  }    case

!     Output
!    --------
!       z     - REAL*4 pressure level height array (m)

! COMMENTS:

!   Dimension of height array may not not be the same as that of the
!     input profile data.

!=======================================================================

!-----------------------------------------------------------------------
!              -- Prevent implicit typing of variables --
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!                           -- Arguments --
!-----------------------------------------------------------------------

! -- Arrays

      real*4    p(*), t(*), w(*), &
               z(*)

! -- Scalars

      integer*4 n_levels, i_dir

      real*4    z_sfc

!-----------------------------------------------------------------------
!                         -- Local variables --
!-----------------------------------------------------------------------

! -- Parameters

      real*4    rog, fac
      parameter ( rog = 29.2898, &
                 fac = 0.5 * rog )

! -- Scalars

      integer*4 i_start, i_end, l

      real*4    v_lower, v_upper, algp_lower, algp_upper, hgt

!***********************************************************************
!                         ** Executable code **
!***********************************************************************

!-----------------------------------------------------------------------
!  -- Calculate virtual temperature adjustment and exponential       --
!  -- pressure height for level above surface.  Also set integration --
!  -- loop bounds                                                    --
!-----------------------------------------------------------------------

      if( i_dir .gt. 0 )then

!       --------------------
!       Data stored top down
!       --------------------

        v_lower = t(n_levels) * ( 1.0 + ( 0.00061 * w(n_levels) ) )

        algp_lower = alog( p(n_levels) )

        i_start = n_levels-1
        i_end   = 1

      else

!       ---------------------
!       Data stored bottom up
!       ---------------------

        v_lower = t(1) * ( 1.0 + ( 0.00061 * w(1) ) )

        algp_lower = alog( p(1) )

        i_start = 2
        i_end   = n_levels

      end if

!-----------------------------------------------------------------------
!                     -- Assign surface height --
!-----------------------------------------------------------------------

      hgt = z_sfc

! .. Following added 18 May 2000 ... previously, z(n_levels) for downward
!       (usual) case was not defined!

      if(i_dir.gt.0) then
         z(n_levels) = z_sfc
      else
         z(1) = z_sfc
      endif

! .. End of addition

!-----------------------------------------------------------------------
!             -- Loop over layers always from sfc -> top --
!-----------------------------------------------------------------------

      do l = i_start, i_end, -1*i_dir

!       ----------------------------------------------------
!       Apply virtual temperature adjustment for upper level
!       ----------------------------------------------------

        v_upper = t(l)
        if( p(l) .ge. 300.0 ) &
         v_upper = v_upper * ( 1.0 + ( 0.00061 * w(l) ) )

!       ----------------------------------------------------- 
!       Calculate exponential pressure height for upper layer
!       ----------------------------------------------------- 

        algp_upper = alog( p(l) )

!       ----------------
!       Calculate height
!       ----------------

        hgt = hgt + ( fac*(v_upper+v_lower)*(algp_lower-algp_upper) )

!       -------------------------------
!       Overwrite values for next layer
!       -------------------------------

        v_lower = v_upper
        algp_lower = algp_upper

!       ---------------------------------------------
!       Store heights in same direction as other data
!       ---------------------------------------------

        z(l) = hgt

      end do

      return
      end subroutine gphite

!>
!!
!! @history  made it f90 ( 02/10/2015 AW)
!
subroutine taudoc(cc,xx,tau)
   
      implicit none
      
      real , intent(in) :: cc(:,:)
      real , intent(in) :: xx(:,:)
      real, intent(out)  :: tau (:)
      
      real :: trap = -999. 
      
! * Strow-Woolf model ... for dry, ozo(ne), and wco (water-vapor continuum)
! .... version of 05.09.02

      integer :: n_lay
      integer :: i , j
     
      integer :: nc , nx
      real ::  taulyr
      real :: yy
      
      n_lay = ubound ( cc, 2)
      nc = ubound ( cc,1 )
      nx = ubound ( xx, 1)

      tau    = 0.
      tau(1) = 1.
      taulyr = 1.
    
  !  -loop over all layers from top 
   do j = 1 , n_lay
      ! assume background stored in last column of coeffecients
      yy = cc(nc,j)
      
      yy = yy + DOT_PRODUCT( cc(:nc-1,j), xx(:nc-1,j) )
      if ( yy > 0. ) taulyr = exp ( -yy)
      tau ( j + 1) = tau ( j ) * taulyr
      
   end do
   
   end subroutine taudoc
   !
   !
   !
   ! * Strow-Woolf model ... for 'wet' (water-vapor other than continuum)
   ! .... version of 05.09.02
subroutine tauwtr( ccs , ccl , xx , tau )
      implicit none
   
      real, intent(in) , target :: ccs (:,:)
      real, intent(in) , target :: ccl (:,:)
      real, intent(in) :: xx(:,:)
      real, intent(out)  :: tau (:)
       
      integer :: nc, nx
      real :: od_sum 
      real, pointer :: cc (:,:) => null()
      real ::  taulyr
      real :: yy
      integer :: i , j , n_lay

      n_lay = ubound ( ccs, 2)
      nc = ubound ( ccs,1 )
      nx = ubound ( xx, 1)
      
      tau    = 0.
      tau(1) = 1.
      taulyr = 1.
      
      od_sum = 0.
      
      cc => ccs
      
      do j = 1 , n_lay


         ! differentiate between ice and water
!---old
!        if ( od_sum >= 5.0 ) cc => ccl
!--- heidinger fix begin
         if ( od_sum >= 5.0 ) then
               cc => ccl
               nc = ubound ( ccl,1 )
         endif
!--- heidinger fix end

         ! assume background stored in last column of coeffecients
         yy = cc(nc,j)
      
         yy = yy + DOT_PRODUCT( cc(:,j), xx(:nc-1,j) )
         od_sum=od_sum + max(yy,0.)
         if ( yy > 0. ) taulyr = exp ( -yy)
         tau ( j + 1) = tau ( j ) * taulyr   
            
      end do
      cc=> null() 
     
      return

  end subroutine tauwtr
   
   
   end module cx_pfaast_tools_mod
