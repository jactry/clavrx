!------------------------------------------------------------------
! module for the SOI model
!
! subroutines:
! 
! truncated_doubling - computes the layer R,T and S
! 
!------------------------------------------------------------------
module SOI
 use CONSTANTS
 implicit none

 public:: FORWARD_MODEL_SOI, INITIALIZE_SOI, FORWARD_MODEL_ABS_APPROX

 private:: TRUNCATED_DOUBLING,COMPUTE_ZERO_ORDER_RAD,COMPUTE_HIGHER_ORDER_RAD,&
          PLEG,GAULEG,PHASE_MATRIX, INVERT

 integer, parameter:: nstream = 8, &
                      nstreams = nstream + 2, &
                      nexp = nstream, &
                      nscat = 200

 logical, parameter:: delta = .false., &
                      accel = .true.

 real, parameter:: delta_tau_scat_iso = 0.001

 real, dimension(nstreams):: mu,w
 real, dimension(nexp,nstreams+2):: Ypleg

 contains
!-----------------------------------------------------------
! 
!
!-----------------------------------------------------------
subroutine INITIALIZE_SOI()

real(kind=real4), dimension(nstreams/2-1):: pts4, wgts4
integer:: i

!----- compute quadrature points
mu( 1 ) = 1.0
w( 1 ) = 0.0
mu( nstreams / 2 + 1 ) = -1.0
w( nstreams / 2 + 1 ) = 0.0

! Use Gauss quadrature for 2-stream
if ( nstream == 2 ) then
  mu( 2 ) = 1. / sqrt( 3.0 )
  w( 2 ) = 1.0
  mu( 4 ) = -mu( 2 )
  w( 4 ) = 1.0
else
  call GAULEG(0.0,1.0,pts4,wgts4,nstreams/2-1)
  do i = 1, nstreams / 2 - 1
    mu( i + 1 ) = pts4( i )
    w( i + 1 ) = wgts4( i )
    mu( i + nstreams / 2 + 1 ) = -pts4( i )
    w( i + nstreams / 2 + 1 ) = wgts4( i )
  enddo
end if   


!------- compute Legendre Polynomials
call PLEG(1.0,mu,0,nexp,nstreams,Ypleg)

end subroutine INITIALIZE_SOI
!--------------------------------------------------------------------------
!
! the forward model - uses SOI to compute TOA radiance
!
!  input: nlayers - number of layers
!         tau_in - vector of layer optical depths
!         wo_in - vector of layer single scatter albedoes
!         g_in - vector of layer asymmetry parameters
!         B - vector of level planck emissions
!         esfc - surface emissivity
!         Bsfc - surface planck emission
!         Bspace - TOA planck emission
!         mu_obs - observation zenith angle cosine
!  output
!         rad_toa - toa radiance along mu_obs
!--------------------------------------------------------------------------
subroutine FORWARD_MODEL_SOI(nlayers,tau_in,wo_in,g_in,B,esfc,Bsfc, &
                             Bspace,mu_obs,rad_toa)
  real(kind=real4), dimension(:), intent(in):: tau_in, wo_in, g_in, B
  real(kind=real4), intent(in)::  esfc,Bsfc,mu_obs, Bspace
  integer, intent(in):: nlayers
  real(kind=real4), intent(out):: rad_toa

  real(kind=real4), dimension(nlayers):: Bo, Bn, tau, wo, g, f
  real, dimension(nstreams/2,nlayers):: Sp,Sm
  real, dimension(nstreams/2,nstreams/2,nlayers):: Pf,Pb,R,Ts,Td
  real, dimension(nstreams/2,nlayers+1,0:nscat):: Tb_d,Tb_u
  real, dimension(nstreams/2):: S,S_prev,Tb_final,Tb_final_prev,dTb_final
  real:: g_previous,dB_dtau,zen,dTb,Tbprc,ddTb,tau_scat_iso
  integer:: lay,ndoub,imethod,iscat,iscat_final,imu_acc
  logical:: accel_conv_prev


 imethod = 0

!----------- adjust quadrature for observation angle
mu( 1 ) = mu_obs
w( 1 ) = 0.0
mu( nstreams / 2 + 1 ) = -mu_obs
w( nstreams / 2 + 1 ) = 0.0
  
!--- apply Delta scaling
 tau = tau_in
 wo = wo_in
 g = g_in

 if (delta .eqv. .true.) then
  f = g_in ** nstream
  tau = ( 1.0 - wo_in * f ) * tau_in
  wo = ( 1.0 - f ) * wo_in / ( 1.0 - f * wo_in )
  g = ( g_in - f ) / ( 1.0 - f )
 endif

!----- make plank emission profile
do lay=1,nlayers
  Bo(lay) = B(lay)
  Bn(lay) = B(lay+1)
enddo 

!---- compute phase matrices
g_previous = 0.0
Pf = 1.0
Pb = 1.0
do lay = 1, nlayers

 if ((g(lay) /= g_previous).or.(lay==1)) then
  call  PHASE_MATRIX(0,g(lay),nstreams,nexp,Ypleg,Pf(:,:,lay),Pb(:,:,lay))
  g_previous = g(lay)
 else
  Pf(:,:,lay) = Pf(:,:,lay-1)
  Pb(:,:,lay) = Pb(:,:,lay-1)
 endif

!---- compute number of doublings to perform
!tau_scat_iso = (1.0 - g(lay)) * wo(lay) * tau(lay)
 tau_scat_iso = wo(lay) * tau(lay)
 if (tau_scat_iso > delta_tau_scat_iso) then
  ndoub = nint(log(tau_scat_iso / delta_tau_scat_iso) / log(2.0))
 else
  ndoub = 0
 endif

 if (ndoub >= 0) then
     call TRUNCATED_DOUBLING(imethod,ndoub,tau(lay),g(lay),wo(lay),Pf(:,:,lay),Pb(:,:,lay),&
                             mu(1:nstreams/2),w(1:nstreams/2), &
                             Bo(lay),Bn(lay),R(:,:,lay),Ts(:,:,lay), &
                             Td(:,:,lay),Sp(:,lay),Sm(:,lay))
 endif

enddo
 
!-----------------------------------------------------------------------------
! compute zeroth order 
!-----------------------------------------------------------------------------
Tb_u = 0.0
Tb_d = 0.0
call COMPUTE_ZERO_ORDER_RAD(Td,Sp,Sm,esfc,Bsfc,Bspace,Tb_d(:,:,0),Tb_u(:,:,0))

!-----------------------------------------------------------------------------
! compute higher orders
!-----------------------------------------------------------------------------
dTb_final = 1.e+10
S_prev = Tb_u(:,1,0)
Tb_final_prev = 0.0
iscat = 1
!Tbprc = 0.1  ! Solution precision in K
Tbprc = 0.005  ! Solution in terms of relative accuracy
accel_conv_prev = .false.

imu_acc = 1

convergence_loop: do 

  call COMPUTE_HIGHER_ORDER_RAD(R,Ts,Td,esfc,Tb_d(:,:,iscat-1),Tb_u(:,:,iscat-1),&
                                Tb_d(:,:,iscat),Tb_u(:,:,iscat))
  S = S_prev + Tb_u(:,1,iscat)

  if (accel .eqv. .true.) then
    Tb_final = S + (Tb_u(:,1,iscat)**2) / (Tb_u(:,1,iscat-1) - Tb_u(:,1,iscat))
    dTb_final = abs(Tb_final - Tb_final_prev)
                                                                                                                                                               
!--- do regular check for convergence
      if (Tb_u(imu_acc,1,iscat) < Tbprc*Tb_final(imu_acc)) then
        iscat_final = iscat
        exit
      endif

!--- check for accelerated convergence
     if ( (iscat > 2) .and. (dTb_final(imu_acc) < 2.0*Tbprc*Tb_final(imu_acc) )) then
      iscat_final = iscat
      if (accel_conv_prev .eqv. .true.) then
        iscat_final = iscat
        exit
      endif
       accel_conv_prev = .true.
     else
       accel_conv_prev = .false.
     endif

    Tb_final_prev = Tb_final

     else
   if (Tb_u(imu_acc,1,iscat) < Tbprc*Tb_final(imu_acc)) then
       iscat_final = iscat
       exit
   endif
  endif

  S_prev = S
  iscat = iscat + 1

end do convergence_loop

rad_toa = Tb_final(1)

end subroutine FORWARD_MODEL_SOI

!----------------------------------------------------------------
subroutine FORWARD_MODEL_ABS_APPROX(nlayers,tau,wo,g,B,esfc,Bsfc,Bspace, &
                                    mu_obs,rad_toa)
 integer, intent(in):: nlayers
 real, dimension(:), intent(in):: tau, wo, g, B
 real, intent(in):: esfc, Bsfc, Bspace, mu_obs
 real, intent(out):: rad_toa
 real:: trans_toa, tau_abs, rad_layer
 integer:: k

trans_toa = 1.0
rad_toa = 0.0
do k = 1, nlayers
  tau_abs = (1.0-wo(k))*tau(k)
  rad_layer = 0.5*( B(k) + B(k+1)) * &
                  (1.0 - exp(-1.0*tau_abs/mu_obs))
  rad_toa = rad_toa + trans_toa * rad_layer
  trans_toa = trans_toa * exp(-1.0*tau_abs/mu_obs)
enddo   
rad_toa = rad_toa + trans_toa * esfc * Bsfc

end subroutine FORWARD_MODEL_ABS_APPROX

!----------------------------------------------------------
! A subroutine to perfrom a truncated doubling process to 
! determine the thermal source, reflection and transmission values
!
!  input:
!      rl = local reflection matrix
!      tl = local transmission matrix
!      tau = optical depth
!      wo = single scatter albedo
!      g = asymmetry parameter
!      Bo = blackbody emission at top of layer
!      Bn = blackbody emission of base of layer
!      imethod - 0 = EIGI, 1 = SS, 2 = IGI, 3 = hybrid, 9 = truth
!      
! output:
!      R = layer reflection matrix for full layer
!      Td = transmission due to direct attenuation (maybe an input)
!      Ts = transmission due to scattering (pseudo source)
!      Sp = Thermal Source upward out of layer top
!      Sm = Thermal Source downward out of layer base
!
! local:
!      Rh = layer reflection matrix for half layer
!      Th,Yh,Zh = half layer values of T,Y,Z
!      Bd = slope of planck function through layer
!      E = the identity matrix
!      m = reflectance airmass matrix
!      n = transmission airmass matrix
!      T = layer tranmission matrix for full layer
!      k - the number of doublings required
!      Gamma, Th_Gamma, Th_Gamma_Rh - arrays used in the doubling process
!
!  Author: Andrew Heidinger, NOAA/NESDIS Office of Research and Applications
!  Date: February 2004
!
!  Remaining Concerns
!   1. picking the right number of doubling steps
!   2. too many exponentials?
!   3. When we don't double, we need to include some accounting for th
!       lapse rate for Sp and Sm (not done yet)
!      
!------------------------------------------------------------------------------
subroutine TRUNCATED_DOUBLING(imethod_input,k,tau,g,wo,Pp,Pm,mu,w,Bo,Bn,R,Ts,Td,Sp,Sm)
integer, intent(in):: imethod_input,k
real, dimension(:,:),intent(in):: Pp,Pm
real, dimension(:),intent(in):: mu,w
real, intent(in):: tau,g,wo,Bo,Bn
real, dimension(:,:),intent(out):: R,Ts,Td
real, dimension(:),intent(out):: Sp,Sm
real, dimension(size(mu),size(mu)):: E,T,Rh,Th,Gamma,Th_Gamma,Th_Gamma_Rh,store
real:: rl,tl,m, n, Bd, delta_tau,exp_delta_tau_mu,gh
real, dimension(size(mu)):: Yh,Zh,Y,Z
real:: x

integer:: i,j,l,nstreams,imethod

nstreams = size(mu)

!---- compute number of doublings to perform
 delta_tau = tau/(2**k)


!--- slope of blackbody emission through layer
Bd = (Bn - Bo)/tau


!---- select method of solution
imethod = imethod_input
 if (imethod_input == 3) then
  if (wo < 0.5) then
     imethod = 1
  else
     imethod =  0
  endif
 endif

if (imethod_input == 9) then
     imethod =  0
endif

!--- split layer into 2^k layers
 do l = 1,max(k,1) 

!--  for the first layer, initialize with single scatter approx.
  if (l == 1) then 
      Yh = 0
      Rh = 0
      Th = 0

i_loop: do i = 1,nstreams     
      x = -delta_tau/mu(i)
!     exp_delta_tau_mu = 1.0 + x*(1.0 + x*(0.5 + x*(0.166667 + (x/24.0))))
      exp_delta_tau_mu = exp(x)

j_loop: do j = 1,nstreams
         E(i,j) = 0.0
         if (i == j) then
               E(i,j) = 1.0
         endif

         m = 1.0/mu(i) + 1.0/mu(j)
         n = 1.0/mu(j) - 1.0/mu(i)

         rl = 0.5 * wo * Pm(i,j) * w(j)
         tl = 0.5 * wo * Pp(i,j) * w(j)

!---- EIGI
       if (imethod == 0) then
         Rh(i,j) = rl * (1.0-exp_delta_tau_mu)
         Th(i,j) = tl * (1.0-exp_delta_tau_mu)
         if (i == j) then
            Th(i,j) = exp_delta_tau_mu + Th(i,j)
         end if

!------ SSI
       elseif (imethod == 1) then
         Rh(i,j) = rl * (1.0 - exp(-m*delta_tau)) / (mu(i) * m )
         if (i == j) then
          Th(i,j) = exp_delta_tau_mu * (E(i,j) + tl * delta_tau / mu(i))
         else 
          Th(i,j) = exp_delta_tau_mu * (tl *  &
                     (1.0 - exp(-n*delta_tau)) / (mu(i) * n ))
         end if

!------- IGI
       elseif (imethod == 2) then
         Rh(i,j) = rl * delta_tau/mu(i)
         Th(i,j) = tl * delta_tau/mu(i)
         if (i == j) then
            Th(i,j) = (1.0-delta_tau/mu(i)) + Th(i,j)
         end if

!----- hybrid
       elseif (imethod == 3) then   !hybrid of ss and odell
         if (wo >= 0.5) then
         Rh(i,j) = rl * (1.0-exp_delta_tau_mu)
         Th(i,j) = tl * (1.0-exp_delta_tau_mu)
         if (i == j) then
            Th(i,j) = exp_delta_tau_mu + Th(i,j)
         end if
         else
          Rh(i,j) = rl * (1.0 - exp(-m*delta_tau)) / (mu(i) * m )
          if (i == j) then
           Th(i,j) = exp_delta_tau_mu * (E(i,j) + tl * delta_tau / mu(i))
          else 
           Th(i,j) = exp_delta_tau_mu * (tl *  &
                      (1.0 - exp(-n*delta_tau)) / (mu(i) * n ))
          end if
        endif
       endif


!--- compute single scatter source vector (component due to scattering)
      if (imethod == 1) then
         Yh(i) = Yh(i) +           &  
                 (1.0-wo)*(rl*(1.0-exp_delta_tau_mu) - Rh(i,j)) 
         if (i == j) then
           Yh(i) = Yh(i) +           &  
                  (1.0-wo)*(tl*(1.0-exp_delta_tau_mu) - (Th(i,j)-exp_delta_tau_mu)) 
         else
            Yh(i) = Yh(i) +           &  
                (1.0-wo)*(tl*(1.0-exp_delta_tau_mu) - (Th(i,j))) 
         endif
      endif

    end do j_loop

!------ initialize the thermal source vectors
      if (imethod == 0) then
          Yh(i) = (1-wo) * ( 1.0 - exp_delta_tau_mu) !     elseif (imethod == 1) then
      elseif (imethod == 1) then
          Yh(i) = Yh(i) + (1.0-wo)*(1.0-exp_delta_tau_mu)   !add in zeroth order emission
      elseif (imethod == 2) then
          Yh(i) = (1-wo) * delta_tau/ mu(i)
      elseif (imethod == 3) then
         if (wo < 0.5) then
          Yh(i) = Yh(i) + (1.0-wo)*(1.0-exp_delta_tau_mu)   !add in zeroth order emission
         else
          Yh(i) = (1-wo) * ( 1.0 - exp_delta_tau_mu)
         endif
      endif

      Zh(i) = 0.0
      gh = delta_tau / 2.0

    end do i_loop

!--- alternate Y initialization based on energy converation
!    if (imethod == 3) then
!     do i = 1,nstreams
!      Yh(i) = 0.0
!      do j = 1, nstreams 
!         Yh(i) = Yh(i) + Rh(i,j) + Th(i,j)
!      enddo
!      Yh(i) = 1.0 - Yh(i)
!    enddo
!   endif

   endif

!--- if this is not the first doubling, start with previous results
   if (l > 1) then 
    Rh = R
    Th = T
    Yh = Y
    Zh = Z
    gh = 2.0 * gh
   endif

!----- if there is no doubling, store final results
   if (k == 0) then
     R = Rh
     T = Th
     Y = Yh
     Z = Zh
   else 

!--- if there is doubling, perform doubling

    store = matmul(Rh,Rh)

!--- truncated doubling
    Gamma = E + store   + matmul(store,store) !approximation for the (E - Rh*Rh)^-1

!--- a full doubling
if (imethod_input == 9) then
    Gamma = E - store
    call INVERT(Gamma,nstreams)
endif

    Th_Gamma = matmul(Th , Gamma)
    Th_Gamma_Rh = matmul(Th_Gamma , Rh)

    R = matmul(Th_Gamma_Rh , Th) + Rh
    T = matmul(Th_Gamma , Th)
    Y = matmul((Th_Gamma + Th_Gamma_Rh + E) , Yh)
    Z = Zh + gh*Yh + matmul((Th_Gamma - Th_Gamma_Rh) , (Zh - gh*Yh))
   end if

   end do

!--------------------------------------------------
!  finish thermal source computation
!-------------------------------------------------
    Sp = 0.5*(Bo + Bn)*Y - Bd*Z
    Sm = 0.5*(Bo + Bn)*Y + Bd*Z

!---------- compute direct transmission
Td = 0.0
do i = 1,nstreams
 Td(i,i) = exp(-tau/mu(i))
end do

!--- subtract off direct from total for scattering transmission
Ts = T - Td

end subroutine TRUNCATED_DOUBLING

!-------------------------------------------------------------------------
! Subroutine to compute to zeroth order intensity using the quantities 
!  computed from the TRUNCATED_DOUBLING routine
!
! layer 1 is at the top of the atmoshere
!--------------------------------------------------------------------------
subroutine COMPUTE_ZERO_ORDER_RAD(Td,Sp,Sm,esfc,Bsfc,Bspace,Tb_d_0,Tb_u_0)
real, dimension(:,:,:), intent(in):: Td
real, dimension(:,:),intent(in):: Sp,Sm
real, intent(in):: esfc,Bsfc,Bspace
real, dimension(:,:),intent(out):: Tb_d_0,Tb_u_0
integer:: nstreams, nlevels, l

nstreams = size(Sp,1)
nlevels = size(Sp,2)  + 1

!--- upwelling radiance - start at surface and
Tb_u_0(:,nlevels) =  esfc * Bsfc
do l = nlevels-1,1, -1 
Tb_u_0(:,l) = Sp(:,l) + matmul(Td(:,:,l),Tb_u_0(:,l+1))
end do

!--- downwelling radiance - start at toa
Tb_d_0(:,1) =  Bspace
do l = 2,nlevels 
Tb_d_0(:,l) = Sm(:,l-1) + matmul(Td(:,:,l-1),Tb_d_0(:,l-1))
end do

end subroutine COMPUTE_ZERO_ORDER_RAD

!------------------------------------------------------------------
! compute radiance due to higher order scatter
!
! input:  R,Ts,Td - described above
!         esfc - surface emissivity
!         Tb_d_0 - the downwelling radiance for previous order of scat.
!         Tb_u_0 - the upwelling radiance for previous order of scat.
! output:
!         Tb_d_1 - the downwelling radiance for next order of scat.
!         Tb_u_1 - the upwelling radiance for next order of scat.
!------------------------------------------------------------------
subroutine COMPUTE_HIGHER_ORDER_RAD(R,Ts,Td,esfc,Tb_d_0,Tb_u_0,Tb_d_1,Tb_u_1)
real, dimension(:,:,:),intent(in):: R,Ts,Td
real, dimension(:,:),intent(in):: Tb_d_0,Tb_u_0
real, intent(in):: esfc
real, dimension(:,:),intent(out):: Tb_d_1,Tb_u_1
integer:: nstreams, nlevels, l
nstreams = size(Tb_d_0,1)
nlevels = size(Tb_d_0,2)


!-----upwelling, start from bottom
Tb_u_1(:,nlevels) = (1.0-esfc)*Tb_d_0(:,nlevels)
do l = nlevels-1,1,-1
  Tb_u_1(:,l) = matmul(R(:,:,l),Tb_d_0(:,l)) +  &     
                matmul(Ts(:,:,l),Tb_u_0(:,l+1)) + &
                matmul(Td(:,:,l),Tb_u_1(:,l+1))

enddo

!-----upwelling, start from top
Tb_d_1(:,1) = 0.0
do l = 2,nlevels
  Tb_d_1(:,l) = matmul(R(:,:,l-1),Tb_u_0(:,l)) +  &     
                matmul(Ts(:,:,l-1),Tb_d_0(:,l-1)) + &
                matmul(Td(:,:,l-1),Tb_d_1(:,l-1))
enddo

end subroutine COMPUTE_HIGHER_ORDER_RAD

subroutine PHASE_MATRIX(m,g,Nzen,Nexp,Ypleg,Pf,Pb)
 integer, intent(in):: m,Nzen,Nexp
 real, intent(in):: g
 real, dimension(:,:), intent(in):: Ypleg
 real, dimension(:,:), intent(out):: Pf, Pb
 real, dimension(Nexp):: chi
 real:: temp,f
 integer:: i,j,l

 chi = 0.0
 chi(1) = 1.0
 do l = 2,Nexp
    if (g /= 0.0) then
         chi(l)=(2.0*(l-1)+1.0)*g**(l-1)
    endif
 enddo

  Pf = 0.0
  do i = 1,Nzen/2
     do j = 1,Nzen/2
        do l = m+1,Nexp
           temp=Ypleg(l,i)*Ypleg(l,j)
           Pf(i,j) = Pf(i,j) + temp*chi(l)
         enddo
      enddo
  enddo
  
  Pb = 0.0
  do i = 1,Nzen/2
     do j = Nzen/2 + 1, Nzen
        do l = m+1,Nexp
           temp=Ypleg(l,i)*Ypleg(l,j)
           Pb(i,j-Nzen/2) = Pb(i,j-Nzen/2) + temp*chi(l)
         enddo
      enddo
  enddo
end subroutine PHASE_MATRIX

!----------------------------------------------------------------------
!**********************************************************
!** gaussian quadrature from Numerical Recipes(gauleg.for)
!***
!** input x1, x2 - the lower and upper bounds of integration
!**       n - the number of quadrature points desired
!** output: x - quadrature abscissa
!**         w - quadrature weights
!**
!** changes from original:
!** 1) made f90 compatible
!** 2) made x,w real (kind=real4)
!** 3) switched order so the x(1) is largest value
!**********************************************************
subroutine GAULEG(x1,x2,x,w,n)
  real, intent(in)::  x1,x2
  integer, intent(in):: n
  real, dimension(:), intent(out):: x,w
  integer:: i,j,m
  double precision:: p1,p2,p3,pp,xl,xm,z,z1
  double precision, parameter:: EPS = 3.0d-14
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do  i=1,m
        z=cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
        do
          p1=1.0d0
          p2=0.0d0
          do j=1,n
            p3=p2
            p2=p1
            p1=((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
          enddo
          pp=n*(z*p1-p2)/(z*z-1.0d0)
          z1=z
          z=z1-p1/pp
         if (abs(z-z1) > EPS) then
!           print *,"cycling"
          cycle
         else
          exit
         end if
        end do
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.0d0*xl/((1.0d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
    enddo

  end subroutine GAULEG

!----------------------------------------------------------------------C
!      Computation of renormalized associated Legendre polynomials     C  
!      of the form:                                                    C
!                                                                      C
!                m                     1/2 m                           C
!               Y (u) = [(l-m)!/(l+m)!]   P (u)                        C
!                l                         l                           C
!                                                                      C
!         where,                                                       C
!                 m                                                    C
!                P (u) = Associated Legendre polynomial of             C
!                 l      order l and degree m                          C
!                                                                      C
!                    u = Cosine of zenith angle                        C
!                                                                      C
!     Reference:                                                       C
!                                                                      C
!             Dave, J. V.,and Armstrong, B. H., 1970: Computations     C
!                  of High-order Associated Legendre Polynomials,      C
!                  JQSRT, 10, 557-562, 1970                            C
!----------------------------------------------------------------------C
!                 I N P U T    V A R I A B L E S :                     C
!----------------------------------------------------------------------C
!         ANGLEO   :    Cosine of solar zenith angle                   C
!         GMU      :    Gaussian quadrature points                     C
!         M        :    Index for degree of Legendre polynomial        C
!         NEXP     :    Tot. no. of polynomial expansion terms         C
!         NZEN     :    Tot. no. of quadrature points                  C
!----------------------------------------------------------------------C
!                O U T P U T    V A R I A B L E S :                    C
!----------------------------------------------------------------------C
!         YPLEG    :    Renormalized associated Legendre polynomials   C
!----------------------------------------------------------------------C
!     REAL  GMU(NZEN),LM1,LM2,MU,YPLEG(NEXP,NZEN+2)
  subroutine PLEG(ANGLEO,GMU,M,NEXP,NZEN,YPLEG)
      integer, intent(in):: NEXP,NZEN,M
      real,intent(in),dimension(:):: GMU
      real,intent(in):: ANGLEO
      real, intent(out),dimension(:,:)::YPLEG
      real:: LM1,LM2,MU
      real, save:: CM1,DM1
      real:: CM2,DM2
      integer:: IL,L,I

      if (M == 0) then
        CM1 = 1.0
        DM1 = 1.0
      else
        CM2 = -sqrt(real( (2.0*M-1.0)/(2.0*M)))*CM1
        DM2 = -sqrt(real( (2.0*M+1.0)/(2.0*M)))*DM1
      end if
      do I=1,NZEN+2
        if (I <= NZEN) then
           MU = GMU(I)
        endif
        if (I == NZEN+1) then
           MU = ANGLEO
        endif
        if (I == NZEN+2) then
           MU = 0.0
        endif
!----------------------------------------------------------------------C
!                  Compute initial values                              C
!----------------------------------------------------------------------C
        if (M == 0) then
          YPLEG(1,I) = 1.0
          if (NEXP >= 2) then
            YPLEG(2,I) = MU       
          endif
        else  
          if(M+1 <= NEXP) then
            YPLEG(M+1,I) = CM2*(1-MU**2)**(M/2.0)
          endif
          if(M+2 <= NEXP) then
            YPLEG(M+2,I) = DM2*MU*(1-MU**2)**(M/2.0)
          endif
        end if
        do IL=1,NEXP-2 
          L = IL-1
          if (IL < M+1) then
            YPLEG(IL,I) = 0.0
            cycle
          else
!----------------------------------------------------------------------C
!          Use varying-degree recurrence formula to compute            C
!          successive values                                           C
!----------------------------------------------------------------------C
            LM2 = (L-M+2)*(L+M+2)
            LM1 = (L+M+1)*(L-M+1)
            YPLEG(IL+2,I) = ((2*L+3)*MU*YPLEG(IL+1,I)- &
                            sqrt(real(LM1))* &
                             YPLEG(IL,I))/sqrt(real(LM2))
          end if
      end do  
    end do 
      if (M > 0)  then
        CM1 = CM2
        DM1 = DM2
      end if 
    end subroutine PLEG

!**********************************************************************
!**** compute matrix inverse of a square matrix
!**********************************************************************
  subroutine INVERT(A,n)

!     MATRIX INVERSION USING GAUSS-JORDAN REDUCTION WITH
!     PARTIAL PIVOTING.
!
!       SEARCH FOR LARGEST PIVOT ELEMENT
!
      integer, intent(in):: n
      real, dimension (:,:), intent(inout)::  A
      integer:: K,JJ,I,KP1,J,L,KROW,IROW
      real:: AB,BIG,STORI,STORJJ,STORK
      integer, dimension(n,2):: INTER
      fourteen: do K=1,n
        JJ = K
        if (K /= n) then
          KP1 = K+1
          BIG = abs(A(K,K))
          five: do I=KP1,n
            AB = abs(A(I,K))
            if (BIG >= AB) then
              cycle
            endif
            BIG = AB
            JJ = I
          end do five
        endif
!
!       STORE NUMBERS OF ROWS INTERCHANGED. if JJ=K, THERE IS
!       NO INTERCHANGE.
!
        INTER(K,1) = K        !6
        INTER(K,2) = JJ
        if (JJ /= K)  then
!
!       ROW INTERCHANGE
!
        eight: do J=1,n
          STORK = A(K,J)
          STORJJ = A(JJ,J)
          A(K,J) = STORJJ
          A(JJ,J) = STORK
        end do eight

        endif
!
!       CALCULATE NEW ELEMENTS OF PIVOT ROW EXCEPT FOR PIVOT
!       ELEMENT.
!
        ten: do J=1,n
          if (J == K) then
           cycle
          endif
          A(K,J) = A(K,J)/A(K,K)
        end do ten
!
!       CALCULATE NEW ELEMENT REPLACING PIVOT TERM.
!
        A(K,K) = 1.0/A(K,K)
!
!       CALCULATE NEW ELEMENTS NOT IN PIVOT ROW OR COLUMN.
!
        twelve: do I=1,n
          if (I == K) then
             cycle
          endif
          eleven: do J=1,n
             if (J == K) then
               cycle
             endif
             A(I,J) = A(I,J)-A(K,J)*A(I,K)
          end do eleven
        end do twelve

!
!       CALCULATE NEW ELEMENTS FOR PIVOT COLUMN, EXCEPT FOR
!       PIVOT ELEMENT.
!
        thirteen: do I=1,n
          if (I == K)  then
            cycle
          endif
          A(I,K) = -A(I,K)*A(K,K)
        end do thirteen

    end do fourteen
!
!       REARRANGE COLUMNS OF FINAL MATRIX.
!
      sixteen: do L=1,n
        K = n-L+1
        KROW = INTER(K,1)
        IROW = INTER(K,2)
        if (KROW == IROW) then
           cycle
        end if
        fifteen: do I=1,n
          STORI = A(I,IROW)
          STORK = A(I,KROW)
          A(I,IROW) = STORK
          A(I,KROW) = STORI
        end do fifteen
      end do sixteen
   end subroutine INVERT
end module SOI
