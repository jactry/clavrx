!$Id:$

!>
!!
!!

module cx_pfaast_mod
   use cx_pfaast_constants_mod, only: &
      pstd, tstd, wstd, ostd, nl
   use cx_pfaast_coef_mod   , only: &
      coef_type, NXD, NXW, NXC, NXO
   use cx_pfaast_tools_mod,only: conpir , calpir, taudoc &
      , tauwtr
   
   type ( coef_type ) :: coef
   
   type profile_type
      real, allocatable :: temp(:)
      real, allocatable :: press(:)
      real, allocatable :: ozon(:)
      real, allocatable :: wvp(:)
   end  type profile_type
   
contains   
   !>
   !!  @param temp Temperature profile
   !!  @kban_native channel number native 
   !!  @todo solve the channel index problem
   !!    each file seems to have own definitions
   !!    encourage Pfaast people to amke better files ..
   !!
   subroutine compute_transmission_pfaast ( &
       ancil_data_path &
       & ,temp &
       & ,wvmr &
       & ,ozmr & 
       & ,theta  &
       & ,rco2 &
       & ,sensor &
       & ,kban_native &
       & ,taut &
       & , modis_chn_eqv )
       
       
       
      implicit none    
     
      character(len = * ) , intent(in) :: ancil_data_path  
      real, intent(in)  :: temp (:)
      real, intent(in)  :: wvmr(:)
      real, intent(in)  :: ozmr(:)
      real, intent(in)  :: theta 
      real, intent(in)  :: rco2 
      character (len =* ), intent(in) :: sensor
      integer, intent(in)  :: kban_native 
      integer, intent(in), optional :: modis_chn_eqv
      
       
      integer :: kban
      real, intent(out)  :: taut (:)  
      integer,parameter :: NM = NL - 1 
      real, save:: pavg(nm),tref(nm),wref(nm),oref(nm)
      logical :: are_ref_profiles_set = .false.
      real :: oamt ( NM )
      real :: wamt ( NM )
      real :: tavg ( NM ) 
      real :: secz ( NM )
      real :: xdry(nxd,nm),xozo(nxo,nm),xwet(nxw,nm),xcon(nxc,nm)
      real, dimension(NL)  ::taud,tauw,tauo,tauc
      
      
      
      taud = 1.
      tauw = 1.
      tauc = 1.
      taut = 1.
      tauo = 1.
      
      ! - start
     
      ! - read coef data if needed      
      call coef % read_it ( trim(sensor) , ancil_data_path ) 
      
      if ( present ( modis_chn_eqv )) then
         !  find kban in coef data 
         if ( .not. any ( modis_chn_eqv ==  coef % modis_channel_eqv)) then
           ! print*,'wrong MODIS channel set for pfaast',modis_chn_eqv
           ! print*,'please use modis equivalent sensor channel number'
           !print*,'sensor ', sensor
           ! print*,'possible channels: ',coef % modis_channel_eqv
           ! print*,'returning...'
            return
         end if
         kban = minloc ( abs ( coef % modis_channel_eqv - modis_chn_eqv ), 1)

      else
      
      
         !  find kban in coef data 
         if ( .not. any ( kban_native ==  coef % native_channel)) then
            print*,'wrong channel set for pfaast',kban_native
            print*,'please use native sensor channel number'
            print*,'sensor ', sensor
            print*,'possible channels: ',coef % native_channel
            print*,'returning...'
            return
         end if
         kban = minloc ( abs ( coef % native_channel - kban_native ), 1)
      end if
      
      ! - computes mid-layer values for reference
      if ( .not. are_ref_profiles_set ) then
         call conpir ( pstd,tstd,wstd,ostd,nl,1,pavg,tref,wref,oref )
         are_ref_profiles_set = .true.
         print*,'PFAAST: compute.. reference profiles'
      end if   
      ! - compute mid-layer for input
      call conpir(pstd,temp,wvmr,ozmr,nl,1,pavg,tavg,wamt,oamt)
            
      ! - get predictors for arch comp group
      secz(:)=( 1./cos(0.01745329 * theta ))
      call calpir( &
            tref, wref, oref &
            ,tavg,wamt,oamt,pavg &
            ,secz &
            , nm,nxd,nxw,nxo,nxc &
            ,xdry,xwet,xozo,xcon)
      ! - matrix multiplication coef x predictors
      
         ! - dry
      call taudoc(coef % dry (:,:,kban),xdry,taud)  
      
         ! - ozone
      call taudoc(coef % ozon (:,:,kban),xozo,tauo)   
         ! - wvp cont
      call taudoc ( coef % wvp_cont (:,:,kban),xcon,tauc)   
         ! - wavp solid+liquid
      call tauwtr(coef % wvp_solid(:,:,kban), &
               coef % wvp_liquid (:,:,kban),xwet,tauw) 
             
        
      ! - build total transmisison
      taut = taud * tauo * tauc * tauw
      ! - done    
       
   
   end subroutine compute_transmission_pfaast

  



end module cx_pfaast_mod
