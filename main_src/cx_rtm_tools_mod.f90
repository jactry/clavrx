!  $Id:$
!
!

module cx_rtm_tools_mod
   use cx_rtm_constants_mod, only : &
           nLevels_RTM
           
   integer , save :: k_Rtm_nwp ( nlevels_rtm)
   real , save :: x_rtm_nwp ( nlevels_rtm)
   
   
   logical , save :: is_mapped = .false.
   
   
   
contains
   !
   !  called by rtm_mod
   !
   subroutine convert_profiles_nwp_rtm  ( &
      &    p_std_nwp &
      &  , sfc_level_nwp &
      &  , p_sfc_nwp &
      &  , wvmr_sfc_nwp &
      &  , tmpair_nwp &
      &  , t_nwp &
      &  , ozmr_nwp &
      &  , wvmr_nwp &
      &  , z_nwp & 
      &  , z_rtm &
      &  , t_prof_rtm &
      &  , wvmr_prof_rtm &
      &  , ozmr_prof_rtm &
       )
      
      use cx_rtm_constants_mod, only : &
           p_std_rtm &
         , t_std_rtm &
         , wvmr_std_rtm &
         , ozmr_std_rtm &
         , nLevels_RTM
         
      implicit none  
      
      real, dimension(:), intent(in) :: p_std_nwp  
      integer , intent(in) :: sfc_level_nwp
      real,  intent(in) :: p_sfc_nwp
      real, intent(in) :: tmpair_nwp
      real, intent(in) :: wvmr_sfc_nwp
      real, dimension(:), intent(in) :: t_nwp
      real, dimension(:), intent(in) :: ozmr_nwp
      real, dimension(:), intent(in) :: wvmr_nwp
      real, dimension(:), intent(in) :: z_nwp
      
      
      integer :: k, Lowest_Level_Rtm_Nwp, Highest_Level_Rtm_Nwp,  &
              Sfc_Level_Rtm
      real:: dT_dP_near_Sfc
      real:: dWvmr_dP_near_Sfc
      real:: dZ_dP_near_Sfc
      real:: P_near_Sfc_Nwp
     
      real:: T_Offset
      
      integer :: next_level
      integer :: nLevels_nwp
      
      
      real, dimension ( NLevels_Rtm ) , intent(out)  ::  z_rtm
      real, dimension ( NLevels_Rtm ) , intent(out)  ::  t_prof_rtm
      real, dimension ( NLevels_Rtm ) , intent(out)  ::  wvmr_prof_rtm
      real, dimension ( NLevels_Rtm ) , intent(out)  ::  ozmr_prof_rtm
      
      !--- initialize indices
      Lowest_Level_Rtm_Nwp = NLevels_Rtm
      Highest_Level_Rtm_Nwp = 1
      Sfc_Level_Rtm = NLevels_Rtm
      P_Near_Sfc_Nwp = P_std_Nwp(Sfc_Level_Nwp)
      
      nLevels_nwp = size(p_std_nwp)
      
      !--- determine some critical Levels in the Rtm profile
      do k = 1, NLevels_Rtm
         if (P_Std_Rtm(k) > P_std_Nwp(1)) then
            Highest_Level_Rtm_Nwp = k
            exit
         end if
      end do

      do k = NLevels_Rtm,1,-1
         if (P_Std_Rtm(k) < P_sfc_nwp) then
            Sfc_Level_Rtm = k
            exit
         endif
      enddo

      do k = NLevels_Rtm,1,-1
         if (P_Std_Rtm(k) < P_near_Sfc_Nwp) then
            Lowest_Level_Rtm_Nwp = k
            exit
         end if
      end do
      
      
      !--- compute T and Wvmr lapse rate near the surface
      dT_dP_near_Sfc = 0.0
      dZ_dP_near_Sfc = 0.0
      dWvmr_dP_near_Sfc = 0.0
      
      next_level = 0
      if (P_sfc_Nwp /= P_Std_Nwp(NLevels_Nwp)) next_level = - 1
      
      dT_dP_near_Sfc =  &
             (T_Nwp(Sfc_Level_Nwp + next_level) - Tmpair_Nwp)/ &
             (P_Std_Nwp(Sfc_Level_Nwp + next_level ) - P_sfc_Nwp)
             
      dWvmr_dP_near_Sfc =  &
             (Wvmr_Nwp(Sfc_Level_Nwp + next_level ) - Wvmr_Sfc_nwp)/ &
             (P_Std_Nwp(Sfc_Level_Nwp + next_level ) - P_sfc_Nwp )
             
      dZ_dP_near_Sfc =  &
             (Z_Nwp(Sfc_Level_Nwp + next_level) - 0.0)/ &
             (P_Std_Nwp(Sfc_Level_Nwp + next_level ) - P_sfc_Nwp)
             
             
             
      !  compute temperature offset between standard and NWP profiles at top
      !--- this will be added to the standard profiles        
      
      T_Offset = T_Nwp(1) - T_Std_Rtm(Highest_Level_Rtm_Nwp)
      
      
       !--- for Rtm Levels above the highest nwp Level, use Rtm standard values
      do k = 1,Highest_Level_Rtm_Nwp-1
         Z_Rtm(k) = Z_Nwp(1)
         T_Prof_Rtm(k) = T_Std_Rtm(k) + T_Offset
         Wvmr_Prof_Rtm(k) = Wvmr_Std_Rtm(k)
         Ozmr_Prof_Rtm(k) = Ozmr_Std_Rtm(k)
      end do
      
      if ( .not. is_mapped ) then
         call map_nwp_rtm ( p_std_nwp )
      
      end if 

      !--- Rtm Levels within standard nwp Levels above the surface
      do k = Highest_Level_Rtm_Nwp, Lowest_Level_Rtm_Nwp
         T_Prof_Rtm(k) = T_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (T_Nwp(k_Rtm_Nwp(k)+1) - T_Nwp(k_Rtm_Nwp(k)))
         Z_Rtm(k) = Z_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (Z_Nwp(k_Rtm_Nwp(k)+1) - Z_Nwp(k_Rtm_Nwp(k)))
         Wvmr_Prof_Rtm(k) = Wvmr_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (Wvmr_Nwp(k_Rtm_Nwp(k)+1) - Wvmr_Nwp(k_Rtm_Nwp(k)))
         Ozmr_Prof_Rtm(k) = Ozmr_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (Ozmr_Nwp(k_Rtm_Nwp(k)+1) - Ozmr_Nwp(k_Rtm_Nwp(k)))
      end do

      !--- Rtm Levels that are below the lowest nwp Level but above the surface
      do k = Lowest_Level_Rtm_Nwp+1, Sfc_Level_Rtm
           T_Prof_Rtm(k) = Tmpair_Nwp + dT_dP_near_Sfc * &
                      (P_Std_Rtm(k) - P_sfc_Nwp)
           Wvmr_Prof_Rtm(k) = Wvmr_Sfc_nwp + dWvmr_dP_near_Sfc * &
                      (P_Std_Rtm(k) - P_sfc_Nwp)
           Z_Rtm(k) = dZ_dP_near_Sfc * &
                      (P_Std_Rtm(k) - P_sfc_Nwp)
           Ozmr_Prof_Rtm(k) = Ozmr_Nwp(NLevels_Nwp)
      end do

      !--- Rtm Levels below the surface
      do k = Sfc_Level_Rtm +1, NLevels_Rtm
         T_Prof_Rtm(k) = Tmpair_Nwp
         Z_Rtm(k) = dZ_dP_near_Sfc * &
                   (P_Std_Rtm(k) - P_sfc_Nwp)
         Wvmr_Prof_Rtm(k) = Wvmr_Sfc_nwp
         Ozmr_Prof_Rtm(k) = Ozmr_Std_Rtm(k)
      end do

      
       
   end subroutine convert_profiles_nwp_rtm
   
   !  -------
   !
   ! ----------
   subroutine map_nwp_rtm (p_std_nwp)
      use cx_tools_science_mod, only : locate
       
      use cx_rtm_constants_mod, only : &
          p_std_rtm 
         
      real, dimension(:), intent(in) :: p_std_nwp  
      integer ::  k, k_temp
      integer :: nLevels_nwp
      integer :: k_nwp_rtm ( nlevels_rtm)
      real :: x_nwp_rtm ( nlevels_rtm)
      
      nLevels_nwp = size(p_std_nwp)
      
      !--- map NWP Levels to RTM Levels
      do k = 1 , nlevels_nwp
         
         call locate(P_Std_Rtm, NLevels_Rtm, P_Std_Nwp(k), k_temp)
         k_Nwp_Rtm(k) = min(NLevels_Rtm-1,max(1,k_temp))

         !-- compute linear weighting factor
         x_Nwp_Rtm(k) = (P_Std_Nwp(k) - P_Std_Rtm(k_Nwp_Rtm(k))) / &
                      (P_Std_Rtm(k_Nwp_Rtm(k)+1) - P_Std_Rtm(k_Nwp_Rtm(k))) 

      end do

      !--- map RTM Levels to NWP Levels
      do k = 1, NLevels_Rtm
      
         !--- locate the Rtm Level within the nwp Levels 
         !--- P_Std_Rtm(k) should fall between nwp Levels k_temp and k_temp +1
         call locate(P_Std_Nwp, NLevels_Nwp, P_Std_Rtm(k), k_temp)
         k_Rtm_Nwp(k) = min(NLevels_Nwp-1,max(1,k_temp))

         !-- compute linear weighting factor
         x_Rtm_Nwp(k) = (P_Std_Rtm(k) - P_Std_Nwp(k_Rtm_Nwp(k))) / &
                      (P_Std_Nwp(k_Rtm_Nwp(k)+1) - P_Std_Nwp(k_Rtm_Nwp(k)))
     
     
      end do
      
      is_mapped = .true.
    
   
   
   end subroutine map_nwp_rtm
   
   
 
   !====================================================================
   ! Function Name: emissivity
   !
   ! Function:
   !  Computes the  effective emissivity
   !
   ! Input:  Rad_Toa - channel radiance at top of atmosphere(toa)
   !         Rad_Clear_Tau - channel radiance at toa for clear skies
   !         Rad_Cloud_BB_Toa - channel radiance at TOA if cloud were a Black-Body
   !
   ! Output: Emiss - the effective cloud emissivity
   !
   !====================================================================
   function EMISSIVITY(Rad_Toa, Rad_Clear_Toa, Rad_Cloud_BB_Toa) result(Emiss)
      real, intent(in) :: Rad_Toa
      real, intent(in) :: Rad_Clear_Toa
      real, intent(in) :: Rad_Cloud_BB_Toa
      real :: Emiss

      Emiss = -999.

      if (Rad_Cloud_BB_Toa /= Rad_Clear_Toa .and. rad_toa > 0 .and. rad_clear_toa > 0 ) then
         Emiss = (Rad_Toa - Rad_Clear_Toa) / (Rad_Cloud_BB_Toa - Rad_Clear_Toa) 
      end if

      return

   end function EMISSIVITY

end module cx_rtm_tools_mod
