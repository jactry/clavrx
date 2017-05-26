!----------------------------------------------------------------------
! Third Include File
! Clavrx version
!----------------------------------------------------------------------

   Sfc_Level_RTM = Rtm(Inwp,Jnwp)%Sfc_Level
   !Inver_Level_RTM = Rtm(Inwp,Jnwp)%Inversion_Level !WCS3
   Tropo_Level_RTM = Rtm(Inwp,Jnwp)%Tropo_Level
   Surface_Temperature = Tsfc_Nwp_Pix(Elem_Idx,Line_Idx)
   Tropopause_Temperature = Ttropo_Nwp_Pix(Elem_Idx,Line_Idx)

   !--- populate height and temperature profiles
   if (Inwp <= 0 .or. Jnwp <= 0) then
     print *, "bad nwp indices in awg"
   endif
   if (Allocated(Rtm(Inwp,Jnwp)%T_Prof) .eqv. .false.) then
      print *, "error, T_prof not allocated"
      print *, Elem_Idx,Line_Idx
      print *, Bad_Pixel_Mask(Elem_Idx,Line_Idx)
      print *, Satzen(Elem_Idx,Line_Idx)
      print *, Bt_Ch31(Elem_Idx,Line_Idx)
   endif
   Temp_Prof_RTM = Rtm(Inwp,Jnwp)%T_prof
   Hght_Prof_RTM = Rtm(Inwp,Jnwp)%Z_prof
   Press_Prof_RTM = P_Std_Rtm

   !--- populate radiance and transmission profiles
   if (Chan_On(Chan_Idx_67um) == sym%YES) then
     Atm_Rad_Prof_67um_RTM = Rtm(Inwp,Jnwp)%ch(27)%d(Ivza)%Rad_Atm_Profile
     Atm_Trans_Prof_67um_RTM = Rtm(Inwp,Jnwp)%d(Ivza)%ch(27)%Trans_Atm_Profile
   endif
   if (Chan_On(Chan_Idx_73um) == sym%YES) then
     Atm_Rad_Prof_73um_RTM = Rtm(Inwp,Jnwp)%ch(28)%d(Ivza)%Rad_Atm_Profile
     Atm_Trans_Prof_73um_RTM = Rtm(Inwp,Jnwp)%ch(28)%d(Ivza)%Trans_Atm_Profile
   endif
   if (Chan_On(Chan_Idx_85um) == sym%YES) then
     Atm_Rad_Prof_85um_RTM = Rtm(Inwp,Jnwp)%d(Ivza)%ch(29)%Rad_Atm_Profile
     Atm_Trans_Prof_85um_RTM = Rtm(Inwp,Jnwp)%d(Ivza)%ch(29)%Trans_Atm_Profile
   endif
   if (Chan_On(Chan_Idx_11um) == sym%YES) then
      Atm_Rad_Prof_11um_RTM = Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Rad_Atm_Profile
      Atm_Trans_Prof_11um_RTM = Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Trans_Atm_Profile
      Black_Body_Rad_Prof_11um_RTM = Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Rad_BB_Cloud_Profile
   endif
   if (Chan_On(Chan_Idx_12um) == sym%YES) then
      Atm_Rad_Prof_12um_RTM = Rtm(Inwp,Jnwp)%d(Ivza)%ch(32)%Rad_Atm_Profile
      Atm_Trans_Prof_12um_RTM = Rtm(Inwp,Jnwp)%d(Ivza)%ch(32)%Trans_Atm_Profile
   endif
   if (Chan_On(Chan_Idx_133um) == sym%YES) then
      Atm_Rad_Prof_133um_RTM = Rtm(Inwp,Jnwp)%d(Ivza)%ch(33)%Rad_Atm_Profile
      Atm_Trans_Prof_133um_RTM = Rtm(Inwp,Jnwp)%d(Ivza)%ch(33)%Trans_Atm_Profile
   endif
