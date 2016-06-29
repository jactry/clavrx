! $Header: https://svn.ssec.wisc.edu/repos/cloud_team_nlcomp/trunk/nlcomp_array.f90 8 2014-01-31 08:14:58Z awalther $
!
!  HISTORY: 2014/01/12
!         : AW first verisob of NLCOMP for arrays
! 
!
!
subroutine nlcomp_array_loop_sub ( input , output, debug_mode_user )
   
   use dncomp_interface_def_mod, only: &
      dncomp_in_type &
      , dncomp_out_type
   

   
   type (dncomp_in_type) , intent(in) :: input
   type (dncomp_out_type), intent(out) :: output
   integer , intent(in) , optional :: debug_mode_user
   
   
   
   

end subroutine nlcomp_array_loop_sub
