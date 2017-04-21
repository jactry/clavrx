! $Id$

module muri_retrieval_mod
   use muri_definitions_mod, only: &
      muri_input_type &
      , muri_output_type
   
   use muri_lut_mod, only: &
      lut
   
   use lib_array
      
   implicit none
   private
   public :: muri_algorithm 
contains

   !
   !
   !
   subroutine muri_algorithm (inp, out)
      type( muri_input_type), intent(in) :: inp
      type( muri_output_type), intent(out) :: out
      
     
      
      integer :: i_fm, i_cm
      real :: fmr
      real :: refl_toa (8,11,6)
      real :: aot_mixed  (8,11,6)
      integer :: i_cha,i_opt, i_fmr
      real :: refl_ch4 (8)
      real :: aot_temp(11)
      real :: refl_corrsp(11,6)
      real :: temp(8)
      real :: aot_allbands(11,6)
      real :: n1 (11)
      real :: n2
      real :: err(11,6)
      real :: err_sqrt (11)
      real :: aod_nleta_temp(4,5)
      real :: err_nleta_temp(4,5)
      real :: fmf_nleta_temp(4,5)
      integer :: idx(1), idx2(2)
      real :: val,val2
      real :: aod_allbands(4,5,6)
      integer:: k,i,j
      !print*,'start muri algorithm'
      ! call inp%info
      
      err_nleta_temp = huge(err_nleta_temp)
      err_sqrt = huge(err_sqrt)
     
      call lut % read_lut
      
      call lut % sub_table (inp % sol,inp%sat,inp%azi,inp%ws)
      
      
      ! simple test
   !   call lut % sub_table (22.,13.,120.,2.)
      
   !   do k = 1,4
    !     print*
   !      print*,'fine ',k
   !      do i= 1,8
   !         write(*,'(f7.5)') lut % refl_fine(i,:,k)
           
  !       end do
 !     end do
      
  !    do k = 1,5
  !    print*
   !      print*,'coarse ',k
   !      do i= 1,8
   !         print*,lut % refl_coarse(i,:,k)
    !     end do
   !   end do
      
      
      
      ! - loop over fine and coarse mode
      do i_fm = 1, 4
         do i_cm = 1, 5
            
            ! loop over fine mode ratio
            do i_fmr = 1,11
               fmr = (i_fmr-1)/10.
               
               do i_opt = 1, 8
                  do i_cha = 1,6
                     refl_toa (i_opt,i_fmr,i_cha) = fmr * lut % refl_fine(i_opt,i_cha,i_fm) &
                        & + (1 - fmr) * lut % refl_coarse(i_opt,i_cha,i_cm)
                     aot_mixed(i_opt,i_fmr,i_cha) = fmr * lut % aot_aer_fine(i_opt,i_cha,i_fm) &
                        & + (1 - fmr) * lut % aot_aer_coarse(i_opt,i_cha,i_cm)   
                     
                  end do
               end do  
               
               !- reference channel is #4
               refl_ch4 = refl_toa(:,i_fmr,4)
               
               
               aot_temp(i_fmr) = interp1d(refl_ch4,lut % aot_550nm, inp % rfl(4),bounds_error = .false., FILL_VALUE = -999.)
               
               
               do i_cha = 1,6
                 temp =  aot_mixed(:,i_fmr,i_cha)
                 aot_allbands (i_fmr,i_cha) = interp1d(refl_ch4,temp,inp % rfl(4),bounds_error = .false., FILL_VALUE = -999.)
               end do
               
                  
            end do
            
             do i_cha = 1,6
                  do i_fmr =1,11
                     refl_corrsp(i_fmr,i_cha) = interp1d(lut % aot_550nm,refl_toa(:,i_fmr,i_cha), aot_temp(i_fmr),bounds_error = .false., FILL_VALUE = -999. )
                     
                  end do
               end do
               
            
            do i_cha = 1,6
               !print*,i_cha,refl_corrsp(:,i_cha)
               n1 = inp % rfl(i_cha) - refl_corrsp(:,i_cha)
               n2 = inp % rfl(i_cha) - refl_toa(1,1,i_cha) + 0.001
               
               err(:,i_cha) = n1/n2
               
            end do
           
            err_sqrt = sqrt ( (err(:,3)**2. + err(:,4)**2.+ err(:,5)**2.+ err(:,6)**2.)/4.)  
            
            
            val= minval(err_sqrt)
            idx=minloc(err_sqrt)
            aod_nleta_temp(i_fm,i_cm) = aot_temp(idx(1))
            err_nleta_temp(i_fm,i_cm) = val
            fmf_nleta_temp(i_fm,i_cm) = (idx(1)) /10.
            
            aod_allbands (i_fm,i_cm,:) = aot_allbands(idx(1),:)
            
            !print*,i_fm,i_cm,    aod_nleta_temp(i_fm,i_cm) ,val,fmf_nleta_temp(i_fm,i_cm)     
            
         end do
      end do
      
      
      val2 = minval(err_nleta_temp)
      idx2 = minloc(err_nleta_temp)
      
      
      !- once find a good match between fwd and measurment give output
      out% aot = aod_nleta_temp(idx2(1),idx2(2))
      out % cm_mode = idx2(1)
      out % fm_mode =  idx2(2)
      out % fmf = fmf_nleta_temp(idx2(1),idx2(2))
      do i_cha =1,6 
         out% aot_channel(i_cha) = aod_allbands(idx2(1),idx2(2),i_cha)
      end do
      
   
   
   end subroutine  muri_algorithm  



end module muri_retrieval_mod
