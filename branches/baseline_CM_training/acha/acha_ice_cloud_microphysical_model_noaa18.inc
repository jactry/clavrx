! Description: Include file that contains ice microphyscial model for MODIS sensor

select case(habit_idx)
   case(1)  !-- ice clouds modeled as droxtals (b=0.10)
        A_BETA_11um_67um_FIT_ICE =    1.31754
        B_BETA_11um_67um_FIT_ICE =   -0.28391
        A_BETA_11um_85um_FIT_ICE =    1.92649
        B_BETA_11um_85um_FIT_ICE =   -0.90985
        A_BETA_11um_133um_FIT_ICE =    0.01653
        B_BETA_11um_133um_FIT_ICE =    1.04672
        A_Re_Beta_FIT_ICE =  -51.86815
        B_Re_Beta_FIT_ICE =  125.68085
        C_Re_Beta_FIT_ICE = -100.35716
        D_Re_Beta_FIT_ICE =   26.94883
        A_Qe_065um_FIT_ICE =    2.37059
        B_Qe_065um_FIT_ICE =   -0.33533
        C_Qe_065um_FIT_ICE =    0.07915
        A_Qe_11um_FIT_ICE =   -0.11737
        B_Qe_11um_FIT_ICE =    2.49136
        C_Qe_11um_FIT_ICE =   -0.71160
        A_wo_11um_FIT_ICE =   -0.07234
        B_wo_11um_FIT_ICE =    0.60712
        C_wo_11um_FIT_ICE =   -0.15519
        A_g_11um_FIT_ICE =    0.32167
        B_g_11um_FIT_ICE =    0.78961
        C_g_11um_FIT_ICE =   -0.23527

   case(2)  !-- ice clouds modeled as solid_bullet_rosettes (b=0.10)


   case(3)  !-- ice clouds modeled as hollow_bullet_rosettes (b=0.10)
   case(4)  !-- ice clouds modeled as solid_columns (b=0.10)
   case(5)  !-- ice clouds modeled as hollow_columns (b=0.10)
   case(6)  !-- ice clouds modeled as plates (b=0.10)
   case(7)  !-- ice clouds modeled as aggregate_columns (b=0.10)
        A_BETA_11um_67um_FIT_ICE =    1.54172
        B_BETA_11um_67um_FIT_ICE =   -0.49635
        A_BETA_11um_85um_FIT_ICE =    2.17033
        B_BETA_11um_85um_FIT_ICE =   -1.14189
        A_BETA_11um_133um_FIT_ICE =    0.11629
        B_BETA_11um_133um_FIT_ICE =    0.93978
        A_Re_Beta_FIT_ICE =  -39.22696
        B_Re_Beta_FIT_ICE =   97.67353
        C_Re_Beta_FIT_ICE =  -81.09035
        D_Re_Beta_FIT_ICE =   23.08110
        A_Qe_065um_FIT_ICE =    2.24408
        B_Qe_065um_FIT_ICE =   -0.23005
        C_Qe_065um_FIT_ICE =    0.05641
        A_Qe_11um_FIT_ICE =   -0.02317
        B_Qe_11um_FIT_ICE =    2.12770
        C_Qe_11um_FIT_ICE =   -0.57337
        A_wo_11um_FIT_ICE =   -0.01382
        B_wo_11um_FIT_ICE =    0.57493
        C_wo_11um_FIT_ICE =   -0.15368
        A_g_11um_FIT_ICE =    0.59489
        B_g_11um_FIT_ICE =    0.46862
        C_g_11um_FIT_ICE =   -0.14131
   case(8)  !-- ice clouds modeled as small_aggregate_plates (b=0.10)
   case(9)  !-- ice clouds modeled as large_aggregate_plates (b=0.10)
   case default  !-- ice clouds modeled as empirical model (b=0.10)

end select
