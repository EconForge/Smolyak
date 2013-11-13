load('num_stab_approxMATLAB.mat', 'X')
load('num_stab_approxMATLAB.mat', 'Y')

b_ols_no = Num_Stab_Approx(X, Y, 1, 0, 0)
b_ols_yes = Num_Stab_Approx(X, Y, 1, 0, 1)
b_svd_yes = Num_Stab_Approx(X, Y, 2, 0, 1)
b_svd_no = Num_Stab_Approx(X, Y, 2, 0, 0)
b_ladpp_no = Num_Stab_Approx(X, Y, 3, 0, 0)
b_ladpp_yes = Num_Stab_Approx(X, Y, 3, 0, 1)
b_laddp_yes = Num_Stab_Approx(X, Y, 4, 0, 1)
b_laddp_no = Num_Stab_Approx(X, Y, 4, 0, 0)
penalty = 7
penalty_tsvd = -7
b_rlsTik_yes = Num_Stab_Approx(X, Y, 5, penalty, 1)
b_rlsSVD_yes = Num_Stab_Approx(X, Y, 6, penalty, 1)
b_rladPP_yes = Num_Stab_Approx(X, Y, 7, penalty, 1)
b_rladDP_yes = Num_Stab_Approx(X, Y, 8, penalty, 1)

