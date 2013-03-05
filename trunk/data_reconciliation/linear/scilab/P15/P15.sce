// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Steam metering system
//Serth, R W, and W A Heenan. 1986.
// Gross error detection and data reconciliation in steam-metering systems. 
//AIChE Journal 32: 733-747.
//Bibtex Citation

//@article{Serth1986,
//author = {Serth, R W and Heenan, W A},
//journal = {AIChE Journal},
//pages = {733--747},
//title = {{Gross error detection and data reconciliation in steam-metering systems}},
//volume = {32},
//year = {1986}
//}

// 28 Streams
// 11 Equipments 
// the measures
clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

xm =[0.875
0.989
108.426
107.799
50.149
110.327
2.201
165.779
0.821
51.779
14.757
67.183
107.734
89.855
58.554
23.258
31.879
15.843
7.920
10.130
91.047
5.555
2.471
44.687
81.656
82.839
70.685
72.913
];
// in original paper the standard deviation is given. so it must be squared.
var = [0.022
0.025
2.796
2.749
1.332
2.807
0.058
4.101
0.021
1.310
0.372
1.682
2.782
2.297
1.500
0.591
0.818
0.406
0.196
0.263
2.183
0.136
0.065
1.166
2.137
2.033
1.770
1.806].^2;
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    25    26    27    28
jac = [ 1   1   -1  1   0   0   0   0    0   0   0   0   0     0     0      0    0      0    0    0     0    0     0     0     0     0    0      0
        0   0   0   0  -1   -1  1   1    -1  0   0   0   0     0     0      0    0      0    0    0     0    0     0     0     0     0    0      0 
        -1  0   0   0   1    0   0   0   0  -1   0   0   0     0     0      0    0      0    0    0     0    0     0     0     0     0    0      0
        0   0   0   0   0   0   0   0    0   1   1   -1  0     0     0      0    0      0    0    0     0    0     0     0     0     0    0      0
        0   0   1   0   0   0   0   0    0   0   -1  0   1     -1    -1     -1   -1     0    0    0     0    0     0     0     0     0    0      0  
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    25    26    27    28
        0   -1  0   0   0   1   0   0    0   0   0   0   -1     0     0     0    0      0    0    0     0    0     0     0     0     0    0      0
        0   0   0   0   0   0   -1  0    0   0   0   0   0     1     0      0    0      1    -1   -1    -1   0     0     0     0     0    0      0
        0   0   0   0   0   0   0   0    0   0   0   0   0     0     1      0    0      -1   0    0     0    1     -1    -1    0     0    0      0
        0   0   0   0   0   0   0   0    0   0   0   1   0     0     0      1    0      0    0    0     0    -1    0     0     -1    0    0      0  
        0   0   0   0   0   0   0   0    0   0   0   0   0     0     0      0    0      0    1    0     0    0     1     0     0     -1    1     0
        0   0   0   0   0   0   0   -1   0   0   0   0   0     0     0      0    0      0    0    1     0    0     0     0     0     1     0     1
        ];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    25    26    27    28

//observability/redundancy tests                  
umeas_P15 = [];
[red_P15, just_measured_P15, observ_P15, non_obs_P15, spec_cand_P15] = qrlinclass(jac,umeas_P15)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_P15 = setdiff([1:length(xm)], umeas_P15);
red = measured_P15;//
// to reconcile with all variables, comment the line above and uncomment bellow
//red = [1:length(xm)];

// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS = 0
// Absolute sum of squares = 1
//Cauchy = 2
//Contamined Normal = 3
//Fair  = 4
//Hampel = 5
//Logistic = 6
//Lorenztian = 7
//Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = 0;
exec ../functions/setup_DR.sce
// to run robust reconciliation, it is also necessary to choose the function to return the problem structure
if obj_function_type > 0 then
[nc_eq, n_non_lin_eq, nv, nnzjac_ineq, nnzjac_eq, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = robust_structure(jac, 0, xm, objfun, res_eq, res_ineq);
else
// for WLS, only the line bellow must be choosen and comment the 3 lines above
[nc, nv, i1, i2, nnzeros, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = wls_structure(jac);
end


params = init_param();
// We use the given Hessian
params = add_param(params,"hessian_approximation","exact");
params = add_param(params,"derivative_test","second-order");
params = add_param(params,"tol",1e-8);
params = add_param(params,"acceptable_tol",1e-8);
params = add_param(params,"mu_strategy","adaptive");
params = add_param(params,"journal_level",5);

[x_sol, f_sol, extra] = ipopt(xm, objfun, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

mprintf("\n\nSolution: , x\n");
for i = 1 : nv
    mprintf("x[%d] = %e\n", i, x_sol(i));
end


mprintf("\n\nObjective value at optimal point\n");
mprintf("f(x*) = %e\n", f_sol);
// Balances mounted based on the visual diagram for final check
// equip 1
x_sol(1)+x_sol(2) + x_sol(4)-(x_sol(3))
// equip 2
x_sol(7)+x_sol(8)-(x_sol(5)+x_sol(6)+x_sol(9))
// equip 3
x_sol(5)-(x_sol(1)+x_sol(10))
// equip 4
x_sol(10)+x_sol(11)-(x_sol(12))
// equip 5
x_sol(3)+x_sol(13)-(x_sol(11)+x_sol(14)+x_sol(15)+x_sol(16)+x_sol(17))
// equip 6
x_sol(6)-x_sol(2)-(x_sol(13))
// equip 7
x_sol(14)+x_sol(18)-(x_sol(7)+x_sol(18)+x_sol(19)+x_sol(20)+x_sol(21))
// equip 8
x_sol(15)+x_sol(22)-(x_sol(18)+x_sol(23)+x_sol(24))
// equip 9
x_sol(16)+x_sol(12)-(x_sol(22)+x_sol(25))
// equip 10
x_sol(27)+x_sol(23)+x_sol(19)-(x_sol(26))
// equip 11
x_sol(20)+x_sol(26)+x_sol(28)-(x_sol(8))


