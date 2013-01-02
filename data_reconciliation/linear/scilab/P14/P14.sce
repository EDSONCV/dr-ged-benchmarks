// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Proposed by author

// 24 Streams
// 14 Equipments 

clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

xm =[54.4
149
157.2
148.1
466
466
510.3
218
292.3
292.3
84.4
13.5
70.9
13.5
208
2.9
2.9
87.3
189.6
204
54.3
149.7
201.7
16.3
];
//the variance 
var = [0.284217
2.132184
2.373325
2.106504
57.932395
57.932395
69.470558
4.564205
8.205589
8.205589
0.684127
0.194481
0.596018
0.194481
4.155075
0.008974
0.008974
8.132763
3.452461
3.996801
0.283173
2.152265
3.907185
0.025517
];
//

//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    
jac = [ 1   1  1    1   -1  0   0   0    0   0   0   0   0     0     0      0    0      0    0    0     0    0     0     0    //M1
        0   0   0   0   1   -1  0   0    0   0   0   0   0     0     0      0    0      0    0    0     0    0     0     0    //F1
        0   0   0   0   0   1  -1   0    0   0   0   0   0     0     0      0    0      0    0    0     0    0     0     0    //T1
        0   0   0   0   0   0   1   -1   -1  0   0   0   0     0     0      0    0      0    0    0     0    0     0     0   //S1  
        0   0   0   0   0   0   0   0    1  -1   0   0   0     0     0      0    0      0    0    0     0    0     0     0     //F2
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    
        0   0   0   0   0   0   0   0    0   1   -1  0   0     0     -1     0    0      0    0    0     0    0     0     0    //M2
        0   0   0   0   0   0    0  0    0   0   1  -1  -1     0     0      0    0      0    0    0     0    0     0     0    //S3
        0   0   0   0   0   0   0   0    0   0   0   1   0     -1    0      0    0      0    0    0     0    0     0     0   //F3
        0   0   0   0   0   0   0   0    0   0   0   0   0     0     1      -1   0      0    -1   0     0    0     0     0    //F4
        0   0   0   0   0   0   0   0    0   0   0   0   0     0     0      1    -1     0    0    0     0    0     0     0    //F5
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    
        0   0   0   0   0   0   0    0   0   0   0   0   1     0     0      0    1     -1    0    0     0    0     0     0    //S5
        0   0   0   0   0   0   0    0   0   0   0   0   0     0     0      0    0      0    1    -1    0    0     0     0    //T2
        0   0   0   0   0   0   0    0   0   0   0   0   0     0     0      0    0      0    0    1     -1   -1    0     0    //S4
        0   0   0   0   0   0   0    1   0   0   0   0   0     0     0      0    0      0    0    0     0    0     -1    -1    //S2
        ];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    

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
// Global Mass balance
x_sol(1)+x_sol(2)+x_sol(3)+x_sol(4)-x_sol(14)-x_sol(18)-x_sol(21)-x_sol(22)-x_sol(23)-x_sol(24)
//M1
x_sol(1)+x_sol(2)+x_sol(3)+x_sol(4)-x_sol(5)
//F1
x_sol(5)-x_sol(6)
//T1
x_sol(6)-x_sol(7)
//S1  
x_sol(7)-x_sol(8)-x_sol(9)
//F2
x_sol(9)-x_sol(10)
//M2
x_sol(10)-x_sol(11)-x_sol(15)
//S3
x_sol(11)-x_sol(12)-x_sol(13)
//F3
x_sol(12)-x_sol(14)
//F4
x_sol(15)-x_sol(16)-x_sol(19)
//F5
x_sol(16)-x_sol(17)
//S5
x_sol(13)+x_sol(17)-x_sol(18)
//T2
x_sol(19)-x_sol(20)
//S4
x_sol(20)-x_sol(21)-x_sol(22)
//S2
x_sol(6)-x_sol(7)
