// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Proposed by author
//10 Streams 
//6 Equipments

clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

// the measures
xm =[49.941
74.890
75.049
47.974
30.016
25.012
5.001
4.950
2.953
1.981
];
//the variance
var=[1
1
1
1
1
1
0.15
0.15
0.1
0.1
];

// gross error
gerror = zeros(length(xm),1);
// to setup gross errors, select the stream and magnitude as the line bellow
//gerror(2) = 9*sqrt(var(2));
xm = xm + gerror;


//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10    
jac = [ 1   -1   0  0   0   1   0   0    0   0      
        0   1    -1  0   0   0   0  0    0   0      
        0   0    1  -1   -1   0   0  0   1   0     
        0   0    0  0   1   -1   -1  0    0   0
        0   0    0  0   0   0   1   -1    0   0
        0   0    0  0   0   0   0   1    -1   -1       ];                                
//      1   2   3   4   5   6   7   8    9   10   

//observability/redundancy tests                  
umeas_P7 = [];
[red_P7, just_measured_P7, observ_P7, non_obs_P7, spec_cand_P7] = qrlinclass(jac,umeas_P7)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_P7 = setdiff([1:length(xm)], umeas_P7);
red = measured_P7;//
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


// global balance
x_sol(1)-x_sol(4) - x_sol(10)
//equipment 1
x_sol(1)+x_sol(6) - x_sol(2)
//equipment 2
x_sol(2)-x_sol(3)
//equipment 3
x_sol(3)+x_sol(9)-x_sol(4)-x_sol(5)
//equipment 4
x_sol(5)-x_sol(6) - x_sol(7)
//equipment 5
x_sol(7)-x_sol(8)
//equipment 6
x_sol(8)-x_sol(9) - x_sol(10)
