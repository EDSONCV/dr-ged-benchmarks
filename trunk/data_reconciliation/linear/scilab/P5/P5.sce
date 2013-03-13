// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Yang, Youqi, Rongbo Ten, and Luiqun Jao. 1995. 
//“A study of gross error detection and data reconciliation in process industries.” Comp. & Chem. 
//Eng 19S:S217-S222.

//Bibtex Citation

//@article{Yang1995,
//author = {Yang, Youqi and Ten, Rongbo and Jao, Luiqun},
//journal = {Comp. \& Chem. Eng},
//keywords = {combinatory approach,data reconciliation,gross error detection},
//pages = {S217--S222},
//title = {{A study of gross error detection and data reconciliation in process industries}},
//volume = {19S},
//year = {1995}
//}
// 8 Streams
// 4 Equipments 

clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

xm =[110.1
41.1
79
30.6
108.3
19.8
55
39
];
//the variance
var = [3.234
1.803
1.665
0.999
3.129
0.418
3.228
1.623
];

// gross error
gerror = zeros(length(xm),1);
// to setup gross errors, select the stream and magnitude as the line bellow
//gerror(2) = 9*sqrt(var(2));
xm = xm + gerror;


//The jacobian of the constraints
//      1    2  3    4    5   6   7    8  
jac = [ 1   -1  0    0    0   0   -1   0
        0   1  -1    0    0   0    0   1
        0   0   1    1   -1   0    0   0 
        0   0   0    0    0   -1   1  -1];
//      1    2  3    4    5   6   7    8

//observability/redundancy tests                  
umeas_P5 = [];
[red_P5, just_measured_P5, observ_P5, non_obs_P5, spec_cand_P5] = qrlinclass(jac,umeas_P5)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_P5 = setdiff([1:length(xm)], umeas_P5);
red = measured_P5;//

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
x_sol(1)+x_sol(4)-x_sol(6)-x_sol(5)
// equip1
x_sol(1)-x_sol(7)-x_sol(2)
// equip2
x_sol(2)-x_sol(3)+x_sol(8)
// equip3
x_sol(3)-x_sol(5)+x_sol(4)
// equip4
x_sol(7)-x_sol(6)-x_sol(8)

