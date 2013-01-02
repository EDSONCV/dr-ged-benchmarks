// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Rosenberg, J, R S H Mah, and C Iordache. 1987. 
//“Evaluation of Schemes for Detecting and Identifying Gross Errors in Process Data.”
// Ind. & Eng. Chem. Proc. Des. Dev, Vol. 26: 555-564.

//Bibtex Citation

//@Article{ Rosenberg1987,
//author = "J Rosenberg and R S H Mah and C Iordache",
//journal = "Ind. \& Eng. Chem. Proc. Des. Dev, Vol.",
//pages = "555--564",
//title = "{Evaluation of Schemes for Detecting and Identifying Gross Errors in Process Data}",
//volume = "26",
//year = "1987"
//}
// 7 Streams
// 4 Equipments 

clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

xm =[5.026
14.722
14.721
4.968
10.088
5.131
4.835
];
//the variance
var=[1
1
1
1
1
1
1];

//The jacobian of the constraints
//      1   2   3   4   5   6   7 
jac = [ 1   -1   0  1   0   1   0  
        0   1    -1  0  0   0   0 
        0   0    1  -1  -1   0   0
        0   0    0  0   1   -1   -1];                                
//      1   2   3   4   5   6   7  
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


