// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Authors
//Dovì, V G, and C Solisio. 2001. Reconciliation of censored measurements in chemical processes: an alternative approach. Chemical Engineering Journal 84, no. 3 (December): 309-314. http://www.sciencedirect.com/science/article/B6TFJ-45KNHB1-F/2/199f358469628f600f10b394d2b55a8b.

//Bibtex Citation

//@article{Dovi2001,
//annote = { The importance of considering the censoring of measured data in the reconciliation of process flow rates has been shown in a previous paper [Chem. Eng. Sci. 52 (17) (1997) 3047]. The purpose of the present paper is to introduce a new technique for carrying out the actual reconciliation procedure and compare its significance and performance with those of previous methods. A numerical example shows how nontrivial differences are to be expected.},
//author = {Dov\`{\i}, V G and Solisio, C},
//isbn = {1385-8947},
//journal = {Chemical Engineering Journal},
//keywords = {Censored data,Data reconciliation,Detection limits},
//month = dec,
//number = {3},
//pages = {309--314},
//title = {{Reconciliation of censored measurements in chemical processes: an alternative approach}},
//url = {http://www.sciencedirect.com/science/article/B6TFJ-45KNHB1-F/2/199f358469628f600f10b394d2b55a8b},
//volume = {84},
//year = {2001}
//}
// 6 Streams
// 3 Equipments 
clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

// the measures
xm =[11
9.9
21.25
11.1
7
3.8
];
//the variance proposed by this work
var = [0.032
0.026
0.120
0.033
0.052
0.015];

// gross error
gerror = zeros(length(xm),1);
// to setup gross errors, select the stream and magnitude as the line bellow
//gerror(2) = 9*sqrt(var(2));
xm = xm + gerror;

//The jacobian of the constraints
jac = [ 1   1  -1    0  0   0
        0   -1  1   -1  0   0
        0   0   0    1  -1  -1 ];
        
//observability/redundancy tests                  
umeas_P2 = [];
[red_P2, just_measured_P2, observ_P2, non_obs_P2, spec_cand_P2] = qrlinclass(jac,umeas_P2)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_P2 = setdiff([1:length(xm)], umeas_P2);
red = measured_P2;//
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


