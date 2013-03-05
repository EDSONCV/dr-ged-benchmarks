// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Martins, Márcio A.F., Carolina A. Amaro, Leonardo S. Souza, Ricardo A. Kalid, and Asher Kiperstok. 2010. 
//New objective function for data reconciliation in water balance from industrial processes. 
//Journal of Cleaner Production (March): 1-6. doi:10.1016/j.jclepro.2010.03.014. http://linkinghub.elsevier.com/retrieve/pii/S0959652610001149.

//Bibtex Citation
//@article{Martins2010,
//author = {Martins, M\'{a}rcio A.F. and Amaro, Carolina A. and Souza, Leonardo S. and Kalid, Ricardo A. and Kiperstok, Asher},
//doi = {10.1016/j.jclepro.2010.03.014},
//file = {::},
//issn = {09596526},
//journal = {Journal of Cleaner Production},
//month = mar,
//pages = {1--6},
//title = {{New objective function for data reconciliation in water balance from industrial processes}},
//url = {http://linkinghub.elsevier.com/retrieve/pii/S0959652610001149},
//year = {2010}
//}

// 13 Streams
// 8 Equipments 

clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

// the measures
xm =[28.06
3.06
5.28
8.89
11.39
3.89
4.17
2.78
5.83
4.44
3.89
15.0
13.33
];
//the variance proposed by the original author with some modifications
var = [0.075618
0.002498
0.029749
0.021084
0.138439
0.016148
0.018556
0.002062
0.009067
0.005259
0.004037
0.021609
0.017065
];
//the variance proposed by this work 
//var = ones(13,1);
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  13
jac = [ 1  -1  -1  -1  -1   0   0   0    0   0   0   0   0
        0   1   0   0   0   0   0  -1    0   0   0   0   0
        0   0   1   0   0   0   0   0   -1   0   0   0   0
        0   0   0   1   0  -1  -1   0    0   0   0   0   0        
        0   0   0   0   1   0   0   0    0   0   1  -1   0
        0   0   0   0   0   1   0   0    0  -1   0   0   0
        0   0   0   0   0   0   1   0    0   0  -1   0   0
        0   0   0   0   0   0   0   1    1   1   0   0  -1];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13

//observability/redundancy tests                  
umeas_P10 = [];
[red_P10, just_measured_P10, observ_P10, non_obs_P10, spec_cand_P10] = qrlinclass(jac,umeas_P10)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_P10 = setdiff([1:length(xm)], umeas_P10);
red = measured_P10;//
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
