// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Fictitious but realistic mineral processing plant
//Alhaj-Dibo, Moustapha, Didier Maquin, and José Ragot. 2008.
//Data reconciliation: A robust approach using a contaminated distribution.
//Control Engineering Practice 16, no. 2 (February): 159-170.
// http://www.sciencedirect.com/science/article/B6V2H-4N4406D-1/2/50cac92b050f160a20a795faec990dc7.

//Bibtex Citation

//@article{Alhaj-Dibo2008,
//author = {Alhaj-Dibo, Moustapha and Maquin, Didier and Ragot, Jos\'{e}},
//isbn = {0967-0661},
//journal = {Control Engineering Practice},
//keywords = {Data reconciliation,Gross error detection,Linear and bilinear mass balances,Robust estimation},
//month = feb,
//number = {2},
//pages = {159--170},
//title = {{Data reconciliation: A robust approach using a contaminated distribution}},
//url = {http://www.sciencedirect.com/science/article/B6V2H-4N4406D-1/2/50cac92b050f160a20a795faec990dc7},
//volume = {16},
//year = {2008}
//}

// 16 Streams
// 9 Equipments 
// the measures
clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

xm =[24.7
26.5
29.2
1.8
18.32
22.02
20.8
9.43
8.01
4.14
6
6.56
1.04
7.38
4.99
7.69]
// in original paper the standard deviation is given. so it must be squared.
var=[1
1.325
1.46
0.20
0.916
1.101
1.04
0.472
0.401
0.207
0.3
0.328
0.052
0.369
0.25
0.385
].^2;
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    
jac = [ 1   -1  0   1   0   0   0   0    0   0   0   0   0     0     0      0    
        0   1   -1  0   0   0   0   0    0   0   -1   0   0     0     0      0    
        0   0   1   -1  -1  0   0   0    0   0   0   0   0     0     0      0    
        0   0   0   0   1   -1  0   0    0   1   0   0   0     0     0      0    
        0   0   0   0   0   1   -1  -1   0   0   0   0   0     0     0      0   
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    
        0   0   0   0   0   0   1   0    -1   -1  0   0   0     0     0      0    
        0   0   0   0   0   0   0   0    0   0   1   -1   -1    0     0      1    
        0   0   0   0   0   0   0   0    0   0   0   1   1     -1     0       0
        0   0   0   0   0   0   0   0    0   0   0   0   0     1     -1      -1];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    

//observability/redundancy tests                  
umeas_P13 = [];
[red_P13, just_measured_P13, observ_P13, non_obs_P13, spec_cand_P13] = qrlinclass(jac,umeas_P13)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_P13 = setdiff([1:length(xm)], umeas_P13);
red = measured_P13;//
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
