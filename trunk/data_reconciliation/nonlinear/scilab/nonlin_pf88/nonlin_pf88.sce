// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//C.C. David Pai, Gary D. Fisher
//Application of Broyden's Method to Reconciliation of Nonlinearly Constrained Data
// AICHE Journal 1988, V 34 No. 5 -p 873-876
// 8 Variables
// 5 Measurements
// 3 Unmeasured

clear xm var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs flow_full flow temp_full temp coef
getd('.');
getd('../functions');
// This example from Pai & Fisher is a general non-linearly constrained data reconciliation
// measured data
//           x1       x2      x3    x4     x5  v1    v2     v3
//meas_full = [ 4.4     5.5     1.7  1.6     5    8    0.7    1.8]';
meas_full = [ 9.4     5.5     1.7  1.6     5    8    0.7    1.8]';
// information of measured, unmeasured (-1) or fixed variables
//           x1       x2      x3    x4     x5  v1    v2     v3
meas =    [ 4.4       5.5      1.7   1.6    5  -1   -1    -1]';

xmfull =[meas_full(:)];    
xm=xmfull;    
// the variance
var = ones(length(xmfull),1);

//The variable classification
                    
[At, umeas, fixed] =  jac_pf88_residuals(meas_full, meas);

[red, just_measured, observ, non_obs, spec_cand] = qrlinclass(At,umeas)

// reconcile with all measured to reconcile with only redundant variables, uncomment the "red" assignments
measured = setdiff([1:8], umeas);
// to reconcile with all variables, uncomment bellow
//measured = [1:30];
nmeasured = length(measured);

red=measured;

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
exec ../functions/setup_DR.sce;
// to run robust reconciliation, it is also necessary to choose the function to return the problem structure


// return the problem structure (jacobian, hessian, number of non-zeros, variable type, etc)
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  =  pf_struct_nl(xmfull);

params = init_param();
// We use the given Hessian
//params = add_param(params,"hessian_constant","yes");
params = add_param(params,"hessian_approximation","exact");
// uncheck bellow to test derivatives
params = add_param(params,"derivative_test","second-order");
//params = add_param(params,"derivative_test","first-order");
params = add_param(params,"tol",1e-8);
params = add_param(params,"acceptable_tol",1e-8);
params = add_param(params,"mu_strategy","monotone");
//params = add_param(params,"mu_strategy","adaptive");
params = add_param(params,"journal_level",5);
disp('begore start ipopt')
tic
// if the user want to use random initial guess, uncomment 2 lines bellow and comment the 3rd line
//xrand = rand(30,1);
//[x_sol, f_sol, extra] = ipopt(xrand, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

[x_sol, f_sol, extra] = ipopt(xmfull, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

toc

mprintf("\n\nSolution: , x\n");
for i = 1 : nv
    mprintf("x[%d] = %e\n", i, x_sol(i));
end

mprintf("\n\nObjective value at optimal point\n");
mprintf("f(x*) = %e\n", f_sol);
