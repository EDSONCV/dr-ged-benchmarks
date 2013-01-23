// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Rod, Vladmir and Hancil, Vladislav
//Iterative Estimation of Model Parameters When Measurements of all Variables are Subject to Error// 
// Case 2: Gas phase catalitic hydrogenation of phenol:
// 2 independent variables
// 1 dependent variable
// 28 data points

clear xm var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs flow_full flow temp_full temp coef
getd('.');
getd('../functions');
x1=[0.0150 0.0300 0.0450 0.1000 0.1800 0.0150 0.0300 0.0450 0.0450 0.1000 0.1430 0.1670 0.2500 0.3330 0.0300 0.0450 0.1000 0.1800 0.2400 0.3000 0.3600 0.0260 0.0500 0.1000 0.2500 0.1500 0.3330 0.5000];
x2=[0.2350 0.2200 0.2050 0.1500 0.0700 0.4850 0.4700 0.4550 0.4550 0.4000 0.3570 0.3330 0.2500 0.1670 0.7200 0.7050 0.6500 0.5700 0.5100 0.4500 0.3900 0.9740 0.9500 0.9000 0.7500 0.8500 0.6670 0.5000];
y1=[6.2500 4.9000 2.9000 1.7500 0.3000 12.3000 14.0000 5.0000 14.2000 10.8100 7.8100 6.4100 3.9000 3.6000 13.0000 20.0000 19.8100 15.1000 8.9000 7.5000 2.0000 13.0000 30.0000 37.5000 25.0000 31.5000 10.0000 4.0000 ];
pars =[7.27 0.681 1602];
xmfull = [x1, x2, y1, pars]';
sd = [0.0075*(ones(1,56)), 2.5*ones(1,28)]';
//var = var';
// return the problem structure (jacobian, hessian, number of non-zeros, variable type, etc)
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  =  wls_structure_reactor(xmfull);

params = init_param();
// We use the given Hessian
//params = add_param(params,"hessian_constant","yes");
params = add_param(params,"hessian_approximation","exact");
// uncheck bellow to test derivatives
//params = add_param(params,"derivative_test","second-order");
//params = add_param(params,"derivative_test","first-order");
params = add_param(params,"tol",1e-8);
params = add_param(params,"acceptable_tol",1e-8);
//params = add_param(params,"mu_strategy","monotone");
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