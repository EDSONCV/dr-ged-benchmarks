// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Heat exchanger with by-pass valve
//Narasimhan, S, and C Jordache. 2000.
//Data Reconciliation and Gross Error Detection: An Intelligent Use of Process Data. 1st ed.
//Houston: Gulf Publishing.

//Bibtex Citation

//@book{Narasimhan2000,
//address = {Houston},
//author = {Narasimhan, S and Jordache, C},
//booktitle = {Process Data. Gulf Professional Publishing, Houston, TX.},
//edition = {1},
//publisher = {Gulf Publishing},
//title = {{Data Reconciliation and Gross Error Detection: An Intelligent Use of Process Data}},
//year = {2000}
//}

// 6 Streams
// 4 Equipments 

clear xm var jac nc nv i1 i2 nnz sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
getd('../functions/wls');
xm =[101.91;64.45;34.65;64.2;36.44;98.88];
//the variance
var = ones(6,1);
//The jacobian of the constraints
jac = [ 1  -1  -1    0  0   0
        0   1   0   -1  0   0
        0   0   1    0 -1   0
        0   0   0    1  1  -1];

[nc, nv, i1, i2, nnz, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = wls_structure(jac);

params = init_param();
// We use the given Hessian
params = add_param(params,"hessian_approximation","exact");
params = add_param(params,"derivative_test","second-order");
params = add_param(params,"tol",1e-8);
params = add_param(params,"acceptable_tol",1e-8);
params = add_param(params,"mu_strategy","adaptive");
params = add_param(params,"journal_level",6);
params = add_param(params,"print_user_options","yes");

//[x_sol, f_sol, extra] = ipopt(xm, objfun, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);
[x_sol, f_sol, extra] = ipopt(zeros(6,1), objfun, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

mprintf("\n\nSolution: , x\n");
for i = 1 : nv
    mprintf("x[%d] = %e\n", i, x_sol(i));
end

mprintf("\n\nObjective value at optimal point\n");
mprintf("f(x*) = %e\n", f_sol);
