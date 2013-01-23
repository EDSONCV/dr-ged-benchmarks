// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Kim, I-W; Liebman, M J;Edgar, T F
//“Robust error- in-variables estimation using nonlinear programming techniques.” 
//AIChE J. 1990, 36 (7), 985-993.
// 

clear xm var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs flow_full flow temp_full temp coef
getd('.');
getd('../functions');
ndata=10;    
// parameters
T_ref = 800;
R_gas = 8.314;
tau = 100;//seconds
itau = 1/tau;
deltaH_r = -4180; //J/mol A
rho = 1; // g/l
Cp = 4.18; // J/(g.K )
// data
CAo_x =[0.9871	1.0003	    1.0039	0.976	1.0129	1.0083	1.0075	0.9994	1.0007	0.9973];
CA_x = [0.8906	0.835	     0.8255	0.802	0.752	0.7193	0.6861	0.6388	0.597	0.558];
CB_x = [0.1157	0.138	     0.185	0.2005	0.242	0.2739	0.3215	0.3741	0.3926	0.4703];
To_x = [547.47	531.77	    512.21	490.59	464.67	438.47	408.04	375.56	340.26	306.55];
T_x =  [663.48	676.04	    684.81	695.47	703.69	714.9	726.09	735.44	745.7	753.94];

teta = [0.1 ;10]

xm =[CAo_x; CA_x; CB_x; To_x; T_x]';

xmfull =[xm(:); teta(:)];    

xmrow = xmfull(:);

// the standard deviation
// 
stdz1 =  0.01*ones(1,ndata);
stdz2 = 1*ones(1,ndata);

// organize the standard deviation according to the order in xm
stdDevAll= [stdz1; stdz1; stdz1; stdz2; stdz2]';

stdDevAllrow = stdDevAll(:);
                    
// return the problem structure (jacobian, hessian, number of non-zeros, variable type, etc)

[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  =  wls_structure_CSTR(xmfull);

params = init_param();
// We use the given Hessian
//params = add_param(params,"hessian_constant","yes");
params = add_param(params,"hessian_approximation","exact");
//params = add_param(params,"hessian_approximation","limited-memory");
// uncheck bellow to test derivatives
//params = add_param(params,"derivative_test","second-order");
//params = add_param(params,"derivative_test","first-order");
params = add_param(params,"tol",1e-8);
params = add_param(params,"acceptable_tol",1e-8);
//params = add_param(params,"mu_strategy","monotone");
params = add_param(params,"mu_strategy","adaptive");
//params = add_param(params,"max_iter",100);
params = add_param(params,"journal_level",5);
disp('begore start ipopt')
tic
// if the user want to use random initial guess, uncomment 2 lines bellow and comment the 3rd line
//xrand = rand(30,1);
//[x_sol, f_sol, extra] = ipopt(xrand, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);
isopt=%T;
[x_sol, f_sol, extra] = ipopt(xmfull*1.001, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

toc

mprintf("\n\nSolution: , x\n");
for i = 1 : nv
    mprintf("x[%d] = %e\n", i, x_sol(i));
end

mprintf("\n\nObjective value at optimal point\n");
mprintf("f(x*) = %e\n", f_sol);
