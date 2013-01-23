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

ndata=5;    
T_ref = 323.15;
R_gas = 8.314;
P_x  = [483.8	493.2	499.9	501.4	469.7];
y1_x = [0.591	0.602	0.612	0.657	0.814];
T_x  = [50		50	50	50	50];
TR_x = (T_x + 273.15)./(T_ref*ones(1,ndata));
x1_x = [0.3		0.4	0.5	0.7	0.9];
teta = [1.5 ;1.5]
AB_x = [8.314*325.15 8.314*325.15];
xm =[P_x;	y1_x;	TR_x;	x1_x]';

xmfull =[xm(:); teta(:)];    
//xmfull =[xm(:); AB_x(:)];    
xmrow = xmfull(:);

// the standard deviation
// 
stdP =  0.75*ones(1,ndata);
stdy = 0.015*ones(1,ndata);
stdT = 0.000309*ones(1,ndata);
//stdT = 0.1*ones(1,ndata);
stdx = 0.005*ones(1,ndata);
// organize the standard deviation according to the order in xm
stdDevAll= [stdP; stdy; stdT; stdx]';

stdDevAllrow = stdDevAll(:);
                    
// return the problem structure (jacobian, hessian, number of non-zeros, variable type, etc)
isopt=%F;
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  =  wls_structure_ELV_VanL(xmfull);

params = init_param();
// We use the given Hessian
//params = add_param(params,"hessian_constant","yes");
params = add_param(params,"hessian_approximation","exact");
//params = add_param(params,"hessian_approximation","limited-memory");
// uncheck bellow to test derivatives
//params = add_param(params,"derivative_test","second-order");
//params = add_param(params,"derivative_test","first-order");
params = add_param(params,"tol",1e-2);
params = add_param(params,"acceptable_tol",1e-2);
//params = add_param(params,"mu_strategy","monotone");
params = add_param(params,"mu_strategy","adaptive");
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
