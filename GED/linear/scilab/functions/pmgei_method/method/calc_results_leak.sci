// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to weighted least squares functions using Lagrange Multipliers 
// for sensitivity analysys
// user can send to this function the parameter list of ipopt or leave this 
// argument blank for the default ones
// 

function f = objfun_leak ( x )
//disp('inside eval f')
	f = sum(((xm-x).^2)./var);

endfunction

function c = confun_leak(x)

	c = jac*x;

endfunction

////////////////////////////////////////////////////////////////////////
// Define gradient and Hessian matrix

function gf = gradf_leak ( x )
//disp('inside eval gradf')


gf=2*(x - xm)./var;

endfunction

function H = hessf_leak ( x )

	H = (2*ones(nv_leak,1)./var);
endfunction

function y = dg1_leak(x)
//disp('inside eval dg1_leak')
for i = 1: nnzeros_leak; 
  y(i)=jac(i1_leak(i),i2_leak(i)); 
end

endfunction

// The Lagrangian
// In fact, IPOPT does not call this line because we say that the constraints are constants
function y = dh_leak(x,lambda,obj_weight)
	y = obj_weight * hessf_leak ( x );
endfunction

// The constraints
function y=dg_leak(x)

	y = dg1_leak(x)
	
endfunction


// thsi function avoids repetition of the structure evaluation
function [nc, nv, i1, i2, nnzeros, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = wls_structure_leak(jac, leak)
    // From here on, the problem generation is automatic
// No need to edit below
//The problem size: nc = number of constraints and an number of variables
[nc,nv] = size(jac);
// index of the non-zero elements of the Jacobian
[i1,i2]=find(jac<>0);

nnzeros = size(i1,2);

// The sparsity structure of the constraints

sparse_dg = [i1', i2'];

// The sparsity structure of the Lagrangian
// the Hessian for this problem is diagonal
sparse_dh = [ [1:nv]', [1:nv]'];

// the variables have lower bounds of 0
lower = zeros(nv,1);
// the variables have upper bounds of 50000
upper = 50000*ones(nv,1);
var_lin_type(1:nv) = 1; // Non-Linear
constr_lin_type (1:nc) = 0; // Linear

// the constraints has lower bound of 0
constr_lhs(1:nc) = 0;
// the constraints has upper bound of 0.
constr_rhs(1:nc) = 0;
endfunction



function [x_sol_var, status]  = calc_results_leak(xfinal, jac, leak, runsize)

params_leak = init_param();
    // We use the given Hessian
params_leak = add_param(params_leak,"hessian_approximation","exact");
//params = add_param(params,"derivative_test","second-order");
params_leak = add_param(params_leak,"tol",1e-8);
params_leak = add_param(params_leak,"acceptable_tol",1e-8);
params_leak = add_param(params_leak,"mu_strategy","monotone");
params_leak = add_param(params_leak,"journal_level",1);

x_sol_var = zeros(size(xfinal,1),size(xfinal,2));
status = -666*ones(runsize*jac_row,1);
[nc_leak, nv_leak, i1_leak, i2_leak, nnzeros_leak, sparse_dg_leak, sparse_dh_leak, lower_leak, upper_leak, var_lin_type_leak, constr_lin_type_leak, constr_lhs_leak, constr_rhs_leak]  = wls_structure_leak(jac, leak);

runsizefinal = size(xfinal,1);

for i = 1 : size(xfinal,1)
//pause
    xm = xfinal(i,:)';
    constr_rhs_leak = leak(i,:)';
    constr_lhs_leak = leak(i,:)';
//    disp('inside opt loop')

    [x_sol, f_sol, extra] = ipopt(xm, objfun_leak, gradf_leak, confun_leak, dg_leak, sparse_dg_leak, dh_leak, sparse_dh_leak, var_lin_type_leak, constr_lin_type_leak, constr_rhs_leak, constr_lhs_leak, lower_leak, upper_leak, params_leak);


    status(i) = extra.status;
    x_sol_var(i,:) = x_sol';

end

endfunction
