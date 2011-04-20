// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Yang, Youqi, Rongbo Ten, and Luiqun Jao. 1995. 
//“A study of gross error detection and data reconciliation in process industries.” Comp. & Chem. 
//Eng 19S:S217-S222.

//Bibtex Citation

//@article{Yang1995,
//author = {Yang, Youqi and Ten, Rongbo and Jao, Luiqun},
//journal = {Comp. \& Chem. Eng},
//keywords = {combinatory approach,data reconciliation,gross error detection},
//pages = {S217--S222},
//title = {{A study of gross error detection and data reconciliation in process industries}},
//volume = {19S},
//year = {1995}
//}
// 8 Streams
// 4 Equipments 

function [x_sol, f_sol, status]=P5(xm, sd)

//The jacobian of the constraints
//      1    2  3    4    5   6   7    8  
jac = [ 1   -1  0    0    0   0   -1   0
        0   1  -1    0    0   0    0   1
        0   0   1    1   -1   0    0   0 
        0   0   0    0    0   -1   1  -1];
//      1    2  3    4    5   6   7    8
// From here on, the problem generation is automatic
// No need to edit below
//The problem size: nc = number of constraints and an number of variables
[nc,nv] = size(jac);
// index of the non-zero elements of the Jacobian
[i1,i2]=find(jac<>0);

nonz = nnz(jac);

function f = objfun ( x )

	f = sum(((x-xm).^2)./sd);

endfunction

function c = confun(x)

	c = jac*x;

endfunction

////////////////////////////////////////////////////////////////////////
// Define gradient and Hessian matrix

function gf = gradf ( x )

gf=2*(x-xm)./sd;

endfunction

function H = hessf ( x )

	H = diag(2*ones(nv,1)./sd);
endfunction

function y = dg1(x)

for i = 1: nonz; 
  y(i)=jac(i1(i),i2(i)); 
end

endfunction

function H = Hg1(x)
H = zeros(nv,nv);
endfunction

// The Lagrangian
function y = dh(x,lambda,obj_weight)
	y = obj_weight * hessf ( x ) + lambda * Hg1(x)
endfunction

// The constraints
function y=dg(x)

	y = dg1(x)
	
endfunction


// The sparsity structure of the constraints

sparse_dg = [i1', i2']

// The sparsity structure of the Lagrangian
// the Hessian for this problem is diagonal
sparse_dh = [ [1:nv]', [1:nv]']

// the variables have lower bounds of 0
lower = zeros(nv,1);
// the variables have upper bounds of 50000
upper = 50000*ones(nv,1);
var_lin_type(1:nv) = 1; // Non-Linear
constr_lin_type (1:nc) = 0; // Non-Linear

// the constraints has lower bound of 0
constr_lhs(1:nc) = 0;
// the constraints has upper bound of 0.
constr_rhs(1:nc) = 0;

params = init_param();
// We use the given Hessian
params = add_param(params,"hessian_approximation","exact");
//params = add_param(params,"derivative_test","first-order");
params = add_param(params,"tol",1e-8);
params = add_param(params,"acceptable_tol",1e-8);
params = add_param(params,"mu_strategy","adaptive");

params = add_param(params,"journal_level",0);

[x_sol, f_sol, extra] = ipopt(xm, objfun, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

status = extra('status');
x_sol = x_sol';
endfunction

function [jac]=jacP5()
//The jacobian of the constraints
//      1    2  3    4    5   6   7    8  
jac = [ 1   -1  0    0    0   0   -1   0
        0   1  -1    0    0   0    0   1
        0   0   1    1   -1   0    0   0 
        0   0   0    0    0   -1   1  -1];
//      1    2  3    4    5   6   7    8
endfunction

