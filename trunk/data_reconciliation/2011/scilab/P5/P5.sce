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

clear xm sd jac nc nv i1 i2 nnz sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
xm =[110.1
41.1
79
30.
108.3
19.8
57.6
37.8
];
//the variance
sd = [0.9870
0.4110
0.7890
0.3020
1.0910
0.1980
0.5760
0.3780];

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

nnz = size(i1,2)

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

for i = 1: nnz; 
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
params = add_param(params,"derivative_test","first-order");
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
// global balance
x_sol(1)+x_sol(4)-x_sol(6)-x_sol(5)
// equip1
x_sol(1)-x_sol(7)-x_sol(2)
// equip2
x_sol(2)-x_sol(3)+x_sol(8)
// equip3
x_sol(3)-x_sol(5)+x_sol(4)
// equip4
x_sol(7)-x_sol(6)-x_sol(8)

