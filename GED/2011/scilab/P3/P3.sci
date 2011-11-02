// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Atmospheric tower example from:
// Zhang, P, G Rong, and Y Wang. 2001.
// A new method of redundancy analysis in data reconciliation and its application.
// Computers and Chemical Engineering 25: 941-949.

//Bibtex Citation

//@article{Zhang2001,
//author = {Zhang, P and Rong, G and Wang, Y},
//journal = {Computers and Chemical Engineering},
//pages = {941--949},
//title = {{A new method of redundancy analysis in data reconciliation and its application}},
//volume = {25},
//year = {2001}
//}

//12 Streams 
//3 Equipments
function [x_sol, f_sol, status]=P3(xm, sd, xr)
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  
jac = [ 1  -1   1   0   0   0   0   0    0   0   -1   0  
        0   1   0  -1  -1  -1  -1  -1    1   0    1   -1  
        0   0   -1  0   0   0   0    0   -1  -1   0   1  ];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  
// No need to edit below
// From here on, the problem generation is automatic
//The problem size: nc = number of constraints and number of variables
[nc,nv] = size(jac);
// index of the non-zero elements of the Jacobian
[i1,i2]=find(jac<>0);
// No need to edit below
nonz = nnz(jac);

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
params = add_param(params,"tol",1e-3);
params = add_param(params,"acceptable_tol",1e-3);
params = add_param(params,"constr_viol_tol",1e-3);
params = add_param(params,"acceptable_constr_viol_tol",1e-3);
params = add_param(params,"mu_oracle","probing");
params = add_param(params,"hessian_constant","yes");
params = add_param(params,"jac_d_constant","yes");
params = add_param(params,"jac_c_constant","yes");
//params = add_param(params,"fast_step_computation","yes");
//params = add_param(params,"mu_oracle","probing");
//params = add_param(params,"mehrotra_algorithm","yes");
params = add_param(params,"mu_strategy","monotone");


params = add_param(params,"journal_level",0);

[x_sol, f_sol, extra] = ipopt(xr, objfun, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

status = extra('status');
x_sol = x_sol';
endfunction

function [jac]=jacP3()
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  
jac = [ 1  -1   1   0   0   0   0   0    0   0   -1   0  
        0   1   0  -1  -1  -1  -1  -1    1   0    1   -1  
        0   0   -1  0   0   0   0    0   -1  -1   0   1  ];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  

endfunction

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