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
function [x_sol, f_sol, status]=P12(xm, sd)
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
        0   0   0   0   0   0   0   0    0   0   0   1   1     -1     0      0
        0   0   0   0   0   0   0   0    0   0   0   0   0     1     -1      -1];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    
// From here on, the problem generation is automatic
// No need to edit below
//The problem size: nc = number of constraints and nv number of variables
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

function [jac]=jacP12()
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
        0   0   0   0   0   0   0   0    0   0   0   1   1     -1     0      0
        0   0   0   0   0   0   0   0    0   0   0   0   0     1     -1      -1];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    
endfunction

