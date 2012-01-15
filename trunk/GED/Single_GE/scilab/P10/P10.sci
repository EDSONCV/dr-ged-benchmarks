// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Martins, Márcio A.F., Carolina A. Amaro, Leonardo S. Souza, Ricardo A. Kalid, and Asher Kiperstok. 2010. 
//New objective function for data reconciliation in water balance from industrial processes. 
//Journal of Cleaner Production (March): 1-6. doi:10.1016/j.jclepro.2010.03.014. http://linkinghub.elsevier.com/retrieve/pii/S0959652610001149.

//Bibtex Citation
//@article{Martins2010,
//author = {Martins, M\'{a}rcio A.F. and Amaro, Carolina A. and Souza, Leonardo S. and Kalid, Ricardo A. and Kiperstok, Asher},
//doi = {10.1016/j.jclepro.2010.03.014},
//file = {::},
//issn = {09596526},
//journal = {Journal of Cleaner Production},
//month = mar,
//pages = {1--6},
//title = {{New objective function for data reconciliation in water balance from industrial processes}},
//url = {http://linkinghub.elsevier.com/retrieve/pii/S0959652610001149},
//year = {2010}
//}

// 13 Streams
// 8 Equipments 
function [x_sol, f_sol, status]=P9(xms, sd)
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  13
jac = [ 1  -1  -1  -1  -1   0   0   0    0   0   0   0   0
        0   1   0   0   0   0   0  -1    0   0   0   0   0
        0   0   1   0   0   0   0   0   -1   0   0   0   0
        0   0   0   1   0  -1  -1   0    0   0   0   0   0        
        0   0   0   0   1   0   0   0    0   0   1  -1   0
        0   0   0   0   0   1   0   0    0  -1   0   0   0
        0   0   0   0   0   0   1   0    0   0  -1   0   0
        0   0   0   0   0   0   0   1    1   1   0   0  -1];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13
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

// The constraints
function y=dg(x)

	y = dg1(x)
	
endfunction

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
params = add_param(params,"mu_oracle","probing");
params = add_param(params,"hessian_constant","yes");
params = add_param(params,"jac_d_constant","yes");
params = add_param(params,"jac_c_constant","yes");
//params = add_param(params,"fast_step_computation","yes");
//params = add_param(params,"mu_oracle","probing");
//params = add_param(params,"mehrotra_algorithm","yes");
params = add_param(params,"mu_strategy","monotone");
tic
for i=1:size(xms,1)
xm=xms(i,:)'
[x_sol(1:nv,i), f_sol(i), extra] = ipopt(xm, objfun, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);
status(i) = extra('status');
end
b=toc;
printf('toc: %f',b());
x_sol = x_sol';
endfunction

function [jac]=jacP10()
//      1   2   3   4   5   6   7   8    9   10  11  12  13
jac = [ 1  -1  -1  -1  -1   0   0   0    0   0   0   0   0
        0   1   0   0   0   0   0  -1    0   0   0   0   0
        0   0   1   0   0   0   0   0   -1   0   0   0   0
        0   0   0   1   0  -1  -1   0    0   0   0   0   0        
        0   0   0   0   1   0   0   0    0   0   1  -1   0
        0   0   0   0   0   1   0   0    0  -1   0   0   0
        0   0   0   0   0   0   1   0    0   0  -1   0   0
        0   0   0   0   0   0   0   1    1   1   0   0  -1];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13   
endfunction

