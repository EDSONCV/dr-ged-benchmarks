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

clear xm sd jac nc nv i1 i2 nnz sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
// the measures
xm =[28.06
3.06
5.28
8.89
11.39
3.89
4.17
2.78
5.83
4.44
3.89
15.0
13.33
];
//the variance proposed by the original author
sd = [0.075618
0.002498
1.673401
0.021084
7.787169
0.908304
1.043769
0.002062
0.009067
0.005259
0.004037
0.021609
0.017065
];
//the variance proposed by this work 
//sd = ones(13,1);
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