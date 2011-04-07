// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Rao, R Ramesh, and Shankar Narasimhan. 1996.
//“Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” 
//Industrial & Engineering Chemistry Research 35:1362-1368. 
//http://dx.doi.org/10.1021/ie940538b.
//Bibtex Citation

//@article{Rao1996,
//author = {Rao, R Ramesh and Narasimhan, Shankar},
//isbn = {0888-5885},
//journal = {Industrial \& Engineering Chemistry Research},
//month = apr,
//number = {4},
//pages = {1362--1368},
//publisher = {American Chemical Society},
//title = {{Comparison of Techniques for Data Reconciliation of Multicomponent Processes}},
//url = {http://dx.doi.org/10.1021/ie940538b},
//volume = {35},
//year = {1996}
//}

// 12 Streams
// 7 Equipments 

clear xm sd jac nc nv i1 i2 nnz sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
// In the original paper, streams 2,3 and 9 are unmeasured 
// in these streams (2,3 and 9), values are estimates givem by the paper's original author.
xm =[3707
1900
1807
2910
737
26
7.6
105
2832
57
668
];
//the variance proposed by the original author
//sd = 0.0001*ones(11,1).^2;
//the variance proposed by this work (must not change the original results, since they are all equal)
sd =(0.0001*ones(11,1)).^2;
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  
jac = [ 1   -1  -1  0   0   0   0   0    0   0  0    
        0   1   1   -1  -1  -1  -1  0    0   0  0    
        0   0   0   0   1   0   0   0    0   -1 -1    
        0   0   0   1   0   0   0   -1   -1  0  0     ];                                
//      1   2   3   4   5   6   7   8    9   10  11  
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



