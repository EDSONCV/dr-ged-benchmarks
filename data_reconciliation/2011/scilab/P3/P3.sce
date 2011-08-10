// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Atmospheric tower example from:
// Zhang, Zhengjiang and Shao, Zhijiang and Chen, Xi and Wang, Kexin and Qian, Jixin, 2010
// Quasi-weighted least squares estimator for data reconciliation.
// Computers and Chemical Engineering 34: 154-162.

//Bibtex Citation
//@article{Zhang2010,
//author = {Zhang, Zhengjiang and Shao, Zhijiang and Chen, Xi and Wang, Kexin and Qian, Jixin},
//isbn = {0098-1354},
//journal = {Computers \& Chemical Engineering},
//keywords = {Data reconciliation,Gross error detection,Quasi-weighted least squares,Robust estimator,Weighted least squares},
//month = feb,
//number = {2},
//pages = {154--162},
//title = {{Quasi-weighted least squares estimator for data reconciliation}},
//url = {http://www.sciencedirect.com/science/article/B6TFT-4XDCHNS-1/2/63a2b79a4cc89a3afb57ff83c4063242},
//volume = {34},
//year = {2010}
//}

//12 Streams 
//3 Equipments
clear xm var jac nc nv i1 i2 nnz sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
// the measures
xm =[190.26000
174.71000
3.11570
32.45800
33.47400
7.35190
0.31392
92.43200
29.11800
24.27500
18.22500
54.83800
];
//the variance proposed by the original author
// in original paper the standard deviation is given, so it must be squared.
var = [2.14120
1.94840
0.03410
0.34510
0.40010
0.08680
0.00356
1.05940
0.36070
0.30010
0.20011
0.64530
].^2;
//var = ones(12,1);
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
nnz = size(i1,2)

function f = objfun ( x )

    f = sum(((x-xm).^2)./var);

endfunction

function c = confun(x)

    c = jac*x;

endfunction

////////////////////////////////////////////////////////////////////////
// Define gradient and Hessian matrix

function gf = gradf ( x )

    gf=2*(x-xm)./var;

endfunction

function H = hessf ( x )

    H = diag(2*ones(nv,1)./var);
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
params = add_param(params,"tol",1e-3);
params = add_param(params,"acceptable_tol",1e-3);
params = add_param(params,"constr_viol_tol",1e-3);
params = add_param(params,"acceptable_constr_viol_tol",1e-3);
params = add_param(params,"mu_strategy","monotone");

params = add_param(params,"journal_level",5);

[x_sol, f_sol, extra] = ipopt(xm, objfun, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

mprintf("\n\nSolution: , x\n");
for i = 1 : nv
    mprintf("x[%d] = %e\n", i, x_sol(i));
end


mprintf("\n\nObjective value at optimal point\n");
mprintf("f(x*) = %e\n", f_sol);

