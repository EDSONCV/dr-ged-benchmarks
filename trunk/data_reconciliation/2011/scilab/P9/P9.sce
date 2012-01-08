// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Mandel, Denis, Ali Abdollahzadeh, Didier Maquin, and Jos� Ragot. 1998. 
//Data reconciliation by inequality balance equilibration: a LMI approach. 
//International Journal of Mineral Processing 53, no. 3 (April): 157-169. 
//http://www.sciencedirect.com/science/article/B6VBN-3VM1X8N-3/2/8bffe94a1153eea8647eed5af0031d36.

//Bibtex Citation
//@article{Mandel1998,
//author = {Mandel, Denis and Abdollahzadeh, Ali and Maquin, Didier and Ragot, Jos�},
//isbn = {0301-7516},
//journal = {International Journal of Mineral Processing},
//keywords = {Linear Matrix Inequality Techniques,data reconciliation,error detection,error isolation},
//month = apr,
//number = {3},
//pages = {157--169},
//title = {{Data reconciliation by inequality balance equilibration: a LMI approach}},
//url = {http://www.sciencedirect.com/science/article/B6VBN-3VM1X8N-3/2/8bffe94a1153eea8647eed5af0031d36},
//volume = {53},
//year = {

// 12 Streams
// 5 Equipments 
// the measures
clear xm var jac nc nv i1 i2 nnz sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
xm =[250.5
21.6
207
36.5
172.5
17.6
144.9
47
212
96
90.5
47.7

];
//the variance proposed by the original author and present work

sd = [37.575
1.08
5
1.825
2
0.88
7.245
1
5
2
18.1
2.385
];
var=sd.^2;
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  
jac = [ 1  -1  -1   0   0   0   0   0    0   0   0   0   
        0   0   1   -1  -1  0   0    0   0   0   0   0   
        0   0   0   0   1   -1  -1  0    0   0   0   0 
        0   0   0   0   0   0   1    1    -1  0   0   0
        0   0   0   0   0   0   0   -1    0  1   0   -1
        0   0  0   0   0   0   0   0    1   -1  -1  0
        ];                              
//      1   2   3   4   5   6   7   8    9   10  11  12  

[nc, nv, i1, i2, nnz, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = wls_structure(jac);

params = init_param();
// We use the given Hessian
params = add_param(params,"hessian_approximation","exact");
params = add_param(params,"derivative_test","second-order");
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