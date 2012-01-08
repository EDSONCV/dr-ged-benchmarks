// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Mitsas, Christos L. 2010. Data reconciliation and variable classification by null space methods. 
//Measurement 43, no. 5 (June): 702-707.
// http://apps.isiknowledge.com/full_record.do?product=UA&search_mode=GeneralSearch&qid=2&SID=2A@bF9dN34I72L1Am9M&page=2&doc=17&colname=WOS.

//Bibtex Citation

//@article{Mitsas2010a,
//author = {Mitsas, Christos L.},
//journal = {Measurement},
//month = jun,
//number = {5},
//pages = {702--707},
//publisher = {ELSEVIER SCI LTD},
//title = {{Data reconciliation and variable classification by null space methods}},
//url = {http://apps.isiknowledge.com/full\_record.do?product=UA\&search\_mode=GeneralSearch\&qid=2\&SID=2A@bF9dN34I72L1Am9M\&page=2\&doc=17\&colname=WOS},
//volume = {43},
//year = {2010}
//}

// 12 Streams
// 6 Equipments 

clear xm var jac nc nv i1 i2 nnz sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

xm =[100.49
29.46
40.41
30.05
15.17
4.98
9.92
70.68
10.03
19.82
81.06
101.57
];
//the variance or uncertainties  proposed by the original author
sd = [1
0.3
0.4
0.3
0.3
0.1
0.2
0.7
0.3
0.4
0.8
1];
var=sd.^2;
//the variance proposed by this work
//var = ones(12,1);
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  
jac = [ 1  -1  -1  -1   0   0   0   0    0   0   0   0   
        0    1  0   0   -1  -1  -1  0    0   0   0   0   
        0    0  1   0   1   1   1   -1   0   0   0   0
        0    0  0   1   0   0   0   0    -1  -1  0   0
        0    0  0   0   0   0   0   1    1   0   -1  0
        0    0  0   0   0   0   0   0    0   1   1   -1
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