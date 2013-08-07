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

// 3 Streams
// 1 Equipments
// 8 compounds

clear xm var jac nc nv i1 i2 nnz sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

//                   1     2      3    
flow_full_rn96_3 =[8.50	4.50	4.00];
flow_rn96_3 =     [8.50	4.50	4.00];
//                    1        2       3     
comp_full_rn96_3 = [0.024	0.000	0.050
                    0.024	0.000	0.050
                   60.987	99.889	17.221
                   17.997	0.019	38.221
                   12.998	0.001	27.620
                   4.998	0.000	10.620
                   2.104	0.000	4.472
                   0.869	0.090	1.745]/100;
//               1       2        3    
comp_rn96_3 =  [0.02   -1        -5; 
                0.02   -1        -5;
                60.07  99.94     6.39; 
                18.88  0.02      40.18;
                13.88  1.0E-5    29.31;
                4.95   -5        14.42;
                2.01   -5        6.55;
                0.37   -5        3.14 ]

xmfull_rn96_3 = [flow_full_rn96_3(:);matrix(comp_rn96_3',-1)];
xm = xmfull_rn96_3;

//the variance proposed by the original author
sd = (0.01*ones(12,1));
//recalc variance
for i=1: length(var)
   if sd(i) >= 0.0001 
       sd(i) = 0.0001;
   end       
end
// to run with equaly weighted standard deviation, uncomment the line below
//sd = ones(size(xmfull_rn96_3,1),size(xmfull_rn96_3,2));

var = sd.^2;
//The jacobian of the constraints
//      1   2   3   
jac = [ 1   -1  -1  ];                                        

ncomp=8 
                    
[At_rn96_3,umeas_rn96_3, fixed_rn96_3] = jac_compound_residuals(jac,ncomp,flow_full_rn96_3,comp_full_rn96_3, flow_rn96_3, comp_rn96_3);
[red_rn96_3, just_measured_rn96_3, observ_rn96_3, non_obs_rn96_3, spec_cand_rn96_3] = qrlinclass(At_rn96_3,umeas_rn96_3)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_rn96_3 = setdiff([1:length(xmfull_rn96_3)], umeas_rn96_3);
// to reconcile with all variables, comment the line above and uncomment bellow
//measured = [1:length(xmfull_rn96_3)];
red=measured_rn96_3;

// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS = 0
// Absolute sum of squares = 1
//Cauchy = 2
//Contamined Normal = 3
//Fair  = 4
//Hampel = 5
//Logistic = 6
//Lorenztian = 7
//Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = 0;
exec ../functions/setup_DR.sce;
// to run robust reconciliation, it is also necessary to choose the function to return the problem structure

pause
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = structure_compound(jac,ncomp, flow_full_rn96_3,comp_full_rn96_3);

params = init_param();
// We use the given Hessian
params = add_param(params,"hessian_constant","yes");
params = add_param(params,"hessian_approximation","exact");
//params = add_param(params,"derivative_test","first-order");
params = add_param(params,"tol",1e-2);
params = add_param(params,"acceptable_tol",1e-2);
params = add_param(params,"mu_strategy","monotone");
params = add_param(params,"journal_level",5);
params = add_param(params,"fixed_variable_treatment", "relax_bounds");
disp('begore start ipopt')
//according to the original paper, we fix the measured total flow
if length(fixed_rn96_3) > 0 then
    lower(fixed_rn96_3) = xmfull(fixed_rn96_3);
    upper(fixed_rn96_3) = xmfull(fixed_rn96_3);
end
tic
[x_sol, f_sol, extra] = ipopt(xmfull_rn96_3, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);
toc

mprintf("\n\nSolution: , x\n");
for i = 1 : nv
    mprintf("x[%d] = %e\n", i, x_sol(i));
end

mprintf("\n\nObjective value at optimal point\n");
mprintf("f(x*) = %e\n", f_sol);

[Aeqp, Astreams] =size(jac)
xx=matrix(x_sol,Astreams,ncomp+1)
TotalFlowMeasured = xx(:,1)'
compoundMeasured  = 100*xx(:, 2:$)


