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

// 11 Streams
// 4 Equipments 
// 7 Compounds
getd('.');
clear flow_full_P8 flow_P8 comp_full_P8 comp_P8 At_P8 umeas_P8 fixed_P8 red_P8 lower_P8 upper_P8 var_lin_type_P8 constr_lin_type_P8 constr_lhs_P8 constr_rhs_P8  just_measured_P8 observ_P8 non_obs_P8 spec_cand_P8 x_sol f_sol lower_P8 upper_P8 extra xmfull ncomp var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

// In the original paper, streams 2,3 and 9 are unmeasured 
// in these streams (2,3 and 9), values are estimates givem by the paper's original author.
//           1     2   3     4   5   6   7   8    9  10  11  
flow_full_P8 =[3707 1900 1807 2910 737 26  7.6 105 2832 57 668];
//information of measured/unmeasured(-1)/fixed(-5)
flow_P8 =     [3707 -1   -1   2910 737 26  7.6 105  -1  57 -1];
//                1        2        3         4        5        6        7        8          9      10      11  
comp_full_P8 = [   1.64     1.64     1.64     0.0052    7.42     0.68     4.19     0.11     0.0013    90.     0.006  
                17.66    17.66    17.66    0.003     90.6     53.67    43.65    0.08     0.00016    0.12    98.75  
                0.38     0.38     0.38     0.19      0.001    30.05    7.91     0.12     0.18      1.0E-3   1.0E-3     
                3.44     3.44     3.44     3.68      1.17     11.05    9.05     99.54    0.25      4.4     0.88   
                61.47    61.47    61.47    77.08     0.79     4.4      30.3     1.0E-3    80.       1.8     0.7    
                15.39    1.0E-3  31.57    18.01     1.0E-3   1.0E-3   1.0E-3   1.0E-3     19.    1.0E-3   1.0E-3     
                0.41     0.03     0.03     0.29      0.87     0.12     4.9      1.0E-3    0.2     1.0E-3    0.88   ]/100;
//information of measured/unmeasured(-1)/fixed(-5)
//            1        2     3         4        5        6        7        8      9        10      11  
comp_P8 =  [   1.64     -1     -1     0.0052    7.42     0.68     4.19     0.11     -1      90.     0.006  
             17.66    -1     -1    0.003      90.6     53.67    43.65    0.08     -1    0.12    98.75  
             0.38     -1     -1     0.19      0.001    30.05    7.91     0.12     -1    -1      -1     
             3.44     -1     -1     3.68      1.17     11.05    9.05     99.54    0.25  -1      0.88   
             61.47    -1     -1    -1         0.79     4.4      30.3     -1       -1.   -1      0.7    
             15.39    -5     -1    18.01      -5       -5       -5       -5       -1.    -5      -5     
             -1      -1       -1    -1        -1       -1       -1       -1        -1    -1      -1  ];
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  
jac = [ 1   -1  -1  0   0   0   0   0    0   0  0    
        0   1   1   -1  -1  -1  -1  0    0   0  0    
        0   0   0   0   1   0   0   0    0   -1 -1    
        0   0   0   1   0   0   0   -1   -1  0  0     ];                                        
//      1   2   3   4   5   6   7   8    9   10  11  
// organizing the vector for the constraints residuals
xmfull=[flow_full_P8(:);matrix(comp_full_P8',-1)];
xm=xmfull;

//the variance proposed by the original author
sd = (0.01*xmfull);
//recalculate sd
for i=1: length(sd)
   if sd(i) <= 0.0001 
       sd(i) = 0.0001;
   end       
end
var = sd.^2;
ncomp = 7 ;
//observability/redundancy tests                  
[At_P8,umeas_P8, fixed_P8] = jac_compound_residuals(jac,ncomp,flow_full_P8,comp_full_P8, flow_P8, comp_P8);
[red_P8, just_measured_P8, observ_P8, non_obs_P8, spec_cand_P8] = qrlinclass(At_P8,umeas_P8)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_P8 = setdiff([1:length(xmfull)], umeas_P8);
red = measured_P8;//
// to reconcile with all variables, comment the line above and uncomment bellow
//red = [1:length(xmfull)];




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

// ipopt needs some information about the problem, such as jacobian and hessian structure, 
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower_P8, upper_P8, var_lin_type_P8, constr_lin_type_P8, constr_lhs_P8, constr_rhs_P8]  = wls_structure_compound(jac,ncomp, flow_full_P8,comp_full_P8);

params = init_param();
// We use the given Hessian
//params = add_param(params,"hessian_constant","yes");
params = add_param(params,"hessian_approximation","exact");
//params = add_param(params,"hessian_approximation","limited-memory");
//params = add_param(params,"derivative_test","first-order");
//params = add_param(params,"derivative_test_print_all","yes");
//params = add_param(params,"derivative_test","second-order");
params = add_param(params,"tol",1e-2);
params = add_param(params,"acceptable_tol",1e-2);
//params = add_param(params,"mu_strategy","adaptive");
params = add_param(params,"mu_strategy","monotone");
params = add_param(params,"journal_level",5);
params = add_param(params,"fixed_variable_treatment", "relax_bounds");
disp('begore start ipopt')
//according to the original paper, we fix some variables
lower_P8(fixed_P8) = xmfull(fixed_P8);
upper_P8(fixed_P8) = xmfull(fixed_P8);
xm_init = xmfull./var;
tic
[x_sol, f_sol, extra] = ipopt(xmfull, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type_P8, constr_lin_type_P8, constr_rhs_P8, constr_lhs_P8, lower_P8, upper_P8, params);
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
