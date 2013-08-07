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
getd('../functions');
clear flow_full_rn96_1 flow_rn96_1 comp_full_rn96_1 comp_rn96_1 At_rn96_1 umeas_rn96_1 fixed_rn96_1 red_rn96_1 lower_rn96_1 upper_rn96_1 var_lin_type_rn96_1 constr_lin_type_rn96_1 constr_lhs_rn96_1 constr_rhs_rn96_1  just_measured_rn96_1 observ_rn96_1 non_obs_rn96_1 spec_cand_rn96_1 x_sol f_sol lower_rn96_1 upper_rn96_1 extra xmfull ncomp var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

// In the original paper, streams 2,3 and 9 are unmeasured 
// in these streams (2,3 and 9), values are estimates givem by the paper's original author.
//               1         2        3     4          5        6     7      8       9        10          11  

flow_full_rn96_1 =[3567.2573    1961.9937    1605.2637    2832.9149    715.5544    13.739003    5.0490039    105.01089    2727.904    55.904904    659.64949];
//information of measured/unmeasured(-1)/fixed(-5)
flow_rn96_1 =     [3600 -1   -1   2836.4  730 26  7.6 105  -1  57 -1];
//                1        2        3         4        5        6        7        8          9      10      11  
comp_full_rn96_1 = [1.74311	1.69954	1.79636	0.34051	7.29291	1.00117	4.21013	0.38987	0.33861	89.63081	0.31507
                    18.23989	18.20923	18.27737	0.32437	88.47788	44.92421	43.48957	0.47360	0.31862	0.63543	95.92221
                    0.64090	0.59376	0.69852	0.50307	0.57290	29.93229	7.90183	0.39511	0.50723	0.61756	0.56912
                    3.33081	3.26804	3.40752	3.83246	1.15872	10.93507	8.99940	97.95141	0.20935	3.72621	0.94113
                    60.90894	60.82690	61.00922	76.45582	0.65070	4.94790	29.99257	0.39515	79.38378	1.58240	0.57175
                    14.47437	14.79000	14.08859	18.07731	0.59000	0.01000	0.00000	0.00000	18.77320	0.00000	0.64000
                    0.66198	0.61253	0.72243	0.46646	1.25689	8.24936	5.40650	0.39486	0.46922	3.80759	1.04073]/100;

//information of measured/unmeasured(-1)/fixed(-5)
//            1        2     3         4     5        6        7        8      9        10      11  
comp_rn96_1 =  [ 1.6     -1     -1    0.212    7.42     0.68     4.19     0.11     -1      90.     0.006  
             17.0    -1     -1    0.000      90.6     53.67    43.65    0.08     -1    0.12    98.75  
             0.4     -1     -1     0.212    0.001    30.05    7.91     0.12     -1    -1      -1     
             3.5     -1     -1     3.93     1.17     11.05    9.05     99.54    0.25  -1      0.88   
             61.5    -1     -1    -1         0.79     4.4      30.3     -1       -1.   -1      0.7    
             15.4    -5     -1    19.546    -5       -5       -5       -5       -1.    -5      -5     
             -1      -1     -1    -1        -1       -1       -1       -1        -1    -1      -1  ];
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  
jac = [ 1   -1  -1  0   0   0   0   0    0   0  0    
        0   1   1   -1  -1  -1  -1  0    0   0  0    
        0   0   0   0   1   0   0   0    0   -1 -1    
        0   0   0   1   0   0   0   -1   -1  0  0     ];                                        
//      1   2   3   4   5   6   7   8    9   10  11  
// organizing the vector for the constraints residuals
xmfull=[flow_full_rn96_1(:);matrix(comp_full_rn96_1',-1)];
xm=xmfull;

//the variance proposed by the original author

sd = (0.01*xmfull);
//recalculate sd
for i=1: length(sd)
   if sd(i) <= 0.0001 
       sd(i) = 0.0001;
   end       
end
// to run with equaly weighted standard deviation, uncomment the line below
//sd = ones(size(xmfull,1),size(xmfull,2));

var = sd.^2;
ncomp = 7 ;
//observability/redundancy tests                  
[At_rn96_1,umeas_rn96_1, fixed_rn96_1] = jac_compound_residuals(jac,ncomp,flow_full_rn96_1,comp_full_rn96_1, flow_rn96_1, comp_rn96_1);
[red_rn96_1, just_measured_rn96_1, observ_rn96_1, non_obs_rn96_1, spec_cand_rn96_1] = qrlinclass(At_rn96_1,umeas_rn96_1)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_rn96_1 = setdiff([1:length(xmfull)], umeas_rn96_1);
red = measured_rn96_1;//
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
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower_rn96_1, upper_rn96_1, var_lin_type_rn96_1, constr_lin_type_rn96_1, constr_lhs_rn96_1, constr_rhs_rn96_1]  = structure_compound(jac,ncomp, flow_full_rn96_1,comp_full_rn96_1);
//pause
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
if length(fixed_rn96_1) > 0 then
    lower_rn96_1(fixed_rn96_1) = xmfull(fixed_rn96_1);
    upper_rn96_1(fixed_rn96_1) = xmfull(fixed_rn96_1);
end
xm_init = xmfull./var;
tic
[x_sol, f_sol, extra] = ipopt(xmfull, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type_rn96_1, constr_lin_type_rn96_1, constr_rhs_rn96_1, constr_lhs_rn96_1, lower_rn96_1, upper_rn96_1, params);
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
matrix([1:88],11,8)';
