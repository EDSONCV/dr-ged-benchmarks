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
//
//
// 16 Streams
// 9 Equipments 
// 3 Compounds
getd('.');
getd('../functions');
clear flow_full_al flow_al comp_full_al comp_al At_al umeas_al fixed_al red_al lower_al upper_al var_lin_type_al constr_lin_type_al constr_lhs_al constr_rhs_al  just_measured_al observ_al non_obs_al spec_cand_al x_sol f_sol lower_al upper_al extra xmfull ncomp var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

//               1  2   3   4   5   6   7   8   9  10  11  12   13  14 15  16  
flow_full_al =[25	27	22	2	20	24	14	10	10	4	5	7	3	10	5	5];
//                1      2      3         4        5     6        7       8      9      10      11      12      13       14     15      16
comp_full_al = [ 3.645	3.715	2.683	4.580	2.493	2.077	2.105	2.039	2.946	0.000	8.255	3.518	11.315	5.858	8.255	3.460
                 3.153	3.273	3.646	4.781	3.533	2.944	2.436	3.655	3.410	0.000	1.633	3.958	0.366	2.881	1.633	4.128
                 93.202	93.012	93.671	90.638	93.975	94.979	95.460	94.306	93.644	100.000	90.111	92.523	88.319	91.262	90.111	92.412 ]/100

//information of measured/unmeasured(-1)/fixed(-5)
//               1  2   3   4   5   6   7   8   9  10  11  12   13  14 15  16  
flow_al =     [-1	27	22	-1	20	24	14	10	10	4	-1	7	3	10	5	5];

//                 1        2     3         4     5        6      7        8      9        10    11        12    13       14     15      16
comp_al =  [    3.50	3.58	2.50	4.60	2.29	2.50	2.80	2.08	2.50	250.00	8.34	3.40	11.80	5.92	8.34	3.50
                3.00	3.13	3.50	4.80	3.37	3.37	3.10	3.74	3.00	300.00	1.52	4.20	0.00	2.76	1.52	4.00
                -1       -1	    -1       -1     -1       -1     -1       -1     -1       -1     -1      -1      -1       -1     -1       -1 ];

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

// organizing the vector for the constraints residuals
xmfull=[flow_full_al(:);matrix(comp_full_al',-1)];
xm=xmfull;

//the variance proposed by the original author
// in original paper the standard deviation is given. so it must be squared.
sd_flow = [1.325	1.325	1.46	0.916	0.916	1.101	1.04	0.472	0.401	0.207	0.328	0.328	0.052	0.369	0.25	0.385];
sd_comp_12 = [0.389	0.27	0.252	0.475	0.209	0.246	0.29	0.201	0.371	0.415	0.349	0.363	0.871	0.424	0.349	0.515
             0.344	0.353	0.355	0.5	0.329	0.351	0.374	0.52	0.329	0.455	0.419	0.434	0.632	0.665	0.41	0.518]/100;
sd_comp_3 = [93.202	93.012	93.671	90.638	93.975	94.979	95.460	94.306	93.644	100.000	90.111	92.523	88.319	91.262	90.111	92.412]/1000;             
sd_123 = [sd_comp_12; sd_comp_3];
sd = [sd_flow(:);matrix(sd_123',-1)];
// to run with equaly weighted standard deviation, uncomment the line below
//sd = ones(size(xmfull,1),size(xmfull,2));

var = sd.^2;
ncomp = 3 ;
//observability/redundancy tests                  
[At_al,umeas_al, fixed_al] = jac_compound_residuals(jac,ncomp,flow_full_al,comp_full_al, flow_al, comp_al);
[red_al, just_measured_al, observ_al, non_obs_al, spec_cand_al] = qrlinclass(At_al,umeas_al)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_al = setdiff([1:length(xmfull)], umeas_al);
red = measured_al;//
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
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower_al, upper_al, var_lin_type_al, constr_lin_type_al, constr_lhs_al, constr_rhs_al]  = structure_compound(jac,ncomp, flow_full_al,comp_full_al);
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
if length(fixed_al) > 0 then
    lower_al(fixed_al) = xmfull(fixed_al);
    upper_al(fixed_al) = xmfull(fixed_al);
end
xm_init = xmfull./var;
tic
[x_sol, f_sol, extra] = ipopt(xmfull, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type_al, constr_lin_type_al, constr_rhs_al, constr_lhs_al, lower_al, upper_al, params);
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

