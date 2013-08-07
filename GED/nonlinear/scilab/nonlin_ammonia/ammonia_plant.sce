// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

// Ammonia plant based on paper from and text book from 

// 11 Streams
// 6 Equipments 
// 5 Compounds
// each stream has the following structure:
//P(bar_abs), F(kgmole/h), T(oC), x_H2 (molar fraction) , x_N2 , x_Ar , x_CH4 , x_NH3
// which are splitted in 2 parts: P, F, T and  x_H2 , x_N2 , x_Ar , x_CH4 , x_NH3


getd('.');
getd('../functions');
clear tstraoflow_full tstraoflow tstraocomp_full tstraocomp At umeas fixed red lower upper var_lin_type constr_lin_type constr_lhs constr_rhs  just_measured observ non_obs spec_cand x_sol f_sol lower upper extra xmfull ncomp var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
// exact simulation data
//                  1      2       3           4         5             6         7          8       9        10         11                
stream_main_full = [10	210	210	     200	    200	    199	    199	    199	199     	210	  199
                    907	907	1565.87	1565.87	1208.40	1208.40	844.71	185.84  658.87	658.87	363.69
                    27	453.9	252.3	 410.0	  410.0	-34.0	  -34.0	  -34.0   -34.0	-30	-34.0];

//Measurement data
//                  1      2       3           4         5          6        7       8       9        10         11                
//stream_main_full = [10	212	      210	     209	    205  	    203	    199	    199	    199   	  210	    199
//                    900	907	     1565.87	1575.87	   1200.40	1218.40  	854.71	189.84	640.87	  668.87	373.69
//                    27	443.9	 250.3	     412.0	    412.0	  -32.0	   -34.0	-34.0	-34.0 	  -30.0	    -34.0];
//

//information of measured/unmeasured(-1)/fixed(-5)
//                  1      2       3     4       5          6          7       8     9       10         11                
stream_main = [     -5	 210	-1	     -1	    205       199	      -1     -1	    -1       -5  	  -1
                    907	-1	    -1	   1575.87	-1	      -1	      -1     189.84 640.87	  -1    	373.69
                    -5	443.9	250.3	 412.0	412.0	  -32.0	      -1	  -1   	-1        -29   	-1];
                  

//                   1       2         3        4          5         6         7        8         9         10         11
stream_comp_full = [ 0.74	0.74   0.7293	0.7293	0.5013	0.5013	0.7146	0.7146	0.7146	0.7146	0.0060
                    0.24	0.24	0.2237	0.2237	0.1420	0.1420	0.2013	0.2013	0.2013	0.2013	0.0042
                    0.01	0.01	0.0229	0.0229	0.0297	0.0297	0.0407	0.0407	0.0407	0.0407	0.0041
                    0.01	0.01	0.0214	0.0214	0.0277	0.0277	0.0371	0.0371	0.0371	0.0371	0.0060
                    1.0e-4   1.0e-4	  0.0026	0.0026	0.2993	0.2993	0.0063	0.0063	0.0063	0.0063	0.9797 ];

//information of measured/unmeasured(-1)/fixed(-5)
//             1       2     3       4      5     6      7        8    9    10       11
stream_comp =[ 0.74	-1	-1	0.729	-1	-1	0.715	-1	-1	-1	0.006
               0.24	-1	-1	0.224	-1	-1	0.201	-1	-1	-1	0.004
               0.01	-1	-1	0.023	-1	-1	0.041	-1	-1	-1	0.004
               0.01	-1	-1	0.021	-1	-1	0.037	-1	-1	-1	0.006
               1.0e-4  -1	-1	0.003	-1	-1	0.006	-1	-1	-1	0.980];
// Q heater 1 , Q heater 2   reaction advance
//x_end = [ 7378754; 25921810; 357.472/2];
x_end = [ 7384260; 26681099; 357.472/2];
// organizing the vector for the constraints residuals
xmfull=[matrix(stream_main_full',-1);matrix(stream_comp_full',-1);x_end];
xm=[matrix(stream_main',-1);matrix(stream_comp',-1);[-1;-1;-1]];
// Ki calculation
// Ki were determined by simulation using iiSE process simulator - valid range: -24 to  -44
// Ki(T(4)) =  a*T(4) + b
//[a; b]
//        [H2     N2,     Ar,      CH4,     NH3]
K_coef = [-3.332	-0.789	-0.048	-0.028	0.000320;
        15.1240	19.0240	6.8367	0.4276	0.076725373]';
// Cp calculation for stream that is heated before the reactor        
// Cp_t = 1.98*(A1 + B1.*tt + C1.*(tt).^2 + D1.* (tt).^-2) (in cal*mol^(-1)*K^(-1))
A1 = [3.249  3.28    2.51  1.702     3.578]
B1 = [0.422  0.593   0      9.081     3.020]*1.0e-3
C1 = [0      0       0     -2.164     0 ]*1.0e-6
D1 = [0.083  0.040   0        0      -0.186]*1.0e5
cp_1 = [A1;B1;C1;D1];
// coeficients for calculation of deltaH in the second heat exchanger 
// These values were adjusted from simulation data
hh = [25.7407  2212.51;
   49.4863 -7630.96]; 

//the variance 
// Pressures are fixed but variances are set to 2 bar
// 5% for flowrates
// 3 oC for temperatures
// 2 % for mole fractions
// for Qheater 1, 2 and reaction advance the var is set to 1 (they do not participate in the reconciliation)
var = zeros(91,1);
var(1:11) = 2*ones(11,1);
var(12:22) = 0.05*xmfull(12:22); 
var(23:33) = 3*ones(11,1);
var(34:$ - 3) = 0.02*xmfull(34:$-3); 
var($-2:$) = ones(3,1);

[At, umeas, fixed] = jac_flowsheet_residuals(xm, xmfull, K_coef, cp_1, hh, 0.22);

[red, just_measured, observ, non_obs, spec_cand] = qrlinclass(At,umeas);

//pause
// reconcile with all measured to reconcile with only redundant variables, uncomment the "red" assignments
measured = setdiff([1:91], umeas);
// to reconcile with all variables, uncomment bellow
//measured = [1:91];
red = measured;
nmeasured = length(measured);

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

// return the problem structure (jacobian, hessian, number of non-zeros, variable type, etc)
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  =  structure_ammonia(xm, xmfull, K_coef, cp_1, hh, 0.22);

lower(fixed) = xmfull(fixed);
upper(fixed) = xmfull(fixed);

params = init_param();
// We use the given Hessian
//params = add_param(params,"hessian_constant","yes");
//params = add_param(params,"hessian_approximation","exact");
// uncheck bellow to test derivatives
//params = add_param(params,"derivative_test","second-order");
//params = add_param(params,"derivative_test","first-order");
params = add_param(params,"tol",1e-1);
params = add_param(params,"acceptable_tol",1e-1);
//params = add_param(params,"mu_strategy","monotone");
//params = add_param(params,"expect_infeasible_problem","yes");
//params = add_param(params,"expect_infeasible_problem","yes");
//params = add_param(params,"mu_strategy","adaptive");
params = add_param(params,"fixed_variable_treatment", "relax_bounds")
params = add_param(params,"journal_level",5);
params = add_param(params,"max_iter",60);
disp('begore start ipopt')
tic
// if the user want to use random initial guess, uncomment 2 lines bellow and comment the 3rd line
//xrand = rand(30,1);
//[x_sol, f_sol, extra] = ipopt(xrand, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

[x_sol, f_sol, extra] = ipopt(xmfull, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

toc

//mprintf("\n\nSolution: , x\n");
//for i = 1 : nv
//    mprintf("x[%d] = %e\n", i, x_sol(i));
//end
//
//mprintf("\n\nObjective value at optimal point\n");
//mprintf("f(x*) = %e\n", f_sol);

[[1:83]' constr_rhs constr_lhs flowsheet_residuals(x_sol, K_coef, cp_1, hh, 0.22)]
xc2 = matrix(xmfull(1:$-3), 11,8)'
xc1 = matrix(x_sol(1:$-3), 11,8)'


