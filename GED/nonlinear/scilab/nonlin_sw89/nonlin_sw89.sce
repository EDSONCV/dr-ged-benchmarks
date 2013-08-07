// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Swartz, C. L. E, 1898
//“Data reconciliation for generalized flowsheet applications.” 
//197th National Meeting of the American Chemical Society, Dallas, TX. 
// Data extracted from example 5.4 from:
// Romagnoli, José, and Mabel C. Sánchez. 1999. Process Systems Engineering - Volume 2 -
// Data Processing and Reconciliation for Chemical Process Operations. Edited by Elsevier. 
// Modified flows and temperatures to obtain exact mass and energy balance
// 30 Variables
// 15 Flow rates
// 15 Temperatures 

clear xm var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs flow_full flow temp_full temp coef
getd('.');
getd('../functions');
// exact flow measurements or estimates
//           A1       A2      A3    A4      A5       A6     A7      A8      B1     B2    B3     C1    C2    D1       D2  
flow_full =[ 954.4     954.4   401.7   401.7  401.7   552.7  552.7   954.4  253.2 253.2  253.2   308.1 308.1 680.1   680.1 ];
// information of measured, unmeasured (-1) or fixed variables
flow =     [ 1000     -1      401.7  -1     -1       552.7  -1     -1      253.2  -1    -1     308.1  -1   -1       680.1];
// exact tempeature measurements or estimates
//           A1        A2      A3    A4      A5       A6     A7      A8      B1     B2     B3     C1      C2     D1       D2  
temp_full =[ 473.8 491.78  491.78   530.09 616.31   491.78  619    617.87  627.71  570.00 504.84  699.99  600.74 662.82  558.34];
// information of measured, unmeasured (-1) or fixed variables
//           A1        A2      A3    A4      A5       A6     A7      A8      B1     B2     B3     C1      C2     D1       D2  
temp =     [ 466.33     -1    481.78  530.09 616.31   -1      619    614.92  618.11  -1     -1    694.99   -1    667.84  558.34];
// Enthalpy calculations:
// Stream in rows and coeficients (nu1 nu2 and nu3) in rows:
//       A1        A2       A3        A4          A5        A6        A7         A8      B1       B2       B3       C1        C2       D1       D2  
coef=[-6.8909 -6.8909   -6.8909   -6.8909   -6.8909   -6.8909    -6.8909   -6.8909   -14.8538  -14.8538 -14.8538  -28.2807  -28.2807 -11.4172 -11.4172;
       0.0991  0.0991   0.0991   0.0991   0.0991    0.0991   0.0991  0.0991    0.1333  0.1333  0.1333  0.1385   0.1385    0.1229 0.1229;
      1.1081e-4 1.1081e-4 1.1081e-4 1.1081e-4 1.1081e-4  1.1081e-4 1.1081e-4  1.1081e-4 7.539e-5 7.539e-5 7.539e-5 9.043e-5 9.043e-5 7.94e-5 7.94e-5];

nflow = length(flow_full);
ntemp = length(temp_full);
//the enthalpies:
temp_exp = [ones(1,ntemp)
            temp_full 
            temp_full.^2];
enthalpy = sum(temp_exp.*coef, 'r');    

xmfull =[flow_full(:); temp_full(:)];    
xm = xmfull;    
// the variance

var = ones(1,nflow + ntemp);
var(:,1:nflow) = (flow_full*0.02).^2;
var(:,nflow + 1: nflow + ntemp) = (0.75).^2;
var = var';
//The variable classification

[At, umeas, fixed] = jac_flowsheet_residuals(flow_full,temp_full, flow, temp, coef);

[red, just_measured, observ, non_obs, spec_cand] = qrlinclass(At,umeas)

// reconcile with all measured to reconcile with only redundant variables, uncomment the "red" assignments
measured = setdiff([1:30], umeas);
// to reconcile with all variables, uncomment bellow
//measured = [1:30];
red = measured;

// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS = 0
// Absolute sum of squares = 1
//Cauchy = 2
//Contamined Normal = 3
//Fair  = 4
//Hampel = 5
//Logistic = 6
//Lorenztian = 7
//Quasi Weighted = 0
// run the configuration functions with the desired objective function type
obj_function_type = 0;
exec ../functions/setup_DR.sce;


// return the problem structure (jacobian, hessian, number of non-zeros, variable type, etc)
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  =  struct_energy_bal(flow_full, temp_full, coef);

params = init_param();
// We use the given Hessian
//params = add_param(params,"hessian_constant","yes");
params = add_param(params,"hessian_approximation","exact");
// uncheck bellow to test derivatives
//params = add_param(params,"derivative_test","second-order");
params = add_param(params,"derivative_test","first-order");
params = add_param(params,"tol",1e-6);
params = add_param(params,"acceptable_tol",1e-6);
//params = add_param(params,"mu_strategy","monotone");
params = add_param(params,"mu_strategy","adaptive");
params = add_param(params,"journal_level",5);
disp('begore start ipopt')
tic
// if the user want to use random initial guess, uncomment 2 lines bellow and comment the 3rd line
//xrand = rand(30,1);
//[x_sol, f_sol, extra] = ipopt(xrand, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

[x_sol, f_sol, extra] = ipopt(xmfull, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

toc

mprintf("\n\nSolution: , x\n");
for i = 1 : nv
    mprintf("x[%d] = %e\n", i, x_sol(i));
end

mprintf("\n\nObjective value at optimal point\n");
mprintf("f(x*) = %e\n", f_sol);
