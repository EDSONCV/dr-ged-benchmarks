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
// 4 Compounds
getd('.');
getd('../functions');
clear tstraoflow_full tstraoflow tstraocomp_full tstraocomp At umeas fixed red lower upper var_lin_type constr_lin_type constr_lhs constr_rhs  just_measured observ non_obs spec_cand x_sol f_sol lower upper extra xmfull ncomp var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

// In the original paper, streams 2 to 12  are unmeasured, 
//theses values are estimates givem by the paper's original author.

tstrao_flow_full = [691.67, 727.54, 699.36, 687.15, 35.87, 12.51, 27.88 , 23.36, 22.67, 4.79, 4.52, 9.31]
//information of measured/unmeasured(-1)/fixed(-5)
tstrao_flow = [-5, -1, -1, -1, -1, -1, -1 , -1, -1, -1, -1, -1] 

tstrao_comp_full = [0.41 0.86 0.32 0.1  9.6 12.4   22   8.1  22.4   23   47.5  34.9;
                    2.58  2.63 2.88 2.94 3.64 4.1  3.52 4.38    4    4.82 2.56   3.5;
                    4     4.38 4.5  4.5  11.8 12.8 12.2  13.2  12.4  14.8 10.2  12.2;
                    93.01   92  92.30  92.46  74.96  70.7 62.28 74.32 66.2  57.38  39.74 49.4]/100                    
//information of measured/unmeasured(-1)/fixed(-5)
tstrao_comp = [0.41 0.86 0.32 0.1  9.6 12.4   22   8.1  22.4   23   47.5  34.9;
              2.58  2.63 2.88 2.94 3.64 4.1  3.52 4.38    4    4.82 2.56   3.5;
              4     4.38 4.5  4.5  11.8 12.8 12.2  13.2  12.4  14.8 10.2  12.2;
              -100   -100     -100   -100    -100   -100  -100    -100    -100    -100   -100    -100]/100   

//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12
jac = [ 1   -1  0   0   1   0   0   0    0   0   0   0 
        0   1   -1  0   0   0   -1  0    0   0   0   0 
        0   0   1   -1  0   -1  0   0    0   0   0   0 
        0   0   0   0   -1  1   0   1    0   0   0   0 
        0   0   0   0   0   0   0   -1   1   0   0  -1  
        0   0   0   0   0   0   1   0    -1  1   0   0
        0   0   0   0   0   0   0   0    0   -1  -1  1 ];                                
//      1   2   3   4   5   6   7   8    9   10  11  12

// organizing the vector for the constraints residuals
xmfull=[tstrao_flow_full(:);matrix(tstrao_comp_full',-1)];
xm=xmfull;
//the variance proposed by the original author
sd = (0.01*xmfull);
//recalculating variance
for i=1: length(sd)
   if sd(i) <= 0.0001 
       sd(i) = 0.0001;
   end       
end
var = sd.^2;
ncomp=4;

//observability/redundancy tests
[At,umeas, fixed] = jac_compound_residuals(jac,ncomp,tstrao_flow_full,tstrao_comp_full, tstrao_flow, tstrao_comp);
[red, just_measured, observ, non_obs, spec_cand] = qrlinclass(At,umeas)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured = setdiff([1:length(xmfull)], umeas);
// to reconcile with all variables, comment the line above and uncomment bellow
//measured = [1:length(xmfull)];
red=measured;


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
obj_function_type = 5;
exec ../functions/setup_DR.sce;
// to run robust reconciliation, it is also necessary to choose the function to return the problem structure

// ipopt needs some information about the problem, such as jacobian and hessian structure, 
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = structure_compound(jac,ncomp, tstrao_flow_full,tstrao_comp_full);

params = init_param();
// We use the given Hessian
params = add_param(params,"hessian_approximation","exact");

// notice that the option bellow must only be set when the Hessians are constant, 
// eg: weighted least squares objective function and linear or bilinear Jacobian
//params = add_param(params,"hessian_constant","yes");
params = add_param(params,"derivative_test","first-order");
//params = add_param(params,"derivative_test","second-order");
//params = add_param(params,"derivative_test_print_all","yes");
params = add_param(params,"tol",1e-4);
params = add_param(params,"acceptable_tol",1e-4);
params = add_param(params,"mu_strategy","monotone");
params = add_param(params,"journal_level",5);
params = add_param(params,"fixed_variable_treatment", "relax_bounds");

disp('begore start ipopt')
//according to the original paper, we fix the measured total flow
lower(fixed) = xmfull(fixed);
upper(fixed) = xmfull(fixed);
tic
//xrd=rand(size(xmfull,1),size(xmfull,2));
//[x_sol, f_sol, extra] = ipopt(xrd, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);
[x_sol, f_sol, extra] = ipopt(xmfull, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);
toc

mprintf("\n\nSolution: , x\n");
for i = 1 : nv
    mprintf("x[%d] = %e\n", i, x_sol(i));
end

mprintf("\n\nObjective value at optimal point\n");
mprintf("f(x*) = %e\n", f_sol);
//printing results
[Aeqp, Astreams] =size(jac)
xx=matrix(x_sol,Astreams,ncomp+1)
TotalFlowMeasured = xx(:,1)'
compoundMeasured  = 100*xx(:, 2:$)
