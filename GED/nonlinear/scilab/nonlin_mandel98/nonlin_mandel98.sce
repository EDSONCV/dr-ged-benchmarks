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
// 2 Compounds
getd('.');
getd('../functions');
clear flow_full_P8 flow_P8 comp_full_P8 comp_P8 At_P8 umeas_P8 fixed_P8 red_P8 lower_P8 upper_P8 var_lin_type_P8 constr_lin_type_P8 constr_lhs_P8 constr_rhs_P8  just_measured_P8 observ_P8 non_obs_P8 spec_cand_P8 x_sol f_sol lower_P8 upper_P8 extra xmfull ncomp var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs

     //               1  2   3   4   5   6   7   8   9  10  11  12  
flow_full_mandel = [230	21	209	35	174	15	159	50	209	94	115	44];

//information of measured/unmeasured(-1)/fixed(-5)
     //               1  2   3   4   5   6   7   8   9  10  11  12  
flow_mandel =     [230	21	-1	35	-1	15	159	-1   -1 -1	115	44];
         //                1      2      3         4        5     6        7       8      9      10      11      12      
comp_full_mandel = [  0.728	    0.708	0.730	0.739	0.728	0.846	0.717	0.541	0.675	0.949	0.452	1.412
                      99.272	99.292	99.270	99.261	99.272	99.154	99.283	99.459	99.325	99.051	99.548	98.588]/100;


//information of measured/unmeasured(-1)/fixed(-5)
 //                1      2      3         4        5     6        7       8      9      10      11      12      
comp_mandel = [  0.0046	0.0045	0.0046	0.0050	-1  	0.0065	-1 	0.0020	0.0038	-1	0.0008	0.0137
                 -1        -1     -1    -1       -1       -1     -1       -1       -1    -1      -1      -1];
//The jacobian of the constraints
jac = [ 1  -1  -1   0   0   0   0   0    0   0   0   0   
        0   0   1   -1  -1  0   0    0   0   0   0   0   
        0   0   0   0   1   -1  -1  0    0   0   0   0 
        0   0   0   0   0   0   1    1    -1  0   0   0
        0   0   0   0   0   0   0   -1    0  1   0   -1
        0   0  0   0   0   0   0   0    1   -1  -1  0];

// organizing the vector for the constraints residuals
xmfull=[flow_full_mandel(:);matrix(comp_full_mandel',-1)];
xm=xmfull;

//the variance proposed by the original author
// in original paper the standard deviation is given. so it must be squared.
//   1     2   3    4     5     6    7   8   9    10   11   12  
sd=[0.15 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05  0.2 0.05;
    0.05 0.02 0.05 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02;
    0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02;]';
sd = matrix(sd', -1);    
// to run with equaly weighted standard deviation, uncomment the line below
//sd = ones(size(xmfull,1),size(xmfull,2));

var = sd.^2;
ncomp = 2 ;
//observability/redundancy tests                  
[At_mandel,umeas_mandel, fixed_mandel] = jac_compound_residuals(jac,ncomp,flow_full_mandel,comp_full_mandel, flow_mandel, comp_mandel);
[red_mandel, just_measured_mandel, observ_mandel, non_obs_mandel, spec_cand_mandel] = qrlinclass(At_mandel,umeas_mandel)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_mandel = setdiff([1:length(xmfull)], umeas_mandel);
red = measured_mandel;//
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
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower_mandel, upper_mandel, var_lin_type_mandel, constr_lin_type_mandel, constr_lhs_mandel, constr_rhs_mandel]  = structure_compound(jac,ncomp, flow_full_mandel,comp_full_mandel);
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
lower_mandel(fixed_mandel) = xmfull(fixed_mandel);
upper_mandel(fixed_mandel) = xmfull(fixed_mandel);
xm_init = xmfull./var;
tic
[x_sol, f_sol, extra] = ipopt(xmfull, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type_mandel, constr_lin_type_mandel, constr_rhs_mandel, constr_lhs_mandel, lower_mandel, upper_mandel, params);
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

