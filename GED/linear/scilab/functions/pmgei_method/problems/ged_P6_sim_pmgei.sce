// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Rosenberg, J and Mah, R S H and Iordache, C
//Evaluation of Schemes for Detecting and Identifying Gross Errors in Process Data
//Ind. & Eng. Chem. Proc. Des. Dev, Vol., V. 26:555--564

//Bibtex Citation
//@article{Rosenberg1987a, 
//author = {Rosenberg, J and Mah, R S H and Iordache, C},
//journal = {Ind. \& Eng. Chem. Proc. Des. Dev, Vol.},
//pages = {555--564},
//title = {{Evaluation of Schemes for Detecting and Identifying Gross Errors in Process Data}},
//volume = {26},
//year = {1987}
//}


//7 Streams 
//4 Equipments
getd('../../');
getd('../../../jacobians/');
getd('../method/');
getd('../method/pls');
cd  '../../'
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt; 

stacksize('max');
//stacksize(149900000)
tic;

// the real values
xr =[5 15 15 5 10 5 5 ]';

szx = size(xr,1);
runsize = 500;
//the variance
sd=[1
1
1
1
1
1
1].^0.5;
//sd = [0.13
//0.38
//0.38
//0.13
//0.25
//0.13
//0.13];
sds = sd;
sds =sd;
var=sd.^2;
jac=jacP6();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);


[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);

//[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.25, 0.35);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 5, 9, 0.07, 0.15);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

//observability/redundancy tests
//user can set unmeasured streams here, if this vector is empty, all streams are measured                  
umeas_P6 = [];
[red_P6, just_measured_P6, observ_P6, non_obs_P6, spec_cand_P6] = qrlinclass(jac,umeas_P6);
measured_P6 = setdiff([1:length(xr)], umeas_P6);
red = measured_P6;//
        
// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS analytical = -1 WLS numerical = 0  ; Absolute sum of squares = 1 ; Cauchy = 2 ;Contamined Normal = 3 ; Fair  = 4
//Hampel = 5 Logistic = 6 ; Lorenztian = 7 ; Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = 2;

[x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);

[res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ] = calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi);

[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.105, 0.105, gamaMeasuremts, runsize, rj, jac_col, jac_row);

[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.0017, 0.011, zadj, runsize, jac_col);

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.033, 0.125, jac_row, runsize, zr_nt_nodal);
pause
runtime=toc();

//[p9train, p9validate]  = generate_trainning(xr, sd, jac, runsize, 2, 3, 1, 7, 0.07, 0.02, 0.15, 0.1, 0.1, 0.134, 0.15, 0.021, 0.116);
// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS analytical = -1 WLS numerical = 0  ; Absolute sum of squares = 1 ; Cauchy = 2 ;Contamined Normal = 3 ; Fair  = 4
//Hampel = 5 Logistic = 6 ; Lorenztian = 7 ; Quasi Weighted = 8
// run the configuration functions with the desired objective function type

//nvalidate = 10; lower_bias = 5; delta_bias = 1; upper_bias = 9; lower_leak = 0.25; delta_leak = 0.2; upper_leak = 0.35; 
nvalidate = 10; lower_bias = 5; delta_bias = 1; upper_bias = 9; lower_leak = 0.07; delta_leak = 0.02; upper_leak = 0.15; 
alfa_gt_mt = 0.105; alfa_gt_nt = 0.105; alfa_mt1 = 0.0017; alfa_mt2 =  0.011; alfa_nt1 = 0.033; alfa_nt2 = 0.125;
//cauchy
//[avt1_mt1_cauchy, avt1_mt2_cauchy, op_mt1_cauchy, op_mt2_cauchy] = measurement_test(0.00084, 0.006, zadj, runsize, jac_col);

is_multiple = 0;
clear res gamaMeasuremts gamaNodal zr_nt_nodal zr_nt_nodal_rand zadj x_sol resGrossErrorNodalRandFi;
[p6_train, p6_validate,xfinal_train,x_sol_train]  = generate_trainning2(xr, sd, jac, runsize, nvalidate, lower_bias, delta_bias, upper_bias, lower_leak,delta_leak,upper_leak, alfa_gt_mt,alfa_gt_nt,alfa_mt1,alfa_mt1, alfa_nt1, alfa_nt2,obj_function_type, is_multiple);

ndatainterval = 5

[list_models_P6, p6_stat] = generate_pls_models_m( 'P6', 7, 4, p6_train, p6_validate, nvalidate,ndatainterval);

[avti_meas, op_meas, selectivity_meas, aee_meas, avti_eqp, op_eqp, selectivity_eqp, aee_eqp] = get_lit_info(p6_stat, jac_col, jac_row)  
cd 'pmgei_method/problems';
