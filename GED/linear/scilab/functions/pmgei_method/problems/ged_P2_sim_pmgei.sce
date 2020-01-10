// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Authors
//Dovì, V G, and C Solisio. 2001. Reconciliation of censored measurements in chemical processes: an alternative approach. Chemical Engineering Journal 84, no. 3 (December): 309-314. http://www.sciencedirect.com/science/article/B6TFJ-45KNHB1-F/2/199f358469628f600f10b394d2b55a8b.

//Bibtex Citation

//@article{Dovi2001,
//annote = { The importance of considering the censoring of measured data in the reconciliation of process flow rates has been shown in a previous paper [Chem. Eng. Sci. 52 (17) (1997) 3047]. The purpose of the present paper is to introduce a new technique for carrying out the actual reconciliation procedure and compare its significance and performance with those of previous methods. A numerical example shows how nontrivial differences are to be expected.},
//author = {Dov\`{\i}, V G and Solisio, C},
//isbn = {1385-8947},
//journal = {Chemical Engineering Journal},
//keywords = {Censored data,Data reconciliation,Detection limits},
//month = dec,
//number = {3},
//pages = {309--314},
//title = {{Reconciliation of censored measurements in chemical processes: an alternative approach}},
//url = {http://www.sciencedirect.com/science/article/B6TFJ-45KNHB1-F/2/199f358469628f600f10b394d2b55a8b},
//volume = {84},
//year = {2001}
//}
// 6 Streams
// 3 Equipments 
getd('../../');
getd('../../../jacobians/');
getd('../method/');
getd('../method/pls');
cd  '../../'
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt; 

//stacksize('max');
tic;
xr =[11;10;21;11;5.5;5.5];

szx = size(xr,1);
runsize = 500;
//the variance proposed by this work
//sd =[0.032
//0.026
//0.120
//0.033
//0.052
//0.015].^(0.5);
sd=ones(6,1);
sds = sd;
var=sd.^2;
jac=jacP2();
rj=rank(jac);
jac_col = size(jac,2);
jac_row = size(jac,1);
sigma=diag(sds.^2);

[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 5, 9, 0.07, 0.15);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

//observability/redundancy tests
//user can set unmeasured streams here, if this vector is empty, all streams are measured                  
umeas_P2 = [];
[red_P2, just_measured_P2, observ_P2, non_obs_P2, spec_cand_P2] = qrlinclass(jac,umeas_P2);
measured_P2 = setdiff([1:length(xr)], umeas_P2);
red = measured_P2;//
        
// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS analytical = -1 WLS numerical = 0  ; Absolute sum of squares = 1 ; Cauchy = 2 ;Contamined Normal = 3 ; Fair  = 4
//Hampel = 5 Logistic = 6 ; Lorenztian = 7 ; Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = 2;

[x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);

[res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ] = calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi);

[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.12, 0.12, gamaMeasuremts, runsize, rj, jac_col, jac_row);

[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.005, 0.029, zadj, runsize, jac_col);

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.037, 0.107, jac_row, runsize, zr_nt_nodal);

pause

nvalidate = 10; lower_bias = 5; delta_bias = 1; upper_bias = 9; lower_leak = 0.07; delta_leak = 0.02; upper_leak = 0.15; 
//cauchy 0.037, 0.107
alfa_gt_mt = 0.12; alfa_gt_nt = 0.12; alfa_mt1 = 0.005; alfa_mt2 =  0.029; alfa_nt1 = 0.037; alfa_nt2 = 0.107;
is_multiple = 0;
//logistic
//alfa_gt_mt = 0.1; alfa_gt_nt = 0.1; alfa_mt1 = 0.01; alfa_mt2 =  0.06; alfa_nt1 = 0.037; alfa_nt2 = 0.107;
//lorenz
//alfa_gt_mt = 0.1; alfa_gt_nt = 0.1; alfa_mt1 = 0.023; alfa_mt2 =  0.128; alfa_nt1 = 0.037; alfa_nt2 = 0.107;

[p2_train, p2_validate]  = generate_trainning2(xr, sd, jac, runsize, nvalidate, lower_bias, delta_bias, upper_bias, lower_leak,delta_leak,upper_leak, alfa_gt_mt,alfa_gt_nt,alfa_mt1,alfa_mt1, alfa_nt1, alfa_nt2,obj_function_type,0);
ndatainterval = 5
[list_models_P2,p2_stat] = generate_pls_models_m( 'P2', 6, 3, p2_train, p2_validate, nvalidate,ndatainterval);

[avti_meas, op_meas, selectivity_meas, aee_meas, avti_eqp, op_eqp, selectivity_eqp, aee_eqp] = get_lit_info(p2_stat, jac_col, jac_row);
cd 'pmgei_method/problems';
