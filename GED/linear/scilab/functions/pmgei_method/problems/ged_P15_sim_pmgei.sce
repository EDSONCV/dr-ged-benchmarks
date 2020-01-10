// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Steam metering system
//Serth, R W, and W A Heenan. 1986.
// Gross error detection and data reconciliation in steam-metering systems. 
//AIChE Journal 32: 733-747.
//Bibtex Citation

//@article{Serth1986,
//author = {Serth, R W and Heenan, W A},
//journal = {AIChE Journal},
//pages = {733--747},
//title = {{Gross error detection and data reconciliation in steam-metering systems}},
//volume = {32},
//year = {1986}
//}

// 28 Streams
// 11 Equipments 
// the measures
getd('../../');
getd('../../../jacobians/');
getd('../method/');
getd('../method/pls');
cd  '../../'
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt; 

stacksize('max');
tic;
// fixed to make all balances exact
xr =[0.86;1;111.82;109.96;53.27;112.27;2.32;164.05;0.83;52.41;14.86;67.27;111.27;91.86;60;23.64;32.73;16.23;7.85;10.6;87.32;5.45;2.59;46.63;85.46;81.23;70.79;72.22];
//xr =[0.86;1;111.82;109.95;53.27;112.27;2.32;164.05;0.86;52.41;14.86;67.27;111.27;91.86;60;23.64;32.73;16.23;7.95;10.5;87.27;5.45;2.59;46.64;85.45;81.32;70.77;72.73];
szx = size(xr,1);
runsize = 500;
//sd = [0.022
//0.025
//2.796
//2.749
//1.332
//2.807
//0.058
//4.101
//0.021
//1.310
//0.372
//1.682
//2.782
//2.297
//1.500
//0.591
//0.818
//0.406
//0.196
//0.263
//2.183
//0.136
//0.065
//1.166
//2.137
//2.033
//1.770
//1.806];
sd = 2*ones(28,1);
sds = sd;
var=sd.^2;
jac=jacP15();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);


[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 5, 9, 0.07, 0.15);
//[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.02, 0.07);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

//observability/redundancy tests
//user can set unmeasured streams here, if this vector is empty, all streams are measured                  
umeas_P15 = [];
[red_P15, just_measured_P15, observ_P15, non_obs_P15, spec_cand_P15] = qrlinclass(jac,umeas_P15);
measured_P15 = setdiff([1:length(xr)], umeas_P15);
red = measured_P15;//
        
// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS analytical = -1 WLS numerical = 0  ; Absolute sum of squares = 1 ; Cauchy = 2 ;Contamined Normal = 3 ; Fair  = 4
//Hampel = 5 Logistic = 6 ; Lorenztian = 7 ; Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = 2;

[x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);

[res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ] = calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi);


//wls
//[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.0048, 0.125, zadj, runsize, jac_col);
//[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.010, 0.105, jac_row, runsize, zr_nt_nodal);
//alfa_gt_mt = 0.1; alfa_gt_nt = 0.1; alfa_mt1 = 0.0048; alfa_mt2 =0.125; alfa_nt1 = 0.01; alfa_nt2 = 0.105;
////cauchy sigma =1
//[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.08, 0.08, gamaMeasuremts, runsize, rj, jac_col, jac_row);
//[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.00000002, 0.000001, zadj, runsize, jac_col);
//[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.0125, 0.13, jac_row, runsize, zr_nt_nodal)
//alfa_gt_mt = 0.08; alfa_gt_nt = 0.08; alfa_mt1 = 0.00000002; alfa_mt2 =0.000001; alfa_nt1 = 0.0125; alfa_nt2 = 0.13;

//cauchy sigma =2
[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.08, 0.08, gamaMeasuremts, runsize, rj, jac_col, jac_row);
[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.00000006, 0.000002, zadj, runsize, jac_col);
[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.0125, 0.13, jac_row, runsize, zr_nt_nodal)

alfa_gt_mt = 0.08; alfa_gt_nt = 0.08; alfa_mt1 = 0.00000006; alfa_mt2 =0.000002; alfa_nt1 = 0.0125; alfa_nt2 = 0.13;

runtime=toc();

nvalidate = 10; lower_bias = 5; delta_bias = 1; upper_bias = 9; lower_leak = 0.07; delta_leak = 0.02; upper_leak = 0.15; 

is_multi = 0;
pause
[p15_train, p15_validate]  = generate_trainning2(xr, sd, jac, runsize, nvalidate, lower_bias, delta_bias, upper_bias, lower_leak,delta_leak,upper_leak, alfa_gt_mt,alfa_gt_nt,alfa_mt1,alfa_mt1, alfa_nt1, alfa_nt2,obj_function_type, is_multi);
ndatainterval = 5;
[list_models_P15_cauchy , p15_stat ] = generate_pls_models_m( 'P15', 28, 11, p15_train, p15_validate, nvalidate,ndatainterval);
[avti_meas, op_meas, selectivity_meas, aee_meas, avti_eqp, op_eqp, selectivity_eqp, aee_eqp] = get_lit_info(p15_stat, jac_col, jac_row) 
runtime=toc();
//streamNames =generateStreamName(szx);
//prettyprinttable([tokens(streamNames), string([xr, rrn(4,sd), rrn(3,adj), rrn(3,detect), rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)])],"latex")
//eqpNames = generateEqpName('', jac_row);
//prettyprinttable([tokens(eqpNames), string([rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)])],"latex")
//[ op_gt_mt avti_gt_mt avt1_mt1 avt1_mt2 avti_glr avt1_nt1 avt1_nt2  avti_glr_nt runtime ]
//prettyprinttable(string([rrn(3,avt1_mt1),  rrn(3,avt1_mt2),  rrn(3,avti_glr),  rrn(3,avt1_nt1),  rrn(3,avt1_nt2),  rrn(3,avti_glr_nt)]))
//[rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)]
//[rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)]
//
////saving results
//aa = clock();
//nowtime = '_' + string(aa(4)) + '-'+ string(aa(5));
//save ('P_resumed_' + date() + nowtime +'.sav', runtime,  adj, detect, op_nt1, op_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, avti_gt_mt, op_glr_mt, op_glr_nt, avti_glr, avti_glr_nt, aee_nt, aee_mt);
cd 'pmgei_method/problems';
