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
getd('../../');
getd('../../../jacobians/');
getd('../method/');
getd('../method/pls');
cd  '../../'
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt; 

stacksize('max');
tic;

xr =[690;725;700;685;35;15;25;20;30;5;5;10];
szx = size(xr,1);
runsize = 500;
//the variance proposed by the original author
//sd = [20.7501
//21.8262
//20.9808
//20.6145
//1.0761
//0.3753
//0.8364
//0.7008
//0.6801
//0.1437
//0.1356
//0.2793];
sd=0.8*ones(12,1);
sds = sd;
var=sd.^2;
jac=jacP11();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);


[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
//[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.1, 0.2);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 5, 9, 0.07, 0.15);
resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

//observability/redundancy tests
//user can set unmeasured streams here, if this vector is empty, all streams are measured                  
umeas_P11 = [];
[red_P11, just_measured_P11, observ_P11, non_obs_P11, spec_cand_P11] = qrlinclass(jac,umeas_P11);
measured_P11 = setdiff([1:length(xr)], umeas_P11);
red = measured_P11;//
        
// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS analytical = -1 WLS numerical = 0  ; Absolute sum of squares = 1 ; Cauchy = 2 ;Contamined Normal = 3 ; Fair  = 4
//Hampel = 5 Logistic = 6 ; Lorenztian = 7 ; Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = 2;

[x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);
[res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ] = calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi);

// Cauchy sd = 1
//
//[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.12, 0.12, gamaMeasuremts, runsize, rj, jac_col, jac_row);
//[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.00033, 0.0039, zadj, runsize, jac_col);
//[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.017, 0.115, jac_row, runsize, zr_nt_nodal);
//
//alfa_gt_mt = 0.12; alfa_gt_nt = 0.12; alfa_mt1 = 0.00033; alfa_mt2 =0.0039; alfa_nt1 = 0.017; alfa_nt2 = 0.115;


//cauchy sd =0.8


[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.1, 0.1, gamaMeasuremts, runsize, rj, jac_col, jac_row);
[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.00035, 0.0045, zadj, runsize, jac_col)
[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.017, 0.113, jac_row, runsize, zr_nt_nodal)

alfa_gt_mt = 0.1; alfa_gt_nt = 0.1; alfa_mt1 = 0.00035; alfa_mt2 =0.0045; alfa_nt1 = 0.017; alfa_nt2 = 0.113;


nvalidate = 10; lower_bias = 5; delta_bias = 1; upper_bias = 9; lower_leak = 0.07; delta_leak = 0.02; upper_leak = 0.15; 
is_multiple = 0;
pause
clear res gamaMeasuremts gamaNodal zr_nt_nodal zr_nt_nodal_rand zadj x_sol resGrossErrorNodalRandFi;
[p11_train, p11_validate]  = generate_trainning2(xr, sd, jac, runsize, nvalidate, lower_bias, delta_bias, upper_bias, lower_leak,delta_leak,upper_leak, alfa_gt_mt,alfa_gt_nt,alfa_mt1,alfa_mt1, alfa_nt1, alfa_nt2,obj_function_type, is_multiple);
ndatainterval = 5;
[list_models_P11, p11_stat] = generate_pls_models_m( 'P11', 12, 7, p11_train, p11_validate, nvalidate,ndatainterval);
[avti_meas, op_meas, selectivity_meas, aee_meas, avti_eqp, op_eqp, selectivity_eqp, aee_eqp] = get_lit_info(p11_stat, jac_col, jac_row)

runtime=toc();




//streamNames =generateStreamName(szx);
//prettyprinttable([tokens(streamNames), string([xr, rrn(4,sd), rrn(3,adj), rrn(3,detect), rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(1,aee_mt)])],"latex")
//eqpNames = generateEqpName('', jac_row);
//prettyprinttable([tokens(eqpNames), string([rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(1,aee_nt)])],"latex")
//
//prettyprinttable([tokens(streamNames), string([xr, rrn(4,sd), rrn(3,adj), rrn(3,detect), rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(4,aee_mt)])],"latex")
//eqpNames = generateEqpName('', jac_row);
//prettyprinttable([tokens(eqpNames), string([rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(4,aee_nt)])],"latex")
//
//[ op_gt_mt avti_gt_mt avt1_mt1 avt1_mt2 avti_glr avt1_nt1 avt1_nt2  avti_glr_nt runtime ]
//prettyprinttable(string([rrn(3,avt1_mt1),  rrn(3,avt1_mt2),  rrn(3,avti_glr),  rrn(3,avt1_nt1),  rrn(3,avt1_nt2),  rrn(3,avti_glr_nt)]))
//[rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(4,aee_mt)]
//[rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(4,aee_nt)]
//
//
cd 'pmgei_method/problems';
