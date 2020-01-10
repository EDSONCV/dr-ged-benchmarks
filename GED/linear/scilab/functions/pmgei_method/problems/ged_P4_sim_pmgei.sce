// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Heat exchanger with by-pass valve
//Narasimhan, S, and C Jordache. 2000.
//Data Reconciliation and Gross Error Detection: An Intelligent Use of Process Data. 1st ed.
//Houston: Gulf Publishing.

//Bibtex Citation

//@book{Narasimhan2000,
//address = {Houston},
//author = {Narasimhan, S and Jordache, C},
//booktitle = {Process Data. Gulf Professional Publishing, Houston, TX.},
//edition = {1},
//publisher = {Gulf Publishing},
//title = {{Data Reconciliation and Gross Error Detection: An Intelligent Use of Process Data}},
//year = {2000}
//}

// 6 Streams
// 4 Equipments getd('../functions/');
getd('../../');
getd('../../../jacobians/');
getd('../method/');
getd('../method/pls');
cd  '../../'
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt; 

//stacksize('max');
tic;
xr=[100;64;36;64;36;100];
szx = size(xr,1);
runsize = 500;
sd = ones(6,1);
sds = ones(6,1);
var=sd.^2;
jac=jacP4();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);

[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize,5, 9, 0.07, 0.15);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

//observability/redundancy tests
//user can set unmeasured streams here, if this vector is empty, all streams are measured                  
umeas_P4 = [];
[red_P4, just_measured_P4, observ_P4, non_obs_P4, spec_cand_P4] = qrlinclass(jac,umeas_P4);
measured_P4 = setdiff([1:length(xr)], umeas_P4);
red = measured_P4;//
        
// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS analytical = -1 WLS numerical = 0  ; Absolute sum of squares = 1 ; Cauchy = 2 ;Contamined Normal = 3 ; Fair  = 4
//Hampel = 5 Logistic = 6 ; Lorenztian = 7 ; Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = 2;

[x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);

[res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ] = calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi);

[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.08, 0.08, gamaMeasuremts, runsize, rj, jac_col, jac_row);

[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.0033, 0.02, zadj, runsize, jac_col)

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.023, 0.088, jac_row, runsize, zr_nt_nodal);


nvalidate = 10; lower_bias = 5; delta_bias = 1; upper_bias = 9; lower_leak = 0.07; delta_leak = 0.02; upper_leak = 0.15; 
alfa_gt_mt = 0.08; alfa_gt_nt = 0.08; alfa_mt1 = 0.0033; alfa_mt2 =  0.02; alfa_nt1 = 0.023; alfa_nt2 = 0.088;
pause
is_multiple = 0;

[p4_train, p4_validate]  = generate_trainning2(xr, sd, jac, runsize, nvalidate, lower_bias, delta_bias, upper_bias, lower_leak,delta_leak,upper_leak, alfa_gt_mt,alfa_gt_nt,alfa_mt1,alfa_mt1, alfa_nt1, alfa_nt2,obj_function_type, is_multiple);
ndatainterval = 5
[list_models_P4,p4_stat] = generate_pls_models_m( 'P4', 6, 4, p4_train, p4_validate, nvalidate,ndatainterval);
[avti_meas, op_meas, selectivity_meas, aee_meas, avti_eqp, op_eqp, selectivity_eqp, aee_eqp] = get_lit_info(p4_stat, jac_col, jac_row)


runtime=toc();
//streamNames =generateStreamName(szx);
//prettyprinttable([tokens(streamNames), string([xr, rrn(4,sd), rrn(3,adj), rrn(3,detect), rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)])],"latex")
//eqpNames = generateEqpName('', jac_row);
//prettyprinttable([tokens(eqpNames), string([rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)])],"latex")
//[ op_gt_mt avti_gt_mt avt1_mt1 avt1_mt2 avti_glr avt1_nt1 avt1_nt2  avti_glr_nt runtime ]
//prettyprinttable(string([rrn(3,avt1_mt1),  rrn(3,avt1_mt2),  rrn(3,avti_glr),  rrn(3,avt1_nt1),  rrn(3,avt1_nt2),  rrn(7,avti_glr_nt)]))
//[rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)]
//[rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)]

cd 'pmgei_method/problems';
