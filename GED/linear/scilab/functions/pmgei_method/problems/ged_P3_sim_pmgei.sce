// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Atmospheric tower example from:
// Zhang, P, G Rong, and Y Wang. 2001.
// A new method of redundancy analysis in data reconciliation and its application.
// Computers and Chemical Engineering 25: 941-949.

//Bibtex Citation

//@article{Zhang2001,
//author = {Zhang, P and Rong, G and Wang, Y},
//journal = {Computers and Chemical Engineering},
//pages = {941--949},
//title = {{A new method of redundancy analysis in data reconciliation and its application}},
//volume = {25},
//year = {2001}
//}

//12 Streams 
//3 Equipments
getd('../../');
getd('../../../jacobians/');
getd('../method/');
getd('../method/pls');
cd  '../../'
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt avti_meas op_meas selectivity_meas aee_meas avti_eqp op_eqp selectivity_eqp aee_eqp list_models_P3 p3_statp3_train p3_validate  alfa_gt_mt alfa_gt_nt alfa_mt1 alfa_mt2 alfa_nt1 alfa_nt2 nvalidate  lower_bias delta_bias upper_bias lower_leak delta_leak  upper_leak ;
is_multiple = 0; ; 
//pause
//stacksize('max');
tic;
xr =[189.98000
174.60000
3.13900
32.77000
33.47000
7.25000
0.31600
92.37600
28.62900
23.80000
18.52600
55.56800
];

szx = size(xr,1);
runsize = 400;
//sd = [2.14120
//1.94840
//0.03410
//0.34510
//0.40010
//0.08680
//0.00356
//1.05940
//0.36070
//0.30010
//0.20011
//0.64530
//];
sd=ones(12,1);
//sd=5*ones(12,1);
sds = sd;
var=sd.^2;
jac=jacP3();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);


[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 5, 9, 0.07, 0.15);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

//observability/redundancy tests
//user can set unmeasured streams here, if this vector is empty, all streams are measured                  
umeas_P3 = [];
[red_P3, just_measured_P3, observ_P3, non_obs_P3, spec_cand_P3] = qrlinclass(jac,umeas_P3);
measured_P3 = setdiff([1:length(xr)], umeas_P3);
red = measured_P3;//
        
// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS analytical = -1 WLS numerical = 0  ; Absolute sum of squares = 1 ; Cauchy = 2 ;Contamined Normal = 3 ; Fair  = 4
//Hampel = 5 Logistic = 6 ; Lorenztian = 7 ; Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = 2;

[x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);

[res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ] = calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi);

[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.105, 0.105, gamaMeasuremts, runsize, rj, jac_col, jac_row);
// cauchy
[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.000015, 0.0002, zadj, runsize, jac_col)

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.038, 0.11, jac_row, runsize, zr_nt_nodal);

nvalidate = 5; lower_bias = 5; delta_bias = 1; upper_bias = 9; lower_leak = 0.07; delta_leak = 0.02; upper_leak = 0.15; 
//cauchy 
//alfa_gt_mt = 0.1; alfa_gt_nt = 0.1; alfa_mt1 = 0.000025; alfa_mt2 =0.00031; alfa_nt1 = 0.038; alfa_nt2 = 0.11;
//sigma = unit*5
//[avt1_mt1, avt1_mt2, op_mt1, op_mt2] =  measurement_test(0.000002, 0.00002, zadj, runsize, jac_col)
//pause
alfa_gt_mt = 0.105; alfa_gt_nt = 0.105; alfa_mt1 = 0.000015; alfa_mt2 =0.0002; alfa_nt1 = 0.038; alfa_nt2 = 0.11;
is_multiple = 0;

[p3_train, p3_validate]  = generate_trainning2(xr, sd, jac, runsize, nvalidate, lower_bias, delta_bias, upper_bias, lower_leak,delta_leak,upper_leak, alfa_gt_mt,alfa_gt_nt,alfa_mt1,alfa_mt1, alfa_nt1, alfa_nt2,obj_function_type,is_multiple);

ndatainterval = 5

[list_models_P3, p3_stat] = generate_pls_models_m( 'P3', 12, 3, p3_train, p3_validate, nvalidate,ndatainterval);
[avti_meas, op_meas, selectivity_meas, aee_meas, avti_eqp, op_eqp, selectivity_eqp, aee_eqp] = get_lit_info(p3_stat, jac_col, jac_row)

runtime=toc();
streamNames =generateStreamName(szx);
prettyprinttable([tokens(streamNames), string([xr, rrn(4,sd), rrn(3,adj), rrn(3,detect), rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)])],"latex")
eqpNames = generateEqpName('', jac_row);
prettyprinttable([tokens(eqpNames), string([rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)])],"latex")
[ op_gt_mt avti_gt_mt avt1_mt1 avt1_mt2 avti_glr avt1_nt1 avt1_nt2  avti_glr_nt runtime ]
prettyprinttable(string([rrn(3,avt1_mt1),  rrn(3,avt1_mt2),  rrn(3,avti_glr),  rrn(3,avt1_nt1),  rrn(3,avt1_nt2),  rrn(7,avti_glr_nt)]))
[rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)]
[rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)]

//saving results
aa = clock();
nowtime = '_' + string(aa(4)) + '-'+ string(aa(5));
save ('P_resumed_' + date() + nowtime +'.sav', runtime,  adj, detect, op_nt1, op_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, avti_gt_mt, op_glr_mt, op_glr_nt, avti_glr, avti_glr_nt, aee_nt, aee_mt);

cd 'pmgei_method/problems';
