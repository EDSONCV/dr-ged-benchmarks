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
// 18 Streams
// 9 Equipments 
getd('../../');
getd('../../../jacobians/');
getd('../method/');
getd('../method/pls');
cd  '../../'
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt; 

stacksize('max');
//stacksize(100000000)
tic;
// changed from 3098.20 to 3098.10 os equipment DA 408, stream S408
// from // DA 402 stream 452 349,07 to 349,04
// from M1 S418 3,84 to 3.80
// from M1 S418 3,84 to 3.83
// from DA401 2872,50 to 2872,47
// from  DA 301 from 225.4 to 225.37
// from DA 301 from 2276,90 to 2276,87

 xr=[582.60; 501.44; 392.10; 33.80; 535.24; 244.46; 147.64; 31.13; 504.11; 7.00; 497.11; 233.72; 247.55; 15.84; 13.37; 2.47; 308.47;310.94]
szx = size(xr,1);
runsize = 350;
//sd =[2.06;1.31;11.12;11.16;18.61;1.22;7.81;4.74;3.32;3.32;3.32;0.85;2.26;1.34;6.35;6.35;6.24;14.97;14.97;5.28;5.28;9.27;23.76;0.03;25.96;25.96;25.96;25.96;25.87;11.46;12.2;0.1;2.98;4.37;3.34;0.29;1;2.15;1.2;1;4.43;0.06;4.18;1.93;2.12;0.12;0.12;0.02;2.86;2.74;];
sd = 3*ones(18,1);
sds = sd;
var=sd.^2;
jac=jacP16_split2_part3();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);

pause
[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
//original
//[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.02, 0.07);
//modified
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 5, 9, 0.07, 0.15);
//
resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];
//
//observability/redundancy tests
//user can set unmeasured streams here, if this vector is empty, all streams are measured                  
umeas_P16 = [];
[red_P16, just_measured_P16, observ_P16, non_obs_P16, spec_cand_P16] = qrlinclass(jac,umeas_P16);
measured_P16 = setdiff([1:length(xr)], umeas_P16);
red = measured_P16;//
        
// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS analytical = -1 WLS numerical = 0  ; Absolute sum of squares = 1 ; Cauchy = 2 ;Contamined Normal = 3 ; Fair  = 4
//Hampel = 5 Logistic = 6 ; Lorenztian = 7 ; Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = 2;

[x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);
//
[res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ] = calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi);
//
//
//[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.000001, 0.000001, gamaMeasuremts, runsize, rj, jac_col, jac_row);
//
////[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.05, 0.13, zadj, runsize, jac_col);
//// cauchy
////[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.000000000000001, 0.000000000000045, zadj, runsize, jac_col)
////[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.0000000025, 0.00000006, jac_row, runsize, zr_nt_nodal);
////



////// cauchy sigma = 1 P_16
//[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.11, 0.11, gamaMeasuremts, runsize, rj, jac_col, jac_row)
//[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.6E-5, 1E-4, zadj, runsize, jac_col)
//[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(1.3E-2, 0.99E-1, jac_row, runsize, zr_nt_nodal)
//
////[avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.135, 0.21, runsize);
////pause
//alfa_gt_mt = 0.11; alfa_gt_nt =0.11; alfa_mt1 = 0.6E-5; alfa_mt2 =1E-4; alfa_nt1 = 1.3E-2; alfa_nt2 =0.99E-1;
//nvalidate = 20; lower_bias = 5; delta_bias = 1; upper_bias = 9; lower_leak = 0.07; delta_leak = 0.02; upper_leak = 0.15; 



//// cauchy sigma = 3 P_16
[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.14, 0.14, gamaMeasuremts, runsize, rj, jac_col, jac_row)
[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(3E-5, 5E-4, zadj, runsize, jac_col)
[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(1.7E-2, 0.13, jac_row, runsize, zr_nt_nodal)
//
////[avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.135, 0.21, runsize);
//pause
alfa_gt_mt = 0.14; alfa_gt_nt =0.14; alfa_mt1 = 3E-5; alfa_mt2 =5E-4; alfa_nt1 = 1.7E-2; alfa_nt2 =0.13;
nvalidate = 10; lower_bias = 5; delta_bias = 1; upper_bias = 9; lower_leak = 0.07; delta_leak = 0.02; upper_leak = 0.15; 
 pause



//streamNames = 'S319 S316 S312 S378 S336 S357 S346 S359P1 S347 S352 S356 S358 S357P S359P2 S359 S338P S338 S341P S341 S414 S502 S411 S401 S415 S402 S404 S405 S407 S408 S453 S460 S456 S452 S511 S503 S384P S52P S592 S581 S525 S524 S536 S527 S549 S550 S537 S598 S599 S267 S538';
//prettyprinttable([tokens(streamNames), string([xr, rrn(4,sd), rrn(3,adj), rrn(3,detect), rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(4,aee_mt)])],"latex")
//eqpNames = "DA-301 Split1 COMP1 EA-32X FA-309 Split2 Mixer1 FSPL H338 H341 DA-401 PTC M1 448-450 DC-401 452-448 DA-408 DA-402 DA-404 DA-405 DC-402 DC-40X DA-407 DA-406 DA-409 Mix2";
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
cd 'pmgei_method/problem
