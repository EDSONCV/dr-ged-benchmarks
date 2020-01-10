// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Yang, Youqi, Rongbo Ten, and Luiqun Jao. 1995. 
//“A study of gross error detection and data reconciliation in process industries.” Comp. & Chem. 
//Eng 19S:S217-S222.

//Bibtex Citation

//@article{Yang1995,
//author = {Yang, Youqi and Ten, Rongbo and Jao, Luiqun},
//journal = {Comp. \& Chem. Eng},
//keywords = {combinatory approach,data reconciliation,gross error detection},
//pages = {S217--S222},
//title = {{A study of gross error detection and data reconciliation in process industries}},
//volume = {19S},
//year = {1995}
//}
// 8 Streams
// 4 Equipments 
getd('../functions/');
getd('../jacobians/');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt; 

stacksize('max');
tic;
xr =[98.7;41.1;78.9;30.2;109.1;19.8;57.6;37.8];

szx = size(xr,1);
runsize = 2500;
sd = [0.9870
0.4110
0.7890
0.3020
1.0910
0.1980
0.5760
0.3780].^(0.5);
sds = sd;
var=sd.^2;
jac=jacP5();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);
sigma_inv=inv(sigma);



[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 5, 9, 0.07, 0.15);
//[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.02, 0.07);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

//observability/redundancy tests
//user can set unmeasured streams here, if this vector is empty, all streams are measured                  
umeas_P5 = [];
[red_P5, just_measured_P5, observ_P5, non_obs_P5, spec_cand_P5] = qrlinclass(jac,umeas_P5);
measured_P5 = setdiff([1:length(xr)], umeas_P5);
red = measured_P5;//
        
// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS analytical = -1 WLS numerical = 0  ; Absolute sum of squares = 1 ; Cauchy = 2 ;Contamined Normal = 3 ; Fair  = 4
//Hampel = 5 Logistic = 6 ; Lorenztian = 7 ; Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = -1;

[x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);

[res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ] = calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi);

[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.11, 0.11, gamaMeasuremts, runsize, rj, jac_col, jac_row);

[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.020, 0.15, zadj, runsize, jac_col);


[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.028, 0.107, jac_row, runsize, zr_nt_nodal);

[avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.15, 0.22, runsize);

//[ avt1_mt1 avt1_mt2 avt1_nt1 avt1_nt2   avti_glr avti_glr_nt  avti_gt_mt avti_gt_nt]
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
//aa = clock();
//nowtime = '_' + string(aa(4)) + '-'+ string(aa(5));
//save ('P_resumed_Logistic' + date() + nowtime +'.sav', runtime,  adj, detect, op_nt1, op_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, avti_gt_mt, op_glr_mt, op_glr_nt, avti_glr, avti_glr_nt, aee_nt, aee_mt);

