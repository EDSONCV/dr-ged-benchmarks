// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Proposed by author

// 24 Streams
// 14 Equipments 
getd('../functions/');
getd('.');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect avt1_mt1 avt1_mt2 resi Qglr betaglr xchiglr ge_glr op_glr ;
stacksize(26840000);
tic;

xr =[50;150;140;140;480;480;480;220;260;260;85.8;15.8;70;15.8;174.2;3.484;3.484;73.484;170.716;170.716;55;115.716;200;20];

szx = size(xr,1);
runsize = 1;
sd = [0.284217
2.132184
2.373325
2.106504
57.932395
57.932395
69.470558
4.564205
8.205589
8.205589
0.684127
0.194481
0.596018
0.194481
4.155075
0.008974
0.008974
8.132763
3.452461
3.996801
0.283173
2.152265
3.907185
0.025517
].^(0.5);

sds = sd;
jac=jacP14();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);


[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.1, 0.2);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

[x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ]=calc_results(xfinal, jac, sigma, resGrossErrorNodalRandFi);

[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.1, 0.1, gamaMeasuremts, runsize, rj, jac_col, jac_row);

[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.05, 0.92, zadj, runsize, jac_col);

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.10, 0.77, jac_row, runsize, zr_nt_nodal_rand);

[avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.14, 0.22, runsize);

//[ avt1_mt1 avt1_mt2 avt1_nt1 avt1_nt2   avti_glr avti_glr_nt  avti_gt_mt avti_gt_nt]
runtime=toc();
//saving results
//aa = clock();
//nowtime = '_' + string(aa(4)) + '-'+ string(aa(5));
//save ('P_resumed_' + date() + nowtime +'.sav', runtime,  adjustability, detect, op_nt1, op_nt2, norm_nt1, norm_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, norm_mt, norm_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, xchi, avti_gt_mt, op_glr, op_glr_nt,  xchiglr,  xchiglr_nt, avti_glr, avti_glr_nt);
