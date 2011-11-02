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
// 50 Streams
// 26 Equipments 
getd('../functions/');
getd('.');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect avt1_mt1 avt1_mt2 resi Qglr betaglr xchiglr ge_glr op_glr ;
stacksize(268400000);
tic;

xr =[225.45;167.89;1332.00;1332.00;2276.90;137.50;917.89;532.38;385.51;385.51;385.51;100.32;285.19;147.69;680.07;683.07;683.07;1593.80;1593.80;582.60;582.60;1178.20;2872.50;3.84;2876.30;2870;2876.30;2876.30;3098.20;1400.00;1337.70;11.36;349.07;501.44;392.10;33.80;535.24;244.46;147.64;31.13;504.11;7.00;497.11;233.72;247.55;15.84;13.37;2.47;308.47;310.94];

szx = size(xr,1);
runsize = 2500;
sd =[2.06;1.31;11.12;11.16;18.61;1.22;7.81;4.74;3.32;3.32;3.32;0.85;2.26;1.34;6.35;6.35;6.24;14.97;14.97;5.28;5.28;9.27;23.76;0.03;25.96;25.96;25.96;25.96;25.87;11.46;12.2;0.1;2.98;4.37;3.34;0.29;1;2.15;1.2;1;4.43;0.06;4.18;1.93;2.12;0.12;0.12;0.02;2.86;2.74;];
sds = sd;
jac=jacP16();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);


[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.1, 0.2);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

[x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ]=calc_results(xfinal, jac, sigma, resGrossErrorNodalRandFi);

[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.1, 0.1, gamaMeasuremts, runsize, rj, jac_col, jac_row)

[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.05, 0.28, zadj, runsize, jac_col)

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.1, 0.1, jac_row, runsize, zr_nt_nodal_rand)

[avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.26, 0.28, runsize)

//[ avt1_mt1 avt1_mt2 avt1_nt1 avt1_nt2   avti_glr avti_glr_nt  avti_gt_mt avti_gt_nt]
runtime=toc();
//saving results
//aa = clock();
//nowtime = '_' + string(aa(4)) + '-'+ string(aa(5));
//save ('P_resumed_' + date() + nowtime +'.sav', runtime,  adjustability, detect, op_nt1, op_nt2, norm_nt1, norm_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, norm_mt, norm_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, xchi, avti_gt_mt, op_glr, op_glr_nt,  xchiglr,  xchiglr_nt, avti_glr, avti_glr_nt);
