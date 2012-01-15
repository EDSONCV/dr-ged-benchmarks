// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Proposed by author
//10 Streams 
//6 Equipments

getd('../functions/');
clear rerror rerror1 xr xrs sd sds x_sol grerror grerrors mySign xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj d_adj zadj_alt Wbar_alt x_chi ge_gt nge_gt ge_nt1_i ge_nt1_j ge_nt2_i ge_nt2_j norm_nt1 norm_nt2 norm_mt norm_mt2 ge_mt_i ge_mt_j ge_mt_alt_i ge_mt1_alt_i ge_mt1_alt_j ge_mt2_alt_i ge_mt2_alt_j ge_mt_alt2_j ge_mt1_i ge_mt1_j ge_mt2_i ge_mt2_j Vinv_r zr_nt_alt ge_nt_alt_i ge_nt_alt_j op_nt1 op_nt2 op_mt1 op_mt2 op_mt1_alt op_mt2_alt op_nt_alt_1 op_nt_alt_2 adjustability detect ge_mt1_i_r ge_mt1_j_r ge_mt21_i_r ge_mt2_j_r  ge_mt1_i_ge ge_mt1_j_ge ge_mt1_alt_i_ge ge_mt1_alt_j_ge ge_mt2_j_ge ge_mt2_alt_i_ge ge_mt2_alt_j_ge e_mt1_ij_ge e_mt2_ij_ge ge_mt1_alt_ij_ge ge_mt2_alt_ij_ge avt1_mt1 avt1_mt2 avt1_mt1_alt avt1_mt2_alt;
tic;
// the real values
xr =[50;75;75;48;30;25;5;5;3;2];
//the variance

szx = size(xr,1);
runsize = 2500;
sd=[1
1
1
1
1
1
0.15
0.15
0.1
0.1
].^0.5;
sds = sd;
jac=jacP7();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);


[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.02, 0.07);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

[x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ]=calc_results(xfinal, jac, sigma, resGrossErrorNodalRandFi);

[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.095, 0.095, gamaMeasuremts, runsize, rj, jac_col, jac_row);

[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.0126, 0.12, zadj, runsize, jac_col);

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.0193, 0.11, jac_row, runsize, zr_nt_nodal);

[avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.12, 0.207, runsize);


runtime=toc();
//streamNames =generateStreamName(szx);
//prettyprinttable([tokens(streamNames), string([xr, rrn(4,sd), rrn(3,adj), rrn(3,detect), rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)])],"latex")
//eqpNames = "A B C D E F";
//prettyprinttable([tokens(eqpNames), string([rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)])],"latex")
//[ op_gt_mt avti_gt_mt avt1_mt1 avt1_mt2 avti_glr avt1_nt1 avt1_nt2  avti_glr_nt runtime ]
//prettyprinttable(string([rrn(3,avt1_mt1),  rrn(3,avt1_mt2),  rrn(3,avti_glr),  rrn(3,avt1_nt1),  rrn(3,avt1_nt2),  rrn(7,avti_glr_nt)]))
//[rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)]
//[rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)]
//
//
////saving results
//aa = clock();
//nowtime = '_' + string(aa(4)) + '-'+ string(aa(5));
//save ('P_resumed_' + date() + nowtime +'.sav', runtime,  adj, detect, op_nt1, op_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, avti_gt_mt, op_glr_mt, op_glr_nt, avti_glr, avti_glr_nt, aee_nt, aee_mt);
//
