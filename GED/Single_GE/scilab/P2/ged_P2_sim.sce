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
getd('../functions/');
getd('.');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect avt1_mt1 avt1_mt2 resi Qglr betaglr xchiglr ge_glr op_glr ;
stacksize(26840000);
tic;
xr =[11;10;21;11;5.5;5.5];

szx = size(xr,1);
runsize = 2500;
//the variance proposed by this work
sd =[0.032
0.026
0.120
0.033
0.052
0.015].^(0.5);
sds = sd;
jac=jacP2();
rj=rank(jac);
jac_col = size(jac,2);
jac_row = size(jac,1);
sigma=diag(sds.^2);


[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.02, 0.07);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

[x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ]=calc_results(xfinal, jac, sigma, resGrossErrorNodalRandFi);

[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.1, 0.1, gamaMeasuremts, runsize, rj, jac_col, jac_row);

[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.028, 0.16, zadj, runsize, jac_col);

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.037, 0.107, jac_row, runsize, zr_nt_nodal);

[avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.16, 0.238, runsize);


runtime=toc();
//streamNames =generateStreamName(szx);
//prettyprinttable([tokens(streamNames), string([xr, rrn(4,sd), rrn(3,adj), rrn(3,detect), rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)])],"latex")
//eqpNames = generateEqpName('', jac_row);
//prettyprinttable([tokens(eqpNames), string([rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)])],"latex")
//[ op_gt_mt avti_gt_mt avt1_mt1 avt1_mt2 avti_glr avt1_nt1 avt1_nt2  avti_glr_nt runtime ]
//prettyprinttable(string([rrn(3,avt1_mt1),  rrn(3,avt1_mt2),  rrn(3,avti_glr),  rrn(3,avt1_nt1),  rrn(3,avt1_nt2),  rrn(7,avti_glr_nt)]))
//[rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)]
//[rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)]
//
////saving results
//aa = clock();
//nowtime = '_' + string(aa(4)) + '-'+ string(aa(5));
//save ('P_resumed_' + date() + nowtime +'.sav', runtime,  adj, detect, op_nt1, op_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, avti_gt_mt, op_glr_mt, op_glr_nt, avti_glr, avti_glr_nt, aee_nt, aee_mt);
//
