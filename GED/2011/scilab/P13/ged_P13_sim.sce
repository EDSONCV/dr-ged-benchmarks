// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Fictitious but realistic mineral processing plant
//Alhaj-Dibo, Moustapha, Didier Maquin, and José Ragot. 2008.
//Data reconciliation: A robust approach using a contaminated distribution.
//Control Engineering Practice 16, no. 2 (February): 159-170.
// http://www.sciencedirect.com/science/article/B6V2H-4N4406D-1/2/50cac92b050f160a20a795faec990dc7.

//Bibtex Citation

//@article{Alhaj-Dibo2008,
//author = {Alhaj-Dibo, Moustapha and Maquin, Didier and Ragot, Jos\'{e}},
//isbn = {0967-0661},
//journal = {Control Engineering Practice},
//keywords = {Data reconciliation,Gross error detection,Linear and bilinear mass balances,Robust estimation},
//month = feb,
//number = {2},
//pages = {159--170},
//title = {{Data reconciliation: A robust approach using a contaminated distribution}},
//url = {http://www.sciencedirect.com/science/article/B6V2H-4N4406D-1/2/50cac92b050f160a20a795faec990dc7},
//volume = {16},
//year = {2008}
//}

// 16 Streams
// 9 Equipments 
getd('../functions/');
getd('.');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect avt1_mt1 avt1_mt2 resi Qglr betaglr xchiglr ge_glr op_glr ;
stacksize(268400000);
tic;
xr =[25;27;22;2;20;24;14;10;10;4;5;7;1;8;5;3];

szx = size(xr,1);
runsize = 2500;
// in original paper the standard deviation is given. so it must be squared.
sd=[1
1.325
1.46
0.20
0.916
1.101
1.04
0.472
0.401
0.207
0.3
0.328
0.052
0.369
0.25
0.385
];

sds = sd;
jac=jacP13();
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
