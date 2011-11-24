// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Rosenberg, J and Mah, R S H and Iordache, C
//Evaluation of Schemes for Detecting and Identifying Gross Errors in Process Data
//Ind. & Eng. Chem. Proc. Des. Dev, Vol., V. 26:555--564

//Bibtex Citation
//@article{Rosenberg1987a, 
//author = {Rosenberg, J and Mah, R S H and Iordache, C},
//journal = {Ind. \& Eng. Chem. Proc. Des. Dev, Vol.},
//pages = {555--564},
//title = {{Evaluation of Schemes for Detecting and Identifying Gross Errors in Process Data}},
//volume = {26},
//year = {1987}
//}


//7 Streams 
//4 Equipments
getd('../functions/');
getd('.');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect avt1_mt1 avt1_mt2 resi Qglr betaglr xchiglr ge_glr op_glr ;
stacksize(26840000);
tic;

// the real values
xr =[5 15 15 5 10 5 5 ]';

szx = size(xr,1);
runsize = 2500;
//the variance
sd=[1
1
1
1
1
1
1].^0.5;
//sd = [0.13
//0.38
//0.38
//0.13
//0.25
//0.13
//0.13];
sds = sd;
sds =sd;
jac=jacP6();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);


[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.1, 0.2);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

[x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ]=calc_results(xfinal, jac, sigma, resGrossErrorNodalRandFi);


[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.095, 0.095, gamaMeasuremts, runsize, rj, jac_col, jac_row);

[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.05, 0.52, zadj, runsize, jac_col);

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.1, 0.34, jac_row, runsize, zr_nt_nodal_rand);

[avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.11, 0.21, runsize);

//[ avt1_mt1 avt1_mt2 avt1_nt1 avt1_nt2   avti_glr avti_glr_nt  avti_gt_mt avti_gt_nt]
runtime=toc();
//saving results
//aa = clock();
//nowtime = '_' + string(aa(4)) + '-'+ string(aa(5));
//save ('P_resumed_' + date() + nowtime +'.sav', runtime,  adjustability, detect, op_nt1, op_nt2, norm_nt1, norm_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, norm_mt, norm_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, xchi, avti_gt_mt, op_glr, op_glr_nt,  xchiglr,  xchiglr_nt, avti_glr, avti_glr_nt);
