// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Martins, Márcio A.F., Carolina A. Amaro, Leonardo S. Souza, Ricardo A. Kalid, and Asher Kiperstok. 2010. 
//New objective function for data reconciliation in water balance from industrial processes. 
//Journal of Cleaner Production (March): 1-6. doi:10.1016/j.jclepro.2010.03.014. http://linkinghub.elsevier.com/retrieve/pii/S0959652610001149.

//Bibtex Citation
//@article{Martins2010,
//author = {Martins, M\'{a}rcio A.F. and Amaro, Carolina A. and Souza, Leonardo S. and Kalid, Ricardo A. and Kiperstok, Asher},
//doi = {10.1016/j.jclepro.2010.03.014},
//file = {::},
//issn = {09596526},
//journal = {Journal of Cleaner Production},
//month = mar,
//pages = {1--6},
//title = {{New objective function for data reconciliation in water balance from industrial processes}},
//url = {http://linkinghub.elsevier.com/retrieve/pii/S0959652610001149},
//year = {2010}
//}

// 13 Streams
// 8 Equipments 
getd('../functions/');
getd('.');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect avt1_mt1 avt1_mt2 resi Qglr betaglr xchiglr ge_glr op_glr ;
stacksize(268400000);
tic;

xr =[28;5;5;7;11;4;3;5;5;4;3;14;14];
szx = size(xr,1);
runsize = 2500;
sd = [0.274988
0.049980
0.172480
0.145203
0.372073
0.127073
0.136220
0.045407
0.095223
0.072520
0.063537
0.147000
0.130634];
sds = sd;
jac=jacP10();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);

[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.1, 0.2);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

[x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ]=calc_results(xfinal, jac, sigma, resGrossErrorNodalRandFi);

[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.1, 0.1, gamaMeasuremts, runsize, rj, jac_col, jac_row);

[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.051, 0.75, zadj, runsize, jac_col);

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.1, 0.57, jac_row, runsize, zr_nt_nodal_rand);

[avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.12, 0.18, runsize);

//[ avt1_mt1 avt1_mt2 avt1_nt1 avt1_nt2   avti_glr avti_glr_nt  avti_gt_mt avti_gt_nt]
runtime=toc();
//saving results
//aa = clock();
//nowtime = '_' + string(aa(4)) + '-'+ string(aa(5));
//save ('P_resumed_' + date() + nowtime +'.sav', runtime,  adjustability, detect, op_nt1, op_nt2, norm_nt1, norm_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, norm_mt, norm_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, xchi, avti_gt_mt, op_glr, op_glr_nt,  xchiglr,  xchiglr_nt, avti_glr, avti_glr_nt);
