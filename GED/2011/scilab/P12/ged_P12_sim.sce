// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Mitsas, Christos L. 2010. Data reconciliation and variable classification by null space methods. 
//Measurement 43, no. 5 (June): 702-707.
// http://apps.isiknowledge.com/full_record.do?product=UA&search_mode=GeneralSearch&qid=2&SID=2A@bF9dN34I72L1Am9M&page=2&doc=17&colname=WOS.

//Bibtex Citation

//@article{Mitsas2010a,
//author = {Mitsas, Christos L.},
//journal = {Measurement},
//month = jun,
//number = {5},
//pages = {702--707},
//publisher = {ELSEVIER SCI LTD},
//title = {{Data reconciliation and variable classification by null space methods}},
//url = {http://apps.isiknowledge.com/full\_record.do?product=UA\&search\_mode=GeneralSearch\&qid=2\&SID=2A@bF9dN34I72L1Am9M\&page=2\&doc=17\&colname=WOS},
//volume = {43},
//year = {2010}
//}

// 12 Streams
// 6 Equipments 
getd('../functions/');
getd('.');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect avt1_mt1 avt1_mt2 resi Qglr betaglr xchiglr ge_glr op_glr ;
stacksize(26840000);
tic;
xr =[100;30;40;30;15;5;10;70;10;20;80;100];
szx = size(xr,1);
runsize = 2500;
//the variance or uncertainties  proposed by the original author
sd = [1
0.3
0.4
0.3
0.3
0.1
0.2
0.7
0.3
0.4
0.8
1];
sds=sd;
jac=jacP12();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);


[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.1, 0.2);

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

[x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ]=calc_results(xfinal, jac, sigma, resGrossErrorNodalRandFi);

[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.1, 0.1, gamaMeasuremts, runsize, rj, jac_col, jac_row);

[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.05, 0.71, zadj, runsize, jac_col);

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.105, 0.48, jac_row, runsize, zr_nt_nodal_rand);

[avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.14, 0.22, runsize);

//[ avt1_mt1 avt1_mt2 avt1_nt1 avt1_nt2   avti_glr avti_glr_nt  avti_gt_mt avti_gt_nt]
runtime=toc();
//saving results
//aa = clock();
//nowtime = '_' + string(aa(4)) + '-'+ string(aa(5));
//save ('P_resumed_' + date() + nowtime +'.sav', runtime,  adjustability, detect, op_nt1, op_nt2, norm_nt1, norm_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, norm_mt, norm_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, xchi, avti_gt_mt, op_glr, op_glr_nt,  xchiglr,  xchiglr_nt, avti_glr, avti_glr_nt);
