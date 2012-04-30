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
getd('../jacobians/');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt; 

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

lb_delta = 2;
ub_delta = 7;
mtsize = ub_delta - lb_delta + 1;
op_gt_mt = zeros(1, mtsize); 
op_mt1 = zeros(jac_col, mtsize); op_mt2 = op_mt1;
op_glr_mt = op_mt1; aee_mt = op_mt1;

// we do the Overall Power Curve in two steps, first for the measurement bias, then for the leakings:
// in order to avoid undesired overwriting , the variables with "_tmp" termination
// is created
for i = lb_delta:ub_delta

    [xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, i, i, 0, 0);

    resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

    [x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ]=calc_results(xfinal, jac, sigma, resGrossErrorNodalRandFi);

    [avti_gt_mt, op_gt_mt(1:jac_col,i), op_gt_nt_tmp] = global_test(0.095, 0.095, gamaMeasuremts, runsize, rj, jac_col, jac_row);

    [avt1_mt1, avt1_mt2, op_mt1(1:jac_col,i), op_mt2(1:jac_col,i)] = measurement_test(0.0166, 0.11, zadj, runsize, jac_col);

    [avti_glr, op_glr_mt(1:jac_col,i), aee_mt(1:jac_col,i), aee_nt_tmp, op_glr_nt_tmp, avti_glr_nt_tmp ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.11, 0.21, runsize);
    clear xfinal  resRand resGrossErrorNodalRand resGrossErrorNodalRandFi;
end

// the leak bounds are 0.01 spaced, in case of other space, it is necessary to 
// set up ntsize and the loop increment
lb_leak = 2;
ub_leak = 7;
ntsize = (ub_leak - lb_leak) + 1
op_gt_nt = zeros(1, ntsize);
aee_nt  = zeros(jac_row, ntsize); op_glr_nt = aee_nt;
op_nt1 = aee_nt; op_nt2 = aee_nt;

// notice that the lb_leak and ub_leak must be divided by 100 in gererate_data function

for j = lb_leak:ub_leak

    [xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 0, 0, j/100, j/100);

    resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

    [x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ]=calc_results(xfinal, jac, sigma, resGrossErrorNodalRandFi);

    [avti_gt_mt_tmp, op_gt_mt_tmp, op_gt_nt(1:jac_row,j)] = global_test(0.095, 0.095, gamaMeasuremts, runsize, rj, jac_col, jac_row);

    [avt1_nt1, avt1_nt2, op_nt1(1:jac_row,j), op_nt2(1:jac_row,j)] = nodal_test(0.0243, 0.0935, jac_row, runsize, zr_nt_nodal);

    [avti_glr_tmp, op_glr_mt_tmp, aee_mt_tmp, aee_nt(1:jac_row,j), op_glr_nt(1:jac_row,j), avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.11, 0.21, runsize);
    clear xfinal  resRand resGrossErrorNodalRand resGrossErrorNodalRandFi;
end
   
clear op_gt_nt_tmp aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp avti_glr_tmp op_glr_mt_tmp aee_mt_tmp;
runtime=toc();
//streamNames =generateStreamName(szx);
//prettyprinttable([tokens(streamNames), string([xr, rrn(4,sd), rrn(3,adj), rrn(3,detect), rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)])],"latex")
//eqpNames = generateEqpName('', jac_row);
//prettyprinttable([tokens(eqpNames), string([rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)])],"latex")
//prettyprinttable([tokens(eqpNames), string([rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)])],"latex")
//[ op_gt_mt avti_gt_mt avt1_mt1 avt1_mt2 avti_glr avt1_nt1 avt1_nt2  avti_glr_nt runtime ]
//prettyprinttable(string([rrn(3,avt1_mt1),  rrn(3,avt1_mt2),  rrn(3,avti_glr),  rrn(3,avt1_nt1),  rrn(3,avt1_nt2),  rrn(7,avti_glr_nt)]))
//[rrn(3,op_mt1), rrn(3,op_mt2), rrn(3,op_glr_mt), rrn(7,aee_mt)]
//[rrn(3,op_nt1), rrn(3,op_nt2), rrn(3,op_glr_nt), rrn(7,aee_nt)]
//
//saving results
aa = clock();
nowtime = '_' + string(aa(4)) + '-'+ string(aa(5)) + '-'+ string(round(aa(6)));
save ('P_resumed_OP_curve' + date() + nowtime +'.sav', runtime,  adj, detect, op_nt1, op_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, avti_gt_mt, op_glr_mt, op_glr_nt, avti_glr, avti_glr_nt, aee_nt, aee_mt);
//
