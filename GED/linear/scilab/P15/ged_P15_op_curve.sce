// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Steam metering system
//Serth, R W, and W A Heenan. 1986.
// Gross error detection and data reconciliation in steam-metering systems. 
//AIChE Journal 32: 733-747.
//Bibtex Citation

//@article{Serth1986,
//author = {Serth, R W and Heenan, W A},
//journal = {AIChE Journal},
//pages = {733--747},
//title = {{Gross error detection and data reconciliation in steam-metering systems}},
//volume = {32},
//year = {1986}
//}

// 28 Streams
// 11 Equipments 
// the measures
getd('../functions/');
getd('../jacobians/');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt; 

stacksize('max');
tic;
xr =[0.86;1;111.82;109.96;53.27;112.27;2.32;164.05;0.83;52.41;14.86;67.27;111.27;91.86;60;23.64;32.73;16.23;7.85;10.5;87.32;5.45;2.59;46.63;85.46;81.32;70.79;72.23];
szx = size(xr,1);
runsize = 2500;
sd = [0.022
0.025
2.796
2.749
1.332
2.807
0.058
4.101
0.021
1.310
0.372
1.682
2.782
2.297
1.500
0.591
0.818
0.406
0.196
0.263
2.183
0.136
0.065
1.166
2.137
2.033
1.770
1.806];
sds = sd;
var=sd.^2;
jac=jacP15();
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

xrand=generate_data_random_err(xr, sd, jac, runsize);

// we do the Overall Power Curve in two steps, first for the measurement bias, then for the leakings:
// in order to avoid undesired overwriting , the variables with "_tmp" termination
// is created

// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS analytical = -1 WLS numerical = 0  ; Absolute sum of squares = 1 ; Cauchy = 2 ;Contamined Normal = 3 ; Fair  = 4
//Hampel = 5 Logistic = 6 ; Lorenztian = 7 ; Quasi Weighted = 8
obj_function_type = -1;

for i = lb_delta:ub_delta

    [xfinal, resRand, resGrossErrorNodalRand]=generate_data_errors(xr, xrand, sd, jac, runsize, i, i, 0, 0);

    resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

    [x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);
    
    [res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ] = calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi);

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

    [xfinal, resRand, resGrossErrorNodalRand]=generate_data_errors(xr, xrand, sd, jac, runsize, 0, 0, j/100, j/100);

    resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

    [x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);
    
    [res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ] = calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi);

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
// plot the op curves
sdx=[1:ub_delta];
scf();
for i = 1:szx
    plot( sdx, op_mt1(i,:));
    p=get("hdl")
    p.children.foreground = i;
end
streamNames =generateStreamName(szx,1);
legend(streamNames,2);


