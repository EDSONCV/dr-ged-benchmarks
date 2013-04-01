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
getd('../jacobians/');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt; 

stacksize('max');
tic;

xr =[225.45;167.89;1332.00;1332.00;2276.90;137.50;917.89;532.38;385.51;385.51;385.51;100.32;285.19;147.69;680.07;683.07;683.07;1593.80;1593.80;582.60;582.60;1178.20;2872.50;3.84;2876.30;2870;2876.30;2876.30;3098.20;1400.00;1337.70;11.36;349.07;501.44;392.10;33.80;535.24;244.46;147.64;31.13;504.11;7.00;497.11;233.72;247.55;15.84;13.37;2.47;308.47;310.94];

szx = size(xr,1);
runsize = 2500;
sd =[2.06;1.31;11.12;11.16;18.61;1.22;7.81;4.74;3.32;3.32;3.32;0.85;2.26;1.34;6.35;6.35;6.24;14.97;14.97;5.28;5.28;9.27;23.76;0.03;25.96;25.96;25.96;25.96;25.87;11.46;12.2;0.1;2.98;4.37;3.34;0.29;1;2.15;1.2;1;4.43;0.06;4.18;1.93;2.12;0.12;0.12;0.02;2.86;2.74;];
sds = sd;
var=sd.^2;
jac=jacP16();
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

