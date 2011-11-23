// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Proposed by author
//10 Streams 
//6 Equipments
getd('../functions/');
getd('.');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect avt1_mt1 avt1_mt2 resi Qglr betaglr xchiglr ge_glr op_glr ;
stacksize(258435454);
tic;
// the real values
xr =[50;75;75;48;30;25;5;5;3;2];
//the variance

szx = size(xr,1);
runsize = 20000;
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
// i = runsize loop starts at "runsize value" and ands at runsize/i
for i =1:5
//j  = loop to evaluate mean values   
    for j=1:4

        [xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.1, 0.2);
        
        resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];
        
        [x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ]=calc_results(xfinal, jac, sigma, resGrossErrorNodalRandFi);
        
        [avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.095, 0.095, gamaMeasuremts, runsize, rj, jac_col, jac_row);
        
        [avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.05, 0.65, zadj, runsize, jac_col);
        
        [avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.1, 0.47, jac_row, runsize, zr_nt_nodal_rand);
        
        [avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.12, 0.17, runsize);
        
        //[ avt1_mt1 avt1_mt2 avt1_nt1 avt1_nt2   avti_glr avti_glr_nt  avti_gt_mt avti_gt_nt]
        runtime=toc();
        //saving results
        aa = clock();
        nowtime = '_' + string(aa(4)) + '-' + string(aa(5))+ '-'+ string(aa(6));
        save ('P_resumed_' + date() + nowtime +'.sav', runtime,  adj, detect, op_nt1, op_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, avti_gt_mt, op_glr_mt, op_glr_nt, avti_glr, avti_glr_nt, aee_mt, aee_nt, runsize);
        
        clear xfinal  resRand  resGrossErrorNodalRand resGrossErrorNodalRandFi x_sol  res  gamaMeasuremts gamaNodal zr_nt_nodal  zr_nt_nodal_rand  zadj avti_gt_mt  op_gt_mt  op_gt_nt avti_gt_mt  op_gt_mt  op_gt_nt avt1_mt1  avt1_mt2  op_mt1  op_mt2 avt1_nt1  avt1_nt2  op_nt1  op_nt2 avti_glr  op_glr_mt  aee_mt  aee_nt  op_glr_nt  avti_glr_nt;
    end   
    runsize = runsize/2; 
end