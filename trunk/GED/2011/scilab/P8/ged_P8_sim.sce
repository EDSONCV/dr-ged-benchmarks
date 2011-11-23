// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Rao, R Ramesh, and Shankar Narasimhan. 1996.
//“Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” 
//Industrial & Engineering Chemistry Research 35:1362-1368. 
//http://dx.doi.org/10.1021/ie940538b.
//Bibtex Citation

//@article{Rao1996,
//author = {Rao, R Ramesh and Narasimhan, Shankar},
//isbn = {0888-5885},
//journal = {Industrial \& Engineering Chemistry Research},
//month = apr,
//number = {4},
//pages = {1362--1368},
//publisher = {American Chemical Society},
//title = {{Comparison of Techniques for Data Reconciliation of Multicomponent Processes}},
//url = {http://dx.doi.org/10.1021/ie940538b},
//volume = {35},
//year = {1996}
//}

// 12 Streams
// 7 Equipments 

getd('../functions/');
getd('.');
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect avt1_mt1 avt1_mt2 resi Qglr betaglr xchiglr ge_glr op_glr ;
stacksize(258435454);
tic;

xr =[3600;1850;1750;2837;730;25;8;137;2700;58;672];
szx = size(xr,1);
runsize = 20000;
sd = [111.2100
57.0000
54.2100
87.3000
22.1100
0.7800
0.2280
3.1500
84.9600
1.7100
20.0400];
sds =sd;
jac=jacP8();
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


    	[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.1, 0.095, gamaMeasuremts, runsize, rj, jac_col, jac_row);

	    [avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(0.05, 0.69, zadj, runsize, jac_col);

	    [avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.1, 0.34, jac_row, runsize, zr_nt_nodal_rand);

	    [avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt ]=calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, 0.2, 0.29, runsize);

        runtime=toc();
        //saving results
        aa = clock();
        nowtime = '_' + string(aa(4)) + '-' + string(aa(5))+ '-'+ string(aa(6));
        save ('P_resumed_' + date() + nowtime +'.sav', runtime,  adj, detect, op_nt1, op_nt2, avt1_nt1, avt1_nt2, op_mt1, op_mt2, avt1_mt1, avt1_mt2, op_gt_mt, op_gt_nt, avti_gt_mt, op_glr_mt, op_glr_nt, avti_glr, avti_glr_nt, aee_mt, aee_nt, runsize);
        
        clear xfinal  resRand  resGrossErrorNodalRand resGrossErrorNodalRandFi x_sol  res  gamaMeasuremts gamaNodal zr_nt_nodal  zr_nt_nodal_rand  zadj avti_gt_mt  op_gt_mt  op_gt_nt avti_gt_mt  op_gt_mt  op_gt_nt avt1_mt1  avt1_mt2  op_mt1  op_mt2 avt1_nt1  avt1_nt2  op_nt1  op_nt2 avti_glr  op_glr_mt  aee_mt  aee_nt  op_glr_nt  avti_glr_nt;
    end   
    runsize = runsize/2; 
end
