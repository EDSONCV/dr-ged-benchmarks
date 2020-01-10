// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Mandel, Denis, Ali Abdollahzadeh, Didier Maquin, and Jos� Ragot. 1998. 
//Data reconciliation by inequality balance equilibration: a LMI approach. 
//International Journal of Mineral Processing 53, no. 3 (April): 157-169. 
//http://www.sciencedirect.com/science/article/B6VBN-3VM1X8N-3/2/8bffe94a1153eea8647eed5af0031d36.

//Bibtex Citation
//@article{Mandel1998,
//author = {Mandel, Denis and Abdollahzadeh, Ali and Maquin, Didier and Ragot, Jos�},
//isbn = {0301-7516},
//journal = {International Journal of Mineral Processing},
//keywords = {Linear Matrix Inequality Techniques,data reconciliation,error detection,error isolation},
//month = apr,
//number = {3},
//pages = {157--169},
//title = {{Data reconciliation by inequality balance equilibration: a LMI approach}},
//url = {http://www.sciencedirect.com/science/article/B6VBN-3VM1X8N-3/2/8bffe94a1153eea8647eed5af0031d36},
//volume = {53},
//year = {

// 12 Streams
// 5 Equipments 
getd('../../');
getd('../../../jacobians/');
getd('../method/');
getd('../method/pls');
cd  '../../'
clear xr sd sds x_sol xfinal jac jac_col jac_col rj sigma sigam_inv res  V V_inv diag_diag_V Wbar gama zr_nt adj zadj   Wbar_alt  adjustability detect resi Qglr betaglr xchiglr ge_glr op_glr;
clear avti_gt_mt op_gt_mt op_gt_nt_tmp avt1_mt1 avt1_mt2 op_mt1 op_mt2 avti_glr op_glr_mt aee_mt aee_nt_tmp op_glr_nt_tmp avti_glr_nt_tmp avti_gt_mt_tmp op_gt_mt_tmp op_gt_nt avt1_nt1 avt1_nt2 op_nt1 op_nt2 avti_glr_tmp op_glr_mt_tmp aee_mt_tmp aee_nt op_glr_nt avti_glr_nt; 

stacksize('max');
tic;

xr =[230;21;209;35;174;15;159;50;209;94;115;44];
//the variance proposed by the original author
//sd = [37.575
//1.08
//5
//1.825
//2
//0.88
//7.245
//1
//5
//2
//18.1
//2.385
//];
szx = size(xr,1);
runsize = 500;
// we are testing equal sigma here
sd=ones(12,1);
sds = sd;
var=sd.^2;
jac=jacP9();
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);
//sigma=eye(12,12);
[adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac);
//[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 2, 7, 0.1, 0.2);
[xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, 5, 9, 0.07, 0.15);


resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

//observability/redundancy tests
//user can set unmeasured streams here, if this vector is empty, all streams are measured                  
umeas_P9 = [];
[red_P9, just_measured_P9, observ_P9, non_obs_P9, spec_cand_P9] = qrlinclass(jac,umeas_P9);
measured_P9 = setdiff([1:length(xr)], umeas_P9);
red = measured_P9;//
        
// to run robust reconciliation,, one must choose between the folowing objective functions to set up the functions path and function parameters:
//WLS analytical = -1 WLS numerical = 0  ; Absolute sum of squares = 1 ; Cauchy = 2 ;Contamined Normal = 3 ; Fair  = 4
//Hampel = 5 Logistic = 6 ; Lorenztian = 7 ; Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = 2;

[x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);

[res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ] = calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi);

// for equal sigma
[avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(0.08, 0.08, gamaMeasuremts, runsize, rj, jac_col, jac_row);
// cauchy
[avt1_mt1, avt1_mt2, op_mt1, op_mt2] =  measurement_test(0.00015, 0.0018, zadj, runsize, jac_col);

[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(0.0175, 0.101, jac_row, runsize, zr_nt_nodal);


nvalidate = 10; lower_bias = 5; delta_bias = 1; upper_bias = 9; lower_leak = 0.07; delta_leak = 0.02; upper_leak = 0.15; 
//cauchy 
alfa_gt_mt = 0.08; alfa_gt_nt = 0.08; alfa_mt1 = 0.00015; alfa_mt2 =0.0018; alfa_nt1 = 0.0175; alfa_nt2 = 0.101;
pause
is_multi = 0;
clear res gamaMeasuremts gamaNodal zr_nt_nodal zr_nt_nodal_rand zadj x_sol resGrossErrorNodalRandFi;
[p9_train, p9_validate]  = generate_trainning2(xr, sd, jac, runsize, nvalidate, lower_bias, delta_bias, upper_bias, lower_leak,delta_leak,upper_leak, alfa_gt_mt,alfa_gt_nt,alfa_mt1,alfa_mt1, alfa_nt1, alfa_nt2,obj_function_type,is_multi);
ndatainterval = 5
//pause
[list_models_P9, p9_stat] = generate_pls_models_m( 'P9', 12, 6, p9_train, p9_validate, nvalidate,ndatainterval);
[avti_meas, op_meas, selectivity_meas, aee_meas, avti_eqp, op_eqp, selectivity_eqp, aee_eqp] = get_lit_info(p9_stat, jac_col, jac_row)  
list_models_P9

runtime=toc();
cd 'pmgei_method/problems';
