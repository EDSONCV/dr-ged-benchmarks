function [final_train, final_validate,varargout]  = generate_trainning2(xr, sd, jac, runsize, n_validate, meas_bias_low, meas_bias_inc, meas_bias_up, leak_low, leak_inc, leak_up, gt_mt_tol, gt_nt_tol, mt_tol_1, mt_tol_2, nt_tol_1, nt_tol_2,obj_function_type, is_multiple)
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// This function generate the trainning dataset based on the problem passed.
//input arguments:
// xr -the real values of the measurements
// sd -the standard deviation of the measurements
// jac -the jacobean of the problem
// runsize - the number of the runs for each measurement bias or leak
// n_validate -  number of samples of each case to be used for validation
// meas_bias_low - the lower limit of the measuremt error bias to be generated
// meas_bias_inc - the increment of the measuremt error bias to be generated
// meas_bias_up - the upper limit of the measuremt error bias to be generated
//              -example: 3, 1, 7 will generate a trainning and validation 
//              data set with errors equals to sigma*[3, 4, 5 , 6, 7]            
// leak_low - the lower limit of the leaking to be generated
// leak_inc - the increment of the leaking to be generated
// leak_up - the upper limit of the leaking to be generated
//              - the genearated leaking is proportional to totalnode flow
//                which is the sum of total flow that enters or leaves the
//                equipment
//              -example: 0.07, 0.02, 0.15 will generate a trainning and validation 
//              data set with erros equals to totalnodeflow*[0.07, 0.09, 0.11, 0.13, 0.15]   
// NOTICE THAT THIS 2 VECTORS MUST HAVE THE SAME NUMBER OF ELEMENTS, EXAMPLE:    
//            [0.07, 0.09, 0.11, 0.13, 0.15] AND
//            [3, 4, 5 , 6, 7] OK! BOTH VECTORS WITH THE SAME NUMBER OF ELEMENTS
// gt_mt_tol - the tolerance of the global test to consider a global balance 
//             with a gross error with the MEASUREMENT BIAS DATASET. It must be previously set in order to keep the 
//             AVT1 (probability of false alarm) as 0.1 within the measurement error
//             range (must be determined outside this routine)
// gt_nt_tol - the tolerance of the global test to consider a global balance 
//             with a gross error with the LEAKING DATASET.
//              It must be previously set in order to keep the 
//             AVT1 (probability of false alarm) as 0.1 within the leakinng 
//             range (must be determined outside this routine)
//mt_tol_1 -  the tolerance of the measurement test for the SIMPLE measuremente test.
//            It must be previously set in order to keep the 
//             AVT1 (probability of false alarm) as 0.1 within the measurement error
//             range (must be determined outside this routine)
//mt_tol_2 -  the tolerance of the measurement test for the MAXIMUM POWER measuremente test.
//            It must be previously set in order to keep the 
//             AVT1 (probability of false alarm) as 0.1 within the measurement error
//             range (must be determined outside this routine)
// nt_tol_1 - the tolerance of the nodal SIMPLE test.
//              It must be previously set in order to keep the 
//             AVT1 (probability of false alarm) as 0.1 within the leakinng 
//             range (must be determined outside this routine)
// nt_tol_2 - the tolerance of the MAXIMUM POWER nodal test.
//              It must be previously set in order to keep the 
//             AVT1 (probability of false alarm) as 0.1 within the leakinng 
//             range (must be determined outside this routine)
//
//
//
//output arguments
// train - the dataset used for trainning the model 
// validate - the dataset for cross validation  (data not included in train dataset)
// varargout(1) - The x used for generation of the trainning data
// varargout(2) - The x used for generation of the validation data

szx = size(xr,1);

sds = sd;
sds =sd;
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);
var =  sd.^2;
jj = [leak_low: leak_inc: leak_up];
jj_counter = 1;
train=[];
validate =[];
final_train=[];
final_validate =[];
final_xfinal=[];
final_x_sol =[];
final_xm_leaks = [];
final_x_sol_leak = [];


final_adj = [];
// generate pure random bias
[xrs] = generate_data_random_err(xr, sd, jac, runsize);

for ii = meas_bias_low:meas_bias_inc:meas_bias_up
    
//generate the random measurement bias, the corresponding residuals, the residuals with only random error and the leaks    

    [xfinal, resRand, resGrossErrorNodalRand, leaks, mysign]=generate_data_errors(xr, xrs, sd, jac, runsize, ii, ii, jj(jj_counter), jj(jj_counter))
// genearate the measuremente from the leakings    
    [xm_leaks] = generate_meas_from_leaks(xr, resRand, leaks, jac, runsize);
// concatenate the vectors    
    resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];
// original    
//    [x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_mt, zr_nt_nodal_rand, zadj, zadj_nodal] = calc_results_tst_pls(xfinal, jac, sigma, resGrossErrorNodalRandFi, leaks, xm_leaks);
//    
//    
//    [avti_gt_mt, avti_gt_nt, op_gt_mt, op_gt_nt, ge_gt_delete_mt ge_gt_delete_nt, ge_gt_delete_low] = global_test2(gt_mt_tol, gt_nt_tol, gamaMeasuremts, gamaNodal, runsize, rj, jac_col, jac_row);
//    
//    [avt1_mt1, avt1_mt2, op_mt1, op_mt2, ge_mt1_indexu, ge_mt2_indexu, ge_mt1_indexu_low, ge_mt2_indexu_low] = measurement_test3(mt_tol_1, mt_tol_2, zadj, runsize, jac_col);
//    
//    [avt1_nt1, avt1_nt2, op_nt1, op_nt2, ge_nt1_indexu, ge_nt2_indexu, ge_nt1_indexu_low, ge_nt2_indexu_low] = nodal_test2(nt_tol_1, nt_tol_2, jac_row, runsize, zr_nt_nodal);
    
// witout abs

// caculate the test statistics:  calc_results_tst_pls_sig is different from  calc_results_tst_pls because the former do not make the absolute operation in zr_nt_nodal and zadj (to get more information for pls)
    [x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_mt, zr_nt_nodal_rand, zadj, zadj_nodal, x_sol_leak] = calc_results_tst_pls_sig(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type, is_multiple, leaks, xm_leaks);
 
    
    [avti_gt_mt, avti_gt_nt, op_gt_mt, op_gt_nt, ge_gt_delete_mt ge_gt_delete_nt, ge_gt_delete_low] = global_test2(gt_mt_tol, gt_nt_tol, abs(gamaMeasuremts), abs(gamaNodal), runsize, rj, jac_col, jac_row);
// the same as measurement_test, but also returns the indices of elements below the test statistics    
    [avt1_mt1, avt1_mt2, op_mt1, op_mt2, ge_mt1_indexu, ge_mt2_indexu, ge_mt1_indexu_low, ge_mt2_indexu_low] = measurement_test3(mt_tol_1, mt_tol_2, abs(zadj), runsize, jac_col, 0);
// the same as nodal_test, but also returns the indices of elements below the test statistics    
    [avt1_nt1, avt1_nt2, op_nt1, op_nt2, ge_nt1_indexu, ge_nt2_indexu, ge_nt1_indexu_low, ge_nt2_indexu_low] = nodal_test2(nt_tol_1, nt_tol_2, jac_row, runsize, abs(zr_nt_nodal));
//disp('inside generate_trainning 2: after caculate the test statistics')
//pause   

//pause

    index_to_clear_low = [];    
    index_to_clear_mt = [];    
    index_to_clear_nt = [];
    // remove unused components (measurement errors which random error were higher 
    //than the test statistics for three tests)    
    
    if ii == meas_bias_low then
        index_to_clear_low = unique([ge_gt_delete_low, ge_mt1_indexu_low, ge_nt1_indexu_low]);
    else
        index_to_clear_low = unique([[1:runsize],ge_gt_delete_low, ge_mt1_indexu_low, ge_nt1_indexu_low]);

    end
    
    
//    index_to_clear_low = unique([ge_gt_delete_low, ge_mt1_indexu_low, ge_nt1_indexu_low]);
    // remove unused components (measurement errors which gross error with lower
    // than the test statistics) for measurement test   
    ge_mt1_indexu_filter =  setdiff([runsize + 1: (jac_col + 1)* runsize],ge_mt1_indexu);
    index_to_clear_mt = unique([ge_mt1_indexu_filter, ge_gt_delete_mt]);
    disp('before index to clear nt')
//pause
//    index_to_clear_mt = unique([ge_mt1_indexu, ge_gt_delete_mt]);
    // remove unused components (leakings which gross error with lower
    // than the test statistics) for nodal test   
    index_to_clear_nt = unique([ge_nt1_indexu, ge_gt_delete_nt]);
    
    //lambdaRes_new = lambdaRes;
    //lambdaRes_new(:, index_to_clear) = [];
    //lambdaVar_new = lambdaVar;
    //lambdaVar_new(:,index_to_clear) = [];
    
    // clear gamaMeasurements vector (gamma for measurement bias dataset)
    gamaMeasuremts_new = gamaMeasuremts;
    gamaMeasuremts_new([index_to_clear_low, index_to_clear_mt]) = [];
    // clear zadj vector    
    zadj_new = zadj;
//    zadj_new = xfinal-x_sol;

    zadj_new([index_to_clear_low, index_to_clear_mt],:) = [];
    
    
    // clear zr_nt_nodal_mt (residuals of the measurement bias dataset) vector        
    zr_nt_nodal_mt_new = zr_nt_nodal_mt;
    zr_nt_nodal_mt_new([index_to_clear_low, index_to_clear_mt], :) = [];
   
    xfinal([index_to_clear_low, index_to_clear_mt], :) = [];
    x_sol([index_to_clear_low, index_to_clear_mt], :) = [];
 //  pause
    // clear gamaNodal (gamma for leaking dataset)
    gamaNodal_new = gamaNodal;    
    gamaNodal_new([index_to_clear_nt]) = [];
    // since gammaNodal(1:runsize) is already included in gamaMeasurements, the
    // first "runsize" elements will be removed
    gamaNodal_new(1:runsize) = [];
    // clear zr_nt_nodal (residuals of the leakings dataset) vector        
    //pause
    zr_nt_nodal_new = zr_nt_nodal;
//    zr_nt_nodal_new = resGrossErrorNodalRandFi;
    zr_nt_nodal_new(index_to_clear_nt, :) = [];
     // since zr_nt_nodal(1:runsize,:) is already included in zr_nt_nodal_mt, the
    // first "runsize" elements will be removed
    zr_nt_nodal_new(1:runsize,:) = [];
    // clear zadj_nodal (adjustments statistics of the leakings dataset) vector        
//    zadj_nodal_new = zadj_nodal;
//pause
//    zadj_nodal_new = [zadj_new(1:runsize,:);xm_leaks - x_sol_leak];
    zadj_nodal_new = [zadj_new(1:runsize,:);zadj_nodal(runsize+1:$,:)];
    zadj_nodal_new(index_to_clear_nt,:) = [];
    // since zadj_nodal(1:runsize,:) is already included in zadj, the
    // first "runsize" elements will be removed
    zadj_nodal_new(1:runsize,:) = [];
    
    
    // arranging MT trainning data
    y=zeros(length(gamaMeasuremts)-runsize,jac_col + jac_row);
    
    y2=ii*ones(runsize ,1);
    
    for i = 1 : jac_col
//            y((i-1)*runsize+1:i*runsize,i) = y2;
            y((i-1)*runsize+1:i*runsize,i) = y2.*mysign(:,i);
//            pause
    end
    
    ytop=zeros(runsize,jac_col+ jac_row);
    
    
    y=[ytop;y];
    
    y([index_to_clear_low, index_to_clear_mt],:) = [];
    
    
    train=[gamaMeasuremts_new, zadj_new, zr_nt_nodal_mt_new, y];

    // arranging NT trainning data
    
    y_nt=zeros(length(gamaNodal), jac_col+ jac_row);
    
    // we multiply for 10 to present a good scalling
    
    y2=10*jj(jj_counter)*ones(runsize ,1);
    
    for i = 1 : jac_row
      
            y_nt((i-1)*runsize + 1 + runsize:(i)*runsize + runsize , i + jac_col) = y2;
    end
    
    y_nt(index_to_clear_nt,:) = [];
    
    y_nt(1:runsize, :) = [];
    disp('before xmleaks');
//pause
if length(index_to_clear_nt) > 0 then
    xm_leaks([index_to_clear_nt - runsize],:) = [];
    x_sol_leak([index_to_clear_nt - runsize],:) = [];    
end

    xfinal = [xfinal;xm_leaks];
    x_sol = [x_sol;x_sol_leak];
//    xm_leaks(1:runsize,:) = [];
//    x_sol_leak(1:runsize,:) = [];
    
    

    train=[train; gamaNodal_new, zadj_nodal_new, zr_nt_nodal_new, y_nt];
    validate_index=[];

// pick some lines for validation
    for i = 0: jac_col+jac_row
//        disp('i: ', i)
        if i == 0 then
            index = find(train(:,jac_row + jac_col + 2 + i) == 0); 
        elseif i <= jac_col then
//pause            
           index = find(train(:,jac_row + jac_col + 2 + i - 1) <> 0); 
           szi_mt = length(index);
           b_meas = jac_row + jac_col + 2 + i - 1;
           printf('b_meas %d, i = %d, ii = %d, jj= %d, index meas: %d \n', b_meas, i, ii, jj_counter, szi_mt)// , %d, %d' ,  ) 
       else
           index = find(train(:,jac_row + jac_col + 2 + i - 1) == 10*jj(jj_counter)); 
           b_eqp = jac_row + jac_col + 2 + i - 1;
           szi_nt = length(index);
           printf('b_eqp %d, i = %d, ii = %d, jj= %d, index eqp: %d \n', b_eqp, i, ii, jj_counter ,szi_nt)// , %d, %d' ,  )    
       end
       // select the first "n_validate" for validation

//       pause
       if length(index) == 0 then
       printf('Warning : No values for training and/or validation found empty column = %d, ii = %d, jj= %d \n', i, ii, jj_counter)// , %d, %d' ,  )    
       elseif length(index) <  n_validate then           
           validate_index = [validate_index, index(1:$)];    
       else
           validate_index = [validate_index, index(1:n_validate)];
       end
       index =[];
       
    end
    validate_index = unique(validate_index);
//    pause
    //add the selected data to validation set
    final_validate = [final_validate; train(validate_index,:)];
    // remove the validation set 
    train(validate_index,:) = [];

    final_train = [final_train; train];

    xfinal(validate_index,:) = [];
    final_xfinal = [final_xfinal; xfinal];

    x_sol(validate_index,:) = [];
    final_x_sol = [final_x_sol; x_sol];

//    final_xm_leaks(validate_index,:) = [];
//    final_xm_leaks = [final_xm_leaks; xm_leaks];
//
//    final_x_sol_leak(validate_index,:) = [];
//    final_x_sol_leak = [final_x_sol_leak; x_sol_leak];



    jj_counter = jj_counter +1;

//    disp('exiting loop  ',  i)
end
varargout(1) = final_xfinal;
varargout(2) = final_x_sol;
//varargout(3) = final_xm_leaks ;
//varargout(4) = final_x_sol_leak ;


endfunction

