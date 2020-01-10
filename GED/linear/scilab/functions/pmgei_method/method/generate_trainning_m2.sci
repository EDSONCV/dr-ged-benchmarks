function [final_train, final_validate] = generate_train_m2(xr, xfinal, sigma, obj_function_type, runsize, sd, jac, gt_mt_tol, gt_nt_tol, mt_tol_1, mt_tol_2, nt_tol_1, nt_tol_2, mult_bias_low, mult_bias_inc, multi_bias_up,multi_leak_low, multi_leak_inc, multi_leak_up ,bias_vector,leak_vector,add_errors, input_train,input_validate, n_validate)
                  
szx = size(xr,1);

sds = sd;
sds =sd;
jac_col = size(jac,2);
jac_row = size(jac,1);
rj=rank(jac);
sigma=diag(sds.^2);
var =  sd.^2;
jj = [multi_leak_low: multi_leak_inc: multi_leak_up];
jj_counter = 1;
validate =[];
final_train=[];
final_validate =[];
final_xfinal=[];
final_x_sol =[];
final_xm_leaks = [];
final_x_sol_leak = [];

length_merrornode = length(find(leak_vector > 0));
ind_merrornode = find(leak_vector > 0);

length_merrorbias = length(find(bias_vector <> 0));
ind_merrorbias = find(bias_vector <> 0);

ind_error_all = [ind_merrorbias, ind_merrornode + jac_col];

for ii = mult_bias_low: mult_bias_inc: multi_bias_up 

    [xfinal2, resRand2, resGrossErrorNodalRand2, grerror, mySign,leaks ]=generate_data_multiple(xfinal(1:runsize,:), sd, jac, ii,ii, multi_leak_low, multi_leak_up ,bias_vector,leak_vector,add_errors);


    [xm_leaks] = generate_meas_from_leaks(xr, resRand2, leaks, jac, runsize, 1, xfinal2(runsize + 1: $, :));

    resGrossErrorNodalR2 = [ resRand2;resGrossErrorNodalRand2];

    is_multiple = 1;

    [x_sol2, res2, gamaMeasuremts2,gamaNodal2,zr_nt_nodal2, zr_nt_nodal_mt2, zr_nt_nodal_rand2, zadj2, zadj_nodal2, x_sol_leak2] = calc_results_tst_pls_sig(xfinal2, jac, sigma, resGrossErrorNodalR2, obj_function_type, is_multiple, leaks, xm_leaks);
disp('after calc_results_tst_pls_sig')    

    [avti_gt_mt, avti_gt_nt, op_gt_mt, op_gt_nt, ge_gt_delete_mt ge_gt_delete_nt, ge_gt_delete_low] = global_test2(gt_mt_tol, gt_nt_tol, abs(gamaMeasuremts2), abs(gamaNodal2), runsize, rj, jac_col, jac_row);

    // the same as measurement_test, but also returns the indices of elements below the test statistics    
    [avt1_mt1, avt1_mt2, op_mt1, op_mt2, ge_mt1_indexu, ge_mt2_indexu, ge_mt1_indexu_low, ge_mt2_indexu_low] = measurement_test3(mt_tol_1, mt_tol_2, abs(zadj2), runsize, jac_col, 1);
    // the same as nodal_test, but also returns the indices of elements below the test statistics    
    [avt1_nt1, avt1_nt2, ge_nt1_indexu, ge_nt2_indexu, ge_nt1_indexu_low, ge_nt2_indexu_low] = nodal_test2_mult(nt_tol_1, nt_tol_2, jac_row, runsize, abs(zr_nt_nodal2));


    index_to_clear_low = [];    
    index_to_clear_mt = [];    
    index_to_clear_nt = [];
    // remove unused components (measurement errors which random error were higher 
    //than the test statistics for three tests)    
    if ii == mult_bias_low then
        index_to_clear_low = unique([ge_gt_delete_low, ge_mt1_indexu_low, ge_nt1_indexu_low]);
    else
        index_to_clear_low = unique([[1:runsize],ge_gt_delete_low, ge_mt1_indexu_low, ge_nt1_indexu_low]);

    end

    // remove unused components (measurement errors which gross error with lower
    // than the test statistics) for measurement test   
    ge_mt1_indexu_filter =  setdiff([runsize + 1: 2* runsize],ge_mt1_indexu);
    index_to_clear_mt = unique([ge_mt1_indexu_filter, ge_gt_delete_mt]);
    // remove unused components (leakings which gross error with lower
    // than the test statistics) for nodal test   
    ge_nt1_indexu_filter =  setdiff([runsize + 1: 2* runsize],ge_nt1_indexu);
    index_to_clear_nt = unique([ge_nt1_indexu_filter, ge_gt_delete_nt]);


    // clear gamaMeasurements vector (gamma for measurement bias dataset)
    gamaMeasuremts_new = gamaMeasuremts2;
    gamaMeasuremts_new([index_to_clear_low, index_to_clear_mt]) = [];
    // clear zadj vector    
    zadj_new = zadj2;
    //    zadj_new = xfinal-x_sol;

    zadj_new([index_to_clear_low, index_to_clear_mt],:) = [];



    // clear zr_nt_nodal_mt (residuals of the measurement bias dataset) vector        
    zr_nt_nodal_mt_new = zr_nt_nodal_mt2;
    zr_nt_nodal_mt_new([index_to_clear_low, index_to_clear_mt], :) = [];

    xfinal2([index_to_clear_low, index_to_clear_mt], :) = [];
    x_sol2([index_to_clear_low, index_to_clear_mt], :) = [];

    // clear gamaNodal (gamma for leaking dataset)
    gamaNodal_new = gamaNodal2;    
    gamaNodal_new([index_to_clear_nt]) = [];
    // since gammaNodal(1:runsize) is already included in gamaMeasurements, the
    // first "runsize" elements will be removed
    gamaNodal_new(1:runsize) = [];
    // clear zr_nt_nodal (residuals of the leakings dataset) vector        
    //pause
    zr_nt_nodal_new = zr_nt_nodal2;
    //    zr_nt_nodal_new = resGrossErrorNodalRandFi;
    zr_nt_nodal_new(index_to_clear_nt, :) = [];
    // since zr_nt_nodal(1:runsize,:) is already included in zr_nt_nodal_mt, the
    // first "runsize" elements will be removed
    zr_nt_nodal_new(1:runsize,:) = [];
    // clear zadj_nodal (adjustments statistics of the leakings dataset) vector        
    //    zadj_nodal_new = zadj_nodal;

    //    zadj_nodal_new = [zadj_new(1:runsize,:);xm_leaks - x_sol_leak];

    zadj_nodal_new = [zadj_nodal2(runsize+1:$,:)];
//        zadj_nodal_new = [zadj_new(1:runsize,:);zadj_nodal2(runsize+1:$,:)];
    zadj_nodal_new(index_to_clear_nt,:) = [];
    // since zadj_nodal(1:runsize,:) is already included in zadj, the
    // first "runsize" elements will be removed

//    zadj_nodal_new(1:runsize,:) = [];
//disp('before y mt moutn');
//pause

    // arranging MT trainning data

   
        // arranging NT trainning data

        y_nt = zeros(runsize, jac_row);

 if length_merrornode > 0 then
     
        // we multiply for 10 to present a good scalling

        y2 = 10*jj(jj_counter)*ones(runsize ,1);

        for i = 1 : length_merrornode

            y_nt(: , ind_merrornode(i) ) = y2;
        end

        //        y_nt(1:runsize, :) = [];
        //    disp('before xmleaks');
       //pause

    end

//pause
    y = [grerror, y_nt];
    
    ytop = zeros(runsize,jac_col+ jac_row);
    
    y=[ytop;y];

    y([index_to_clear_low, index_to_clear_mt],:) = [];

    train=[gamaMeasuremts_new, zadj_new, zr_nt_nodal_mt_new, y];
    
    
    xfinal = [xfinal;xm_leaks];
    x_sol = [x_sol2;x_sol_leak2];

    validate_index=[];

    // pick some lines for validation
    for i = 0: length(ind_error_all)
//                disp('i: ', i)
//        disp('inside generate train_m2')
//        pause
        if i == 0 then
            // FIXME ALWAYS PICK 10 FIRST POINTS TO VALIDATION OF RANDOM NUMBERS
//            index = find(train(:,jac_row + jac_col + 2 + i) == 0); 
            index = [1:10]; 
        elseif i <= jac_col then

            index = find(train(:,jac_row + jac_col + 2 + ind_error_all(i) - 1) <> 0); 
        else
            index = find(train(:,jac_row + jac_col + 2 + ind_error_all(i) - 1) == 10*jj(jj_counter)); 
//            index = find(train(:,jac_row + jac_col + 2 + i - 1) == 10*jj(jj_counter)); 
        end
        // select the first "n_validate" for validation
        //       pause
        if length(index) == 0 then

        elseif length(index) <  n_validate then           
            validate_index = [validate_index, index(1:$)];    
        else
            validate_index = [validate_index, index(1:n_validate)];
        end
//        pause
        index =[];
    end

    //add the selected data to validation set
    final_validate = [final_validate; train(validate_index,:)];
    // remove the validation set 
    train(validate_index,:) = [];

    final_train = [final_train; train];

    xfinal2(validate_index,:) = [];
    final_xfinal = [final_xfinal; xfinal2];

    x_sol2(validate_index,:) = [];
    final_x_sol = [final_x_sol; x_sol2];

    //    final_xm_leaks(validate_index,:) = [];
    //    final_xm_leaks = [final_xm_leaks; xm_leaks];
    //
    //    final_x_sol_leak(validate_index,:) = [];
    //    final_x_sol_leak = [final_x_sol_leak; x_sol_leak];



    jj_counter = jj_counter +1;

    //    disp('exiting loop  ',  i)
end
//varargout(1) = final_xfinal;
//varargout(2) = final_x_sol;
//varargout(3) = final_xm_leaks ;
//varargout(4) = final_x_sol_leak ;


endfunction
