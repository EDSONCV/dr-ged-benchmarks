function [model,models_stat] = generate_pls_models_m( string_problem, nofstreams, nofeqps, train, validate, ndatarand,ndatainterval)
//    streams model
models_stat = list();
total_train = nofstreams + nofeqps + 1; // total number of inputs
length_train = size(validate,'r');  
total_validate = total_train + nofstreams + nofeqps;  //total number of variables (input + outputs)
//pause
 // create the PLS model with the 'train' dataset the inputs are: 
 // train(:,1:total_train) and the outputs are:
 //train(:,total_train+1:$)
model = qpls(train(:,1:total_train),train(:,total_train+1:$),total_train,1,0); // create the PLS model with the 'train' dataset the inputs are: 
//pause
// simulate the outputs from the generated mode, 'model', using the 'validate' dataset as input (validate(:,1:total_train) 
//and compare with the outputs, (validate(:,(total_train + 1):total_validate))
[Y0hat, pred_error,abs_pred_error,mean_absolute_error,imax_error,max_error] =predictq(validate(:,1:total_train),model,model.stat.bestdirs, validate(:,(total_train + 1):total_validate));
//Yhat holds the predictions
//find the max value  in the  random error dataset which belongs to the measured streams
find_max_rand_meas = max(abs(Y0hat(1:ndatarand,1:nofstreams )));
//find the max value  in the  random error dataset which belongs to the leakings. They are separeted because of the magnitude of the random errorsi
find_max_rand_eqp = max(abs(Y0hat(1:ndatarand,nofstreams+1:$ )));
printf('\n find_max_rand_meas: %f', find_max_rand_meas)
printf('\n find_max_rand_eqp: %f \n' , find_max_rand_eqp)
    for i = (total_train + 1):total_validate
       // pause
        model_s = struct();
        // prediction
        Y0hat_int = Y0hat(:,i - total_train);
        predict_model.Y0hat = Y0hat_int;
        Y0_int =  validate(:,i);
        predict_model.Y0 = Y0_int;
        predict_model.pred_error = Y0hat_int - Y0_int;
        predict_model.abs_pred_error = abs(Y0hat_int - Y0_int);
        predict_model.mean_absolute_error = sum(abs((Y0hat_int - Y0_int)))/length_train;
        [ predict_model.max_error,  predict_model.imax_error] = max(abs(Y0hat_int - Y0_int));
        // FIXME - CHECK FOR MULTIPLE ERRORS
//        find_max_rand = max(abs(validate(1:ndatatrain,i )));
//pause        
        [find_simulated_error_i,find_simulated_error_j] = find(validate(:,i ) <> 0);
        [not_simulated_error_i,not_simulated_error_j]  = find(validate(:,i ) == 0);
        find_simulated_error_j = unique(find_simulated_error_j);
        not_simulated_error_j = unique (not_simulated_error_j);
        min_lower = min(abs(Y0hat_int(find_simulated_error_i,find_simulated_error_j)));


//disp('         inside generate_pls_models:   ')     
//disp(i);


//        [find_upper_i,find_upper_j] = find(abs(Y0hat_int(find_simulated_error_i,find_simulated_error_j)) >= min_lower) ;
//  if i == 73 then
//      pause
//  else
      
//  end
        
        if length(find_simulated_error_i) > 0 then
            if( i - (total_train) <= nofstreams  ) then
                [find_upper_i,find_upper_j] = find(abs(Y0hat_int(find_simulated_error_i,find_simulated_error_j)) >= find_max_rand_meas) ;
            else
                [find_upper_i,find_upper_j] = find(abs(Y0hat_int(find_simulated_error_i,find_simulated_error_j)) >= find_max_rand_eqp) ;
            end
            find_upper = length(find_upper_i )/length(find_simulated_error_i);
//            aee = sum(abs(validate(find_simulated_error_i,i)) - abs(Y0hat_int(find_simulated_error_i,find_simulated_error_j)))/length(find_simulated_error_i);
            aee = 100*abs(sum((validate(find_simulated_error_i,i) -Y0hat_int(find_simulated_error_i,find_simulated_error_j))./validate(find_simulated_error_i,i))/length(find_simulated_error_i));

       else     
            find_upper = 0;
            aee = -666;
        end


//  pause                
        if length(not_simulated_error_i) > 0 & length(min_lower) > 0 then

           [find_lower_i,find_lower_j] = find(Y0hat_int(not_simulated_error_i,not_simulated_error_j) >= min_lower)
            avti = length(find(abs(validate(1:ndatarand,i ))) >= min_lower);    
            find_lower = length(find_lower_i)/length(not_simulated_error_i);
        
        else
        
        find_lower = 1;
        avti = -6666    
        end
        
//        disp('before aee')
//        pause

//
        predict_model.false_detection = find_lower;
        predict_model.selectivity = 1 - find_lower;
        predict_model.correct_detection = find_upper;
        predict_model.avti = avti;
        predict_model.aee = aee;
//        disp(i-total_train);
        model_s.prediction = predict_model;
        if (i - total_train) <= nofstreams then
            model_s.name = ('model_'  + string_problem + '_stream_' + string(i - total_train));
        else
            model_s.name = ('model_'  + string_problem + '_leak_' + string(i - total_train - nofstreams));
        end
        models_stat(i - total_train) = model_s;
    end
endfunction
