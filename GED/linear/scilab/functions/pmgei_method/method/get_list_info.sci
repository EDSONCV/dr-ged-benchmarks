function [avti_meas, op_meas, selectivity_meas, aee_meas, avti_eqp, op_eqp, selectivity_eqp, aee_eqp] = get_lit_info(list_data, n_streams, n_equipments)
    
avti_meas = []; op_meas= []; selectivity_meas= []; aee_meas= []; avti_eqp= []; op_eqp= []; selectivity_eqp= []; aee_eqp = [];    
for i = 1: n_streams
    avti_meas(i,1) = list_data(i).prediction.avti;
    op_meas(i,1) = list_data(i).prediction.correct_detection;
    selectivity_meas(i,1) = list_data(i).prediction.selectivity;
    aee_meas(i,1) = list_data(i).prediction.aee;
end    
       
for i = 1 + n_streams: n_equipments + n_streams;
    avti_eqp(i - n_streams,1) = list_data(i).prediction.avti;
    op_eqp(i - n_streams,1) = list_data(i).prediction.correct_detection;
    selectivity_eqp(i - n_streams,1) = list_data(i).prediction.selectivity;
    aee_eqp(i - n_streams,1) = list_data(i).prediction.aee;
end        
//varargout(1) = list_data.nconverged;
//varargout(2) = list_data.nnotconverged;
    
endfunction
