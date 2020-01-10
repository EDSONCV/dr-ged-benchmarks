function [aee, op, selectivity,avti] = get_list_info_mult(model_stat)

size_struct = length(model_stat.aee);
avti = model_stat.avti(1);

for i = 1: size_struct
    aee(i,1) =  model_stat.aee(i)
    op(i,1) =  model_stat.correct_detection(i);
    selectivity(i,1) =  model_stat.selectivity(i);
end
    
    
endfunction
