// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// This funcion calculates a Measurement corresponding to a leak in a certain node
function [xm_leaks] = generate_meas_from_leaks(xreal, resRand, leaks, jac, runsize, varargin)

[lhs, rhs] = argn(0);
if rhs > 5 then
    is_multiple = varargin(1);
    x_add = varargin(2)
else
    is_multiple = 0;
end
    
jac_row = size(jac,1);
jac_col = size(jac,2);
resRand_plus_leaks = ones(size(leaks,1),size(leaks,2));
xm_leaks = ones(size(leaks,1),jac_col);
//disp('inside generate_meas_from_leaks')
//pause
if is_multiple == 1  then
    for j=0:jac_row -1
//        resRand_plus_leaks = resRand - leaks;
        resRand_plus_leaks = + leaks;
    end
    
    for i=1:size(leaks,1)
    xm_leaks(i,:) = (pinv(jac)*resRand_plus_leaks(i,:)' + x_add(i,:)')';

 //   xm_leaks(i,:) = (sigma*jac'*inv(jac*sigma*jac')*resRand_plus_leaks(1,:)' + xr)'
    end

else
    for j=0:jac_row -1

        resRand_plus_leaks((j*runsize) +1 : runsize*(j+1),:) = resRand + leaks((j*runsize) +1 : runsize*(j+1),:);
    end
    
    for i=1:size(leaks,1)

    xm_leaks(i,:) = (pinv(jac)*resRand_plus_leaks(i,:)' + xr)';
    // the line above have the same results as
 //   xm_leaks(i,:) = (sigma*jac'*inv(jac*sigma*jac')*resRand_plus_leaks(1,:)' + xr)'
    end

end


endfunction