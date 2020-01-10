// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
function [resRand, resGrossErrorNodalRand, varargout]=generate_leaks(xr, sd, jac, runsize, lbres, ubres)
jac_row = size(jac,1);
xrs = zeros(szx,runsize);
rerror1 = zeros(runsize,szx);
grerrornodal = zeros(runsize*jac_row, jac_row);
resGrossErrorNodalRand = zeros(runsize*jac_row, jac_row);
//rerror=grand(runsize,szx,'nor',0,1);
// random number generators: rerror1 prefered as rerror
for i=1:szx
    rerror1(:,i)=grand(runsize,1,'nor',0,sd(i));        
//the above line must be commeted and the beelow line uncommented if the user wants to test just one measurement 
//    rerror1(:,i)=0;
end

// adding random error to exact x
for i=1:runsize
//    xrs(:,i) = xr + sd.*rerror(i,:)';
    xrs(:,i) = xr + rerror1(i,:)';
end

// gross errors for leakings (Nodal Test)

// the total flow that passes through the node (streams entering or leaving the node)
totalNodeFlow = abs(jac)*xr;

//gross error to nodes (leaking)
k = 0;
for i=1:jac_row
    for j=1:runsize
        grerrornodal(j+k*runsize,i) = grand(1,1,'unf',totalNodeFlow(i)*lbres,totalNodeFlow(i)*ubres);
//        leaks(j+k*runsize) = grerrornodal(j+k*runsize,i);
    end
    k=k+1;
end
// residuals when pure random noise is added to measurements
resRand = zeros(runsize,jac_row);
for i = 1: runsize
    resRand(i,:) = (jac*(xrs(:,i)))';
end
//sum of leaking plus randon error residuals
//Adding gross error plus random error

for i=0:(size(grerrornodal,1)-1)
            resGrossErrorNodalRand(i+1,1:jac_row) = (resRand(modulo(i,runsize)+1,:) -grerrornodal(i+1,1:jac_row) );

end

varargout = list(grerrornodal);
//disp('before leave generate data');
//pause
endfunction


