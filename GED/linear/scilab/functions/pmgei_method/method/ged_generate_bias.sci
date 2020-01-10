// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
function [xfinal]=generate_bias(xr, sd, jac, runsize, lbm, ubm)
jac_row = size(jac,1);
xrs = zeros(szx,runsize);
rerror1 = zeros(runsize,szx);
grerror = zeros(szx*runsize,szx);
grerrors = zeros(szx*runsize,szx);

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
// gross errors for measurements (Measurement Test)
// random sign generator to add to gross error
// for some tests, like the MT it doesn't matter because the absolute value of the adjustments are considered
mySign=sign(grand(runsize,szx,'unf',-1,1));

k=0;
for i=1:szx
    for j=1:runsize
//        grerror(j+k*runsize,i) = grand(1,1,'unf',xrs(i,j)*0.05,xrs(i,j)*0.1)*mySign(j,i);
        grerror(j+k*runsize,i) = grand(1,1,'unf',lbm*sd(i),ubm*sd(i))*mySign(j,i);
    end
    k=k+1;
end
//Adding gross error plus random error
for i=0:(size(grerror,1)-1)
    grerrors(i+1,1:szx) = grerror(i+1,1:szx) + xrs(:,modulo(i,runsize)+1)';
end
// xfinal concatenates pure random measurement in sensors with gross errors in sensors
xfinal =[xrs';grerrors]; 

//disp('before leave generate data');
//pause
endfunction


