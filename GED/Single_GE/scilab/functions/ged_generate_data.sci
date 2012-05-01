// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
function [xfinal, resRand, resGrossErrorNodalRand]=generate_data(xr, sd, jac, runsize, lbm, ubm, lbres, ubres)
jac_row = size(jac,1);
xrs = zeros(szx,runsize);
rerror1 = zeros(runsize,szx);
grerror = zeros(szx*runsize,szx);
grerrors = zeros(szx*runsize,szx);
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
// gross errors for leakings (Nodal Test)

// the total flow that passes through the node (streams entering or leaving the node)
totalNodeFlow = abs(jac)*xr;
// the total flow that passes through the node divided by the total number of streams 
// entering or leaving it
meanNodeFlow = totalNodeFlow./sum(abs(jac),2);
//gross error to nodes (leaking)
k = 0;
for i=1:jac_row
    for j=1:runsize
        grerrornodal(j+k*runsize,i) = grand(1,1,'unf',totalNodeFlow(i)*lbres,totalNodeFlow(i)*ubres);
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

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

endfunction


// the function was divided to make a better implementation of OP curves generation
function [xrs] = generate_data_random_err(xr, sd, jac, runsize)
jac_row = size(jac,1);
xrs = zeros(szx,runsize);
rerror1 = zeros(runsize,szx);
grerror = zeros(szx*runsize,szx);
grerrors = zeros(szx*runsize,szx);
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
endfunction
function [xfinal, resRand, resGrossErrorNodalRand]=generate_data_errors(xr, xrandom, sd, jac, runsize, lbm, ubm, lbres, ubres)

jac_row = size(jac,1);

rerror1 = zeros(runsize,szx);
grerror = zeros(szx*runsize,szx);
grerrors = zeros(szx*runsize,szx);
grerrornodal = zeros(runsize*jac_row, jac_row);
resGrossErrorNodalRand = zeros(runsize*jac_row, jac_row);

// gross errors for measurements (Measurement Test)
// random sign generator to add to gross error
// for some tests, like the MT it doesn't matter because the absolute value of the adjustments are considered
mySign=sign(grand(runsize,szx,'unf',-1,1));

k=0;
for i=1:szx
    for j=1:runsize
//        grerror(j+k*runsize,i) = grand(1,1,'unf',xrandom(i,j)*0.05,xrandom(i,j)*0.1)*mySign(j,i);
        grerror(j+k*runsize,i) = grand(1,1,'unf',lbm*sd(i),ubm*sd(i))*mySign(j,i);
    end
    k=k+1;
end
//Adding gross error plus random error
for i=0:(size(grerror,1)-1)
    grerrors(i+1,1:szx) = grerror(i+1,1:szx) + xrandom(:,modulo(i,runsize)+1)';
end
// xfinal concatenates pure random measurement in sensors with gross errors in sensors
xfinal =[xrandom';grerrors]; 
// gross errors for leakings (Nodal Test)

// the total flow that passes through the node (streams entering or leaving the node)
totalNodeFlow = abs(jac)*xr;
// the total flow that passes through the node divided by the total number of streams 
// entering or leaving it
meanNodeFlow = totalNodeFlow./sum(abs(jac),2);
//gross error to nodes (leaking)
k = 0;
for i=1:jac_row
    for j=1:runsize
        grerrornodal(j+k*runsize,i) = grand(1,1,'unf',totalNodeFlow(i)*lbres,totalNodeFlow(i)*ubres);
    end
    k=k+1;
end
// residuals when pure random noise is added to measurements
resRand = zeros(runsize,jac_row);
for i = 1: runsize
    resRand(i,:) = (jac*(xrandom(:,i)))';
end
//sum of leaking plus randon error residuals
//Adding gross error plus random error

for i=0:(size(grerrornodal,1)-1)
            resGrossErrorNodalRand(i+1,1:jac_row) = (resRand(modulo(i,runsize)+1,:) -grerrornodal(i+1,1:jac_row) );

end

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];

endfunction



