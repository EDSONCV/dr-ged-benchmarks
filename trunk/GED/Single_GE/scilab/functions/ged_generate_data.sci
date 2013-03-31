// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// we have 4 functions in this file:
// -generate_data
// -generate_data_multiple
// -generate_data_random_err
// -generate_data_errors
function [xfinal, resRand, resGrossErrorNodalRand, varargout]=generate_data(xr, sd, jac, runsize, lbm, ubm, lbres, ubres)
// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// generate single gross error vector according to user input
//***************************************************************
//This function receives the users argument and generate a single gross error for measurement bias and leakings. 

//Outputs:
//    xfinal: a vector with [xrs; gross errors in measurements]
//    resRand: residuals when pure random error was added
//    resGrossErrorNodalRand: residuals with pure random error - leaking
//    varargout(1) = grerrornodal =  only the leaking
//
//Inputs:

//    xrs:                measurements with pure random error
//    sd:                 standard deviations of measurements
//    jac                 Jacobian  matrix
//    lbm               lower bound of measurement error
//    ubm               upper bound of measurement error
//    lbres               lower bound residuals leaking
//    ubres               upper bound residuals leaking
  
    
    
jac_row = size(jac,1);
xrs = zeros(szx,runsize);
rerror1 = zeros(runsize,szx);
grerror = zeros(szx*runsize,szx);
grerrors = zeros(szx*runsize,szx);
grerrornodal = zeros(runsize*jac_row, jac_row);
//leaks = zeros(runsize*jac_row);
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

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];
varargout(1) = list(grerrornodal);

//disp('before leave generate data');
//pause

endfunction

function [xfinal, resRand, varargout]=generate_data_multiple(xrs, sd, jac, lbm, ubm, lbres, ubres,vec_bias_error,vec_nodal_error, flag_sum)

// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// generate multiple gross error vector according to user input
//***************************************************************
//This function receives the users argument and generate MULTIPLE gross error for measurement bias and leakings. 
//User have the option to sum these errors at the end for the residuals
//Outputs:
//    xfinal: a vector with [xrs; gross errors in measurements]
//    resRand: residuals when pure random error was added
//    varargout(1) = resGrossErrorNodalRand = residuals with pure random error - leaking
//    varargout(2) = grerrornodal =  only the leaking
//    varargout(3) = residuals calculated based on the measurement bias (with gross error) - grerrornodal 
//                 only calculated if "flag_sum" = 1
//
//Inputs:

//    xrs:                measurements with pure random error
//    sd:                 standard deviations of measurements
//    jac                 Jacobian  matrix
//    lbm               lower bound of measurement error
//    ubm               upper bound of measurement error
//    lbres               lower bound residuals leaking
//    ubres               upper bound residuals leaking
//    vec_bias_error      a gross error signature vector for measurement bias. A vector os size of xrs with one on streas
//                        that a gross error must be added and zero elsewere.
//                        ex: if we have 6 streams and want a gross error in streams 2 and 5:
//                        vec_bias_error = [0 1 0 0 1 0]    
//    vec_nodal_error     a gross error signature vector for leakings. A vector os size of jac_row (number of residuals) 
//                        with one on residuals that an error will be added and zero elsewere.
//                        ex: if we have 4 balances (equipments) and want a gross error in equipment 2:
//                        vec_bias_error = [0 1 0 0]    
//    flag_sum        sum residuals of measurement bias and leakings

[lhs ,rhs]=argn();

runsize = size(xrs,1);
// finding the sizes
length_merrorbias = length(find(vec_bias_error > 0));
ind_merrorbias = find(vec_bias_error > 0);
length_merrornode = length(find(vec_nodal_error > 0));
ind_merrornode = find(vec_nodal_error > 0);

jac_row = size(jac,1);

grerror = zeros(runsize,szx);
grerrors = zeros(runsize,szx);
grerrornodal = zeros(runsize, jac_row);
mySign = zeros(runsize,1);

resGrossErrorNodalRand = zeros(runsize*jac_row, jac_row);

// gross errors for measurements (Measurement Test)

k=0;
// pure gross error
//pause
for i=1:length_merrorbias
     // random sign generator to add to gross error
    // for some tests, like the MT it doesn't matter because the absolute value of the adjustments are considered
    mySign=sign(grand(runsize,1,'unf',-1,1));
    for j=1:runsize
//            disp([i,j]);
            grerror(j,ind_merrorbias(i)) = grand(1,1,'unf',lbm*sd(ind_merrorbias(i)),ubm*sd(ind_merrorbias(i)))*mySign(j,1);
    end
end

//Adding gross error plus random error
grerrors = grerror + xrs;

// xfinal concatenates pure random measurement in sensors with gross errors in sensors
xfinal =[xrs;grerrors]; 
// gross errors for leakings (Nodal Test)

// the total flow that passes through the node (streams entering or leaving the node)
totalNodeFlow = abs(jac)*xr;
// the total flow that passes through the node divided by the total number of streams 
// entering or leaving it
meanNodeFlow = totalNodeFlow./sum(abs(jac),2);
//gross error to nodes (leaking)

k = 0;
for i=1:length_merrornode
    for j=1:runsize
        grerrornodal(j,ind_merrornode(i)) = grand(1,1,'unf',totalNodeFlow(ind_merrornode(i))*lbres,totalNodeFlow(ind_merrornode(i))*ubres);
//        leaks(j+k*runsize) = grerrornodal(j+k*runsize,i);
    end
end

// residuals when pure random noise is added to measurements
resRand = zeros(runsize,jac_row);
for i = 1: runsize
    resRand(i,:) = (jac*(xrs(i,:)'))';
end
if length_merrornode == 0 
    return
end

//sum of leaking plus randon error residuals
//Adding gross error plus random error

resGrossErrorNodalRand = resRand - grerrornodal;

resGrossErrorNodalRandFi = [ resRand;resGrossErrorNodalRand];
varargout(1) = (resGrossErrorNodalRand);
varargout(2) = (grerrornodal);
// sum the residuals
if flag_sum == 1 & lhs > 4 then
  sumres = zeros(runsize,jac_row); 
  
  for i = 1: runsize
    sumres(i,:) = (jac*(grerrors(i,:)'))' - grerrornodal(i,:) ;
  end  
    varargout(3) = sumres;    
end

//disp('before leave generate data');
//pause

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



