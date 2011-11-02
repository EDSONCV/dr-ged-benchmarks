// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
function [x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj ]=calc_results(xfinal, jac, sigma, resGrossErrorNodalRandFi)
V=jac*sigma*jac';
V_inv = inv(V);
// covariance matrix of adjustments: narasimham pg. 183 eq. 7-13
Wbar=sigma*jac'*inv(V)*jac*sigma;
// variance-covariance matrix: narasimham pg. 178 eq. 7-3
diag_diag_V = diag(diag(V));
sigma_inv=inv(sigma);
rj=rank(jac);
jac_col = size(jac,2);
jac_row = size(jac,1);

runsizefinal = size(xfinal,1);
x_sol= zeros(szx,runsizefinal);
BB=(eye(szx,szx)-sigma*jac'*V_inv*jac);
for i=1:runsizefinal
//    xrs(:,i) = xr + sd.*rerror(i,:)';
    x_sol(:,i) = BB*xfinal(i,:)';
//    [x_sol(i,:), f_sol(i), status(i)] = P5(xfinal(i,:)',sds.^2,xr');
end

x_sol = x_sol';

res = zeros(runsizefinal, jac_row);
gamaMeasuremts = zeros(runsizefinal);
gamaNodal = zeros(runsizefinal);
zr_nt_measurement = zeros(runsizefinal, jac_row);
zr_nt_nodal = zeros(runsize, jac_row);
zr_nt_nodal_rand = zeros(runsize, jac_row);
adj = zeros(runsizefinal,szx);
zadj = zeros(runsizefinal,szx);
d_adj = zeros(szx,runsizefinal);
runsizenodal = size(resGrossErrorNodalRandFi,1);
//adjustments  narasimham pg. 183 eq. 7-11     
 adj = xfinal - x_sol;  
for i=1:runsizefinal
// residuals when random noise and gross errors are
// added to the measurements: narasimham pg. 178 eq. 7-2    
// residuals of res(1:runsize,:) is related with measurements
// randon noise, while res(runsize + 1: end,:) to gross errors
    res(i,1:jac_row)=(jac*xfinal(i,:)')';
// global test statistics: narasimham pg. 178 eq. 7-4       
// gamaMeasuremts is used in global test
// gamaMeasuremts(1:runsize) is is related with measurements
// randon noise, while gamaMeasuremts(runsize + 1: end) to gross errors
    gamaMeasuremts(i) = res(i,:)*V_inv*res(i,:)';
// nodal test statistics: narasimham pg. 180 eq. 7-5
// the residuals of resGrossErrorNodalRand has a different structure of 'res'
// each row, 'i', has a gross error in balance 'i', added 'runsize' times.
    for k=1:jac_row
        if i <=  runsizenodal then
// gamaNodal is calculated as  gamaMeasuremts but the residuals used are the
// residuals from gross errors in the balances (leakings)          
            gamaNodal(i) = resGrossErrorNodalRandFi(i,:)*V_inv*resGrossErrorNodalRandFi(i,:)';
// zr_nt_nodal is calculated as Narasimhan pg. 180 eq 7.5 also but the residuals used are the
// residuals from gross errors in the balances (leakings)                  
            zr_nt_nodal(i,k)=abs(resGrossErrorNodalRandFi(i,k))./(diag_diag_V(k,k).^0.5);  
            
            if i <= runsize then
                zr_nt_nodal_rand(i,k)=abs(res(i,k))./(diag_diag_V(k,k).^0.5);                    
            end //if
         end
    end //for
    for j=1:szx
// measurements test statistics: narasimham pg. 183 eq. 7-14        
        zadj(i,j)=abs(adj(i,j))/sqrt(Wbar(j,j));
    end
end

endfunction