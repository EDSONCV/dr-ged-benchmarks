// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
function [x_sol, res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_mt, zr_nt_nodal_rand, zadj, varargout] = calc_results_tst_pls_sig(xfinal, jac, sigma, resGrossErrorNodalRandFi ,obj_function_type, is_multiple, varargin)
[lhs,rhs]=argn(0);
if is_multiple == 0  then
    runsize =  size(xfinal,1)/(size(jac,2) + 1);
else
    runsize = size(xfinal,1)/2;
end

if rhs > 6 then
    leak = varargin(1);
    xfinal_leak = varargin(2);
    [x_sol_leak, status]  = calc_results_leak(xfinal_leak, jac, leak, runsize);
end

V=jac*sigma*jac';
V_inv = inv(V);
//disp('inv-V',V_inv)
// covariance matrix of adjustments: narasimham pg. 183 eq. 7-13
Wbar=sigma*jac'*inv(V)*jac*sigma;

// variance-covariance matrix: narasimham pg. 178 eq. 7-3
diag_diag_V = diag(diag(V));
diag_diag_inv_V = diag(diag(V_inv));

sigma_inv=inv(sigma);
rj=rank(jac);
jac_col = size(jac,2);
jac_row = size(jac,1);

runsizefinal = size(xfinal,1);
x_sol= zeros(szx,runsizefinal);

[x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, obj_function_type);
//sum leak + measurement bias for x_sol also
if is_multiple == 1 then
     x_sol(runsize + 1: $,:)   = x_sol_leak;
end
    
//disp('inside ged_calc_results_pls')
//pause
//BB=(eye(szx,szx)-sigma*jac'*V_inv*jac);
//for i=1:runsizefinal
////    xrs(:,i) = xr + sd.*rerror(i,:)';
//    x_sol(:,i) = BB*xfinal(i,:)';
////    [x_sol(i,:), f_sol(i), status(i)] = P5(xfinal(i,:)',sds.^2,xr');
//endgamaMeasuremts2,gamaNodal2,zr_nt_nodal2, zr_nt_nodal_mt2, zr_nt_nodal_rand2, zadj2, zadj_nodal2
//pause
//x_sol = x_sol';

res = zeros(runsizefinal, jac_row);
runsizenodal = size(resGrossErrorNodalRandFi,1);
gamaMeasuremts = zeros(runsizefinal);
gamaNodal = zeros(runsizenodal);
zr_nt_measurement = zeros(runsizefinal, jac_row);
zr_nt_nodal = zeros(runsizenodal, jac_row);
zr_nt_nodal_mt = zeros(runsizefinal, jac_row);
zr_nt_nodal_rand = zeros(runsize, jac_row);
adj = zeros(runsizefinal,szx);
zadj = zeros(runsizefinal,szx);
d_adj = zeros(szx,runsizefinal);

//adjustments  narasimham pg. 183 eq. 7-11     

 adj = xfinal - x_sol;
 


if rhs > 5 then 
    //    adj_nodal = zeros(jac_row*runsize,szx);
    if is_multiple == 1 then
        zadj_nodal = zeros(runsize,szx);
    else
        zadj_nodal = zeros(jac_row*runsize,szx);
    end
    adj_nodal = xfinal_leak - x_sol_leak;
end
//pause
disp('befor zad calcl');
//pause
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
//originals
//            zr_nt_nodal(i,k)=abs(resGrossErrorNodalRandFi(i,k))./(diag_diag_V(k,k).^0.5);  
//            zr_nt_nodal_mt(i,k)=abs(res(i,k))./(diag_diag_V(k,k).^0.5);  
//without abs
           
           // pause
            zr_nt_nodal(i,k)=(resGrossErrorNodalRandFi(i,k))./(diag_diag_V(k,k).^0.5);  
            zr_nt_nodal_mt(i,k)=(res(i,k))./(diag_diag_V(k,k).^0.5);  
// for maximum power constraint test (Narasimhan pg. 180 eq 7.9), comment the line above and uncomment the 2 lines bellow
//pause
//            V_inv_res = (V_inv*resGrossErrorNodalRandFi(i,:)')';
//            zr_nt_nodal(i,k)=abs(V_inv_res(k))./(diag_diag_inv_V(k,k).^0.5);  
            if i <= runsize then
// original
// zr_nt_nodal_rand(i,k)=abs(res(i,k))./(diag_diag_V(k,k).^0.5);                 
//without abs
                 zr_nt_nodal_rand(i,k)=(res(i,k))./(diag_diag_V(k,k).^0.5);                  
// for maximum power constraint test (Narasimhan pg. 180 eq 7.9), comment the line above and uncomment the 2 lines bellow
//                 V_inv_res = (V_inv*res(i,:)')';
//                 zr_nt_nodal_rand(i,k)=abs(V_inv_res(k))./(diag_diag_inv_V(k,k).^0.5);                    
            end //if
        else
            //original
//         zr_nt_nodal_mt(i,k)=abs(res(i,k))./(diag_diag_V(k,k).^0.5);  
        //without abs
         zr_nt_nodal_mt(i,k)=(res(i,k))./(diag_diag_V(k,k).^0.5);  
         end
    end //for
//    disp('inside calc_results_tst_sig')
//    pause
    for j=1:szx
// measurements test statistics: narasimham pg. 183 eq. 7-14  
// original
//        zadj(i,j)=abs(adj(i,j))/sqrt(Wbar(j,j));   
// without abs        
        zadj(i,j)=(adj(i,j))/sqrt(Wbar(j,j));
// zadj_nodal has the size of runsize*jac_row, so this if below is necessary
         if rhs > 5 & i <= runsize*jac_row then
//            printf('i: %d j: %d \n', i, j ) 
//orignial
//            zadj_nodal(i,j)=abs(adj_nodal(i,j))/sqrt(Wbar(j,j));
//withou abs
            if is_multiple ==1 then
                if i <=runsize then
                     zadj_nodal(i,j)=(adj_nodal(i,j))/sqrt(Wbar(j,j));
                end

            else
                zadj_nodal(i,j)=(adj_nodal(i,j))/sqrt(Wbar(j,j));
            end
        
            
        end
        
    end
end
//
disp('inside calc_results_tst_pls_sig')
//pause
if lhs > 8 then
    varargout(1) = [zadj(1:runsize,:);zadj_nodal];
    varargout(2) = x_sol_leak;
end
endfunction
