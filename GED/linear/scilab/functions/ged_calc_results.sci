// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

function [x_sol] = calc_results_DR(xfinal, jac, sigma, resGrossErrorNodalRandFi, opt_type)
obj_function_type = opt_type;
exec ../functions/setup_DR.sce

// to run robust reconciliation, it is also necessary to choose the function to return the problem structure
if obj_function_type > 0 then
    xm = xr;
    [nc_eq, n_non_lin_eq, nv, nnzjac_ineq, nnzjac_eq, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = robust_structure(jac, 0, xr, objfun, res_eq, res_ineq);
    //opt_type = 1; // linear
    opt_type = 2; // nonlinear : robust functions or nonlinear constraits
    x_sol = calc_results_opt(xfinal, jac, sigma, resGrossErrorNodalRandFi, opt_type);

else // WLS
    if obj_function_type < 0 then // WLS analytical
        x_sol = calc_results_analytical(xfinal, jac, sigma, resGrossErrorNodalRandFi);
    else // WLS numeric
        [nc, nv, i1, i2, nnzeros, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = wls_structure(jac);
        opt_type = 1; // linear
        //opt_type = 2; // nonlinear : robust functions or nonlinear constraits
        x_sol = calc_results_opt(xfinal, jac, sigma, resGrossErrorNodalRandFi, opt_type);
    end
end

    
endfunction    
    

function [x_sol]=calc_results_analytical(xfinal, jac, sigma, resGrossErrorNodalRandFi, varargin)
[lhs,rhs]=argn(0);
V=jac*sigma*jac';
V_inv = inv(V);
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
BB=(eye(szx,szx)-sigma*jac'*V_inv*jac);
for i=1:runsizefinal
//    xrs(:,i) = xr + sd.*rerror(i,:)';
    x_sol(:,i) = BB*xfinal(i,:)';
//    [x_sol(i,:), f_sol(i), status(i)] = P5(xfinal(i,:)',sds.^2,xr');
end

x_sol = x_sol';
endfunction


function [x_sol]=calc_results_opt(xfinal, jac, sigma, resGrossErrorNodalRandFi, opt_type)
[lhs,rhs]=argn(0);
rj=rank(jac);
jac_col = size(jac,2);
jac_row = size(jac,1);

runsizefinal = size(xfinal,1);
x_sol= zeros(szx,runsizefinal);
f_sol=zeros(runsizefinal);
params = init_param();
// We use the given Hessian
params = add_param(params,"hessian_approximation","exact");
if  opt_type == 1 then
    params = add_param(params,"hessian_constant","yes");
end
params = add_param(params,"tol",1e-8);
params = add_param(params,"acceptable_tol",1e-8);
params = add_param(params,"mu_strategy","adaptive");
params = add_param(params,"journal_level",0);
params = add_param(params,"print_level",0);



for i=1:runsizefinal
//    xrs(:,i) = xr + sd.*rerror(i,:)';
//    x_sol(:,i) = BB*xfinal(i,:)';
    xm = xfinal(i,:)';
    [x_sol(:,i), f_sol(i), extra] = ipopt(xfinal(i,:)', objfun, gradf, confun, dg, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);
    
//    [x_sol(i,:), f_sol(i), status(i)] = P5(xfinal(i,:)',sds.^2,xr');
//disp(i);
end
x_sol = x_sol';

endfunction

function [res, gamaMeasuremts,gamaNodal,zr_nt_nodal, zr_nt_nodal_rand, zadj, varargout ]=calc_results_index(x_sol, jac, sigma, resGrossErrorNodalRandFi, varargin)
[lhs,rhs]=argn(0);
V=jac*sigma*jac';
V_inv = inv(V);
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

if rhs > 4 then 
    x_sol_leak= varargin(1);
    adj_nodal = zeros(runsizefinal,szx);
    zadj_nodal = zeros(runsizefinal,szx);
    adj_leak = xfinal - x_sol_leak;

end
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
// for maximum power constraint test (Narasimhan pg. 180 eq 7.9), comment the line above and uncomment the 2 lines bellow
//pause
//            V_inv_res = (V_inv*resGrossErrorNodalRandFi(i,:)')';
//            zr_nt_nodal(i,k)=abs(V_inv_res(k))./(diag_diag_inv_V(k,k).^0.5);  
            if i <= runsize then
                 zr_nt_nodal_rand(i,k)=abs(res(i,k))./(diag_diag_V(k,k).^0.5); 
// for maximum power constraint test (Narasimhan pg. 180 eq 7.9), comment the line above and uncomment the 2 lines bellow
//                 V_inv_res = (V_inv*res(i,:)')';
//                 zr_nt_nodal_rand(i,k)=abs(V_inv_res(k))./(diag_diag_inv_V(k,k).^0.5);                    
            end //if
         end
    end //for
    for j=1:szx
// measurements test statistics: narasimham pg. 183 eq. 7-14        
        zadj(i,j)=abs(adj(i,j))/sqrt(Wbar(j,j));
        if rhs > 4 then
            zadj_nodal(i,j)=abs(adj_nodal(i,j))/sqrt(Wbar(j,j));
        end
    end
end

if lhs > 7 then
    varargout = zadj_nodal;
end

endfunction
