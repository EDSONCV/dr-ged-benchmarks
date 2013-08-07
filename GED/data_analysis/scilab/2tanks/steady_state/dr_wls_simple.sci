// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
function [x_sol]=dr_wls_simple(xfinal, jac, sigma)
[lhs,rhs]=argn(0);
//xfinal = xfinal';
szx = size(xfinal,1);
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

runsizefinal = size(xfinal,2);
x_sol= zeros(szx,runsizefinal);

BB=(eye(szx,szx)-sigma*jac'*V_inv*jac);
for i=1:runsizefinal
//    xrs(:,i) = xr + sd.*rerror(i,:)';
    x_sol(:,i) = BB*xfinal(:,i);
//    [x_sol(i,:), f_sol(i), status(i)] = P5(xfinal(i,:)',sds.^2,xr');
end
//x_sol = x_sol';

endfunction
