// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// adjustability and detectability
function [adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac)
    sigma_inv=inv(sigma);
// variance-covariance matrix: narasimham pg. 178 eq. 7-3
    V=jac*sigma*jac';
    V_inv= inv(V);
    diag_diag_V = diag(diag(V));
// covariance matrix of adjustments: narasimham pg. 183 eq. 7-13
    Wbar=sigma*jac'*V_inv*jac*sigma;
    varx=diag(sigma);
    rowjac=size(jac,1);
    coljac=size(jac,2);
    W=inv(sigma);
    M=[W, jac';jac,zeros(rowjac,rowjac)];
    MI= inv(M);
    // variance of the reconciled variables calculated according to:
    //Heyen, G, E. Marechál, and B. Kalitventzeff. 1996. Sensitivity calculations and variance analysis in plant measurement reconciliation. Computers & Chemical Engineering 20, no. 96: S539-S544. doi:10.1016/0098-1354(96)00099-3. http://linkinghub.elsevier.com/retrieve/pii/0098135496000993.
    for i=1:coljac
        a(i)=sum((MI(i,1:coljac).^2)./varx')
//Narasimhan pg 211 eq. 7.73
        adj(i)= 1-sqrt(a(i))/sqrt(varx(i))
//Narasimhan pg 211 eq. 7.73
        detect(i)=sqrt(1-a(i)/varx(i))
    end

endfunction