function y = residuals(x)
//*********************************************************************        
// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//*********************************************************************
// this function is prepared to use the automatic derivatives toolbox of 
// Scilab. This toolbox can be instaled using the ATOMS installer (package name: diffcode).
// This function evaluates the residuals of the flowsheet of the problem proposed by Swartz, 1989
//
// Outputs:
// y,:        the contraints residuals
// Inputs:
// x:         the column vector of the variables: x = [flow, temperatures]

// all the equations are coded "by hand":
e1 = x(1:28);
e2 = x(29:56);
e3 = x(57:84);
//e3 = xmfull(57:84);
teta1 = x(85);
teta2 = x(86);
teta3 = x(87);
y = e3 - ((teta1.*(teta2.^2).*teta3.*e1.*(e2.^2))./((ones(1,28)'+teta1.*e1+teta2.*e2).^3));

endfunction
