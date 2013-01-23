function y = residuals_cstr(x)
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
// y:        the contraints residuals. The number of constraints depends on the number of the 
//            data sets choosen 
// Inputs:
// x:         the column vector of the variables: x - it is reorganized inside this function


xdata = matrix(x(1:$-2),ndata,5)';
//pause
CAo = xdata(1,:)';
CA  = xdata(2,:)';
CB  = xdata(3,:)';
To  = xdata(4,:)';
T   = xdata(5,:)';

//
teta1 = x($-1);
teta2 = x($);
k = teta1*exp(-teta2.*((T_ref*ones(ndata,1)./T) - 1));

y=zeros(3*ndata,1);

y(1:ndata,1) = itau.*CAo - itau.*CA - k.*CA;
y(ndata+1:2*ndata,1) = -itau.*CB + k.*CA;
y(2*ndata+1:3*ndata,1) = itau.*To - itau.*T - (deltaH_r/(rho*Cp)).*k.*CA;



endfunction
