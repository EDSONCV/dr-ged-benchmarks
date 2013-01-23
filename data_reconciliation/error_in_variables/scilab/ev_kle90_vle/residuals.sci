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
// y:        the contraints residuals. The number of constraints depends on the number of the 
//            data sets choosen 
// Inputs:
// x:         the column vector of the variables: x - it is reorganized inside this function


xdata = matrix(x(1:$-2),ndata,4)';
//pause
P = xdata(1,:)';
y1 = xdata(2,:)';
T = xdata(3,:)';
x1 = xdata(4,:)';
//Aa = x($-1);
//Bb = x($);
//gamma1 = exp((Aa./(R_gas*T)).*(1 + (Aa./Bb).*(x1./(1-x1))).^(-2) );
//gamma2 = exp((Bb./(R_gas*T)).*(1 + (Bb./Aa).*((1-x1)./x1)).^(-2) );
//p1o = exp(18.5875 - (3626.55./(T -34.29)) );
//p2o = exp(16.1764 - (2927.17./(T - 50.22)) );
//
teta1 = x($-1);
teta2 = x($);
gamma1 = exp((teta1./T).*(1 + (teta1/teta2).*(x1./(1-x1))).^(-2) );
gamma2 = exp((teta2./T).*(1 + (teta2/teta1).*((1-x1)./x1)).^(-2) );
p1o = exp(18.5875 - (3626.55./((T.*T_ref) -34.29)) );
p2o = exp(16.1764 - (2927.17./((T.*T_ref) - 50.22)) );


y=zeros(2*ndata,1);

y(1:ndata,1) = gamma1.*x1.*p1o -y1.*P;
y(ndata+1:2*ndata,1) = gamma2.*(1 - x1).*p2o  - (1-y1).*P;

endfunction
