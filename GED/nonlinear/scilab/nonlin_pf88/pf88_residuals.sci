function y = pf88_residuals(x)
//*********************************************************************        
// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//*********************************************************************
// this function is prepared to use the automatic derivatives toolbox of 
// Scilab. This toolbox can be instaled using the ATOMS installer (package name: diffcode).
// This function evaluates the residuals of the flowsheet of the problem proposed by Pai and Fisher 1988:
// C.C. David Pai, Gary D. Fisher
// Application of Broyden's Method to Reconciliation of Nonlinearly Constrained Data
// AICHE Journal 1988, V 34 No. 5 -p 873-876
//
// Outputs:
// y:        the contraints residuals
// Inputs:
// x:         the column vector of the variables: 
// 8 variables and 6 equations
//
xi = x(1:5);
vi = x(6:$);
y(1) = 0.5*xi(1).^2 - 0.7*xi(2) + xi(3)*vi(1) + vi(1)*vi(2)*xi(2).^2  + 2*xi(3)*vi(3).^2;
y(2) = xi(1) - 2*xi(2) + 3*xi(1)*xi(3) - 2*xi(2)*vi(1) - xi(2)*vi(2)*vi(3);
y(3) = xi(3)*vi(1) - xi(1) + 3*xi(2) + xi(1)*vi(2) - xi(3)*vi(3).^0.5;
y(4) = xi(4) - xi(1) - xi(3).^2 + vi(2) + 3*vi(3);
y(5) = xi(5) - 2*xi(3)*vi(2)*vi(3);
y(6) = 2*xi(1) + xi(2)*xi(3)*vi(1) + vi(2) - vi(3);

endfunction
