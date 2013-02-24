function y = flowsheet_residuals(x, Nflow, coef)
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
// Nflow:     number of flow streams in the x vector
//            the enthalpy will be evaluated in the following way:
// coef:                the enthalpy coeficients: it is a nflow x ncoeficients (in the paper case, 3 coeficients)
//                      they will be calculated according to the following formula:
//                      h_{i} = nu_{1,i} + nu_{2,i}*temp_{i}  + nu_{3,i}*temp_{i}^{2}  
//                      in the matrix form (simplified example for 4 streams):
//                      temp_exp= [ 1          1            1        1;
//                                 t_{1}       t_{2}      t_{3}     t_{4};
//                                 t_{1}^{2}] t_{2i}^{2}] t_{3}^{2} t_{4}^{2} ]
//
//                      coef = [ nu_{1,1} nu_{1,2} nu_{1,3} nu_{1,4} ;
//                               nu_{2,1} nu_{2,2} nu_{2,3} nu_{2,4} ;
//                               nu_{3,1} nu_{3,2} nu_{3,3} nu_{3,4} ]
//
//                      enthalpy_internal = sum(temp_exp.*coef, 'r'); 
// The resulting system is according to the papper from Swartz (1989):
//(1) = fA(1) - fA(2);
//y(2) = fB(2) - fB(3);
//y(3) = fA(2) - fA(3) - fA(6);
//y(4) = fA(3) - fA(4);
//y(5) = fB(1) - fB(2);
//y(6) = fA(4) - fA(5);
//y(7) = fC(1) - fC(2);
//y(8) = fA(5) + fA(7) - fA(8);
//y(9) = fA(6) - fA(7);
//y(10) = fD(1) - fD(2);
//
//// now the enthalpy balances
//
//y(11) = fA(1)*hA(1) - fA(2)*hA(2) + fB(2)*hB(2) - fB(3)*hB(3);
//y(12) = fA(3)*hA(3) - fA(4)*hA(4) + fB(1)*hB(1) - fB(2)*hB(2);
//y(13) = fA(4)*hA(4) - fA(5)*hA(5) + fC(1)*hC(1) - fC(2)*hC(2);
//y(14) = fA(5)*hA(5) + fA(7)*hA(7) - fA(8)*hA(8);
//y(15) = fA(6)*hA(6) - fA(7)*hA(7) + fD(1)*hD(1) - fD(2)*hD(2);
//// temperature equality in mixers
//y(16) = tempA(2)-tempA(3) ;
//y(17) = tempA(2)-tempA(6) ;
// all the equations are coded "by hand":

flow_internal = x(1:Nflow);
fA = flow_internal(1:8);
fB = flow_internal(9:11);
fC = flow_internal(12:13);
fD = flow_internal(14:15);

temp_internal = x(Nflow + 1: $)';

tempA = temp_internal(1:8);
tempB = temp_internal(9:11);
tempB = temp_internal(12:13);
tempD = temp_internal(14:15);

// relationships between h and t

temp_exp = [ones(1,Nflow)
            temp_internal 
            temp_internal.^2];

enthalpy_internal = sum(temp_exp.*coef, 'r');    


hA = enthalpy_internal(1:8);
hB = enthalpy_internal(9:11);
hC = enthalpy_internal(12:13);
hD = enthalpy_internal(14:15);


// first we build the mass balance

y(1) = fA(1) - fA(2);
y(2) = fB(2) - fB(3);
y(3) = fA(2) - fA(3) - fA(6);
y(4) = fA(3) - fA(4);
y(5) = fB(1) - fB(2);
y(6) = fA(4) - fA(5);
y(7) = fC(1) - fC(2);
y(8) = fA(5) + fA(7) - fA(8);
y(9) = fA(6) - fA(7);
y(10) = fD(1) - fD(2);

// now the enthalpy balances

y(11) = fA(1)*hA(1) - fA(2)*hA(2) + fB(2)*hB(2) - fB(3)*hB(3);
y(12) = fA(3)*hA(3) - fA(4)*hA(4) + fB(1)*hB(1) - fB(2)*hB(2);
y(13) = fA(4)*hA(4) - fA(5)*hA(5) + fC(1)*hC(1) - fC(2)*hC(2);
y(14) = fA(5)*hA(5) + fA(7)*hA(7) - fA(8)*hA(8);
y(15) = fA(6)*hA(6) - fA(7)*hA(7) + fD(1)*hD(1) - fD(2)*hD(2);
// temperature equality in mixers
y(16) = tempA(2)-tempA(3) ;
y(17) = tempA(2)-tempA(6) ;

endfunction
