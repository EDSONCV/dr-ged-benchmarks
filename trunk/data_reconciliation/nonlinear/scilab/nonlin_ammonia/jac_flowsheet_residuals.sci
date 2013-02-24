function [At, varargout] = jac_flowsheet_residuals(x, xfull, K_coef, cp1_coef, h_hx_coef, frac )
//******************************************************************************
// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//******************************************************************************
// Builds a Jacobian matrix of a heat exchanger network  problem proposed by Swartz, 1989,
// considering the measured and unmeasured flows and temperatures entered by user.
// This function is prepared to use the automatic derivatives toolbox of 
// Scilab. This toolbox can be instaled using the ATOMS installer (package name: diffcode).
// The objective of this jacobian is to use it further in variable classification routines
// Do not use it for optimization purposes, since it make some costly computational operations

// Outputs:
// At:    The full Jacobian of the system (as if all variables were measured).
// Ax:   The Jacobian of the system of the measured part.
// Au:   Jacobian of the system (as if all variables were measured).
// varargout(1):   unmeasured streams (identified by -1).
// varargout(2):   fixed streams (identified by -5).
// Inputs
//
// x_full:            The  measurements. It is a row vector with
//                        the values measured. In the case where a variable is unmeasured, an estimate must be given, eg:
//                       [ 300 200 100] (variable 3 is unmeasured, but ans estimate is given )
// flow:                The variable with unmeasured information. It is a row vector with
//                        either the values measured, -1 if the variable is not measured or -5 if the variable is fixed
//                        (fixed variables arrise, for exemple, if we have a flow controler in the process)
//                       It is a matrix with the form:
//                       [ 300 -5 -1] 
//                      means variable 1 = 3000 ; variable 2 is fixed and variable 3 unmeasured 
//                      the fixed value information comes from x_full input
//
// K_coef               coeficient for calculation of partition coeficient K in the flash
//                      Ki = yi/xi. Since it is not constant, it was adjusted from simulation data
//                      in a range for temperature from -24:-44 oC  
//cp1_coef              coeficient for calculation of Cp in Heat exchanger 1 considering the mixture an ideal gas
//h_hx_coef            coeficient for calculation of h in the output of the reactor and in the output of second heat
//                     exchanger. It was also adjusted from simulation data
//frac                 split fraction in the splitter
//******************************************************************************



A = diffcode_jacobian(list(flowsheet_residuals,K_coef, cp1_coef, h_hx_coef, frac),xfull)
//pause
At=A';

varargout = list(find(x == -1), find(x == -5));
   
endfunction
