function [At, varargout] = jac_pf88_residuals(xm_full, measured)
//******************************************************************************
// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//******************************************************************************
// Builds a Jacobian matrix of the generic non-linear problem from Pai & Fisher (88):
// C.C. David Pai, Gary D. Fisher
// Application of Broyden's Method to Reconciliation of Nonlinearly Constrained Data
// AICHE Journal 1988, V 34 No. 5 -p 873-876
//
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
// xm_full:              The measured variable information plus the unmeasured [(x1:x5) (v1:v3)]' 
//                        In the case where a stream is unmeasured, an estimate must be given, eg:
//                       [ 300 200 100] (variable 3 is unmeasured, but ans estimate is given )
// measured:             The vector of measrued/unmeasured information. It is a row vector with
//                        either the values measured, -1 if the flow is not measured or -5 if the stream is fixed
//                        (fixed variables arrise, for exemple, if we have a flow controler in the process)
//                       It is a matrix with the form:
//                       [ 300 -5 -1] 
//                      means variable 1 = 3000 ; variable2 is fixed and variable3 unmeasured 
//                      the fixed value information comes from xm_full input
//
//
//******************************************************************************

xi_full = [xm_full(:)];
// x for variables classification:

xi_class = [measured(:)];


A = diffcode_jacobian(pf88_residuals,xi_full)
//pause
At=A';

varargout = list(find(xi_class == -1), find(xi_class == -5));
   
endfunction
