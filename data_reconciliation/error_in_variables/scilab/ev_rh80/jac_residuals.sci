function [At, varargout] = jac_flowsheet_residuals(flow_full,temp_full, flow, temp, coef )
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
// flow_full:            The total flow measurements. It is a row vector with
//                        the values measured. In the case where a stream is unmeasured, an estimate must be given, eg:
//                       [ 300 200 100] (stream 3 is unmeasured, but ans estimate is given )
// temp_full:            The flow temperature measurements, same as flow_full: It is a row vector with
//                        the temperature values measured. In the case where a stream temperature is unmeasured, an estimate must be given, eg:
//                       [ 353 279 387] (the temperature of stream 3 is unmeasured, but ans estimate is given )
// flow:                The total flow with unmeasured information. It is a row vector with
//                        either the values measured, -1 if the flow is not measured or -5 if the stream is fixed
//                        (fixed variables arrise, for exemple, if we have a flow controler in the process)
//                       It is a matrix with the form:
//                       [ 300 -5 -1] 
//                      means Stream1 = 3000 ; Stream2 is fixed and Stream3 unmeasured 
//                      the fixed value information comes from flow_full input
//
// temp:                same as flow: 
//                        The values are filled with the values measured, -1 if the temperatures is unmeasured or
//                         -5 if the temperature is fixed. The fixed value information comes from temp_full input//                        
//                        Ex: for a simple splitter (one input and 2 output streams) (stream 2 with unmeasured temperature and
//                        stream's temperature 3, fixed:
//                        [ 353 -1 -5]
// coef:                the enthalpy coeficients, it is a nflow x ncoeficients (in the paper case, 3 coeficient)
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
//
//******************************************************************************

xi_full = [flow_full(:); temp_full(:)];
// x for variables classification:

xi_class = [flow(:); temp(:)];

nflow = length(flow_full);

A = diffcode_jacobian(list(flowsheet_residuals,nflow,coef),xi_full)
//pause
At=A';

varargout = list(find(xi_class == -1), find(xi_class == -5));
   
endfunction
