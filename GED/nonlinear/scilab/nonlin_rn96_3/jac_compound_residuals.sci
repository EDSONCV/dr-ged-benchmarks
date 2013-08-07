function [At, varargout] = jac_compound_residuals(Ainput,Ncomp,totalFlowMeasured,compoundMeasured, totalFlowUnmeasured, compoundUnmeasured)
//******************************************************************************
// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//******************************************************************************
// Builds a Jacobian matrix of a multicompound balance considering the measured and
// unmeasured flows and compositions entered by the user
// Normalization equations: \sum_{xi}^n x_{i,j} -1 = 0  are appended at the end 
// of A
// Outputs:
// At:    The full Jacobian of the system (as if all variables were measured).
// varargout(1):   unmeasured streams (identified by -1).
// varargout(2):   fixed streams (identified by -5).
// Inputs
//
// Ainput:               The incidence matrix considering only the total flow
//                       Example: For a simple splitter:   Ainput = [1 -1 -1]
// Ncomp:                Number of compounds 
// totalFlowMeasured:    The total flow measurements. It is a row vector with
//                        either the values measured. In the case where a stream is unmeasured, an estimate must be given, eg:
//                       [ 300 200 100] (stream 3 is unmeasured, but ans estimate is given )
// compoundMeasured:     mass or molar fraction of the compounds in the stream, it is a matrix with 
//                        compounds in lines and streams in columns.
//                        x_{j,i}, where i= stream and j = compounds 
//                        The values are filled either the values measured or -1 if the composition is not measured.
//                        [ x11  x12   x13
//                          x21  x22  x23  
//                         x31  x32   x33]
//                        In the case where a stream is unmeasured, an estimate must be given.
//                        Ex: for a simple splitter (one input and 2 output streams) with 3 compounds:
//                        [ 0.5 0.55   0.75
//                           0.4  0.15  .05 
//                          0.1   0.3  0.2]
// totalFlowUnmeasured:    The total flow with unmeasured information. It is a row vector with
//                        either the values measured, -1 if the flow is not measured or -5 if the stream is fixed
//                        ()fixed variables arrise, for exemple, if we have a flow controler in the process)
//                       It is a matrix with the form:
//                       [ 300 -5 -1] 
//                      means Stream1 = 3000 ; Stream2 is fixed and Stream3 unmeasured 
//                      the fixed value information comes from totalFlowMeasured input
//
// compoundUnmeasured:     mass or molar fraction of the compounds in the stream, it is a matrix with 
//                        compounds in lines and streams in columns.
//                        x_{j,i}, where i= stream and j = compounds 
//                        The values are filled with the values measured, -1 if the composition is not measured or
//                         -5 if the composition is fixed. The fixed value information comes from compoundMeasured input
//                        [ x11  x12   x13
//                          x21  x22  x23  
//                         x31  x32   x33]
//                        Ex: for a simple splitter (one input and 2 output streams) with 3 compounds (compound 2 unmeasured
//                        in all streams) and compound 3 fixed in stream 3:
//                        [ 0.5 0.5 0.5;
//                           -1  -1  -1 ;
//                          0.1 .1  5]
// Example for a 3 stream system with 3 compounds
// for a simple splitter where the incidence matrix is Ainput = [a11 a12 a13]
// the resulting system is:
//      eq1 = a11.F1.x11 + a12.F2.x12 +a13.F3.x13
//      eq2 = a12.F1.x21 + a12.F2.x22 +a13.F3.x23
//      eq3 = a12.F1.x31 + a12.F2.x32 +a13.F3.x33
// next we'll build the derivatives of the normalization equations:
// for each stream, we have \sum_{xi}^n x_{i,j} -1 = 0, , resulting in 3 equaitons  
//      eq4 = x11 + x21 + x31 = 1
//      eq5 = x12 + x22 + x32 = 1
//      eq6 = x13 + x23 + x33 = 1

//  the resulting matrix "A" is :
//
//           Ftot1     Ftot2        Ftot3     x11   x12    x13     x21      x22     x23     x31    x32    x33
//
//  eq1 |    a11.x11     a12.x12     a13.x13  a11.F1 a12.F2 a13.F3  0       0          0      0      0      0  
//  eq2 |    a11.x21     a12.x22     a13.x23    0    0       0      a11.F1 a12.F2   a13.F3    0       0      0
//  eq3 |    a11.x31     a12.x32     a13.x33   0     0       0      0       0          0    a11.F1 a12.F2 a13.F3
//  eq4 |     0             0          0       1     0      0       1       0          0     1        0      0
//  eq5 |     0             0          0       0     1      0       0       1          0     0        1      0
//  eq6 |     0             0          0       0     0      1       0       0          1     0        0      1
// Then Ax and Au is built appropriatelly based in the totalFlowMeasured or compoundMeasured
// Notice that even if a measurement is zero (either totalflow or composition) do not set it as zero, because it may lead to 
// a misleading classification. It is suggested to set a very small value , say 1.0e-2.
// Another point is the linear dependent columns or rown, because it may lead to a missleading classification:
// Avoid compounds structures like the one below:
// [ 0.1 0.1 0.1
//   0.1 0.1 0.1
//   0.1 0.1 0.1]
//
//******************************************************************************
[Aeqp, Astreams] =size(Ainput);
// organize the order of the elements of compoundMeasured because the (:) operator mess things up

x=[totalFlowMeasured(:);matrix(compoundMeasured',-1)];


A = diffcode_jacobian(list(compound_residuals,Ainput,Ncomp),x);
At=A';
 
compoundUnmeasured_t = compoundUnmeasured';
//cc = [totalFlowUnmeasured(:); matrix(compoundUnmeasured',-1)];
cc = [totalFlowUnmeasured(:); compoundUnmeasured_t(:)];

varargout = list(find(cc == -1), find(cc == -5));
 
    
endfunction
