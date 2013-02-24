// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to sum of absolute errors
// it is necessary to install the "diffcode" package using ATOMS in Scilab
// Logistic Robust function, according to Ozyurt and Pike - Comp.  & Chem. Eng.
// 28, p. 381-402, (2004) 
function f = objfun ( x )

    e1 = (xm(red)-x(red))./(var(red).^(0.5));
    f = sum( 2*(log( ones(length(red),1) + exp(e1/const_logist))) -e1/const_logist);

endfunction

// gradient of the objetive function
function gf = gradf ( x )
// in the future we can express this function analytically
    gf = diffcode_jacobian(objfun,x)';

endfunction

function H = hessf ( x )
// For the robust functions, the lagrangean of the objective function is not constant
// as in weigthed least squares.

// in the future we can express this function analytically
    H = diffcode_hessian(objfun,x);

endfunction
