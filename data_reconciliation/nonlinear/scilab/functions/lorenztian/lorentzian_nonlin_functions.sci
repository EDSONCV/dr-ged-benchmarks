// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to sum of absolute errors
// it is necessary to install the "diffcode" package using ATOMS in Scilab
// Lorentzian Robust function, according to Ozyurt and Pike - Comp.  & Chem. Eng.
// 28, p. 381-402, (2004) 
function f = objfun ( x )

    e1 = (xm(red)-x(red))./(var(red).^(0.5));
    f = sum( -ones(length(red),1)./(ones(length(red),1) + (e1.^2/(2*const_lor^2)))  );

endfunction

// gradient of the objetive function

function gf = gradf ( x )
//    gf = diffcode_jacobian(objfun,x)';
// analytical
gf = zeros (nv,1)
gf(red,1) = -1*(xm(red)-x(red))./(((1 + ((xm(red)-x(red)).^(2))./(2*var(red)*const_lor^2))^2).*(var(red)*const_lor^2));

endfunction

function H = hessf ( x )
// For the robust functions, the lagrangean of the objective function is not constant
// as in weigthed least squares.
//    H = diffcode_hessian(objfun,x);
// analytical
// TODO CHECK THIS HESSIAN WITH UNMEASURE VARIABLES
    t1 = zeros (nv,1);
	t1(red,1) = -2*((xm(red)-x(red))^2)./(((1 + ((xm(red)-x(red)).^(2))./(2*var(red)*const_lor^2)).^3).*((var(red).^2)*const_lor.^4)) + (((1 + ((xm(red)-x(red)).^(2))./(2*var(red)*const_lor^2))^2).*(var(red)*const_lor^2))^(-1);
    H = diag(t1);

endfunction
