// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to weighted least squares functions

function f = objfun ( x )

	f = sum(((xm(red)-x(red)).^2)./var(red));

endfunction

////////////////////////////////////////////////////////////////////////
// Define gradient and Hessian matrix

function gf = gradf ( x )
gf = zeros (nv,1)
gf(red,1) =2*(x(red) - xm(red))./var(red);

endfunction

function H = hessf ( x )
    t1 = zeros (nv,1);
	t1(red,1) = (2*ones(length(red),1)./var(red));
    H=diag(t1);
endfunction

