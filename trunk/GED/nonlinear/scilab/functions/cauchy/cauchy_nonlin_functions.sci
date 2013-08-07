// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to sum of absolute errors
// it is necessary to install the "diffcode" package using ATOMS in Scilab
// Cauchy Robust function, according to Ozyurt and Pike - Comp.  & Chem. Eng.
// 28, p. 381-402, (2004) 
function f = objfun ( x )

    e1 = (xm(red)-x(red))./(var(red).^(0.5)); 
    f = sum( (const_cauchy^2)*(log( ones(length(red),1) + (e1.^2/const_cauchy^2))));

endfunction

// gradient of the objetive function
function gf = gradf ( x )
//    gf = diffcode_jacobian(objfun,x)';
// analytical
    gf = zeros (nv,1)
    gf(red,1)  = -2*(xm(red) - x(red))./(var(red).*(1 + ((xm(red)-x(red)).^2)./(var(red).*const_cauchy^2)));

endfunction

function H = hessf ( x )
//    H = diffcode_hessian(objfun,x);
  t1 = zeros (nv,1);
	t1(red,1) = 2 ./(var(red).*(1 +  ((xm(red)-x(red)).^2)./((const_cauchy.^2).*var(red)) ))  -  (4*(xm(red)-x(red)).^2)./((const_cauchy.^2).*(var(red).^2).*(1 + ((xm(red)-x(red)).^2)./((const_cauchy.^2).*var(red))).^2);
    H=diag(t1);

endfunction
