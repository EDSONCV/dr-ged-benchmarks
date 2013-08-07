// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to sum of absolute errors
// it is necessary to install the "diffcode" package using ATOMS in Scilab
// smooth functions according to Gopal and Biegler
// AICHE Journal 45(7) 1535-1547 - July 1999
// Quasi Weighted Robust function, according to Zhang et al. - Comp.  & Chem. Eng.
// 34, p. 154-162-402, (2010) 
function f = objfun ( x )

    e1 = (xm(red)-x(red))./(var(red).^(0.5)); 
// smoothing functions
// for sigmoidal function (Eq. 24 from paper)
//    abs_error = sum(sig1=1./alpha_smooth*log(2+exp(alpha_smooth*e1)+exp(-alpha_smooth*e1)));
    // for interior point function (Eq 25 from paper)
    abs_error = (e1.^2 + beta_smooth.^2).^0.5;
	// sigmoidal, but based in max operator property (Eq 28 from paper)
    // this one leads to a small error when e1 = 0  
    
//    abs_error = e1 + beta_smooth*log(1+exp(-2*alpha_smooth*e1));
    f = sum( ((e1.^(2))./(2 + const_qw*abs_error))  );

endfunction

// gradient of the objetive function
function gf = gradf ( x )
// in the future we can express this function analytically
//    gf = diffcode_jacobian(objfun,x)';
    gf = zeros(nv,1);
    sqrarg = const_qw.*sqrt(((xm(red)-x(red)).^2)./(var(red)) + beta_smooth.^2);
    sqrargdiv = sqrarg + 2;
    gf(red,1) = (sqrarg.*(xm(red)-x(red)))./(var(red).*sqrargdiv.^2)  - (2*(xm(red)-x(red)))./(var(red).*sqrargdiv);
endfunction

function H = hessf ( x )
// For the robust functions, the lagrangean of the objective function is not constant
// as in weigthed least squares.

// in the future we can express this function analytically
//    H = diffcode_hessian(objfun,x);

    onesqw = ones(length(red),1);
    sqrarg = const_qw.*sqrt(((xm(red)-x(red)).^2)./(var(red)) + beta_smooth.^2);
    sqrarg2 =  sqrarg./const_qw;
    sqrargdiv = sqrarg + 2;
    t1 = zeros (nv,1);
	t1(red,1) = -3.*const_qw.*((xm(red) - x(red)).^2)./((var(red).^2).*(sqrargdiv.^2).*sqrarg2) + 2*(const_qw.^2).*((xm(red) - x(red)).^2)./(var(red).*sqrargdiv.^3) - (sqrarg)./(var(red).*sqrargdiv.^2) + (2*onesqw)./(var(red).*sqrargdiv);
    H=diag(t1);
endfunction
