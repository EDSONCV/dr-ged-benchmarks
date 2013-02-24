// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to sum of absolute errors
// it is necessary to install the "diffcode" package using ATOMS in Scilab
// Smooth functions according to Gopal and Biegler
// AICHE Journal 45(7) 1535-1547 - July 1999
// Hampel function, according to Ozyurt and Pike - Comp.  & Chem. Eng.
// 28, p. 381-402, (2004) 
function f = objfun ( x )

    e1 = (xm(red)-x(red))./(var(red).^(0.5));
    
// smoothing functions
// for sigmoidal function (Eq. 24 from paper)
//    abs_error = sig1=1./alpha_smooth*log(2+exp(alpha_smooth*e1)+exp(-alpha_smooth*e1));
// for interior point function (Eq 25 from paper)
    abs_error = (e1.^2 + beta_smooth.^2).^0.5;
	// sigmoidal, but based in max operator property (Eq 28 from paper)
    // this one leads to a small error when e1 = 0  
    
//    abs_error = e1 + beta_smooth*log(1+exp(-2*alpha_smooth*e1));
    // the logic of the Hampel robust estimator is long.
    // TODO IMPROVE THIS LOGIC
    
//    from -a to a    
    th1 = (tanh(100*(e1 + ones_a)) + ones_xm)./2;
    th2 = -(tanh(100*(e1 - ones_a)) - ones_xm)./2;
//
    hap1 = 0.5*e1.^2; 
//
    hap11 = hap1.*(th1 + th2 - ones_xm) ;
//    // from a < abs(e1) < b
//    
    th3 = (tanh(100*(e1 - ones_a)) + ones_xm)/2;
    th4 = -(tanh(100*(e1 - ones_b)) - ones_xm)/2;
    th3n = (tanh(100*(-e1 - ones_a)) + ones_xm)/2;
    th4n = -(tanh(100*(-e1 - ones_b)) - ones_xm)/2;
//    
    hap2a=a*abs_error - 0.5*ones_a.^2; 
    hap22 = hap2a.*(th3+th4-1) + hap2a.*(th3n+th4n - ones_xm);
        // from b < abs(e1) < c

    th5 = (tanh(100*(e1  - ones_b)) + ones_xm)/2;
    th6 = -(tanh(100*(e1 - ones_c)) - ones_xm)/2;

    th5n = (tanh(100*(-e1  - ones_b)) + ones_xm)/2;
    th6n = -(tanh(100*(-e1 - ones_c)) - ones_xm)/2;
//pause

    hap3 = ones_a.*ones_b - 0.5*ones_a + 0.5*(c-b)*(1-((c*ones_xm-abs(e1))/(c-b)).^2);
    hap33 = hap3.*(th5 + th6 - ones_xm) + hap3.*(th5n + th6n - ones_xm);

    // from c < abs(e1) 
    
    th7 = (tanh(100*(e1 - ones_c)) + ones_xm)/2;
    th7n = (tanh(100*(-e1 - ones_c)) + ones_xm)/2;

    hap4 = ones_a.*ones_b - 0.5*ones_a + 0.5*(ones_c-ones_b); // > c
    hap44 = hap4.*th7 + hap4.*th7n;

    f = sum(hap11 + hap22+ hap33 + hap44);
//    f = (hap11 + hap22+ hap33 + hap44);
//    f = sum(hap11);

endfunction

// gradient of the objetive function
function gf = gradf ( x )
// in the future we can express this function analytically
// For Hampel, we are using the finite difference formula due to a limitation of
// diffcode when providing the exact differences of tanh
//    gf = diffcode_jacobian(objfun,x)';
//    gf = derivative(objfun, x, 1.0e-2, order = 4)';
    gf = derivative(objfun, x, order = 4)';
    

endfunction

function H = hessf ( x )
// For the robust functions, the lagrangean of the objective function is not constant
// as in weigthed least squares.

// in the future we can express this function analytically
// For Hampel, we are using the finite difference formula due to a limitation of
// diffcode when providing the exact differences of tanh
//  H = diffcode_hessian(objfun,x);
//[J,H] = derivative(objfun, x, H_form = "hypermat");
[J,H] = derivative(objfun, x, H_form = "hypermat");
//[J,H] = derivative(objfun, x, 1.0e-2, order =2, H_form = "hypermat");

endfunction
