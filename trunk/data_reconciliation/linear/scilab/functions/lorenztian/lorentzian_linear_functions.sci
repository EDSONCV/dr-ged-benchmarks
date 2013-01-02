// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to sum of absolute errors
// it is necessary to install the "diffcode" package using ATOMS in Scilab
// Lorentzian Robust function, according to Ozyurt and Pike - Comp.  & Chem. Eng.
// 28, p. 381-402, (2004) 
function f = objfun ( x )

    e1 = (xm-x)./(var.^(0.5));
    f = sum( -ones(nv,1)./(ones(nv,1) + (e1.^2/(2*const_lor^2)))  );

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

////////////////////////////////////////////////////////////////////////
// Define constraints, gradient and Hessian matrix


// The constraints function, Jacobian and Hessian

// First the vector of inequalyties and equalyties 
// We generallu don't use inequalyties in classical DR problems
// but it is up to the user use it or not

function c = confun(x)
    

    if   nnzjac_ineq <> 0 then
        c1 =  res_ineq(x);
        c =[ c1; res_eq(x)];
    else
        c =[ res_eq(x)']
    end

endfunction

function c = res_ineq(x)

	c = [];

endfunction


function c = res_eq(x)

	c = jac*(x);

endfunction

function y=dg(x)

    if   nnzjac_ineq <> 0 then

        y = [dg_ineq(x);dg_eq(x)];

    else

        y = [dg_eq(x)];

    end

endfunction

function y = dg_ineq(x)

    if   nnzjac_ineq <> 0 then
        
        ytmp = diffcode_jacobian(res_ineq,x)';

        for i = 1: nnzjac_ineq; 
            y(i)=ytmp(sparse_dg(i,1),sparse_dg(i,2)); 
        end
        
    else 
        y =[];
    end
endfunction


function y = dg_eq(x)
// we use the jacobian here, if use wants to use a different Jacobian , comment and 
// uncomment the lines approprieatelly
//   ytmp = diffcode_jacobian(res_eq,x)';

    for i = nnzjac_ineq + 1: nnzjac_ineq + nnzjac_eq; 
//        y(i - nnzjac_ineq)=ytmp(sparse_dg(i-nnzjac_ineq,1),sparse_dg(i-nnzjac_ineq,2)); 
        y(i - nnzjac_ineq)=jac(sparse_dg(i-nnzjac_ineq,1),sparse_dg(i-nnzjac_ineq,2)); 
    end

endfunction

// The Hessian of the Lagrangian

function y = dh(x,lambda,obj_weight)

    ysum = zeros(nv,nv);
    if obj_weight <> 0 then        
        yobj = obj_weight * hessf ( x );
    else
        yobj = zeros(nv,nv);
    end
    if sum(abs(lambda)) <> 0  & n_non_lin_eq > 2  then
        // the hessian of the constraints
        ytmpconstr = diffcode_hessian(confun,x);

        for i = 1: nc; 
            if lambda(i) <> 0 then

                ysum = ysum + lambda(i)*ytmpconstr(:,:,i); 

            end
        end

    else
        ysum = zeros(nv,nv);

    end

    ysumall = ysum + yobj;


    for i = 1: nnz_hess
        y(i) = ysumall(sparse_dh(i,1),sparse_dh(i,2));
    end

endfunction


