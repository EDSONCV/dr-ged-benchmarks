// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to weighted least squares functions

function f = objfun ( x )

	f = sum(((xm(measured)-x(measured)).^2)./var(measured));

endfunction

function c = confun(x)

	c = jac*x;

endfunction

////////////////////////////////////////////////////////////////////////
// Define gradient and Hessian matrix

function gf = gradf ( x )
gf = zeros(nv,1);
tmp2=2*(x(measured) - xm(measured))./var(measured);
gf(measured)=tmp2;

endfunction

function H = hessf ( x )
    tmp1 = zeros(nv,1);
    tmp1(measured) = 1;
	H = 2*tmp1./var;
endfunction

function y = dg1(x)

for i = 1: nnzeros; 
  y(i)=jac(i1(i),i2(i)); 
end

endfunction

// The Lagrangian
// In fact, IPOPT does not call this line because we say that the constraints are constants
function y = dh(x,lambda,obj_weight)
	y = obj_weight * hessf ( x );
endfunction

// The constraints
function y=dg(x)

	y = dg1(x)
	
endfunction

function [nc, nv, i1, i2, nnzeros, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = wls_structurel(jac)
    // From here on, the problem generation is automatic
// No need to edit below
//The problem size: nc = number of constraints and an number of variables
[nc,nv] = size(jac);
// index of the non-zero elements of the Jacobian
[i1,i2]=find(jac<>0);

nnzeros = size(i1,2);

// The sparsity structure of the constraints

sparse_dg = [i1', i2'];

// The sparsity structure of the Lagrangian
// the Hessian for this problem is diagonal
sparse_dh = [ [1:nv]', [1:nv]'];

// the variables have lower bounds of 0
lower = (-8*var.^(0.5)).*ones(nv,1);
// the variables have upper bounds of 50000
upper = 50000*ones(nv,1);
var_lin_type(1:nv) = 1; // Non-Linear
constr_lin_type (1:nc) = 0; // Linear

// the constraints has lower bound of 0
constr_lhs(1:nc) = 0;
// the constraints has upper bound of 0.
constr_rhs(1:nc) = 0;
endfunction

