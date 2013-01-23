// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to weighted least squares functions

function f = objfun ( x )

f = obj_user(x);

endfunction

function c = confun(x)
//flowsheet_residuals is it's own file
	c =  residuals(x);

endfunction

////////////////////////////////////////////////////////////////////////
// Define gradient and Hessian matrix
// gradient of the objetive function
function gf = gradf ( x )

// it is encouraged to use the code bellow to become more generic 
//(eg. in case the user want's to use a objective function form robust
// statistics ).
//gf = diffcode_jacobian(obj_user,x)';
// If one wants to use the classical weighted 
// least squares approach, the evaluation of the hessian using the function
// bellow is less lime consuming    
gf=[2*(x(1:84) - xmfull(1:84))./(sd(1:84).^2) ;0;0;0];

endfunction

// hessian of the objetive function

function H = hessf ( x )
// it is encouraged to use the code bellow to become more generic 
//(eg. in case the user want's to use a objective function form robust
// statistics ).
//    H = diffcode_hessian(obj_user,x);
// If one wants to use the classical weighted 
// least squares approach, the evaluation of the hessian using the function
// bellow is less time consuming    
	H = diag([2*ones(84,1)./(sd(1:84).^2);0 ;0;0]);

endfunction

// gradient of the constraints

function y = dg1(x)
//flowsheet_residuals is it's own file
ytmp = diffcode_jacobian(residuals,x)';
//disp('inside dg')
//pause
for i = 1: nnzjac; 
  y(i)=ytmp(sparse_dg(i,1),sparse_dg(i,2)); 
end
//pause
endfunction


// The Hessian of the Lagrangian

function y = dh(x,lambda,obj_weight)
//  disp('inside dh')
    ysum = zeros(nv,nv);
    if obj_weight <> 0 then
         yobj = obj_weight * hessf ( x );
     else
         yobj = zeros(nv,nv);
    end
//    pause
    if sum(abs(lambda)) <> 0           then
        // the hessian of the constraints
        //flowsheet_residuals is it's own file
//        ytmpconstr = diffcode_hessian(residuals,x);
        [J,ytmpconstr] = derivative(residuals, x);

        for i = 1: nc; 
            if lambda(i) <> 0 then

//                ysum = ysum + lambda(i)*ytmpconstr(:,:,i); 
                ytmp(i,:) = lambda(i)*ytmpconstr(i,:); 

            end
        end
        ysum = matrix(sum(ytmp,'r'),nv,nv);
        
    else
         ysum = zeros(nv,nv);

    end

    ysumall = ysum + yobj;
//   pause    
    for i = 1: nnz_hess
        y(i) = ysumall(sparse_dh(i,1),sparse_dh(i,2));
    end

endfunction

function y = obj_user( x )

y = sum(((xmfull(1:84)-x(1:84))./sd(1:84)).^2 );

endfunction    


function [nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = wls_structure_reactor(x)
// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to ipopt solver
//***************************************************************
//This function analyses the structure of the problem 
//and return vectors and matrices that will be used by ipopt solver 
//Outputs:
//
//    nc:                 number of constraints
//    nv:                 number of variables
//    nnzjac              number of non zero elements in the Jacobian  of the constraints
//    nnz_hess            number of non zero elements in the Lagrangean Hessian's
//    sparse_dg           sparsity structure of the Jacobian matrix of the constraints
//    sparse_dh           sparsity structure of the Lagrangean Hessian's
//    lower               lower bound of the variables
//    upper               upper bound of the variables
//    var_lin_type        type of the variable (linear or non-linear)
//    constr_lin_type     type of the constraints (linear or non-linear)
//    constr_lhs         lower bound of the constraints residuals
//    constr_rhs         upper bound of the constraints residuals
//
//Inputs:
//    x                    The variables
// From here on, the problem generation is automatic
// No need to edit below
// Arrange the vectors

xlocal=x;

// Jacobian and its structure
jactst = diffcode_jacobian(residuals,xlocal)';

jactstSparse=sparse(jactst);
[ij,v,mn]=spget(jactstSparse);
// index of the non-zero elements of the Jacobian
nnzjac = size(ij,1);
// The sparsity structure of the constraints
sparse_dg = ij;

//The problem size: nc = number of constraints and an number of variables
[nc,nv] = size(jactst);


// The sparsity structure of the Hessian Lagrangian
// the Hessian of the objective function is diagonal but the hessian of the constraints not!

//first retrieve the constraints Hessian structure, notice that with diffcode_hessian, the
// Hessian has the following dimensions: nvar x nvar x nconstr , so it is in fact a 
// 3 dimensional matrix.

hess_constr_tst = diffcode_hessian(residuals,xlocal);

// cumulative sums the constraints in one hessian
hess_constr = zeros(nv,nv);
//pause
for i = 1: nc
    hess_constr = hess_constr + abs(hess_constr_tst(:,:,i));
end
// the Hessian of the objective function

hess_f = diffcode_hessian(objfun,xlocal)
//pause
// sum both of the Hessians
hess_Sparse=sparse(hess_constr + hess_f);
// get the hessian structure
[ij_hess,v_hess,mn_hess]=spget(hess_Sparse);

//filters the hessian to remove symmetric indexes
ij_hess_filtered = filter_symmetric(ij_hess);


// index of the non-zero elements of the Hessian
sparse_dh = ij_hess_filtered;
nnz_hess = length(ij_hess_filtered(:,1));

// the variables have lower bounds of 0
lower = zeros(nv,1);
lower (85) = 5;
lower (86) = 0.4;
lower (87) = 1000;
lower ($) = 0;
// the variables have upper bounds of 50000
upper = 50000*ones(nv,1);
upper (85) = 10;
upper (86) = 1;
upper (87) = 2000;
upper ($) = 20000;

// in the non-linear case, all constraints and variables are non-linear (bilinear)
var_lin_type(1:nv) = 1; // Non-Linear
constr_lin_type (1:nc) = 1; // Non-Linear


// the constraints has lower and upper bound of 0
constr_lhs(1:nc) = 0;
constr_rhs(1:nc) = 0;

endfunction

