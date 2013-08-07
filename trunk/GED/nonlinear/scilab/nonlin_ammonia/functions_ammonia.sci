// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to weighted least squares functions


function c = confun(x)
//flowsheet_residuals is it's own file
	c =  flowsheet_residuals(x,K_coef, cp_1, hh, 0.22);

endfunction

////////////////////////////////////////////////////////////////////////

// gradient of the constraints

function y = dg1(x)
//flowsheet_residuals is it's own file
ytmp =  derivative(list(flowsheet_residuals,K_coef, cp_1, hh, 0.22), x, order = 2);
//ytmp = diffcode_jacobian_par(list(flowsheet_residuals,K_coef, cp_1, hh, 0.22),x,nc)';
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
        [J,ytmpconstr] = derivative(list(flowsheet_residuals,K_coef, cp_1, hh, 0.22), x, H_form = "hypermat");
//        ytmpconstr = diffcode_hessian_par(list(flowsheet_residuals,K_coef, cp_1, hh, 0.22),x,nc);

        for i = 1: nc; 
            if lambda(i) <> 0 then

                ysum = ysum + lambda(i)*ytmpconstr(:,:,i); 

            end
        end

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
// these 2 evaluations bellow are equivalent, the first one is supposed to be faster
y = sum(((xmfull(measured)-x(measured)).^2)./var(measured)');
//    y = (xmfull(measured)-x(measured))'*diag(ones(1,length(var(measured)))./var(measured))*(xmfull(measured)-x(measured));
// if user wants to make a reconciliation with all variables, it is necessary to set in the SC89.sce
// red = ones(1,length(xm)) 
endfunction    


function [nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = structure_ammonia(x, xfull, K_coef, cp1_coef, h_hx_coef, frac)
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

//    flow_full    Total flow measurements (or estimates, in case of unmeasured stream)
//    temp_full    Temperature measurements (or estimates, in case of unmeasured compound mass fraction)
//    coef:        Coeficients to enthalpy calculations

// From here on, the problem generation is automatic
// No need to edit below
// Arrange the vectors

xlocal=xfull;

// Jacobian and its structure
//jactst = diffcode_jacobian_par(list(flowsheet_residuals,K_coef, cp1_coef, h_hx_coef, frac),xlocal,83)';
jactst = diffcode_jacobian(list(flowsheet_residuals,K_coef, cp1_coef, h_hx_coef, frac),xlocal)';

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

//hess_constr_tst = diffcode_hessian_par(list(flowsheet_residuals,K_coef, cp1_coef, h_hx_coef, frac),xlocal,83);
hess_constr_tst = diffcode_hessian(list(flowsheet_residuals,K_coef, cp1_coef, h_hx_coef, frac),xlocal);

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
// the variables have upper bounds of 50000
upper = ones(nv,1);
// pressure upper bounds
upper(1:11) = 220*ones(11,1);
// flow upper bounds
upper(12:22) = 1600*ones(11,1);
// temperature upper bounds
upper(23:33) = 480*ones(11,1);
// compounds upper bounds
upper(34:$ - 3) = ones(5*11,1); 
// Q upper bounds
upper($-2:$ -1) = (1.0e10)*ones(2,1);
// eps upper bounds
upper($) = 2000;

// pressure  lower bounds
lower(2:11) = 180*ones(10,1);
lower(1) = -60;

// temperature lower bounds
lower(23:27) = 15*ones(5,1);
lower(28:33) = -60*ones(6,1);


// make some more tight bound constraints for the reactor and for the flash:

lower(27) = 400;
upper(27) = 425;
lower(28) = -45;
upper(28) = -20;

// in the non-linear case, all constraints and variables are non-linear (bilinear)
var_lin_type(1:nv) = 1; // Non-Linear
constr_lin_type (1:nc) = 1; // Non-Linear


// the constraints has lower and upper bound of 0
constr_lhs(1:nc) = 0;
constr_rhs(1:nc) = 0;
// pressure drop in heater 1 (H 1)
constr_lhs(15) = 10;
constr_rhs(15) = 10;
// pressure drop in heater 1 (H 2)
constr_lhs(28) = -1;
constr_rhs(28) = -1;
// flash : sum zi = 1 in top and botton
constr_lhs(39) = 1;
constr_rhs(39) = 1;

constr_lhs(40) = 1;
constr_rhs(40) = 1;
// temperature and pressure in compressor C 2
constr_lhs(64) = -69.2746;
constr_rhs(64) = -69.2746;

constr_lhs(65) = -11;
constr_rhs(65) = -11;
// temperature and pressure in compressor C 1
constr_lhs(71) = 222.0027;
constr_rhs(71) = 222.0027;

constr_lhs(72) = -200;
constr_rhs(72) = -200;
// all other streams : sum zi = 1 
constr_lhs(73:81) = 1;
constr_rhs(73:81) = 1;
//constr_lhs(84) = 4;
//constr_rhs(84) = 4;

endfunction
// This function removes the symmetric coeficients of the Hessian matrix.
// Since ipopt uses only the upper triangular part of the Hessian matrix
// it is necessary to remove the lower triangular part of the matrix, 
// which is the purpose of this function.
// inputs:
// ij :  The 2 column matrix of the hessian structure (first columns is the row indices)
//       and the second column is the column indices)
// outputs
// ijnew: The "filtered" Hessian structure where only the upper triangular part of the Hessian is 
//        considered
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

function [ijnew] = filter_symmetric(ij)
    [isize, jsize] = size(ij);
    ijnew = [];
    count = 1;
    for i =1: isize
            if ij(i,1) <= ij(i,2) then
                ijnew(count,:) = ij(i,:);
                count = count + 1;
              
            end
    end
endfunction

