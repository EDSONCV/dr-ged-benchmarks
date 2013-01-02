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
//    nc_eq               number of equality constraints
//    n_non_lin_eq        number of non-linear equality constraints
//    nv:                 number of variables
//    nnzjac_ineq         number of non zero elements in the Jacobian  of the inequality constraints
//    nnzjac_eq           number of non zero elements in the Jacobian  of the equality constraints
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
//    jac_eq:                is the flowsheet jacobian matrix (regarding only total flows), 
//                            also known as the incidence matrix.
//    n_inequality_contraints since in classical reconciliation problems, no extra inequality constraints are used
//                            it is necessary to tell the number of them
//    x_init                 measurements (used to calculate the hessian and jacobian structures in an initial point
//    iobjfun                objective function 
//    eqfun                  equality constraint function
//    eqfun                  inequality constraint function
//*****************************************************************************

function [nc_eq, n_non_lin_eq, nv, nnzjac_ineq, nnzjac_eq, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = robust_structure(jac_eq, n_inequality_constraints,x_init, iobjfun, eqfun, ineqfun)
    
    // From here on, the problem generation is automatic
// No need to edit below
//The problem size: nc = number of constraints and an number of variables
[nc_eq,nv] = size(jac_eq);

if n_inequality_constraints > 0 then
    // Call Jacobian of inequality constraints to get its structure
    jac_ineqconstraints = diffcode_jacobian(ineqfun,x_init)';
    jac_ineqconstr_sparse = sparse(jac_ineqconstraints);
    [ij_ineq,v_ineq,mn_ineq]=spget(jac_ineqconstr_sparse);
else
    ij_ineq = [];
end

// Jacobian of equality constraints to get its structure
[i1_eq,i2_eq]=find(jac<>0);
// The sparsity structure of the equality constraints

ij_eq = [i1_eq', i2_eq'];


// organizing the indices of the non-zero elements
ij_eq(:,1) = ij_eq(:,1) + n_inequality_constraints;
// The sparsity structure of the constraints
sparse_dg = [ij_ineq; ij_eq];
nnzjac_ineq = size(ij_ineq,1);
nnzjac_eq = size(ij_eq,1);

// The sparsity structure of the Hessian Lagrangian

//first retrieve the constraints Hessian structure, notice that with diffcode_hessian, the
// Hessian has the following dimensions: nvar x nvar x nconstr , so it is in fact a 
// 3 dimensional matrix.
if n_inequality_constraints > 0 then
    hess_ineq_constr = diffcode_hessian(ineqfun,x_init + 10*rand(nv,1 ));
end

//since our jacobian is constant, the hessian is null, in case of a non-linear jacobian
//the lines bellow must be uncommented/commented appropriatelly
//hess_eq_constr = diffcode_hessian(eqfun,x_init);
hess_eq_constr = [];

// cumulative sums the constraints in one hessian
hess_constr = zeros(nv,nv);

//pause
for i = 1: n_inequality_constraints
    // we sum the absolute values to avoid zero cancelation
    hess_constr = hess_constr + abs(hess_ineq_constr(:,:,i));
end

size_hess_eq_constr =  size(hess_eq_constr);

if length(size_hess_eq_constr) > 2  then
    n_non_lin_eq = size_hess_eq_constr(1,3);

    for i = 1: nc_eq
        hess_constr = hess_constr + abs(hess_eq_constr(:,:,i));
    end
else
    n_non_lin_eq = 0;

end

// the Hessian of the objective function
// For Hampel, we are using the finite difference formula due to a limitation of
// diffcode when providing the exact differences of tanh

if obj_function_type == 5 then
    [J,hs_f] = derivative(objfun, x_init , H_form = "hypermat");
    diaghess = diag(hs_f);
    hess_f =  diag(diaghess);
else
//    hess_f = diffcode_hessian(objfun, 100*rand(nv,1 ));
    hess_f = diffcode_hessian(objfun, x_init + 10*rand(nv,1 ));
end

//pause
// sum both of the Hessians
if length(size_hess_eq_constr) > 2  then
    hess_Sparse=sparse(hess_constr + hess_f);
else
    hess_Sparse=sparse(hess_f);
end
// get the hessian structure
[ij_hess,v_hess,mn_hess]=spget(hess_Sparse);

//filters the hessian to remove symmetric indexes
ij_hess_filtered = filter_symmetric(ij_hess);
// index of the non-zero elements of the Hessian
sparse_dh = ij_hess_filtered;
nnz_hess = length(ij_hess_filtered(:,1));

// in case of energy and/or compound balance added, user needs to check these limits!
// the variables have lower bounds of 0
lower = zeros(nv,1);
// the variables have upper bounds of 50000
upper = 50000*ones(nv,1);
// if the user added extra constraints, these lines also need review
var_lin_type(1:nv) = 1; // Non-Linear
constr_lin_type (1:nc_eq) = 0; // Linear
// These bounds must be changed in case of extra constraints are added to the problem 
// (eg. compound balances, energy balances etc)
// the constraints has lower bound of 0
constr_lhs(1:nc_eq) = 0;
// the constraints has upper bound of 0.
constr_rhs(1:nc_eq) = 0;

endfunction

// filter_symmetric removes the symmetric coeficients of the Hessian matrix.
// Since ipopt uses only the upper triangular part of the Hessian matrix
// it is necessary to remove the lower triangular part of the matrix, 
// that is the purpose of this function.
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



