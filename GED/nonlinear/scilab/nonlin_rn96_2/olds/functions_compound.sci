// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to weighted least squares functions

// gradient of the constraints

function y = dg1(x)
    //compound_residuals is it's own file
    ytmp = diffcode_jacobian(list(compound_residuals,jac,ncomp),x)';

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
        //compound_residuals is it's own file
        ytmpconstr = diffcode_hessian(list(compound_residuals,jac,ncomp),x);

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

function y=obj_full_compounds(x)
// if user wants to make a reconciliation with all variables, it is necessary 
// to set in the te "red" vector in the PXX_comp.sce
    y =  sum(((xmfull(red) - x(red)).^2)./var(red));
// The same written in matrix notation
// y = (xmfull(red)-x(red))'*diag(ones(length(var(red)),1)./var(red))*(xmfull(red)-x(red))

endfunction    

function [nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  = structure_compound(jac,ncompounds, totalFlowMeasured,compoundMeasured)
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
//    jac:                is the flowsheet jacobian matrix (regarding only total flows), 
//                            also known as the incidence matrix.
//    ncompounds          number of compounds
//    totalFlowMeasured   Total flow measurements (or estimates, in case of unmeasured stream)
//    compoundMeasured    Compounds measuremente (or estimates, in case of unmeasured compound mass fraction)
//    

nstreams = length(totalFlowMeasured);    
// From here on, the problem generation is automatic
// No need to edit below
// Arrange the vectors
xlocal=[totalFlowMeasured(:);matrix(compoundMeasured',-1)];

// Jacobian and its structure
jactst = diffcode_jacobian(list(compound_residuals,jac,ncompounds),xlocal)';

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

hess_constr_tst = diffcode_hessian(list(compound_residuals,jac,ncompounds),xlocal);

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
// the flow variables have upper bounds of 50000
upperf = 50000*ones(length(totalFlowMeasured(:)),1);
// the compounds variables have upper bounds of 1
upperc = ones(length(matrix(compoundMeasured',-1)),1);
upper=[ upperf ; upperc];
// in the non-linear case, all constraints and variables are non-linear (bilinear)
var_lin_type(1:nv) = 1; // Non-Linear
constr_lin_type (1:nc) = 1; // Non-Linear

// Mass and compounds balance has a lower and upper bound of zero
// the constraints has lower and upper bound of 0
constr_lhs(1:nc-nstreams) = 0;
constr_rhs(1:nc-nstreams) = 0;

// Compositions in each streams must sum 1
// the constraints has lower and upper bound of 1
constr_lhs(nc-nstreams+1:nc) = 1;
constr_rhs(nc-nstreams+1:nc) = 1;

endfunction
// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// aux functions to weighted least squares functions   
//***************************************************************
// This function removes the symmetric coeficients of the Hessian matrix.
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

