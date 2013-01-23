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
