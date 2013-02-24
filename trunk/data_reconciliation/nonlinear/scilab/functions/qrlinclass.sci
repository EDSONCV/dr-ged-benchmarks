// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Classification of a constrained
// input:
// A         :matrix with the jacobian of the constraints
// umeas     :vector with the columns of the A matrix that have the unmeasured 
//       streams
// output:    
// redund         :redundant streams (in the Ax set)
// just_measured  :streams that are just measured but not redundant (in the Ax set)
// observ         :streams that are not measured but are observable by the balance
//                 property
// non_obs        :streams that are not measured and cannot be determined
//spec_candidates :streams that, once measured, will make the non observables, 
//                 observables, if your are planning a simulation, you must specify
//                 this variables  
function [redund, just_measured, observ, non_obs, varargout] =qrlinclass(A,umeas)
    [lhs ,rhs]=argn();
    [row_A, col_A] = size(A);
    Au = A(:,umeas);
    // separate measured, Ax, and unmeasured, Au
    col_vec = [1:col_A];
  
    meas = setdiff(col_vec,umeas);
    Ax = A(:,meas)
    rankAu=rank(Au);
    // perform QR decomposition    
    [Q,R,E]=qr(Au);
    // separating the matrix
    Qu1 = Q(:, 1:rankAu);
    Qu2 = Q(:, rankAu + 1:$);    
    Ru1 = R(1:rankAu, 1:rankAu);
    
    // Gx is the  row vector with information of the redundants:
    //the redundants measurements are the non-zeros of the matrix Gx
    Gx = Qu2'*Ax;
    Gx = sum(abs(Gx), 1);
    // we use a tolerance here    
    //  pause
    if length(Qu2) == 0 then
            redund = [];
            just_measured = gsort(meas, 'r', 'i');
        else
             redund = gsort(meas(find(Gx < -1.0E-7 | Gx > 1.0E-7)), 'r', 'i');
             just_measured = gsort(meas(find(Gx >= -1.0E-7 & Gx <= 1.0E-7)), 'r', 'i');

    end
    // we need to reorder the umeas vector due to the QR factorization which 
    // changed the order of the R matrix
    ordered_umeas = E'*umeas';
    // The first rankAu are candidates to observable (will be tested latter)
    ordered_umeas_up = ordered_umeas(1:rankAu);
    // The last rankAu +1:$ are candidates to specification and are non observable
    ordered_umeas_down = ordered_umeas(rankAu + 1:$);
        
    if rankAu == size(R,2) then
// all unmeasured are observables         
       Ru2 = 0;
       non_obs = [];
       observ = umeas;
       spec_candidates = [];
    else
// we need to find the non observables in the unmeasured set        

        Ru2 = R(1:rankAu, rankAu + 1:$);        
        class_non_obs = sum(abs(inv(Ru1)*Ru2),2);
        // we use a tolerance here

        observ = gsort(ordered_umeas_up(find(class_non_obs >= -1.0E-7 & class_non_obs <= 1.0E-7)), 'r', 'i')';
        non_obs = ordered_umeas_up(find(class_non_obs < -1.0E-7 | class_non_obs > 1.0E-7));
        non_obs = gsort([non_obs; ordered_umeas_down], 'r', 'i')';
        spec_candidates = gsort(ordered_umeas_down, 'r', 'i')';
    end
    if lhs > 4 then
        varargout(1) = gsort(spec_candidates, 'r', 'i');
    end

endfunction
