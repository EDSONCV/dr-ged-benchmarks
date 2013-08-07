function y=compound_residuals(x, Ainput,Ncomp )
//******************************************************************************
// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//*********************************************************************        
// This function is prepared to use the automatic derivatives toolbox of 
// Scilab. This toolbox can be instaled using the ATOMS installer (package diffcode).
// This function evaluates the residuals of the compound balances for each stream:
// the 'A' matrix, and then concatenates with the normalization equations:
// x_{j,i}, where i= stream and j = compounds 
// Example for a 3 stream system with 3 compounds
// for a simple splitter where the incidence matrix is Ainput = [a11 a12 a13]
// the resulting system is:
//      eq1 = a11.F1.x11 + a12.F2.x12 +a13.F3.x13
//      eq2 = a12.F1.x21 + a12.F2.x22 +a13.F3.x23
//      eq3 = a12.F1.x31 + a12.F2.x32 +a13.F3.x33
// for each stream, we have \sum_{xi}^n x_{i,j} -1 = 0, , resulting in 3 equaitons  
//      eq4 = x11 + x21 + x31 = 1
//      eq5 = x12 + x22 + x32 = 1
//      eq6 = x13 + x23 + x33 = 1
// Notice that the x is a column vector and must pe previoulsy organized
// Outputs:
// y,:        the constraints residuals
// Inputs:
// x:         the column vector of the variables, after the x = [flow, compounds]
//            and x = x(:) operation
// Ainput:    the incidence matrix of the total flow 
// Ncomp:     number of compounds
//  
// get the sizes
//    
[Aeqp, Astreams] =size(Ainput);
// resize the x vector
xx=matrix(x,Astreams,Ncomp+1)
// organize the variables appropriately
TotalFlowMeasured = xx(:,1)';
compoundMeasured  = xx(:, 2:$);  

// Build the functions for the constraints residuals
    for i = 0:Aeqp - 1     
        
//              pause
               A(i*(Ncomp) + 1:(i + 1)*Ncomp ) = sum((ones(Ncomp,1)*Ainput(i+1,:)).*(ones(Ncomp,1)*TotalFlowMeasured.*compoundMeasured'),'c')
    end      
    
// next we'll build the derivatives of the normalization equations:
// for each stream, we have \sum_{j}^streams x_{j,i} -1 = 0, , resulting in "Astreams" equations  
//pause
    Asum = sum (compoundMeasured','r')' ;
y=[A;Asum];
endfunction
