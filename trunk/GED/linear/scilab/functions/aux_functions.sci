// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// adjustability and detectability
function [adj, detect, V, V_inv, sigma_inv, diag_diag_V, Wbar] = adjust(sigma, jac)
    sigma_inv=inv(sigma);
// variance-covariance matrix: narasimham pg. 178 eq. 7-3
    V=jac*sigma*jac';
    V_inv= inv(V);
    diag_diag_V = diag(diag(V));
// covariance matrix of adjustments: narasimham pg. 183 eq. 7-13
    Wbar=sigma*jac'*V_inv*jac*sigma;
    varx=diag(sigma);
    rowjac=size(jac,1);
    coljac=size(jac,2);
    W=inv(sigma);
    M=[W, jac';jac,zeros(rowjac,rowjac)];
    MI= inv(M);
    // variance of the reconciled variables calculated according to:
    //Heyen, G, E. Marechál, and B. Kalitventzeff. 1996. Sensitivity calculations and variance analysis in plant measurement reconciliation. Computers & Chemical Engineering 20, no. 96: S539-S544. doi:10.1016/0098-1354(96)00099-3. http://linkinghub.elsevier.com/retrieve/pii/0098135496000993.
    for i=1:coljac
        a(i)=sum((MI(i,1:coljac).^2)./varx')
//Narasimhan pg 211 eq. 7.73
        adj(i)= 1-sqrt(a(i))/sqrt(varx(i))
//Narasimhan pg 211 eq. 7.73
        detect(i)=sqrt(1-a(i)/varx(i))
    end

endfunction
function [a,  b] = between(m, tol)
 // return a vector of indexes of vector elements between max values within tolerance

[a, b] = find (m > max(m) - tol & m < max(m) + tol)    
endfunction
function out=rrn(Ncifras,x)
// Source http://stackoverflow.com/questions/202302/rounding-to-an-arbitrary-number-of-significant-digits
// with some modifications to allow vectors as inputs    
// Return Ncifras significative digits for the given x 
[a b] = size(x);
out = zeros(a,b);
for i = 1:a
    for j=1:b
        if (x(i,j) == 0) then,
            out(i,j)=0;
            continue;
        end
        if (x(i,j) < 0 ) then
            x(i,j) = -x(i,j);        
        end
        d = ceil(log(x(i,j)));
        power = Ncifras - int(d);
        magnitude = 10.^ power;
        shifted = round(x(i,j)*magnitude);
        out(i,j) =  shifted/magnitude;
//        signo=sign(x(i,j));
//        x(i,j)=signo*x(i,j);                 // x is now positive
//        pot=round(log(x(i,j))./log(10)); //log base 10 rounded to the nearest integer
////        mprintf('%f\n',pot);
//        out(i,j)=signo*round(x(i,j)./10.^(pot-Ncifras+1))*10.^(pot-Ncifras+1); // include sign
    end
end
return;
endfunction
function a = generateStreamName(sizex, varargin)
    [lhs,rhs]=argn(0);
  // generate a string with the stream names
    a=string('');
    for i = 1:sizex,
        if  rhs >=2 then // vector form
            a(i) = string(i);
        else // names are concatenated in 1 string
            a = a + string(' ') + string(i);               
        end
    end

endfunction

function a = generateEqpName(sufix, sizex)
    // generate a string with the equipment names
    // Do not separate sufix with spaces
    a = string('');
    for i = 1: sizex
        a = a + string(' ') + sufix + string(i);
    end
endfunction
    
