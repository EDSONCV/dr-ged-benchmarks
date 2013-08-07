function [y]=moving(x,m)
//MOVING will compute moving averages of order n (best taken as odd)
//
//Usage: y=moving(x,n[,fun])
//where x 	is the input vector (or matrix) to be smoothed. 
//      m 	is number of points to average over (best odd, but even works)
//      y 	is output vector of same length as x
//      fun  (optional) is a custom function rather than moving averages
//
// Note:if x is a matrix then the smoothing will be done 'vertically'.
// 
//
// Example:
//
// x=randn(300,1);
// plot(x,'g.'); 
// hold on;
// plot(moving(x,7),'k'); 
// plot(moving(x,7,'median'),'r');
// plot(moving(x,7,@(x)max(x)),'b'); 
// legend('x','7pt moving mean','7pt moving median','7pt moving max','location','best')
//
// optimized Aslak Grinsted jan2004
// enhanced Aslak Grinsted Apr2007


if m==1
    y=x;
    return
end
if size(x,1)==1
    x=x';
end
    f=zeros(m,1)+1/m;
    n=size(x,1);
    isodd=bitand(m,1);
    m2=floor(m/2);


    if (size(x,2)==1) then
        y=filter(f,1,x);
        y=y([zeros(1,m2-1+isodd)+m,m:n,zeros(1,m2)+n]);
    else
        disp('error')
    end

    

endfunction


