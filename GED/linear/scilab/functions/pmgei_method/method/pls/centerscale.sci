function [Xsc,Xbar,Xstd] = centerscale(X,Xbar,Xstd)

// Ouput variables initialisation (not found in input variables)
Xsc=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


if %nargin==1 then
  Xbar = mean(X,"r");
  Xstd = stdev(X,'r');
end;
[n,k] = size(X);
[zi,zj]=find((ones(n,k)*diag(Xstd)) == 0);
//pause
Xsc=(X-ones(n,k)*diag(Xbar))./(ones(n,k)*diag(Xstd));
Xsc(zi,zj) = 0;
endfunction
