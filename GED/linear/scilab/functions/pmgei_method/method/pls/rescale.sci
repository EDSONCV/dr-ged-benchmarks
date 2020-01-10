function [X] = rescale(Xsc,Xbar,Xstd)

// Ouput variables initialisation (not found in input variables)
X=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


[n,k] = size(Xsc);

X=(Xsc.*[ones(n,1)*Xstd])+(ones(n,1)*Xbar);
endfunction
