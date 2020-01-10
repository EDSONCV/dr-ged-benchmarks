function [digse] = digse_eval(xk0, xk)
	digse=-(0.3+log10((1e-12)+abs(xk-xk0)/abs(xk)));
	if digse<0
		digse=0;
    end
endfunction
	
