function [Y0hat, varargout] =predictq(X0,model,a, Y0)

// Ouput variables initialisation (not found in input variables)
Y0hat=[];

[%nargout,%nargin] = argn(0)


//Extraindo as matrizes da estrutura model:
W=model.arrays.W;
P=model.arrays.P;
Q=model.arrays.Q;
C=model.pars.C;
X0bar=model.scale.X0bar;
X0std=model.scale.X0std;
Y0bar=model.scale.Y0bar;
Y0std=model.scale.Y0std;

//Dimensões das matrizes de dados
[n,k]=size(X0);
[a2,m]=size(Q);

X=centerscale(X0,X0bar,X0std);

if %nargin==2 then
   a=a2;
end

Yhat=zeros(n,m);

for i=1:a
   w=W(i,:)';
   c=C(i,:)';
   q=Q(i,:)';
   p=P(i,:)';   
   
   //extraindo LV1 de X
   t=X*w;
   
   //modelo
   T2=[ones(size(t,1),1) t t.^2];
   r=T2*c;
      
   //atualizando Y
   Yhat=Yhat+r*q';
   
   //resíduos
   X=X-t*p';
   
end;

Y0hat=rescale(Yhat,Y0bar,Y0std);
//pause
if %nargout > 1 then
    pred_error = Y0hat - Y0;
    abs_pred_error = abs(Y0hat - Y0);
    mean_absolute_error = sum(abs((Y0hat - Y0)))/n;
    [max_error, imax_error] = max(abs(Y0hat - Y0));
    varargout(1) = (pred_error);
    varargout(2) = (abs_pred_error);
    varargout(3) = (mean_absolute_error);
    varargout(4) = (imax_error);
    varargout(5) = (max_error);
//    pause
    if %nargout > 6 then
         varargout(6) = struct('pred_error',pred_error,'Y0hat',Y0hat,'mean_absolute_error',mean_absolute_error,'max_error',max_error,'index_of_max_error',imax_error);
        end
        
end
endfunction

   
