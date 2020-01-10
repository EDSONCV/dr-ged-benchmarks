function [model]=qpls(X0,Y0,comps,predict_info, compact_info)

//
//Carry out the PLS decomposition of the matrices X and Y using the NIPALS algorithm
//
//INPUTS
//X     	   - inputs matrix to be decomposed 						(required)
//Y          - response matrix to be decomposed						(required)
//comps 	   - numbers of components to be used						(default = m)
//tol   	   - precision of the calculations 							(default = 1e-15)
//maxiter  	- maximum numbers of iterations allowed				(default = 1000)
//
//OUTPUTS
//T - score matrix of X 	(n times comp)
//P - loading matriz of X	(comp times m)
//
//---> explicar os outros no futuro <---
nconverged=0;
nnotconverged=0;
model.method='QPLS';
[%nargout,%nargin] = argn(0)
//extracting information from data and pars structures
[n,k]=mtlb_size(X0);
[n,m]=mtlb_size(Y0);
if %nargin < 3 then
   comps=k;
end

have_best_dirs = %F;
//setting defaults values
//digse_tol=5;
digse_tol=3;
maxiter=3000;

//centering and scaling data
[X0sc,X0bar,X0std]=centerscale(X0);
[Y0sc,Y0bar,Y0std]=centerscale(Y0);

//saving scale constants
model.scale.X0bar=X0bar;
model.scale.X0std=X0std;
model.scale.Y0bar=Y0bar;
model.scale.Y0std=Y0std;
//pause
//#####################################
//STEP A
X=X0sc; Y=Y0sc;
//#####################################

//#####################################
//STEP B
a=0;
//#####################################

u=zeros(size(Y,'r'),1);

for i=1:comps

   a=a+1;
//   pause
   //(1)Set the output scores u as some Y column:
   if  sum(abs(Y(:,1)),'r') == 0 then
       if sum(abs(Y(:,2)),'r') > 0 then
           u=Y(:,2);
           printf('initial output scores at 2, %d \n ', sum(abs(Y(:,2)),'r'));
       end
       if sum(abs(Y(:,3)),'r') > 0 then
           u=Y(:,3);
           printf('initial output scores at 3, %d \n ', sum(abs(Y(:,3)),'r'));
       end
       if sum(abs(Y(:,4)),'r') > 0 then
           u=Y(:,4);
           printf('initial output scores at 4, %d \n ', sum(abs(Y(:,4)),'r'));
       end
       if sum(abs(Y(:,5)),'r') > 0 then
           u=Y(:,5);
           printf('initial output scores at 5, %d \n ', sum(abs(Y(:,5)),'r'));
       end
       if sum(abs(Y(:,6)),'r') > 0 then
           u=Y(:,6);
           printf('initial output scores at 6, %d \n ', sum(abs(Y(:,6)),'r'));
       end  
        if sum(abs(Y(:,7)),'r') > 0 then
           u=Y(:,7);
           printf('initial output scores at 7, %d \n ', sum(abs(Y(:,7)),'r'));
       end  
         if sum(abs(Y(:,8)),'r') > 0 then
           u=Y(:,8);
           printf('initial output scores at 8, %d \n ', sum(abs(Y(:,8)),'r'));
       end 
         if sum(abs(Y(:,9)),'r') > 0 then
           u=Y(:,9);
           printf('initial output scores at 9, %d \n ', sum(abs(Y(:,9)),'r'));
       end   
       if sum(abs(Y(:,10)),'r') > 0 then
           u=Y(:,10);
           printf('initial output scores at 10, %d \n ', sum(abs(Y(:,10)),'r'));
       end 
         if sum(abs(Y(:,11)),'r') > 0 then
           u=Y(:,11);
           printf('initial output scores at 11, %d \n ', sum(abs(Y(:,11)),'r'));
       end  
        if sum(abs(Y(:,12)),'r') > 0 then
           u=Y(:,12);
           printf('initial output scores at 12, %d \n ', sum(abs(Y(:,12)),'r'));
       end  
        if sum(abs(Y(:,13)),'r') > 0 then
           u=Y(:,13);
           printf('initial output scores at 13, %d \n ', sum(abs(Y(:,13)),'r'));
       end  
        if sum(abs(Y(:,14)),'r') > 0 then
           u=Y(:,14);
           printf('initial output scores at 14, %d \n ', sum(abs(Y(:,14)),'r'));
       end  
       
   else     
       u=Y(:,1);
          printf('initial output scores at 1, %d \n ', sum(abs(Y(:,1)),'r') );
   end
//       pause
       
       
   
   
   //(2)Regress columns of X on u:
   w=u'*X/(u'*u);
//pause   
   //(3)Normalise w to unit legth:
   w=w'/norm(w');
   
   //(4)Calculate the input scores  
   t0=X*w/(w'*w);
   
   for iter=1:maxiter
      //(5)Fit the quadratic relation:

      T2=[ones(size(t0,1),1) t0 t0.^2];
      c=inv(T2'*T2)*(T2'*u);
      
      //(6)Calculate the quadratic prediction r of u:
      r=T2*c;
      
      //(7)Regress the columns of Y on r:
      q=r'*Y/(r'*r);
      
      //(8)Normalise q to unit legth:
      q=q'/norm(q');
      
      //(9)Calculate the new output scores:
      u=Y*q/(q'*q);
      
      //(10)Update the input weights:
      Z=[];
      for j=1:k
         Z=[Z (c(2)+2*c(3).*t0).*X(:,j)];
      end
      e=u-r;
//pause
      dw=pinv(Z'*Z)*(Z'*e);
      
      //dw=0.15*dw;
      w=w+dw;
      
      //(11)Normalise w to unit legth:
      w=w/norm(w);
      
      //(12)Calculate new input scores:
      t=X*w/(w'*w);
      
      //(13)check convergence on t:
      //TEST(iter)=norm(t);
      for j=1:k
         digsei(i,1)=digse_eval(w(i)-dw(i),w(i));
      end
      digse=min(digsei); //disp(digse)
      dif=norm(t-t0)/norm(t0);
      DIF(iter)=dif;
      if digse>digse_tol  then //covergrence achived goto 14
         clear DIF
         printf('convergence achieved, %d, %d \n' , digse, iter )
         nconverged = nconverged + 1;
         //pause
         break
      else //covergrence not achived goto 5         
         t0=t;      
  //       printf('convergence NOT achieved, %d, %d \n' , digse, iter)
         nnotconverged=nnotconverged + 1;
      end
      
   end
   
   if iter==maxiter then
     
      clear DIF;
//      printf('the iteration  %d did not converge  \n ', iter);
//      nnotconverged=nnotconverged + 1;
  else 
//      printf('converged at iteration  %d  \n ', iter);
//      nconverged = nconverged + 1;
   end   
   
   //(14)Fit the quadratic relation:
   T2=[ones(size(t,1),1) t t.^2];
   c=inv(T2'*T2)*(T2'*u);
   
   //(15)Calculate the quadratic prediction r of u:
   r=T2*c;
   
   q=r'*Y/(r'*r);
   q=q'/norm(q');
   u=Y*q/(q'*q);
   
   //(16)Calculate the X loadings p
   p=(t'*X)/(t'*t); p=p';
   
   //(17)Calculate the input residual matrix
   X=X-t*p';
   
   //(18)Calculate the output residual matrix   
   Y=Y-r*q';
   
   //#####################################
   //SAVING RESULTS FOR PREDICTION:
   W(i,:)=w'; model.arrays.W=W;
   Q(i,:)=q'; model.arrays.Q=Q;
   P(i,:)=p'; model.arrays.P=P;
   if compact_info == 0 then

        U(:,i)=u; model.arrays.U=U;
        T(:,i)=t; model.arrays.T=T;
        R(:,i)=r; model.arrays.R=R;
   end
   //#####################################
   
   ITER(i,:)=iter; model.stat.ITER=ITER;
   model.pars.C(i,:)=c';
   // calculate best directions
   if minreq(Y-r*q', Y, X, i) == %T then
       model.stat.bestdirs = i-1;
       have_best_dirs = %T;
   end
   

   //axis([-3 3 -3 3])
end

for i=1:comps
   it=model.stat.ITER(i,:);
   par=model.pars.C(i,:);
   //disp([it par])
end
//pause
// calculate R^2
//1-((yhat-y).SumSquare())/((y-ymeanVector).SumSquare());
//if %nargin == 4 & predict_info ==1 then
if predict_info ==1 then
    if have_best_dirs then
        y_hat = predictq(X0, model, model.stat.bestdirs );          
    else
        model.stat.bestdirs = comps;
        y_hat = predictq(X0, model, comps );  
    end
      r_SQR = 1- sum((y_hat - Y0).^2)/sum((Y0 - ones(n,m)*diag(model.scale.Y0bar)).^2);
//  r_SQR = 1- sum((y_hat - Y0).^2)/sum((Y0 - model.scale.Y0bar).^2);
  model.stat.r2 = r_SQR;
end
// calculate mean absolute error
model.stat.mean_absolute_error = sum(abs((y_hat - Y0)))/n;
// calculate max error
model.stat.max_error = max(abs(y_hat - Y0));
model.nconverged = nconverged;
model.nnotconverged = nnotconverged;
endfunction

function[is_best] = minreq(Y1, Y, X, i)
//minreq(Y-r*q.t(), Y, X, i)
	N = size(Y,1);
    M = size(Y,2);
	Yd = Y1 - Y;	// Reduction obtained
	YY = Y'*Y;
	C0 = 2*(YY)./(N-i);	// 2 * residual variance
    bim = [];
	is_best = %F;
	k = %F;
	// ------------ Rule 6 for each response var
	for j = 1: M 
		bim = Yd(:,j)'*Yd(:,j); // Variation obtained
		if  bim(1,1) > C0(j,j) then
			k = %T;
        end
	end
	if (k == %F) then
		is_best = %T;
    end

endfunction    

