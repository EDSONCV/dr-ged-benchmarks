The aim of these problems is to test DR in more complex and non-linear problems. Some features are tested such as: 

-How DR behaves with observable non-measured variables.
-How DR behaves with non-observable non-measured variables. 
-How DR behaves with models which are approximated by model simplifications. 

In these problems, since we use non-linear constraints and robust objective functions, and the derivatives are complicate to express, we use Ipopt optimizer to run the optimization. We also use automatic differentiation toolbox from Scilab (diffcode) to evaluate the Hessian and Jacobian structure of objective function. The constraints Jacobian and Hessian are evaluated by finite differences method using "derivative" function form Scilab. 

"non_linear" folder have 3 sub-folders: 

-iise: .ise files for ammonia plant simulation with iiSE process simulator (http://www.vrtech.com.br/en_us/simulador-ise/simulador-de-processo-qu-micos-e-petroqu-micos-ise-3.html)
-diagrams: a diagram of the process in the .dia and .png format.
-scilab: the Scilab implementation.
-spreadsheet: a Excel spreadsheet of the data, results or some data handling necessary to run the model. 
