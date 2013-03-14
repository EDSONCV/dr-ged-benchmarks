Error in Variables (EVM) is a problem where a DR is carried out together with parameter estimation, the aim of this problem is to test algorithms for this problem class. Since all the residuals (constraints) are non-linear, we use automatic differentiation toolbox from Scilab (diffcode) to evaluate the Hessian and Jacobian structure. The Jacobian and Hessian evaluation are done by finite difference (Scilab function "derivative"), and for this reason, some problems takes long to converge. 

We have 4 cases, collected from literature review, to test algorithms and solve strategies. 

-Heat exchanger network from Biegler and Tjoa, 1993. The idea is to reconcile data and estimate heat exchangers parameters (U.A)
-Kinetic constant estimation from a CSTR from Kim, Liebman and Edgar, 1990.
-Vapour-Liquid equilibria data estimation and reconciliation from Kim, Liebman and Edgar, 1990.
-Gas phase catalytic hydrogenation of phenol, kinetic constant estimation and DR from Rod and Hancil. 

Each folder has 3 subfolders: 

-diagrams: a diagram of the process in the .dia and .png format.
-scilab: the scilab implementation.
-spreadsheet: a Excel spreadsheet of the data, results or some data handling necessary to run the model. 

Inside "scilab" folder we find:

-ev_bt93: Heat exchanger network from Biegler and Tjoa.
-ev_kle90_cstr: Kinetic constant estimation from a CSTR from Kim, Liebman and Edgar, 1990.
-ev_kle90_vle: Vapour-Liquid equilibria data estimation and reconciliation from Kim, Liebman and Edgar, 1990
-ev_rh80: Gas phase catalytic hydrogenation of phenol, kinetic constant estimation and DR from Rod and Hancil
-functions: auxiliary functions for the optimizer.

References

Heat Exchanger Network Biegler & Tjoa, 1993 “A parallel implementation for parameter estimation with implicit models” Annals of Operations Research, V42, Issue 1, pp 1-23 , 1993

CSTR and VLE data

Kim, I-W; Liebman, M J;Edgar, T F “Robust error- in-variables estimation using nonlinear programming techniques.” AIChE J. 1990, 36 (7), 985-993.

Gas phase catalytic hydrogenation of phenol

Rod, Vladmir and Hancil, Vladislav Iterative Estimation of Model Parameters When Measurements of all Variables are Subject to Error. Computers & Chemical Engineering, 1980, v.4, 33-38.  
