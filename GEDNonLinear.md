## Non-Linear Problems ##

> In this type of problem, we deal with non-linear data reconciliation/GED (with compound balance, enthalpy balance and other generic non-linearity), with all variables measured/unmeasured and known variance/standard deviation data.
> The aim of this GED benchmark is:

  * Non-linear method for DR/GED solving.
  * Test GED challenge problems where the issues of each problem is presented.
  * Test GED with measured and unmeasured streams.
  * Test Robust functions in GED and the appropriate test parameters.
  * Test/present the Overall Power of each problem.

> In this problem, since we use non-linear constraints and robust objective functions, and the derivatives are complicate to express, we use Ipopt optimizer to run the optimization.  We also use automatic differentiation toolbox from Scilab (diffcode) to evaluate the Hessian and Jacobian structure of objective function. The constraints Jacobian and Hessian are evaluated by finite differences method using "derivative" function form Scilab.

Inside "non\_linear" folder, we have 2 subfolders:

  * scilab: the scilab implementation.
  * diagrams: the process flowsheet.

"scilab" folder
  * jacobians
  * functions
  * nonlin\_alhaj08
  * nonlin\_mandel98
  * nonlin\_pf88
  * nonlin\_rn96\_1
  * nonlin\_rn96\_2
  * nonlin\_rn96\_3
  * nonlin\_sw89
  * nonlin\_ammonia

"functions" folder

This folder has implementation of robust and other type of objective functions, the explanation about them can be found the the references

  * absolute
  * cauchy
  * contamined\_normal
  * fair
  * hampel
  * logistic
  * lorenztian
  * quasi\_weighted
  * wls

There are also 2 files used to setup the functions, these files must not be edited by final users.

  * setup\_DR.sce: method selection

Other functions used for GED:

  * qrlinclass.sci: variable classification

Inside each problem folder, we find 4 files:

  * compound\_residuals.sci : evaluate the residuals of the constraints
  * functions\_compound.sci : evaluate the objective function and problem structure
  * jac\_compound\_residuals.sci : evaluate the jacobians of the residuals
  * nonlin\_XXXXX.sce :  the problem itself




## Step by step file description ##

Initial setup: no need to modify.
```
getd('.');
getd('../functions');
clear flow_full_al flow_al comp_full_al comp_al At_al umeas_al fixed_al red_al lower_al upper_al var_lin_type_al constr_lin_type_al constr_lhs_al constr_rhs_al  just_measured_al observ_al non_obs_al spec_cand_al x_sol f_sol lower_al upper_al extra xmfull ncomp var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
```


Problem setup, in this case a nonlinear problem with compoound balance, first we setup the total flowrates
```
flow_full_al =[25	27	22	2	20	24	14	10	10	4	5	7	3	10	5	5];
```

then the compound flow rates
```
comp_full_al = [ 3.645	3.715	2.683	4.580	2.493	2.077	2.105	2.039	2.946	0.000	8.255	3.518	11.315	5.858	8.255	3.460
                 3.153	3.273	3.646	4.781	3.533	2.944	2.436	3.655	3.410	0.000	1.633	3.958	0.366	2.881	1.633	4.128
                 93.202	93.012	93.671	90.638	93.975	94.979	95.460	94.306	93.644	100.000	90.111	92.523	88.319	91.262	90.111	92.412 ]/100
```

Then, we place information about measured and unmeasured variables:
measured the variable value; unmeasured(-1) or fixed(-5)

```
//               1  2   3   4   5   6   7   8   9  10  11  12   13  14 15  16  
flow_al =     [-1	27	22	-1	20	24	14	10	10	4	-1	7	3	10	5	5];

//                 1        2     3         4     5        6      7        8      9        10    11        12    13       14     15      16
comp_al =  [    3.50	3.58	2.50	4.60	2.29	2.50	2.80	2.08	2.50	250.00	8.34	3.40	11.80	5.92	8.34	3.50
                3.00	3.13	3.50	4.80	3.37	3.37	3.10	3.74	3.00	300.00	1.52	4.20	0.00	2.76	1.52	4.00
                -1       -1	    -1       -1     -1       -1     -1       -1     -1       -1     -1      -1      -1       -1     -1       -1 ];

```


The Jacobian of the constraints, in this problem we made a function to deal with the Jacobian generation for the compound balances. For generic Jacobian, user must implement its own function and its its derivatives.

```

//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    
jac = [ 1   -1  0   1   0   0   0   0    0   0   0   0   0     0     0      0    
        0   1   -1  0   0   0   0   0    0   0   -1   0   0     0     0      0    
        0   0   1   -1  -1  0   0   0    0   0   0   0   0     0     0      0    
        0   0   0   0   1   -1  0   0    0   1   0   0   0     0     0      0    
        0   0   0   0   0   1   -1  -1   0   0   0   0   0     0     0      0   
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    
        0   0   0   0   0   0   1   0    -1   -1  0   0   0     0     0      0    
        0   0   0   0   0   0   0   0    0   0   1   -1   -1    0     0      1    
        0   0   0   0   0   0   0   0    0   0   0   1   1     -1     0      0
        0   0   0   0   0   0   0   0    0   0   0   0   0     1     -1      -1];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    


```


Organizing the vector for the constraints residuals
```
xmfull=[flow_full_al(:);matrix(comp_full_al',-1)];
xm=xmfull;
```


The variance configuration
```
sd_flow = [1.325	1.325	1.46	0.916	0.916	1.101	1.04	0.472	0.401	0.207	0.328	0.328	0.052	0.369	0.25	0.385];
sd_comp_12 = [0.389	0.27	0.252	0.475	0.209	0.246	0.29	0.201	0.371	0.415	0.349	0.363	0.871	0.424	0.349	0.515
             0.344	0.353	0.355	0.5	0.329	0.351	0.374	0.52	0.329	0.455	0.419	0.434	0.632	0.665	0.41	0.518]/100;
sd_comp_3 = [93.202	93.012	93.671	90.638	93.975	94.979	95.460	94.306	93.644	100.000	90.111	92.523	88.319	91.262	90.111	92.412]/1000;             
sd_123 = [sd_comp_12; sd_comp_3];
sd = [sd_flow(:);matrix(sd_123',-1)];
// to run with equaly weighted standard deviation, uncomment the line below
//sd = ones(size(xmfull,1),size(xmfull,2));
var = sd.^2;
ncomp = 3 ;
```

Observability/redundancy tests based on information previously given
```
[At_al,umeas_al, fixed_al] = jac_compound_residuals(jac,ncomp,flow_full_al,comp_full_al, flow_al, comp_al);
[red_al, just_measured_al, observ_al, non_obs_al, spec_cand_al] = qrlinclass(At_al,umeas_al)
```

To reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
```
measured_al = setdiff([1:length(xmfull)], umeas_al);
red = measured_al;//
// to reconcile with all variables, comment the line above and uncomment bellow
//red = [1:length(xmfull)];
```

To select the appropriate objective function the user must select the problem folder (e.g. cd nonlin\_ammonia) and edit the selected script file (edit nonlin\_ammonia.sce). After it, the user may choose the objective function type by selection the right objective function number, e.g., to select WLS function:

```
//WLS = 0
// Absolute sum of squares = 1
//Cauchy = 2
//Contamined Normal = 3
//Fair  = 4
//Hampel = 5
//Logistic = 6
//Lorenztian = 7
//Quasi Weighted = 8
// run the configuration functions with the desired objective function type
obj_function_type = 0;
exec ../functions/setup_DR.sce;
// to run robust reconciliation, it is also necessary to choose the function to return the problem structure
```

The ipopt solver needs some information about the problem, such as jacobian and hessian structure
```
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower_al, upper_al, var_lin_type_al, constr_lin_type_al, constr_lhs_al, constr_rhs_al]  = structure_compound(jac,ncomp, flow_full_al,comp_full_al);
//pause
params = init_param();
params = add_param(params,"hessian_approximation","exact");
params = add_param(params,"tol",1e-2);
params = add_param(params,"acceptable_tol",1e-2);

params = add_param(params,"mu_strategy","monotone");
params = add_param(params,"journal_level",5);
params = add_param(params,"fixed_variable_treatment", "relax_bounds");

//according to the original paper, we fix some variables
lower_al(fixed_al) = xmfull(fixed_al);
upper_al(fixed_al) = xmfull(fixed_al);
xm_init = xmfull./var;
```

An now, the optimization itself

```
[x_sol, f_sol, extra] = ipopt(xmfull, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type_al, constr_lin_type_al, constr_rhs_al, constr_lhs_al, lower_al, upper_al, params);
```

Notice that in linear problems, only the measurements, standard deviations and Jacobians must be given, but here, since we can be dealing with generic nonlinear problems, the problem creation must be implemented individually .


## The Problems ##

  * nonlin\_alhaj08 folder
    * nonlin\_alhaj08.sce: the running script
    * compound\_residuals.sci mass and compound balance (constraints)
    * functions\_compound.sci: problem structure (Jacobian, Hessian, etc.)
    * jac\_compound\_residuals.sci: Jacobian evaluation

> Problem flowsheet:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/GED/nonlinear/diagrams/png/nonlin_alhaj08.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/GED/nonlinear/diagrams/png/nonlin_alhaj08.png)


  * nonlin\_amonia folder
    * ammonia\_plant.sce: the running script
    * flowsheet\_res\_s.sci: mass and energy balance (constraints)
    * functions\_ammonia.sci: problem structure (Jacobian, Hessian, etc.)
    * jac\_flowsheet\_residuals.sci: Jacobian evaluation

> Problem flowsheet:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/nonlinear/diagrams/png/flow_ammonia2.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/nonlinear/diagrams/png/flow_ammonia2.png)


  * nonlin\_mandel98 folder
    * nonlin\_mandel98.sce: the running script
    * compound\_residuals.sci mass and compound balance (constraints)
    * functions\_compound.sci: problem structure (Jacobian, Hessian, etc.)
    * jac\_compound\_residuals.sci: Jacobian evaluation

> Problem flowsheet:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/GED/nonlinear/diagrams/png/nonlin_mandel98_.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/GED/nonlinear/diagrams/png/nonlin_mandel98_.png)


  * nonlin\_pf88 folder
    * nonlin\_pf88.sce: the running script
    * pf88\_residuals.sci:constraints
    * pf\_functions\_nl.sci: problem structure (Jacobian, Hessian, etc.)
    * jac\_pf88\_residuals.sci: Jacobian evaluation

> No flowsheet available.

  * nonlin\_rn96\_1 folder
    * nonlin\_rn96\_1.sce: the running script
    * compound\_residuals.sci mass and compound balance (constraints)
    * functions\_compound.sci: problem structure (Jacobian, Hessian, etc.)
    * jac\_compound\_residuals.sci: Jacobian evaluation
    * old: unused, old implementations

> Problem flowsheet:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/nonlinear/diagrams/png/nonlin_rn96_1.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/nonlinear/diagrams/png/nonlin_rn96_1.png)

  * nonlin\_rn96\_2 folder
    * nonlin\_rn96\_2: the running script
    * compound\_residuals.sci mass and compound balance (constraints)
    * functions\_compound.sci: problem structure (Jacobian, Hessian, etc.)
    * jac\_compound\_residuals.sci: Jacobian evaluation
    * old: unused, old implementations

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/nonlinear/diagrams/png/nonlin_rn96_2.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/nonlinear/diagrams/png/nonlin_rn96_2.png)

  * nonlin\_sw89
    * nonlin\_sw89.sce: the running script
    * flowsheet\_residuals.sci: mass and energy balance (constraints)
    * functions\_energy\_balance.sci: problem structure (Jacobian, Hessian, etc.)
    * jac\_flowsheet\_residuals.sci: Jacobian evaluation

> Problem flowsheet

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/nonlinear/diagrams/png/nonlin_sw89.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/nonlinear/diagrams/png/nonlin_sw89.png)


References

  * Alhaj-Dibo, Moustapha, Didier Maquin, and José Ragot. 2008. Data reconciliation: A robust approach using a contaminated distribution. Control Engineering Practice 16, no. 2 (February): 159-170.

  * Mandel, Denis, Ali Abdollahzadeh, Didier Maquin, and José Ragot. 1998. Data reconciliation by inequality balance equilibration: a LMI approach. International Journal of Mineral Processing 53, no. 3 (April): 157-169.

  * C.C. David Pai, Gary D. Fisher.Application of Broyden's Method to Reconciliation of Nonlinearly Constrained Data. AICHE Journal 1988, V 34 No. 5 -p 873-876

  * Rao, R Ramesh, and Shankar Narasimhan. 1996.“Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” Industrial & Engineering Chemistry Research 35:1362-1368.

  * Swartz, C. L. E, 1898. “Data reconciliation for generalized flowsheet applications.” 197th National Meeting of the American Chemical Society, Dallas, TX.  Data extracted from example 5.4 from: Romagnoli, José, and Mabel C. Sánchez. 1999. Process Systems Engineering - Volume 2 - Data Processing and Reconciliation for Chemical Process Operations. Edited by Elsevier.

  * Rao, R Ramesh, and Shankar Narasimhan. 1996.“Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” Industrial & Engineering Chemistry Research 35:1362-1368.

  * Robust functions:
    * Ozyurt and Pike - Comp.  & Chem. Eng. 28, p. 381-402, (2004)
    * Smooth functions according to Gopal and Biegler. AICHE Journal 45(7) 1535-1547 - July 1999
    * Quasi Weighted Robust function, according to Zhang et al. - Comp.  & Chem. Eng.
34, p. 154-162-402, (2010)