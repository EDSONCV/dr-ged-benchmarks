# Data Reconciliation Benchmark Problems #
Data reconciliation are divided in 4 categories, which are the sub folders contained in /trunk/data\_reconciliation folder:

> ## 1. data\_analysis ##

> The aim of this problems is to test data filter (filter on the variables before DR is applied), data filter window selection and influence on the dynamics on the DR results. These problems were built using the xcos (simulink version of scilab, **compatible with scilab version 5.4**).

Inside the data\_analysis folder, we find an example called 2\_tanks. The problem is a 2 tanks in parallel with measured (or unmeasured) total flows. Each tank has its own dynamic behaviour and a bypass valve to close the tanks feeding. The flowsheet is presented below:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2tanks.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2tanks.png)

When the dynamic system starts:
  * At time t0 where all flow is by-passed to streams 4 and 7.
  * At time t1, the bypass valve is closed and the inlet valve of tank 1 is opened.
  * At time t2, the valve to tank 1 is closed and inlet valve to tank 2 is opened.
  * At time t3, the inlet valves of tank 1 and 2 are closed and by-pass valve is opened.
  * Depending on the case, the user can add gross error to the measurement. The time and location can be set by the user.

There are 2 scripts in each folder. The first one is to test the influence of the size of the moving average window filter in the data reconciliation. The second one, uses a default windows size of 30 data points to filter data and perform reconciliation of all measured streams with and without filtering. The results are plotted in several windows (one for each stream).

We have four implementations of this version:

1. The dynamic behaviour is **neglected** and **all streams are measured** (in folder ss\_all\_measured). The block diagram can be checked below:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2_tanks_ss_all_measured_block_diagram.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2_tanks_ss_all_measured_block_diagram.png)

2. The dynamic behaviour is **neglected** and **some streams are unmeasured** (in folder ss\_unmeasured). The block diagram can be checked below:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2_tanks_ss_all_measured_block_diagram.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2_tanks_ss_all_measured_block_diagram.png)

3. The dynamic behaviour is **neglected** and **all streams are measured** but a **gross error is added** at a given time, which can be set-up by user (in folder gross\_error). The user can select the time of gross error and the location by turning on-off the respective switch. The block diagram can be checked below:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2_tanks_ss_all_measured_block_diagram.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2_tanks_ss_all_measured_block_diagram.png)

4. The dynamic behaviour **must be considered** and **all streams are measured** (in folder dynamic). The block diagram can be checked below:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2_tanks_dynamic_block_diagram.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2_tanks_dynamic_block_diagram.png)

5. The dynamic behaviour **must be considered** and **all streams are measured** and a **gross error** is added at a given time, which can be setup by user (in folder dynamic\_gross\_error).  The user can select the time of gross error and the location by turning on-off the respective switch. The block diagram can be checked below:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2_tanks_dynamic_gross_error_block_diagram.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/data_analysis/diagrams/png/2_tanks_dynamic_gross_error_block_diagram.png)

### File Structure ###

At the root folder of data\_analysis folder we find the following files and folders:

  * diagrams: a diagram of the process in the .dia and .png format.
  * scilab: the scilab implementation.
  * spreadsheet: a Excel spreadsheet of the data, results or some data handling necessary to run the model.

**Files inside "scilab" folder
  * dr\_wls\_simple.sci: Data reconciliation using weighted least squares (WLS) formulation.
  * moving.sci: Moving average method for data filtering.**

  * Folders
    * ss\_all\_measured: Steady state with all streams measured.
    * ss\_unmeasured: Steady state with some streams measured and some unmeasured.
    * gross\_error: Steady state with all streams measured and a gross error (position and time can be chosen).
    * dynamic: Different dynamic behaviour for each tank with all streams measured.
    * dynamic\_gross\_error: Different dynamic behaviour for each tank with all streams measured and a gross error (position and time can be chosen).

> Files inside each folder

Since the block diagrams were generated with Scilab 5.4 and Xcos (Simulink similar for Scilab), they are not compatible with older Scilab versions, in this case, a ".sav" file is available and can be loaded instead of running Xcos simulation. If you want to load a previously saved simulation data, set _run\_new_ flag to 0, in the respective ".sce" script, as in the code below:
```
clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
getd('../');
getd('.');
run_new = 0;
```

> #### ss\_all\_measured folder ####

  * tanks\_ss\_noerror.sce: script file to open Xcos diagram, simulate, filter, reconcile data and plot results.
  * tanks\_tst\_window.sce: script to test moving window size.
  * ss\_no\_error\_mat.sav: saved simulation file, in case the user have problems with Xcos.
  * steady\_state\_no\_ge.zcos: Xcos diagram
  * wls\_linear\_functions.sce: unused by now (will be used in the future)
  * wls\_linear\_functions.sci:unused by now (will be used in the future)

> #### ss\_unmasured folder ####

In this case a simulation with all variables is carried out and a filter is used just like the case before. The difference is that in DR routine, the user can set tu unmeasured streams (in the base case 1, 2, 3 and 4). And a data reconciliation is performed (a little bit slow because we use an optimizer). The files are listed below:

  * tanks\_ss\_umeas.sce: script file to open Xcos diagram, simulate, filter, reconcile data and plot results.
  * tanks\_tst\_window.sce: script to test moving window size.
  * ss\_no\_error\_mat.sav: saved simulation file, in case the user have problems with Xcos.
  * steady\_state\_no\_ge.zcos: Xcos diagram
  * wls\_linear\_functions.sce: used by optimizer (weighted least squares formulation)
  * wls\_linear\_functions.sci:used by optimizer (weighted least squares formulation)

> all the other folders have basically the same structure and will not be described.



> ## 2. error\_in\_variables ##

> Error in Variables (EVM) is a problem where a DR is carried out together with parameter estimation, the aim of this problem is to test algorithms for this problem class. Since all the residuals (constraints) are non-linear, we use automatic differentiation toolbox from Scilab (diffcode) to evaluate the Hessian and Jacobian structure. The Jacobian and Hessian evaluation are done by finite difference (Scilab function "derivative"), and for this reason, some problems takes long to converge.

> We have 4 cases, collected from literature review, to test algorithms and solve strategies.

  * Heat exchanger network from Biegler and Tjoa, 1993. The idea is to reconcile data and estimate heat exchangers parameters (_U.A_)
  * Kinetic constant estimation from a CSTR from Kim, Liebman and Edgar, 1990.
  * Vapour-Liquid equilibria data estimation and reconciliation from Kim, Liebman and Edgar, 1990.
  * Gas phase catalytic hydrogenation of phenol, kinetic constant estimation and DR from Rod and Hancil.

### Folder structure ###


> Each folder has 3 sobfolders:

  * diagrams: a diagram of the process in the .dia and .png format.
  * scilab: the scilab implementation.
  * spreadsheet: a Excel spreadsheet of the data, results or some data handling necessary to run the model.

Inside "scilab" folder we find:

  * ev\_bt93: Heat exchanger network from Biegler and Tjoa.
  * ev\_kle90\_cstr: Kinetic constant estimation from a CSTR from Kim, Liebman and Edgar, 1990.
  * ev\_kle90\_vle: Vapour-Liquid equilibria data estimation and reconciliation from Kim, Liebman and Edgar, 1990
  * ev\_rh80: Gas phase catalytic hydrogenation of phenol, kinetic constant estimation and DR from Rod and Hancil
  * functions:  auxiliary functions for the optimizer.

Folder : ev\_bt93
  * bt93.sce: running script.
  * flowsheet\_residuals.sci: mass and energy balance evaluation (constraints)
  * wls\_functions\_energy\_balance.sci:  optimizer data (objective functions, residuals Jacobian and Hessian, structures, etc).

Since some results take long to solve, we saved some results:

  * x\_sol\_ok\_new\_var\_objfun\_12p\_bfgs.sav
  * x\_sol\_ok\_new\_var\_objfun\_20p\_bfgs.sav
  * x\_sol\_ok\_new\_var\_objfun\_20p\_const\_hess.sav
  * x\_sol\_ok\_new\_var\_objfun\_20p\_exact\_Hess.sav

Folder ev\_kle90\_cstr

  * kle90\_cstr.sce: running script.
  * residuals\_cstr.sci: mass and energy balance evaluation (constraints)
  * wls\_functions\_energy\_CSTR.sci:
optimizer data (objective functions, residuals Jacobian and Hessian, structures, etc).

Folder ev\_kle90\_vle

  * kle90\_vle.sce: running script.
  * residuals.sci: constraints evaluation.
  * wls\_functions\_energy\_ELV\_VanL.sci:
optimizer data (objective functions, residuals Jacobian and Hessian, structures, etc).

Folder ev\_rh80

  * rh80.sce: running script.
  * residuals.sci: constraints evaluation.
  * jac\_residuals.sci: unused by now.
  * wls\_functions\_reactor.sci:optimizer data (objective functions, residuals Jacobian and Hessian, structures, etc).
  * x\_sol\_ok\_new\_var\_objfun.sav: saved results.



### References ###

Heat Exchanger Network

Biegler & Tjoa, 1993
“A parallel implementation for parameter estimation with implicit models”
Annals of Operations Research, V42, [Issue 1](https://code.google.com/p/dr-ged-benchmarks/issues/detail?id=1), pp 1-23 , 1993

CSTR and VLE data

Kim, I-W; Liebman, M J;Edgar, T F
“Robust error- in-variables estimation using nonlinear programming techniques.”
AIChE J. 1990, 36 (7), 985-993.

Gas phase catalytic hydrogenation of phenol

Rod, Vladmir and Hancil, Vladislav
Iterative Estimation of Model Parameters When Measurements of all Variables are Subject to Error. Computers & Chemical Engineering, 1980, v.4, 33-38.



> ## 3. linear ##

> In this type of problem, we deal with linear data reconciliation (total flow measurements), with all streams measured and known variance/standard deviation data.  Although linear problems are easy to solve in DR, since it have analytical solution, the aim of this problem is to test different objective function types, called robust functions. These robust functions can provide results with gross errors removed, however, to do this, some parameters must be tuned in these function. The aim of this class problems is to test these parameters set up.

In this problem, since we use robust functions, and the derivatives are complicate to express, we use Ipopt optimizer to run the optimizations.  We also use automatic differentiation toolbox from Scilab (diffcode) to evaluate the Hessian and Jacobian structure of objective function.

inside "linear" folder, we have 3 sobfolders:

  * diagrams: a diagram of the process in the .dia and .png format.
  * scilab: the scilab implementation.
  * spreadsheet: a Excel spreadsheet of the data, results or some data handling necessary to run the model.

"scilab" folder

  * functions
  * P1
  * P2
  * P3
  * P4
  * P5
  * P6
  * P7
  * P8
  * P9
  * P10
  * P11
  * P12
  * P13
  * P14
  * P15
  * P16

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

  * robust\_structure.sci: method setup
  * setup\_DR.sce: method selection

To select the appropriate objective function the user must select the problem folder (e.g. cd P1) and edit the selected script file (edit P1.sce). After it, the user may choose the objective function type by selection the right objective function number, e.g., to select WLS function:

```
// to run robust reconciliation,, one must choose between the following objective functions to set up the functions path and function parameters:
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

```

Then the file must be saved and executed (exec P1.sce).

The user can also add gross error to a specific Stream. For example, to add a gross error of 9\*sigma to Stream 2 of problem P1, edit P1.sce and remove comments of line "gerror(2) = 9\*sqrt(var(2));" as below:

```
// gross error
gerror = zeros(length(xm),1);
// to setup gross errors, select the stream and magnitude as the line bellow
gerror(2) = 9*sqrt(var(2));
xm = xm + gerror;
```

Then the file must be saved and executed (exec P1.sce).

It is also possible to classify and remove measurements from data reconciliation.
To do this, for example, for problem P1:

-Edit the

`umeas_P1 = [];`

line and place the unmeasured streams inside the brackets . E.g. to make Stream 1 unmeasured:
```
umeas_P1 = [1];
```
```
//observability/redundancy tests                  
umeas_P1 = [];
[red_P1, just_measured_P1, observ_P1, non_obs_P1, spec_cand_P1] = qrlinclass(jac,umeas_P1)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_P1 = setdiff([1:length(xm)], umeas_P1);
red = measured_P1;//
// to reconcile with all variables, comment the line above and uncomment bellow
//red = [1:length(xm)];
```

Then the file must be saved and executed (exec P1.sce).


Problems:

P1:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P1.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P1.png)

P2

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P2.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P2.png)

P3

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P3.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P3.png)

P4

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P4.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P4.png)

P5

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P5.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P5.png)

P6

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P6.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P6.png)

P7

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P7.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P7.png)

P8

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P8.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P8.png)

P9

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P9.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P9.png)

P10

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P10.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P10.png)

P11

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P11.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P11.png)

P12

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P12.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P12.png)

P13

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P13.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P13.png)

P14

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P14.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P14.png)

P15

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P15.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P15.png)

P16

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P16.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/linear/diagrams/png/P16.png)


References:
Reference Index:

P1	Rao, R Ramesh, and Shankar Narasimhan. 1996. “Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” Industrial & Engineering Chemistry Research 35:1362-1368. http://dx.doi.org/10.1021/ie940538b.

P2	Dovì, V G, and C Solisio. 2001. “Reconciliation of censored measurements in chemical processes: an alternative approach.” Chemical Engineering Journal 84:309-314. http://www.sciencedirect.com/science/article/B6TFJ-45KNHB1-F/2/199f358469628f600f10b394d2b55a8b.

P3	Zhang, Zhengjiang, Zhijiang Shao, Xi Chen, Kexin Wang, and Jixin Qian. 2010. “Quasi-weighted least squares estimator for data reconciliation.” Computers & Chemical Engineering 34:154-162. http://www.sciencedirect.com/science/article/B6TFT-4XDCHNS-1/2/63a2b79a4cc89a3afb57ff83c4063242.

P4	Narasimhan, S, and C Jordache. 2000. Data Reconciliation and Gross Error Detection: An Intelligent Use of Process Data. 1st ed. Houston: Gulf Publishing.

P5	Yang, Youqi, Rongbo Ten, and Luiqun Jao. 1995. “A study of gross error detection and data reconciliation in process industries.” Comp. & Chem. Eng 19S:S217-S222.

P6	﻿Rosenberg, J, R S H Mah, and C Iordache. 1987. “Evaluation of Schemes for Detecting and Identifying Gross Errors in Process Data.” Ind. & Eng. Chem. Proc. Des. Dev, Vol. 26: 555-564.

P7	Proposed

P8	Rao, R Ramesh, and Shankar Narasimhan. 1996. “Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” Industrial & Engineering Chemistry Research 35:1362-1368. http://dx.doi.org/10.1021/ie940538b.

P9	Mandel, Denis, Ali Abdollahzadeh, Didier Maquin, and José Ragot. 1998. “Data reconciliation by inequality balance equilibration: a LMI approach.” International Journal of Mineral Processing 53:157-169. http://www.sciencedirect.com/science/article/B6VBN-3VM1X8N-3/2/8bffe94a1153eea8647eed5af0031d36.

P10	Martins, Márcio A.F., Carolina A. Amaro, Leonardo S. Souza, Ricardo A. Kalid, and Asher Kiperstok. 2010. “New objective function for data reconciliation in water balance from industrial processes.” Journal of Cleaner Production 1-6. http://linkinghub.elsevier.com/retrieve/pii/S0959652610001149.

P11	Rao, R Ramesh, and Shankar Narasimhan. 1996. “Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” Industrial & Engineering Chemistry Research 35:1362-1368. http://dx.doi.org/10.1021/ie940538b.

P12	Mitsas, Christos L. 2010. “Data reconciliation and variable classification by null space methods.” Measurement 43:702-707. http://apps.isiknowledge.com/full_record.do?product=UA&search_mode=GeneralSearch&qid=2&SID=2A@bF9dN34I72L1Am9M&page=2&doc=17&colname=WOS (Accessed July 22, 2010).

P13	Alhaj-Dibo, Moustapha, Didier Maquin, and José Ragot. 2008. “Data reconciliation: A robust approach using a contaminated distribution.” Control Engineering Practice 16:159-170. http://www.sciencedirect.com/science/article/B6V2H-4N4406D-1/2/50cac92b050f160a20a795faec990dc7.

P14	Proposed

P15	Serth, R W, and W A Heenan. 1986. “Gross error detection and data reconciliation in steam-metering systems.” AIChE Journal 32:733-747.

P16	Zhang, Zhengjiang, Zhijiang Shao, Xi Chen, Kexin Wang, and Jixin Qian. 2010. “Quasi-weighted least squares estimator for data reconciliation.” Computers & Chemical Engineering 34:154-162. http://www.sciencedirect.com/science/article/B6TFT-4XDCHNS-1/2/63a2b79a4cc89a3afb57ff83c4063242.


Robust functions:

Ozyurt and Pike - Comp.  & Chem. Eng.
> 28, p. 381-402, (2004)

Smooth functions according to Gopal and Biegler
> AICHE Journal 45(7) 1535-1547 - July 1999

Quasi Weighted Robust function, according to Zhang et al. - Comp.  & Chem. Eng.
34, p. 154-162-402, (2010)



> ## 4 non-linear ##

> The aim of these problems is to test DR in more complex and non-linear problems. Some features are tested such as:

> How DR behaves with observable non-measured variables.
> How DR behaves with non-observable non-measured variables.
> How DR behaves with models which are approximated by model simplifications.

> In this problem, since we use non-linear constraints and robust objective functions, and the derivatives are complicate to express, we use Ipopt optimizer to run the optimization.  We also use automatic differentiation toolbox from Scilab (diffcode) to evaluate the Hessian and Jacobian structure of objective function. The constraints Jacobian and Hessian are evaluated by finite differences method using "derivative" function form Scilab.

> "non\_linear" folder have 3 sub-folders:

  * diagrams: a diagram of the process in the .dia and .png format.
  * scilab: the Scilab implementation.
  * spreadsheet: a Excel spreadsheet of the data, results or some data handling necessary to run the model.

> The scilab folder have the following folders:

  * functions: Objective functions and problems structure to communicate with Ipopt
  * nonlin\_ammonia: Ammonia plant
  * nonlin\_pf88:  Pai & Fisher non-linear example
  * nonlin\_rn96\_1: Mineral processing unit
  * nonlin\_rn96\_2: Juice extraction unit
  * nonlin\_sw89: Heat exchanger network

Inside "functions" folder there are the robust function implementation, constraints, Jacobian and gradient evaluation.

  * absolute
  * cauchy
  * contamined\_normal
  * fair
  * hampel
  * logistic
  * lorenztian
  * qrlinclass.sci
  * quasi\_weighted
  * setup\_DR.sce
  * wls

  * nonlin\_amonia folder
    * ammonia\_plant.sce: the running script
    * flowsheet\_res\_s.sci: mass and energy balance (constraints)
    * functions\_ammonia.sci: problem structure (Jacobian, Hessian, etc.)
    * jac\_flowsheet\_residuals.sci: Jacobian evaluation

> Problem flowsheet:

![http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/nonlinear/diagrams/png/flow_ammonia2.png](http://dr-ged-benchmarks.googlecode.com/svn/trunk/data_reconciliation/nonlinear/diagrams/png/flow_ammonia2.png)


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

To select the appropriate objective function the user must select the problem folder (e.g. cd nonlin\_ammonia) and edit the selected script file (edit nonlin\_ammonia.sce). After it, the user may choose the objective function type by selection the right objective function number, e.g., to select WLS function:

```
// to run robust reconciliation,, one must choose between the following objective functions to set up the functions path and function parameters:
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

```

References:

C.C. David Pai, Gary D. Fisher - Application of Broyden's Method to Reconciliation of Nonlinearly Constrained Data - AICHE Journal 1988, V 34 No. 5 -p 873-876

Rao, R Ramesh, and Shankar Narasimhan. 1996.
“Comparison of Techniques for Data Reconciliation of Multicomponent Processes.”
Industrial & Engineering Chemistry Research 35:1362-1368.

Swartz, C. L. E, 1898
“Data reconciliation for generalized flowsheet applications.”
197th National Meeting of the American Chemical Society, Dallas, TX.

Heyen, Marechál, Kalitventzeff - Sensitivity calculations and variance analysis in plant measurement reconciliation - Computers & Chemical Engineering - V. 20 I 96 - 1996 - S539-S544.