Data Reconciliation Benchmark Problems From Lietrature Review Release 2011
Author: Edson Cordeiro do Valle
Contact - edsoncv@{gmail.com}{vrtech.com.br}
Skype: edson.cv

Linear Problems


"scilab" folder structure

    functions
    P1
    P2
    P3
    P4
    P5
    P6
    P7
    P8
    P9
    P10
    P11
    P12
    P13
    P14
    P15
    P16 

"functions" folder

This folder has implementation of robust and other type of objective functions, the explanation about them can be found the the references

-absolute
-cauchy
-contamined_normal
-fair
-hampel
-logistic
-lorenztian
-quasi_weighted
-wls 

There are also 2 files used to setup the functions, these files must not be edited by final users.

-robust_structure.sci: method setup
-setup_DR.sce: method selection 

To select the appropriate objective function the user must select the problem folder (e.g. cd P1) and edit the selected script file (edit P1.sce). After it, the user may choose the objective function type by selection the right objective function number, e.g., to select WLS function:
{{{
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
}}}
Then the file must be saved and executed (exec P1.sce) 


The user can also add gross error to a specific Stream. For example, to add a gross error of 9*sigma to Stream 2 of problem P1, edit P1.sce and remove comments of line "gerror(2) = 9*sqrt(var(2));" as below:

{{{
// gross error
gerror = zeros(length(xm),1);
// to setup gross errors, select the stream and magnitude as the line bellow
gerror(2) = 9*sqrt(var(2));
xm = xm + gerror;
}}}

It is also possible to classify and remove measurements from data reconciliation.
To do this, for example, for problem P1:
-Edit the "umeas_P1 = [];" line and place the unmeasured streams inside the brackets . E.g. to make Stream 1 unmeasured:

umeas_P1 = [1];

{{{
//observability/redundancy tests                  
umeas_P1 = [];
[red_P1, just_measured_P1, observ_P1, non_obs_P1, spec_cand_P1] = qrlinclass(jac,umeas_P1)

// reconcile with all measured. To reconcile with only redundant variables, uncomment the "red" assignments
measured_P1 = setdiff([1:length(xm)], umeas_P1);
red = measured_P1;//
// to reconcile with all variables, comment the line above and uncomment bellow
//red = [1:length(xm)];
}}}


