File structure:

-dr_wls_simple.sci: Data reconciliation using weighted least squares (WLS) formulation.
-moving.sci: Moving average method for data filtering. 

sub-folders
  -ss_all_measured: Steady state with all streams measured.
  -ss_unmeasured: Steady state with some streams measured and some unmeasured.
  -gross_error: Steady state with all streams measured and a gross error (position and time can be chosen).
  -dynamic: Different dynamic behaviour for each tank with all streams measured.
  -dynamic_gross_error: Different dynamic behaviour for each tank with all streams measured and a gross error (position and time can be chosen). 

Files inside each folder 

Since the block diagrams were generated with Scilab 5.4 and Xcos (Simulink similar for Scilab), they are not compatible with older Scilab versions, in this case, a ".sav" file is available and can be loaded instead of running Xcos simulation. If you want to load a previously saved simulation data, set run_new flag to 0, in the respective ".sce" script, as in the code below:
{{{
clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
getd('../');
getd('.');
run_new = 0;
}}}

   Files in folder ss_all_measured folder

    -tanks_ss_noerror.sce: script file to open Xcos diagram, simulate, filter, reconcile data and plot results.
    -tanks_tst_window.sce: script to test moving window size.
    -ss_no_error_mat.sav: saved simulation file, in case the user have problems with Xcos.
    -steady_state_no_ge.zcos: Xcos diagram
    -wls_linear_functions.sce: unused by now (will be used in the future)
    -wls_linear_functions.sci:unused by now (will be used in the future) 

    Files in folder ss_unmasured folder

In this case a simulation with all variables is carried out and a filter is used just like the case before. The difference is that in DR routine, the user can set tu unmeasured streams (in the base case 1, 2, 3 and 4). And a data reconciliation is performed (a little bit slow because we use an optimizer). The files are listed below:

    -tanks_ss_umeas.sce: script file to open Xcos diagram, simulate, filter, reconcile data and plot results.
    -tanks_tst_window.sce: script to test moving window size.
    -ss_no_error_mat.sav: saved simulation file, in case the user have problems with Xcos.
    -steady_state_no_ge.zcos: Xcos diagram
    -wls_linear_functions.sce: used by optimizer (weighted least squares formulation)
    -wls_linear_functions.sci:used by optimizer (weighted least squares formulation) 

    All the other folders have basically the same structure and will not be described. 
