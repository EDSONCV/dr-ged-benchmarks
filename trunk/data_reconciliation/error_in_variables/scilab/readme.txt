File structure

Folder : ev_bt93

-bt93.sce: running script.
-flowsheet_residuals.sci: mass and energy balance evaluation (constraints)
-wls_functions_energy_balance.sci: optimizer data (objective functions, residuals Jacobian and Hessian, structures, etc). 

Since some results take long to solve, we saved some results:

-x_sol_ok_new_var_objfun_12p_bfgs.sav
-x_sol_ok_new_var_objfun_20p_bfgs.sav
-x_sol_ok_new_var_objfun_20p_const_hess.sav
-x_sol_ok_new_var_objfun_20p_exact_Hess.sav 

Folder ev_kle90_cstr

-kle90_cstr.sce: running script.
-residuals_cstr.sci: mass and energy balance evaluation (constraints)
-wls_functions_energy_CSTR.sci: 

optimizer data (objective functions, residuals Jacobian and Hessian, structures, etc).

Folder ev_kle90_vle

-kle90_vle.sce: running script.
-residuals.sci: constraints evaluation.
-wls_functions_energy_ELV_VanL.sci: 

optimizer data (objective functions, residuals Jacobian and Hessian, structures, etc).

Folder ev_rh80

-rh80.sce: running script.
-residuals.sci: constraints evaluation.
-jac_residuals.sci: unused by now.
-wls_functions_reactor.sci:optimizer data (objective functions, residuals Jacobian and Hessian, structures, etc).
-x_sol_ok_new_var_objfun.sav: saved results. 
