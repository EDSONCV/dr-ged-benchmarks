// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
getd('.');
getd('../');
run_new = 0;
// if you face problems with scilab 5.4, load a previoulsy saved result
if run_new == 1 then

    // run the dynamic_no_ge.zcos
    xcos('dynamic_no_ge.zcos');
    importXcosDiagram('dynamic_no_ge.zcos');
    scicos_simulate(scs_m);
    savematfile('-mat','dynamic_no_error_mat.sav', 'dyn_sum_inlet', 'dyn_outlet_tanks','dyn_inlet_tanks','-v6');
else
    loadmatfile('-mat','dynamic_no_error_mat.sav', 'dyn_sum_inlet', 'dyn_outlet_tanks','dyn_inlet_tanks','-v6');
    
end

wsize = 30;

filtered_out_sum_1 = moving(dyn_outlet_tanks.values(:,1),wsize);
filtered_out_sum_2 = moving(dyn_outlet_tanks.values(:,2),wsize);
filtered_out_sum_3 = moving(dyn_outlet_tanks.values(:,3),wsize);
filtered_out_sum_4 = moving(dyn_outlet_tanks.values(:,4),wsize);
filtered_in_sum_1 = moving(dyn_sum_inlet.values,wsize);
filtered_in_sum_2 = moving(dyn_inlet_tanks.values(:,1),wsize);
filtered_in_sum_3 = moving(dyn_inlet_tanks.values(:,2),wsize);
filtered_in_sum_4 = moving(dyn_inlet_tanks.values(:,3),wsize);

xm_full_unfiltered1 = [dyn_sum_inlet.values, dyn_inlet_tanks.values(:,1), dyn_inlet_tanks.values(:,2), dyn_inlet_tanks.values(:,3), dyn_outlet_tanks.values(:,2), dyn_outlet_tanks.values(:,3), dyn_outlet_tanks.values(:,4), dyn_outlet_tanks.values(:,1)];

xm_full_filtered1=[filtered_in_sum_1, filtered_in_sum_2, filtered_in_sum_3, filtered_in_sum_4,filtered_out_sum_2, filtered_out_sum_3, filtered_out_sum_4,  filtered_out_sum_1];

[xm_f_r1, xm_f_col1] = size(xm_full_filtered1);
[xm_u_r1, xm_u_col1] = size(xm_full_unfiltered1);

xm_full_filtered1 = xm_full_filtered1';
xm_full_unfiltered1 = xm_full_unfiltered1';


var = [0.5 0.5 0.3 1 0.25 0.5 1.5 1].^2;
var = var';
sigma=diag(var);
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    
jac = [ 1  -1  -1  -1   0   0   0   0   
0   1   0   0  -1   0   0   0   
0   0   1   0   0   -1  0   0   
0   0   0   1   0   0  -1   0   
0   0   0   0   1   1   1   -1   ];                                
//      1   2   3   4   5   6   7   8

// to solve with analytic WLS
[x_sol_unfiltered1]=dr_wls_simple(xm_full_unfiltered1, jac, sigma);
[x_sol_filtered1]=dr_wls_simple(xm_full_filtered1, jac, sigma);
use_subplot = 1;
if use_subplot == 1 then

    i=0;
    j=0;
    dx=800;
    dy=445;
    vsize=800;
    hsize=800

    a1=scf(1);
    subplot(2,2,1);
    set(a1,"figure_size",[vsize,hsize]);
    set(a1,"figure_position",[i*dx,j*dy]);
    //i=i+1;
    plot(dyn_outlet_tanks.time,xm_full_filtered1(1,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(1,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(1,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(1,:)','yel');
    title("m =" + string(wsize) + " - Input - Stream 1");
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",4);
    subplot(2,2,2);
    title('m =' + string(wsize) + ' - Tank 1 Input - Stream 2');
    plot(dyn_outlet_tanks.time,xm_full_filtered1(2,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(2,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(2,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(2,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);
    subplot(2,2,3);
    title('m =' + string(wsize) + ' - Tank 2 Input - Stream 3');
    plot(dyn_outlet_tanks.time,xm_full_filtered1(3,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(3,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(3,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(3,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);
    subplot(2,2,4);
    title('m =' + string(wsize) + ' - By-pass Input - Stream 4')
    plot(dyn_outlet_tanks.time,xm_full_filtered1(4,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(4,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(4,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(4,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",2);    

    i=i+1
    a2=scf(2);

    subplot(2,2,1);
    set(a2,"figure_size",[vsize,hsize]);
    set(a2,"figure_position",[i*dx,j*dy]);
    title('m =' + string(wsize) + ' - Tank 1 Output - Stream 5');
    plot(dyn_outlet_tanks.time,xm_full_filtered1(5,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(5,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(5,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(5,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);
    subplot(2,2,2);
    title('m =' + string(wsize) + ' - Tank 2 Output - Stream 6');
    plot(dyn_outlet_tanks.time,xm_full_filtered1(6,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(6,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(6,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(6,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);        
    subplot(2,2,3);
    title('m =' + string(wsize) + ' - By-pass Input - Stream 7');
    plot(dyn_outlet_tanks.time,xm_full_filtered1(7,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(7,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(7,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(7,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",2);        
    subplot(2,2,4);
    title('m =' + string(wsize) + ' - Outnput - Stream 8');
    plot(dyn_outlet_tanks.time,xm_full_filtered1(8,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(8,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(8,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(8,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",4);
    // IF USER FACE PROBLEMS WITH THE SUBPLOT FUNCTION (SUCH AS SCILAB CRASH) IT IS POSSIBLE TO ANALYSE THE CHARTS INDIVIDUALLY, IN THIS CASE, COMMENT THE LINES ABOVE AND UNCOMMENT BELLOW
else
    //
    i=0;
    j=0;
    dx=400;
    dy=445;
    vsize=400;
    hsize=440

    a3=scf(3);

    set(a3,"figure_size",[vsize,hsize]);
    set(a3,"figure_position",[i*dx,j*dy]);
    i=i+1;
    title("m =" + string(wsize) + " - Input - Stream 1");
    plot(dyn_outlet_tanks.time,xm_full_filtered1(1,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(1,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(1,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(1,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",4);

    a4=scf(4);
    set(a4,"figure_size",[vsize,hsize]);
    set(a4,"figure_position",[i*dx,j*dy]);
    i=i+1;
    title("m =" + string(wsize) + " - Input - Stream 2");    
    plot(dyn_outlet_tanks.time,xm_full_filtered1(2,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(2,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(2,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(2,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);

    a5=scf(5);
    set(a5,"figure_size",[vsize,hsize]);
    set(a5,"figure_position",[i*dx,j*dy]);
    i=i+1;
    title("m =" + string(wsize) + " - Input - Stream 3");    
    plot(dyn_outlet_tanks.time,xm_full_filtered1(3,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(3,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(3,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(3,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);

    a6=scf(6);
    set(a6,"figure_size",[vsize,hsize]);
    set(a6,"figure_position",[i*dx,j*dy]);
    i=0;
    j=j+1;
    title("m =" + string(wsize) + " - Input - Stream 4");    
    plot(dyn_outlet_tanks.time,xm_full_filtered1(4,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(4,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(4,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(4,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",2);    

    a7=scf(7);
    set(a7,"figure_size",[vsize,hsize]);
    set(a7,"figure_position",[i*dx,j*dy]);
    i=i+1;
    title("m =" + string(wsize) + " - Input - Stream 5");    
    plot(dyn_outlet_tanks.time,xm_full_filtered1(5,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(5,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(5,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(5,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);

    a8=scf(8);
    set(a8,"figure_size",[vsize,hsize]);
    set(a8,"figure_position",[i*dx,j*dy]);
    i=i+1;
    title("m =" + string(wsize) + " - Input - Stream 6");    
    plot(dyn_outlet_tanks.time,xm_full_filtered1(6,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(6,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(6,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(6,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);

    a9=scf(9);
    set(a9,"figure_size",[vsize,hsize]);
    set(a9,"figure_position",[i*dx,j*dy]);
    i=i+1;
    title("m =" + string(wsize) + " - Input - Stream 7");    
    plot(dyn_outlet_tanks.time,xm_full_filtered1(7,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(7,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(7,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(7,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",2);

    a10=scf(10);
    set(a10,"figure_size",[vsize,hsize]);
    set(a10,"figure_position",[i*dx,j*dy]);
    title("m =" + string(wsize) + " - Input - Stream 8");    
    plot(dyn_outlet_tanks.time,xm_full_filtered1(8,:)','r', dyn_outlet_tanks.time,x_sol_filtered1(8,:)','blu',dyn_outlet_tanks.time,xm_full_unfiltered1(8,:)','g', dyn_outlet_tanks.time,x_sol_unfiltered1(8,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",4);
end
//
//
