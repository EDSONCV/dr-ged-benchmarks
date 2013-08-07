// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

clear xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
getd('.');
getd('../');
run_new = 1;

// setup the constant bias stream
const_stream = 1
    if const_stream == 1 then
    
        xcos('serial_steady_miscal_s1.zcos');
        importXcosDiagram('serial_steady_miscal_s1.zcos');
    
    elseif const_stream == 2  then
        
        xcos('serial_steady_miscal_s2.zcos');
        importXcosDiagram('serial_steady_miscal_s2.zcos');
       

    elseif const_stream == 3 then
        xcos('serial_steady_miscal_s3.zcos');
        importXcosDiagram('serial_steady_miscal_s3.zcos');

        
    elseif const_stream == 4 then
        xcos('serial_steady_miscal_s4.zcos');
        importXcosDiagram('serial_steady_miscal_s4.zcos');

    elseif const_stream == 5 then
        xcos('serial_steady_miscal_s5.zcos');
        importXcosDiagram('serial_steady_miscal_s5.zcos');

    else
        xcos('serial_steady_miscal_s6.zcos');
        importXcosDiagram('serial_steady_miscal_s6.zcos');
                           
    end        
    scicos_simulate(scs_m);

wsize = 30;

filtered_out_sum_1 = moving(serial_outlet_tanks.values(:,1),wsize);
filtered_out_sum_2 = moving(serial_outlet_tanks.values(:,2),wsize);
filtered_out_sum_3 = moving(serial_outlet_tanks.values(:,3),wsize);
filtered_out_sum_4 = moving(serial_outlet_tanks.values(:,4),wsize);
filtered_out_sum_5 = moving(serial_outlet_tanks.values(:,5),wsize);
filtered_in_sum_1 = moving(serial_inlet_tanks.values(:,1),wsize);


xm_full_unfiltered1 = [serial_inlet_tanks.values(:,1), serial_outlet_tanks.values(:,1), serial_outlet_tanks.values(:,2), serial_outlet_tanks.values(:,3), serial_outlet_tanks.values(:,4), serial_outlet_tanks.values(:,5)];

xm_full_filtered1=[filtered_in_sum_1, filtered_out_sum_1, filtered_out_sum_2, filtered_out_sum_3, filtered_out_sum_4,  filtered_out_sum_5];

[xm_f_r1, xm_f_col1] = size(xm_full_filtered1);
[xm_u_r1, xm_u_col1] = size(xm_full_unfiltered1);

xm_full_filtered1 = xm_full_filtered1';
xm_full_unfiltered1 = xm_full_unfiltered1';


sd =[0.5
0.3
1
0.3
0.055
0.055].^(0.5);
sds = sd;
var=sd.^2;
sigma=diag(var);
//The jacobian of the constraints
jac = [ 1   1  -1    0  0   0
        0   -1  1   -1  0   0
        0   0   0    1  -1  -1 ];
// to solve with analytic WLS
[x_sol_unfiltered1]=dr_wls_simple(xm_full_unfiltered1, jac, sigma);
[x_sol_filtered1]=dr_wls_simple(xm_full_filtered1, jac, sigma);
use_subplot = 1;
if use_subplot == 0 then

    i=0;
    j=0;
    dx=800;
    dy=445;
    vsize=800;
    hsize=800

    a1=scf(1);
    subplot(2,3,1);
    set(a1,"figure_size",[vsize,hsize]);
    set(a1,"figure_position",[i*dx,j*dy]);
    //i=i+1;
    plot(serial_outlet_tanks.time,xm_full_filtered1(1,:)','r', serial_outlet_tanks.time,x_sol_filtered1(1,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(1,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(1,:)','yel');
    title("m =" + string(wsize) + " - Reactor Input - Stream 1");
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);
    subplot(2,3,2);
    title('m =' + string(wsize) + ' -Reactor Output - Stream 2');
    plot(serial_outlet_tanks.time,xm_full_filtered1(2,:)','r', serial_outlet_tanks.time,x_sol_filtered1(2,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(2,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(2,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);
    subplot(2,3,3);
    title('m =' + string(wsize) + ' - Reactor Recycle- Stream 3');
    plot(serial_outlet_tanks.time,xm_full_filtered1(3,:)','r', serial_outlet_tanks.time,x_sol_filtered1(3,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(3,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(3,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);
    subplot(2,3,4);
    title('m =' + string(wsize) + ' - Separator input - Stream 4');
    plot(serial_outlet_tanks.time,xm_full_filtered1(4,:)','r', serial_outlet_tanks.time,x_sol_filtered1(4,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(4,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(4,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);    

    subplot(2,3,5);
    title('m =' + string(wsize) + ' - Separator Top - Stream 5');
    plot(serial_outlet_tanks.time,xm_full_filtered1(5,:)','r', serial_outlet_tanks.time,x_sol_filtered1(5,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(5,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(5,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);    

    subplot(2,3,6);
    title('m =' + string(wsize) + ' - Separator Botton - Stream 6')
    plot(serial_outlet_tanks.time,xm_full_filtered1(6,:)','r', serial_outlet_tanks.time,x_sol_filtered1(6,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(6,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(6,:)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);    


    // IF USER FACE PROBLEMS WITH THE SUBPLOT FUNCTION (SUCH AS SCILAB CRASH) IT IS POSSIBLE TO ANALYSE THE CHARTS INDIVIDUALLY, IN THIS CASE, SET THE "USE_SUBPLOT" FLAG TO 0    
else

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
    title("m =" + string(wsize) + " - Reactor Input - Stream 1");
    plot(serial_outlet_tanks.time,xm_full_filtered1(1,:)','r', serial_outlet_tanks.time,x_sol_filtered1(1,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(1,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(1,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);

    a4=scf(4);
    set(a4,"figure_size",[vsize,hsize]);
    set(a4,"figure_position",[i*dx,j*dy]);
  title('m =' + string(wsize) + ' -Reactor Output - Stream 2');
    i=i+1;
    plot(serial_outlet_tanks.time,xm_full_filtered1(2,:)','r', serial_outlet_tanks.time,x_sol_filtered1(2,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(2,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(2,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);

    a5=scf(5);
    set(a5,"figure_size",[vsize,hsize]);
    set(a5,"figure_position",[i*dx,j*dy]);
    i=i+1;
    title('m =' + string(wsize) + ' - Reactor Recycle- Stream 3');
    plot(serial_outlet_tanks.time,xm_full_filtered1(3,:)','r', serial_outlet_tanks.time,x_sol_filtered1(3,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(3,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(3,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);

    a6=scf(6);
    set(a6,"figure_size",[vsize,hsize]);
    set(a6,"figure_position",[i*dx,j*dy]);
    i=0;
    j=j+1;
    title('m =' + string(wsize) + ' - Separator input - Stream 4') 
    plot(serial_outlet_tanks.time,xm_full_filtered1(4,:)','r', serial_outlet_tanks.time,x_sol_filtered1(4,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(4,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(4,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);    

    a7=scf(7);
    set(a7,"figure_size",[vsize,hsize]);
    set(a7,"figure_position",[i*dx,j*dy]);
    i=i+1;
    title('m =' + string(wsize) + ' - Separator Top - Stream 5')
    plot(serial_outlet_tanks.time,xm_full_filtered1(5,:)','r', serial_outlet_tanks.time,x_sol_filtered1(5,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(5,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(5,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);

    a8=scf(8);
    set(a8,"figure_size",[vsize,hsize]);
    set(a8,"figure_position",[i*dx,j*dy]);
    i=i+1;
    title('m =' + string(wsize) + ' - Separator Botton - Stream 6')
    plot(serial_outlet_tanks.time,xm_full_filtered1(6,:)','r', serial_outlet_tanks.time,x_sol_filtered1(6,:)','blu',serial_outlet_tanks.time,xm_full_unfiltered1(6,:)','g', serial_outlet_tanks.time,x_sol_unfiltered1(6,:)','yel')
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",3);

end
//
//
