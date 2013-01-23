// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

clear a1 a2 xm var jac nc nv i1 i2 nnzeros sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs
getd('.');
getd('../');
run_new = 0;
// if you face problems with scilab 5.4, load a previoulsy saved result
if run_new ==1 then

    // run the steady_state_no_ge.zcos
    xcos('steady_state_no_ge.zcos');
    importXcosDiagram('steady_state_no_ge.zcos');
    scicos_simulate(scs_m);
    savematfile('-mat','ss_no_error_mat.sav', 'ss_sum_inlet', 'ss_outlet_tanks','ss_inlet_tanks','-v6');
else
    loadmatfile('-mat','ss_no_error_mat.sav', 'ss_sum_inlet', 'ss_outlet_tanks','ss_inlet_tanks','-v6');
    
end

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

wsize = 5;
for k = 1:4

    filtered_out_sum_1 = moving(ss_outlet_tanks.values(:,1),wsize);
    filtered_out_sum_2 = moving(ss_outlet_tanks.values(:,2),wsize);
    filtered_out_sum_3 = moving(ss_outlet_tanks.values(:,3),wsize);
    filtered_out_sum_4 = moving(ss_outlet_tanks.values(:,4),wsize);
    filtered_in_sum_1 = moving(ss_sum_inlet.values,wsize);
    filtered_in_sum_2 = moving(ss_inlet_tanks.values(:,1),wsize);
    filtered_in_sum_3 = moving(ss_inlet_tanks.values(:,2),wsize);
    filtered_in_sum_4 = moving(ss_inlet_tanks.values(:,3),wsize);
   if k ==1  then
    xm_full_unfiltered1 = [ss_sum_inlet.values, ss_inlet_tanks.values(:,1), ss_inlet_tanks.values(:,2), ss_inlet_tanks.values(:,3), ss_outlet_tanks.values(:,2), ss_outlet_tanks.values(:,3), ss_outlet_tanks.values(:,4), ss_outlet_tanks.values(:,1)];
    xm_full_filtered1 = [filtered_in_sum_1, filtered_in_sum_2, filtered_in_sum_3, filtered_in_sum_4,filtered_out_sum_2, filtered_out_sum_3, filtered_out_sum_4,  filtered_out_sum_1];
    [xm_f_r1, xm_f_col1] = size(xm_full_filtered1);
    [xm_u_r1, xm_u_col1] = size(xm_full_unfiltered1);

    [x_sol_unfiltered1]=dr_wls_simple(xm_full_unfiltered1', jac, sigma)';
    [x_sol_filtered1]=dr_wls_simple(xm_full_filtered1', jac, sigma)';
    
       
   else
    xm_full_unfiltered1(:,:,k) = [ss_sum_inlet.values, ss_inlet_tanks.values(:,1), ss_inlet_tanks.values(:,2), ss_inlet_tanks.values(:,3), ss_outlet_tanks.values(:,2), ss_outlet_tanks.values(:,3), ss_outlet_tanks.values(:,4), ss_outlet_tanks.values(:,1)];
    xm_full_filtered1(:,:,k) = [filtered_in_sum_1, filtered_in_sum_2, filtered_in_sum_3, filtered_in_sum_4,filtered_out_sum_2, filtered_out_sum_3, filtered_out_sum_4,  filtered_out_sum_1];
    [x_sol_unfiltered1(:,:,k)]=dr_wls_simple(xm_full_unfiltered1(:,:,k)', jac, sigma)';
    [x_sol_filtered1(:,:,k)]=dr_wls_simple(xm_full_filtered1(:,:,k)', jac, sigma)';
      
    end

//    xm_full_filtered1 = xm_full_filtered1';
//    xm_full_unfiltered1 = xm_full_unfiltered1';

    //plot figures

    i=0;
    j=0;
    dx=800;
    dy=445;
    vsize=800;
    hsize=900;

    a1(k)=scf(2*k);
    subplot(2,2,1);
    set(a1,"figure_size",[vsize,hsize]);
    set(a1,"figure_position",[i*dx,j*dy]);
    //i=i+1;

    plot(ss_outlet_tanks.time,xm_full_filtered1(:,1,k)','r', ss_outlet_tanks.time,x_sol_filtered1(:,1,k)','blu',ss_outlet_tanks.time,xm_full_unfiltered1(:,1,k)','g', ss_outlet_tanks.time,x_sol_unfiltered1(:,1,k)','yel');
    title("m =" + string(wsize) + " - Input - Stream 1");
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",4);
    subplot(2,2,2);
    title('m =' + string(wsize) + ' - Tank 1 Input - Stream 2');
    plot(ss_outlet_tanks.time,xm_full_filtered1(:,2,k)','r', ss_outlet_tanks.time,x_sol_filtered1(:,2,k)','blu',ss_outlet_tanks.time,xm_full_unfiltered1(:,2,k)','g', ss_outlet_tanks.time,x_sol_unfiltered1(:,2,k)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);
    subplot(2,2,3);
    title('m =' + string(wsize) + ' - Tank 2 Input - Stream 3');
    plot(ss_outlet_tanks.time,xm_full_filtered1(:,3,k)','r', ss_outlet_tanks.time,x_sol_filtered1(:,3,k)','blu',ss_outlet_tanks.time,xm_full_unfiltered1(:,3,k)','g', ss_outlet_tanks.time,x_sol_unfiltered1(:,3,k)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);
    subplot(2,2,4);
    title('m =' + string(wsize) + ' - By-pass Input - Stream 4');
    plot(ss_outlet_tanks.time,xm_full_filtered1(:,4,k)','r', ss_outlet_tanks.time,x_sol_filtered1(:,4,k)','blu',ss_outlet_tanks.time,xm_full_unfiltered1(:,4,k)','g', ss_outlet_tanks.time,x_sol_unfiltered1(:,4,k)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",2);    
    i=i+1;
    a2(k)=scf();
    subplot(2,2,1);
    set(a2,"figure_size",[vsize,hsize]);
    set(a2,"figure_position",[i*dx,j*dy]);
    title('m =' + string(wsize) + ' - Tank 1 Output - Stream 5');
    plot(ss_outlet_tanks.time,xm_full_filtered1(:,5,k)','r', ss_outlet_tanks.time,x_sol_filtered1(:,5,k)','blu',ss_outlet_tanks.time,xm_full_unfiltered1(:,5,k)','g', ss_outlet_tanks.time,x_sol_unfiltered1(:,5,k)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);
    subplot(2,2,2);
    title('m =' + string(wsize) + ' - Tank 2 Output - Stream 6');
    plot(ss_outlet_tanks.time,xm_full_filtered1(:,6,k)','r', ss_outlet_tanks.time,x_sol_filtered1(:,6,k)','blu',ss_outlet_tanks.time,xm_full_unfiltered1(:,6,k)','g', ss_outlet_tanks.time,x_sol_unfiltered1(:,6,k)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",1);        
    subplot(2,2,3);
    title('m =' + string(wsize) + ' - By-pass Input - Stream 7');
    plot(ss_outlet_tanks.time,xm_full_filtered1(:,7,k)','r', ss_outlet_tanks.time,x_sol_filtered1(:,7,k)','blu',ss_outlet_tanks.time,xm_full_unfiltered1(:,7,k)','g', ss_outlet_tanks.time,x_sol_unfiltered1(:,7,k)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",2);        
    subplot(2,2,4);
    title('m =' + string(wsize) + ' - Output - Stream 8');
    plot(ss_outlet_tanks.time,xm_full_filtered1(:,8,k)','r', ss_outlet_tanks.time,x_sol_filtered1(:,8,k)','blu',ss_outlet_tanks.time,xm_full_unfiltered1(:,8,k)','g', ss_outlet_tanks.time,x_sol_unfiltered1(:,8,k)','yel');
    legend("meas_filtered", "reconc_filtered", "meas_unfiltered", "reconc_unfiltered",4);
    wsize = 2*wsize;
end;
pause
scf();
plot(ss_outlet_tanks.time,x_sol_filtered1(:,1,1)',ss_outlet_tanks.time,x_sol_filtered1(:,1,2),ss_outlet_tanks.time,x_sol_filtered1(:,1,3),ss_outlet_tanks.time,x_sol_filtered1(:,1,4));
title('Filtered data - Reconciled - Stream 1');
legend("w=5", "w=10", "w=20", "w=40",4);

scf();
plot(ss_outlet_tanks.time,x_sol_filtered1(:,2,1)',ss_outlet_tanks.time,x_sol_filtered1(:,2,2),ss_outlet_tanks.time,x_sol_filtered1(:,2,3),ss_outlet_tanks.time,x_sol_filtered1(:,2,4));
title('Filtered data - Reconciled - Stream 2');
legend("w=5", "w=10", "w=20", "w=40",1);

scf();
plot(ss_outlet_tanks.time,x_sol_filtered1(:,3,1)',ss_outlet_tanks.time,x_sol_filtered1(:,3,2),ss_outlet_tanks.time,x_sol_filtered1(:,3,3),ss_outlet_tanks.time,x_sol_filtered1(:,3,4));
title('Filtered data - Reconciled - Stream 3');
legend("w=5", "w=10", "w=20", "w=40",1);

scf();
plot(ss_outlet_tanks.time,x_sol_filtered1(:,4,1)',ss_outlet_tanks.time,x_sol_filtered1(:,4,2),ss_outlet_tanks.time,x_sol_filtered1(:,4,3),ss_outlet_tanks.time,x_sol_filtered1(:,4,4));
title('Filtered data - Reconciled - Stream 4');
legend("w=5", "w=10", "w=20", "w=40",1);

scf();
plot(ss_outlet_tanks.time,x_sol_filtered1(:,5,1)',ss_outlet_tanks.time,x_sol_filtered1(:,5,2),ss_outlet_tanks.time,x_sol_filtered1(:,5,3),ss_outlet_tanks.time,x_sol_filtered1(:,5,4));
title('Filtered data - Reconciled - Stream 5');
legend("w=5", "w=10", "w=20", "w=40",1);

scf();
plot(ss_outlet_tanks.time,x_sol_filtered1(:,6,1)',ss_outlet_tanks.time,x_sol_filtered1(:,6,2),ss_outlet_tanks.time,x_sol_filtered1(:,6,3),ss_outlet_tanks.time,x_sol_filtered1(:,6,4));
title('Filtered data - Reconciled - Stream 6');
legend("w=5", "w=10", "w=20", "w=40",1);

scf();
plot(ss_outlet_tanks.time,x_sol_filtered1(:,7,1)',ss_outlet_tanks.time,x_sol_filtered1(:,7,2),ss_outlet_tanks.time,x_sol_filtered1(:,7,3),ss_outlet_tanks.time,x_sol_filtered1(:,7,4));
title('Filtered data - Reconciled - Stream 7');
legend("w=5", "w=10", "w=20", "w=40",1);

scf();
plot(ss_outlet_tanks.time,x_sol_filtered1(:,8,1)',ss_outlet_tanks.time,x_sol_filtered1(:,8,2),ss_outlet_tanks.time,x_sol_filtered1(:,8,3),ss_outlet_tanks.time,x_sol_filtered1(:,8,4));
title('Filtered data - Reconciled - Stream 8');
legend("w=5", "w=10", "w=20", "w=40",4);


