The aim of this problems is to test data filter (filter on the variables before DR is applied), data filter window selection and influence on the dynamics on the DR results. These problems were built using the xcos (simulink version of scilab, compatible with scilab version 5.4). 

Inside the data_analysis folder, we find an example called 2_tanks. The problem is a 2 tanks in parallel with measured (or unmeasured) total flows. Each tank has its own dynamic behaviour and a bypass valve to close the tanks feeding.

When the dynamic system starts:

-At time t0 where all flow is by-passed to streams 4 and 7.
-At time t1, the bypass valve is closed and the inlet valve of tank 1 is opened.
-At time t2, the valve to tank 1 is closed and inlet valve to tank 2 is opened.
-At time t3, the inlet valves of tank 1 and 2 are closed and by-pass valve is opened.
-Depending on the case, the user can add gross error to the measurement. The time and location can be set by the user. 

There are 2 scripts in each folder. The first one is to test the influence of the size of the moving average window filter in the data reconciliation. The second one, uses a default windows size of 30 data points to filter data and perform reconciliation of all measured streams with and without filtering. The results are plotted in several windows (one for each stream). 

We have 5 implementations of this version:

1. The dynamic behaviour is neglected and all streams are measured (in folder ss_all_measured). Diagram 2_tanks_ss_all_measured_block_diagram.png.
2. The dynamic behaviour is neglected and some streams are unmeasured (in folder ss_unmeasured). Diagram 2_tanks_ss_all_measured_block_diagram.png.
3. The dynamic behaviour is neglected and all streams are measured but a gross error is added at a given time, which can be set-up by user (in folder gross_error). The user can select the time of gross error and the location by turning on-off the respective switch. Diagram 2_tanks_ss_gross_error_block_diagram.png.
4. The dynamic behaviour must be considered and all streams are measured (in folder dynamic). Diagram 2_tanks_dynamic_block_diagram.png
5. The dynamic behaviour must be considered and all streams are measured and a gross error is added at a given time, which can be setup by user (in folder dynamic_gross_error). The user can select the time of gross error and the location by turning on-off the respective switch. Diagram 2_tanks_dynamic_gross_error_block_diagram.png
