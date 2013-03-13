// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Biegler & Tjoa, 1993
//“A parallel implementation for parameter estimation with implicit models” 
//Annals of Operations Research, V42, Issue 1, pp 1-23 , 1993
// 

clear xm var jac nc nv nnzjac nnz_hess sparse_dg sparse_dh lower upper var_lin_type constr_lin_type constr_lhs constr_rhs flow_full flow temp_full temp coef
getd('.');
getd('../functions');
stacksize('max');
FA6x=[5.56	5.61	5.51	5.56	5.56	5.56	5.56	5.56	5.56	5.56	5.51	5.61	5.56	5.56	5.56	5.56	5.56	5.56	5.56	5.56];
TA2x=[4.826	4.826	4.826	4.826	4.826	4.822	4.824	4.848	4.826	4.826	4.827	4.826	4.826	4.826	4.829	4.828	4.804	4.826	4.826	4.835];
TA4x=[5.655	5.656	5.65	5.655	5.655	5.608	5.542	5.788	5.655	5.655	5.653	5.66	5.655	5.655	5.696	5.768	5.522	5.655	5.655	5.655];
TA5x=[6.02	6.018	6.025	6.013	6.02	6.077	6.156	5.859	6.064	6.02	6.022	6.015	6.027	6.02	5.97	5.883	6.181	5.976	6.02	6.041];
TA7x=[5.489	5.485	5.492	5.489	5.487	5.488	5.488	5.491	5.489	5.506	5.492	5.485	5.489	5.49	5.489	5.489	5.486	5.489	5.471	5.498];
TA8x=[5.713	5.71	5.72	5.71	5.712	5.737	5.771	5.647	5.723	5.724	5.717	5.707	5.717	5.714	5.693	5.656	5.78	5.695	5.703	5.728];
TD1x=[6.115	6.118	6.113	6.115	6.11	6.118	6.117	6.1 	6.115	6.15	6.113	6.118	6.115	6.121	6.113	6.114	6.131	6.115	6.081	6.125];
TC1x=[6.304	6.299	6.323	6.287	6.304	6.442	6.634	5.915	6.382	6.304	6.308	6.284	6.321	6.304	6.183	5.974	6.692	6.225	6.304	6.342];
TB1x=[6.78	6.784	6.788	6.78	6.78	6.659	6.517	7.062	6.78	6.78	6.775	6.771	6.78	6.78	6.886	7.042	6.497	6.78	6.78	6.767];
TB2x=[5.443	5.444	5.443	5.443	5.443	5.417	5.358	5.547	5.443	5.443	5.442	5.443	5.443	5.443	5.465	5.527	5.338	5.443	5.443	5.448];
FA1x=[9.64	9.69	9.64	9.64	9.64	9.64	9.64	9.64	9.64	9.64	9.59	9.64	9.64	9.64	9.64	9.64	9.64	9.64	9.64	9.65];
FA3x=[4.08	4.08	4.13	4.08	4.08	4.08	4.08	4.08	4.08	4.08	4.08	4.03	4.08	4.08	4.08	4.08	4.08	4.08	4.08	4.09];
FC1x=[3.08	3.08	3.08	3.13	3.08	3.08	3.08	3.08	3.08	3.08	3.08	3.08	3.03	3.08	3.08	3.08	3.08	3.08	3.08	3.09];
FD1x=[6.88	6.88	6.88	6.88	6.93	6.88	6.88	6.88	6.88	6.88	6.88	6.88	6.88	6.83	6.88	6.88	6.88	6.88	6.88	6.89];
FB1x=[2.53	2.53	2.53	2.53	2.53	2.58	2.53	2.53	2.53	2.53	2.53	2.53	2.53	2.53	2.49	2.53	2.53	2.53	2.53	2.54];
TA1x=[4.66	4.66	4.66	4.66	4.66	4.66	4.68	4.66	4.66	4.66	4.66	4.66	4.66	4.66	4.66	4.64	4.66	4.66	4.66	4.67];
TB3x=[4.81	4.81	4.81	4.81	4.81	4.81	4.81	4.83	4.81	4.81	4.81	4.81	4.81	4.81	4.81	4.81	4.79	4.81	4.81	4.82];
TC2x=[5.82	5.82	5.82	5.82	5.82	5.82	5.82	5.82	5.84	5.82	5.82	5.82	5.82	5.82	5.82	5.82	5.82	5.8	    5.82	5.83];
TD2x=[5.58	5.58	5.58	5.58	5.58	5.58	5.58	5.58	5.58	5.6	    5.58	5.58	5.58	5.58	5.58	5.58	5.58	5.58	5.56	5.59];
// unmeasured variables TA3 and TA6, initialized as TA2
TA3x=[4.826	4.826	4.826	4.826	4.826	4.822	4.824	4.848	4.826	4.826	4.827	4.826	4.826	4.826	4.829	4.828	4.804	4.826	4.826	4.835];
TA6x=[4.826	4.826	4.826	4.826	4.826	4.822	4.824	4.848	4.826	4.826	4.827	4.826	4.826	4.826	4.829	4.828	4.804	4.826	4.826	4.835];

ndata = 20

xm =[FA6x(1:ndata);	TA2x(1:ndata);	TA4x(1:ndata);	TA5x(1:ndata);	TA7x(1:ndata);	TA8x(1:ndata);	TD1x(1:ndata);	TC1x(1:ndata);	TB1x(1:ndata);	TB2x(1:ndata);	FA1x(1:ndata);	FA3x(1:ndata);	FC1x(1:ndata);	FD1x(1:ndata);	FB1x(1:ndata);	TA1x(1:ndata);	TB3x(1:ndata);	TC2x(1:ndata);	TD2x(1:ndata); TA3x(1:ndata); TA6x(1:ndata)];
xmrow = xm(:);

UA = [5 4 7 5];

xmfull =[xm(:); UA(:)];    
    
// the standard deviation
// for flow = 0.5, for temperature = 0.01
sflow =  0.5*ones(1,ndata);
stemp = 0.01*ones(1,ndata);
// organize the standard deviation according to the order in xm
stdDevAll= [sflow; stemp; stemp; stemp; stemp; stemp; stemp; stemp; stemp; stemp; sflow; sflow; sflow; sflow; sflow; stemp; stemp; stemp; stemp; stemp; stemp];
//var = var';
stdDevAllrow = stdDevAll(:);
                    
// return the problem structure (jacobian, hessian, number of non-zeros, variable type, etc)
[nc, nv, nnzjac, nnz_hess, sparse_dg, sparse_dh, lower, upper, var_lin_type, constr_lin_type, constr_lhs, constr_rhs]  =  wls_structure_energy_bal(xmfull);

params = init_param();
// We use the given Hessian
//params = add_param(params,"hessian_constant","yes");
//params = add_param(params,"hessian_approximation","exact");
//params = add_param(params,"hessian_approximation","limited-memory");
// uncheck bellow to test derivatives
//params = add_param(params,"derivative_test","second-order");
//params = add_param(params,"derivative_test","first-order");
params = add_param(params,"tol",1e-4);
params = add_param(params,"acceptable_tol",1e-4);
//params = add_param(params,"mu_strategy","monotone");
//params = add_param(params,"limited_memory_max_history",20);
params = add_param(params,"mu_strategy","adaptive");
params = add_param(params,"journal_level",5);
disp('begore start ipopt')
tin = tic();
// if the user want to use random initial guess, uncomment 2 lines bellow and comment the 3rd line
//xrand = rand(30,1);
//[x_sol, f_sol, extra] = ipopt(xrand, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

[x_sol, f_sol, extra] = ipopt(xmfull*1.001, objfun, gradf, confun, dg1, sparse_dg, dh, sparse_dh, var_lin_type, constr_lin_type, constr_rhs, constr_lhs, lower, upper, params);

tend = toc();

mprintf("\n\nSolution: , x\n");
for i = 1 : nv
    mprintf("x[%d] = %e\n", i, x_sol(i));
end

mprintf("\n\nObjective value at optimal point\n");
mprintf("f(x*) = %e\n", f_sol);
