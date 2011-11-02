function[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(Q1mt, Q1mtbeta, zadj, runsize, jac_col)
//Measurement test
// The choice of Q1/P1 must be choosen to guarantee that all methods has the same AVTI in order to compare the 
// methods with the same overall power basis.

//Q1mt = 0.05;
P1 = 1 - Q1mt;

norm_mt=cdfnor("X",0,1,P1,Q1mt);

//Q1mtbeta = 0.52;
beta_m = (1-((1-Q1mtbeta).^(1/jac_col)));
Q2=beta_m/2;
// The choice of Q2/P2 must be choosen to guarantee that all methods has the same AVTI in order to compare the 
// methods with the same overall power basis.
P2=1-Q2;

norm_mt2=cdfnor("X",0,1,P2,Q2);

//only with random errors
[ge_mt1_i_r,ge_mt1_j_r]=find(zadj(1:runsize,:)>norm_mt);
[ge_mt2_i_r,ge_mt2_j_r]=find(zadj(1:runsize,:)>norm_mt2);
//with gross errors
[ge_mt1_i_ge,ge_mt1_j_ge]=find(zadj(runsize+1:$,:)>norm_mt);
ge_mt1_ij_ge = [ge_mt1_i_ge',ge_mt1_j_ge'];
[ge_mt2_i_ge,ge_mt2_j_ge]=find(zadj(runsize+1:$,:)>norm_mt2);
ge_mt2_ij_ge = [ge_mt2_i_ge',ge_mt2_j_ge'];

avt1_mt1=length(ge_mt1_i_r)./(runsize*jac_col);
avt1_mt2=length(ge_mt2_i_r)./(runsize*jac_col);

op_mt1=[];
op_mt2=[];
for i=0:jac_col-1
    op_mt1(i+1) = length(find(ge_mt1_ij_ge(:,1) >=(i*runsize+1) & ge_mt1_ij_ge(:,1) <=((i+1)*runsize) & ge_mt1_ij_ge(:,2) == (i+1)))./runsize;
    op_mt2(i+1) =  length(find(ge_mt2_ij_ge(:,1) >=(i*runsize+1) & ge_mt2_ij_ge(:,1) <=((i+1)*runsize) & ge_mt2_ij_ge(:,2) == (i+1)))./runsize;
end   
// MEASUREMENT TEST ENDS HERE

endfunction