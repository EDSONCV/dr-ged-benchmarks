function[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(Qnt1, Q2nt1, jac_row, runsize, zr_nt_nodal_rand)
    //Nodal test
    // The choice of Q1/P1 must be choosen to guarantee that all methods has the same AVTI in order to compare the 
    // methods with the same overall power basis.
    //Qnt1=0.1;

    Q1=Qnt1/2;
    //Q2nt1 = 0.34
    beta_r = (1-((1-Q2nt1).^(1/jac_row)));
    Q2=beta_r/2;

    P1=1-Q1;
    P2=1-Q2;
    norm_nt1=cdfnor("X",0,1,P1,Q1);
    norm_nt2=cdfnor("X",0,1,P2,Q2);
    //
    [ge_nt1_i_r , ge_nt1_j_r] = find(zr_nt_nodal_rand(1:runsize,:) > norm_nt1);
    [ge_nt2_i_r,  ge_nt2_j_r]=  find(zr_nt_nodal_rand(1:runsize,:) > norm_nt2);

    avt1_nt1=length(ge_nt1_i_r)./(runsize*jac_row);
    avt1_nt2=length(ge_nt2_i_r)./(runsize*jac_row);

    [ge_nt1_i_nt,ge_nt1_j_nt]=find(zr_nt_nodal(runsize+1:$,:) > norm_nt1);
    [ge_nt2_i_nt,ge_nt2_j_nt]=find(zr_nt_nodal(runsize+1:$,:) > norm_nt2);

    ge_nt1_ij_nt = [ge_nt1_i_nt', ge_nt1_j_nt'];
    ge_nt2_ij_nt = [ge_nt2_i_nt', ge_nt2_j_nt'];

    for i=0:jac_row-1
        op_nt1(i+1) = length(find(ge_nt1_ij_nt(:,1) >=(i*runsize+1) & ge_nt1_ij_nt(:,1) <=((i+1)*runsize) & ge_nt1_ij_nt(:,2) == (i+1)))./runsize;
        op_nt2(i+1) =  length(find(ge_nt2_ij_nt(:,1) >=(i*runsize+1) & ge_nt2_ij_nt(:,1) <=((i+1)*runsize) & ge_nt2_ij_nt(:,2) == (i+1)))./runsize;
    end   
endfunction