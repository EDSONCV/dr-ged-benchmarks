function[avt1_nt1, avt1_nt2, op_nt1, op_nt2] = nodal_test(Qnt1, Q2nt1, jac_row, runsize, zr_nt_nodal)
    //Nodal test
    // The choice of Q1/P1 must be choosen to guarantee that all methods has the same AVTI in order to compare the 
    // methods with the same overall power basis.
    //Qnt1=0.1;

    Q1=Qnt1/2;
    beta_r = (1-((1-Q2nt1).^(1/jac_row)));
    Q2=beta_r/2;

    P1=1-Q1;
    P2=1-Q2;
    norm_nt1=cdfnor("X",0,1,P1,Q1);
    norm_nt2=cdfnor("X",0,1,P2,Q2);
    runsizenodal = size(zr_nt_nodal,1)    
    //
//    [ge_nt1_i_r , ge_nt1_j_r] = find(zr_nt_nodal(1:runsize,:) > norm_nt1);
//    [ge_nt2_i_r,  ge_nt2_j_r]=  find(zr_nt_nodal(1:runsize,:) > norm_nt2);
//
//    avt1_nt1=length(ge_nt1_i_r)./(runsize*jac_row);
//    avt1_nt2=length(ge_nt2_i_r)./(runsize*jac_row);
//
//    [ge_nt1_i_nt,ge_nt1_j_nt]=find(zr_nt_nodal(runsize+1:$,:) > norm_nt1);
//    [ge_nt2_i_nt,ge_nt2_j_nt]=find(zr_nt_nodal(runsize+1:$,:) > norm_nt2);
//
//    ge_nt1_ij_nt = [ge_nt1_i_nt', ge_nt1_j_nt'];
//    ge_nt2_ij_nt = [ge_nt2_i_nt', ge_nt2_j_nt'];
//
//    for i=0:jac_row-1
//        op_nt1(i+1) = length(find(ge_nt1_ij_nt(:,1) >=(i*runsize+1) & ge_nt1_ij_nt(:,1) <=((i+1)*runsize) & ge_nt1_ij_nt(:,2) == (i+1)))./runsize;
//        op_nt2(i+1) =  length(find(ge_nt2_ij_nt(:,1) >=(i*runsize+1) & ge_nt2_ij_nt(:,1) <=((i+1)*runsize) & ge_nt2_ij_nt(:,2) == (i+1)))./runsize;
//    end   


Tsup_nt= zeros(runsize*jac_row,1);
Tsupindex_nt= zeros(runsize*jac_row,2);
// The 'low' termination means pure random error
Tsup_nt_low= zeros(runsize,1);
Tsupindex_nt_low= zeros(runsize,2);

zr_nt_nodalT = zr_nt_nodal';

for j = 1 : runsizenodal
    if j <= runsize then
        //catches Tsup for random noise, used in AVTI          
        [Tsup_nt_low(j), Tsupindex_nt_low(j,:)] = max(zr_nt_nodalT(:,j));
    else
        //catches Tsup for gross errors, used in OP
        [Tsup_nt(j - runsize), Tsupindex_nt(j - runsize,:)] = max(zr_nt_nodalT(:,j));
    end
    //     end
end

ge_glr_nt1 = zeros(jac_row);
avt1_nt1 = zeros(1,1);
op_nt1 = zeros(jac_row);
correct_index_nt1 = [];

ge_glr_nt2 = zeros(jac_row);
avt1_nt2 = zeros(1,1);
op_nt2 = zeros(jac_row);
correct_index_nt2 = [];

if length(Tsupindex_nt) > 0 then
    for i = 1 : jac_row
        correct_index_nt1 = [];
        correct_index_nt1 = intersect(find(Tsupindex_nt((i-1)*runsize+1:i*runsize,1) == i), find(Tsup_nt((i-1)*runsize+1:i*runsize,1) >= norm_nt1));       
        correct_index_nt2 = [];
        correct_index_nt2 = intersect(find(Tsupindex_nt((i-1)*runsize+1:i*runsize,1) == i), find(Tsup_nt((i-1)*runsize+1:i*runsize,1) >= norm_nt2));
        // find number of gross errors correctly identified
        ge_glr_nt1(i) =length(correct_index_nt1);
        ge_glr_nt2(i) =length(correct_index_nt2);
        // calculate the mean of gross errors correctly identified
        // Overall Power
        op_nt1(i) = ge_glr_nt1(i)/runsize;
        op_nt2(i) = ge_glr_nt2(i)/runsize;


    end // for
end // if

// Average of Type I error
avt1_nt1 = length(find(Tsup_nt_low >=  norm_nt1))/runsize;
avt1_nt2 = length(find(Tsup_nt_low >=  norm_nt2))/runsize;

endfunction