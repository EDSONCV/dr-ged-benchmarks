function[avt1_mt1, avt1_mt2, op_mt1, op_mt2] = measurement_test(Q1mt, Q1mtbeta, zadj, runsize, jac_col)
//Measurement test
// The choice of Q1/P1 must be choosen to guarantee that all methods has the same AVTI in order to compare the 
// methods with the same overall power basis.

//Q1mt = 0.05;
//P1 = 1 - Q1mt;
P1 = (1 - Q1mt/2);

norm_mt=cdfnor("X",0,1,P1,Q1mt/2);

//Q1mtbeta = 0.52;
beta_m = (1-((1-Q1mtbeta).^(1/jac_col)));
Q2=beta_m/2;
// The choice of Q2/P2 must be choosen to guarantee that all methods has the same AVTI in order to compare the 
// methods with the same overall power basis.
P2=1-Q2;

norm_mt2=cdfnor("X",0,1,P2,Q2);

op_mt1=[];
op_mt2=[];

//This is a code snippet from ged_GLR
//The definition of Power is according to Iordache, Mah, Tamhane, AICHE JOURNAL, V 31 No. 7, 1985
//Eq 42 Pai = P[|zi| >= |zj|, for all j<>i and |zi| > k]
runsizefinal = size(zadj,1);
Tsup_mt= zeros(runsizefinal-runsize);
Tsupindex_mt= zeros(runsizefinal-runsize,2);
Tsuplow_mt= zeros(runsize);
Tsupindexlow_mt= zeros(runsize,2);
runsizefinal = size(zadj,1);

zadjT = zadj';
for j = 1 : runsizefinal
    if j <= runsize  then
        //catches Tsup for random noise, used in AVTI
        // here, any location of the max is registered
        [Tsuplow_mt(j), Tsupindexlow_mt(j,:)] = max(zadjT(:,j));

    else

        //catches Tsup for gross errors, used in OP
        // here we are not using the definition of power according to IORDACHE 1985, defined as:
        // according to the definition of power (Iordache, Mah, Tamhane, AICHE JOURNAL, V 31 No. 7, 1985)
        //if 2 or more streams have the same zadj, zadj(i) is identified as gross error, where i is the
        //stream where gross error was added

        [a, b] = between(zadjT(:,j)', 1.0e-6);
        if length(b) > 1 then
            c = -1;
//            To use the definition of Iordache, comment the line above and uncomment the line bellow
//            c = int((j-1)/runsize);    
            Tsup_mt(j - runsize) =  zadjT(a(1),j);
            Tsupindex_mt(j - runsize,:) = c;

        else

            [Tsup_mt(j - runsize), Tsupindex_mt(j - runsize,:)] = max(zadjT(:,j));

        end

    end
end
if length(Tsupindex_mt) > 0 then
    for i = 1 : jac_col
        // find number of gross errors correctly identified
        ge_glr_mt1(i) = length(intersect(find(Tsupindex_mt((i-1)*runsize+1:i*runsize,1) == i), find(Tsup_mt((i-1)*runsize+1:i*runsize,1) >= norm_mt)) );
        ge_glr_mt2(i) = length(intersect(find(Tsupindex_mt((i-1)*runsize+1:i*runsize,1) == i), find(Tsup_mt((i-1)*runsize+1:i*runsize,1) >= norm_mt2)) );
        // Overall Power
        op_mt1(i) = ge_glr_mt1(i)/runsize;
        op_mt2(i) = ge_glr_mt2(i)/runsize;
    end
end
avt1_mt1 = length(find(Tsuplow_mt >=  norm_mt))/runsize;
avt1_mt2 = length(find(Tsuplow_mt >=  norm_mt2))/runsize;

    
endfunction