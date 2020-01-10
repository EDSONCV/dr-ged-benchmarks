function[avt1_mt1, avt1_mt2, op_mt1, op_mt2, ge_mt1_indexu, ge_mt2_indexu, ge_mt1_indexu_low, ge_mt2_indexu_low] = measurement_test3(Q1mt, Q1mtbeta, zadj, runsize, jac_col, varargin)
    [lhs, rhs] = argn(0);
    if rhs > 5 then
        is_multiple = varargin(1);
    else
        is_multiple = 0;
    end
    disp('inside MT3');
    //Measurement test
    // The choice of Q1/P1 must be choosen to guarantee that all methods has the same AVTI in order to compare the 
    // methods with the same overall power basis.

    //Q1mt = 0.05;
    //P1 = 1 - Q1mt;
    Q1 = Q1mt/2;
    P1 = (1 - Q1);
    
    if Q1 == 0 then
        Q1 = %eps/1000;
        P1 = 1 - Q1;
    end
    
    norm_mt=cdfnor("X",0,1,P1,Q1);

    //Q1mtbeta = 0.52;
    beta_m = (1-((1-Q1mtbeta).^(1/jac_col)));
    Q2=beta_m/2;
    // The choice of Q2/P2 must be choosen to guarantee that all methods has the same AVTI in order to compare the 
    // methods with the same overall power basis.
    P2=1-Q2;
    if Q2 == 0 then
        Q2 = %eps/1000;
        P2 = 1 - Q2;
    end
    norm_mt2=cdfnor("X",0,1,P2,Q2);

//    printf('xchi MT1: %f \n', norm_mt);
//    printf('xchi MT2: %f \n', norm_mt2);

    op_mt1=[];
    op_mt2=[];

    //This is a code snippet from ged_GLR
    //The definition of Power is according to Iordache, Mah, Tamhane, AICHE JOURNAL, V 31 No. 7, 1985
    //Eq 42 Pai = P[|zi| >= |zj|, for all j<>i and |zi| > k]
    runsizefinal = size(zadj,1);
    nrun = (runsizefinal/runsize) -1; 
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
    ge_mt1_index=[];
    ge_mt2_index=[];
    if length(Tsupindex_mt) > 0 then
        for i = 1 : nrun
            // find number of gross errors correctly identified
            //pause
            ge_glr_mt1(i) = length(intersect(find(Tsupindex_mt((i-1)*runsize+1:i*runsize,1) == i), find(Tsup_mt((i-1)*runsize+1:i*runsize,1) >= norm_mt)) );
            ge_glr_mt2(i) = length(intersect(find(Tsupindex_mt((i-1)*runsize+1:i*runsize,1) == i), find(Tsup_mt((i-1)*runsize+1:i*runsize,1) >= norm_mt2)) );    

            //find the indexes of measurement error vector which is bellow the test statistics  
            //for multiple gross error
            if is_multiple == 0 then

                //for single gross error
//                 [found_mt1_i,found_mt1_j] = find(zadj((i-1)*runsize+1 + runsize:i*runsize + runsize,:) > norm_mt);      
//                 [found_mt2_i,found_mt2_j] = find(zadj((i-1)*runsize+1 + runsize:i*runsize + runsize,:) > norm_mt2);
//                if length(found_mt1_i) > 0 then
//                    ge_mt1_index =[ge_mt1_index , runsize*i + unique(found_mt1_i)]; 
//                end
//
//                if length(found_mt2_i) > 0 then
//                    ge_mt2_index =[ge_mt2_index , runsize*i + unique(found_mt2_i)]; 
//                end

                found_mt1 = find(zadj((i-1)*runsize+1 + runsize:i*runsize + runsize,i) > norm_mt);      
                found_mt2 = find(zadj((i-1)*runsize+1 + runsize:i*runsize + runsize,i) > norm_mt2);
//                disp('inside MT3')
//                pause
                if length(found_mt1) > 0 then
                    ge_mt1_index =[ge_mt1_index , runsize*i + found_mt1]; 
                end

                if length(found_mt2) > 0 then
                    ge_mt2_index =[ge_mt2_index , runsize*i + found_mt2]; 
                end
                //       
            end  
            // Overall Power
            op_mt1(i) = ge_glr_mt1(i)/runsize;
            op_mt2(i) = ge_glr_mt2(i)/runsize;
        end
        //    pause
    end
    // notice that for multiple errors the ge_mt1_index is just the opposite
    if is_multiple == 1 then
        [found_mt1_i,found_mt1_j] = find(zadj(runsize + 1 :$,:) > norm_mt);      
        [found_mt2_i,found_mt2_j]=  find(zadj(runsize + 1 :$,:) > norm_mt2);
        if length(found_mt1_i) > 0 then
            ge_mt1_index =[runsize + unique(found_mt1_i)]; 
        end

        if length(found_mt2_i) > 0 then
            ge_mt2_index =[runsize + unique(found_mt2_i)]; 
        end
    end
    
    avt1_mt1 = length(find(Tsuplow_mt >=  norm_mt))/runsize;
    avt1_mt2 = length(find(Tsuplow_mt >=  norm_mt2))/runsize;
    //find the indexes of random error vector which exceeds the test statistics and
    // the indexes of gross errors vector which is bellow the test statistics for 
    // measurement bias vector, they need to be removed later
    ge_mt1_indexu = unique(ge_mt1_index);
    ge_mt1_indexu_low = unique(find(Tsuplow_mt >=  norm_mt));
    ge_mt2_indexu = unique(ge_mt2_index);
    ge_mt2_indexu_low = unique(find(Tsuplow_mt >=  norm_mt2));
    disp('before end MT3')
//    pause

endfunction
