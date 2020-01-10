// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
function [avti_glr, op_glr_mt, aee_mt, aee_nt, op_glr_nt, avti_glr_nt,varargout ] = calc_GLR(res, V_inv, xfinal, jac, sigma, resGrossErrorNodalRandFi, Qglr_mt, Qglr_nt, runsize)
    jac_row = size(jac,1);
    jac_col= size(jac,2); 
    runsizefinal = size(res,1);
    runsizenodal = size(resGrossErrorNodalRandFi,1);

    // GLR for measurement error dataset
    // the gross error signature for measurement bias
    fi = jac';
    di = zeros(size(fi,1), runsizefinal);
    Ci = zeros(size(fi,1), runsizefinal);
    bi = zeros(size(fi,1), runsizefinal);
    Ti = zeros(size(fi,1), runsizefinal);
    Tsup= zeros(runsizefinal-runsize);
    Tsupindex= zeros(runsizefinal-runsize,2);
    Tsuplow= zeros(runsize);
    Tsupindexlow= zeros(runsize,2);
//    disp(V_inv)
    for j = 1 : runsizefinal
        for i = 1: size(fi,1)
            // Narasimhan pg. 188 eq. 7-31
            di(i,j) = fi(i,:)*V_inv*res(j,:)';
            // Narasimhan pg. 188 eq. 7-32        
            Ci(i,j) = fi(i,:)*V_inv*fi(i,:)';
            // Narasimhan pg. 188 eq. 7-29        
            bi(i,j) = inv(Ci(i,j))*di(i,j);
            // Narasimhan pg. 188 eq. 7-30        
            Ti(i,j) = di(i,j).^2/Ci(i,j) ;
        end
        if j <= runsize  then
            //catches Tsup for random noise, used in AVTI
            [Tsuplow(j), Tsupindexlow(j,:)] = max(Ti(:,j));
        else
 //catches Tsup for gross errors, used in OP

        [a, b] = between(Ti(:,j)', 1.0e-6);
        if length(b) > 1 then
//            c = int((j-1)/runsize);    
            Tsup(j - runsize) =  Ti(a(1),j);
            Tsupindex(j - runsize,:) = -1;

        else

            [Tsup(j - runsize), Tsupindex(j - runsize,:)] = max(Ti(:,j));

        end        end
    end
//    pause
//    disp(Ti)
    //Qglr=0.145;
    //Qglr=0.14;
    // The choice of Qglr must be choosen to guarantee that all methods has the same AVTI in order to compare the 
    // methods with the same overall power basis.
    //Qglr_mt=0.12;
    betaglr = 1 - (1 - Qglr_mt).^(1/size(fi,1));
    xchiglr=cdfchi("X",1, 1-betaglr , betaglr);
    ge_glr = zeros(jac_col);
    ge_avg = zeros(jac_col);
    avti_glr = zeros(jac_col);
    op_glr = zeros(jac_col);
    xfinalT = xfinal';
    aee_mt = zeros(jac_col,1);
    estimate_ge_mt = zeros(jac_col,1);
    ge_avg_glr_abs_mt = zeros(jac_col);
    ge_avg_glr_abs_mt2 = zeros(jac_col);
    ge_avg_glr_mt = zeros(jac_col);
    ge_avg_glr_mt2 = zeros(jac_col);
//    disp(length(Tsupindex))
//    disp(Tsupindex)
    if length(Tsupindex) > 0 then
        for i = 1 : jac_col
            correct_index_mt = [];
            correct_index_mt = intersect(find(Tsupindex((i-1)*runsize+1:i*runsize,1) == i), find(Tsup((i-1)*runsize+1:i*runsize,1) >= xchiglr));
            // find number of gross errors correctly identified
            ge_glr_mt(i) = length(intersect(find(Tsupindex((i-1)*runsize+1:i*runsize,1) == i), find(Tsup((i-1)*runsize+1:i*runsize,1) >= xchiglr)) );
            // calculate the mean of gross errors correctly identified
            if  length(correct_index_mt) <> 0 then

                //ge_avg_glr_mt(i) = mean(bi(i,correct_index_nt + runsize));
                ge_avg_glr_abs_mt(i) = mean(abs(bi(i,correct_index_mt + runsize*i)));            
                ge_avg_glr_abs_mt2(i) = mean(abs(xfinalT(i,correct_index_mt + runsize*i) - xr(i)));          

                // Average error of estimation
                estimate_ge_mt(i) = mean(abs(bi(i,correct_index_mt + runsize*i) ));
                aee_mt(i) = 100*mean(abs((bi(i,correct_index_mt + runsize*i) - (xfinalT(i,correct_index_mt + runsize*i) - xr(i)))./(xfinalT(i,correct_index_mt + runsize*i) - xr(i))));
                // calculate the standrd deviation of gross errors correctly identified
                ge_stdv_glr_mt(i) = stdev(bi(i,correct_index_mt + runsize*i));
            else
                // if there are no elements correctly identified set all elements to "-666"
                ge_avg_glr_mt(i) = -666;
                ge_avg_glr_mt2(i) = -666;
                aee_mt(i) = -666;
                ge_stdv_glr_mt(i) = -666;

            end

            // Overall Power
            op_glr_mt(i) = ge_glr_mt(i)/runsize;
//            clear correct_index_nt;
        end // for
    end //if
//    pause
    // Average of Type I error for measurement random noise
    avti_glr = length(find(Tsuplow >=  xchiglr))/runsize;
//    pause
    // To reduce memory use, can be commented for debug or small problems
    clear di Ci bi Ti Tsup Tsupindex Tsuplow Tsupindexlow;

    // GLR for leaking dataset:
    // It have a different structure from GLR test above
    // The fi2 is the error signature for leakings which is different from 
    // measurement bias
    fi2 = eye(jac_row,jac_row);
    di_nt = zeros(size(fi2,1), runsizenodal);
    Ci_nt = zeros(size(fi2,1), runsizenodal);
    bi_nt = zeros(size(fi2,1), runsizenodal);
    Ti_nt = zeros(size(fi2,1), runsizenodal);

    Tsup_nt= zeros(runsize*jac_row,1);
    Tsupindex_nt= zeros(runsize*jac_row,2);
    // The 'low' termination means pure random error
    Tsup_nt_low= zeros(runsize,1);
    Tsupindex_nt_low= zeros(runsize,2);


    for j = 1 : runsizenodal
        for i = 1: size(fi2,1)
            // Narasimhan pg. 188 eq. 7-31        
            di_nt(i,j) = fi2(i,:)*V_inv*resGrossErrorNodalRandFi(j,:)';        

            // Narasimhan pg. 188 eq. 7-32 
            Ci_nt(i,j) = fi2(i,:)*V_inv*fi2(i,:)';
            // Narasimhan pg. 188 eq. 7-29          
            bi_nt(i,j) = inv(Ci_nt(i,j))*di_nt(i,j);
            // Narasimhan pg. 188 eq. 7-30               
            Ti_nt(i,j) = di_nt(i,j).^2/Ci_nt(i,j) ;        
        end
        if j <= runsize then
            //catches Tsup for random noise, used in AVTI          
            [Tsup_nt_low(j), Tsupindex_nt_low(j,:)] = max(Ti_nt(:,j));
        else
            //catches Tsup for gross errors, used in OP
            [Tsup_nt(j - runsize), Tsupindex_nt(j - runsize,:)] = max(Ti_nt(:,j));
        end
        //     end
    end

    //Qglr=0.145;
    //Qglr=0.14;
    // The choice of Qglr must be choosen to guarantee that all methods has the same AVTI in order to compare the 
    // methods with the same overall power basis.
    //Qglr_nt=0.23;
    betaglr = 1 - (1 - Qglr_nt).^(1/size(fi,1));
    xchiglr_nt=cdfchi("X",1, 1-betaglr , betaglr);
    ge_glr_nt = zeros(jac_row);
    ge_avg_nt = zeros(jac_row);
    avti_glr_nt = zeros(1,1);
    op_glr_nt = zeros(jac_row);
    resGrossErrorNodalRandFT = resGrossErrorNodalRandFi';
    correct_index_nt = [];
    estimate_leak = zeros(jac_row,1);
    aee_nt = zeros(jac_row,1);
    ge_avg_glr_abs_nt = zeros(jac_row,1);
    ge_avg_glr_abs_nt_2 = zeros(jac_row,1);
    if length(Tsupindex_nt) > 0 then
        for i = 1 : jac_row
            correct_index_nt = [];
            correct_index_nt = intersect(find(Tsupindex_nt((i-1)*runsize+1:i*runsize,1) == i), find(Tsup_nt((i-1)*runsize+1:i*runsize,1) >= xchiglr_nt));       
            // find number of gross errors correctly identified
            ge_glr_nt(i) =length(correct_index_nt);
            // calculate the mean of gross errors correctly identified
            if  length(correct_index_nt) <> 0 then
                ge_avg_glr_abs_nt(i) = mean(abs(bi_nt(i,correct_index_nt + i* runsize)));
                ge_avg_glr_abs_nt_2(i) = mean(abs(resGrossErrorNodalRandFT(i,correct_index_nt + i*runsize)));        
                estimate_leak(i) = mean(abs((bi_nt(i,correct_index_nt + runsize*i))))
                aee_nt(i) = 100*mean(abs((bi_nt(i,correct_index_nt + runsize*i) - (resGrossErrorNodalRandFT(i,correct_index_nt + runsize*i)))./(resGrossErrorNodalRandFT(i,correct_index_nt + runsize*i))));

                // calculate the standrd deviation of gross errors correctly identified
                ge_stdv_glr_nt(i) = stdev(bi_nt(i,correct_index_nt + i*runsize));
            else
                ge_avg_glr_nt(i) = -666;
                // Average error of estimation
                aee_nt(i) = -666;
                ge_avg_glr_nt_2(i) = -666;
                // calculate the standard deviation of gross errors correctly identified
                ge_stdv_glr(i) = -666;
            end
            // Overall Power
            op_glr_nt(i) = ge_glr_nt(i)/runsize;
            // Average error of estimation

        end // for
    end // if

    // Average of Type I error
    avti_glr_nt = length(find(Tsup_nt_low >=  xchiglr_nt))/runsize;
    //GLR ENDS HERE
    varargout(1) = estimate_ge_mt;
    varargout(2) = estimate_leak;
endfunction
