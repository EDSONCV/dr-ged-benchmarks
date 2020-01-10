function [avti_gt_mt, avti_gt_nt, op_gt_mt, op_gt_nt, ge_gt_delete_mt, ge_gt_delete_nt, ge_gt_delete_low] = global_test2(Q_mt, Q_nt, gamaMeasuremts, gamaNodal, runsize, rj, jac_col, jac_row)
    // GLOBAL TEST MEASUREMENT ERRORS
    //Q=0.1;
    P=1-Q_mt;
    xchi=cdfchi("X",rj,P,Q_mt);
    printf('xchi MT: %f \n', xchi);
    ge_gt_delete_nt = [];
    ge_gt_delete_mt = [];
    ge_gt_mt=find(gamaMeasuremts(runsize+1:$) >= xchi);
    ge_gt_delete_low=(find(gamaMeasuremts(1:runsize) >= xchi));
    //find the indexes of random error vector which exceeds the test statistics and
    // the indexes of gross errors vector which is bellow the test statistics for 
    // measurement bias vector , they need to be removed later
    found_mt = find(gamaMeasuremts(runsize+1:$) < xchi);
    if length(found_mt > 0) then
            ge_gt_delete_mt=[found_mt + runsize];
    end

    avti_gt_mt = length(ge_gt_delete_low)/runsize;

    // Overall Power
    op_gt_mt = length(ge_gt_mt)./(runsize*jac_col);

    // GLOBAL TEST LEAKING
    //Q=0.1;
    P=1-Q_nt;
    xchi=cdfchi("X",rj,P,Q_nt);
    printf('xchi NT: %d \n', xchi);
    ge_gt_nt=find(gamaNodal(runsize+1:$) >= xchi);
    ge_gt_low_nt=(find(gamaNodal(1:runsize) >= xchi));    
     //find the indexes of random error vector which exceeds the test statistics and
    // the indexes of gross errors vector which is bellow the test statistics for 
    // leaking vector , they need to be removed later.
    // For leaking data set the ge_gt_low_nt will not be concatenated since this 
    //data is already in the gamaMeasuremts dataset
//    ge_gt_delete_nt=[ge_gt_low_nt, find(gamaNodal(runsize+1:$) < xchi) + runsize];
     found_nt = find(gamaNodal(runsize+1:$) < xchi);
     if length(found_nt) > 0 then
              ge_gt_delete_nt=[found_nt + runsize];
     end

//    ge_gt_all_nt=(find(gamaNodal > xchi));
    avti_gt_nt = length(ge_gt_low_nt)/runsize;
    // Overall Power
    op_gt_nt = length(ge_gt_nt)./(runsize*jac_row);
disp('before end GT');
//pause
    // GLOBAL TEST ENDS HERE

endfunction
