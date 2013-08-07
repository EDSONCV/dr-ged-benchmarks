function [avti_gt_mt, op_gt_mt, op_gt_nt] = global_test(Q_mt, Q_nt, gamaMeasuremts, runsize, rj, jac_col, jac_row)
    // GLOBAL TEST MEASUREMENT ERRORS
    //Q=0.1;
    P=1-Q_mt;
    xchi=cdfchi("X",rj,P,Q_mt);
//    printf('xchi MT: %d \n', xchi);
    ge_gt_mt=find(gamaMeasuremts(runsize+1:$)>xchi);
    ge_gt_low_mt=length(find(gamaMeasuremts(1:runsize)>xchi));
    avti_gt_mt = ge_gt_low_mt/runsize;

    // Overall Power
    op_gt_mt = length(ge_gt_mt)./(runsize*jac_col);

    // GLOBAL TEST LEAKING
    //Q=0.1;
    P=1-Q_nt;
    xchi=cdfchi("X",rj,P,Q_nt);
//    printf('xchi NT: %d \n', xchi);
    ge_gt_nt=find(gamaNodal(runsize+1:$)> xchi);
    ge_gt_low_nt=length(find(gamaNodal(1:runsize)>xchi))    ;
    avti_gt_nt = ge_gt_low_nt/runsize;
    nge_gt_nt=length(ge_gt_nt)
    // Overall Power
    op_gt_nt = nge_gt_nt./(runsize*jac_row);


    // GLOBAL TEST ENDS HERE

endfunction
function [avti_gt_mt, op_gt_mt, op_gt_nt] = global_test_multi(Q_mt, Q_nt, gamaMeasuremts, runsize, rj, jac_col, jac_row)
    // GLOBAL TEST MEASUREMENT ERRORS
    //Q=0.1;
    P=1-Q_mt;
    xchi=cdfchi("X",rj,P,Q_mt);
//    printf('xchi MT: %d \n', xchi);
    ge_gt_mt=find(gamaMeasuremts(runsize+1:$)>xchi);
    ge_gt_low_mt=length(find(gamaMeasuremts(1:runsize)>xchi));
    avti_gt_mt = ge_gt_low_mt/runsize;

    // Overall Power
    op_gt_mt = length(ge_gt_mt)./(runsize);

    // GLOBAL TEST LEAKING
    //Q=0.1;
    P=1-Q_nt;
    xchi=cdfchi("X",rj,P,Q_nt);
//    printf('xchi NT: %d \n', xchi);
    ge_gt_nt=find(gamaNodal(runsize+1:$)> xchi);
    ge_gt_low_nt=length(find(gamaNodal(1:runsize)>xchi))    ;
    avti_gt_nt = ge_gt_low_nt/runsize;
    nge_gt_nt=length(ge_gt_nt)
    // Overall Power
    op_gt_nt = nge_gt_nt./(runsize);


    // GLOBAL TEST ENDS HERE

endfunction
