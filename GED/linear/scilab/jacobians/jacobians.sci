// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Rao, R Ramesh, and Shankar Narasimhan. 1996.
//“Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” 
//Industrial & Engineering Chemistry Research 35:1362-1368. 
//http://dx.doi.org/10.1021/ie940538b.
//Bibtex Citation

//@article{Rao1996,
//author = {Rao, R Ramesh and Narasimhan, Shankar},
//isbn = {0888-5885},
//journal = {Industrial \& Engineering Chemistry Research},s
//month = apr,
//number = {4},
//pages = {1362--1368},
//publisher = {American Chemical Society},
//title = {{Comparison of Techniques for Data Reconciliation of Multicomponent Processes}},
//url = {http://dx.doi.org/10.1021/ie940538b},
//volume = {35},
//year = {1996}
//}

// 3 Streams
// 1 Equipment
function [jac]=jacP1()
    jac = [ 1  -1   -1  ];  
endfunction

// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Authors
//Dovì, V G, and C Solisio. 2001. Reconciliation of censored measurements in chemical processes: an alternative approach. Chemical Engineering Journal 84, no. 3 (December): 309-314. http://www.sciencedirect.com/science/article/B6TFJ-45KNHB1-F/2/199f358469628f600f10b394d2b55a8b.

//Bibtex Citation

//@article{Dovi2001,
//annote = { The importance of considering the censoring of measured data in the reconciliation of process flow rates has been shown in a previous paper [Chem. Eng. Sci. 52 (17) (1997) 3047]. The purpose of the present paper is to introduce a new technique for carrying out the actual reconciliation procedure and compare its significance and performance with those of previous methods. A numerical example shows how nontrivial differences are to be expected.},
//author = {Dov\`{\i}, V G and Solisio, C},
//isbn = {1385-8947},
//journal = {Chemical Engineering Journal},
//keywords = {Censored data,Data reconciliation,Detection limits},
//month = dec,
//number = {3},
//pages = {309--314},
//title = {{Reconciliation of censored measurements in chemical processes: an alternative approach}},
//url = {http://www.sciencedirect.com/science/article/B6TFJ-45KNHB1-F/2/199f358469628f600f10b394d2b55a8b},
//volume = {84},
//year = {2001}
//}
// 6 Streams
// 3 Equipments 
function [jac]=jacP2()
jac = [ 1   1  -1    0  0   0
        0   -1  1   -1  0   0
        0   0   0    1  -1  -1 ];
endfunction
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Atmospheric tower example from:
// Zhang, P, G Rong, and Y Wang. 2001.
// A new method of redundancy analysis in data reconciliation and its application.
// Computers and Chemical Engineering 25: 941-949.

//Bibtex Citation

//@article{Zhang2001,
//author = {Zhang, P and Rong, G and Wang, Y},
//journal = {Computers and Chemical Engineering},
//pages = {941--949},
//title = {{A new method of redundancy analysis in data reconciliation and its application}},
//volume = {25},
//year = {2001}
//}

//12 Streams 
//3 Equipments
function [jac]=jacP3()
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  
jac = [ 1  -1   1   0   0   0   0   0    0   0   -1   0  
        0   1   0  -1  -1  -1  -1  -1    1   0    1   -1  
        0   0   -1  0   0   0   0    0   -1  -1   0   1  ];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  

endfunction
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Heat exchanger with by-pass valve
//Narasimhan, S, and C Jordache. 2000.
//Data Reconciliation and Gross Error Detection: An Intelligent Use of Process Data. 1st ed.
//Houston: Gulf Publishing.

//Bibtex Citation

//@book{Narasimhan2000,
//address = {Houston},
//author = {Narasimhan, S and Jordache, C},
//booktitle = {Process Data. Gulf Professional Publishing, Houston, TX.},
//edition = {1},
//publisher = {Gulf Publishing},
//title = {{Data Reconciliation and Gross Error Detection: An Intelligent Use of Process Data}},
//year = {2000}
//}

// 6 Streams
// 4 Equipments 
function [jac]=jacP4()
jac = [ 1  -1  -1    0  0   0
        0   1   0   -1  0   0
        0   0   1    0 -1   0
        0   0   0    1  1  -1]; 
endfunction
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Yang, Youqi, Rongbo Ten, and Luiqun Jao. 1995. 
//“A study of gross error detection and data reconciliation in process industries.” Comp. & Chem. 
//Eng 19S:S217-S222.

//Bibtex Citation

//@article{Yang1995,
//author = {Yang, Youqi and Ten, Rongbo and Jao, Luiqun},
//journal = {Comp. \& Chem. Eng},
//keywords = {combinatory approach,data reconciliation,gross error detection},
//pages = {S217--S222},
//title = {{A study of gross error detection and data reconciliation in process industries}},
//volume = {19S},
//year = {1995}
//}
// 8 Streams
// 4 Equipments 
function [jac]=jacP5()
//The jacobian of the constraints
//      1    2  3    4    5   6   7    8  
jac = [ 1   -1  0    0    0   0   -1   0
        0   1  -1    0    0   0    0   1
        0   0   1    1   -1   0    0   0 
        0   0   0    0    0   -1   1  -1];
//      1    2  3    4    5   6   7    8
endfunction
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Rosenberg, J and Mah, R S H and Iordache, C
//Evaluation of Schemes for Detecting and Identifying Gross Errors in Process Data
//Ind. & Eng. Chem. Proc. Des. Dev, Vol., V. 26:555--564

//Bibtex Citation
//@article{Rosenberg1987a, 
//author = {Rosenberg, J and Mah, R S H and Iordache, C},
//journal = {Ind. \& Eng. Chem. Proc. Des. Dev, Vol.},
//pages = {555--564},
//title = {{Evaluation of Schemes for Detecting and Identifying Gross Errors in Process Data}},
//volume = {26},
//year = {1987}
//}


//7 Streams 
//4 Equipments
//      1   2     3   4   5   6   7   
function [jac]=jacP6()
jac = [ 1   -1    0   1   0   1   0  
        0    1   -1   0   0   0   0 
        0    0    1  -1  -1   0   0
        0    0    0   0   1   -1   -1];                                 
//      1   2   3   4   5   6   7   
endfunction
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Proposed by author
//10 Streams 
//6 Equipments
function [jac]=jacP7()
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10    
jac = [ 1   -1   0  0   0   1   0   0    0   0      
        0   1    -1  0   0   0   0  0    0   0      
        0   0    1  -1   -1   0   0  0   1   0     
        0   0    0  0   1   -1   -1  0    0   0
        0   0    0  0   0   0   1   -1    0   0
        0   0    0  0   0   0   0   1    -1   -1       ];                                
//      1   2   3   4   5   6   7   8    9   10    
endfunction
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Rao, R Ramesh, and Shankar Narasimhan. 1996.
//“Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” 
//Industrial & Engineering Chemistry Research 35:1362-1368. 
//http://dx.doi.org/10.1021/ie940538b.
//Bibtex Citation

//@article{Rao1996,
//author = {Rao, R Ramesh and Narasimhan, Shankar},
//isbn = {0888-5885},
//journal = {Industrial \& Engineering Chemistry Research},
//month = apr,
//number = {4},
//pages = {1362--1368},
//publisher = {American Chemical Society},
//title = {{Comparison of Techniques for Data Reconciliation of Multicomponent Processes}},
//url = {http://dx.doi.org/10.1021/ie940538b},
//volume = {35},
//year = {1996}
//}

// 12 Streams
// 7 Equipments 
function [jac]=jacP8()
//      1   2   3   4   5   6   7   8    9   10  11  
jac = [ 1   -1  -1  0   0   0   0   0    0   0  0    
        0   1   1   -1  -1  -1  -1  0    0   0  0    
        0   0   0   0   1   0   0   0    0   -1 -1    
        0   0   0   1   0   0   0   -1   -1  0  0     ];                                
//      1   2   3   4   5   6   7   8    9   10  11    
endfunction
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Mandel, Denis, Ali Abdollahzadeh, Didier Maquin, and Jos� Ragot. 1998. 
//Data reconciliation by inequality balance equilibration: a LMI approach. 
//International Journal of Mineral Processing 53, no. 3 (April): 157-169. 
//http://www.sciencedirect.com/science/article/B6VBN-3VM1X8N-3/2/8bffe94a1153eea8647eed5af0031d36.

//Bibtex Citation
//@article{Mandel1998,
//author = {Mandel, Denis and Abdollahzadeh, Ali and Maquin, Didier and Ragot, Jos�},
//isbn = {0301-7516},
//journal = {International Journal of Mineral Processing},
//keywords = {Linear Matrix Inequality Techniques,data reconciliation,error detection,error isolation},
//month = apr,
//number = {3},
//pages = {157--169},
//title = {{Data reconciliation by inequality balance equilibration: a LMI approach}},
//url = {http://www.sciencedirect.com/science/article/B6VBN-3VM1X8N-3/2/8bffe94a1153eea8647eed5af0031d36},
//volume = {53},
//year = {

// 12 Streams
// 5 Equipments 
// the measures
function [jac]=jacP9()
jac = [ 1  -1  -1   0   0   0   0   0    0   0   0   0   
        0   0   1   -1  -1  0   0    0   0   0   0   0   
        0   0   0   0   1   -1  -1  0    0   0   0   0 
        0   0   0   0   0   0   1    1    -1  0   0   0
        0   0   0   0   0   0   0   -1    0  1   0   -1
        0   0  0   0   0   0   0   0    1   -1  -1  0
        ];    
endfunction
function [jac]=jacP9_u_s10()
jac = [ 1  -1  -1   0   0   0   0   0    0   0   0   
        0   0   1   -1  -1  0   0    0   0   0   0   
        0   0   0   0   1   -1  -1  0    0   0   0 
        0   0   0   0   0   0   1    1   -1  0   0
        0   0   0   0   0   0   0   -1   1   -1  -1 ];    
endfunction

function [jac]=jacP9_u_s7()
jac = [ 1  -1  -1   0   0   0   0    0    0   0   0   
        0   0   1   -1  -1  0   0    0    0   0   0   
        0   0   0   0   1   -1  1    -1   0   0   0 
        0   0   0   0   0   0  -1    0    1   0   -1
        0   0  0   0   0   0   0     1   -1  -1   0];    
endfunction
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//Martins, Márcio A.F., Carolina A. Amaro, Leonardo S. Souza, Ricardo A. Kalid, and Asher Kiperstok. 2010. 
//New objective function for data reconciliation in water balance from industrial processes. 
//Journal of Cleaner Production (March): 1-6. doi:10.1016/j.jclepro.2010.03.014. http://linkinghub.elsevier.com/retrieve/pii/S0959652610001149.

//Bibtex Citation
//@article{Martins2010,
//author = {Martins, M\'{a}rcio A.F. and Amaro, Carolina A. and Souza, Leonardo S. and Kalid, Ricardo A. and Kiperstok, Asher},
//doi = {10.1016/j.jclepro.2010.03.014},
//file = {::},
//issn = {09596526},
//journal = {Journal of Cleaner Production},
//month = mar,
//pages = {1--6},
//title = {{New objective function for data reconciliation in water balance from industrial processes}},
//url = {http://linkinghub.elsevier.com/retrieve/pii/S0959652610001149},
//year = {2010}
//}

// 13 Streams
// 8 Equipments 
function [jac]=jacP10()
//      1   2   3   4   5   6   7   8    9   10  11  12  13
jac = [ 1  -1  -1  -1  -1   0   0   0    0   0   0   0   0
        0   1   0   0   0   0   0  -1    0   0   0   0   0
        0   0   1   0   0   0   0   0   -1   0   0   0   0
        0   0   0   1   0  -1  -1   0    0   0   0   0   0        
        0   0   0   0   1   0   0   0    0   0   1  -1   0
        0   0   0   0   0   1   0   0    0  -1   0   0   0
        0   0   0   0   0   0   1   0    0   0  -1   0   0
        0   0   0   0   0   0   0   1    1   1   0   0  -1];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13   
endfunction
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Rao, R Ramesh, and Shankar Narasimhan. 1996.
//“Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” 
//Industrial & Engineering Chemistry Research 35:1362-1368. 
//http://dx.doi.org/10.1021/ie940538b.
//Bibtex Citation

//@article{Rao1996,
//author = {Rao, R Ramesh and Narasimhan, Shankar},
//isbn = {0888-5885},
//journal = {Industrial \& Engineering Chemistry Research},
//month = apr,
//number = {4},
//pages = {1362--1368},
//publisher = {American Chemical Society},
//title = {{Comparison of Techniques for Data Reconciliation of Multicomponent Processes}},
//url = {http://dx.doi.org/10.1021/ie940538b},
//volume = {35},
//year = {1996}
//}

// 12 Streams
// 7 Equipments 
function [jac]=jacP11()
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12
jac = [ 1   -1  0   0   1   0   0   0    0   0   0   0 
        0   1   -1  0   0   0   -1  0    0   0   0   0 
        0   0   1   -1  0   -1  0   0    0   0   0   0 
        0   0   0   0   -1  1   0   1    0   0   0   0 
        0   0   0   0   0   0   0   -1   1   0   0  -1  
        0   0   0   0   0   0   1   0    -1  1   0   0
        0   0   0   0   0   0   0   0    0   -1  -1  1 ];                                
//      1   2   3   4   5   6   7   8    9   10  11  12
endfunction
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Mitsas, Christos L. 2010. Data reconciliation and variable classification by null space methods. 
//Measurement 43, no. 5 (June): 702-707.
// http://apps.isiknowledge.com/full_record.do?product=UA&search_mode=GeneralSearch&qid=2&SID=2A@bF9dN34I72L1Am9M&page=2&doc=17&colname=WOS.

//Bibtex Citation

//@article{Mitsas2010a,
//author = {Mitsas, Christos L.},
//journal = {Measurement},
//month = jun,
//number = {5},
//pages = {702--707},
//publisher = {ELSEVIER SCI LTD},
//title = {{Data reconciliation and variable classification by null space methods}},
//url = {http://apps.isiknowledge.com/full\_record.do?product=UA\&search\_mode=GeneralSearch\&qid=2\&SID=2A@bF9dN34I72L1Am9M\&page=2\&doc=17\&colname=WOS},
//volume = {43},
//year = {2010}
//}

// 12 Streams
// 6 Equipments 
function [jac]=jacP12()
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12
jac = [ 1  -1  -1  -1   0   0   0   0    0   0   0   0   
        0    1  0   0   -1  -1  -1  0    0   0   0   0   
        0    0  1   0   1   1   1   -1   0   0   0   0
        0    0  0   1   0   0   0   0    -1  -1  0   0
        0    0  0   0   0   0   0   1    1   0   -1  0
        0    0  0   0   0   0   0   0    0   1   1   -1
        ];  
//      1   2   3   4   5   6   7   8    9   10  11  12
endfunction
// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Fictitious but realistic mineral processing plant
//Alhaj-Dibo, Moustapha, Didier Maquin, and José Ragot. 2008.
//Data reconciliation: A robust approach using a contaminated distribution.
//Control Engineering Practice 16, no. 2 (February): 159-170.
// http://www.sciencedirect.com/science/article/B6V2H-4N4406D-1/2/50cac92b050f160a20a795faec990dc7.

//Bibtex Citation

//@article{Alhaj-Dibo2008,
//author = {Alhaj-Dibo, Moustapha and Maquin, Didier and Ragot, Jos\'{e}},
//isbn = {0967-0661},
//journal = {Control Engineering Practice},
//keywords = {Data reconciliation,Gross error detection,Linear and bilinear mass balances,Robust estimation},
//month = feb,
//number = {2},
//pages = {159--170},
//title = {{Data reconciliation: A robust approach using a contaminated distribution}},
//url = {http://www.sciencedirect.com/science/article/B6V2H-4N4406D-1/2/50cac92b050f160a20a795faec990dc7},
//volume = {16},
//year = {2008}
//}

// 16 Streams
// 9 Equipments 
// the measures
function [jac]=jacP13()
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    
jac = [ 1   -1  0   1   0   0   0   0    0   0   0   0   0     0     0      0    
        0   1   -1  0   0   0   0   0    0   0   -1   0   0     0     0      0    
        0   0   1   -1  -1  0   0   0    0   0   0   0   0     0     0      0    
        0   0   0   0   1   -1  0   0    0   1   0   0   0     0     0      0    
        0   0   0   0   0   1   -1  -1   0   0   0   0   0     0     0      0   
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    
        0   0   0   0   0   0   1   0    -1   -1  0   0   0     0     0      0    
        0   0   0   0   0   0   0   0    0   0   1   -1   -1    0     0      1    
        0   0   0   0   0   0   0   0    0   0   0   1   1     -1     0      0
        0   0   0   0   0   0   0   0    0   0   0   0   0     1     -1      -1];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    
endfunction

// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv

//Proposed by author

// 24 Streams
// 14 Equipments 

function [jac]=jacP14()
//The jacobian of the constraints
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    
jac = [ 1   1  1    1   -1  0   0   0    0   0   0   0   0     0     0      0    0      0    0    0     0    0     0     0    //M1
        0   0   0   0   1   -1  0   0    0   0   0   0   0     0     0      0    0      0    0    0     0    0     0     0    //F1
        0   0   0   0   0   1  -1   0    0   0   0   0   0     0     0      0    0      0    0    0     0    0     0     0    //T1
        0   0   0   0   0   0   1   -1   -1  0   0   0   0     0     0      0    0      0    0    0     0    0     0     0   //S1  
        0   0   0   0   0   0   0   0    1  -1   0   0   0     0     0      0    0      0    0    0     0    0     0     0     //F2
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    
        0   0   0   0   0   0   0   0    0   1   -1  0   0     0     -1     0    0      0    0    0     0    0     0     0    //M2
        0   0   0   0   0   0    0  0    0   0   1  -1  -1     0     0      0    0      0    0    0     0    0     0     0    //S3
        0   0   0   0   0   0   0   0    0   0   0   1   0     -1    0      0    0      0    0    0     0    0     0     0   //F3
        0   0   0   0   0   0   0   0    0   0   0   0   0     0     1      -1   0      0    -1   0     0    0     0     0    //F4
        0   0   0   0   0   0   0   0    0   0   0   0   0     0     0      1    -1     0    0    0     0    0     0     0    //F5
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    
        0   0   0   0   0   0   0    0   0   0   0   0   1     0     0      0    1     -1    0    0     0    0     0     0    //S5
        0   0   0   0   0   0   0    0   0   0   0   0   0     0     0      0    0      0    1    -1    0    0     0     0    //T2
        0   0   0   0   0   0   0    0   0   0   0   0   0     0     0      0    0      0    0    1     -1   -1    0     0    //S4
        0   0   0   0   0   0   0    1   0   0   0   0   0     0     0      0    0      0    0    0     0    0     -1    -1    //S2
        ];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24   
endfunction

// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Steam metering system
//Serth, R W, and W A Heenan. 1986.
// Gross error detection and data reconciliation in steam-metering systems. 
//AIChE Journal 32: 733-747.
//Bibtex Citation

//@article{Serth1986,
//author = {Serth, R W and Heenan, W A},
//journal = {AIChE Journal},
//pages = {733--747},
//title = {{Gross error detection and data reconciliation in steam-metering systems}},
//volume = {32},
//year = {1986}
//}

// 28 Streams
// 11 Equipments 
// the measures
function [jac]=jacP15()
//The jacobian of the constraints
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    25    26    27    28
jac = [ 1   1   -1  1   0   0   0   0    0   0   0   0   0     0     0      0    0      0    0    0     0    0     0     0     0     0    0      0
        0   0   0   0  -1   -1  1   1    -1  0   0   0   0     0     0      0    0      0    0    0     0    0     0     0     0     0    0      0 
        -1  0   0   0   1    0   0   0   0  -1   0   0   0     0     0      0    0      0    0    0     0    0     0     0     0     0    0      0
        0   0   0   0   0   0   0   0    0   1   1   -1  0     0     0      0    0      0    0    0     0    0     0     0     0     0    0      0
        0   0   1   0   0   0   0   0    0   0   -1  0   1     -1    -1     -1   -1     0    0    0     0    0     0     0     0     0    0      0  
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    25    26    27    28
        0   -1  0   0   0   1   0   0    0   0   0   0   -1     0     0     0    0      0    0    0     0    0     0     0     0     0    0      0
        0   0   0   0   0   0   -1  0    0   0   0   0   0     1     0      0    0      1    -1   -1    -1   0     0     0     0     0    0      0
        0   0   0   0   0   0   0   0    0   0   0   0   0     0     1      0    0      -1   0    0     0    1     -1    -1    0     0    0      0
        0   0   0   0   0   0   0   0    0   0   0   1   0     0     0      1    0      0    0    0     0    -1    0     0     -1    0    0      0  
        0   0   0   0   0   0   0   0    0   0   0   0   0     0     0      0    0      0    1    0     0    0     1     0     0     -1    1     0
        0   0   0   0   0   0   0   -1   0   0   0   0   0     0     0      0    0      0    0    1     0    0     0     0     0     1     0     1
        ];                                
//      1   2   3   4   5   6   7   8    9   10  11  12  13    14    15    16    17    18   19   20    21   22    23    24    25    26    27    28

endfunction


// Data Reconciliation Benchmark Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
// Atmospheric tower example from:
// Zhang, P, G Rong, and Y Wang. 2001.
// A new method of redundancy analysis in data reconciliation and its application.
// Computers and Chemical Engineering 25: 941-949.

//Bibtex Citation

//@article{Zhang2001,
//author = {Zhang, P and Rong, G and Wang, Y},
//journal = {Computers and Chemical Engineering},
//pages = {941--949},
//title = {{A new method of redundancy analysis in data reconciliation and its application}},
//volume = {25},
//year = {2001}
//}
// 50 Streams
// 26 Equipments 

function [jac]=jacP16()
//The jacobian of the constraints

//S319	S316	S312	S378	S336	S357	S346	S359P1	S347	S352	S356	S358	S357P	S359P2	S359	S338P	S338	S341P	S341	S414	S502	S411	S401	S415	S402	S404	S405	S407	S408	S453	S460	S456	S452	S511	S503	S384P	S52P	S592	S581	S525	S524	S536	S527	S549	S550	S537	S598	S599	S267	S538

jac = [1	1	1	1	-1	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	-1	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	1	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	-1	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	-1	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	1	-1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	-1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	-1	-1	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	-1	-1	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	-1	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	-1	-1	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	-1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	-1];

endfunction


function [jac]=jacP16_P1()
// P16 part 1
//The jacobian of the constraints
// P16 part 1 // firs 35 streams plus the  lats stream (S538) , only 19 first equipment
//S319	S316	S312	S378	S336	S357	S346	S359P1	S347	S352	
//S356	S358	S357P	S359P2	S359	S338P	S338	S341P	S341	S414	
//S502	S411	S401	S415	S402	S404	S405	S407	S408	S453	
//S460	S456	S452	S511	S503	S538

jac = [
1	1	1	1	-1	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 // DA301
0	0	0	0	0	0	1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 // Split 1
0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 // Comp1
0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 // EA 32X
0	0	0	0	0	0	0	0	0	0	1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 // FA 309
0	0	0	0	0	-1	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 // Split 2
0	0	0	0	0	0	0	1	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 // Mixer 1
0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	-1	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 // FSPL
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 // H 338
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 // H 341
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	-1	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0 // DA 401
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 // PTC
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	-1	0	0	0	0	0	0	0	0	0	0	0 // M1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0 // 448-450
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0 // DC 401
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0 // 552-448
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	1	-1	1	0	0	0	0	0	0 // DA 408
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	-1	-1	-1	0	0	0 // DA 402
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	-1	-1	1]; // DA 404


endfunction 
 
 
 function [jac]=jacP16_P2()
// P16 part 2 // last 17 streams , only last 7 equipment
//	S511	S503	S384P	S52P	S592	S581	S525	
//S524	S536	S527	S549	S550	S537	S598	S599	S267	S538

jac = [	0	1	0	0	-1	-1	0	0	0	0	0	0	0	0	0	0	0 // DA 405
	1	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0 // DC 402
	0	0	0	1	0	0	-1	-1	0	0	0	0	0	0	0	0	0 // DC 40x
	0	0	0	0	0	0	0	1	-1	-1	0	0	0	0	0	0	0 // DA 407
	0	0	0	0	0	0	0	0	0	1	-1	-1	-1	0	0	0	0 // DA 406
	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	-1	0	0 // DA 409
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	-1]; // MIX 2

endfunction

 
