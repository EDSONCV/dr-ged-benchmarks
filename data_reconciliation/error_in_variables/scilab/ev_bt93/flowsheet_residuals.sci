function y = flowsheet_residuals(x)
//*********************************************************************        
// Data Reconciliation Benchmark Problems From Literature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
//*********************************************************************
// this function is prepared to use the automatic derivatives toolbox of 
// Scilab. This toolbox can be instaled using the ATOMS installer (package name: diffcode).
// This function evaluates the residuals of the flowsheet of the problem proposed by Swartz, 1989
//
// Outputs:
// y,:        the contraints residuals. The number of constraints depends on the number of the 
//            data sets choosen (see 'bt93.sec' file)
// Inputs:
// x:         the column vector of the variables: x - it is reorganized inside this function
// The resulting system is according to the papper from Biegler and Tjoa (1993):

xdata = matrix(x(1:$-4),21, ndata);
//pause
FA6=xdata(1,:);
TA2=xdata(2,:);
TA4=xdata(3,:);
TA5=xdata(4,:);
TA7=xdata(5,:);
TA8=xdata(6,:);
TD1=xdata(7,:);
TC1=xdata(8,:);
TB1=xdata(9,:);
TB2=xdata(10,:);
FA1=xdata(11,:);
FA3=xdata(12,:);
FC1=xdata(13,:);
FD1=xdata(14,:);
FB1=xdata(15,:);
TA1=xdata(16,:);
TB3=xdata(17,:);
TC2=xdata(18,:);
TD2=xdata(19,:);
TA3=xdata(20,:);
TA6=xdata(21,:);
UAint = x($-3:$);
dt1_1 = (TB2 - TA2);
dt1_2 = (TB3 - TA1);
dt2_1 = (TB1 - TA4);
dt2_2 = (TB2 - TA3);
dt3_1 = (TC1 - TA5);
dt3_2 = (TC2 - TA4);
dt4_1 = (TD1 - TA7);
dt4_2 = (TD2 - TA6);
//mldt1 =(dt1_1.*dt1_2.*((dt1_1+dt1_2)/2)).^(1/3);
//mldt2 =(dt2_1.*dt2_2.*((dt2_1+dt2_2)/2)).^(1/3);
//mldt3 =(dt3_1.*dt3_2.*((dt3_1+dt3_2)/2)).^(1/3);
//mldt4 =(dt4_1.*dt4_2.*((dt4_1+dt4_2)/2)).^(1/3);
//
mldt1 = ((TB2 - TA2) - (TB3 - TA1))./log((TB2 - TA2)./(TB3 - TA1)); 
mldt2 = ((TB1 - TA4) - (TB2 - TA3))./log((TB1 - TA4)./(TB2 - TA3));
mldt3 = ((TC1 - TA5) - (TC2 - TA4))./log((TC1 - TA5)./(TC2 - TA4));
mldt4 = ((TD1 - TA7) - (TD2 - TA6))./log((TD1 - TA7)./(TD2 - TA6));
//pause
// mass balance
//y = zeros (1:14*ndata);
y(1:ndata) = FA1 - FA3 - FA6;
// energy balance
y(ndata+1: 2*ndata) = FA1.*(TA2 - TA1) - FB1.*(TB2 - TB3) ;
y(2*ndata+1: 3*ndata) = FA3.*(TA4 - TA3) - FB1.*(TB1 - TB2);
y(3*ndata+1: 4*ndata) = FA3.*(TA5 - TA4) - FC1.*(TC1 - TC2) ;
y(4*ndata+1: 5*ndata) = FA6.*(TA7 - TA6) - FD1.*(TD1 - TD2);
y(5*ndata+1: 6*ndata) = FA1.*TA8 - FA3.*TA5 - FA6.*TA7;
// equations related with parameters
y(6*ndata+1: 7*ndata) = UAint(1).*mldt1 - FB1.*(TB2 - TB3); 
y(7*ndata+1: 8*ndata) = UAint(2).*mldt2 - FB1.*(TB1 - TB2);
y(8*ndata+1: 9*ndata) = UAint(3).*mldt3 - FC1.*(TC1 - TC2);
y(9*ndata+1: 10*ndata) =UAint(4).*mldt4 - FD1.*(TD1 - TD2); 

//equations relating unmeasured streams
//inequality constraints
y(10*ndata+1: 11*ndata) = TB2 - TA2;
y(11*ndata+1: 12*ndata) = TB1 - TA4;
y(12*ndata+1: 13*ndata) = TC1 - TA5;
y(13*ndata+1: 14*ndata) = TD1 - TA7;
//y(14*ndata+1: 15*ndata) = TD2 - TA6;
//y(15*ndata+1: 16*ndata) = TB2 - TA3;
y=y';
//pause

endfunction
