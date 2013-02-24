function res = flowsheet_residuals(x, K_coef, cp1_coef, h_hx_coef, frac)
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
// y:        the contraints residuals
// Inputs:
// x:         the column vector of the variables: x = [flow, temperatures]
// all the equations are coded "by hand":

xc = matrix(x(1:$-3), 11,8)';
Q_heater1 = x($-2);
Q_heater2 = x($-1);
advance = x($);
//eq_const = x($);

press = xc(1,:);
flows = xc(2,:);
temper = xc(3,:);
comps = xc(4:$,:);
res = zeros(83,1);

A1 = cp1_coef(1,:);
B1 = cp1_coef(2,:);
C1 = cp1_coef(3,:);
D1 = cp1_coef(4,:);

tt = (temper(3) + 273.15 + temper(4) + 273.15)/2

Cp_t = 8.314*(A1 + B1.*tt + C1.*(tt).^2 + D1.* (tt).^-2); // (in J/(mol.K)
//pause
//  Mixer (M)
// - mixing heat = 0 (in this consideration we are assuming an error of ~ 0.2% )
// mass compound balance (5 equations)
res(1:5) = flows(10)*comps(:,10) + flows(2)*comps(:,2) - flows(3)*comps(:,3);
// energy balance need a scaling factor here
//res(6)  = ((flows(10)*temper(10) + flows(2)*temper(2)) -flows(3)*temper(3))/10000;
res(6)  = ((flows(10)*temper(10) + flows(2)*temper(2)) -flows(3)*temper(3));
// pressures
res(7)  = press(10) - press(3);
res(8)  = press(2) - press(3);
// Heater (H)
// mass balance
res(9:13) = flows(3)*comps(:,3) - flows(4)*comps(:,4); 
// energy balance this may result in a large error, but the relative error is small, we use a scalling factor here
res(14) = (Q_heater1 - flows(3)*(comps(1,3)*Cp_t(1) + comps(2,3)*Cp_t(2) + comps(3,3)*Cp_t(3) + comps(4,3)*Cp_t(4) + comps(5,3)*Cp_t(5) )*(temper(4)-temper(3)))/10000;
// pressure -  pressure drop of 1 bar in heater H1 // this will be accounted in constraints limits
res(15) = press(3) - press(4);
// Reactor (R)
// pressure are fixed (controled)
// equilibrium constant
k_eq = (10^(-2.691122*log10(temper(5) + 273.15) - (temper(5) +273.15)*5.519265e-5 +((temper(5) + 273.15)^2)*1.848863e-7 +(2001.6/(temper(5)+273.15) + 2.6899)))^2;
// conversion
//res(16) = comps($,5) - (comps(1,5)^(0.5)*comps(2,5)^(1.5))*press(5)*k_eq;
res(16) = comps($,5)^2 - ((comps(1,5)^(3)*comps(2,5))*press(5)^2)*k_eq;
// Ammonia
// no: total moles in the reactor inlet
//no =flows(4)*(comps(1,4) + comps(2,4) + comps($,4));
//res(17) = (flows(4) - advance)*comps($,5) - (flows(4)*(comps($,4)) + advance); 
res(17) = (flows(4) - 2*advance)*comps($,5) - (flows(4)*(comps($,4)) + 2*advance); 
// Hydrogen
res(18) = (flows(4) - 2*advance)*comps(1,5) - (flows(4)*(comps(1,4)) - 3*advance);
// Nitrogen
res(19) = (flows(4) - 2*advance)*comps(2,5) - (flows(4)*(comps(2,4)) - 1*advance);
// Inherts (CH4 and Ar)
res(20) = flows(4)*comps(3,4) - flows(5)*comps(3,5);
res(21) = flows(4)*comps(4,4) - flows(5)*comps(4,5);

// Heater 2 (H 2)
// Mass balance
res(22:26) = flows(5)*comps(:,5) - flows(6)*comps(:,6);
// Energy Balance

h1 = h_hx_coef(1,1)*temper(5) + h_hx_coef(1,2);        // adjusted from simulation data
h2 = h_hx_coef(2,1)*temper(6) + h_hx_coef(2,2);     // adjusted from simulation data
// we need to add an scalling factor here
res(27) = (Q_heater2 - flows(5)*(-h2+h1))/10000;
// Pressure -  pressure drop of 1 bar in heater H2 // this will be accounted in constraints limits
res(28) = press(6) - press(5);
// Flash (F)
// Ki calculation
res(29:33) = comps(:,11).*(K_coef(:,1)*temper(6) + K_coef(:,2)) - comps(:,7); 
//res(14) = comps(2,9)*(coef(1,2)*temper(4) + coef(2,2)) - comps(2,5);
//res(15) = comps(3,9)*(coef(1,3)*temper(4) + coef(2,3)) - comps(3,5);
//res(16) = comps(4,9)*(coef(1,4)*temper(4) + coef(2,4)) - comps(4,5);
//res(17) = comps(5,9)*(coef(1,5)*temper(4) + coef(2,5)) - comps(5,5);
// mass balance
res(34:38) = flows(6)*comps(:,6) - flows(7)*comps(:,7) - flows(11)*comps(:,11);
// sum xi = 1 must be set in the upper and lower bounds
res(39) = sum(comps(:,7)) ;
res(40) = sum(comps(:,11)) ;
// energy and pressure
res(41) = temper(6) -  temper(7) ;
res(42) = temper(6) -  temper(11) ;

res(43) = press(6) -  press(7) ;
res(44) = press(6) -  press(11)  ;
// Splitter (S)
// mass balance
res(45:49) = flows(8)*comps(:,8) - flows(7)*comps(:,7) * frac ;
res(50:54) = flows(9)*comps(:,9) - flows(7)*comps(:,7) * (1- frac) ;
// energy and pressure
res(55) = temper(7) -  temper(8) ;
res(56) = temper(7) -  temper(9) ;

res(57) = press(7) -  press(8) ;
res(58) = press(7) -  press(9) ; 

// Compressor 2 (C2)
// mass balance
res(59:63) = flows(9)*comps(:,9)  - flows(10)*comps(:,10);

// energy and pressure // need to add to independent terms
//(-69.27) and -11 respectivelly
res(64) = temper(10) - temper(9) - ( 0.348713*press(10));
//res(64) = temper(10) - ( 0.348713*press(10));
res(65) = press(9) -  press(10)  ;

// Compressor 1 (C1)
// mass balance
res(66:70) = flows(1)*comps(:,1)  - flows(2)*comps(:,2);

// energy and pressure // need to add to independent terms (222) and -200 respectivelly
res(71) = temper(2) - temper(1) - ( 0.975429*press(2));
//res(71) = temper(2) - ( 0.975429*press(2));
res(72) = press(1) -  press(2)  ;

// The last equations helps to add redundancy to the system
// for stream 7 and 11  the equations were already added in the flash
res(73) = sum(comps(:,1)) ;
res(74) = sum(comps(:,2)) ;
res(75) = sum(comps(:,3)) ;
res(76) = sum(comps(:,4)) ;
res(77) = sum(comps(:,5)) ;
res(78) = sum(comps(:,6)) ;
res(79) = sum(comps(:,8)) ;
res(80) = sum(comps(:,9)) ;
res(81) = sum(comps(:,10)) ;
res(82) = temper(4) - temper(5);
res(83) = press(4) - press(5);
//res(84) = temper(10) - temper(9);

endfunction
