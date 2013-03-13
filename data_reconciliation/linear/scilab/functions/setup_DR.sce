mode(-1);
select obj_function_type 
	case 0 then
	getd('../functions/wls');
 
	case 1 then
		getd('../functions/absolute');
		getd('../functions');
		beta_smooth = 0.25;
		alpha_smooth = 1/beta_smooth;
	case 2 then
		getd('../functions/cauchy');
		getd('../functions');
		beta_smooth = 0.25;
		alpha_smooth = 1/beta_smooth;
		// cauchy parameter
		const_cauchy = 2.3849;
	case 3 then
		getd('../functions/contamined_normal');
		getd('../functions');
		beta_smooth = 0.25;
		alpha_smooth = 1/beta_smooth;
		// contamined_normal parameters
		const1_cont_nor = 0.235;
		const2_cont_nor = 10;

	case 4 then
		getd('../functions/fair');
		getd('../functions');
		beta_smooth = 0.25;
		alpha_smooth = 1/beta_smooth;
		// fair function parameter
		const_fair = 1.3998;


	case 5 then
		getd('../functions/hampel');
		getd('../functions');
		beta_smooth = 0.25;
		alpha_smooth = 1/beta_smooth;
		// Hampel parameters
		a=0.5;
		b=1;
		c=2;
		//a=1;
		//b=2;
		//c=4;
		//a=2;
		//b=4;
		//c=8;
		//a=4;
		//b=8;
		//c=16;
		ones_a = a*ones(length(red),1);
		ones_b = b*ones(length(red),1);
		ones_c = c*ones(length(red),1);
		ones_xm = ones(length(red),1);


	case 6 then
		getd('../functions/logistic');
		getd('../functions');
		beta_smooth = 0.25;
		alpha_smooth = 1/beta_smooth;
		// logistic parameter
//		const_logist = 0.602;
		const_logist = 0.10;

	case 7 then
		getd('../functions/lorenztian');
		getd('../functions');
		beta_smooth = 0.25;
		alpha_smooth = 1/beta_smooth;
		// Lorentzian parameter
		const_lor = 2.6;

	case 8 then
		getd('../functions/quasi_weighted');
		getd('../functions');
		beta_smooth = 0.25;
		alpha_smooth = 1/beta_smooth;
		// quasi-weighted parameter
//		const_qw = 0.15;
        ones_xm = ones(length(red),1);
		const_qw = 0.1;

end
