Data Reconciliation Benchmark Problems From Lietrature Review Release 2011
Author: Edson Cordeiro do Valle
Contact - edsoncv@{gmail.com}{vrtech.com.br}
Skype: edson.cv


In this type of problem, we deal with linear data reconciliation (total flow measurements), with all streams measured and known variance/standard deviation data. Although linear problems are easy to solve in DR, since it have analytical solution, the aim of this problem is to test different objective function types, called robust functions. These robust functions can provide results with gross errors removed, however, to do this, some parameters must be tuned in these function. The aim of this class problems is to test these parameters set up. 

In this problem, since we use robust functions, and the derivatives are complicate to express, we use Ipopt optimizer to run the optimizations. We also use automatic differentiation toolbox from Scilab (diffcode) to evaluate the Hessian and Jacobian structure of objective function.

inside "linear" folder, we have 3 sobfolders:

-diagrams: a diagram of the process in the .dia and .png format.
-scilab: the scilab implementation.
-spreadsheet: a Excel spreadsheet of the data, results or some data handling necessary to run the model. 


Reference Index:

P1	Rao, R Ramesh, and Shankar Narasimhan. 1996. “Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” Industrial & Engineering Chemistry Research 35:1362-1368. http://dx.doi.org/10.1021/ie940538b.

P2	Dovì, V G, and C Solisio. 2001. “Reconciliation of censored measurements in chemical processes: an alternative approach.” Chemical Engineering Journal 84:309-314. http://www.sciencedirect.com/science/article/B6TFJ-45KNHB1-F/2/199f358469628f600f10b394d2b55a8b.

P3	Zhang, Zhengjiang, Zhijiang Shao, Xi Chen, Kexin Wang, and Jixin Qian. 2010. “Quasi-weighted least squares estimator for data reconciliation.” Computers & Chemical Engineering 34:154-162. http://www.sciencedirect.com/science/article/B6TFT-4XDCHNS-1/2/63a2b79a4cc89a3afb57ff83c4063242.

P4	Narasimhan, S, and C Jordache. 2000. Data Reconciliation and Gross Error Detection: An Intelligent Use of Process Data. 1st ed. Houston: Gulf Publishing.

P5	Yang, Youqi, Rongbo Ten, and Luiqun Jao. 1995. “A study of gross error detection and data reconciliation in process industries.” Comp. & Chem. Eng 19S:S217-S222.

P6	﻿Rosenberg, J, R S H Mah, and C Iordache. 1987. “Evaluation of Schemes for Detecting and Identifying Gross Errors in Process Data.” Ind. & Eng. Chem. Proc. Des. Dev, Vol. 26: 555-564.

P7	Proposed

P8	Rao, R Ramesh, and Shankar Narasimhan. 1996. “Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” Industrial & Engineering Chemistry Research 35:1362-1368. http://dx.doi.org/10.1021/ie940538b.

P9	Mandel, Denis, Ali Abdollahzadeh, Didier Maquin, and José Ragot. 1998. “Data reconciliation by inequality balance equilibration: a LMI approach.” International Journal of Mineral Processing 53:157-169. http://www.sciencedirect.com/science/article/B6VBN-3VM1X8N-3/2/8bffe94a1153eea8647eed5af0031d36.

P10	Martins, Márcio A.F., Carolina A. Amaro, Leonardo S. Souza, Ricardo A. Kalid, and Asher Kiperstok. 2010. “New objective function for data reconciliation in water balance from industrial processes.” Journal of Cleaner Production 1-6. http://linkinghub.elsevier.com/retrieve/pii/S0959652610001149.

P11	Rao, R Ramesh, and Shankar Narasimhan. 1996. “Comparison of Techniques for Data Reconciliation of Multicomponent Processes.” Industrial & Engineering Chemistry Research 35:1362-1368. http://dx.doi.org/10.1021/ie940538b.

P12	Mitsas, Christos L. 2010. “Data reconciliation and variable classification by null space methods.” Measurement 43:702-707. http://apps.isiknowledge.com/full_record.do?product=UA&search_mode=GeneralSearch&qid=2&SID=2A@bF9dN34I72L1Am9M&page=2&doc=17&colname=WOS (Accessed July 22, 2010).

P13	Alhaj-Dibo, Moustapha, Didier Maquin, and José Ragot. 2008. “Data reconciliation: A robust approach using a contaminated distribution.” Control Engineering Practice 16:159-170. http://www.sciencedirect.com/science/article/B6V2H-4N4406D-1/2/50cac92b050f160a20a795faec990dc7.

P14	Proposed

P15	Serth, R W, and W A Heenan. 1986. “Gross error detection and data reconciliation in steam-metering systems.” AIChE Journal 32:733-747.

P16	Zhang, Zhengjiang, Zhijiang Shao, Xi Chen, Kexin Wang, and Jixin Qian. 2010. “Quasi-weighted least squares estimator for data reconciliation.” Computers & Chemical Engineering 34:154-162. http://www.sciencedirect.com/science/article/B6TFT-4XDCHNS-1/2/63a2b79a4cc89a3afb57ff83c4063242.
