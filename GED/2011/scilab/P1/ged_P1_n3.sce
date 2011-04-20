// Data Reconciliation Benchmark and GED Problems From Lietrature Review
// Author: Edson Cordeiro do Valle
// Contact - edsoncv@{gmail.com}{vrtech.com.br}
// Skype: edson.cv
getd('.');
clear rerror rerror1 xr xrs sd sds x_sol f_sol status grerror grerrors mySign xfinal jac jac_col jac_col rj sigma sigam_inv res opt_error V V_inv diag_diag_V Wbar gama zr_nt adj zadj d_adj zadj_alt Wbar_alt x_chi ge_gt nge_gt ge_nt1_i ge_nt1_j ge_nt2_i ge_nt2_j norm_nt1 norm_nt2 norm_mt norm_mt2 ge_mt_i ge_mt_j ge_mt_alt_i ge_mt1_alt_i ge_mt1_alt_j ge_mt2_alt_i ge_mt2_alt_j ge_mt_alt2_j ge_mt1_i ge_mt1_j ge_mt2_i ge_mt2_j Vinv_r zr_nt_alt ge_nt_alt_i ge_nt_alt_j op_nt1 op_nt2 op_mt1 op_mt2 op_mt1_alt op_mt2_alt op_nt_alt_1 op_nt_alt_2 adjustability detect;
xr=[8.5;4.5;4];
szx = size(xr,1);
runsize = 100;
sd = 0.0001*xr;
sds = 0.0001*xr;
//rerror=grand(runsize,szx,'nor',0,1);
// random number generators: rerror1 prefered as rerror
for i=1:szx
    rerror1(:,i)=grand(runsize,1,'nor',0,sd(i));
end

// adding random error to exact x
for i=1:runsize
//    xrs(:,i) = xr + sd.*rerror(i,:)';
    xrs(:,i) = xr + rerror1(i,:)';
end

// random sign generator to add to gross error
mySign=sign(grand(runsize,szx,'unf',-1,1));
k=0;
for i=1:szx
    for j=1:runsize
        grerror(j+k*runsize,i) = grand(1,1,'unf',xrs(i,j)*0.05,xrs(i,j)*0.1)*mySign(j,i);
    end
    k=k+1;
end

for i=0:(size(grerror,1)-1)
    grerrors(i+1,1:szx) = grerror(i+1,1:szx) + xrs(:,modulo(i,runsize)+1)';
end
xfinal =[xrs';grerrors]; 
// dataset generation ends here
tic;
runsizefinal = size(xfinal,1);
for i=1:runsizefinal
//    xrs(:,i) = xr + sd.*rerror(i,:)';
    [x_sol(i,:), f_sol(i), status(i)] = P1(xfinal(i,:)',sds.^2);
end

toc
opt_error=find(status<>0);
// Global test
jac=jacP1();
rj=rank(jac);
jac_col = size(jac,2);
jac_row = size(jac,1);
sigma=diag(sds.^2);
sigma_inv=inv(sigma);
// adjustability and detectability
[adjustability, detect] = adjust(sigma, jac);
// variance-covariance matrix: narasimham pg. 178 eq. 7-3
V=jac*sigma*jac';
V_inv= inv(V);
diag_diag_V = diag(diag(V));

// covariance matrix of adjustments: narasimham pg. 183 eq. 7-13
Wbar=sigma*jac'*inv(V)*jac*sigma;
for i=1:runsizefinal
// residuals: narasimham pg. 178 eq. 7-2    
    res(i,1:jac_row)=(jac*xfinal(i,:)')';
// global test statistics: narasimham pg. 178 eq. 7-4       
    gama(i) = res(i,:)*V_inv*res(i,:)';
// nodal test statistics: narasimham pg. 180 eq. 7-5
    for k=1:jac_row

        zr_nt(i,k)=abs(res(i,k))./(diag_diag_V(k,k).^0.5);
        
    end
    for j=1:szx
//adjustments  narasimham pg. 183 eq. 7-11       
        adj(i,j) = xfinal(i,j)-x_sol(i,j);
// measurements test statistics: narasimham pg. 183 eq. 7-14        
        zadj(i,j)=abs(adj(i,j))/sqrt(Wbar(j,j));
    end
// alternative measurements test statistics: narasimham pg. 183 eq. 7-15            
    d_adj(:,i)=sigma_inv*adj(i,:)';
end
// Alternative residual statistics narasimham pg. 183 eq. 7-9


for i=1:runsizefinal
    Vinv_r(i,1:jac_row) = (inv(V)*res(i,:)')';;
    for j=1:jac_row
        zr_nt_alt(i,k)=abs(Vinv_r(i,j))./(diag_diag_V(k,k).^0.5);       
    end
end

// covariance matrix of adjustments, alternative formulation: narasimham pg. 183 eq. 7-16
Wbar_alt=jac'*inv(jac*sigma*jac')*jac;
for i=1:runsizefinal
    for j=1:szx
// alternative measurements test statistics: narasimham pg. 183 eq. 7-17        
        zadj_alt(i,j) = d_adj(j,i)/sqrt(Wbar_alt(j,j));
    end
end

Q=0.05;
P=1-Q;
xchi=cdfchi("X",rj,P,Q);
ge_gt=find(gama>xchi);
ge_gt_low=length(find(gama(1:runsize)>xchi));
nge_gt=length(ge_gt)
//find elements wrongly idendified
//ge_gt_remove_index =[];
//for i=0:jac_col-1
//    ge_gt_tmp = (i*runsize+ge_gt_low+(ge_gt(find(ge_gt <=runsize))))';
//    if ge_gt_tmp == 0 then         
//             break
//             else ge_gt_remove_index = [ge_gt_remove_index;ge_gt_tmp];
//    end,
//    
//end
//
//ge_gt_remove_index_all = [(find(ge_gt <=runsize))';ge_gt_remove_index];
//ge_gt_new_removed = ge_gt;
//ge_gt_new_removed(ge_gt_remove_index_all) = [];

// Overall Power
op_gt = length(find(ge_gt>runsize))./(runsize*jac_col);


//Nodal test
beta_r = (1-((1-Q).^(1/jac_row)));

Q1=Q/2;
Q2=beta_r/2;

P1=1-Q1;
P2=1-Q2;
norm_nt1=cdfnor("X",0,1,P1,Q1);
norm_nt2=cdfnor("X",0,1,P2,Q2);

[ge_nt1_i,ge_nt1_j]=find(zr_nt>norm_nt1);
[ge_nt2_i,ge_nt2_j]=find(zr_nt>norm_nt2);

[ge_nt1_alt_i,ge_nt1_alt_j]=find(zr_nt_alt>norm_nt1);
[ge_nt2_alt_i,ge_nt2_alt_j]=find(zr_nt>norm_nt2);

op_nt1 = length(find(ge_nt1_i>runsize))./(runsize*jac_col);
op_nt2 = length(find(ge_nt2_i>runsize))./(runsize*jac_col);
op_nt_alt_1 = length(find(ge_nt1_alt_i>runsize))./(runsize*jac_col);
op_nt_alt_2 = length(find(ge_nt2_alt_i>runsize))./(runsize*jac_col);


//Measurement test
norm_mt=cdfnor("X",0,1,P1,Q1);

beta_m = (1-((1-Q).^(1/jac_col)));
Q2=beta_m/2;

P2=1-Q2;

norm_mt2=cdfnor("X",0,1,P2,Q2);

[ge_mt1_i,ge_mt1_j]=find(zadj>norm_mt);
[ge_mt2_i,ge_mt1_j]=find(zadj>norm_mt2);


[ge_mt1_alt_i,ge_mt1_alt_j]=find(zadj_alt>norm_mt);
[ge_mt2_alt_i,ge_mt2_alt_j]=find(zadj_alt>norm_mt2);

op_mt1 = length(find(ge_mt1_i>runsize))./(runsize*jac_col);
op_mt2 = length(find(ge_mt2_i>runsize))./(runsize*jac_col);
op_mt1_alt = length(find(ge_mt1_alt_i>runsize))./(runsize*jac_col);
op_mt2_alt = length(find(ge_mt2_alt_i>runsize))./(runsize*jac_col);
//measurement
//for i=1:runsize
//    for j=1:size(xr,1)
//        adj(i,j) = xrs(j,j)-x_sol(i,j);
//    end
//end


//nge_nt=length(ge_nt)


