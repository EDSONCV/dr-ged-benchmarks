getd('.');
clear rerror rerror1 xr xrs sd sds x_sol f_sol status jac jac_col jac_col rj sigma sigam_inv res opt_error V V_inv diag_diag_V Wbar gama zr_nt adj zadj d_adj zadj_alt Wbar_alt x_chi ge_gt nge_gt ge_nt_i ge_nt_j ge_mt_i ge_mt_j ge_mt_alt_i ge_mt_alt_j;
//xr=[100;99;1;99;1;100];
xr=[101.91;68.45;34.65;64.2;36.44;98.88];
szx = size(xr,1);
runsize=1;
sds = ones(6,1);
//sds = (0.01*xr).^2;
x_sol = zeros(runsize,szx);
f_sol=zeros(runsize);
status = -20*ones(runsize)
rerror=grand(runsize,szx,'nor',0,1);

tic;
for i=1:runsize
//    xrs(:,i) = xr + sd.*rerror(i,:)';
    [x_sol(i,:), f_sol(i), status(i)] = P4(xr,sds);
end
toc
opt_error=find(status<>0);
// Global test
jac=jacP4();
rj=rank(jac);
jac_col = size(jac,2);
jac_row = size(jac,1);
sigma=diag(sds.^2);
sigma_inv=inv(sigma);
// variance-covariance matrix: narasimham pg. 178 eq. 7-3
V=jac*sigma*jac';
V_inv= inv(V);
diag_diag_V = diag(diag(V));
// covariance matrix of adjustments: narasimham pg. 183 eq. 7-13
Wbar=sigma*jac'*inv(V)*jac*sigma;
for i=1:runsize
// residuals: narasimham pg. 178 eq. 7-2    
    res(i,1:jac_row)=(jac*xr(:,i))';
// global test statistics: narasimham pg. 178 eq. 7-4       
    gama(i) = res(i,:)*V_inv*res(i,:)';
// nodal test statistics: narasimham pg. 180 eq. 7-5
    for k=1:jac_row

        zr_nt(i,k)=abs(res(i,k))./(diag_diag_V(k,k).^0.5);
    
    end
    for j=1:szx
//adjustments  narasimham pg. 183 eq. 7-11       
        adj(i,j) = xr(j,i)-x_sol(i,j);
// measurements test statistics: narasimham pg. 183 eq. 7-14        
        zadj(i,j)=abs(adj(i,j))/sqrt(Wbar(j,j));
    end
// alternative measurements test statistics: narasimham pg. 183 eq. 7-15            
    d_adj(:,i)=sigma_inv*adj(i,:)';
end
// covariance matrix of adjustments, alternative formulation: narasimham pg. 183 eq. 7-16
Wbar_alt=jac'*inv(jac*sigma*jac')*jac;
for i=1:runsize
    for j=1:szx

// alternative measurements test statistics: narasimham pg. 183 eq. 7-17        
        zadj_alt(i,j) = d_adj(j,i)/sqrt(Wbar_alt(j,j));
    end
end


Q=0.05;
P=1-Q;
xchi=cdfchi("X",rj,P,Q);
ge_gt=find(gama>xchi)
nge_gt=length(ge_gt)
Q=Q/2;
P=1-Q;
//Nodal test
norm_nt=cdfnor("X",0,1,P,Q);
[ge_nt_i,ge_nt_j]=find(zr_nt>norm_nt);
//Measurement test
[ge_mt_i,ge_mt_j]=find(zadj>norm_nt);
beta_m = (1-((1-Q).^(1/3)));
Q=beta_m/2
P=1-Q;
norm_mt_alt=cdfnor("X",0,1,P,Q);
[ge_mt_alt_i,ge_mt_alt_j]=find(zadj_alt>norm_mt_alt);
//measurement
//for i=1:runsize
//    for j=1:size(xr,1)
//        adj(i,j) = xrs(j,j)-x_sol(i,j);
//    end
//end


//nge_nt=length(ge_nt)


