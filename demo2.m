clc,clear

load("UMIST.mat")
Y=X;
M=[20 20 20];
R=15;
ROU=[1e1 1e3 1e0];
mu=1e1;
maxiter=5000;
num_center=10;
% % % 
lam=0.2;

[ES_Y,iter_time,re_list,total_time,a,idx,C] = JTKM_THALS(Y,M,R,ROU,mu,lam,num_center,maxiter);


% sir=SIR(Y,ES_Y);
% r_e=norm(Y(:)-ES_Y(:))/norm(Y(:));
% psnr=PSNR(Y,ES_Y);

sub=10;
n=10;


ac=zeros(1,sub);
% 
% 
for i=1:sub
    il=idx((i-1)*n+1:(i-1)*n+n);
    mil=mode(il);
    ac(i)=sum(il==mil)/n;
end
sum(ac)/sub
% imshow(Y1(:,:,:,1))
real_idx=zeros(1,100);
for i = 1:10
    real_idx((i-1)*10 + 1:i*10) = i; % 将每 10 个元素赋值为 i
end
nmi=NMI(real_idx,idx)
