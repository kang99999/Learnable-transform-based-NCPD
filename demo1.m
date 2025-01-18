clc,clear
% sz=[70 70 70 70];
% R=4;
% S= 0.9;
% SNR=10;
%     Y=gendata_ncp(sz,R,S,SNR);
A=0;

load("chicago.mat")
%     Y=Y1;
M=[40 20 20];
R=10;
ROU=[1e-5 0.1 1e-2];
mu=1e5;
maxiter=5000;
% [ES_Y,iter_time,re_list,total_time,a,U] = T_HALS(Y,M,R,ROU,mu,maxiter);
% [ES_Y,iter_time,re_list,total_time,a,U] = iT_HALS(Y,M,R,ROU,mu,maxiter);


sir=SIR(Y,ES_Y);
r_e=norm(Y(:)-ES_Y(:))/norm(Y(:));
% psnr=PSNR(Y,ES_Y);

% sub=10;
% n=10;
% 
% 
% ac=zeros(1,sub);
% % 
% % 
% for i=1:sub
%     il=idx((i-1)*n+1:(i-1)*n+n);
%     mil=mode(il);
%     ac(i)=sum(il==mil)/n;
% end
% sum(ac)/sub
% % imshow(Y1(:,:,:,1))
% real_idx=zeros(1,100);
% for i = 1:10
%     real_idx((i-1)*10 + 1:i*10) = i; % 将每 10 个元素赋值为 i
% end
% nmi=NMI(real_idx,idx)
