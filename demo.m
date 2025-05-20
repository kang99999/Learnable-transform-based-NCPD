clc,clear
sz=[50 50 50];
R=10;
S= 0.9;
SNR=10;
for i=1:1
%     [B,Y]=gendata_ncp(sz,R,S,SNR);
    load("F:/NCPD/实验记录/Indian_pines.mat");
    Y=indian_pines;
    M=[30 30 30];
    M1=[30 30 30];
    ROU=[0.01 0.1 0.01];
    mu=10;
    maxiter=5000;
    [ES_Y,iter_time,re_list,total_time,a,hatB] = T_HALS(Y,M,R,ROU,mu,maxiter);
%     [ES_Y,iter_tirme,re_list,total_time,a,hatB] = iT_HALS(Y,M,R,ROU,mu,maxiter);
%    [ES_Y,iter_time,re_list,total_time,a,hatB] = AIBCD(Y,R,maxiter);
%     [ES_Y,iter_time,re_list,total_time,a,hatB] = FHALS(Y,R,maxiter);
%     [ES_Y,iter_time,re_list,total_time,a,hatB] = RHALS(Y,M1,R,maxiter);

%     sir=MeanSIR(B,hatB);
%     sirlist(i)=sir;
    r_e=norm(Y(:)-ES_Y(:))/norm(Y(:));
    relist(i)=r_e;

    psnrlist(i)=PSNR(Y,ES_Y);
    tmlist(i)=total_time;
end

