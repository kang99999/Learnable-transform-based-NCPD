function [ES_Y,iter_time,re_list,total_time,a,U] = AIBCD(Y,R,maxiter)
%NTD_CP_L 此处显示有关此函数的摘要
%   此处显示详细说明
%initial
iniTotalTime = tic;
sz=size(Y);
N=length(sz);
epsilon=1e-14;
beta=0.4;
U=cell(N,1);
oldU=cell(N,1);
for i=1:N
    oldU{i}=zeros(sz(i),R);
end
T1=1;
for n=1:N
    U{n}=rand(sz(n),R);
    if n<N
        for r=1:R
            U{n}(:,r)=U{n}(:,r)/norm(U{n}(:,r),'fro');
        end
    end
    T1=T1.*(U{n}.'*U{n});
end
%main part

re_list = zeros(1, maxiter);
iter_time = zeros(1, maxiter+1);
iter_time(1) = 0;
idx=sort(2:N,"descend");
YB=U{idx(1)};
for j=2:N-1
    YB=khatri_rao(YB,U{idx(j)});
end
YB=U{1}*YB.';
ES_Y=reshape(YB,sz);
re=norm(Y(:)-ES_Y(:),"fro")/norm(Y(:),"fro");
%     fit=1-re;
re_list(1)=re;
for it =1:maxiter
    iniIterTime=tic;
    %update U
    lam=diag(U{N}.'*U{N});
    for n=1:N
        if n==N
            lam=1;
        end
        idx=[sort(n+1:N,"descend"),sort(1:n-1,"descend")];
        kr_U=U{idx(1)};
        for j=2:N-1
            kr_U=khatri_rao(kr_U,U{idx(j)});
        end
        Yn=classical_mode_unfolding(Y,n);
        T2=Yn*kr_U;
        T3=T1./(U{n}.'*U{n});
%         T4=T2/T3;
     
        for m=1:15
        
        if n<N
            for r=1:R
                U{n}(:,r)=U{n}(:,r)+0.05*(T2(:,r)-U{n}*T3(:,r))/lam(n);
                U{n}(U{n}(:,r)<0,r)=epsilon;

                 s=norm(U{n}(:,r),"fro");
                U{n}(:,r)=U{n}(:,r)/s;
                
                if all(U{n}(:,r)<0)
                    c=U{n}(:,r);
                    maxc=max(c);
                    U{n}(:,r)=0;
                    U{n}(c==maxc,r)=1;
                else
                    U{n}(U{n}(:,r)<0,r)=-1e-14;
                    U{n}(:,r)=U{n}(:,r)/norm(U{n}(:,r),"fro");
                end
            end
        else
            for r=1:R
                U{n}(:,r)=U{n}(:,r)+0.05*(T2(:,r)-U{n}*T3(:,r));
                U{n}(U{n}(:,r)<0,r)=epsilon;

            end
        end
        

        end
        U{n}=U{n}+beta*(U{n}-oldU{n});
        oldU{n}=U{n};
        T1=T3.*(U{n}.'*U{n});
    end
    
    

    iter_time(it+1) = iter_time(it) + toc(iniIterTime);
    % check convergence
    idx=sort(2:N,"descend");
    YB=U{idx(1)};
    for j=2:N-1
        YB=khatri_rao(YB,U{idx(j)});
    end
    YB=U{1}*YB.';
    ES_Y=reshape(YB,sz);
    re=norm(Y(:)-ES_Y(:),"fro")/norm(Y(:),"fro");
    fit=1-re;
    re_list(it+1)=re;
%     if it==1
%         oldES_Y=ES_Y;
%     else
%         re=norm(oldES_Y(:)-ES_Y(:),"fro")/norm(oldES_Y(:),"fro");
%         re_list(it)=re;
%         if re<1e-4
%             fprintf('\tRelative error below tol \n');
%             
%             break
%         else
%             fprintf('\tRelative error after iteration %d: %.8f\n', it, re);
%             oldES_Y=ES_Y;
%         end
%     end
    if it==1
        oldfit=fit;
    else
        fitchange=fit-oldfit;
        
        if abs(fitchange)<1e-6 && fit>0
            fprintf('\tRelative error below tol \n');
            
            break
        else
            fprintf('\tRelative error after iteration %d: %.8f\n', it, abs(fitchange));
            oldfit=fit;
            if it==maxiter
                a=0;
            else
                a=1;
            end
        end
    end
end
total_time = toc(iniTotalTime);
