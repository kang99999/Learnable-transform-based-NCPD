function [ES_Y,iter_time,re_list,total_time,a,B] = RHALS(Y,M,R,maxiter)
%NTD_CP_L 此处显示有关此函数的摘要
%   此处显示详细说明
%initial
iniTotalTime = tic;
sz=size(Y);
N=length(sz);
w=randn(sz(1),M(1));
Yw=tensorprod(Y,w,1,1);
Yw=permute(Yw,[N,1:N-1]);
D=cell(N,1);
for n=2:N
    Z=classical_mode_unfolding(Yw,n);
%     [D{n},~,~]=svd(Z,"econ");
    [D{n},~,~]=eig(Z*Z');
    D{n}=D{n}(:,sz(n)-M(n)+1:sz(n));
end
X=Y;
for n=2:N
    X=tensorprod(X,D{n},2,1);
end
U=cell(N,1);
B=cell(N,1);
T1=1;
for n=1:N
    if n==1
        U{n}=rand(sz(n),R);
    elseif n>1
        U{n}=rand(M(n),R);
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
for n=1:N
    if n>1
            for r=1:R
                b=D{n}*U{n}(:,r);
                b(b<0)=1e-14;
                B{n}(:,r)=b;
            end
     else
            for r=1:R
                U{n}(U{n}(:,r)<0,r)=1e-14;
                B{n}=U{n};
            end
    end
end
idx=sort(2:N,"descend");
YB=B{idx(1)};
for j=2:N-1
    YB=khatri_rao(YB,B{idx(j)});
end
YB=B{1}*YB.';
ES_X=reshape(YB,size(Y));

re=norm(Y(:)-ES_X(:),"fro")/norm(Y(:),"fro");
re_list(1)=re;

for it =1:maxiter
    iniIterTime=tic;
    %update U
    lam=diag(U{1}.'*U{1});
    for n=1:N
%         if n==1
%             lam=1;
%         end
        idx=[sort(n+1:N,"descend"),sort(1:n-1,"descend")];
        kr_U=U{idx(1)};
        for j=2:N-1
            kr_U=khatri_rao(kr_U,U{idx(j)});
        end
        Xn=classical_mode_unfolding(X,n);
        T2=Xn*kr_U;
        T3=T1./(U{n}.'*U{n});
%         T4=T2/T3;
        if n>1
            for r=1:R
                U{n}(:,r)=U{n}(:,r)+(T2(:,r)-U{n}*T3(:,r))/lam(n);
                s=norm(U{n}(:,r),"fro");
                U{n}(:,r)=U{n}(:,r)/s;
                b=D{n}*U{n}(:,r);
                b(b<0)=1e-14;
                B{n}(:,r)=b;
                U{n}(:,r)=D{n}'*b;
%                 U{n}(:,r)=T4(:,r);
            end
%             B=D{n}*U{n};
%             B(B<0)=-1e-14;
%             U{n}=D{n}'*B;

        else
            for r=1:R
                U{n}(:,r)=U{n}(:,r)+T2(:,r)-U{n}*T3(:,r);
                U{n}(U{n}(:,r)<0,r)=1e-14;
%                 U{n}(:,r)=T4(:,r);
                B{n}=U{n};
            end
        end
        
%         if n<N
%             for r=1:R
%                 s=norm(U{n}(:,r),"fro");
% %                 if s~=0
%                  U{n}(:,r)=U{n}(:,r)/s;
% %                 end
%             end
%         end
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
    ES_X=reshape(YB,size(X));
%     ES_Y=ES_X;
%     for n=2:N
%         ES_Y=tensorprod(ES_Y,D{n},2,2);
%     end
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
    re=norm(X(:)-ES_X(:),"fro")/norm(X(:),"fro");
%     ES_Y=ES_X;
%     for n=2:N
%         ES_Y=tensorprod(ES_Y,D{n},2,2);
%     end
%     re2=norm(Y(:)-ES_Y(:),"fro")/norm(Y(:),"fro");
    re_list(it+1)=re;
    fit=1-re;
    if it==1
        oldfit=fit;
    else
        fitchange=fit-oldfit;
        
        if (abs(fitchange)<1e-6 && fit>0)||(it==maxiter)
            fprintf('\tRelative error below tol \n');
            ES_Y=ES_X;
            for n=2:N
                ES_Y=tensorprod(ES_Y,D{n},2,2);
            end
            if it==maxiter
                a=0;
            else
                a=1;
            end
            break
        else
            fprintf('\tRelative error after iteration %d: %.8f\n', it, abs(fitchange));
            oldfit=fit;
            
        end
    end
end
total_time = toc(iniTotalTime);
