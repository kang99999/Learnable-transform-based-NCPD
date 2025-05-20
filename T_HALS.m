function [ES_Y,iter_time,re_list,total_time,a,B] = T_HALS(Y,M,R,ROU,mu,maxiter)
%NTD_CP_L 此处显示有关此函数的摘要
%   此处显示详细说明
%initial
iniTotalTime = tic;
sz=size(Y);
N=length(sz);
% D=cell(N,1);
X=Y;
D=R_HOSVD(Y,M);
% for n=1:N
%     D{n}=randn(sz(n),M(n));
% end
% for it=1:5
% for n=1:N
%     XX=Y;
%     for j =1:N
%         if j<n
%             XX=tensorprod(XX,D{n},1,1);
%         elseif j>n
%             XX=tensorprod(XX,D{n},2,1);
%         end
%     end
%     XX=classical_mode_unfolding(XX,1);
%     [D{n},~,~]=svd(XX);
%     D{n}=D{n}(:,1:M(n));
% %     X=tensorprod(X,D{n},1,1);
% %     Yn=classical_mode_unfolding(Y,n);
% %     [D{n},~,~]=svd(Yn);
% %     D{n}=D{n}(:,1:M(n));
% %     X=tensorprod(X,D{n},1,1);
% end
% end
for n=1:N
    X=tensorprod(X,D{n},1,1);
end
U=cell(N,1);
B=cell(N,1);
T1=1;
bestfit=-Inf;
bestB=0;
for n=1:N
    U{n}=rand(M(n),R);
    
    if n<N
        for r=1:R
            U{n}(:,r)=U{n}(:,r)/norm(U{n}(:,r),'fro');
        end
    end
    T1=T1.*(U{n}.'*U{n});
end
% oldU=cell(N,1);
% for n=1:N
%     oldU{n}=zeros(size(U{n}));
% end
%main part

re_list = zeros(1, maxiter);
iter_time = zeros(1, maxiter+1);
iter_time(1) = 0;
idx=sort(2:N,"descend");
for n=1:N
    if n~=N
        for r=1:R
            B{n}(:,r)=D{n}*U{n}(:,r);
            if all(B{n}(:,r)<0)
                c=B{n}(:,r);
                maxc=max(c);
                B{n}(:,r)=0;
                B{n}(c==maxc,r)=sqrt(1/diag(D{n}(c==maxc,:)*D{n}(c==maxc,:)'));
            else
                B{n}(B{n}(:,r)<0,r)=-1e-14;
                B{n}(:,r)=B{n}(:,r)/norm(B{n}(:,r),"fro");
            end 
            
        end
    else
        for r=1:R
            B{n}(:,r)=D{n}*U{n}(:,r);
            B{n}((B{n}(:,r)<0),r)=-1e-15;              

        end
    end
end
YB=B{idx(1)};
for j=2:N-1
    YB=khatri_rao(YB,B{idx(j)});
end
YB=B{1}*YB.';
ES_Y=reshape(YB,sz);
re=norm(Y(:)-ES_Y(:),"fro")/norm(Y(:),"fro");
re_list(1)=re;

for it =1:maxiter
    iniIterTime=tic;
    %update U
    lam=diag(U{N}.'*U{N});
%     print(lam)
%     print('\n')
    for n=1:N
%         oldU{n}=U{n};
%         w=0.01;
        if n==N
            lam=1;
        end
        idx=[sort(n+1:N,"descend"),sort(1:n-1,"descend")];
        kr_U=U{idx(1)};
        for j=2:N-1
            kr_U=khatri_rao(kr_U,U{idx(j)});
        end
        Xn=classical_mode_unfolding(X,n);
        T2=Xn*kr_U;
        T3=T1./(U{n}.'*U{n});
%         T4=T2/T3;
%         t1=1;
        for m=1:1
%         t2=1/2*(1+sqrt(1+4*t1^2));
%         w=(t1-1)/t2;
        if n<N
            for r=1:R
                U{n}(:,r)=U{n}(:,r)+(T2(:,r)-U{n}*T3(:,r))/(ROU(1)+lam(r));
%                  U{n}(:,r)=U{n}(:,r)+w*(U{n}(:,r)-oldU{n}(:,r));
%                 U{n}(:,r)=U{n}(:,r)/norm(U{n}(:,r),"fro");
%                 B{n}(:,r)=D{n}*U{n}(:,r);
%                 B{n}((B{n}(:,r)<0),r)=-1e-14;
%                 U{n}(:,r)=D{n}'*B{n}(:,r);
                B{n}(:,r)=D{n}*U{n}(:,r);
                if all(B{n}(:,r)<0)
                    c=B{n}(:,r);
                    maxc=max(c);
                    B{n}(:,r)=0;
                    B{n}(c==maxc,r)=sqrt(1/diag(D{n}(c==maxc,:)*D{n}(c==maxc,:)'));
                else
                    B{n}(B{n}(:,r)<0,r)=-1e-14;
                    B{n}(:,r)=B{n}(:,r)/norm(B{n}(:,r),"fro");
                end
                
                U{n}(:,r)=D{n}'*B{n}(:,r);
%                 oldU{n}(:,r)=U{n}(:,r);
%                 uu=U{n}(:,r);
%                 U{n}(:,r)=U{n}(:,r)+w*(U{n}(:,r)-oldU{n}(:,r));
%                 oldU{n}(:,r)=uu;
            end
        else
            for r=1:R
                U{n}(:,r)=U{n}(:,r)+(T2(:,r)-U{n}*T3(:,r))/(ROU(1)+1);
%                  U{n}(:,r)=U{n}(:,r)+w*(U{n}(:,r)-oldU{n}(:,r));
                B{n}(:,r)=D{n}*U{n}(:,r);
                B{n}((B{n}(:,r)<0),r)=-1e-15;
%                 
                U{n}(:,r)=D{n}'*B{n}(:,r);
%                 oldU{n}(:,r)=U{n}(:,r);
%                 uu=U{n}(:,r);
%                 U{n}(:,r)=U{n}(:,r)+w*(U{n}(:,r)-oldU{n}(:,r));
%                 oldU{n}(:,r)=uu;
            end
        end
        
%         B{n}=D{n}*U{n};
%         B{n}(B{n}<0)=-1e-14*rand(length(B{n}(B{n}<0)),1);
%         B{n}(B{n}<0)=-1e-14;
%         if n<N
%             for r=1:R
%                 if all(B{n}(:,r)<0)
%                     c=B{n}(:,r);
%                     maxc=max(c);
%                     B{n}(:,r)=0;
%                     B{n}(c==maxc,r)=sqrt(1/diag(D{n}(c==maxc,:)*D{n}(c==maxc,:)'));
%                 else
%                     B{n}(B{n}<0)=-1e-14;
%                     B{n}(:,r)=B{n}(:,r)/norm(B{n}(:,r),"fro");
%                 end
%             end
%         else
%             B{n}(B{n}<0)=-1e-14;
%         end
%         U{n}=D{n}'*B{n};
        end
%         U{n}=U{n}+w*(U{n}-oldU{n});
%         if n<N
%             for r=1:R
%                 U{n}(:,r)=U{n}(:,r)/norm(U{n}(:,r),"fro");
%             end
%         end
        T1=T3.*(U{n}.'*U{n});
    end
    
    
    %update X
    idx=sort(2:N,"descend");
    XU=U{idx(1)};
    for j=2:N-1
        XU=khatri_rao(XU,U{idx(j)});
    end
    XU=U{1}*XU.';
    XU=reshape(XU,M);
    YD1=Y;
    for n=1:N
        YD1=tensorprod(YD1,D{n},1,1);
    end
    X=(XU+mu*YD1+ROU(2)*X)/(1+mu+ROU(2));
    %update D
    for n=1:N
    %     idx=[sort(n+1:N,"descend"),sort(1:n-1,"descend")];
    %     Xn=classical_mode_unfolding(X,n);
%         XD=X;
%         for j=1:n-1
%             XD=tensorprod(XD,D{j},1,2);
%         end
%         for j=n+1:N
%             XD=tensorprod(XD,D{j},2,2);
%         end
%         XD=reshape(XD,[M(n),prod(sz)/sz(n)]);
%         Yn=classical_mode_unfolding(Y,n);
%         f=mu/2*Yn*XD.'+ROU(3)/2*D{n};
        YD1=tensorprod(YD1,D{n},n,2);
        YD1=permute(YD1,[1:n-1,N,n:N-1]);
        Yn=classical_mode_unfolding(YD1,n);
        Xn=classical_mode_unfolding(X,n);
        f=mu/2*Yn*Xn.'+ROU(3)/2*D{n};
        [P,~,V]=svd(f);
        D{n}=P(:,1:M(n))*V';
%         D{n}=
        YD1=tensorprod(YD1,D{n},n,1);
        YD1=permute(YD1,[1:n-1,N,n:N-1]);
    end

    iter_time(it+1) = iter_time(it) + toc(iniIterTime);
    % check convergence
    idx=sort(2:N,"descend");
    YB=B{idx(1)};
    for j=2:N-1
        YB=khatri_rao(YB,B{idx(j)});
    end
    YB=B{1}*YB.';
    ES_Y=reshape(YB,sz);
    re=norm(Y(:)-ES_Y(:),"fro")/norm(Y(:),"fro");
    re_list(it+1)=re;
    fit=1-re;

    if fit>=bestfit
        bestfit=fit;
        bestB=B;
    end
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
%     idx=sort(2:N,"descend");
%     X1=U{idx(1)};
%     for j=2:N-1
%         X1=khatri_rao(X1,U{idx(j)});
%     end
%     X1=U{1}*X1.';
%     ES_X=reshape(X1,M);
%     re=norm(X(:)-ES_X(:),"fro")/norm(X(:),"fro");
%     if re<1e-14
%         fprintf('\tRelative error below tol \n');
%         
%         break
%     else
%         fprintf('\tRelative error after iteration %d: %.8f\n', it, re);
% %         oldES_Y=ES_Y;
%     end
end
total_time = toc(iniTotalTime);
