function [ES_Y,iter_time,re_list,total_time,a,S,C] = JTKM(Y,R,eta,lam,mu,num_center,maxiter)
%NTD_CP_L 此处显示有关此函数的摘要
%   此处显示详细说明
%initial
iniTotalTime = tic;
sz=size(Y);
N=length(sz);

U=cell(N,1);
for i=1:N
    U{i}=rand(sz(i),R);
    
end
% D=zeros(sz(1));
for i=1:sz(1)
%     D(i,i)=norm(U{1}(i,:));
    U{1}(i,:)=U{1}(i,:)/norm(U{1}(i,:));
end
Z=U{1};

W=U{N};
WTW=U{N}'*U{N};
idx=sort(2:N-1,"descend");
for i=1:N-2
    W=khatri_rao(W,U{idx(i)});
    WTW=WTW.*(U{idx(i)}'*U{idx(i)});
end
Y1=reshape(Y,[sz(1),prod(sz(2:N))]);
Y1W=Y1*W;
% U{1}=(D\Y1W+lam*SC+mu*Z)/(WTW+(lam+mu)*eye(R));
% U{1}(U{1}<0)=-1e-14;
D=diag(diag(Y1W*U{1}')./diag(U{1}*WTW*(U{1}')));
% 规范化
% H=zeros(sz(N));
% F=zeros(size(B{N}));
% for i=1:sz(N)
%     H(i,i)=norm(B{N}(i,:));
%     F(i,:)=B{N}(i,:)/H(i,i);
% end
[S,C]=kmeans(U{1}',num_center);
SC=C(:,S);
SC=SC';
% oldU=cell(N,1);
% for n=1:N
%     oldU{n}=zeros(size(U{n}));
% end
%main part

re_list = zeros(1, maxiter);
iter_time = zeros(1, maxiter+1);
iter_time(1) = 0;

for it =1:maxiter
    iniIterTime=tic;
    %update U
    %U1
    Y1=reshape(Y,[sz(1),prod(sz(2:N))]);
    W=U{N};
    WTW=U{N}'*U{N};
    idx=sort(2:N-1,"descend");
    for i=1:N-2
        W=khatri_rao(W,U{idx(i)});
        WTW=WTW.*(U{idx(i)}'*U{idx(i)});
    end
    Y1W=Y1*W;
    U{1}=(D\Y1W+lam*SC+mu*Z)/(WTW+(lam+mu)*eye(R));
    U{1}(U{1}<0)=-1e-15;
    D=diag(diag(Y1W*U{1}')./diag(U{1}*WTW*(U{1}')));
    for i=1:sz(1)
        Z(i,:)=U{1}(i,:)/norm(U{1}(i,:));
    end
    [S,C]=kmeans(U{1}',num_center);
    SC=C(:,S);
    SC=SC';
%     SC=C(S,:);
    DU1=D*U{1};
%     update u
    for n=2:N
        idx=[sort(n+1:N,"descend"),sort(2:n-1,"descend")];
        G=U{idx(1)};
        GTG=U{idx(1)}'*U{idx(1)};
        for i=2:length(idx)
            G=khatri_rao(G,U{idx(i)});
            GTG=GTG.*(U{idx(i)}'*U{idx(i)});
        end
        G=khatri_rao(G,DU1);
        GTG=GTG.*(DU1'*DU1);
        Yn=classical_mode_unfolding(Y,n);
        U{n}=Yn*G/(GTG+eta*eye(R));
        U{n}(U{n}<0)=-1e-15;
    end

    iter_time(it+1) = iter_time(it) + toc(iniIterTime);
    % check convergence
    idx=sort(2:N,"descend");
    YB=U{idx(1)};
    for j=2:N-1
        YB=khatri_rao(YB,U{idx(j)});
    end
    YB=DU1*YB.';
    ES_Y=reshape(YB,sz);
    re=norm(Y(:)-ES_Y(:),"fro")/norm(Y(:),"fro");
    re_list(it)=re;
    fit=1-re;
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
