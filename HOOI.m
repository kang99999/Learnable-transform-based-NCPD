function [U] = HOOI(Y,R,maxiter)
%HOOI 此处显示有关此函数的摘要
%   此处显示详细说明
sz=size(Y);
N=length(sz);
U=cell(N,1);
for n=1:N
    U{n}=randn(sz(n),R(n));
end
for it=1:maxiter
    for n=1:N
        X=Y;
        for j=1:n-1
            X=tensorprod(X,U{j},1,1);
        end
        for j=n+1:N
            X=tensorprod(X,U{j},2,1);
        end
        Xn=classical_mode_unfolding(X,1);
        [u,~,~]=eig(Xn*Xn');
        U{n}=u(:,sz(n)-R(n)+1:sz(n));
%         [u,~,~]=svd(Xn);
%         U{n}=u(:,1:R(n));
    end
%     X=Y;
%     for n=1:N
%         X=tensorprod(X,U{n}*U{n}',1,1);
%     end
% %     e=X-Y;
%     re=norm(X(:)-Y(:))/norm(Y(:));
%     if re<1e-14
%         break
%     else
%         fprintf('\tRelative error after iteration %d: %.8f\n', it, re);
%     end
end

end

