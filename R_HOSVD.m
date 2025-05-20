function [U] = R_HOSVD(Y,R)
%R_HOSVD 此处显示有关此函数的摘要
%   此处显示详细说明
sz=size(Y);
U=cell(ndims(Y),1);
for i=1:ndims(Y)
X=classical_mode_unfolding(Y,i);
% [U{i},~,~]=svd(X*randn(prod(sz)/sz(i),R(i)),"econ");
[u,~,~]=eig(X*X');
U{i}=u(:,sz(i)-R(i)+1:sz(i));
end

end
