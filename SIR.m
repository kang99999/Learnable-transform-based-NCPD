function [sir] = SIR(Y,e_Y)
%SIR 此处显示有关此函数的摘要
%   此处显示详细说明
E=Y-e_Y;
sir=norm(e_Y,'fro')^2/norm(E(:),'fro')^2;
end

