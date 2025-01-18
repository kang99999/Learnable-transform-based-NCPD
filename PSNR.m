function [p] = PSNR(Y,E_Y)
%PSNR 此处显示有关此函数的摘要
%   此处显示详细说明
mse=(norm(Y(:)-E_Y(:)))^2/numel(Y);
m=max(E_Y(:));
p=10*log10(m^2/mse);
end

