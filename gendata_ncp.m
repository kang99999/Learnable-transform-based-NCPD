function [B,Y] = gendata_ncp(sz,R,S,SNR)
%GENDATA_NCP 此处显示有关此函数的摘要
%   此处显示详细说明
N=length(sz);
B=cell(N,1);
for n=1:N
    B{n}=rand(sz(n),R);
    B{n}(B{n}<S)=0;
end

idx=sort(2:N,"descend");
y1=B{idx(1)};
for n=2:N-1
    y1=khatri_rao(y1,B{idx(n)});
end
Y1=B{1}*y1.';
signal_power = 1; % 信号功率，可以根据实际情况调整
noise_power = signal_power / (10^(SNR/10)); % 噪声功率
% 计算噪声矩阵 E
E=randn(size(Y1));
rate=norm(Y1(:),'fro')/norm(E(:),'fro');
E = rate*sqrt(noise_power) * E;
Y=E+Y1;
% snr_check = 10 * log10(norm(Y1(:), 'fro')^2 / norm(E(:), 'fro')^2);
Y=reshape(Y,sz);
end

% SNR_dB = 20; % 信噪比，单位为分贝
% signal_power = 1; % 信号功率，可以根据实际情况调整
% noise_power = signal_power / (10^(SNR_dB/10)); % 噪声功率
% 
% % 随机生成非负矩阵 A 和 X
% A = rand(5, 3); % 示例维度，根据需要调整
% X = rand(3, 4);
% 
% % 计算信号矩阵 AX
% AX = A * X;
% 
% % 计算噪声矩阵 E
% E = sqrt(noise_power) * randn(size(AX));
% 
% % 生成最终的 Y
% Y = AX + E;
% 
% % 检查SNR是否接近20dB
% snr_check = 10 * log10(norm(AX, 'fro')^2 / norm(E, 'fro')^2);
% disp(['Calculated SNR: ', num2str(snr_check), ' dB']);


