function kr = khatri_rao(F, G)
    nR_F = size(F, 1); % F 的行数
    nR_G = size(G, 1); % G 的行数
    mul = ones(nR_G, 1);
    FF = kron(F, mul); % 通过 kron 函数实现对 F 矩阵的扩充
    GG = repmat(G, nR_F, 1); % 通过 repmat 函数实现对 G 矩阵的扩充
    kr = FF .* GG; % Hadamard 积（逐元素乘积）
end