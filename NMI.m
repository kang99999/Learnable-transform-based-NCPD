function nmi = NMI(X, Y)
    % 计算两个聚类结果的归一化互信息 (NMI)
    % 输入:
    %   X - 第一个聚类结果的标签向量
    %   Y - 第二个聚类结果的标签向量
    % 输出:
    %   nmi - 归一化互信息值

    % 计算联合概率分布
    [uniqueX, ~, idxX] = unique(X);
    [uniqueY, ~, idxY] = unique(Y);
    jointHist = accumarray([idxX, idxY], 1, [length(uniqueX), length(uniqueY)]);
    jointProb = jointHist / sum(jointHist(:));

    % 计算边缘概率分布
    marginalX = sum(jointProb, 2);
    marginalY = sum(jointProb, 1);

    % 计算互信息
    mutualInfo = 0;
    for i = 1:length(uniqueX)
        for j = 1:length(uniqueY)
            if jointProb(i, j) > 0
                mutualInfo = mutualInfo + jointProb(i, j) * log2(jointProb(i, j) / (marginalX(i) * marginalY(j)));
            end
        end
    end

    % 计算熵
    entropyX = -sum(marginalX .* log2(marginalX + eps));
    entropyY = -sum(marginalY .* log2(marginalY + eps));

    % 计算归一化互信息
    nmi = 2 * mutualInfo / (entropyX + entropyY);
end

