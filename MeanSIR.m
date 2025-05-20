% function [sir] = SIR(Y,e_Y)
% %SIR 此处显示有关此函数的摘要
% %   此处显示详细说明
% E=Y-e_Y;
% sir=norm(e_Y,'fro')^2/norm(E(:),'fro')^2;
% end
%  function SIR = SIRNMF2(XH,X)
%   X = bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));
%  
%  XH = bsxfun(@rdivide,XH,sqrt(sum(XH.^2,2)));
%  
%  G = X*XH'/ (X*X'); %X*XH'*pinv(X*X');
%  
%  [Gmax,maxind] = max(abs(G),[],2);
%  
%  SIR = 10*log10(Gmax.^2./sum((XH -bsxfun(@times,X(maxind,:),Gmax)).^2,2));
%  end
function sir=MeanSIR(A,Ahat)
sir=cell2mat(cellfun(@CalcSIR,A,Ahat,'uni',0));
sir=mean(mean(sir));
end
function [SIR maps] = CalcSIR(A,Aest)
% Sergio Cruces & Andrzej Cichocki
% A(:,maps) matches Aest properly.  --Added by Guoxu Zhou

% mean value should be extracted first --Added by Guoxu Zhou
A=real(A);
A=bsxfun(@minus,A,mean(A));
Aest=bsxfun(@minus,Aest,mean(Aest));

A=A*diag(1./(sqrt(sum(A.^2))+eps));
Aest=Aest*diag(1./(sqrt(sum(Aest.^2))+eps));

col=size(A,2);flag=zeros(1,col);
MSE=inf*ones(1,col);
for i=1:size(Aest,2)
    temp=min(sum(bsxfun(@minus,Aest(:,i),A).^2,1),...
        sum(bsxfun(@plus,Aest(:,i),A).^2,1));
    temp=max(temp,flag);
    [MSE(i),maps(i)]=min(temp);
    flag(maps(i))=inf;
end
SIR=-10*log10(MSE);
end

% function [meansir]=MeanSIR(A,hatA)
% meansir=0;
% for i=1:length(hatA)
%     [MSIR, ~] = minimize_F_norm(A{i}, hatA{i});
%     meansir=meansir+MSIR;
% end
% meansir=meansir/length(hatA);
% end
% 
% function [MSIR, permuted_B] = minimize_F_norm(A, B)
%     A=real(A);
%     A=bsxfun(@minus,A,mean(A));
%     B=bsxfun(@minus,B,mean(B));
%     A=A*diag(1./(sqrt(sum(A.^2))+eps));
%     B=B*diag(1./(sqrt(sum(B.^2))+eps));
%     % 计算成本矩阵
%     cost_matrix = zeros(size(A,2), size(B,2));
%     for i = 1:size(A,2)
%         for j = 1:size(B,2)
%             cost_matrix(i,j) = sum((A(:,i) - B(:,j)).^2);
%         end
%     end
%     
%     % 使用匈牙利算法求解指派问题
%     [cols] = solve_assignment(cost_matrix);
% %     根据匹配结果排列B的列
%     permuted_B = B(:, cols);
%     % 计算MEAN SIR
%     MSIR=0;
%     for i=1:size(A,2)
%         MSIR=MSIR+10*log10(norm(A(:,i))/norm(A(:,i)-permuted_B(:,i)));
%     end
%     MSIR=MSIR;
% 
% end
% 
% function [cols] = solve_assignment(cost_matrix)
%     % 使用匈牙利算法找到最小成本匹配
%     [row, col] = size(cost_matrix);
%     cost = zeros(row, col);
%     for i = 1:row
%         for j = 1:col
%             cost(i,j) = cost_matrix(i,j);
%         end
%     end
%     
%     % 使用匈牙利算法
%     [~, cols] = matchpairs(cost_matrix, 0);
%     
%     % 根据匹配结果排列B的列
% %     permuted_B = B(:, assignment);
% end