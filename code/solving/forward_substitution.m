function x = forward_substitution(L, b)
    % L: 下三角矩阵 (n x n)
    % b: 右端向量 (n x 1)
    % x: 解向量 (n x 1)
    
    n = length(b);
    x = zeros(n, 1); % 初始化解向量
    
    for i = 1:n
        x(i) = (b(i) - L(i, 1:i-1) * x(1:i-1)) / L(i, i);
    end
end
