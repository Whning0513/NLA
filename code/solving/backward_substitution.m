function x = backward_substitution(U, b)
    % U: 上三角矩阵 (n x n)
    % b: 右端向量 (n x 1)
    % x: 解向量 (n x 1)
    
    n = length(b);
    x = zeros(n, 1); % 初始化解向量
    
    for i = n:-1:1
        x(i) = (b(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
end