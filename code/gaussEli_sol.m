function [L, U, x] = gaussEli_sol(A, b)
    % A: 系数矩阵 (n x n)
    % b: 右端向量 (n x 1)
    % L: 下三角矩阵 (n x n)
    % U: 上三角矩阵 (n x n)
    % x: 解向量 (n x 1)
    
    n = length(b);
    L = eye(n); % 初始化 L 为单位矩阵
    U = A;      % 初始化 U 为 A
    
    % 高斯消元过程（带部分选主元）
    for k = 1:n-1
        % 部分选主元：找到第 k 列中绝对值最大的行
        [~, pivot_row] = max(abs(U(k:n, k)));
        pivot_row = pivot_row + k - 1; % 调整行索引
        
        % 交换当前行和主元行
        U([k, pivot_row], :) = U([pivot_row, k], :);
        b([k, pivot_row]) = b([pivot_row, k]);
        
        % 更新 L 矩阵
        if k > 1
            L([k, pivot_row], 1:k-1) = L([pivot_row, k], 1:k-1);
        end
        
        % 消元过程
        for i = k+1:n
            L(i, k) = U(i, k) / U(k, k); % 更新 L 矩阵
            U(i, k:n) = U(i, k:n) - L(i, k) * U(k, k:n); % 更新 U 矩阵
        end
    end
    
    % 求解 Ly = b
    y = forward_substitution(L, b);
    
    % 求解 Ux = y
    x = backward_substitution(U, y);
end