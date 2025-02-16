function [L, U, P, Q, x] = partial_pivoting_gaussian(A, b)
    % 列主元高斯消元法
    % A: 系数矩阵 (n x n)
    % b: 右端向量 (n x 1)
    % L: 每个元绝对值都小于等于1的单位下三角矩阵 (n x n)
    % U: 对角线非零元的数目等于rankA的上三角矩阵 (n x n)
    % P: 含有n个車的1的行置换矩阵 (n x n)
    % Q: 含有n个車的1列置换矩阵 (n x n)
    % x: 解向量 (n x 1)
    
    n = size(A, 1); % 矩阵的阶数
    
    % 初始化
    L = eye(n); % 单位下三角矩阵
    U = A;      % 初始化为 A
    P = eye(n); % 行置换矩阵
    Q = eye(n); % 列置换矩阵
    
    % 列主元高斯消元过程
    for k = 1:n-1
        % 找到第 k 列中绝对值最大的元素所在的行
        [~, pivot_row] = max(abs(U(k:n, k)));
        pivot_row = pivot_row + k - 1; % 调整行索引
        
        % 交换kth row和pivot所在那一行
        U([k, pivot_row], :) = U([pivot_row, k], :);
        P([k, pivot_row], :) = P([pivot_row, k], :);
        if k > 1
            L([k, pivot_row], 1:k-1) = L([pivot_row, k], 1:k-1);
        end
        
        % 消元过程,把a_kk下方的元素换元作0
        for i = k+1:n
            L(i, k) = U(i, k) / U(k, k); % 更新 L 矩阵
            U(i, k:n) = U(i, k:n) - L(i, k) * U(k, k:n); % 更新 U 矩阵
        end
    end
    
    % 求解 Ly = Pb
    y = forward_substitution(L, P * b);
    
    % 求解 Ux = y
    x = backward_substitution(U, y);
end
