function [L, U, P, Q, x] = complete_pivoting_gaussian(A, b)
    % 全主元高斯消元法
    % A: 系数矩阵 (n x n)
    % b: 右端向量 (n x 1)
    % L: 单位下三角矩阵 (n x n)
    % U: 上三角矩阵 (n x n)
    % P: 行置换矩阵 (n x n)
    % Q: 列置换矩阵 (n x n)
    % x: 解向量 (n x 1)

    n = size(A, 1); % 矩阵的阶数

    % 初始化
    L = eye(n); % 单位下三角矩阵
    U = A;      % 初始化为 A
    P = eye(n); % 行置换矩阵
    Q = eye(n); % 列置换矩阵

    % 全主元高斯消元过程
    for k = 1:n-1
        % 找到当前子矩阵中绝对值最大的元素及其位置
        [max_val, max_index] = max(abs(U(k:n, k:n)), [], 'all', 'linear');
        [pivot_row, pivot_col] = ind2sub([n-k+1, n-k+1], max_index);
        pivot_row = pivot_row + k - 1; % 调整行索引
        pivot_col = pivot_col + k - 1; % 调整列索引

        % 交换行
        U([k, pivot_row], :) = U([pivot_row, k], :);
        P([k, pivot_row], :) = P([pivot_row, k], :);
        if k > 1
            L([k, pivot_row], 1:k-1) = L([pivot_row, k], 1:k-1);
        end

        % 交换列
        U(:, [k, pivot_col]) = U(:, [pivot_col, k]);
        Q(:, [k, pivot_col]) = Q(:, [pivot_col, k]);

        % 消元过程
        for i = k+1:n
            L(i, k) = U(i, k) / U(k, k); % 更新 L 矩阵
            U(i, k:n) = U(i, k:n) - L(i, k) * U(k, k:n); % 更新 U 矩阵
        end
    end

    % 求解 Ly = Pb
    y = forward_substitution(L, P * b);

    % 求解 Uz = y
    z = backward_substitution(U, y);

    % 恢复解向量 x = Qz
    x = Q * z;
end