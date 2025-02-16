function [L,D] = ldlt_cholesky_decomposition(A)
    % LDLT 分解
    % A: 对称正定矩阵 (n x n)
    % L: 单位下三角矩阵 (n x n)
    % D: 对角矩阵 (n x n)

    % 检查矩阵是否为方阵
    [n, m] = size(A);
    if n ~= m
        error('输入矩阵必须为方阵');
    end

    % 检查矩阵是否对称
    if ~isequal(A, A')
        error('输入矩阵必须对称');
    end

    % 检查矩阵是否正定
    % 正定矩阵的所有特征值必须大于零
    eigenvalues = eig(A);
    if any(eigenvalues <= 0)
        error('输入矩阵必须正定');
    end

    % 初始化
    L = eye(n); % 单位下三角矩阵
    D = zeros(n, n); % 对角矩阵

    % LDLT 分解
    for j = 1:n
        % 计算 D(j, j)
        D(j, j) = A(j, j) - L(j, 1:j-1) * D(1:j-1, 1:j-1) * L(j, 1:j-1)';

        % 计算 L(i, j) for i = j+1:n
        for i = j+1:n
            L(i, j) = (A(i, j) - L(i, 1:j-1) * D(1:j-1, 1:j-1) * L(j, 1:j-1)') / D(j, j);
        end
    end
end