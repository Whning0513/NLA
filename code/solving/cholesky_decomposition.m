function L = cholesky_decomposition(A)
    % Cholesky 分解
    % A: 对称正定矩阵 (n x n)
    % L: 下三角矩阵，满足 A = LL^T

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

    % 初始化下三角矩阵 L
    L = zeros(n, n);

    % Cholesky 分解
    for i = 1:n
        for j = 1:i
            % 计算 L(i, j)
            if i == j
                % 对角线元素
                L(i, j) = sqrt(A(i, j) - sum(L(i, 1:j-1).^2));
            else
                % 非对角线元素
                L(i, j) = (A(i, j) - L(i, 1:j-1) * L(j, 1:j-1)') / L(j, j);
            end
        end
    end
end