function result = allLeaPriMinorNot0(A)
    % 判断矩阵 A 是否所有顺序主子式非零
    % A: 输入的 n x n 矩阵
    % result: 逻辑值，true 表示所有顺序主子式非零，false 表示至少有一个顺序主子式为零
    
    n = size(A, 1); % 获取矩阵的阶数
    
    % 检查是否为方阵
    if size(A, 2) ~= n
        error('输入矩阵必须是方阵');
    end
    
    % 逐个计算顺序主子式
    for k = 1:n
        % 提取前 k 行和前 k 列的子矩阵
        submatrix = A(1:k, 1:k);
        
        % 计算子矩阵的行列式
        det_submatrix = det(submatrix);
        
        % 如果行列式为零，返回 false
        if abs(det_submatrix) < eps
            result = false;
            return;
        end
    end
    
    % 如果所有顺序主子式都非零，返回 true
    result = true;
end