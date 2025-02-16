% test_lu.m
% 测试前代法、后代法和高斯消去法

% 清空工作区
clear;
clc;

% 1. 生成一个满秩随机矩阵 A 和向量 b
n = 5; % 矩阵大小
A = randn(n, n); % 随机生成 n x n 矩阵

% 确保 A 是满秩的
while rank(A) < n
    A = randn(n, n); % 如果秩不足，重新生成
end

b = randn(n, 1); % 随机生成 n x 1 向量

% 2. 使用高斯消去法求解 Ax = b，并保留 L 和 U
[L, U, x] = gaussEli_sol(A, b);

% 3. 显示结果
disp('随机矩阵 A:');
disp(A);
disp('下三角矩阵 L:');
disp(L);
disp('上三角矩阵 U:');
disp(U);
disp('解向量 x:');
disp(x);

% 4. 验证结果
disp('验证 Ax - b 的误差:');
disp(norm(A * x - b)); % 计算误差

% 5. 测试秩不足的情况
disp('测试秩不足的情况:');
A_singular = A; % 复制 A
A_singular(:, end) = A_singular(:, 1); % 使最后一列与第一列相同，制造秩不足
disp('秩不足矩阵 A_singular:');
disp(A_singular);

try
    [L_singular, U_singular, x_singular] = gaussEli_sol(A_singular, b);
    disp('解向量 x_singular:');
    disp(x_singular);
catch ME
    disp('高斯消元法失败，矩阵秩不足:');
    disp(ME.message);
end