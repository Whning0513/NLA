% test_lu.m
% 测试前代法、后代法和高斯消去法

% 清空工作区
clear;
clc;

% 1. 添加 solving and checking 文件夹到 MATLAB 路径
addpath('./solving');
addpath('./checking');

A = [4, 12, -16; 12, 37, -43; -16, -43, 1118]; % 对称正定矩阵
[L,D] =ldlt_cholesky_decomposition(A);
disp('L:'); disp(L);
disp('D:'); disp(D);
disp('L * L'':'); disp(L * L'); % 验证 A = LL^T