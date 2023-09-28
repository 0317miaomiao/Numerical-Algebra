%% 实验一
h = 0.01;                                   %步长为0.01
N = 100;                                    %划分为100份
n = 99;
e = ones(n, 1);
A = spdiags([-e, 1.9996*e, -e], -1:1, n, n);  %用于生成三对角型矩阵
A = full(A)                                 %转化为常规形式
% 生成b()
B = zeros(n, 1);
k = 1;
for i = 0.01:0.01:0.99
    B(k) = -h*h*(2*cos(pi*i) + 3*sin(pi*i))*(4 - pi*pi); 
    k = k+1;
end
B(1) = B(1) + 2;
B(99) = B(99) - 2;
%% 使用Cholesky分解计算
tic
x1 = Cholesky(A, B);
toc
%% 使用Thomas分解计算
a = -1*ones(n-1, 1);
b = 1.9996*ones(n, 1);
c = -1*ones(n-1, 1);
d = B;
tic
x2 = Thomas(a, b, c, d);
toc
x2 = vpa(x2, 5);                %保留5位有效数字
%% 使用Doolitle分解计算
tic
x3 = Doolitle(A, B);
toc
%% 使用QR分解计算
tic
x4 = qrfact(A, B);
toc
%% 计算运行速度
tic
x3 = Doolitle(A, B);
toc
%% 计算精确值
R = zeros(n, 1);
k = 1;
for i = 0.01:0.01:0.99
    R(k) = 2*cos(pi*i) + 3*sin(pi*i);
    k = k + 1;
end
%% 计算A的各阶顺序主子式，为Doolitle分解做准备
dett = zeros(99, 1);
for i = 1:99
    dett(i) = det(A(1:i, 1:i));
end
min(abs(dett))
%% 计算误差
 err = max(abs(R - x));