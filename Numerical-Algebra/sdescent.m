function [X, n] = sdescent(A, b, p, N)
%n 为迭代次数， X 为最终的迭代值
%A 为线性方程组，b为方程式右端值， p为收敛条件 
%N为矩阵范数的选择，输入1，2，Inf分别代表1范数，2范数与无穷范数
    n = length(b);
    X0 = zeros(n, 1);
    r0 = b - A*X0;
    n = 0;
    while norm(r0, N) >= p
        n = n + 1;
        alpha0 = r0'*r0/(r0'*A*r0);
        X1 = X0 + alpha0*r0;
        r0 = b - A*X1;
        X0 = X1;
    end
    X = X0;
end

