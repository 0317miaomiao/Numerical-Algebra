function [X, n] = precongrad(A, b, M, p, N)
%n 为迭代次数， X 为最终的迭代值
%A 为线性方程组，b为方程式右端值， p为收敛条件 M为预优矩阵
%N为矩阵范数的选择，输入1，2，Inf分别代表1范数，2范数与无穷范数
    n = length(b);
    X0 = zeros(n, 1);
    r0 = b - A*X0;
    zeta0 = M\r0;
    rho0 = r0'*zeta0;
    P0 = zeta0;
    while norm(r0, N) >= p
        n = n + 1;
        omega = A*P0;
        alpha = rho0/(P0'*omega);
        X0 = X0 + alpha*P0;
        r0 = r0 - alpha*omega;
        zeta1 = M\r0;
        rho1 = r0'*zeta1;
        beta = rho1/rho0;
        P0 = zeta1 + beta*P0;
        rho0 = rho1;
    end
    X = X0;
end

