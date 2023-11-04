function [X, k] = Newton(F, N, n, p, xnew, ynew, znew)
% 表示最终输出的解 k表示迭代次数
%F表示给定的非线性方程组,N表示使用的矩阵范数（在最开始检验时均使用矩阵无穷范数）p为误差 n为有几个未知变量
%L表示定义域的范围,其中定义域要为闭凸集，若原点不在集合中，L的值要随之修改（多加几个变量）
%在输入参数时，由于矩阵A中带有x，因此需要先输入"syms xnew real" 从而对x进行定义，否则报错，其余参数同理(real 防止转置后出现复数)
%在输入变量时，格外建议参数使用形如"xnew"的形式，避免后续程序因为变量重复而出错
%由于本人能力有限，在变量数目增多时需要手动更改程序的几处与变量数目有关的地方（已经在程序中注明）
%如果后面的同学有更好的思路欢迎分享
%Newton迭代法的收敛条件判定涉及的参数过多，比较繁琐，因此没有给出代码，建议手动验证
%syms xnew real; syms ynew real; syms znew real; F = [6*xnew-2*cos(ynew*znew)-1;sqrt(xnew^2+sin(znew)+1.06)-9*(ynew+0.1);3*exp(-xnew*ynew)+60*znew+10*pi-3]
%N = Inf; L = 1; p = 10^(-8); n=3;
%求Jacobi矩阵
J = sym(zeros(1, n)); 
for i = 1:n
    J(1, i) = sum(F(i,:));
end
FJocobi = jacobian(J, [xnew;ynew;znew]);          %注意xyz的数量要随着方程更改！！！
k = 0;                                %记录迭代次数
x1 = zeros(n, 1);
x2 = ones(n, 1);
X = zeros(n, 1);
delta = zeros(n, 1);
while norm(x2 - x1, N) >= p
     x2 = x1;
     xx = x2(1); yy = x2(2); zz = x2(3);   %注意修改变量个数！！！
     B = F;
     F = subs(F,xnew,xx); F = subs(F,ynew,yy);F = subs(F,znew,zz);
     BJocobi = FJocobi;
     FJocobi = subs(FJocobi,xnew,xx); FJocobi = subs(FJocobi,ynew,yy);FJocobi = subs(FJocobi,znew,zz);
     FJocobi = double(FJocobi); F = double(-F);
     [delta, c] = congrad(FJocobi, F, p, N);
     x1 = x2 + delta;
     x1 = double(x1);
     F = B;
     FJocobi= BJocobi;
     k = k + 1;
end
X = x1;
X
k
end
%tic; Newton(F, N, n, p, xnew, ynew, znew);toc;
