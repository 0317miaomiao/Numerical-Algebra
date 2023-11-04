function [X, k] = Picard(A, b, N, L, p, xnew, ynew, znew)
%X 表示最终输出的解 k表示迭代次数
%A表示非线性方程组左端，b表示非线性方程组左端， N表示使用的矩阵范数（在最开始检验时均使用矩阵无穷范数）p为误差
%L表示定义域的范围,其中定义域要为闭凸集，若原点不在集合中，L的值要随之修改（多加几个变量）
%在输入参数时，由于矩阵A中带有x，因此需要先输入"syms xnew real" 从而对x进行定义，否则报错，其余参数同理(real 防止转置后出现复数)
%在输入变量时，格外建议参数使用形如"xnew"的形式，避免后续程序因为变量重复而出错
%由于本人能力有限，在变量数目增多时需要手动更改程序的几处与变量数目有关的地方（已经在程序中注明）
%如果后面的同学有更好的思路欢迎分享
%syms xnew real; syms ynew real; syms znew real; b = [xnew,ynew,znew]'; N = Inf; L = 1; 
%A = [(1/3)*cos(ynew*znew) + 1/6; sqrt(xnew^2 + sin(znew) + 1.06)/9 - 0.1; -exp(-xnew*ynew)/20 - pi/6 + 1/20]
%p = 10^(-8);
n = length(b);     %确定方程有几个参数
%求Jacobi矩阵
J = sym(zeros(1, n)); 
F1 = sym(zeros(1, n)); 
F2 = sym(zeros(1, n)); 
S1 = zeros(1, n);        %构造第一个储存矩阵，用于储存A的每一行元素绝对值之和的最大值
S2 = zeros(1, n);        %构造第二个储存矩阵，用于储存A的Jacobi矩阵的每一行元素绝对值之和的最大值
for i = 1:n
    J(1, i) = sum(A(i,:));
    F1(1, i) = sum(abs(A(i,:)));        %构造对应的无穷范数矩阵，便于求无穷范数
end
AJocobi = jacobian(J, [xnew;ynew;znew]);          %注意xyz的数量要随着方程更改！！！
for i = 1:n
    F2(1, i) = sum(abs(AJocobi(i,:)));
end
%验证是否收敛
%A矩阵范数
for i = 1:n
    f = -F1(i);                      %取负值便于求最大值
    xm = sym(ones(1,n));               %给符号变量赋初值
    for ii=1:n
       xm(ii)=['x' num2str(ii)];       %定义符号变量的形式为x1,...,xn
    end
    f = char(f);                     %转化为字符变量
    f = strrep(f,'xnew','x(1)');        %将x，y转化为向量的形式，便于后面求极值
    f = strrep(f,'ynew','x(2)');        %若有更多的变量使用相同的方法操作！！！
    f = strrep(f,'znew','x(3)'); 
    for iii = n:-1:1
       f = replace(f,['x' num2str(iii)],['x(' num2str(iii) ')']);
    end
    f = eval(['@(x)' char(f) ';']);
    f(1:n);                                                          %注意此处也要修改！！！
    lb = [-L,-L,-L]; ub = [L,L,L]; T = [];b = [];Aeq = [];beq = [];  %求极值函数的对应参数设置
    x0 = (lb + ub)/2;
    [xd, fval] = fmincon(f,x0,T,b,Aeq,beq,lb,ub);
    S1(i) = fval;                                                %记录每一列绝对值之和的值
    f = 0;
end
%Ajacobi矩阵范数
for i = 1:n
    f = -F2(i);                      %取负值便于求最大值
    xm = sym(ones(1,n));               %给符号变量赋初值
    for ii=1:n
       xm(ii)=['x' num2str(ii)];       %定义符号变量的形式为x1,...,xn
    end
    f = char(f);                     %转化为字符变量
    f = strrep(f,'xnew','x(1)');        %将x，y转化为向量的形式，便于后面求极值
    f = strrep(f,'ynew','x(2)');        %若有更多的变量使用相同的方法操作!！！
    f = strrep(f,'znew','x(3)');
    for iii = n:-1:1
       f = replace(f,['x' num2str(iii)],['x(' num2str(iii) ')']);
    end
    f = eval(['@(x)' char(f) ';']);
    f(1:n);                                                          %注意此处也要修改！！！
    lb = [-L,-L,-L]; ub = [L,L,L]; T = [];b = [];Aeq = [];beq = [];  %求极值函数的对应参数设置
    x0 = (lb + ub)/2;
    [xd, fval] = fmincon(f,x0,T,b,Aeq,beq,lb,ub);
    S2(i) = fval;                                                %记录每一列绝对值之和的值
    f = 0;
end
%开始计算
if  max(abs(S1)) <= L && max(abs(S2)) <= 1
    disp('Picard方法收敛')
    k = 0;                                %记录迭代次数
    x1 = zeros(n, 1);
    x2 = ones(n, 1);
    X = zeros(n, 1);
    while norm(x2 - x1, N) >= p
        x2 = x1;
        xx = x2(1); yy = x2(2); zz = x2(3);   %注意修改变量个数！！！
        B = A;
        A = subs(A,xnew,xx);
        A = subs(A,ynew,yy);
        A = subs(A,znew,zz);
        for i = 1:n
           x1(i) = sum(A(i,:));
        end
        A = B;
        k = k + 1;
    end
    X = x1;
else
    disp('Picard方法不收敛')
    return
end
X
k
end




