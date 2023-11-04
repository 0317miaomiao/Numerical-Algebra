%生成矩阵A
A = zeros(50,50);
for i = 1:50
    for j = 1:50
        if i==j
            A(i,j) = 50*i;
        else
            A(i,j) = max(i, j);
        end
    end
end
%生成向量b
b = zeros(50, 1);
for i = 1:50
    for j = 1:50
        b(i) = b(i) + A(i,j)*(51 - j);
    end
end
x = diag(A);
M = diag(x);
