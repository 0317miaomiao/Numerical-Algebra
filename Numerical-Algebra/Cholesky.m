function x = Cholesky(A, B)    %高阶情形快于LU分解
  n = length(B);
  v = zeros(n);
  x = zeros(n, 1);
  y = zeros(n, 1);
  for j = 1:n
      for i = 1:j-1
          v(j, i) = A(j, i)*A(i, i);
      end
      A(j, j) = A(j, j) - A(j, 1:j-1)*v(j, 1:j-1)';  %因为在后续计算中a这个矩阵前面的数字已经没有用了，所以为了简化计算时间，直接将L和D矩阵的值填到a中
      A(j+1:n, j) = (A(j+1:n, j) - A(j+1:n, 1:j-1)*v(j, 1:j-1)')/A(j, j); %A(j, j)代表了d(j), A(j+1:n, j)代表了L（i，j）
  end
  L = tril(A, -1) + eye(n);, U = diag(diag(A))*L';
  for i = 1:n
      y(i) = B(i) - L(i, 1:i-1)*y(1: i-1);
  end
  for j = n:-1:1
      x(j) = (y(j) - U(j, j+1:n)*x(j+1:n))/U(j, j);
  end
end
      