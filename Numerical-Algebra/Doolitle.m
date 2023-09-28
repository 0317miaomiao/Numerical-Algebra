function x = Doolitle(A, B)
 [n, n] = size(A);
 L = zeros(n);
 U = zeros(n);
 x = zeros(n, 1);
 y = zeros(n, 1);
 for r = 1:n
     for i = r:n
         U(r, i) = A(r, i) - sum(L(r, 1:r-1).*U(1:r-1, i)');
         L(i, r) = (A(i, r)- sum(L(i, 1:r-1).*U(1:r-1, r)'))/U(r, r);
     end
 end;
 L;, U;                                             %即为所得的L与U
 for i= 1:n
     y(i) = B(i) - sum(L(i, 1:i-1).*y(1:i-1)');     %直接除的方式要慢
 end
 for j = n:-1:1
     x(j) = (y(j) - sum(U(j, j+1:n).*x(j+1:n)'))/U(j, j);
 end
end

