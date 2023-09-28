function x = qrfactzj(A, B)
  [m n] =size(A);
  Q = eye(m);   %生成单位矩阵
  R = zeros(n);
  for j = 1:n      
      [v, beta] = Householder(A(j:m, j));   %不断的对A矩阵做householder变换，但不断向右下角缩小                                              %这里是第一次缩短计算时间，因为将v(1)化1
      A(j:m, j:n) = (eye(m -j +1) - beta*v*v')*A(j:m, j:n); %储存R，将R储存在A的右上角,并且由于H矩阵的性质（左上角1：j-1为单位阵）只用储存j：n的部分即可，其余部分在后续计算中不改变
      d(j) = beta;
      if j < m
          A(j + 1:m, j) = v(2: m - j + 1);     %储存Q,通过储存v的方式，由于v对应的向量第一位一直为1，因此只储存后面的，刚好储存在A的左下角（不算对角线）
      end
  end
  R = triu(A);    %从储存的A中获取R
  for k = 1:n
      H = eye(m);
      H1 = eye(m - k +1) - d(k)*[1, A(k+1:m, k)']'*[1, A(k+1:m, k)'];    %计算Q 因为Q = H1*H2...*Hn，因此首先计算每个H，再乘到一起
      H(k:m, k:m) = H1;                                                %每个H1只是n*n矩阵的右下角的一部分（详细见Householder变换）
      Q = Q*H;
  end
  x = zeros(n, 1);                                                     % Hn*...H2*H1*A = R, 所以Q = H1*H2...*Hn（H为正交矩阵）
  B = Q'*B;
  for j = n:-1:1
    x(j) = (B(j) - sum(R(j, j+1:n).*x(j+1:n)'))/R(j, j);
  end
end
