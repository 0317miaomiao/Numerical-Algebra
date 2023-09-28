function x = qrfact(A, B)
  n =length(A);
  Q = eye(n);   %生成单位矩阵
  R = zeros(n);
  for j = 1:n
      if j < n
          [v, beta] = Householder(A(j:n, j));   %不断的对A矩阵做householder变换，但不断向右下角缩小
      else                                        %这里是第一次缩短计算时间，因为将v(1)化1
          v = 1; 
          beta = 2 - 2*mod(n, 2);
      end
      A(j:n, j:n) = (eye(n -j +1) - beta*v*v')*A(j:n, j:n); %储存R，将R储存在A的右上角,并且由于H矩阵的性质（左上角1：j-1为单位阵）只用储存j：n的部分即可，其余部分在后续计算中不改变
      d(j) = beta;
      if j < n
          A(j + 1:n, j) = v(2: n - j + 1);     %储存Q,通过储存v的方式，由于v对应的向量第一位一直为1，因此只储存后面的，刚好储存在A的左下角（不算对角线）
      end
  end
  R = triu(A);    %从储存的A中获取R
  for k = 1:n
      H = eye(n);
      H1 = eye(n - k +1) - d(k)*[1, A(k+1:n, k)']'*[1, A(k+1:n, k)'];    %计算Q 因为Q = H1*H2...*Hn，因此首先计算每个H，再乘到一起
      H(k:n, k:n)=H1;                                                %每个H1只是n*n矩阵的右下角的一部分（详细见Householder变换）
      Q = Q*H;
  end
  x = zeros(n, 1);                                                     % Hn*...H2*H1*A = R, 所以Q = H1*H2...*Hn（H为正交矩阵）
  B = Q'*B;
  for j = n:-1:1
    x(j) = (B(j) - sum(R(j, j+1:n).*x(j+1:n)'))'/R(j, j);
  end
end

  
      
  

