function [x, n, err] = Jacobi(A, b, p, k, N, wp, r)
%n 为迭代次数， x 为最终的迭代值 err为误差值
%A 为线性方程组，b为方程式右端值， p为收敛条件 
%k为迭代方法的选择，如果k=1则使用Jacobi方法，若k = 2则使用JOR迭代法
%N为矩阵范数的选择，输入1，2，Inf分别代表1范数，2范数与无穷范数
%r为精确解，用来计算误差,如果没有精确值则不输入
%wp为如果没有最佳松弛因子时自己选择输入的松弛因子
% 初始解始终设为 x1 = 0
  U = triu(A, 1);     %提出上三角
  L = tril(A, -1);     %提出下三角
  d = diag(A);     %提取对角线
  D = diag(d);     %组成矩阵
  n = length(A);
  B = -D\(L + U);    %求出迭代矩阵B
  TB = max(abs(eig(B)));   %求出谱半径
  if k == 1
      disp('使用Jacobi方法')
      if TB < 1
         disp('Jacobi方法收敛')
         x1 = zeros(n,1);     %构造初始值
         x2 = zeros(n, 1) + 0.5;
         n = 0;               %记录迭代次数
         if exist('r', 'var')              %根据是否有精确解来计算误差，如果没有则根据后验误差计算
            while norm((r - x1), N) >= p   %判断收敛条件
               x2 = x1;
               x1 = B*x2 + D\b;
               n = n + 1;
            end
            err = norm(x1 - r);     %计算误差值，统一采用矩阵二范数计算
         else
             while norm((x2 - x1), N) >= p   %判断收敛条件
               x2 = x1;
               x1 = B*x2 + D\b;
               n = n + 1;
            end
            err = norm(x1 - x2);     %计算误差值，统一采用矩阵二范数计算
         end
         x = x1;
      else
          disp('Jacobi方法不收敛')
      end
  elseif k == 2
      disp('使用JOR方法')
      if isreal(eig(B)) && TB < 1
          disp('JOR方法存在最佳松弛因子')
          MB = max(eig(B));     %计算松弛因子
          mB = min(eig(B));
          w = 2/(2 - MB - mB);
          BJ = eye(n) - w*(D\A);  %计算JOR迭代中的B矩阵
          if max(abs(eig(BJ))) < 1   %判断JOR方法是否收敛
              disp('JOR迭代法收敛')
              x1 = zeros(n,1);     %构造初始值
              x2 = zeros(n, 1) + 1;
              n = 0;               %记录迭代次数
              if exist('r', 'var')              %根据是否有精确解来计算误差，如果没有则根据后验误差计算
                  while norm((r - x1), N) >= p   %判断收敛条件         
                     x2 = x1;
                     x1 = BJ*x2 + w*(D\b);
                     n = n + 1;
                  end
                  err = norm(x1 - r);     %计算误差值，统一采用矩阵二范数计算
              else
                  while norm((x2 - x1), N) >= p        
                     x2 = x1;
                     x1 = BJ*x2 + w*(D\b);
                     n = n + 1;
                  end
                  err = norm(x1 - x2);    
              end
              x = x1;
          end
      elseif isreal(eig(B)) == 0
          disp('请使用正交化对矩阵进行处理')
      elseif  isreal(eig(B)) && TB > 1                           %需要注意的是对于使用A进行迭代还是A'*A进行迭代迭代次数会有显著的变化
          disp('不存在最佳松弛因子，请输入松弛因子使JOR方法收敛')
          if all(eig(A)>0) && all(eig(2*D/wp-A)>0)
              disp('松弛因子选择合理，JOR方法收敛')
              BJ = eye(n) - wp*(D\A);  %计算无最佳松弛因子的JOR迭代中的迭代矩阵
              x1 = zeros(n,1);     %构造初始值
              x2 = zeros(n, 1) + 1;
              n = 0;               %记录迭代次数
              if exist('r', 'var')              %根据是否有精确解来计算误差，如果没有则根据后验误差计算
                  while norm((r - x1), N) >= p   %判断收敛条件         
                     x2 = x1;
                     x1 = BJ*x2 + wp*(D\b);
                     n = n + 1;
                  end
                  err = norm(x1 - r);     %计算误差值，统一采用矩阵二范数计算
              else
                  while norm((x2 - x1), N) >= p        
                     x2 = x1;
                     x1 = BJ*x2 + wp*(D\b);
                     n = n + 1;
                  end
                  err = norm(x1 - x2);    
              end
              x = x1;
          else
              disp('松弛因子选择不合理，JOR方法不收敛,请从新选择')
              return
          end
      end
  end
end
      
  
%A = [5 -1 -1 -1; -1 10 -1 -1; -1 -1 5 -1; -1 -1 -1 10], b = [-4 12 8 34]'
%A = [1 2 3 4 5; -2 3 4 5 6; -3 -4 5 6 7; -4 -5 -6 7 8; -5 -6 -7 -8 9], 
%b = [55 66 63 36 -25]', r =[1 2 3 4 5]', p = 0.00001, k = 2, N=Inf,wp=0.07


