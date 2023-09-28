function [x, n, err] = Gauss(A, b, p, k, N, wp, r)
%n 为迭代次数， x 为最终的迭代值 err为误差值
%A 为线性方程组，b为方程式右端值， p为收敛条件 
%k为迭代方法的选择，如果k=1则使用Gauss方法，若k = 2则使用SOR迭代法
%N为矩阵范数的选择，输入1，2，Inf分别代表1范数，2范数与无穷范数
%r为精确解，用来计算误差,如果没有精确值则不输入
%wp为如果没有最佳松弛因子时自己选择输入的松弛因子
% 初始解始终设为 x1 = 0
  U = triu(A, 1);     %提出上三角
  L = tril(A, -1);     %提出下三角
  d = diag(A);     %提取对角线
  D = diag(d);     %组成矩阵
  n = length(A);
  Bjacobi = -D\(L + U);    %求出Jacobi迭代矩阵B
  TB = max(abs(eig(Bjacobi)));   %求出Jacobi谱半径
  Bgauss = -(D + L)\U;       %求出gauss迭代矩阵B
  GB = max(abs(eig(Bgauss)));   %求出gauss谱半径
  if k == 1
      disp('使用Gauss方法')
      if GB < 1
         disp('Gauss方法收敛')
         x1 = zeros(n,1);     %构造初始值
         x2 = zeros(n, 1) + 0.5;
         n = 0;               %记录迭代次数
         if exist('r', 'var')              %根据是否有精确解来计算误差，如果没有则根据后验误差计算
            while norm((r - x1), N) >= p   %判断收敛条件
               x2 = x1;
               x1 = Bgauss*x2 + (D + L)\b;
               n = n + 1;
            end
            err = norm(x1 - r);     %计算误差值，统一采用矩阵二范数计算
         else
             while norm((x2 - x1), N) >= p   %判断收敛条件
               x2 = x1;
               x1 = Bgauss*x2 + (D + L)\b;
               n = n + 1;
            end
            err = norm(x1 - x2);     %计算误差值，统一采用矩阵二范数计算
         end
         x = x1;
      else
          disp('Gauss方法不收敛')
      end
  elseif k == 2
      disp('使用SOR方法')
      if isreal(eig(Bjacobi)) && TB < 1
          disp('SOR方法存在最佳松弛因子')                        
          w = 2/(1 + sqrt(1 - max(abs(eig(-inv(D)*(L+U))))^2));  %计算松弛因子
          BS = (D + w*L)\((1 - w)*D - w*U);  %计算SOR迭代中的B矩阵
          if max(abs(eig(BS))) < 1   %判断SOR方法是否收敛
              disp('SOR迭代法收敛')
              x1 = zeros(n,1);     %构造初始值
              x2 = zeros(n, 1) + 1;
              n = 0;               %记录迭代次数
              if exist('r', 'var')              %根据是否有精确解来计算误差，如果没有则根据后验误差计算
                  while norm((r - x1), N) >= p   %判断收敛条件         
                     x2 = x1;
                     x1 = BS*x2 + (D + w*L)\(w*b);
                     n = n + 1;
                  end
                  err = norm(x1 - r);     %计算误差值，统一采用矩阵二范数计算
              else
                  while norm((x2 - x1), N) >= p        
                     x2 = x1;
                     x1 = BS*x2 + (D + w*L)\(w*b);
                     n = n + 1;
                  end
                  err = norm(x1 - x2);    
              end
              x = x1;
          end
      elseif isreal(eig(Bjacobi)) == 0
          disp('请使用正交化对矩阵进行处理')
      elseif isreal(eig(Bjacobi)) &&  TB > 1                    %需要注意的是对于使用A进行迭代还是A'*A进行迭代迭代次数会有显著的变化
          disp('不存在最佳松弛因子，请输入松弛因子使SOR方法收敛')
          if wp > 0 && wp < 2
              disp('松弛因子选择合理，SOR方法收敛')
              BS = (D + wp*L)\((1 - wp)*D - wp*U);  %计算无最佳松弛因子的SOR迭代中的迭代矩阵
              x1 = zeros(n,1);     %构造初始值
              x2 = zeros(n, 1) + 1;
              n = 0;               %记录迭代次数
              if exist('r', 'var')              %根据是否有精确解来计算误差，如果没有则根据后验误差计算
                  while norm((r - x1), N) >= p   %判断收敛条件         
                     x2 = x1;
                     x1 = BS*x2 + (D + wp*L)\(wp*b);
                     n = n + 1;
                  end
                  err = norm(x1 - r);     %计算误差值，统一采用矩阵二范数计算
              else
                  while norm((x2 - x1), N) >= p        
                     x2 = x1;
                     x1 = BS*x2 + (D + wp*L)\(wp*b);
                     n = n + 1;
                  end
                  err = norm(x1 - x2);    
              end
              x = x1;
          else
              disp('松弛因子选择不合理，SOR方法不收敛,请从新选择')
              return
          end
      end
  end
end