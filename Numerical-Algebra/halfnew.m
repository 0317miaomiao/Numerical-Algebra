function [c, k] = halfnew(a, b, tol)
%c最终近似解，k迭代次数
%a左边值，b右边值， tol误差值
c = (a + b)/2;
k = 1;
m = 1 + round((log(b - a) - log(2*tol))/log(2));
function y = f(x)
   y = x^3 - 48;
end
while k <= m
    if f(c) == 0
        c = c;
    return;
    elseif f(a)*f(c) < 0
        b = (a + b)/2;
    else
        a = (a + b)/2;
    end
    c = (a + b)/2;
    k = k+1;
end
end
