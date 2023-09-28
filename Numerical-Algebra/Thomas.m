function x = Thomas(a, b, c, d)
  n = length(b);
  f(1) = c(1)/b(1);
  g(1) = d(1)/b(1);
  for i = 2:n-1     %由于f只到n-1，所以在这里范围到n-1， g（n）单独算
      h(i) = b(i) -f(i-1)*a(i - 1);
      f(i) = c(i)/h(i);
      g(i) = (d(i) - g(i-1)*a(i-1))/h(i);
  end
  g(n) = (d(n) - g(n-1)*a(n-1))/(b(n) - f(n-1)*a(n-1));
  x(n) = g(n);
  for i = n-1:-1:1
      x(i) = g(i) - f(i)*x(i+1);
  end
end

