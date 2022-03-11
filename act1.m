clear; clc

% coefficients of quadratic polynomial p(x) = a*x^2 + b*x + c
load('coefficients.mat');

discriminant = b^2 - 4*a*c;
root1 = (-b - sqrt(discriminant))/(2*a)
root2 = (-b + sqrt(discriminant))/(2*a)

save('roots','root1','root2')
