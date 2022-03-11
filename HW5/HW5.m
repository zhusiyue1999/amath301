%% Problem1
clear; close all; clc

f = @(t) 1/3*exp(-t/2)+3*t*exp(-t/2);
f1 = @(t) -1/3*exp(-t/2)-3*t*exp(-t/2);
fprime = @(t) (-17/6)*exp(-t/2) + (3*exp(-t/2)*t)/2;
fprime2 = @(t) (-3*exp(-t/2)*t)/4 + (35/12)*exp(-t/2);
ans1 = [fminbnd(f1, 0, 10), f(fminbnd(f1, 0, 10))];

a = 0;
b = 10;
c = (-1+sqrt(5))/2;
x1 = c*a + (1-c)*b;
x2 = (1-c)*a + c*b;

f1 = -1/3*exp(-x1/2)-3*x1*exp(-x1/2);
f2 = -1/3*exp(-x2/2)-3*x2*exp(-x2/2);

while b-a >= 10^(-3)
    if f1 < f2
        b = x2;
        x2 = x1;
        f2 = f1;
        x1 = c*a + (1-c)*b;
        f1 = -1/3*exp(-x1/2)-3*x1*exp(-x1/2);
    else
        a = x1;
        x1 = x2;
        f1 = f2;
        x2 = (1-c)*a + c*b;
        f2 = -1/3*exp(-x2/2)-3*x2*exp(-x2/2);
    end
end
ans2 = [a, b];


tol = 10^(-3);
xk = 0;
maxIter = 10^5;
[ans3, ans4] = Newton(fprime, fprime2, xk, tol, maxIter);


%% Problem2
clear; close all; clc

f = @(p) (2-p(1))^2 + (p(2) - p(1)^2)^2;
ans5 = fminsearch(f, [0; 5]);

p = [0; 5]; % The current guess
fgrad = @(p) [-2*(2-p(1))-4*p(1)*(p(2)-p(1)^2); 2*(p(2)-p(1)^2)];
iteration = 0;
for k = 1:10000
    grad = fgrad(p); % Find which direction to go
    phi = @ (t) p - t*grad; % Define the "path"
    f_phi = @ (t) f(phi(t)); 
    tmin = fminbnd(f_phi,0,1);
    p = phi(tmin); 
    iteration = iteration + 1;
    if norm(grad, inf) < 10^(-4)
        break
    end
end
ans6 = p;
ans7 = iteration;

%% Problem3
clear; close all; clc

f = @(p) (2-p(1))^2 + (p(2) - p(1)^2)^2;
p = [0; 5]; % The current guess
fgrad = @(p) [-2*(2-p(1))-4*p(1)*(p(2)-p(1)^2); 2*(p(2)-p(1)^2)];
tstep = 0.01;
iteration = 0;
for k = 1:10000
    grad = fgrad(p); % Find which direction to go
    p = p - tstep*grad;
    iteration = iteration + 1;
    if norm(grad, inf) < 10^(-4)
        break
    end
end
ans8 = iteration;

f = @(p) (2-p(1))^2 + (p(2) - p(1)^2)^2;
p = [0; 5]; % The current guess
fgrad = @(p) [-2*(2-p(1))-4*p(1)*(p(2)-p(1)^2); 2*(p(2)-p(1)^2)];
tstep = 0.03;
iteration = 0;
for k = 1:10000
    grad = fgrad(p); % Find which direction to go
    p = p - tstep*grad;
    iteration = iteration + 1;
    if norm(grad, inf) < 10^(-4)
        break
    end
end
ans9 = iteration;

f = @(p) (2-p(1))^2 + (p(2) - p(1)^2)^2;
p = [0; 5]; % The current guess
fgrad = @(p) [-2*(2-p(1))-4*p(1)*(p(2)-p(1)^2); 2*(p(2)-p(1)^2)];
tstep = 0.05;
iteration = 0;
for k = 1:10000
    grad = fgrad(p); % Find which direction to go
    p = p - tstep*grad;
    iteration = iteration + 1;
    if norm(grad, inf) < 10^(-4)
        break
    end
end
ans10 = iteration;

%% Problem4
clear; close all; clc

f = @(p) (2-p(1))^2 + (p(2) - p(1)^2)^2;
ans5 = fminsearch(f, [0; 5]);

p = [0; 5]; % The current guess
fgrad = @(p) [-2*(2-p(1))-4*p(1)*(p(2)-p(1)^2); 2*(p(2)-p(1)^2)];
iteration = 0;
tic
for k = 1:10000
    grad = fgrad(p); % Find which direction to go
    phi = @ (t) p - t*grad; % Define the "path"
    f_phi = @ (t) f(phi(t)); 
    tmin = fminbnd(f_phi,0,1);
    p = phi(tmin); 
    iteration = iteration + 1;
    if norm(grad, inf) < 10^(-4)
        break
    end
end
time_fmin = toc
ans6 = p;
ans7 = iteration;


f = @(p) (2-p(1))^2 + (p(2) - p(1)^2)^2;
p = [0; 5]; % The current guess
fgrad = @(p) [-2*(2-p(1))-4*p(1)*(p(2)-p(1)^2); 2*(p(2)-p(1)^2)];
tstep = 0.01;
iteration = 0;
tic
for k = 1:10000
    grad = fgrad(p); % Find which direction to go
    p = p - tstep*grad;
    iteration = iteration + 1;
    if norm(grad, inf) < 10^(-4)
        break
    end
end
time_grad = toc
ans8 = iteration;
%%
function [xk,iteration] = Newton(f,fprime,xk,tol,maxIter)
% Performs Newton's method
% Inputs:
%   f = function to find the roots of
%   fprime = derivative of f
%   xk = initial guess
%   tol = tolerance for stopping criteria
%   maxIter = maximum number of iterations

change = 2*tol;
iteration = 0;
while change > tol && iteration < maxIter
    xkplus1 = xk - f(xk)/fprime(xk); 
    change = abs(xkplus1-xk);
    xk = xkplus1; 
    iteration = iteration+1;
end

end