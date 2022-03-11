%% Problem1
clear; close all; clc

a1 = [1,80];
for k = 1:80
    a1(k) = 2;
end
a2 = [1,79];
for k = 1:79
    a2(k) = -1;
end
A = diag(a1) + diag(a2,1) + diag(a2,-1);

b = [80;1];
for k = 1:80
    b(k) = exp(-15*pi/10)*sin(15*pi*k/81);
end


ans1 = A\b;

D = diag(diag(A));
T = A - D;
M = -D\T;
lambda = eig(M);
ans2 = max(abs(lambda));


x0 = 2*ones(80,1);
tol = 1e-4;
maxIter = 10000000;
[x_Jac, num_Jac] = Jacobi(A,b,x0,tol,maxIter);
ans3 = [num_Jac, norm(ans1 - x_Jac)];


S = tril(A);
T = A - S;
M = -S\T;
lambda = eig(M);
ans4 = max(abs(lambda));


[x_Gau, num_Gau] = GaussSeidel(A,b,x0,tol,maxIter);
ans5 = [num_Gau, norm(ans1 - x_Gau)];


%% Problem2 
clear; close all; clc

a1 = [1,80];
for k = 1:80
    a1(k) = 2;
end
a2 = [1,79];
for k = 1:79
    a2(k) = -1;
end
A = diag(a1) + diag(a2,1) + diag(a2,-1);

b = [80;1];
for k = 1:80
    b(k) = exp(-15*pi/10)*sin(15*pi*k/81);
end

D = diag(diag(A));
L = tril(A, -1);
U = triu(A, 1);

w = 1.2;
P = 1/w*D + L;
T = (w - 1)/w * D + U;
M = -P\T;
lambda = eig(M);
ans6 = max(abs(lambda));

x0 = 2*ones(80,1);
tol = 1e-4;
maxIter = 10000000;
[x_SOR, num_SOR] = SOR(A,b,x0,tol,maxIter,w);
ans7 = [num_SOR, norm(A\b - x_SOR)];

max_lambda = 1000;
best_w = 1;
for k = 1:0.01:2
    w = k;
    P = 1/w*D + L;
    T = (w - 1)/w * D + U;
    M = -P\T;
    lambda = max(abs(eig(M)));
    if max_lambda >= lambda
        max_lambda = lambda;
        best_w = k;
    end
end
ans8 = best_w;
ans9 = max_lambda;

[x_SOR, num_SOR] = SOR(A,b,x0,tol,maxIter,best_w);
ans10 = [num_SOR, norm(A\b - x_SOR)];

%% Problem4
clear; close all; clc

a1 = [1,80];
for k = 1:80
    a1(k) = 2;
end
a2 = [1,79];
for k = 1:79
    a2(k) = -1;
end
A = diag(a1) + diag(a2,1) + diag(a2,-1);

b = [80;1];
for k = 1:80
    b(k) = exp(-15*pi/10)*sin(15*pi*k/81);
end

x_exact = A\b;

x0 = 2*ones(80,1);
tol = 1e-4;
maxIter = 10000000;

D = diag(diag(A));
T = A - D;
M = -D\T;
change = 2*tol;
c = D\b;
xk = x0; 
iteration = 0;

error_vector = [norm(x_exact - x0)];
k = [0];

while change > tol && iteration < maxIter
    xkplus1 = M*xk + c;
    change = norm(xkplus1-xk,inf);
    iteration = iteration + 1;
    k(iteration + 1) = iteration;
    error_vector(iteration + 1) = norm(x_exact - xkplus1);
    xk = xkplus1;
end
semilogy(k, error_vector, 'r', 'LineWidth', 2)
hold on

lambda = max(abs(eig(M)));
lambda_vector = [];
for j = 1:length(k)
    lambda_vector(j) = lambda^(j-1);
end
semilogy(k, lambda_vector, 'k--', 'LineWidth', 2)
title('Jacobi Method Error', 'fontsize', [20])
xlabel('Iterations, k', 'fontsize', [15])
ylabel('Error (log scale)', 'fontsize', [15])
legend('Jacobi Error', 'λ^k_m_a_x', 'location', 'best', 'fontsize', [15])
print('HW4_fig1.png','-dpng')


%%
function [x, numIter] = Jacobi(A,b,x0,tol,maxIter)
D = diag(diag(A));
T = A - D;
M = -D\T;
change = 2*tol;
c = D\b;
xk = x0; 
iteration = 0;
while change > tol && iteration < maxIter
    xkplus1 = M*xk + c;
    change = norm(xkplus1-xk,inf);
    iteration = iteration + 1;
    xk = xkplus1;
end
x = xk;
numIter = iteration;
end

function [x, numIter] = GaussSeidel(A,b,x0,tol,maxIter)
S = tril(A);
T = A - S;
M = -S\T;
change = 2*tol;
c = S\b;
xk = x0; 
iteration = 0;
while change > tol && iteration < maxIter
    xkplus1 = M*xk + c;
    change = norm(xkplus1-xk,inf);
    iteration = iteration + 1;
    xk = xkplus1;
end
x = xk;
numIter = iteration;
end

function [x, numIter] = SOR(A,b,x0,tol,maxIter,w)
L = tril(A, -1);
U = triu(A, 1);
P = 1/w*(diag(diag(A))) + L;
T = (w - 1)/w * (diag(diag(A))) + U;
M = -P\T;
change = 2*tol;
c = P\b;
xk = x0; 
iteration = 0;
while change > tol && iteration < maxIter
    xkplus1 = M*xk + c;
    change = norm(xkplus1-xk,inf);
    iteration = iteration + 1;
    xk = xkplus1;
    error = xkplus1-xk;
end
x = xk;
numIter = iteration;
end