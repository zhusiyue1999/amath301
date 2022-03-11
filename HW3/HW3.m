%% Problem1
clear; close all; clc

% x = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13;f14;f15];
A = [-1/2,1,0,0,0,0,0,0,0,1/2,0,0,0,0,0;
     -sqrt(3)/2,0,0,0,0,0,0,0,0,-sqrt(3)/2,0,0,0,0,0;
     0,-1,1,0,0,0,0,0,0,0,-1/2,1/2,0,0,0;
     0,0,0,0,0,0,0,0,0,0,-sqrt(3)/2,-sqrt(3)/2,0,0,0;
     0,0,-1,1,0,0,0,0,0,0,0,0,-1/2,1/2,0;
     0,0,0,0,0,0,0,0,0,0,0,0,-sqrt(3)/2,-sqrt(3)/2,0;
     0,0,0,-1,1/2,0,0,0,0,0,0,0,0,0,-1/2;
     0,0,0,0,-sqrt(3)/2,0,0,0,0,0,0,0,0,0,-sqrt(3)/2;
     0,0,0,0,-1/2,-1,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,1,-1,0,0,0,0,0,0,-1/2,1/2;
     0,0,0,0,0,0,0,0,0,0,0,0,0,sqrt(3)/2,sqrt(3)/2;
     0,0,0,0,0,0,1,-1,0,0,0,-1/2,1/2,0,0;
     0,0,0,0,0,0,0,0,0,0,0,sqrt(3)/2,sqrt(3)/2,0,0;
     0,0,0,0,0,0,0,1,-1,-1/2,1/2,0,0,0,0;
     0,0,0,0,0,0,0,0,0,sqrt(3)/2,sqrt(3)/2,0,0,0,0];
 b = [0;0;0;0;0;0;0;0;0;0;10000;0;11500;0;15000];
 forces = A\b;
 max_force = 0;
 max_force = max(abs(forces));
 
 x = A\b;
 while max(abs(x)) < 30000
     b(11) = b(11) + 100;
     x = A\b;
 end
 W7 = b(11);
 collapsed_beam = 3;
%% Problem2
clear; close all; clc

A = [75,-30,-20;
     -30,50,-15;
     -20,-15,45];
[L,U,P] = lu(A);
ULP = U*L*P;
b = [20;30;15];
y = L\(P*b);
currents = U\y;

b = [50;10;1];
I2_vec = [];
for k = 1:100
    y = L\b;
    x = U\y;
    I2_vec(k) = x(2);
    b(3) = b(3) + 1;
end
I2_vec

%% Problem3
clear; close all; clc

A = hilb(11);
v = [1;1;1;1;1;1;1;1;1;1;1];
b = A*v;
x = A\b;
e = v - x;
r = A*x - b;
backslash_result = [norm(e) norm(r)];

[L,U,P] = lu(A);
y = L\(P*b);
x = U\y;
e = v - x;
r = A*x - b;
LU_result = [norm(e) norm(r)];

x = inv(A)*b;
e = v - x;
r = A*x - b;
inv_result = [norm(e) norm(r)];

%% Problem4
clear; close all; clc

nvec = 1000:500:4000;

for n = nvec
    A = rand(n);
    b = rand(n,1);
    
    tic
    x = A\b;
    bstime = toc; 
    loglog(n,bstime, 'ob', 'markersize', 10)
    hold on
    
    tic
    [L,U,P] = lu(A);
    y = L\(P*b);
    x = U\y;
    lutime = toc; 
    loglog(n,lutime, 'r+', 'markersize', 10)
    hold on
    
    tic
    x = inv(A)*b;
    intime = toc; 
    loglog(n,intime, 'kd', 'markersize', 10)
    hold on
end
xlabel('n')
ylabel('time')
set(gca, 'fontsize', 14)
loglog(nvec,3e-8*nvec.^2,'m') 
loglog(nvec,3e-11*nvec.^3,'g')
print('HW3_fig1.png','-dpng')

%% Problem5

clear; close all; clc
A = toeplitz(1:100);
b = rand(100,1);
[L,U,P] = lu(A);
tic
for k = 1:1000
    x = A\b;
end
backslashtime1 = toc; 
tic
for k = 1:10000
    x = A\b;
end
backslashtime2 = toc; 
tic
for k = 1:1000
    y = L\(P*b);
    x = U\y;
end
lutime1 = toc; 
tic
for k = 1:10000
    y = L\(P*b);
    x = U\y; 
end
lutime2 = toc;
tic
for k = 1:1000
    x = inv(A)*b;
end
intime1 = toc; 
tic
for k = 1:10000
    x = inv(A)*b;
end
intime2 = toc; 
backslashtime1
lutime1
intime1
backslashtime2
lutime2
intime2