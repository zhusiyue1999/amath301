clear; close all; clc;


% Problem 1
A = [4 2; 1 -2];
B = [1 -3; 3 -2];
C = [1 -3 -2; 1 7 4];
D = [-1 -2; 4 1; 2 -3];

x = [1; -1];
y = [2; 1];
z = [1; -2; 4];

A1 = A + B;

A2 = 2 * x - 3 * y;

A3 = C.';

A4 = C * z + y;

A5 = A * B;

A6 = B * A;

A7 = D * A;

% Problem 2
A = zeros(16);
for j = 1 : 16
    for k = 1 : 16
        A(j, k) = 1/(j + k - 1);
    end
end
disp(A)
x = [16:-1:1].'

A8 = A * x;

A9 = A(1:3, 13:16);

A10 = [A(:, 4), A(:, 8)];

%Problem 3
sum1 = 0;
for k = 1:8000
    sum1 = sum1 + 0.25;
end
x1 = abs(2000-sum1);
disp(x1);

sum2 = 0;
for k = 1:10000
    sum2 = sum2 + 0.2;
end
x2 = abs(2000-sum2);
disp(x2);

sum3 = 0;
for k = 1:16000
    sum3 = sum3 + 0.125;
end
x3 = abs(2000-sum3);
disp(x3);

sum4 = 0;
for k = 1:20000
    sum4 = sum4 + 0.1;
end
x4 = abs(2000-sum4);
disp(x4);

%Problem 4
p = [1:101];
p(1) = 0.8;

for j = 1:3
    subplot(3,1,j)
    r = 3 + 0.5 * (j-1);
    for k = 2:101
        p(k) = r * p(k-1)*(1-p(k-1));
    end
    time = [0:100];
    plot(time, p, 'ko-', 'Linewidth', [2])

    xlabel('time')
    ylabel('population')
    title(['r = ' num2str(r)])
end
print('HW1_fig1.png','-dpng')



% subplot(3,1,1)
% r1 = 3;
% for k = 2:101
%     p(k) = r1 * p(k-1)*(1-p(k-1));
% end
% disp(p);
% time = [0:100];
% plot(time, p, 'ko-', 'Linewidth', [2])
% 
% xlabel('time')
% ylabel('population')
% title('r = 3')
% 
% subplot(3,1,2)
% r2 = 3.5
% for k = 2:101
%     p(k) = r2 * p(k-1)*(1-p(k-1));
% end
% disp(p);
% time = [0:100];
% plot(time, p, 'ko-', 'Linewidth', [2])
% 
% xlabel('time')
% ylabel('population')
% title('r = 3.5')
% 
% subplot(3,1,3)
% r3 = 4
% for k = 2:101
%     p(k) = r3 * p(k-1)*(1-p(k-1));
% end
% disp(p);
% time = [0:100];
% plot(time, p, 'ko-', 'Linewidth', [2])
% 
% xlabel('time')
% ylabel('population')
% title('r = 4')

xx = 0:.1:1
yy= linspace(0,1,10)