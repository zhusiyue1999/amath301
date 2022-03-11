%% Problem 1
clear; close all; clc

load('Plutonium.mat')

ans1 = (P(12) - P(10)) / (2*1);

ans2 = (-3*P(1) + 4*P(2) - P(3)) / (2*1);

ans3 = (3*P(41) - 4*P(40) + P(39)) / (2*1);

P_prime = zeros([1 41]);
P_prime(1) = ans2;
P_prime(end) = ans3;
P_prime(2:end-1) = (P(3:end)-P(1:end-2))/(2*1);
ans4 = zeros([1 41]);
ans4(1:end) = -1 * P_prime / P(1:end);


ans5 = log(2) / (sum(ans4) / 41);

%% Problem 2
clear; close all; clc

load('BloodFlow.mat')

ans6 = 0.001/2*(v(1) .* r(1) + 2*sum(v(2:end-1) .* r(2:end-1)) + v(end) .* r(end)) * 2 * pi;

ans7 = 0.001 * 2 * pi * sum(r(2:end));

ans8 = ans6 / ans7;

%% Problem 3
clear; close all; clc

load('Dye.mat')

dx = 0.1;
ans9 = dx/3 * (c(1) + 4*sum(c(2:2:end-1)) + 2* sum(c(3:2:end-2)) + c(end));

ans10 = 3 / ans9;

%% Problem 4
clear; close all; clc

f = @(x) exp(-(x-85).^2/50) / sqrt(50*pi);
exact = integral(f, 76, 86);

n = 0:-1:-16;
dx = 2.^n;
x = 76;
error_left = zeros([1 17]);
error_right = zeros([1 17]);
error_midpoint = zeros([1 17]);
error_trap = zeros([1 17]);
error_simpson = zeros([1 17]);

midpoint = zeros([1 17]);
for k = 1:length(dx)
    x = 76:dx(k):86;
    
    left = dx(k) * sum(f(x(1:end-1)));
    error_left(k) = abs(exact - left);
    
    right = dx(k) * sum(f(x(2:end)));
    error_right(k) = abs(exact - right);
    
    mid = f((x(1:end-1) + x(2:end)) / 2);
    midpoint(k) = dx(k) * sum(mid);
    error_midpoint = abs(exact - midpoint);
    
    trap(k) = dx(k)/2*(f(x(1)) + 2*sum(f(x(2:end-1))) + f(x(end)));
    error_trap = abs(exact - trap);
    
    simpson(k) = dx(k)/3 * (f(x(1)) + 4*sum(f(x(2:2:end-1))) + 2*sum(f(x(3:2:end-2))) + f(x(end)));
    error_simpson = abs(exact - simpson);
end


loglog(dx, error_left, 'ko','MarkerSize',12)
hold on
loglog(dx, error_right, 'ro', 'LineWidth', 3)
loglog(dx, error_midpoint, 'cd', 'LineWidth', 3)
loglog(dx, error_trap, 'b+', 'LineWidth', 3)
loglog(dx, error_simpson, 'mo', 'LineWidth', 3)

loglog(dx, 0.03*dx, 'k-', 'LineWidth', 3) 
loglog(dx, 0.0005*dx.^2, 'r-.', 'LineWidth', 3)
loglog(dx, 0.000002*dx.^4, 'k--', 'LineWidth', 3)

line = 10^-16 + zeros([1 length(dx)]);
loglog(dx, line, 'r:', 'LineWidth', 3)

title('Convergence of Numerical Integration Schemes', 'fontsize', [20])
xlabel('Grid Spacing, h, k', 'fontsize', [15])
ylabel('Error', 'fontsize', [15])
legend('Left Rectangle', 'Right Rectangle', 'Midpoint', 'Trapezoidal', 'Simpson’s', 'O(h)', 'O(h^2)', 'O(h^4)', 'Machine Precision', 'fontsize', [15], 'Location', 'Northeastoutside')
ylim([10^-18 1])
print('HW8_fig1.png','-dpng')
