clear; close all; clc;


% % Problem 1
% left = 0;
% right = 1.2;
% tolerance = 1e-9;  
% max_iter = 1000;
% f = @ (x) x.^5-2*x.^4+4*x.^3-2*x.^2+x-1-cos(30*x);
% 
% f_mid = 1;
% iterations = 0; 
% while abs(f_mid) > tolerance && iterations < max_iter
%     
%     root = (left+right)/2;
%     f_mid = f(root);
%     
%     if f_mid*f(left) < 0
%         right = root; 
%     elseif f_mid*f(right) < 0
%         left = root; 
%     else 
%         disp('No root found') 
%         break 
%     end
%     
%     iterations = iterations + 1; 
%     
% end
% 
% 
% % Problem 2
% f = @ (x) x.^2-3;
% f1 = @ (x) 2*x;
% x = 10;
% for k = 1:50
%     x = x - f(x)/f1(x);
% end
% error = abs(sqrt(3) - x);
% 
% 
% tolerance = 1e-7;
% iterations = 0;
% y = 10;
% guesses = [10];
% error_vector = [abs(sqrt(3) - 10)];
% while iterations <= 50 && abs(f(y)) > tolerance
%     y = y - f(y)/f1(y);
%     iterations = iterations + 1;
%     guesses(iterations + 1) = y;
%     error_vector (iterations + 1) = abs(sqrt(3) - y);
% end
% 
% Problem 3
% fun = @ (x) 9-x-0.5*sin(2*x);
% fun1 = @ (x) -cos(2*x)-1;
% exact = fzero(fun,0);
% 
% iterations = 0;
% tolerance = 1e-9;
% x = 0;
% while iterations <100 && abs(fun(x)) > tolerance
%     x = x - fun(x)/fun1(x);
%     iterations = iterations + 1;
% end
% Newton1 = [abs(exact-x), iterations];
% 
% iterations = 0;
% tolerance = 1e-9;
% x = 2;
% while iterations <100 && abs(fun(x)) > tolerance
%     x = x - fun(x)/fun1(x);
%     iterations = iterations + 1;
% end
% Newton2 = [abs(exact-x), iterations];
% 
% iterations = 0;
% tolerance = 1e-9;
% x = 50;
% while iterations <100 && abs(fun(x)) > tolerance
%     x = x - fun(x)/fun1(x);
%     iterations = iterations + 1;
% end
% Newton3 = [abs(exact-x), iterations];
%
% Problem 4
f = @ (x) x^3-x^2+2*x+3;
f1 = @ (x) 3*x^2-2*x+2;
exact = fzero(f,0);

left = -100;
right = 100;
tolerance = 10 ^ (-8);  

f_mid = 1;
iterations = 0; 
errors1 = [];
n1 = [];
while abs(f_mid) > tolerance
    
    mid = (left+right)/2;
    f_mid = f(mid);
    
    if f_mid*f(left) < 0
        right = mid; 
    elseif f_mid*f(right) < 0
        left = mid; 
    else 
        disp('No root found') 
        break 
    end
    errors1(iterations + 1) = abs(exact - mid);
    n1(iterations + 1) = iterations + 1;
    iterations = iterations + 1; 
end

iterations = 0;
tolerance = 10 ^ (-8);
x = 100;
errors2 = [abs(exact - x)];
n2 = [0];
while iterations <100 && abs(f(x)) > tolerance
    x = x - f(x)/f1(x);
    
    errors2(iterations + 2) = abs(exact - x);
    n2(iterations + 2) = iterations + 1;
    iterations = iterations + 1;
end

plot(n1, errors1, 'ko'), hold on
plot(n2, errors2, 'r*')
set(gca, 'fontsize', [15])
legend('bisection method', 'Newton’s method', 'location', 'best', 'fontsize', [15])
print('HW2_fig1.png','-dpng')

% Problem 5
semilogy(n1, errors1, 'ko'), hold on
semilogy(n2, errors2, 'r*')
legend('bisection method', 'Newton’s method', 'location', 'best', 'fontsize', [15])
print('HW2_fig2.png','-dpng')

% Problem 6
f = @ (x) x^5-3*x^4+5*x^3-7*x^2+6*x-2;
f1 = @ (x) 5*x^4-12*x^3+15*x^2-14*x+6;
exact = fzero(f,0);

left = -100;
right = 100;
tolerance = 10 ^ (-8);  

f_mid = 1;
iterations = 0; 
errors1 = [];
n1 = [];
while abs(f_mid) > tolerance
    
    mid = (left+right)/2;
    f_mid = f(mid);
    
    if f_mid*f(left) < 0
        right = mid; 
    elseif f_mid*f(right) < 0
        left = mid; 
    else 
        disp('No root found') 
        break 
    end
    errors1(iterations + 1) = abs(exact - mid);
    n1(iterations + 1) = iterations + 1;
    iterations = iterations + 1; 
end

iterations = 0;
tolerance = 10 ^ (-8);
x = 100;
errors2 = [abs(exact - x)];
n2 = [0];
while iterations <100 && abs(f(x)) > tolerance
    x = x - f(x)/f1(x);
    
    errors2(iterations + 2) = abs(exact - x);
    n2(iterations + 2) = iterations + 1;
    iterations = iterations + 1;
end

semilogy(n1, errors1, 'ko'), hold on
semilogy(n2, errors2, 'r*')
legend('bisection method', 'Newton’s method', 'location', 'best', 'fontsize', [15])
print('HW2_fig3.png','-dpng')