%% Problem2
clear; close all; clc

load('SeaPopData.mat')
plot(t, Seattle_Pop, 'ko', 'LineWidth', 2)
hold on

lin_fit = polyfit(t, Seattle_Pop, 1);
ans1 = lin_fit(1);
ans2 = lin_fit(1) * 159 + lin_fit(2);
x_plot = 0:0.1:160;
y_plot = polyval(lin_fit, x_plot);
plot(x_plot, y_plot, 'b', 'LineWidth', 2)
hold on


poly_fit = polyfit(t, Seattle_Pop, 2);
y_plot = polyval(poly_fit, x_plot);
plot(x_plot, y_plot, 'r', 'LineWidth', 2)
hold on
ans3 = poly_fit(1) * 159^2 + poly_fit(2) * 159 + poly_fit(3);

poly_fit5 = polyfit(t, Seattle_Pop, 5);
y_plot = polyval(poly_fit5, x_plot);
plot(x_plot, y_plot, 'm', 'LineWidth', 2);
hold on

sum1 = 0;
for k = 1:length(poly_fit5)
    sum1 = sum1 + poly_fit5(k) * 159^(6-k);
end
sum1;

poly_fit9 = polyfit(t, Seattle_Pop, 9);
y_plot = polyval(poly_fit9, x_plot);
plot(x_plot, y_plot, 'g', 'LineWidth', 2)
hold on
axis([0 160 0 800000])
xlabel('Years since 1860', 'fontsize', [20])
ylabel('Seattle Population', 'fontsize', [20])
legend('data', 'deg 1', 'deg 2', 'deg 5', 'deg 9', 'fontsize', [20], 'location', 'northwest')
print('HW6_fig1.png','-dpng')

sum2 = 0;
for k = 1:length(poly_fit9)
    sum2 = sum2 + poly_fit9(k) * 159^(10-k);
end
sum2;

ans4 = [sum1, sum2];

%% Problem4
clear; close all; clc

load('CO2_data.mat')

% lin_fit = polyfit(t, log(CO2), 1);
lin_fit = polyfit(t, (CO2), 1);
ans5 = lin_fit(1);
ans6 = exp(lin_fit(2));

x_plot = 0:0.1:500;
y_plot = polyval(lin_fit, x_plot);

% plot(t, log(CO2), '-k.', 'LineWidth', 2)
plot(t, (CO2), '-k.', 'LineWidth', 2)

hold on
plot(x_plot, y_plot, 'r', 'LineWidth', 2)
xlim([0 65])
xlabel('Years since 1958', 'fontsize', [20])
ylabel('Atmospheric CO2', 'fontsize', [20])
legend('data', 'fit curve', 'fontsize', [20], 'location', 'northwest')
print('HW6_fig2.png','-dpng')

%% Problem6
clear; close all; clc

load('CO2_data.mat')

a = 30;
r = 0.03;
b = 300;

% error = exponential_shift(a, r, b);
guess = [30, 0.03, 300];
[best, ans8] = fminsearch(@(arb) exponential_shift(arb), guess);
ans7 = best.';


% plot(t, log(CO2), '-k.', 'LineWidth', 2)
plot(t, (CO2), '-k.', 'LineWidth', 2)
hold on
x_plot = 0:0.1:70;
% y_plot = log(best(1)*exp(best(2) * x_plot) + best(3));
y_plot = (best(1)*exp(best(2) * x_plot) + best(3));
plot(x_plot, y_plot, 'r', 'LineWidth', 2)
xlim([0 65])
xlabel('Years since 1958', 'fontsize', [20])
ylabel('Atmospheric CO2', 'fontsize', [20])
legend('data', 'fit curve', 'fontsize', [20], 'location', 'northwest')
print('HW6_fig3.png','-dpng')


%  Problem8

guess2 = [ans7(1), ans7(2), ans7(3), 3, 6, 0];
[best2, ans10] = fminsearch(@(arbcde) exponential_shift_sin(arbcde), guess2);
ans9 = best2.';


% plot(t, log(CO2), '-k.', 'LineWidth', 2)
plot(t, (CO2), '-k.', 'LineWidth', 2)
hold on
x_plot = 0:0.1:70;
% y_plot = log(best2(1)*exp(best2(2) * x_plot) + best2(3) + best2(4)*sin(best2(5)*(x_plot - best2(6))));
y_plot = (best2(1)*exp(best2(2) * x_plot) + best2(3) + best2(4)*sin(best2(5)*(x_plot - best2(6))));
plot(x_plot, y_plot, 'r', 'LineWidth', 2)
xlim([0 65])
xlabel('Years since 1958', 'fontsize', [20])
ylabel('Atmospheric CO2', 'fontsize', [20])
legend('data', 'fit curve', 'fontsize', [20], 'location', 'northwest')
print('HW6_fig4.png','-dpng')

function [sum2, error2] = exponential_shift_sin(arbcde)
a = arbcde(1);
r = arbcde(2);
b = arbcde(3);
c = arbcde(4);
d = arbcde(5);
e = arbcde(6);
load('CO2_data.mat')
sum2 = 0;
for k = 1:length(t)
    sum2 = sum2 + (a*exp(r*t(k)) + b + c*sin(d*(t(k)-e)) - CO2(k))^2;
end 
sum2;
error2 = sqrt(1/746*sum2);
end
%%


function [sum, error] = exponential_shift(arb)
a = arb(1);
r = arb(2);
b = arb(3);
load('CO2_data.mat')
sum = 0;
for k = 1:length(t)
    sum = sum + (a*exp(r*t(k)) + b - CO2(k))^2;
end
sum;
error = sqrt(1/746*sum);
end
