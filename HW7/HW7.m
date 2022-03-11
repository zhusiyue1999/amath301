%% Problem 1
clear; close all; clc

load('particle_position.mat')

A = A - mean(A,2);

[U,S,V] = svd(A, 'econ');
ans1 = diag(S);
Arank1 = U(:,1:1)*S(1:1,1:1)*V(:,1:1)';
ans2 = norm(A-Arank1);

Arank2 = U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
ans3 = norm(A-Arank2);

%% Problem 2

plot3( A(1,:),A(2,:),A(3,:), 'k.');
hold on
plot3( Arank1(1,:),Arank1(2,:),Arank1(3,:), 'r.');
axis vis3d
xlabel('x', 'fontsize', [20]);
ylabel('y', 'fontsize', [20]);
zlabel('z', 'fontsize', [20]);
legend('Positions of the particle', 'Rank-1 approximation', 'fontsize', [15]);
print('HW7_fig1.png','-dpng')

%%
plot3( A(1,:),A(2,:),A(3,:), 'k.');
hold on
plot3( Arank2(1,:),Arank2(2,:),Arank2(3,:), 'r.');
axis vis3d
xlabel('x', 'fontsize', [20]);
ylabel('y', 'fontsize', [20]);
zlabel('z', 'fontsize', [20]);
legend('Positions of the particle', 'Rank-2 approximation', 'fontsize', [15]);
print('HW7_fig2.png','-dpng')


%% Problem3
clear; close all; clc

A = imread('llama.jpg');
A = im2double(rgb2gray(A));

[U,S,V] = svd(A, 'econ');
ans4 = diag(S(1:20, 1:20));

ans5 = S(1,1)/sum(sum(S));

ans6 = sum(ans4)/sum(sum(S));

count = 0;
sum = 0;
diagonal = diag(S);
sum_diagonal = 0;
for k = 1:876
    sum_diagonal = sum_diagonal + diagonal(k,:);
end

while sum / sum_diagonal < 0.9
    count = count + 1;
    sum = sum + diagonal(count,:);
end
ans7 = count;

%% Problem 4
[U,S,V] = svd(A, 'econ');
Arank1 = U(:,1:1)*S(1:1,1:1)*V(:,1:1)';
Arank20 = U(:,1:20)*S(1:20,1:20)*V(:,1:20)';
Arank206 = U(:,1:206)*S(1:206,1:206)*V(:,1:206)';

subplot(2,2,1)
imshow(A)
subplot(2,2,2)
imshow(Arank1)
subplot(2,2,3)
imshow(Arank20)
subplot(2,2,4)
imshow(Arank206)

print('HW7_fig3.png','-dpng')

%% Problem 5
clear; close all; clc

load('NoisyImage.mat')

[U,S,V] = svd(A_noise, 'econ');

diagonal = diag(S);
sum_diagonal = 0;
for k = 1:length(S)
    sum_diagonal = sum_diagonal + diagonal(k,:);
end

ans8 = (S(1,1) + S(2,2))/sum_diagonal;

A_noise_rank2 = U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
ans9 = norm(A_noise-A_noise_rank2);
ans10 = norm(A-A_noise_rank2);

%% Problem6
 
semilogy(diag(S), 'bo')
 
% energies = cumsum(diagonal) / sum_diagonal
% semilogy(energies)
print('HW7_fig4.png','-dpng')
 
subplot(1,3,1)
imshow(A)
title('True Image', 'fontsize', [20])
subplot(1,3,2)
imshow(A_noise)
title('Noisy Image', 'fontsize', [20])
subplot(1,3,3)
imshow(A_noise_rank2)
title('Rank-2 Approx', 'fontsize', [20])
print('HW7_fig5.png','-dpng')

