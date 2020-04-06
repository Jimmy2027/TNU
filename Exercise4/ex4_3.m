%% 4.3
t = linspace(0,0.2,200);
ka1 = 80;
ka2 = 50;
f1 = 50;
f2 = 50;
af = 3000;
ab = 1000;
c = 1;
mu = 0.05;
sig = 0.01;

u = normpdf(t,mu,sig);
A = [0,1,0,0;-ka1^2,-f1,af,0;0,0,0,1;ab,0,-ka2^2,-f2];
C = [0,0,0,c]';
x = eulerintegrate(A,C,u);
load('tn20_ex4.mat')
x_true = x_condition_1;

figure;
plot(x_true(1,:))
hold on;
plot(x_true(2,:))
plot(x_true(3,:))
plot(x_true(4,:))

figure;
plot(x(1,:))
hold on;
plot(x(2,:))
plot(x(3,:))
plot(x(4,:))
%% 4.3

t = linspace(0,0.2,200);
ka1 = 80;
ka2 = 50;
f1 = 50;
f2 = 50;
af = 3000;
ab = 1000;
c = 1;
mu = 0.05;
sig = 0.01;
u = normpdf(t,mu,sig);
number_of_tries = 50;

%% Testing model 1 where ka1 was changed
explained_variances1 = zeros(1,number_of_tries);
i=1;
searchspace = linspace(40,120,number_of_tries);
for ka1 = linspace(40,120,number_of_tries)
A = [0,1,0,0;-ka1^2,-f1,af,0;0,0,0,1;ab,0,-ka2^2,-f2];
C = [0,0,0,c]';
x = eulerintegrate(A,C,u);
vE = 1 - var(x_condition_2-x)/var(x_condition_2);
explained_variances1(1,i) = vE;
i = i+1;
end
[m,i] = max(explained_variances);
max1 = m
best_param1 = searchspace(i)

%% Testing model 2, where ka2 was changed
explained_variances2 = zeros(1,number_of_tries);
i=1;
searchspace = linspace(10,90,number_of_tries);
for ka2 = linspace(10,90,number_of_tries)
A = [0,1,0,0;-ka1^2,-f1,af,0;0,0,0,1;ab,0,-ka2^2,-f2];
C = [0,0,0,c]';
x = eulerintegrate(A,C,u);
vE = 1 - var(x_condition_2-x)/var(x_condition_2);
explained_variances2(1,i) = vE;
i = i+1;
end
[m,i] = max(explained_variances2);
max2 = m
best_param2 = searchspace(i)
%% Testing model 3, where af was changed
explained_variances3 = zeros(1,number_of_tries);
i=1;
searchspace = linspace(2800,3200,number_of_tries);
for af = linspace(3000-3000/number_of_tries,3000+3000/number_of_tries,number_of_tries)
A = [0,1,0,0;-ka1^2,-f1,af,0;0,0,0,1;ab,0,-ka2^2,-f2];
C = [0,0,0,c]';
x = eulerintegrate(A,C,u);
vE = 1 - var(x_condition_2-x)/var(x_condition_2);
explained_variances3(1,i) = vE;
i = i+1;
end
[m,i] = max(explained_variances3);
max3 = m
best_param3 = searchspace(i)
%% Testin model 4, where ab was changed
explained_variances4 = zeros(1,number_of_tries);
i=1;
searchspace = linspace(900,1100,number_of_tries);
for ab = linspace(900,1100,number_of_tries)
A = [0,1,0,0;-ka1^2,-f1,af,0;0,0,0,1;ab,0,-ka2^2,-f2];
C = [0,0,0,c]';
x = eulerintegrate(A,C,u);
vE = 1 - var(x_condition_2-x)/var(x_condition_2);
explained_variances4(1,i) = vE;
i = i+1;
end
[m,i] = max(explained_variances4);
max4 = m
best_param4 = searchspace(i)

best_params = [best_param1,best_param2,best_param3,best_param4]
[m,i] = max([max1,max2, max3, max4])
max_all = m
fprintf('best model is model%i with best param = %i',i,best_params(i))


