clear all;
close all;
%% Model 2 with parameters from 3.1
% neural model construction

%construct struct P of matrices A and B and vector C
P.A = [-0.5,0;1,-0.5];
P.B = [0,0;-0.5,0];
P.C = [1;0];

%define x0 at t = 0
x0 = [0;0];

%construct u
u_vector = zeros(2,800);
u_vector(2,301:601)= 1;
u_vector(1,70:70:631)=5;

% hrf model construction

%hemodynamic state vector at t=0 (s,f,v,q)
h0 = [0;1;1;1];

% prameters for hrf : kappa, gamma, tau, alpha and E_0
Phrf=[0.64,0.32,2,0.32,0.4];

% compute dcm
t = linspace(0,80,800);
[y,h,x] = euler_integrate_dcm(u_vector,P,Phrf,x0,h0);

%% a) generate noisy BOLD trace \hat{y}
sigma = 0.005;
y_hat = zeros(2,800);
for i = 1:800
    noise = normrnd(0,sigma);
    y_hat(:,i) = y(:,i) + noise;
end

% plot trace
figure;
plot(t,y(:,:));
hold on;
plot(t,y_hat(:,:));
legend('$y_{x1}$','$y_{x2}$','$\hat{y}_{x1}$','$\hat{y}_{x2}$','Interpreter','Latex')
title('noisy BOLD signal');

%% b) compute log-likelihood 

llh_y1 = -0.5 * log(2 *pi * sigma^2)-0.5*(y_hat(1,:)-y(1,:))*(y_hat(1,:)-y(1,:))'/(sigma^2);
llh_y2 = -0.5 * log(2 *pi * sigma^2)-0.5*(y_hat(2,:)-y(2,:))*(y_hat(2,:)-y(2,:))'/(sigma^2);

%% c) compute log joint distribution



