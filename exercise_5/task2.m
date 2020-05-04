clear all;
close all;

%% a) 10 samples from gaussian y= mu+eps, mu=tau=1
%%
% *p(eps)= 1/sqrt(2pi)* exp(-eps/2)*

obs = zeros(10,1);

for i = 1:10
    obs(i)= 1 + normrnd(0,1);
end


%% b) update equations 
%%
% *mu0=0,lambda0=3, a0=2,b0=2* 
mu =1;
tau =1;
mu0 = 0;
lambda0=3;
a0=2;
b0=2;
N=10;
mu_vec = ones(10,1);


%% c) negative free energy F

threshold = 0.0001;

prior_inizialize = false;
f_vec = program(a0,b0,mu0,lambda0,mu_vec,obs,N,threshold,prior_inizialize);

prior_inizialize = true;
f_vec_prior = program(a0,b0,mu0,lambda0,mu_vec,obs,N,threshold,prior_inizialize);

x = [0:length(f_vec)-1];
plot(x,f_vec)
hold on
x_prior =  [0:length(f_vec_prior)-1];
plot(x_prior,f_vec_prior)
legend('F with random inizialized values','F with values inizialized from prior')
title('Negative free energy, threshold = 0.0001')
xlabel('iteration')
ylabel('F')
hold off


