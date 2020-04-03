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
sigma = 0.01;
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


%% d) evaluate MAP

joint =@param -post(param);
[map, feval] = fminsearch(log_post,zeros(1,5));

%% c) compute log joint distribution

%prior mit zero mean und sigma = 0.1
%log(llh*p) = log(llh)+log(p);

function log_post=post(param)
sigmap =0.1;
cov= sigmap^2 * eye(5);
invcovp= inv(cov);

lp =-0.5 * log(det(2*pi*cov))-0.5*param'*invcovp *param ;

log_post =  likeli(param)+lp;

end

%% b) compute log-likelihood 

function llh = likeli(param)

    P.A(1,1)=param(1);
    P.A(2,1)=param(2);
    P.A(2,2)=param(3);
    P.B(2.1)=param(4);
    P.C(1)= param(5);

    [y,h,x] = euler_integrate_dcm(u_vector,P,Phrf,x0,h0);

    cov_llh = eye(2).* sigma^2;
    determinant= det(2*pi*cov_llh);
    inv_cov= inv(cov_llh);

    for i=1:length(y(1))
            llh=llh -0.5* log(determinant)-0.5*(([y_hat(1,i);y_hat(2,i)]-[y(1,i);y(2,i)]).'*inv_cov*([y_hat(1,i);y_hat(2,i)]-[y(1,i);y(2,i)]));
    end
end


