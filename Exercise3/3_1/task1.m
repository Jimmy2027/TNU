
clear all;
close all;

%% neural model construction

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



%% hrf model construction

%hemodynamic state vector at t=0 (s,f,v,q)
h0 = [0;1;1;1];

% prameters for hrf : kappa, gamma, tau, alpha and E_0
Phrf=[0.64,0.32,2,0.32,0.4];

%% compute dcm
[y,h,x] = euler_integrate_dcm(u_vector,P,Phrf,x0,h0);


%%plot results
t = linspace(0,80,800);

figure;

subplot(3,2,1);
plot(t,u_vector(:,:));
ylim([0,6]);
title('Inputs')
legend('u1','u2')

subplot(3,2,3);
plot(t,x(:,:));
ylim([0,0.6]);
title('Neural activity')
legend('x1','x2')

subplot(3,2,5);
plot(t,y(:,:));
ylim([0,0.08]);
title('BOLD Signal')
legend('y_{x1}','y_{x2}')
xlabel('time[s]')

subplot(3,2,2);
plot(t,h(1:4,:));
title('Hemodynamic states for x1')
legend('s','f','v','q')

subplot(3,2,4);
plot(t,h(5:8,:));
title('Hemodynamic states for x2')
legend('s','f','v','q')
xlabel('time[s]')

savefig('task1_figure')
