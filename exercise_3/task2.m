clear all;
close all;

%If B = B1 = 0, the dynamics of the neural system are independent of the
%modulatory input. Descriptively we see that x1 and x2 are periodic and
%u2=1 has no influence on them.
%If B = B2, the modulatory input influences the dynamics of the neural
%state x2, linked with the previous state x1. Descriptively we see that if b>0, u2>0 has an positive impact on
%x2 (x2 increases) and b<0, u2>0 has a negative impact on x2 (x2
%decreases).
%If B = B3, the modulatory input also influences the dynamics of the state
%x2, linked with the previous state x2. Descriptively we see that if b<0,
%u2>0 has a negative impact on x2 (x2 decreases). If b>0, the positive impact of 
%u2>0 is changed by setting the diagonal terms of (A+u2B) to -0.5exp(aii+u2bii), for the sability of the system.
%This explains the negative impact of b=1 and the reduced positive effect
%of b=0.5.


%% define which B-matrix is used
index = 3;
b_vector= linspace(-1,1,5);

for j = 1:5
    %% construct struct P of matrices A and B and vector C
    P.A = [-0.5,0;1,-0.5];
    P.B = createB(index,b_vector(j));
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


    %% plot results
    t = linspace(0,80,800);
    figure;
    

    subplot(3,1,1);
    plot(t,u_vector(:,:));
    ylim([0,6]);
    title('Inputs')
    legend('u1','u2')

    subplot(3,1,2);
    plot(t,x(:,:));
    %ylim([0,0.6]);
    title('Neural activity')
    legend('x1','x2')

    subplot(3,1,3);
    plot(t,y(:,:));
    %ylim([0,0.08]);
    title('BOLD Signal')
    legend('y_{x1}','y_{x2}')
    xlabel('time[s]')

    
    if(index == 1)
        sgtitle(['B',num2str(index)])
        break;
    else
        sgtitle(['B',num2str(index), ' with b = ',num2str(b_vector(j))])
    end
end
