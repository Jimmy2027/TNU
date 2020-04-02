%IMPULSE TRAINS
u = zeros(2, 800);
u(1,70:70:end) = 5;
u(2,300:1:600) = 1;

%INPUT STRUCT
U.one = u(1,:);
U.two = u(2,:);

%NEURONAL DYNAMIC STRUCT
P.A = [-0.5, 0; 1, -0.5];
P.B = [0, 0; -0.5, 0];
P.C = [1, 0; 0, 0];

%HEMODYNAMIC STRUCT
Phrf.k = 0.64;
Phrf.gamma = 0.32;
Phrf.tao = 2;
Phrf.alpha = 0.32;
Phrf.E = 0.4;

%INITIALI CONDITIONS
x0 = [0;0];
h10 =  [0;1;1;1];
h20 = [0;1;1;1];

axis = linspace(1,800,800);

%3.A
[yo,ho,xo] = euler_integrate_dcm(U,P,Phrf,x0,h10);
sigmah = 0.004; %RANDOM VARIANCE THAT LOOKS GOOD ON THE PLOT 
noise = normrnd(0,sigmah,[2,800]);
nn.one = yo.one + noise(1,:); %noisy bold1
nn.two = yo.two + noise(2,:); %noisy bold2

figure;
plot(axis, nn.one , 'g-','LineWidth', 2); hold on
plot(axis,yo.one, 'r-', 'LineWidth' ,2);
legend(['$\hat{y_{1}}$'], ['$y_{1}$'], 'Interpreter', 'latex', 'FontSize', 18);
title(['BOLD signal 1 with AWGN'],'Interpreter', 'latex', 'FontSize', 15);

figure;
plot(axis, nn.two , 'g-','LineWidth', 2); hold on
plot(axis,yo.two, 'r-', 'LineWidth' ,2);
legend(['$\hat{y_{2}}$'], ['$y_{2}$'], 'Interpreter', 'latex', 'FontSize', 18);
title(['BOLD signal 2 with AWGN'],'Interpreter', 'latex', 'FontSize', 15);





%3.D
lk = @(A)-likelihoo(nn,A); %FIX NOISY INPUT, BECOMES FUNCTION OF PARAMETERSS ONLY
[ml, mval] = fminsearch(lk, [0,0,0,0,0]); %ML ESTIMATION
ml_t = table(ml, mval);
ml_t.Properties.VariableNames = {'Parameters', 'log-likelihood'};
ml_t



map_t = table(["a11";"a21"; "a22"; "b21"; "c11";"log-L"], NaN(6,1), NaN(6,1), NaN(6,1), NaN(6,1));
map_t.Properties.VariableNames = {'Parameters', 'sigma_1^2 =10', 'sigma _2^2 =0.1', 'sigma_3^2 =0.001', 'sigma_1^2 =10^-5'};

%MAP FOR DIFFERENT VAR OF PRIORS
vars = [10,0.1,0.001, 0.00001]; %VARIOUS VARIANCES OF CONNECTIVITY PARM.
for i=1:length(vars)
    
    
    joi = @(A)-joint(nn,A,vars(i)); %anonymous function where we fix noisy input and variance of parameter, becomes function of parameters vector A only.
    [map, feval] = fminsearch(joi, [0,0,0,0,0]); %MAP ESTIMATE
    map_t{:,i+1} = [map.';likelihoo(nn,map)];
    
    
    
end
map_t

%3.B LOG LIKELIHOOD

function L = likelihoo(yn, A) %yn is vector of noisy observations y_hat, sigma noise std., a is vector of connectivity parameters
    var = 0.004^2;                  %STANDARD DEVIATION OF NOISE KNOWN (AS SAID IN PROBLEM STATEMENT) %U,pHRF,x0.h10 FIXED AS IN PROB:STAT:
    %%FIXED PARAMETERS OF BOLD SIGNAL
    Phrf.k = 0.64;
    Phrf.gamma = 0.32;
    Phrf.tao = 2;
    Phrf.alpha = 0.32;
    Phrf.E = 0.4;
    %IMPULSE TRAINS
    u = zeros(2, 800);
    u(1,70:70:end) = 5;
    u(2,300:1:600) = 1;
    %INPUT STRUCT
    U.one = u(1,:);
    U.two = u(2,:);
    x0 = [0;0];
    h10 =  [0;1;1;1];
    P.A(1,2)=0;
    P.B = zeros(2,2);
    P.C = zeros(2,2);
    %VARIABLE CONNECTIVITY PARAMETERS
    P.A(1,1)= A(1);
    P.A(2,1) = A(2);
    P.A(2,2) = A(3);
    P.B(2,1) = A(4);
    P.C(1,1) = A(5);
    L=0; %initial value of likelihood
    [ko,ho,xo] = euler_integrate_dcm(U,P,Phrf,x0,h10); %NON NOISY BOLD SIGNAL CREATED WITH PARAMETERS IN VECTOR A, WE NEED ONLY KO (BOLD)
    
    for i=1:length(ko.one)
        L=L + log(1/(2*pi*var))+ (-1/2*([yn.one(i);yn.two(i)]-[ko.one(i);ko.two(i)]).'*[1/var,0;0,1/var]*([yn.one(i);yn.two(i)]-[ko.one(i);ko.two(i)]));
    end

end

%3.C PROBABILITY OF JOINT
function JP = joint(yn, A, var_p) %A is a vector whose entries correspond to a11,a21,a22,b21,c11 in this order, var_p is the cariance of the parameters
                                      
    l_prior_p = (2*pi)*(-5/2)*(var_p^5)^(-0.5)*exp((-1/2*[A(1);A(2);A(3);A(4);A(5)].'*eye(5)*1/var_p*[A(1);A(2);A(3);A(4);A(5)]));
    
    JP = likelihoo(yn,A) + log(l_prior_p);

end



























%FUNCTIONS FORM TASKS 1 AND 2
function dxdt = single_step_neural(x,u,P)
    endo = (P.A + u.two*P.B);
    %replacement of diagonal elements 
    if endo(1,1)>=0 || endo(2,2) >=0
        endo(1,1) = -0.5*exp(endo(1,1));
        endo(2,2) = -0.5*exp(endo(2,2));
    end
    dxdt = endo*x + P.C*[u.one;0];
end

%FUNCTION FOR POINT 2

function dhdt = single_step_hrf(h,x,Phrf)
    svar = x-Phrf.k*h(1)-Phrf.gamma*(h(2)-1);
    fvar = h(1);
    vvar = 1/Phrf.tao *(h(2)-h(3)^(1/Phrf.alpha));
    qvar = 1/Phrf.tao * (h(2)*((1-(1-Phrf.E)^(1/h(2)))/Phrf.E)- h(3)^(1/Phrf.alpha)*h(4)/h(3));
    dhdt = [svar; fvar; vvar; qvar];
end
    


%FUNCTION FOR POINT 3

function y = compute_bold_signal(h, Phrf)
    %Parameters for BOLD at 3T
    theta0 = 80.6;
    r0 = 110;
    TE = 0.035;
    epsilon = 0.47;
    V0 = 0.04;
    k1 = 4.3*theta0*Phrf.E*TE;
    k2 = epsilon*r0*Phrf.E*TE;
    k3 = 1-epsilon;
    y = V0*(k1*(1-h(4)) + k2*(1-h(4)/h(3)) + k3*(1-h(3)));
end
       
       
    
function [y, h, x] = euler_integrate_dcm(U,P,pHRF,x0,h0)
    %FIRST INTEGRATE AKA SUM X AND H
    ht1 = h0;
    ht2 = h0;
    xt = x0;
    x = [x0];
    h1 = [h0];
    h2 = [h0];
    for i=1:800
        u.one = U.one(i);
        u.two = U.two(i);
        dxdt = single_step_neural(xt,u,P);
        xt = xt+ 0.1*dxdt;
        x = [x, xt];
        dhdt_neuron1 = single_step_hrf(ht1,xt(1),pHRF);
        dhdt_neuron2 = single_step_hrf(ht2,xt(2),pHRF);
        ht1 = ht1 + 0.1*dhdt_neuron1;
        ht2 = ht2 + 0.1*dhdt_neuron2;
        h1 = [h1, ht1];
        h2 = [h2, ht2];
    end
    
    y1 = [];
    y2 = [];
    for k = 1:800
        dy1 = compute_bold_signal(h1(:,k),pHRF);
        dy2 = compute_bold_signal(h2(:,k),pHRF);
        y1 = [y1, dy1];
        y2 = [y2, dy2];
    end
    y.one = y1;
    y.two = y2;
    h.one = h1;
    h.two = h2;
    
end





