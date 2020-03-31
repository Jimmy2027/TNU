
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



[yo,ho,xo] = euler_integrate_dcm(U,P,Phrf,x0,h10);

figure; 
plot(axis, U.one, '-k','LineWidth', 1); hold on
plot(axis,U.two, '-', 'Color', 'g', 'LineWidth',2);
legend(['$u_{1}$'], ['$u_{2}$'], 'Interpreter', 'latex');

figure;
plot(axis, xo(1,2:end) , 'r-','LineWidth', 2); hold on
plot(axis,xo(2,2:end), 'g-', 'LineWidth' ,2);
legend(['$x_{1}$'], ['$x_{2}$'], 'Interpreter', 'latex');

figure;
plot(axis, yo.one , 'r-','LineWidth', 2); hold on
plot(axis,yo.two, 'g-', 'LineWidth' ,2);
legend(['$y_{1}$'], ['$y_{2}$'], 'Interpreter', 'latex');













%FUNCTION FOR POINT 1

function dxdt = single_step_neural(x,u,P)
    endo = (P.A + u.two*P.B);
    %replacement of diagonal elements 
    endo(1,1) = -0.5*exp(endo(1,1));
    endo(2,2) = -0.5*exp(endo(2,2));
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







