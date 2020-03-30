function [y,h1,h2,x] = euler_integrate_dcm(U,P,pHRF,x0,h0)
y = zeros(800,2);
h1 = zeros(800,4);
h2 = zeros(800,4);
x = zeros(800,2);
h1(1,:) = h0;
h2(1,:) = h0;
x(1,:) = x0;
y(1,1) = compute_bold_signal(h1(1,:), pHRF);
y(1,2) = compute_bold_signal(h2(1,:), pHRF);

for i= 2:800
    dxdt = single_step_neural(x(i-1,:), U(i,:), P);
    x(i,:) = [dxdt(1,1)*x(i-1,1) + dxdt(1,2)*x(i-1,2), (dxdt(2,1)*x(i-1,1) + dxdt(2,2)*x(i-1,2))];
    h1(i,:) = single_step_hrf(h1(i-1,:), x(i-1,1), pHRF);
    h2(i,:) = single_step_hrf(h2(i-1,:), x(i-1,2), pHRF);
    y(i,1) = compute_bold_signal(h1(i-1,:), pHRF);
    y(i,2) = compute_bold_signal(h2(i-1,:), pHRF);
end

return
    