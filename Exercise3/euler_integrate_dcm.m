function [y,h,x] = euler_integrate_dcm(U,P,pHRF,x0,h0)
y = zeros(800,1);
h = zeros(800,4);
x = zeros(800,2);
h(1,:) = h0;
x(1,:) = x0;
y(1,:) = compute_bold_signal(h(1,:), pHRF);

for i= 2:800
    dxdt = single_step_neural(x(i-1,:), U(i,:), P)
    x(i,:) = [exp(dxdt(0,0) + exp(dxdt(0,1), exp(1,0) + exp(1,1)]
    h(i,:) = single_step_hrf(h(i-1,:), x(i-1,:), pHRF)
    y(i,:) = compute_bold_signal(h(i-1,:), pHRF)
end

return
    