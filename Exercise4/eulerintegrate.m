function x = eulerintegrate(A,C,u)
x = zeros(4,200);

for i = 1:200
    dxdt = A*x(:,i)+ C*u(i);
    x(:,i+1) = x(:,i) + 0.001*dxdt;
end
end