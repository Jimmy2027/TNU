function [x] = EulerIntegration(A,C,u)
    x = zeros(4,200);
    x_start = [0;0;0;0];
    x(:,1)=x_start;
    ts= 0.001;
    
    for i = 1:200
        dxdt = A*x(:,i)+ C*u(i);
        x(:,i+1)= x(:,i)+ ts*dxdt;
    end
end

