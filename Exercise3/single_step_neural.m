function dxdt = single_step_neural(x,u,P)
A = P.A;
B = P.B;
C = P.C;

dxdt = (A + u(2)*B)*x' + C*u(1)
end