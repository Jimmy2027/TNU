
function dxdt = single_step_neural(x,u,P)
    %r1 = (-0.5)*exp(P.A(1,1)+u(2)*P.B(1,1));
    %r2 =  (-0.5)*exp(P.A(2,2)+u(2)*P.B(2,2));
    con = (P.A + u(2)*P.B);
    %con_new = [r1,con(1,2);con(2,1),r2];
    
    dxdt = con*x + P.C*u(1);
return