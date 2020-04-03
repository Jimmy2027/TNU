
function dxdt = single_step_neural(x,u,P)
    con = (P.A + u(2)*P.B);
    r1 = (-0.5)*exp(con(1,1));
    r2 =  (-0.5)*exp(con(2,2));
    con_new = [r1,con(1,2);con(2,1),r2];
    
    % replace diagonal entries of connectivity with -0.5exp(aii+u2bii) to ensure
    % stability if entries are >0
    
    if con(1,1)>=0 || con(2,2)>=0
        dxdt = con_new*x + P.C*u(1);
    else
        dxdt = con*x + P.C*u(1);
    end
return