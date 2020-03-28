%Generate inputs



function inputs = generate(K2,W2,W3,x3_init,x2_init,u_init)
close all;
addpath('tapas/HGF')

x3 = zeros(101,1)
x2 = zeros(101,1)
x1 = zeros(101,1);
u = zeros(101,1);

    for i=1:100
        x3(i+1) = normrnd(x3_init,exp(W3));
        x2(i+1) = normrnd(x2_init,exp(K2*x3(i+1)+W2));
        x1(i+1) = binornd(1,1/(1+exp(x2(i+1))));
        u(i+1) = u_init^x1(i+1)*(1-u_init)^(1-x1(i+1));
        x3_init=x3(i+1);
        x2_init=x2(i+1);
        u_init = u(i+1);
    end
    inputs=[u,x2,x3];
    plot(u); 
    hold on;
    plot(x2);
    plot(x3);
    legend('u','x2','x3')
    hold off;
    
    r =tapas_fitModel([x2,x3],u);
    sim = tapas_hgf_binary(r,[NaN 0 1 NaN 1 1 NaN 0 0 0 1 NaN W2 W3]);
    tapas_hgf_binary_plotTraj(r)

    
end


