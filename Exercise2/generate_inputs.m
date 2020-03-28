%Generate inputs

function inputs = generate(K2,W2,W3,x3_init,x2_init,u_init)
x3 = zeros(101,1);
x2 = zeros(101,1);
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
end

