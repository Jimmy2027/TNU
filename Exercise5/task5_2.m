%(a)
mu = 1;
tao = 1;
obs = 1 + randn(1,10);

% (b) (c) (d)
mu0 = 0;
lbd0 = 3;
a0 = 2;
b0 = 2;

%function defined at the end of the  script

% (e)
%run with initial values as priors
[s1,m1,a1,b1,F1,c1] = update_eqs(obs,length(obs),lbd0,a0,b0,mu0,mu0,(b0/(a0*lbd0)),a0,b0,1e-5)

%runs with random values
[s2,m2,a2,b2,F2,c2] = update_eqs(obs,length(obs),lbd0,a0,b0,mu0,100*rand,100*rand,100*rand,100*rand,1e-5)
[s3,m3,a3,b3,F3,c3] = update_eqs(obs,length(obs),lbd0,a0,b0,mu0,100*rand,100*rand,100*rand,100*rand,1e-5)


%ALL FUNCTIONS

function [s_2,m,a,b,F1,c] = update_eqs(y,N,lbd0,a0,b0,mu0,m_i,s_i,a_i,b_i,threshold)
    c = 0;
    y_m = mean(y);
    F0 =-a_i*log(b_i) + gammaln(a_i) - gammaln(a0) + a0*log(b0) ...
    +0.5*log(lbd0) + log(sqrt(s_i)) -N/2*log(2*pi) +0.5;
    s_2 = 1/((a_i/b_i)*(N+lbd0));
    m = (lbd0*mu0+N*y_m)/(lbd0+N);
    a = a0 + (N+1)/2;
    b = b0 + 0.5*(sum((y-m).^2)+lbd0*(mu0-m)^2+(N+lbd0)*s_2);
    while 1
    F1 = -a*log(b) + gammaln(a) - gammaln(a0) + a0*log(b0) ...
    +0.5*log(lbd0) + log(sqrt(s_2)) -N/2*log(2*pi) +0.5;   
    if (abs(F1-F0) <= threshold)
        break
    end
    F0 = F1;
    s_2 = 1/((a/b)*(N+lbd0));
    m = (lbd0*mu0+N*y_m)/(lbd0+N);
    a = a0 + (N+1)/2;
    b = b0 + 0.5*(sum((y-m).^2)+lbd0*(mu0-m)^2+(N+lbd0)*s_2);
    c = c +1;
    end

end

