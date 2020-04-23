dt = 1e-3; %integration step

%constants
k1 = 80;
k2 = 50;
f1 = 50;
f2 = 50;
af = 3000;
ab = 1000;
c = 1;
mu = 0.05;
sigma = 0.01;

t = linspace(0,0.2,201);
u = normpdf(t,mu,sigma);

A = [0,1,0,0; -(k1^2),-f1, af,0;0,0,0,1;ab,0,-(k2)^2,-f2];
C = [0;0;0;c];

X = euler_int(A,C,u,dt);


%a 
figure(5);
plot(t,X(1,:),'-r');
hold on;
plot(t, x_condition_1(1,:),'--g');
legend('Euler integration', 'Given data');
title("Comparison between integrated data and given data");

%b
%plot difference between condtion 1 and condition 2 for x1
figure(6);
plot(t,x_condition_1(1,:),'-b');
hold on;
plot(t,x_condition_2(1,:), '-y');
legend("condition 1", "condition 2");
title('Comparison of x1 between condition 1 and condition 2');


%find the parameter that has been changed to get from x_condition_1 to 
%x_condition_2 between k1,k2,af,ab, and find its value
%}

res = zeros(2,length(pam));
%define array of functions to minimize
func_diff = {@(k1)(norm(x_condition_2 - euler_int([0,1,0,0; -(k1^2),-f1, af,0;0,0,0,1;ab,0,-(k2)^2,-f2],C,u,dt))),
    @(k2)(norm(x_condition_2 - euler_int([0,1,0,0; -(k1^2),-f1, af,0;0,0,0,1;ab,0,-(k2)^2,-f2],C,u,dt))),
    @(af)(norm(x_condition_2 - euler_int([0,1,0,0; -(k1^2),-f1, af,0;0,0,0,1;ab,0,-(k2)^2,-f2],C,u,dt))),
    @(ab)(norm(x_condition_2 - euler_int([0,1,0,0; -(k1^2),-f1, af,0;0,0,0,1;ab,0,-(k2)^2,-f2],C,u,dt)))};

for i = 1:length(pam)
    [pam_est,fval] = fminsearch(func_diff{i},0);
    res(:,i) = [pam_est;fval];
end 
res = array2table(res, 'VariableNames', {'k1','k2','af','ab'},'RowNames',{'Value of parameter','Difference'});
%form table is clear that parameter k1 was changed, in fact it has lowest
%difference in value from the orginal matrix x_condition_2

%plot the x's, original x_condition_2 and with best estiamated parameter
A = [0,1,0,0; -(res{1,1}^2),-f1, af,0;0,0,0,1;ab,0,-(k2)^2,-f2];

for i = 1:4
x = euler_int(A,C,u,dt);
figure(i);
plot(t,x_condition_2(i,:), '-g'); hold on
plot(t,x(i,:), '--r');
legend("Original data","Estimated data");
end











function [X] = euler_int(A,C,u,dt)
    X = zeros(4,length(u));
    for i = 1:length(u)
        dx = A*X(:,i) + C*u(i);
        X(:,i+1) = X(:,i) + dx*dt;
    end
    X = X(:,1:201);
end