T = readtable("NN2emv.csv");




%FIXED PARAMETERS
theta0 = 0.3;
theta1 = -0.1;
theta2 = 0.5;
thetas = [theta0, theta1, theta2];
sigma = sqrt(0.001);

x = linspace(-0.5,0.5, 101);
est_parm = NaN(100, 8, 3); %allocate estimated parameter tensor before to speed up
P = [1 2 7]; %grade of the polynomial

for p = 1:3
for k = 1:100
eps = normrnd(0, sigma, [1,length(x)]); %generate white gaussian noise

y = nan(size(x));
for n = 1:length(x)
    y(n) = theta0 + theta1*x(n) + theta2*x(n)^2 + eps(n);
end

%create the data matrix X
X = NaN(length(x), P(p)+1);

for n = 1:length(x)     %fill data matrix X
    for j = 1:P(p)+1
        X(n, j) = x(n)^(j-1);
    end
end
    
mlest = inv(X.'*X)*X.'*y.';

est_parm(k,[1:P(p)+1], p) = mlest.';
end
end


%histogram for theta0
figure;hold on;
histogram(est_parm(:,1,1), 15, 'facecolor', 'red', 'facealpha', 0.4 );
histogram(est_parm(:,1,3), 15, 'facecolor', 'yellow', 'facealpha', 0.4 );
histogram(est_parm(:,1,2), 15, 'facecolor', 'blue', 'facealpha', 0.4 );
line([theta0, theta0], ylim, 'LineWidth', 1, 'Color', 'r')
xlabel('\theta_{0}');
ylabel('Number of instances');
legend('P=1', 'P=7', 'P=2', 'True \theta_{0}', 'location', 'northwest'); 
title('ML estimation of \theta_{0}');


%histogram for theta1
figure;hold on;
histogram(est_parm(:,2,3), 15, 'facecolor', 'yellow', 'facealpha', 0.4 );
histogram(est_parm(:,2,1), 15, 'facecolor', 'red', 'facealpha', 0.4 );
histogram(est_parm(:,2,2), 15, 'facecolor', 'blue', 'facealpha', 0.4 );
line([theta1, theta1], ylim, 'LineWidth', 1, 'Color', 'r')
xlabel('\theta_{1}');
ylabel('Number of instances');
legend('P=7', 'P=1', 'P=2', 'true \theta_{1}', 'location', 'northwest');
title('ML estimation of \theta_{1}');


%histogram for theta0
figure;hold on;
histogram(est_parm(:,3,3), 15, 'facecolor', 'yellow', 'facealpha', 0.4 );
histogram(est_parm(:,3,2), 15, 'facecolor', 'blue', 'facealpha', 0.4 );
line([theta2, theta2], ylim, 'LineWidth', 1, 'Color', 'r')
xlabel('\theta_{2}');
ylabel('Number of instances');
legend('P=7', 'P=2', 'true \theta_{2}', 'location', 'northwest');
title('ML estimation of \theta_{2}');


%ll = likelihood(X,y,mlest, sigma ,x)


function L = likelihood(X, y, theta, sigma, x)
    L = -length(x)*log(sqrt(2*pi)) - 1/(2*sigma^2)*norm(y.'-X*theta)^2
    
end
