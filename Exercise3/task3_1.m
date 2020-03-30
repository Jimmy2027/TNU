%% Generating inputs


u1 = zeros(800,1);
u2 = zeros(800,1);

for i= 1:800
    if rem(i,70) == 0
        u1(i) = 5;
    else
        u1(i) = 0;
    end
    if i > 299 && i < 601
        u2(i) = 1;
    else
        u2(i) = 0;
    end
end

U = [u1, u2];
Phrf = struct('ka', 0.64, 'ga', 0.32, 'ta', 2, 'al', 0.32, 'E0', 0.4);
A = [-0.5, 0; 1, -0.5];
B = [0,0;-0.5, 0];
C = [1, 0;0,0];
P = struct('A',A, 'B',B, 'C',C);
h0 = [0,1,1,1];
x0 = [0, 0];
[y, h1, h2, x] = euler_integrate_dcm(U, P, Phrf, x0, h0);

figure(1)
plot(u1)
title('Input')
hold on;
plot(u2)
legend('u1', 'u2')
hold off;

figure(2)
plot(x(:,1))
hold on;
plot(x(:,2))
title('Neural activity')
legend('x1','x2')
hold off;

figure(3)
plot(y(:,1))
hold on;
plot(y(:,2))
title('BOLD signal')
legend('yBold1','yBold2')
hold off;
