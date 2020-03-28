%% Generating inputs


u1 = zeros(800,1);
u2 = zeros(800,1);

for i= 1:800
    if i%7 == 0
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
[y, h, x] = euler_integrate_dcm(U, P, Phrf, x0, h0)