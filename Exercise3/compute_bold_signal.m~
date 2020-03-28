function y = compute_bold_signal(h,Phrf)
eps = 0.47;
E0 = Phrf.E0;
TE = 0.035;
Th0 = 80.6;
r0 = 110;
V0 = 0.04;

k1 = 4.3 * Th0 * E0 * TE;
k2 = eps * r0 * E0 * TE;
k3 = 1 - eps;

v = h(3);
q = h(4);
y = V0 * (k1*(1-q)+k2*(1-(q/v))+k3*(1-v));
return
