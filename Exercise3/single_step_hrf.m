function dhdt = single_step_hrf(h,x,Phrf)
s = h(1);
f = h(2);
v = h(3);
q = h(4);
ka = Phrf(1);
ga = Phrf(2);
ta = Phrf(3);
al = Phrf(4);
E0 = Phrf(5);

ds = x - ka*s- ga*(f-1);
df = s;
dv = 1/ta * (f - v^(1/al));
dq = 1/ta*(f * (1-(1-E0)^(1/f))/E0 - v^(1/al)*(q/v));

dhdt = [ds;df;dv;dq];

return


