function dhdt = single_step_hrf(h,x,Phrf)
s = h(1);
f = h(2);
v = h(3);
q = h(4);
ka = Phrf.ka;
ga = Phrf.ga;
ta = Phrf.ta;
al = Phrf.al;
E0 = Phrf.E0;

ds = x - ka*s- ga*(f-1);
df = s;
dv = 1/ta * (f - v^(1/al));
dq = 1/ta*(f * (1-(1-E0)^(1/f))/E0 - v^(1/al)*(q/v));

dhdt = [ds,df,dv,dq];

return

