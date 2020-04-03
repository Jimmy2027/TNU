function y = compute_bold_signal(h,Phrf)
    eps= 0.47;

    k1= 4.3*80.6*Phrf(5)*0.035;
    k2= eps*110*Phrf(5)*0.035;
    k3= 1- eps;

    v= h(3);
    q= h(4);
    y= 0.04 * (k1*(1-q)+k2*(1-q/v)+k3*(1-v));
return

