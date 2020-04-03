function [y,h,x] = euler_integrate_dcm(u,P,Phrf,x0,h0)
    %construct matrices x and h and y 
    x = zeros(2,800);
    x(:,1) = x0;

    h1 = zeros(4,800);
    h2 = zeros(4,800);
    h1(:,1)= h0;
    h2(:,1)= h0;

    y = zeros(2,800);
    

    % call single_step_neural to calculate the temporal derivative at time
    % point t and perform euler integration to calculate x
    for i = 1:799
        dxdt = single_step_neural(x(:,i),u(:,i),P);
        x(:,i+1) = x(:,i) + 0.1*dxdt;
    end

    
    % call single_step_hrf to calculate the temporal derivative at time
    % point t and perform euler integration to calculate h
    for i = 1:799
        dhdt_x1 = single_step_hrf(h1(:,i),x(1,i),Phrf); %x1
        h1(:,i+1)=h1(:,i)+ 0.1*dhdt_x1;
        dhdt_x2 = single_step_hrf(h2(:,i),x(2,i),Phrf); %x2
        h2(:,i+1)=h2(:,i)+ 0.1*dhdt_x2;
    end

    h = [h1;h2];

    %compute bold signal
    for i = 1:800
        y(1,i) = compute_bold_signal(h1(:,i),Phrf);
        y(2,i) = compute_bold_signal(h2(:,i),Phrf);
    end


return

