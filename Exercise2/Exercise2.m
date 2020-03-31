%1.d Posterior for y=5
%{
[X1, X2] = meshgrid(-7:0.1:7, -7:0.1:7);
X = [X1(:) X2(:)];
Z = 2/2*pi*exp(-2*(sqrt(X1.^2+X2.^2)-5).^2);
surf(X1,X2,Z)
title('Posterior for y=5')
xlabel('x1');
ylabel('x2');
zlabel(['$P(\vec{x} \mid y=5 )$'],'Interpreter','latex');
%}



%1.e Prior
%{
[X1, X2] = meshgrid(0:0.1:10, -5:0.1:5);
X = [X1(:) X2(:)];
Z = 1/2*pi*exp(-0.5*((X1-6).^2+X2.^2));
surf(X1,X2,Z)
title(['Prior $\mathcal{N}(0,6)$'], 'Interpreter','latex')
xlabel('x1');
ylabel('x2');
zlabel(['$P(\vec{x})$'], 'Interpreter','latex');
hold on
plot3(ones(1,101)*6, [-5:0.1:5], zeros(1,101), 'r', 'LineWidth',2);
plot3([0:0.1:10], zeros(1,101), zeros(1,101), 'r', 'LineWidth',2);
%}


%1.e likelihood
%{
[X1, X2] = meshgrid(-7:0.1:7, -7:0.1:7);
X = [X1(:) X2(:)];
Z = 2/2*pi*exp(-2*(sqrt(X1.^2+X2.^2)-5).^2);
surf(X1,X2,Z)
title('Likelihood for y=5')
xlabel('x1');
ylabel('x2');
zlabel(['$P(y=5 \mid \vec{x})$'],'Interpreter','latex');
%}

%1.e posterior
%{
[X1, X2] = meshgrid(-15:0.1:15, -15:0.1:15);
X = [X1(:) X2(:)];
Z = 2/2*pi*1/2*pi*exp(-2*(sqrt(X1.^2+X2.^2)-5).^2)*exp(-0.5*((X1-6).^2+X2.^2)) ;
surf(X1,X2,Z)
title('Unnormalized posterior for y=5')
xlabel('x1');
ylabel('x2');
zlabel(['$P(y=5 \mid \vec{x})$'],'Interpreter','latex')
hold on 
plot3(ones(1,301)*6, [-15:0.1:15], zeros(1,301), 'r', 'LineWidth', 3); 
plot3([-15:0.1:15],ones(1,301)*5, zeros(1,301), 'r', 'LineWidth', 3); 
plot3([-15:0.1:15] ,ones(1,301)*-5, zeros(1,301), 'r', 'LineWidth', 3); 
%}


