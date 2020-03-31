%CONSTANTS
k2 = 1;
w2 = -4;
w3 = -6; % aka w3 = -6

addpath('C:\Users\PC\Google Drive\ETH\Master\2 semester\Translational Neuromodelling\tapas\HGF')

[u1,w22,w33] = binary(w3,k2,w2,320,0,0);


function [U, x2s, x3s] = binary(w3, k2, w2, n, x2_0, x3_0)
    x3s = [x3_0, zeros(1,n)];
    x2s = [x2_0, zeros(1,n)];
    U = zeros(1,n);
    for k = 1:n
        x3 = normrnd(x3s(k), exp(w3));
        x2 = normrnd(x2s(k), exp(k2*x3+w2));
        x1 = 1/(1+exp(-x2));
        u = binornd(1,x1);
        x3s(k+1) = x3;
        x2s(k+1) = x2;
        U(k) = u;
    end
    
    trials = linspace(1,n,n);
    hold on;
    plot(trials, U, 'o', 'MarkerSize', 2, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
    plot(trials, x3s(2:n+1), '-', 'LineWidth', 2, 'Color', 'r');
    plot(trials, x2s(2:n+1), '-', 'LineWidth', 2, 'Color', 'b');
    xlabel('Number of trials');
    legend(['$u^{(k)}$'], ['$x_{3}^{(k)}$'], ['$x_{(2)}^{(k)}$'], 'Interpreter', 'latex');
    
    
    est = tapas_fitModel([],U.','tapas_hgf_binary_config','tapas_bayes_optimal_binary_config','tapas_quasinewton_optim_config');
    sim = tapas_simModel(U.','tapas_hgf_binary',[NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN est.optim.final(13) est.optim.final(14)]);
    tapas_hgf_binary_plotTraj(sim);
                    
                    
       
    
    
end


