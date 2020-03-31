addpath('C:\Users\PC\Google Drive\ETH\Master\2 semester\Translational Neuromodelling\tapas\HGF')

mu3_0 = 1;
sigma3_0 = 1;
k2 = 2.5; w2 = -4;
w3 = -6;
zeta = 5; % response parameter

%BINARY INPUT FROM TAPAS
f = fopen("example_binary_input.txt");
data = textscan(f, '%f');
fclose(f);
input = cell2mat(data);

%SIMULATE MODEL
sim = tapas_simModel(input,'tapas_hgf_binary', [NaN 0 1 NaN 1 1 NaN 0 0 1 2.5 NaN -4 -6], 'tapas_unitsq_sgm', [zeta]);

est1 = tapas_fitModel(sim.y,sim.u,'tapas_hgf_binary_config_2','tapas_bayes_optimal_binary_config','tapas_quasinewton_optim_config');


%tapas_hgf_binary_plotTraj(sim);

tapas_fit_plotCorr(est1)