%% Task 2.3
%% 
addpath('../tapas/HGF')
clear all;
close all;

u = load('example_binary_input.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% a)k2=2.5, w2=-4, w3=-6, mu3=1, sa3=1, ze =5

% The \omega_3 is very accurate in all our results.
% The estimated parameters are not as accurate as the set ones since the
% variance is not zero. For example the \kappa_2 is not accurate when we
% estimate it but nearly 2.5 when not estimating this parameter.

% The first covariance plot shows large correlation between \kappa_2 and
% \mu_0.
% The second covariance plot shows large correlation between \omega_2 and
% \mu_0.


%simulate model
sim = tapas_simModel(u,...
'tapas_hgf_binary',...
[NaN 0 1 NaN 1 1 NaN 0 0 1 2.5 NaN -4 -6],...
'tapas_unitsq_sgm',...
5);

%estimate param: (ze,mu3,K2,exp(w3))
est1 = tapas_fitModel(sim.y,...
                     sim.u,...
                     'tapas_hgf_binary_config_2',...
                     'tapas_unitsq_sgm_config',...
                     'tapas_quasinewton_optim_config')

                 
%plot posterior correlation
tapas_fit_plotCorr(est1)

%plot trajectories
tapas_hgf_binary_plotTraj(est1)


%estimate param: (ze,mu3,w2,exp(w3))
est2 = tapas_fitModel(sim.y,...
                     sim.u,...
                     'tapas_hgf_binary_config_3',...
                     'tapas_unitsq_sgm_config',...
                     'tapas_quasinewton_optim_config')

                 
%plot posterior correlation
tapas_fit_plotCorr(est2)

%plot trajectories
tapas_hgf_binary_plotTraj(est2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% b)k2=1, w2=-4, w3=-4.1674, mu3=2.5, sa3=6.25, ze = 5

% The \omega_3 is very accurate in all our results.
% The estimated parameters are not as accurate as the set ones since the
% variance is not zero. For example the \kappa_2 is not accurate when we
% estimate it but equal to 1 when not estimating it.

% The estimates for \zeta and \theta are more accurate in b) than in a).
% We think a possible explanation for the better estimates is that the
% volatility coeffitient \theta is higher and therefore the model is more
% flexible

% The first covariance plot shows large correlation between \kappa_2 and
% \mu_0.
% The second covariance plot shows large correlation between \omega_2 and
% \mu_0.

%simulate model
sim2 = tapas_simModel(u,...
'tapas_hgf_binary',...
[NaN 0 2.5 NaN 1 6.25 NaN 0 0 1 1 NaN -4 -4.1674],...
'tapas_unitsq_sgm',...
5);

%estimate model
est3 = tapas_fitModel(sim2.y,...
                     sim2.u,...
                     'tapas_hgf_binary_config_4',...
                     'tapas_unitsq_sgm_config',...
                     'tapas_quasinewton_optim_config')

                 
%plot posterior correlation
tapas_fit_plotCorr(est3)

%plot trajectories
tapas_hgf_binary_plotTraj(est3)

%estimate model
est4 = tapas_fitModel(sim2.y,...
                     sim2.u,...
                     'tapas_hgf_binary_config_5',...
                     'tapas_unitsq_sgm_config',...
                     'tapas_quasinewton_optim_config')

                 
%plot posterior correlation
tapas_fit_plotCorr(est4)

%plot trajectories
tapas_hgf_binary_plotTraj(est4)




%% c)

%The \mu_3 can only be changed without changing the other beliefs by
%changing \sigma_3 or \kappa_2.
%% d) tapas_unitsq_sgm_mu3 as response model
%Including a readout of $\mu _3$ in the response model, we now see much
%less correlation between the observed parameters than in a).
%The estimate for $\mu _3$ is -1.0342 while the actual value is 1.
%The estimate for $\kappa _2$ is much more accurate: it is 0.3246 while the actual value is 1.
%The estimate for $\omega _2$ is for off the original value: the estimate
%is 0.0820 while the true value is -4.

%simulate model
sim = tapas_simModel(u,...
'tapas_hgf_binary',...
[NaN 0 1 NaN 1 1 NaN 0 0 1 2.5 NaN -4 -6],...
'tapas_unitsq_sgm',...
5);

%estimate param: (ze,mu3,K2,exp(w3))
est1 = tapas_fitModel(sim.y,...
                     sim.u,...
                     'tapas_hgf_binary_config_2',...
                     'tapas_unitsq_sgm_mu3_config',...
                     'tapas_quasinewton_optim_config')

                
%plot posterior correlation
tapas_fit_plotCorr(est1)

%plot trajectories
tapas_hgf_binary_plotTraj(est1)


%estimate param: (ze,mu3,w2,exp(w3))
est2 = tapas_fitModel(sim.y,...
                     sim.u,...
                     'tapas_hgf_binary_config_3',...
                     'tapas_unitsq_sgm_mu3_config',...
                     'tapas_quasinewton_optim_config')

                 
%plot posterior correlation
tapas_fit_plotCorr(est2)

%plot trajectories
tapas_hgf_binary_plotTraj(est2)
%% 
close all


