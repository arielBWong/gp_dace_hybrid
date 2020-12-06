%% 
 % prob = Forrestor();
 clc;
 clear all;

% prob = 'levy(2,2)';
% prob = 'Shekel(1, 1)';
% prob = 'ackley(1, 1)';
% prob = 'tp3(2,2)';
%  prob = 'rastrigin(1,1)';
% prob = 'Zakharov(2,2)';
prob = 'Shekel_curve()';
% rng(1, 'twister');
% prob =  'rosenbrock(2, 2)';
% [match_xl, n_fev, flag] = ego_bumpfindingGP3([0, 0], prob,1);
[match_xl, n_fev, flag] = ego_trustregiontest([0, 0], prob,1);
% [match_xl, n_fev, flag] = ego_BelieverGP([0, 0], prob, 1);
% [match_xl, n_fev, flag] = ego_Believer([0, 0], prob, 1);
% [match_xl, n_fev, flag] = ego_EI([0, 0], prob, 1);
% [match_xl, n_fev, flag] = ego_EIgp([0, 0], prob, 2);
% [match_xl, n_fev, flag] = ego_EIdace([0, 0], prob, 1);

% [match_xl, n_fev, flag] = ego_HybridGP([0, 0], prob, 1);

disp(match_xl);

