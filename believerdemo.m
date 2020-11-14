%% 
 % prob = Forrestor();
 clc;
 clear all;

prob = 'Shekel_curve()';
 % prob= 'ackley(2, 2)';
% prob = 'tp3(1,1)';
% prob = 'rastrigin(1,1)';
rng(1, 'twister');
% [match_xl, n_fev, flag] = ego_BelieverGP([0, 0], prob, 1);
% [match_xl, n_fev, flag] = ego_Believer([0, 0], prob, 1);
% [match_xl, n_fev, flag] = ego_EI([0, 0], prob, 1);
% [match_xl, n_fev, flag] = ego_EIgp([0, 0], prob, 1);
% [match_xl, n_fev, flag] = ego_EIdace([0, 0], prob, 1);

[match_xl, n_fev, flag] = ego_HybridGP([0, 0], prob, 1);

disp(match_xl);

