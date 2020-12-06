%%
problems ={'ackley(3, 3)', 'levy(3, 3)','rastrigin(3, 3)','dsm1(3, 3)', ... %  multimodal global structure  heavy modality and weak modality
            'tp3(3, 3)', 'tp5(3, 3)', 'tp7(3, 3)','Shekel(3, 3)', ... % multimodal no global structure
            'Zakharov(3, 3)', 'smd2(3, 3)',  'rosenbrock(3, 3)'}; % unimodal
% 'tp3(3, 3)', 'tp5(3, 3)', 'tp7(3, 3)','Shekel(3, 3)','rastrigin(3, 3)', 'rosenbrock(3, 3)', 'Zakharov(3, 3)', 'smd2(3, 3)'

problems = {'ackley(3, 3)', 'Zakharov(3, 3)', 'Shekel(3, 3)' ,  'rosenbrock(3, 3)'};
seeds = linspace(12, 21, 10);
match_methods = {'ego_bumpfindingGP3'}; % 'ego_EI', 'ego_Believer', 'ego_BelieverGP', 'ego_EIgp'


np = length(problems);  
ns = length(seeds);
na = length(match_methods);

parameterfile = strcat(pwd, '\parameter_ble.m');
run(parameterfile) 
parameterfile = strcat(pwd, '\parameter_ei.m');
run(parameterfile)
parameterfile = strcat(pwd, '\parameter_hyb.m');
run(parameterfile)

paras=cell(1, np * ns * na);
nn = 1; 
for i = 1:np
    for j = 1:ns
        for k = 1:na
            paras{ nn } = {problems{i}, seeds(j), str2func(match_methods{k})};
            nn = nn + 1;
        end
    end
end

nrun = length(paras);
parfor i =1:nrun
    
tic;
paras{i}{3}([0, 0, 0], paras{i}{1},  paras{i}{2});
toc
end
