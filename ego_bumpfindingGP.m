function[match_xl, n_fev, flag] = ego_bumpfindingGP(xu, prob, seed)

% parameterfile = strcat(pwd, '\parameter_ble.m');
% run(parameterfile)
load('param_hyb');
rng(seed, 'twister');


num_pop = param_hyb.num_pop;
num_gen = param_hyb.num_gen;
init_size = param_hyb.init_size;
iter_size = param_hyb.iter_size; 
norm_str = param_hyb.norm_str;
par = param_hyb.par;

prob = eval(prob);
l_nvar = prob.n_lvar;


upper_bound = prob.xl_bu;
lower_bound = prob.xl_bl;

% init_size = 11 * l_nvar - 1;
xu_init = repmat(xu, init_size, 1);
train_xl = lhsdesign(init_size,l_nvar,'criterion','maximin','iterations',1000);
train_xl = repmat(lower_bound, init_size, 1) ...
    + repmat((upper_bound - lower_bound), init_size, 1) .* train_xl;

initx = train_xl;

% evaluate/get training fl from xu_init and train_xl
% compatible with non-constriant problem
[train_fl, train_fc] = prob.evaluate_l(xu_init, train_xl);

if size(train_xl, 2) == 1
    fighn = figure(1);
else
    fighn = [];
end

normhn = str2func(norm_str);
mdl_save = cell(1, iter_size);
fithn = str2func('EIM_eval');
switch_record = [];
daceflag = false;

for g = 1: iter_size   
    tic;
    [krg, krgc] = surrogate_create(train_xl, train_fl, train_fc, normhn, daceflag);
    disconnected_bump = bump_detection_disconnectedRegion(train_xl, train_fl, prob.xl_bl, prob.xl_bu);
    toc;
    if disconnected_bump   % multimodal
        fprintf('multi modal mode iteration %d \n', g);
        switch_record = [switch_record, 1];
        [new_xl, ~] = switchEIM_gpr(train_xl, train_fl, prob.xl_bu, prob.xl_bl, ...
            num_pop, num_gen, train_fc, fithn, normhn, krg, krgc);
    else                    % single modal
        fprintf('uni modal mode iteration %d \n', g);
        switch_record = [switch_record, 0];
        [new_xl, ~] = switchBeliever_gpr(train_xl, krg, krgc, ...
                                 num_pop, num_gen, daceflag, prob);
    end
    
    if g==8
        a = 0;
    end
    
    if size(train_xl, 2) == 1
        processplot1d(fighn, train_xl, train_fl, krg, prob, initx, new_xl, daceflag)
    end  
 
    mdl_save{g} = krg;
    
    % fprintf('problem %s, iter %d, seed %d, dace flag is %d \n', prob.name,g,  seed, daceflag);
    [new_fl, new_fc] = prob.evaluate_l(xu, new_xl);
    
    %--- closeness check---
    check = abs(train_xl - new_xl);
    check = round(check, 5);
    if sum(check < 1e-20) >= 1     
        fprintf('fail unique check %d \n', g);
        % continue;
    end
    
    train_xl = [train_xl; new_xl];
    train_fl = [train_fl; new_fl];
    train_fc = [train_fc; new_fc];  % compatible with nonconstraint
 
end

[match_xl, n_fev, flag] = post_infillsaveprocess(xu, train_xl, train_fl, train_fc, 'lladp_gp', seed, prob, init_size);
save_eachmodel(seed, 'lladp_gp', mdl_save, init_size, prob, switch_record);
end


