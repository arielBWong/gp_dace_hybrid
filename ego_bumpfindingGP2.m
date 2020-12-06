function[match_xl, n_fev, flag] = ego_bumpfindingGP2(xu, prob, seed)
% based on EI
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
local_ax = [];
local_ay = [];

for g = 1: iter_size
    fprintf('generation %d \n', g);
    if g==20
        a = 0;
    end
    
    [krg, krgc] = surrogate_create(train_xl, train_fl, train_fc, normhn, daceflag);
    % disconnected_bump = bump_detection_disconnectedRegion(train_xl, train_fl, prob.xl_bl, prob.xl_bu);
    
    % switch_record = [switch_record, 1];
    [new_xl, ~] = switchEIM_gpr(train_xl, train_fl, prob.xl_bu, prob.xl_bl, ...
        num_pop, num_gen, train_fc, fithn, normhn, krg, krgc);
    
    
    % plot 2d with bump hunting
    if size(train_xl, 2) == 2
        processplot2d_bumpstudy(fighn, train_xl, train_fl, krg, prob, initx, new_xl, daceflag)
    end
    
    
    
    % local search
    trainx_all = [train_xl; local_ax];
    trainy_all = [train_fl; local_ay];
    
    prim =  bump_detection_4test(trainx_all, trainy_all, prob.xl_bl, prob.xl_bu);
    nboxes = length(prim.boxes);
    if nboxes > 2
        upper_stack = [];
        lower_stack = [];
        for ii = 1:2
            fprintf('box %d, mean %f, support %d \n', ii, prim.boxes{ii}.mean, prim.boxes{ii}.supportN);
            % adjust boundary
            mentioned_var = prim.boxes{ii}.vars;
            bump_lb = prob.xl_bl;
            bump_ub = prob.xl_bu;
            
            nb = length(mentioned_var);
            if ~isempty(mentioned_var)
                for jj = 1:nb % varible indicator
                    if ~isnan(prim.boxes{ii}.min(jj))
                        bump_lb(mentioned_var(jj)) =  prim.boxes{ii}.min(jj);
                    end
                    
                    if ~isnan(prim.boxes{ii}.max(jj))
                        bump_ub(mentioned_var(jj)) =  prim.boxes{ii}.max(jj);
                    end
                    
                end
            end
            upper_stack = [upper_stack; bump_ub];
            lower_stack = [lower_stack; bump_lb];
        end
        
        upper_bound = max(upper_stack, [], 1);
        lower_bound = min(lower_stack,[],  1);
        
        % run local search
        [localx, localy] = pickup_localxy(upper_bound, lower_bound, trainx_all, trainy_all);
        kk = size(localx, 1);
        fprintf('Local training data size %d \n', kk);
        % build local surrogate
        [local_krg, local_krgc] =  surrogate_create(localx, localy, [] , normhn, daceflag);
        
        % run believer to find next point
        [local_newx, ~] = localBeliever_gpr(localx, local_krg, local_krgc, ...
            num_pop, num_gen, daceflag, prob, lower_bound, upper_bound);
        % plot local search
        plot1d_global_withbump(fighn, prim, prob, train_xl,train_fl);
        plotlocal_search(fighn, upper_bound, lower_bound, localy, local_krg);
        
        [local_newy,~] = prob.evaluate_l(xu, local_newx);
        
        local_ax = [local_ax; local_newx];
        local_ay = [local_ay; local_newy];
        
    end
    
     if size(train_xl, 2) == 1
        processplot1d_bumpstudy(fighn, trainx_all, trainy_all, krg, prob, initx, new_xl, daceflag)
        % processplot1d(fighn, train_xl, train_fl, krg, prob, initx, new_xl, daceflag)
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


function [localx, localy] = pickup_localxy(upper_bound, lower_bound, xx, yy)
localx = [];
localy = [];
numx = size(xx, 1);
nx = size(xx, 2);
for ii = 1: numx
    x = xx(ii, :);
    y = yy(ii, :);
    
    % with in boundary
    if sum(upper_bound - x >= 0 ) == nx && sum(x - lower_bound >= 0) == nx
        localx = [localx; x];
        localy = [localy; y];
    end
end
end

function plotlocal_search(fighn, upper_bound, lower_bound, trainy, krg)
testdata = linspace(upper_bound, lower_bound, 100);
testdata = testdata';

[pred, ~] = surrogate_predict(testdata, krg, false);
pred = denormzscore(trainy, pred);
fighn;
plot(testdata, pred);

pause(1);

end

function plot1d_global_withbump(fighn, prim, prob, trainx,trainy)

clf(fighn);
testdata = linspace(prob.xl_bl, prob.xl_bu, 100);
testdata = testdata';

[freal, ~]= prob.evaluate_l([], testdata);
minf = min(freal);
plot(testdata, freal); hold on;

scatter(trainx, trainy, 30, 'ro', 'filled');

for ii = 1:2
    
    fprintf('box %d, mean %f, support %d \n', ii, prim.boxes{ii}.mean, prim.boxes{ii}.supportN);
    % adjust boundary
    mentioned_var = prim.boxes{ii}.vars;
    bump_lb = prob.xl_bl;
    bump_ub = prob.xl_bu;
    
    nb = length(mentioned_var);
    
    if ~isempty(mentioned_var)
        for jj = 1:nb % varible indicator
            if ~isnan(prim.boxes{ii}.min(jj))
                bump_lb(mentioned_var(jj)) =  prim.boxes{ii}.min(jj);
            end
            
            if ~isnan(prim.boxes{ii}.max(jj))
                bump_ub(mentioned_var(jj)) =  prim.boxes{ii}.max(jj);
            end
            
        end
    else
        continue;
    end
    
    
    % find corresponding y
    bump_ly = prob.evaluate_l([], bump_lb);
    bump_uy = prob.evaluate_l([], bump_ub);
    
    if bump_ly > bump_uy
        tmp = bump_uy;
        bump_uy = bump_ly;
        bump_ly = tmp;
    end
    
    xx = [bump_lb, bump_ub,bump_ub, bump_lb];
    yy = [minf, minf, bump_uy, bump_uy];
    
    
    %
    fill(xx, yy, 'k', 'FaceAlpha', 0.1, 'EdgeColor','k');
    pause(1);
end
end

function[best_x, info] = localBeliever_gpr(train_x, krg_obj, krg_con, ...
    num_pop, num_gen, varargin)
% fitness is an useless position, but cannot be deleted at this stage

daceflag = varargin{1};
prob =  varargin{2};
lb = varargin{3};
ub = varargin{4};


funh_obj = @(x)surrogate_predict(x, krg_obj, daceflag);
funh_con = @(x)surrogate_predict(x, krg_con, daceflag);

param.gen=num_gen;
param.popsize = num_pop;

l_nvar = size(train_x, 2);
[~,~,~, archive] = gsolver(funh_obj, l_nvar,  lb, ub, [], funh_con, param);
[newx, growflag] = believer_select(archive.pop_last.X, train_x, prob, false, true);

if growflag % there is unseen data in evolution
    [best_x] = surrogate_localsearch(newx, lb, ub, krg_obj, krg_con, daceflag);
    % inprocess_plotsearch(fighn, prob, cons_hn, new_xl, train_xl);
    
else % there is no unseen data in evolution
    % re-introduce random individual
    fprintf('In believer, no unseen data in last population, introduce randomness \n');
    best_x = [];
end


info = [];
end

function[newx] = surrogate_localsearch(newx, lb, ub, krg_obj, krg_con, daceflag)
%----------------------------

if length(krg_obj) > 1
    % no local search for MO
    return
end

funh_obj = @(x)surrogate_predict(x, krg_obj, daceflag);
funh_con = @(x)surrogate_predict(x, krg_con, daceflag);

opts = optimset('fmincon');
opts.Algorithm = 'sqp';
opts.Display = 'off';
opts.MaxFunctionEvaluations = 100;
[newx, newf, ~, output] = fmincon(funh_obj, newx, [], [],[], [],  ...
    lb, ub, [],opts);

end

