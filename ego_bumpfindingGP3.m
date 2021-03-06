function[match_xl, n_fev, flag] = ego_bumpfindingGP3(xu, prob, seed)
% based on EI 
% local believer considers that local range should include minf

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

fprintf('problem  %s, seed %d \n', prob.name, seed);
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
    % fprintf('generation %d \n', g);
    
    [krg, krgc] = surrogate_create(train_xl, train_fl, train_fc, normhn, daceflag);
    % disconnected_bump = bump_detection_disconnectedRegion(train_xl, train_fl, prob.xl_bl, prob.xl_bu);
    
    % switch_record = [switch_record, 1];
    [new_xl, ~] = switchEIM_gpr(train_xl, train_fl, prob.xl_bu, prob.xl_bl, ...
        num_pop, num_gen, train_fc, fithn, normhn, krg, krgc);
       
    % plot 2d with bump hunting
    if size(train_xl, 2) == 2
        % processplot2d_bumpstudy(fighn, train_xl, train_fl, krg, prob, initx, new_xl, daceflag)
    end

    % local search
    trainx_all = [train_xl; local_ax];
    trainy_all = [train_fl; local_ay];
    
    % identify minf
    [miny, id_min] = min(trainy_all, [], 1);
     minx = trainx_all(id_min, :);
    
    prim =  bump_detection_4test(trainx_all, trainy_all, prob.xl_bl, prob.xl_bu);
    nboxes = length(prim.boxes);
    localsurrogate = false;
    % select local region, and make sure fmin is included in the region
    if nboxes > 2 && mod(g, 5) == 0  % local search frequency
        upper_stack = [];
        lower_stack = [];
        
        % default select two regions 
        for ii = 1:2
            % fprintf('box %d, mean %f, support %d \n', ii, prim.boxes{ii}.mean, prim.boxes{ii}.supportN);
            % adjust boundary
            [bump_lb, bump_ub] = boxboundary(prim, ii, prob);
            upper_stack = [upper_stack; bump_ub];
            lower_stack = [lower_stack; bump_lb];
        end
        
        upper_bound = max(upper_stack, [], 1);
        lower_bound = min(lower_stack,[],  1);
        
        % check whether fmin is in the region
        % inboxflag = inboxcheck(minx, upper_bound, lower_bound);
        if ~inboxcheck(minx, upper_bound, lower_bound)
            
            localsurrogate = false;
            % expand upper bound and lower bound until fmin is within range
            for ii = 3:nboxes-1               
                [bump_lb, bump_ub] = boxboundary(prim, ii, prob);
                % include new boundary
                upper_stack = [upper_stack; bump_ub];
                lower_stack = [lower_stack; bump_lb];
                upper_bound = max(upper_stack, [], 1);
                lower_bound = min(lower_stack,[],  1);
                % check whether new boundary cover fmin 
                if inboxcheck(minx, upper_bound, lower_bound)
                    localsurrogate = true;
                    break
                end
            end
        else
            localsurrogate = true;
            
        end
       
         
        if localsurrogate
            % run local search
            [localx, localy] = pickup_localxy(upper_bound, lower_bound, trainx_all, trainy_all);
            kk = size(localx, 1);
            % fprintf('Local training data size %d \n', kk);
            % fprintf('Archive min is %f \n', miny);
            
            % build local surrogate
            [local_krg, local_krgc] =  surrogate_create(localx, localy, [] , normhn, daceflag);
            
            % run believer to find next point
            [local_newx, ~] = localBeliever_gpr(localx, local_krg, local_krgc, ...
                num_pop, num_gen, daceflag, prob, lower_bound, upper_bound);
            % plot local search
            if size(train_xl, 2) == 1
                plot1d_global_withbump(fighn, prim, prob, train_xl,train_fl);
                plotlocal_search(fighn, upper_bound, lower_bound, localy, local_krg);
            end
            
            if size(train_xl, 2) == 2
                plot2d_global_withbump( prim, prob, local_krg, localx, localy, minx, local_newx);
            end
            
            [local_newy,~] = prob.evaluate_l(xu, local_newx);
            
            local_ax = [local_ax; local_newx];
            local_ay = [local_ay; local_newy];
        end
        
    end
    
%      if size(train_xl, 2) == 1
%         processplot1d_bumpstudy(fighn, trainx_all, trainy_all, krg, prob, initx, new_xl, daceflag)
%      end
%     
    
    
    
    mdl_save{g} = krg;
    
    % fprintf('problem %s, iter %d, seed %d, dace flag is %d \n', prob.name,g,  seed, daceflag);
    [new_fl, new_fc] = prob.evaluate_l(xu, new_xl);
    
    %--- closeness check---
%     check = abs(train_xl - new_xl);
%     check = round(check, 5);
%     if sum(check < 1e-20) >= 1
%         fprintf('fail unique check %d \n', g);
%         % continue;
%     end
    
    train_xl = [train_xl; new_xl];
    train_fl = [train_fl; new_fl];
    train_fc = [train_fc; new_fc];  % compatible with nonconstraint
    
   
    
end



% organise local search and EI search
train_xl = [train_xl; local_ax];
train_fl = [train_fl; local_ay];
train_fc = [train_fc; []];

[match_xl, n_fev, flag] = post_infillsaveprocess(xu, train_xl, train_fl, train_fc, 'lladp_bh', seed, prob, init_size);
% save_eachmodel(seed, 'lladp_gp', mdl_save, init_size, prob, switch_record);
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


function plot2d_global_withbump( prim, prob, local_krg, localx, localy, minx, newx)

% draw local search region box/ copy from main process
% default select two regions 
upper_stack = [];
lower_stack = [];

for ii = 1:2
    % fprintf('box %d, mean %f, support %d \n', ii, prim.boxes{ii}.mean, prim.boxes{ii}.supportN);
    % adjust boundary
    [bump_lb, bump_ub] = boxboundary(prim, ii, prob);
    upper_stack = [upper_stack; bump_ub];
    lower_stack = [lower_stack; bump_lb];
end

upper_bound = max(upper_stack, [], 1);
lower_bound = min(lower_stack,[],  1);

% check whether fmin is in the region 
% inboxflag = inboxcheck(minx, upper_bound, lower_bound);
if ~inboxcheck(minx, upper_bound, lower_bound)
    
    localsurrogate = false;
    % expand upper bound and lower bound until fmin is within range
    for ii = 3:nboxes-1
        [bump_lb, bump_ub] = boxboundary(prim, ii, prob);
        % include new boundary
        upper_stack = [upper_stack; bump_ub];
        lower_stack = [lower_stack; bump_lb];
        upper_bound = max(upper_stack, [], 1);
        lower_bound = min(lower_stack,[],  1);
        % check whether new boundary cover fmin
        if inboxcheck(minx, upper_bound, lower_bound)
            localsurrogate = true;
            break
        end
    end
else
    localsurrogate = true;    
end
% draw local search boundary, two plots,one for real function, one for
% local surrogate
% real plot
lb1 = lower_bound(1);
ub1 = upper_bound(1);

lb2 = lower_bound(2);
ub2 = upper_bound(2);

num_points = 101;
x1 = linspace(lb1, ub1, num_points);
x2 = linspace(lb2, ub2, num_points);
[x1, x2] = meshgrid(x1, x2);
testdata = zeros(num_points, num_points);

xx = [lb1, ub1, ub1, lb1, lb1];
yy = [lb2, lb2, ub2, ub2, ub1];


for i = 1:num_points
    for j = 1:num_points
        [f(i, j), ~] = prob.evaluate_l([], [x1(i, j), x2(i, j)]);
        [fs(i, j), ~] = surrogate_predict([x1(i, j), x2(i, j)], local_krg, false);
        fs(i, j) = denormzscore( localy,  fs(i, j));
    end
end

f1 = figure(1);
surfc(x1, x2, f); hold on;
% contour(x1, x2, f); hold on;
xlim([prob. xl_bl(1), prob.xl_bu(1)]);
ylim([prob.xl_bl(2), prob.xl_bu(2)]);
colormap jet
shading interp
fill(xx, yy, 'k', 'FaceAlpha', 0.1, 'EdgeColor','k');

% surrogate plot
f2 = figure(2);
surfc(x1, x2, fs); hold on;
% contour(x1, x2, fs); hold on;
xlim([prob.xl_bl(1), prob.xl_bu(1)]);
ylim([prob.xl_bl(2), prob.xl_bu(2)]);
colormap jet
shading interp
fill(xx, yy, 'k', 'FaceAlpha', 0.1, 'EdgeColor','k');

newf = surrogate_predict(newx, local_krg, false);
newf = denormzscore( localy,  newf);

minf = surrogate_predict(minx, local_krg, false);
minf = denormzscore( localy,  minf);

scatter3(newx(1), newx(2), newf, 80, 'g', 'filled');
scatter3(minx(1), minx(2), minf, 80, 'r', 'filled');
scatter3(localx(:, 1), localx(:, 2), localy, 80, 'kx');

pause(2);
close(f1);
close(f2);
a = 0;
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

plottrue = 1;
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

function [inbox_flag] = inboxcheck(x,  bu, lb)
    inbox_flag = false;
    nx = size(x,2); 
   if sum(bu - x >= 0 ) == nx && sum(x - lb >= 0)
       inbox_flag = true;
   end
    
end

function[bump_lb, bump_ub] = boxboundary(prim, ii, prob)
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
end
