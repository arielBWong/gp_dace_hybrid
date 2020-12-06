function[match_xl, n_fev, flag] = ego_trustregiontest(xu, prob, seed)
% visualize trust region range
% local believer considers that local range should include minf

load('param_hyb');
rng(seed, 'twister');


num_pop = param_hyb.num_pop;
num_gen = param_hyb.num_gen;
init_size = param_hyb.init_size;
iter_size = param_hyb.iter_size;
norm_str = param_hyb.norm_str;
par = param_hyb.par;


max_eval = init_size + iter_size;
num_regions = 2;

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

% divided into num_regions
[indr, train_xl, train_fl] = init_regions(num_regions, train_xl, train_fl, lower_bound, upper_bound);

if size(train_xl, 2) == 1
    fighn = figure(1);
    testdata = linspace(prob.xl_bl, prob.xl_bu, 100);
    testdata = testdata';   
    [freal, ~]= prob.evaluate_l([], testdata);
    
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

length = zeros(num_regions, 1);
global regionlenmax
regionlenmax = min(prob.xl_bu - prob.xl_bl);
regionlen = ones(num_regions, l_nvar) * -1;


%--------------- 
num_eval = init_size;
while num_eval <= max_eval
    % fprintf('generation %d \n', g);
    
    candidate_x = [];
    candidate_y = [];
    candidate_id = [];
    for ii = 1:num_regions
        % extract regions
        regional_id = (indr(:, ii) == ii);
        
        regional_x = train_xl(regional_id, :);
        regional_f = train_fl(regional_id, :);
        
        subplot(1, 2, ii);
        plot(regional_x, regional_f, 'ro' ); hold on;
        plot(testdata, freal); hold on;
        plot(train_xl, train_fl, 'k+');
        
        ub = max(regional_x, [], 1);
        lb = min(regional_x, [], 1);
        length = ub - lb;
        
        for jj = 1:l_nvar 
         if regionlen(ii, jj) == -1
             regionlen(ii, jj) = length(jj);
         end
        end
        
        length = regionlen(ii);
        
        
%          n = 5;
%         [regional_x, regional_f, krg, krgc] = localEIM(regional_x, regional_f,  [], normhn, daceflag, n,  lb, ub, fithn, prob);
%         
%         
%         new_indr = ones(n, num_regions) * -1;
%         new_indr(:, ii) = ii;
%         indr = [indr; new_indr];
%         train_xl = [train_xl; regional_x(end - n + 1 : end, :)];
%         train_fl = [train_fl; regional_f(end - n + 1 : end, :)];
        


        % adjust regional index 
        [regional_fmin, fminid] = min(regional_f, [], 1);
        regional_xmin = regional_x(fminid, :); 
        regional_lb = regional_xmin - length ./ 2;
        regional_ub = regional_xmin + length ./ 2;
        [localx, localy, indr] = pickup_localxy(regional_ub, regional_lb, train_xl, train_fl, indr, ii);
        
        plot(localx, localy, 'b*'); hold on;
        
        
   
%         % min f and find search region
%         lengthscale = krg{1}.KernelInformation.KernelParameters(1:end-1);
%         nn = size(lengthscale, 1);
%         weights = lengthscale ./ mean(lengthscale);
%         weights = weights ./ prod( weights.^ (1.0/nn) ); %  We now have weights.prod() = 1
%         
%         % identify regional min
%         [regional_fmin, minid] = min(regional_f);
%         regional_xmin = regional_x(minid, :); 
%         
%         % identify local search region      
%         regional_lb = regional_xmin - length .* weights / 2;
%         regional_ub = regional_xmin + length .* weights / 2;
%         
%         % recreate local surrogate 
%         [localx, localy] = pickup_localxy(regional_ub, regional_lb, train_xl, train_fl);
%         plot(localx, localy, 'bo'); hold on;
        [krg, krgc] = surrogate_create(localx, localy, [], normhn, daceflag);
        
        if size(train_xl, 2) == 1
            fmax = max(regional_f);
            fmin = min(regional_f);
            xx = [regional_lb, regional_ub, regional_ub, regional_lb, regional_lb];
            yy = [fmin, fmin, fmax, fmax, fmin];
            fill(xx, yy, 'k', 'FaceAlpha', 0.05, 'EdgeColor','k'); hold on;
            
            %--prediction
            ynorm = surrogate_predict(testdata, krg, false);
            ydenorm = denormzscore(localy, ynorm);
            plot(testdata, ydenorm, '--');
                 
        end
        
        % apply KB to add to regional group
        [regional_newx, ~] = localBeliever_gpr(localx, krg, krgc, ...
                                                num_pop, num_gen, daceflag, prob,regional_lb, regional_ub);
         
        [regional_newyPred, ~] = surrogate_predict(regional_newx,krg, daceflag);  
        regional_newyPred = denormzscore(localy, regional_newyPred);
        
        plot(regional_newx, regional_newyPred, 'ko');
        
        candidate_x = [candidate_x; regional_newx];
        candidate_y = [candidate_y; regional_newyPred];    
        candidate_id = [candidate_id; ii];
    end
    
    
    % select infill from candidate
    [ymin, minid] = min(candidate_y);
    infillx = candidate_x(minid, :);
    infilly = prob.evaluate_l(xu, infillx);
    
    
    % adjust local region
    for ii = 1:num_regions
        if ii == minid % region to be adjusted
            regionlen = length_adjust(ii, infilly, train_xl, train_fl, indr, regionlen);
        end 
    end
    
    % expand archive
    train_xl = [train_xl; infillx];
    train_fl = [train_fl; infilly];
    
    % expand indr
    new_indr = ones(1, num_regions) * -1;
    new_indr(minid) = minid;
    indr = [indr; new_indr];
   
    
    % restart if necessary    
    for ii = 1:num_regions
        % if regionlen(ii) < lengthmin
         if restart_check(indr, trainx,regionlen , region_id)  
            % delete previous region index
            reg2del = indr(:, ii) == ii;
            indr(reg2del) = -1;
            % restart somewhere centered somewhere by EIM 
            % build global model, eim propose next point
            [krg_global, krgc_global] = surrogate_create(train_xl, train_fl, [], normhn, daceflag);
            [new_xk, ~] = switchEIM_gpr(train_xl, train_fl, prob.xl_bu, prob.xl_bl, ...
                                num_pop, num_gen, train_fc, fithn, normhn, krg, krgc);
    
              
            new_yk = prob.evaluate_l(xu, new_xk);
            
            % create new region around new_xk
            trainx_tmp = []; 
            trainy_tmp = [];
            
            newregionk = k; %
            dist = pdist2(new_xk, train_xl);
            [~, id] = sort(dist);
            id_newr = id(1:newregionk);
            
            % assign id_newr to deleted 
            indr(idd_newr, ii) = ii;   
            
            % add newxk newfk into train_xl
            train_xl = [train_xl; new_xk];
            train_fl = [train_fl; new_fk];
            tmp = ones(1, num_regions) * -1;
            tmp(1, ii) = ii;
            indr = [indr; tmp];
            
            
        end
    end
end


[match_xl, n_fev, flag] = post_infillsaveprocess(xu, train_xl, train_fl, [], 'llregional', seed, prob, init_size);
end

function restart_flag = restart_check(indr, trainx,regionlen , region_id)
restart_flag = false;
regionx_id = indr(:, region_id) == region_id;
regionx = trainx(regionx_id, :);
ub = regionx + regionlen(region_id, :)./2;
lb = regionx - regionlen(region_id, :) ./2;

count = 0;
[numx, nx] = size(trainx);
for ii = 1: numx
    x = trainx(ii, :);
    % with in boundary
    if sum(ub - x >= 0 ) == nx && sum(x - lb >= 0) == nx
        count = count + 1;
    end
end

if count <= 4
    restart_flag = true;   
end

end

function regionlen =  length_adjust(regionid, infilly, trainx, trainf, indr, regionlen)
% length_adjust(ii, infilly, train_xl, train_fl, indr)
% condition on adjust length
regionalid = indr(:, regionid)==regionid;
regionalx = trainx(regionalid, :);
regionalf = trainf(regionalid, :);

if infilly <= min(regionalf)- abs(min(regionalf)) * 1e-3
    regionlen(regionid) = min([regionlen(regionid)*2, regionlenmax]);
else
    regionlen(regionid) = regionlen(regionid)/2;
end


end


function[regional_x, regional_f, krg,krgc] = localEIM(regional_x, regional_f,  regional_c, ...
                                                      normhn, daceflag, n, bl, ul, fithn, ...
                                                      prob)
load('param_hyb');
num_pop = param_hyb.num_pop;
num_gen = param_hyb.num_gen;


for i = 1:n
    [krg, krgc] = surrogate_create(regional_x, regional_f, regional_c, normhn, daceflag);
    [new_xl, ~] = switchEIM_gpr(regional_x, regional_f, ul, bl, ...
        num_pop, num_gen, regional_c, fithn, normhn, krg, krgc);
    [new_fl, new_fc] = prob.evaluate_l([], new_xl);
    
    regional_x = [regional_x; new_xl];
    regional_f = [regional_f; new_fl];

end

end


function[indr, train_x, train_y] = init_regions(num_regions, train_x, train_y, lb, ub)

[numx, nvar] = size(train_x);
indr = ones(numx, num_regions) * -1;

k = numx/num_regions;
centers = lhsdesign(num_regions,nvar,'criterion','maximin','iterations',1000);
centers = repmat(lb, num_regions, 1) ...
    + repmat((ub - lb), num_regions, 1) .* centers;

for i = 1: num_regions
    % re-order trainx
    c = centers(i, :);
    dist = pdist2(c, train_x);
    [~, id] = sort(dist);
    
    % use index 
    indr(id(1:k), i) = i;
end
a = 0;

end


function [localx, localy, indr] = pickup_localxy(upper_bound, lower_bound, xx, yy, indr, region_id)
localx = [];
localy = [];
numx = size(xx, 1);
nx = size(xx, 2);

% adjust region
indr(:, region_id) = -1;
for ii = 1: numx
    x = xx(ii, :);
    y = yy(ii, :);
    
    % with in boundary
    if sum(upper_bound - x >= 0 ) == nx && sum(x - lower_bound >= 0) == nx
        localx = [localx; x];
        localy = [localy; y];
        indr(ii, region_id) = region_id;
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
