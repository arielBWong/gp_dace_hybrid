function[match_xl, n_fev, flag] = ego_HybridGP(xu, prob, seed)

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

global eps_dist 
eps_dist = sqrt(l_nvar) * par;  %  1% of max normalizated distance



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
    fighn = figure(2);
else
    fighn = [];
end

normhn = str2func(norm_str);
mdl_save = cell(1, iter_size);
fithn = str2func('EIM_eval');
switchrecord = []; % 1- KB, 0 -EI

for g = 1: iter_size   
     
    daceflag = false;
    [new_xl, infor] = Believer_next(train_xl, train_fl, upper_bound, lower_bound, ...
        num_pop, num_gen, train_fc, [], normhn, daceflag, prob);
    
    switchrecord(g) = 1;
 
    tooclose = archive_check(new_xl, train_xl, prob);
    if tooclose
        % ---- determine whether to switch to eim -----
        fprintf('adopt eim at iter: %d\n', g); 
        [new_xl, infor] = switchEIM_gpr(train_xl, train_fl, upper_bound, lower_bound, ...
            num_pop, num_gen, train_fc, fithn, normhn, infor.krg, infor.krgc);
        switchrecord(g) = 0;
             
    end
    
    
   if size(train_xl, 2) == 1
        processplot1d(fighn, train_xl, train_fl, infor.krg, prob, initx, new_xl, daceflag)
    end 
    mdl_save{g} = infor.krg;
    
    
    
    % fprintf('problem %s, iter %d, seed %d, dace flag is %d \n', prob.name,g,  seed, daceflag);
    [new_fl, new_fc] = prob.evaluate_l(xu, new_xl);
    
    %--- closeness check---
    check = abs(train_xl - new_xl);
    check = round(check, 5);
    if length(unique(check, 'rows')) < size(train_xl, 1)
        fprintf('fail unique check \n');
        % continue;
    end
    
    train_xl = [train_xl; new_xl];
    train_fl = [train_fl; new_fl];
    train_fc = [train_fc; new_fc];  % compatible with nonconstraint
 
end

[match_xl, n_fev, flag] = post_infillsaveprocess(xu, train_xl, train_fl, train_fc, 'lladp_gp', seed, prob, init_size);
save_eachmodel(seed, 'lladp_gp', mdl_save, init_size, prob, switchrecord);
end


%-----------paper demo-------------


% ---paper demo-
function[] = processplot(fighn, trainx, trainy, krg, prob, initx)
clf(fighn);

% (1) create test
testdata = linspace(prob.xl_bl, prob.xl_bu, 2000);
testdata = testdata';

% (2) predict
[fpred, sig] = surrogate_predict(testdata, krg);
fpred = denormzscore(trainy, fpred);
plot(testdata, fpred, 'r--', 'LineWidth', 2); hold on;

y1 = fpred + sig * 1.5;
y2 = fpred - sig * 1.5;
y = [y1', fliplr(y2')];
x = [testdata', fliplr(testdata')];
fill(x, y, 'r', 'FaceAlpha', 0.1, 'EdgeColor','none'); hold on;

% (3) real
[freal, sig]= prob.evaluate_l([], testdata);
plot(testdata, freal, 'b', 'LineWidth', 2);hold on;


% (4) scatter train
scatter(trainx, trainy, 80, 'ko', 'LineWidth', 2);

inity = prob.evaluate_l([], initx);
scatter(initx, inity, 40, 'ro', 'filled');


pause(1);

end

 function f = denormzscore(trainy, fnorm)
[train_y_norm, y_mean, y_std] = zscore(trainy);
f = fnorm * y_std + y_mean;
 end
 
 
 
 
function tooclose = archive_check(newx, trainx, prob)
% ---check newx whether it is
tooclose = false;
global eps_dist

upper_bound = prob.xl_bu;
lower_bound = prob.xl_bl;

trainx_norm = (trainx - lower_bound) ./ (upper_bound - lower_bound);
newx_norm = (newx - lower_bound) ./ (upper_bound - lower_bound);

%--- 
mindistance = min(pdist2(newx_norm,trainx_norm));

if mindistance < eps_dist
     tooclose =  true;
end
end