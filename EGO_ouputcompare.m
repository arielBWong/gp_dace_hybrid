seedmax = 21;
median_num = 11;

problems = {'ackley(3, 3)', 'levy(3, 3)','rastrigin(3, 3)','dsm1(3, 3)', ... %  multimodal global structure  heavy modality and weak modality
            'tp3(3, 3)', 'tp5(3, 3)', 'tp7(3, 3)','Shekel(3, 3)', ... % multimodal no global structure
            'Zakharov(3, 3)', 'smd2(3, 3)',  'rosenbrock(3, 3)',...
            'smd1(3, 3)', 'smd3(3, 3)', 'smd4(3, 3)' ,  'tp6(3, 3)', 'tp9(3, 3)'};
        
algos = { 'llble_gp', 'lleim_gp','lladp_bh4'}; % , 'lleim_gp', 'lladp_gp'
legs = {'BLE', 'EIM', 'BHEIM'};
colors = {'b', 'r', 'k'};
legend_algos = {'BLE', 'EIM', 'BHEIM'};


np = length(problems);  
% ns = length(seeds);
na = length(algos);

%--construct result save path
plotpath = strcat(pwd, '\plots');


best_fl = cell(1,np);

% process each problem
for ii = 1: np
    prob = eval(problems{ii});
    nvar = prob.n_lvar;
    ninit = 11 * nvar - 1;
    
    %collect results for each algo
    algos_out = cell(1, na);
    best_fl{ii} = zeros(na, seedmax);

    for kk = 1:na
      
        algos_out{kk} = [];
        for jj = 1:seedmax
            % save folder and save name
            fout_folder = strcat(pwd, '\resultfolder_gp\', prob.name, '_', num2str(nvar), '_', algos{kk}, '_init_', num2str(ninit));
            fout_file = strcat(fout_folder, '\fl_', num2str(jj), '.csv' )
            
            fl = csvread(fout_file);
            fl_size = length(fl);
            fl_part = [];
            
            % min f so far
            for mm = ninit : fl_size
                tmp = min(fl(1: mm));
                fl_part = [fl_part, tmp];
            end           
            algos_out{kk} = [algos_out{kk}; fl_part];
            best_fl{ii}(kk, jj) = fl_part(end);
        end
    end
    
    
    % generate plot: median fl from start from 32 to last 
    figure(ii);
    % plot each algo output
    for kk = 1:na
        algo_mean = mean(algos_out{kk}, 1);
        plot(algo_mean, colors{kk}, 'LineWidth', 2); hold on;
        
        % std
        algo_std = std(algos_out{kk},1);
        
        x = 1:length(algo_mean);
        x = [x, fliplr(x)];
        y1 = algo_mean + algo_std;
        y2 = algo_mean - algo_std;
        y = [y1, fliplr(y2)];
        fill(x, y, colors{kk}, 'FaceAlpha', 0.1, 'EdgeColor','none');
        
    end
    title(prob.name, 'FontSize', 16);
    legend(legs{1}, strcat(legs{1}, ' std'),  legs{2}, strcat(legs{2}, ' std'), legs{3}, strcat(legs{3}, ' std'),  'FontSize', 14);
    
    % save to plot folder   
    savename = strcat(plotpath, '\', prob.name,'_', num2str(nvar), '_converge.fig');
    saveas(figure(ii), savename);
    close(figure(ii));
end

%--- median to csv
median_overproblems = zeros(np, na);
median_seed = zeros(np, na);
for ii = 1:np
    for jj = 1:na
        median_overproblems(ii, jj) = median(best_fl{ii}(jj, :));
        one_record = best_fl{ii}(jj, :);
        % ----
        [~, id] = sort(one_record);
        median_seed(ii, jj) = id(median_num);        
    end
end
% 
% %--create median plot---
for ii = 1:np
    figure(ii);
    prob = eval(problems{ii});
    nvar = prob.n_lvar;
    ninit = 11 * nvar - 1;
    
    for jj =  1:na
        % retrieve median run
        fout_folder = strcat(pwd, '\resultfolder_gp\', prob.name, '_', num2str(nvar), '_', algos{jj}, '_init_', num2str(ninit));
        seed = median_seed(ii, jj);
        fout_file = strcat(fout_folder, '\fl_', num2str(seed), '.csv' );
        fl = csvread(fout_file);
        fl_part = [];
        for mm = ninit : fl_size
           tmp = min(fl(1: mm));
           fl_part = [fl_part, tmp];
        end 
            
        
        % plot this run
        plot(fl_part, colors{jj}, 'LineWidth', 3); hold on;      
    end
    
    % add rest legend things
    % save and close
    title(strcat(prob.name,' k ', num2str(nvar)), 'FontSize', 16);
    legend(legend_algos{1}, legend_algos{2}, legend_algos{3},  'FontSize', 14);
    
    % save to plot folder    
    savename = strcat(plotpath,'\', prob.name,'_', num2str(nvar), '_meidan_converge.fig');
    saveas(figure(ii), savename);
    close(figure(ii));
    
end
% %-----------------------


savename = strcat(plotpath,'\median_bestf.csv');
fp=fopen(savename,'w');
fprintf(fp, 'problem_method, ');
for jj = 1:na
    fprintf(fp,'%s, ', algos{jj} );
end
fprintf(fp,'\n');
for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
    name = strcat( prob.name, '_', num2str(prob.n_lvar));
    fprintf(fp, '%s,', name);
    for jj = 1:na
       st = num2str(median_overproblems(ii, jj));
       fprintf(fp, '%s,', st);
    end
    fprintf(fp, '\n');
end
fclose(fp);


% %---- generate profiling
min_alg = min(median_overproblems, [], 2);
rangelist = max(median_overproblems,  [], 2) ./ min(median_overproblems, [], 2);
range = max(rangelist);

Data1 = median_overproblems';
PerformanceProfile(Data1,legs);

n = 100; % 
tao_list = linspace(1, range, n);
% plot rows
plotrows = zeros(na, n);
for ii = 1:n
    for jj = 1:na
        count = 0;
        for kk = 1:np
            % ---------
            performance_p = median_overproblems(kk, jj);
            if performance_p < 0
                r = (1/performance_p)/(1/min_alg(kk));
            else
                r = performance_p/min_alg(kk);
            end
            
            if r <= tao_list(ii)
                count = count + 1;
            end           
        end       
        plotrows(jj, ii) = count/np;
    end
end
mkr = {'--s', '--o', '--d'};
for ii = 1:na
    plot(log(tao_list), plotrows(ii, :), mkr{ii}, 'MarkerSize',10,  'LineWidth',2); hold on;
end
 legend( legend_algos{1}, legend_algos{2}, legend_algos{3},'FontSize', 14);
 title('Profiling',  'FontSize', 16);
 xlabel('log(\tau)',  'FontSize', 14)
 ylabel('p(\tau)',  'FontSize', 14)
 savename = strcat(plotpath,'\', prob.name,'_', num2str(nvar), '_profiling.fig');
saveas(figure(ii), savename);

