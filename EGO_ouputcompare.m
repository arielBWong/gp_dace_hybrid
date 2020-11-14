seedmax = 21;
median_num = 11;
problems ={'ackley(3, 3)', 'levy(3, 3)','rastrigin(3, 3)','dsm1(3, 3)', ... %  multimodal global structure  heavy modality and weak modality
            'tp3(3, 3)', 'tp5(3, 3)', 'tp7(3, 3)','Shekel(3, 3)', ... % multimodal no global structure
            'Zakharov(3, 3)', 'smd2(3, 3)',  'rosenbrock(3, 3)'};
algos = { 'llble_gp', 'lleim_gp'};
colors = {'b', 'r'};

np = length(problems);  
% ns = length(seeds);
na = length(algos);

%--construct result save path
plotpath = strcat(pwd, '\plots');



% process each problem
for ii = 1: np
    prob = eval(problems{ii});
    nvar = prob.n_lvar;
    ninit = 11 * nvar - 1;
    
    %collect results for each algo
    algos_out = cell(1, na);
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
    title(strcat(prob.name,' k ', num2str(nvar)), 'FontSize', 16);
    legend(algos{1}, strcat(algos{1}, ' std'),  algos{2}, strcat(algos{2}, ' std'),  'FontSize', 14);
    
    % save to plot folder
    
    savename = strcat(plotpath,'\', prob.name,'_', num2str(nvar), '_converge.fig');
    saveas(figure(ii), savename);
    close(figure(ii));
end
