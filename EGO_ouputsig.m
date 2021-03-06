seedmax = 21;
median_num = 11;

problems = { 'smd1(3, 3)',  'smd2(3, 3)','rosenbrock(3, 3)', 'Zakharov(3, 3)', ...
             'levy(3, 3)','ackley(3, 3)', 'smd3(3, 3)', 'smd4(3, 3)' , ...
               'dsm1(3, 3)', 'tp7(3, 3)','tp9(3, 3)', 'Shekel(3, 3)',...
                'tp3(3, 3)', 'tp6(3, 3)','rastrigin(3, 3)',  'tp5(3, 3)'};    
        
algos = { 'llble_gp', 'lleim_gp','lladp_bh4'}; % , 'lleim_gp', 'lladp_gp'
legs = {'BLE', 'EIM', 'BHEIM'};
colors = {'b', 'r', 'k'};
legend_algos = {'BLE', 'EIM', 'BHEIM'};
%--construct result save path
plotpath = strcat(pwd, '\plots');
np = length(problems);  
% ns = length(seeds);
na = length(algos);
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
            fout_file = strcat(fout_folder, '\fl_', num2str(jj), '.csv' );
            
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

end

%---significant test on KB and BH
median_overproblems = zeros(np, na);
median_seed = zeros(np, na);

sigTest = zeros(np, 3);
for ii = 1:np
 % KB-EI
  KB = best_fl{ii}(1, :);
  EI = best_fl{ii}(2, :);
  [p,h] = ranksum(KB,EI, 'tail','left');
  
  sigTest(ii, 1) = h;
  
  [p,h] = ranksum(KB,EI, 'tail','right');
  sigTest(ii, 2) = h;
  
  [p,h] = ranksum(KB,EI);
  sigTest(ii, 3) = h;
end


savename = strcat(plotpath,'\median_sigKB2EI.csv');
fp=fopen(savename,'w');
fprintf(fp, 'problem_method, ');
testtype = {'left', 'right', 'both'};
for jj = 1:na
    fprintf(fp,'%s, ', testtype{jj} );
end
fprintf(fp,'\n');
for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
    name = strcat( prob.name, '_', num2str(prob.n_lvar));
    fprintf(fp, '%s,', name);
    for jj = 1:3
       st = num2str(sigTest(ii, jj));
       fprintf(fp, '%s,', st);
    end
    fprintf(fp, '\n');
end
fclose(fp);


sigTest = zeros(np, 3);
for ii = 1:np
 % KB-EI
  KB = best_fl{ii}(1, :);
  BH = best_fl{ii}(3, :);
  [p,h] = ranksum(KB,BH, 'tail','left');
  
  sigTest(ii, 1) = h;
  
  [p,h] = ranksum(KB,BH, 'tail','right');
  sigTest(ii, 2) = h;
  
  [p,h] = ranksum(KB,BH);
  sigTest(ii, 3) = h;
end

savename = strcat(plotpath,'\median_sigKB2BH.csv');
fp=fopen(savename,'w');
fprintf(fp, 'problem_method, ');
testtype = {'left', 'right', 'both'};
for jj = 1:na
    fprintf(fp,'%s, ', testtype{jj} );
end
fprintf(fp,'\n');
for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
    name = strcat( prob.name, '_', num2str(prob.n_lvar));
    fprintf(fp, '%s,', name);
    for jj = 1:3
       st = num2str(sigTest(ii, jj));
       fprintf(fp, '%s,', st);
    end
    fprintf(fp, '\n');
end
fclose(fp);

%--------------

sigTest = zeros(np, 3);
for ii = 1:np
 % EI-BH
  EI = best_fl{ii}(2, :);
  BH = best_fl{ii}(3, :);
  [p,h] = ranksum(EI,BH, 'tail','left');
  
  sigTest(ii, 1) = h;
  
  [p,h] = ranksum(EI,BH, 'tail','right');
  sigTest(ii, 2) = h;
  
  [p,h] = ranksum(EI,BH);
  sigTest(ii, 3) = h;
end

savename = strcat(plotpath,'\median_sigEI2BH.csv');
fp=fopen(savename,'w');
fprintf(fp, 'problem_method, ');
testtype = {'left', 'right', 'both'};
for jj = 1:na
    fprintf(fp,'%s, ', testtype{jj} );
end
fprintf(fp,'\n');
for ii = 1:np
    prob = problems{ii};
    prob = eval(prob);
    name = strcat( prob.name, '_', num2str(prob.n_lvar));
    fprintf(fp, '%s,', name);
    for jj = 1:3
       st = num2str(sigTest(ii, jj));
       fprintf(fp, '%s,', st);
    end
    fprintf(fp, '\n');
end
fclose(fp);




