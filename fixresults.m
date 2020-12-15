%% 
seedmax = 11;
problems = {'smd1(3, 3)', 'smd3(3, 3)', 'smd4(3, 3)' ,  'tp6(3, 3)', 'tp9(3, 3)'};
algos = {'lladp_bh'};

np = length(problems);  
% ns = length(seeds);
na = length(algos);
ei_end = 140 + 32;
middle = 201;
bk_end = 28;
% 
% for ii = 1: np
%     prob = eval(problems{ii});  
%      for kk = 1:na
%       
%         algos_out{kk} = [];
%         for jj = 1:seedmax
%             % save folder and save name
%             fout_folder = strcat(pwd, '\resultfolder_gp\', prob.name, '_', num2str(nvar), '_', algos{kk}, '_init_', num2str(ninit));
%             fout_file = strcat(fout_folder, '\fl_', num2str(jj), '.csv' );
%             
%             fl = csvread(fout_file);
%             fl_size = length(fl);            
%             
%             % output folder 
%             fout_folder = strcat(pwd, '\resultfolder_gp\', prob.name, '_', num2str(nvar), '_', algos{kk}, '_fixed_init_', num2str(ninit));
%             n = exist(fout_folder);
%             if n ~= 7
%                 mkdir(fout_folder)
%             end
% 
%             
%             fl_newei = fl(1: ei_end, :);
%             fl_newbk = fl(middle: middle + bk_end - 1, :);
%             
%             fl_new = [fl_newei; fl_newbk];
%             flmin = min(fl_new);
%             fl_new = [fl_new; flmin];
%             
%             fl_rearrange = [];
%             count = 1;
%             for mm = 33: 140 + 32
%                 indiv = fl_new(mm);
%                 fl_rearrange = [fl_rearrange; indiv];
%                 
%                 if mod(mm-32, 5) == 0 
%                     indiv = fl_new(168 + count);
%                     fl_rearrange = [fl_rearrange; indiv];
%                     count = count + 1;
%                 end
%             end 
%             
%             
%            %  arrange head and tail
%            indiv = fl_new(1: 32);
%            fl_rearrange = [indiv;fl_rearrange];
%            fl_rearrange = [fl_rearrange; fl_new(end)];
%            
%             
%             
%             
%             
%             
%             fout_file = strcat(fout_folder, '\fl_', num2str(jj), '.csv' );
%             csvwrite(fout_file, fl_rearrange);
%             
%         end
%     end
% end
% 
% for ii = 1: np
%     prob = eval(problems{ii});  
%      for kk = 1:na
%       
%         algos_out{kk} = [];
%         for jj = 12:21
%             % save folder and save name
%             fout_folder = strcat(pwd, '\resultfolder_gp\', prob.name, '_', num2str(nvar), '_', algos{kk}, '_init_', num2str(ninit));
%             fout_file = strcat(fout_folder, '\fl_', num2str(jj), '.csv' );
%             
%             fl = csvread(fout_file);
%             fl_size = length(fl);            
%             
%             % output folder 
%             fout_folder = strcat(pwd, '\resultfolder_gp\', prob.name, '_', num2str(nvar), '_', algos{kk}, '_fixed_init_', num2str(ninit));
%             n = exist(fout_folder);
%             if n ~= 7
%                 mkdir(fout_folder)
%             end
% 
%             fl_rearrange = [];
%             count = 1;
%             for mm = 33: 140 + 32
%                 indiv = fl(mm);
%                 fl_rearrange = [fl_rearrange; indiv];
%                 
%                 if mod(mm-32, 5) == 0 
%                     indiv = fl(168 + count);
%                     fl_rearrange = [fl_rearrange; indiv];
%                     count = count + 1;
%                 end
%             end 
%             
%            %  arrange head and tail
%            indiv = fl(1: 32);
%            fl_rearrange = [indiv;fl_rearrange];
%            fl_rearrange = [fl_rearrange; fl(end)];
%            
%            % save to 
%            fout_file = strcat(fout_folder, '\fl_', num2str(jj), '.csv' );
%            csvwrite(fout_file, fl_rearrange)
%         end
%     end
% end

for ii = 1: np
    prob = eval(problems{ii});
    for jj = 1:21
        
        % save folder and save name
        fout_folder = strcat(pwd, '\resultfolder_gp\', prob.name, '_', num2str(nvar), '_', algos{kk}, '_init_', num2str(ninit));
        fout_file = strcat(fout_folder, '\fl_', num2str(jj), '.csv' );   
        delete(fout_file);
        
        % output folder
        from_folder = strcat(pwd, '\resultfolder_gp\', prob.name, '_', num2str(nvar), '_', algos{kk}, '_fixed_init_', num2str(ninit));
        from_file = strcat(from_folder, '\fl_', num2str(jj), '.csv');
        copyfile(from_file, fout_file);
    end
end