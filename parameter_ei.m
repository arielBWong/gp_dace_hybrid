param_ei= struct();
param_ei.num_pop = 100;
param_ei.num_gen = 100;
param_ei.init_size = 32;
param_ei.iter_size = 168; 
param_ei.propose_nextx = 'EIMnext';
param_ei.norm_str = 'normalization_z';
param_ei.llfit_hn = 'EIM_evaldace';
save('param_ei', 'param_ei');