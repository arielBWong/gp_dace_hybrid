para_blem = struct();
para_blem.num_pop = 100;
para_blem.num_gen = 100;
para_blem.init_size = 32;
para_blem.iter_size = 168; 
para_blem.propose_nextx = 'Believer_next';
para_blem.norm_str = 'normalization_z';
para_blem.llfit_hn = 'anything';
save('para_blem', 'para_blem');