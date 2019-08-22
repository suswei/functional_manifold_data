% takes all sim_functional_data outputs and runs automatic tuning for
% pisomap. records how many of the chosen deltas = 0
% first run createSimTextFiles_pisomap.R to create data directories, than
% run_isomap_p in matlab

S=dir('2019*')
for k = 1:numel(S)
    try
        fnm = fullfile(S(k).name,'n_K.txt');  
        file = fopen(fnm,'r') ;
        temp = fscanf(file,'%f');
        fclose(file);
        n=temp(1);
        K=temp(2);

        fnm = fullfile(S(k).name,'grid.txt');  
        file = fopen(fnm,'r') ;
        grid_temp = fscanf(file,'%f');
        fclose(file);
        grid_temp2=repmat(grid_temp',n,1);
        gridt=mat2cell(grid_temp2,repmat([1],1,n));
        T=gridt';

        fnm = fullfile(S(k).name,'discrete_data.txt');  
        file = fopen(fnm,'r') ;
        X_temp = fscanf(file,'%f');
        fclose(file);
        X_temp_2 = reshape(X_temp,K,n)';
        X_t = mat2cell(X_temp_2,repmat([1],1,n));
        X = X_t';
        
        
        [hat_delta,estim_pI] = mani_p_iso(X,T);
               
    
        sprintf('delta_chosen=%d',hat_delta)

        
    catch
        
    end
end