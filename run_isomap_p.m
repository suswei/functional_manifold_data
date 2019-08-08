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

        fnm = fullfile(S(k).name,'possible_K.txt');  
        file3 = fopen(fnm,'r') ;
        num_neigh = fscanf(file3,'%f');
        fclose(file3);

        fnm = fullfile(S(k).name,'true_geo.txt');  
        file2 = fopen(fnm,'r') ;
        geo_temp = fscanf(file2,'%f');
        fclose(file2);
        true_geo=reshape(geo_temp,n,n);

        fnm = fullfile(S(k).name,'discrete_data.txt');  
        file = fopen(fnm,'r') ;
        X_temp = fscanf(file,'%f');
        fclose(file);
        data = reshape(X_temp,K,n)';

        fnm = fullfile(S(k).name,'possible_delta.txt');  
        file = fopen(fnm,'r') ;
        delta_can = fscanf(file,'%f');
        fclose(file);

        Dis= L2_distance(data',data',1); %*sqrt(range(grid)/(K-1));
        norm_true_geo = norm(true_geo,'fro');

        err=zeros(length(delta_can),length(num_neigh));
        for del=1:length(delta_can)
            for nei=1:length(num_neigh)
                Manidis = penalized_Isomap(Dis,num_neigh(nei),delta_can(del));
                err(del,nei)=norm(Manidis -true_geo,'fro')/norm_true_geo;
            end
        end

        [row_min,col_min]=find(err == min(err(:)));
        row_min = min(row_min);
        col_min = min(col_min);
        sprintf('delta_chosen=%d',delta_can(row_min))

        % Manidis = penalized_Isomap (Dis,num_neigh(col_min),delta_can(row_min));
        %
        % fichier = fopen('./err_delta_neigh.txt','w');
        % fprintf(fichier,'%12.8f\n',err);
        % fclose(fichier);
        %
        % fichier = fopen('./manidis.txt','w');
        % fprintf(fichier,'%12.8f\n',Manidis);
        % fclose(fichier);
    catch
        
    end
end