

file = fopen('./n_K.txt','r') ;
temp = fscanf(file,'%f');
fclose(file)
n=temp(1);
K=temp(2);

file3 = fopen('./possible_K.txt','r') ;
num_neigh = fscanf(file3,'%f');
fclose(file3)

file2 = fopen('./true_geo.txt','r') ;
geo_temp = fscanf(file2,'%f');
fclose(file2)
true_geo=reshape(geo_temp,n,n);

file = fopen('./discrete_data.txt','r') ;
X_temp = fscanf(file,'%f');
fclose(file)
data = reshape(X_temp,K,n)';

file = fopen('./possible_delta.txt','r') ;
delta_can = fscanf(file,'%f');
fclose(file)

Dis= L2_distance(data',data',1); %*sqrt(range(grid)/(K-1));
norm_true_geo = norm(true_geo,'fro');

err=zeros(length(delta_can),length(num_neigh));
for del=1:length(delta_can)
   for nei=1:length(num_neigh)
       Manidis = penalized_Isomap (Dis,num_neigh(nei),delta_can(del));
       err(del,nei)=norm(Manidis -true_geo,'fro')/norm_true_geo;
   end
end

[row_min,col_min]=find(err == min(err(:)));
row_min = min(row_min);
col_min = min(col_min);
Manidis = penalized_Isomap (Dis,num_neigh(col_min),delta_can(row_min));

fichier = fopen('./err_delta_neigh.txt','w');
fprintf(fichier,'%12.8f\n',err);
fclose(fichier);

fichier = fopen('./manidis.txt','w');
fprintf(fichier,'%12.8f\n',Manidis);
fclose(fichier);