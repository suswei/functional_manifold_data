function [Manidis] = penalized_Isomap (Dis,K,delta)
    N = size(Dis,1);    
    % initial step: connect K nearest points
    [Sort_dis,Ord] = sort(Dis);
    radius = median(Sort_dis(K+1,:));
    Counts = [];
    Indmat = zeros(N,N);
    for i=1:N
        Counts(i) = sum(Dis(:,i)<=radius)-1;
        Indmat(Ord(1:K+1,i),i)=ones(K+1,1);
    end
    Indmat = (Indmat+Indmat')>=1;
    % fprintf([num2str(sum(Counts<0.5)) ' points have no neighbor. \n ']);
    T1 = quantile(Counts,delta)+1;
    Pen = @(x)((x+1)./T1).^-2.*(x<T1)+1;
    Pval = max(repmat(Pen(Counts)',[1 N]),repmat(Pen(Counts),[N 1]));
    Manipath = Pval.*Dis;
    Manipath(~Indmat) = inf;
    Manidis = Dis;
    Manidis(~Indmat) = inf;
    for k=1:N
        for i=1:N
            for j=1:N
                if Manipath(i,k)+Manipath(k,j)<Manipath(i,j)
                    Manipath(i,j) = Manipath(i,k)+Manipath(k,j);
                    Manidis(i,j) = Manidis(i,k)+Manidis(k,j);
                end
            end
        end
    end
    Manidis = (Manidis+Manidis')/2;
    % check unconnected points
    tempid = sum(Manidis==inf);
    if sum(tempid)==0
        % fprintf('All points are connected. \n');
        index = 1:N;
    elseif min(tempid)==sum(tempid==N-1)
        % fprintf([ num2str(min(tempid)) ' points are not connected to others. \n ']);
        index = find(tempid==min(tempid));
        Manidis = Manidis(index,index);
        Manipath = Manipath(index,index);
    else
        index = find(tempid==min(tempid));
        % fprintf(['There are more than one clusterings. The largest clustering with ' num2str(length(index))  ' points is embeded. \n '])
        Manidis = Manidis(index,index);
        Manipath = Manipath(index,index);
    end
end