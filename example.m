%Take forest as an example.
load('data.mat');
t=3;%time scale=3 months
p_for=Precipitation.*repmat(land_forest,[1,1,240]);
e_for=EVI.*repmat(land_forest,[1,1,240]);
l_for=LST.*repmat(land_forest,[1,1,240]);
[m,n,z] = size(p_for);
DI_for=zeros(m,n,(z-t+1));
for i=1:m
    for j=1:n
        if any(isnan(p_for(i,j,:)))
            DI_for(i,j,:)=NaN;
        else
            p=p_for(i,j,:);
            P = permute( p,[3 2 1]);
            e=e_for(i,j,:);
            E = permute( e,[3 2 1]);
            l=l_for(i,j,:);
            L = permute( l,[3 2 1]);
            data=[P E L];
            D=index_for(data,t);
            DI_for(i,j,:) = permute(D,[1 3 2]);
        end
    end
end