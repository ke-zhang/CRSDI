%Take forest land as an example.
load('data.mat')
t=3;%time scale=3 month
[m,n,z] = size(Precipitation);
for i = 1:z
    p_for(:,:,i)=Precipitation(:,:,i).*land_forest;
    e_for(:,:,i)=EVI(:,:,i).*land_forest;
    l_for(:,:,i)=LST(:,:,i).*land_forest;
end
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