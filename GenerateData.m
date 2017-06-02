function data = GenerateData(V,NC,stD,N_Samples,seed)
rng(seed);
c=zeros(NC,V);
for k=1:NC
    if k==1
        c(k,:) = rand(1,V);
    else
        d=0;
        while d < 1*stD
            x=rand(1,V);
            d=2;
            for j=1:k-1  
                d=min(sqrt((c(j,:)-x)*(c(j,:)-x)'),d);
            end
        end
        c(k,:) = x;
    end
end

k=1;

data = repmat(c(k,:),N_Samples,1) + stD.*randn(N_Samples,V);
for k=2:NC
    data = [data; repmat(c(k,:),N_Samples,1) + stD.*randn(N_Samples,V)];
end