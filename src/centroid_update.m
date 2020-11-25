function [c, nk] = centroid_update(cid, data, N )
% N: number of centroids
% K: number of data samples
% D: number of visible units
% cid can have values from 1 to 2^H

[K M] =size(cid);
[K D] =size(data);

c = zeros(N,D);
nk = zeros(N,1);

for k=1:K
    c(cid(k),:) = c(cid(k),:) + data(k,:);
    nk(cid(k)) = nk(cid(k))+1;
end

for n=1:N
    c(n,:) = c(n,:)./nk(n);
end
    
end