function [cid] = classification( c, data )
% N = 2^j
% D: number of visible units
[N D] =size(c);
[K D1] =size(data);
cid=zeros(K,1);

if D ~= D1
    error('Dimension mismatch!');
else
    dist = zeros(N,1);
    for k=1:K
        for n=1:N
            dist(n) = norm(c(n,:)-data(k,:));
        end
        [~, cid(k)] = min(dist);
    end
end
end