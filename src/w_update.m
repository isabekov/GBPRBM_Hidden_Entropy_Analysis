function [c, w] = w_update(c2n, cn, nk )

[M, D] =size(c2n);
[N, D] =size(cn);

if M ~= 2*N
    error('Dimension mismatch!');
else
    dist = zeros(M,1);
    for m=1:M
        ww = c2n(m,:)- cn(floor((m+1)/2),:);
        for n=1:N
            c_o_diff = c2n(2*n-1,:) - cn(n,:);
            c_e_diff = c2n(2*n,:)   - cn(n,:);
            d1 = norm(c_o_diff - ww ) + norm(c_e_diff + ww );
            d2 = norm(c_o_diff + ww ) + norm(c_e_diff - ww );
            dist(m) = dist(m) + min(d1, d2);
        end
    end
    nknorm = (nk+1) ./ sum(nk);
    rnk = 1./nk;
    rnknorm = rnk ./ sum(rnk);
    distnorm = dist ./ sum(dist);
    alpha=0.7;
    dist = alpha*distnorm + (1-alpha)*rnknorm;
    
    [~, m] = min(dist);
    w = c2n(m,:)- cn(floor((m+1)/2),:);
%     plot(dist);
%     pause(3);
    
    for n=1:N
        c(2*n-1,:) = cn(n,:) + w;
        c(2*n,:) = cn(n,:) - w;
    end    
end    
end