function [c2n] = centroid_split(cn, D, stD )

[N D] =size(cn);
    
    ww = [0.01*stD.*randn(1,D)];
    for n=1:N
        c2n(2*n-1,:) = cn(n,:) + ww;
        c2n(2*n,:) = cn(n,:) - ww;
    end    
end