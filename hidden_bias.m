function [bh] = hidden_bias(c2n, w, data, stD )

[L, ~] =size(data);

cid = classification(c2n, data);
pdata = c2n(1,:);
ndata = c2n(2,:);
for k=1:L
    if mod(cid(k),2) == 1
        pdata = [pdata; data(k,:)];
    else
        ndata = [ndata; data(k,:)];
    end
end
hp=(2*pdata*w'/stD);
hn=(2*ndata*w'/stD);
bh = -0.25*(min(hp)+max(hn));

end