clear;
rng(10);
NC=8;     % Number of centroids
stD=0.01; % STD of Gaussians
N = 2000; % Number of random samples per centroid
D=3;      % Observation dimension
seed = 1;
data = GenerateData(D,NC,stD*4,N,seed);
NH=9;     % Number of hidden nodes
c1 = mean(data);
c0 = c1;
b_v = c1';
clear c;

for j=1:NH
    jj=pow2(j);
    %c(j,1:jj,1:D) = centroid_split(c0, D, stD );
    c(j,1:jj,1:D) = centroid_split(c0, D, stD );
    clear ct;
    
    for k=1:5
        for kk=1:1
            ct(1:jj,1:D) = c(j,:,:);
            cid = classification(ct, data);
            [c(j,1:jj,1:D), nk] = centroid_update(cid, data, jj);
        end
        ct(1:jj,1:D) = c(j,:,:);
        [c(j,1:jj,1:D), w(j,:)] = w_update(ct, c0, nk);
    end
    ct(1:jj,1:D) = c(j,:,:);
    bh(j) = hidden_bias(ct, w(j,:), data, stD );
    
    c0(1:jj,1:D) = c(j,:,:);
    %c0 = c(j,:,:);
    fprintf('Hidden unit %i\n', j);
end


for j=1:NH
    h(j,1:N*NC)=sign(2*data*w(j,:)'/stD + 2*bh(j));
end

figure('NumberTitle', 'Off', 'Name', 'GBPRBM LBG-like pretraining');
scatter3(data(:,1),data(:,2),data(:,3),'.y');
hold on;
for k=1:10:N*NC
    if h(1,k)>0
        scatter3(data(k,1),data(k,2),data(k,3),'.g');
    else
        scatter3(data(k,1),data(k,2),data(k,3),'.b');
    end
end
% for k=2:100:N*NC
%     if h(2,k)>0
%         plot(data(k,1),data(k,2),'*g');
%     else
%         plot(data(k,1),data(k,2),'+b');
%     end
% end
       
plot3([0;c1(1)],[0;c1(2)],[0;c1(3)],'b', 'LineWidth',2);
pp(1,:) = c1;
for j=1:NH
    jj=pow2(j);
    for k=1:2:jj
        kk = ceil(k/2);
        p(k,:) = pp(kk,:) + (-1)^(k-1)*w(j,:);
        plot3([pp(kk,1);p(k,1)], [pp(kk,2);p(k,2)], [pp(kk,3);p(k,3)], 'b', 'LineWidth',2);
        p(k+1,:) = pp(kk,:) + (-1)^(k)*w(j,:);
        plot3([pp(kk,1);p(k+1,1)],[pp(kk,2);p(k+1,2)],[pp(kk,3);p(k+1,3)],'r', 'LineWidth',2);
    end
    pp = p;
end
set(gca, 'Projection', 'Perspective')
clear codev;
clear hd;

for k=1:pow2(NH)
    y(1:D)=c(NH,k,:);
    fprintf('c%d: [',k);
    fprintf('%4.2f ',y);
    fprintf('] with h: ');
    str=sprintf('h: ');
    for j=1:NH
        hd(j)=sign(2*y*w(j,:)'/stD + 2*bh(j));
        fprintf('%3d',hd(j));
        ss=sprintf('%d',hd(j));
        str=strcat(str,ss);
    end
    fprintf('\n');
%    str = sprintf('%d%d%d',hd(1),hd(2),hd(3));
%    str = string(hd);
 %   text(y(1)+0.01,y(2)+0.01,y(3)+0.01,str);
    codev(k) = (hd+1)./2 * pow2([NH-1:-1:0])';
end

fprintf('Out of %d codevectors %d of them are in use.\n',pow2(NH), max(size(unique(codev))));

hold off;

sigma_v = stD*ones(D,1);
b_h = bh';
W = w';
save('GeometryLBG.mat', 'W','b_h', 'b_v', 'sigma_v');
