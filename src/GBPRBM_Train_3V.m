function [opt]  = GBPRBM_Train_3V(data, opt)
% Inputs:
% Name            Size             Explanation
% data.features:  (N_Samples x V)  "N_Samples" instances of feature vector of size (1 x V)
% data.feat_mean: (1 x V)          Mean of the data features 
% data.feat_std:  (1 x V)          Standard deviation of features
% opt.mBatch_Size: (integer)      Mini-batch size
% opt.H:         (positive integer)        Number of hidden units
% opt.k:         (positive integer)        Contrastive divergence order (number of Gibbs steps)
% opt.nu:        (positive real)  Learning rate
% opt.lambda:    (positive real)  Regularization constant
% opt.mu:        (positive real)  Momentum
% opt.epoch:     (positive integer)  Momentum

% Output
% RMSE_Evolution
% Centroid_Evolution
% b_v_Evolution

[N_Samples, V] = size(data.features);

% Number of mini-batches
M = ceil(N_Samples/opt.mBatch_Size);
opt.M = M;
% Starting and ending indices of the batch
BI = cell(1,M);
% Partition the batch into M mini-batches
for m=1:M-1
   BI{m}.idx1 = opt.mBatch_Size*(m-1) + 1;
   BI{m}.idx2 = opt.mBatch_Size*m;
end
BI{end}.idx1 = opt.mBatch_Size*(M-1) + 1;
BI{end}.idx2 = N_Samples;

H = opt.H;
opt.Centroid_Evolution = cell(H,1);
for j=1:H
    opt.Centroid_Evolution{j} = zeros(V, 2^j, M);
end
opt.b_v_Evolution = zeros(V,M);


H = opt.H;

optim_HE = opt.optim_HE;

% Initialize delta updates
D_W   = zeros(V,H);
D_b_v = zeros(V,1);
D_b_h = zeros(H,1);
D_z_v = zeros(V,1);

RMSE_Evolution = zeros(1,M);

opt.m = 0;
if opt.Enable_Visual_Debugging==true
    %ContrastiveDivergenceVisualDebugging(opt);
end

% For every mini-batch "m" do:
for m=1:M
    opt.m = m;
    % Clear derivatives of the log-likelihood function:
    d_W = zeros(V,H);
    d_b_v = zeros(V,1);
    d_b_h = zeros(H,1);
    d_z_v = zeros(V,1);       
    % Assign mini-batch
    mBatch = data.features(BI{m}.idx1:BI{m}.idx2,:)';
    RMSE = 0;
    for v=mBatch
       opt.v_0 = v;
       [out, v_cd ]= GBPRBM_Contrastive_Divergence(opt);       
       % Contrastive divergence updates
       d_W = d_W + out.d_W;
       d_b_v = d_b_v + out.d_b_v;
       d_b_h = d_b_h + out.d_b_h;
       d_z_v = d_z_v + out.d_z_v;
       % Root Mean Square Error  
       RMSE = RMSE + sqrt(1/V*sum((v - v_cd(:,opt.CD_order+1)).^2));
    end
    % Normalize by mini-batch size
    sz = size(mBatch,2);
    d_W = d_W/sz;
    d_b_v = d_b_v/sz;
    d_b_h = d_b_h/sz;
    d_z_v = d_z_v/sz;
    if optim_HE == true
        [p_h TT] = Calculate_Probability_of_Hidden_Units(opt);      
        Entropy = Calculate_Entropy(p_h);
        
        [d_W_HE, d_b_v_HE, d_b_h_HE]= Derivative_of_HEntropy(opt, p_h, TT);
        
        d_W = d_W + 1/Entropy*d_W_HE;
        d_b_v = d_b_v + 1/Entropy*d_b_v_HE;        
        d_b_h = d_b_h + 1/Entropy*d_b_h_HE;
    end

    RMSE_Evolution(m) = RMSE/sz;
    
    % Update all parameters
    D_W = opt.mu*(opt.nu*d_W - opt.lambda*opt.W) + (1-opt.mu)*D_W;
    opt.W = opt.W + D_W;
    
    D_b_v = opt.mu*(opt.nu*d_b_v - opt.lambda*opt.b_v) + (1-opt.mu)*D_b_v;
    opt.b_v = opt.b_v + D_b_v;
    
    D_b_h = opt.mu*(opt.nu*d_b_h - opt.lambda*opt.b_h) + (1-opt.mu)*D_b_h;
    opt.b_h = opt.b_h + D_b_h;
%     opt.b_h = -opt.W'*opt.Sinv*opt.b_v;
    
    D_z_v = opt.mu*(opt.nu*d_z_v - opt.lambda*opt.z_v) + (1-opt.mu)*D_z_v;
    %opt.z_v = opt.z_v + D_z_v;
              
    if (rem(m,20)==0) || (m==1) || (m==M)
        str = ['Epoch %0' num2str(opt.N_Epoch_Digits) 'd/' num2str(opt.N_Epochs) ...
               ', iteration %0' num2str(opt.N_Iter_Digits) 'd/' num2str(opt.M) ', RMSE=%f\n'];
        fprintf(str, opt.epoch, opt.m, RMSE_Evolution(m));
    end
    
    if opt.Enable_Visual_Debugging==true
        ContrastiveDivergenceVisualDebugging(opt);
    end
    [opt] = GBPRBM_Visible_Units_Span(opt);
    for j=1:H
        opt.Centroid_Evolution{j}(:,:,m) = opt.Centroid{j};
    end
    opt.b_v_Evolution(:,m) = opt.b_v;
    %opt.nu = opt.nu*0.999;
end
opt.RMSE_Evolution = RMSE_Evolution;
disp('======================================================');



function ContrastiveDivergenceVisualDebugging(opt)
set(0, 'CurrentFigure', opt.h_fig);
cla(opt.h_axes);
%opt = GBPRBM_Visible_Units_Span(opt);
%GBPRBM_Plot_PDF_of_Visible_Units_2V(opt);
GBPRBM_Visible_Units_Plot_Geometry(opt);
opt.v(1).AxesLimMin = opt.XLimLow;
opt.v(1).AxesLimMax = opt.XLimHigh;
opt.v(2).AxesLimMin = opt.YLimLow;
opt.v(2).AxesLimMax = opt.YLimHigh;
set(opt.h_axes, 'XLim', [0 1]);
set(opt.h_axes, 'YLim', [0 1]);
set(opt.h_axes, 'ZLim', [0 1]);
view( -37.5000,30);
%GBPRBM_Plot_Decision_Region_for_Visible_Units(opt);
str = ['p(v_1,v_2) at epoch #%0' num2str(opt.N_Epoch_Digits) 'd, iteration #%0' num2str(opt.N_Iter_Digits) 'd'];
title( sprintf(str, opt.epoch, opt.m));
colorbar('off');
F = getframe(gcf);
[im,map] = frame2im(F);    % Return associated image data 
if isempty(map)            % Truecolor system
  rgb = im;
else                       % Indexed system
  rgb = ind2rgb(im,map);   % Convert image data
end
str = ['CD_%0' num2str(opt.N_Epoch_Digits)  'de_%0' num2str(opt.N_Iter_Digits)  'di.png'];
FileName = fullfile(opt.DirSave, sprintf(str,opt.epoch,opt.m));
imwrite(rgb,FileName,'png');  



function [d_W, d_b_v, d_b_h] = Derivative_of_HEntropy(opt, p_h, TT)
Sinv = opt.Sinv;
V = opt.V;
H = opt.H;
W = opt.W;
b_v = opt.b_v;
E_W = zeros(V,H);
E_b_v = zeros(V,1);
E_b_h = zeros(H,1);

for i = 1:2^H       
    E_W   = E_W + p_h(i)*Sinv*(b_v*TT(:,i)' + (W*TT(:,i))*TT(:,i)');
    E_b_v = E_b_v + p_h(i)*Sinv*W*TT(:,i);
    E_b_h = E_b_h + p_h(i)*TT(:,i);
end

d_W   = zeros(V,H);
d_b_v = zeros(V,1);
d_b_h = zeros(H,1);
for i = 1:2^H       
    d_W   = d_W   + p_h(i)*(Sinv*(b_v*TT(:,i)' + (W*TT(:,i))*TT(:,i)') - E_W)*(log(p_h(i)) + 1);
    d_b_v = d_b_v + p_h(i)*(Sinv*W*TT(:,i) - E_b_v)*(log(p_h(i)) + 1);
    d_b_h = d_b_h + p_h(i)*(TT(:,i) - E_b_h)*(log(p_h(i)) + 1);
end
d_W   = -1/log(2)*d_W;
d_b_v = -1/log(2)*d_b_v;
d_b_h = -1/log(2)*d_b_h;
