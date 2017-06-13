function opt = GBPRBM_Train(data, opt)
%                       Inputs
% =================================================================
% Name                   Size             Explanation
% data:            (N_Samples x V)     "N_Samples" instances of feature vector of size (1 x V)
% opt.mBatch_Size: (positive integer)  Mini-batch size
% opt.H:           (positive integer)  Number of hidden units
% opt.k:           (positive integer)  Contrastive divergence order (number of Gibbs steps)
% opt.nu:          (positive real)     Learning rate
% opt.lambda:      (positive real)     Regularization constant
% opt.mu:          (positive real < 1) Momentum
% opt.epoch:       (positive integer)  Epoch number to display

%                       Outputs
% =================================================================
% saved *.mat file with model parameters ("opt" structure)
% opt.RMSE_Evolution

[N_Samples, V] = size(data);

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

% Initialize delta updates
D_W   = zeros(V,H);
D_b_v = zeros(V,1);
D_b_h = zeros(H,1);
%D_z_v = zeros(V,1);

RMSE_Evolution = zeros(1,M);
disp('===========================================');
fprintf('             GBPRBM, V = %d, H = %d           \n', V, H);
% For every mini-batch "m" do:
for m=1:M
    opt.m = m;
    % Clear derivatives of the log-likelihood function:
    d_W = zeros(V,H);
    d_b_v = zeros(V,1);
    d_b_h = zeros(H,1);
    %d_z_v = zeros(V,1);       
    % Assign mini-batch
    mBatch = data(BI{m}.idx1:BI{m}.idx2,:)';
    RMSE = 0;
    for v=mBatch
       opt.v_0 = v;
       [out, v_cd ]= GBPRBM_Contrastive_Divergence(opt);       
       % Contrastive divergence updates
       d_W = d_W + out.d_W;
       d_b_v = d_b_v + out.d_b_v;
       d_b_h = d_b_h + out.d_b_h;
       %d_z_v = d_z_v + out.d_z_v;
       % Root Mean Square Error  
       RMSE = RMSE + sqrt(1/V*sum((v - v_cd(:,opt.CD_order+1)).^2));
    end
    % Normalize by mini-batch size
    mBatch_Size = size(mBatch,2);
    d_W = d_W/mBatch_Size;
    d_b_v = d_b_v/mBatch_Size;
    d_b_h = d_b_h/mBatch_Size;
    %d_z_v = d_z_v/mBatch_Size;
    
    RMSE_Evolution(m) = RMSE/mBatch_Size;
    
    % Update all parameters
    D_W = (1-opt.mu)*opt.nu*(d_W - opt.lambda*opt.W) + opt.mu*D_W;
    opt.W = opt.W + D_W;
    
    D_b_v = (1-opt.mu)*opt.nu*(d_b_v - opt.lambda*opt.b_v) + opt.mu*D_b_v;
    opt.b_v = opt.b_v + D_b_v;
    
    D_b_h = (1-opt.mu)*opt.nu*(d_b_h - opt.lambda*opt.b_h) + opt.mu*D_b_h;
    opt.b_h = opt.b_h + D_b_h;
    
    %D_z_v = (1-opt.mu)*(opt.nu*d_z_v - opt.lambda*opt.z_v) + opt.mu*D_z_v;
    %opt.z_v = opt.z_v + D_z_v;
              
    if (rem(m,20)==0) || (m==1) || (m==M)
        str = ['Epoch %0' num2str(opt.N_Epoch_Digits) 'd/' num2str(opt.N_Epochs) ...
               ', iteration %0' num2str(opt.N_Iter_Digits) 'd/' num2str(opt.M) ', RMSE=%f\n'];
        fprintf(str, opt.epoch, opt.m, RMSE_Evolution(m));
    end            
end
opt.RMSE_Evolution = RMSE_Evolution;