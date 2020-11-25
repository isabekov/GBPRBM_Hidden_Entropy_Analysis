function [opt, RMSE_Evolution] = MNIST_GBPRBM_Train(varargin)
% This function trains a Gaussian-Bipolar RBM model using default 
% parameters. If number of hidden units is not specified as the first 
% argument then default value is used. Training data is loaded from 
% the file or supplied as the second argument.
%                    Optional Inputs
% =================================================================
%         H:     positive integer            number of hidden units
% trainData:     N_Samples x V               training data
%
%                        Outputs
% =================================================================
%            opt:  a structure with model parameters
% RMSE_Evolution:  a vector with RMSE values for all mini-batch updates

if isempty(varargin)
    % Number of hidden units (default)
    opt.H = 1024;
elseif length(varargin)>=1
    opt.H = varargin{1};
end
if length(varargin)==2
    trainData = varargin{2};
else
    % Load MNIST training data
    load(fullfile('MNIST', 'MNIST_Train_Medal_Normalized.mat'));
end
% Set random number generator's seed value
rng(1);
% trainData = trainData(1:5000,:);
[N_Samples_Train, V] = size(trainData); 

% Visible variances
opt.z_v = -5*ones(V,1);

% Mini-batch size (samples)
mBatch_Size = 20;
% Number of epochs
N_Epochs = 2;
% Contrastive divergence order
opt.CD_order = 1;
% Learning rate
opt.nu = 1e-5;
% Regularization constant
opt.lambda = 0.1;
% Momentum
opt.mu = 0.5;
% Number of mini-batches
M = ceil(N_Samples_Train/mBatch_Size);
% Number of decimal digits for printing iterations
N_Iter_Digits = log10(M);
if rem(N_Iter_Digits,1)==0
    % Consider logarithm of 1,10,100,1000 etc.
    N_Iter_Digits = N_Iter_Digits + 1;
else
    N_Iter_Digits = ceil(N_Iter_Digits);
end
opt.N_Iter_Digits = N_Iter_Digits;

% Number of decimal digits for printing epochs
N_Epoch_Digits = log10(N_Epochs);
if rem(N_Epoch_Digits,1)==0
    % Consider logarithm of 1,10,100,1000 etc.
    N_Epoch_Digits = N_Epoch_Digits + 1;
else
    N_Epoch_Digits = ceil(N_Epoch_Digits);
end
opt.N_Epoch_Digits = N_Epoch_Digits;

RMSE_Evolution = zeros(1,N_Epochs*M);

opt.V = V;
opt.mBatch_Size = mBatch_Size;
opt.N_Samples_Train = N_Samples_Train;
opt.N_Epochs = N_Epochs;

% Initialize weights
opt.W = 0.01*randn(opt.V,opt.H);
% Initialize visible bias
opt.b_v = (mean(trainData,1))';
% Initialize hidden bias
opt.b_h = - opt.W'*diag(1./exp(opt.z_v))*opt.b_v;    
%% Training
for x = 1:N_Epochs
    opt.epoch = x;
    idx = randperm(N_Samples_Train);
    % Shuffle samples
    trainData = trainData(idx,:);   
    % Train the GRBM model using Contrastive Divergence
    opt = GBPRBM_Train(trainData, opt);   
    RMSE_Evolution((x-1)*M+1:x*M) = opt.RMSE_Evolution;
end
opt.DirSave = 'MNIST_GBPRBM_Model';
if exist(opt.DirSave,'dir') ~=7   
   mkdir(opt.DirSave); 
end
% A *.mat file to save the GBPRBM model
FileName = fullfile(opt.DirSave, sprintf('MNIST_GBPRBM_Model,H=%i.mat', opt.H));
save(FileName, '-struct', 'opt');
save(FileName, 'RMSE_Evolution', '-append');

h_fig = figure('Name', sprintf('MNIST GBPRBM, H=%d: Evolution of RMSE (train)', opt.H), ...
               'NumberTitle', 'Off'); 
plot(1:M*N_Epochs,RMSE_Evolution);
set(gca, 'YLim', [0 0.5]);
%axis tight;
title(sprintf('Evolution of RMSE (train), final RMSE=%0.3f', mean(RMSE_Evolution(end-mBatch_Size+1:end))));
xlabel('Iteration');
ylabel('RMSE');
FileName = fullfile(opt.DirSave, sprintf( 'Evolution_of_RMSE,H=%d', opt.H));
saveas(h_fig, [FileName '.png'], 'png');
hgsave(h_fig, [FileName '.fig']);