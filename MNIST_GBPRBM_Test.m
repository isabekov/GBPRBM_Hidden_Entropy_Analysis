function [RMSE_Test, NEEoHU] = MNIST_GBPRBM_Test(varargin)
% This function tests a Gaussian-Bipolar RBM model, which is passed in a
% form of a structure as the second argument or loaded from a *.mat file
% based on the number of hidden units specified as a the first argument. 
% Test data is loaded from the file or supplied as the third argument.
%                    Optional Inputs
% =================================================================
%         H:     positive integer            number of hidden units
%       opt:     structure                   model parameters
%  testData:     N_Samples x V               test data (V = 784)
%  testLabels    N_Samples x N_Classes       test labels (N_Classes = 10)
%
%                        Outputs
% =================================================================
%  RMSE_Test :  RMSE value of the test data
%    NEEoIHU :  Normalized Empirical Entropy of Individual Hidden Units

% A model to load
if isempty(varargin)
    % Number of hidden units (default)
    H = 1024;
elseif length(varargin)>=1
    H = varargin{1};
end
if length(varargin)>=2
    opt = varargin{2};
else
    FileName = fullfile('MNIST_GBPRBM_Model', sprintf('MNIST_GBPRBM_Model,H=%i.mat', H));
    opt = load(FileName);
end
if length(varargin)==4
    testData = varargin{3};
    testLabels = varargin{4};
else
    % Load MNIST test data
    load(fullfile('MNIST', 'MNIST_Test_Medal_Normalized.mat'));
end

%testData = testData(1:1000,:);
N_Samples_Test = size(testData,1); 
%% Testing
H_Test = zeros(opt.H, N_Samples_Test);
% Preallocate memory for the reconstructed images
v_mean = zeros(opt.V, N_Samples_Test);
RMSE = 0;
for i=1:N_Samples_Test
    [H_Test(:,i), ~, ~] = GBPRBM_Gibbs_Sampling_of_Hidden_Units(opt, testData(i,:)');
    [v_mean(:,i), ~]    = GBPRBM_Gibbs_Sampling_of_Visible_Units(opt, H_Test(:,i));  
    RMSE = RMSE + sqrt(1/opt.V*sum((testData(i,:)' - v_mean(:,i)).^2));
end
RMSE_Test = RMSE/N_Samples_Test;

% Plot Normalized Empirical Entropy of Individual Hidden Units, RMSE
[h_fig, NEEoHU] = GBPRBM_NEEoIHU(H_Test, sprintf('MNIST: GBPRBM, H=%d: Normalized Empirical Entropy of Individual Hidden Units', opt.H), RMSE_Test);
FileName = fullfile(opt.DirSave, sprintf( 'NEEoHU,H=%02d', opt.H));
saveas(h_fig, [FileName '.png'], 'png');
hgsave(h_fig, [FileName '.fig']);

% Number of images per digit to display
N_IPD = 10;
% Number of unique digits
N_UD = 10;
idx = zeros(N_UD, N_IPD);
for j=1:N_UD
    cnt = 1;
    for k=1:N_Samples_Test
       Digit = find(testLabels(k,:));
       if Digit == j
          idx(j,cnt) = k;
          cnt = cnt + 1;
       end
       if cnt == N_IPD+1
          break
       end
    end
end
idx = reshape(idx, N_IPD*N_UD,1);
% Plot original and reconstructed images
h_fig = MNIST_Plot_Pairs_of_Images(testData(idx,:)', v_mean(:,idx),...
      sprintf('MNIST: GBPRBM, H=%d, RMSE (Test)=%f, Orig. & Recon. Images', opt.H, RMSE_Test));
FileName = fullfile(opt.DirSave, sprintf( 'Orig_Recon_Images,H=%02d', opt.H));
saveas(h_fig,[FileName '.png'], 'png');
hgsave(h_fig,[FileName '.fig']);

% Plot weights of the GBPRBM model
N_Plots = min(opt.H,100);    
h_fig = MNIST_Visualize_Weights(opt.W(:,1:N_Plots), sprintf('MNIST: GBPRBM (H=%d) Model Weights', opt.H));
FileName = fullfile(opt.DirSave, sprintf('Weights,H=%d', opt.H));
saveas(h_fig,[FileName '.png'], 'png');
hgsave(h_fig,[FileName '.fig']);