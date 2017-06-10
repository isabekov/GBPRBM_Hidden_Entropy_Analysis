function [RMSE_Test, NEEoHU] = MNIST_GBPRBM_Test(varargin)
% This function tests a Gaussian-Bipolar RBM model, which is passed in a
% form of a structure as the second argument or loaded from a *.mat file
% based on the number of hidden units specified as a the first argument. 
% Test data is loaded from the file or supplied as the third rgument.
%                    Optional Inputs
% =================================================================
%         H:     positive integer            number of hidden units
%       opt:     structure                   model parameters
%  testData:     N_Samples x V               test data
%
%                        Outputs
% =================================================================
%  RMSE_Test :  RMSE value of the test data
%     NEEoHU :  Normalized Empirical Entropy of Hidde Units

% A model to load
if isempty(varargin)
    % Number of hidden units (default)
    H = 49;
elseif length(varargin)>=1
    H = varargin{1};
end
if length(varargin)>=2
    opt = varargin{2};
else
    FileName = fullfile('MNIST_GBPRBM_Model', sprintf('MNIST_GBPRBM_Model,H=%i.mat', H));
    opt = load(FileName);
end
if length(varargin)==3
    testData = varargin{3};
else
    % Load MNIST test data
    load(fullfile('MNIST', 'MNIST_Test_Medal_Normalized.mat'));
end

%testData = testData(1:1000,:);
N_Samples_Test = size(testData,1); 
% Number of images to reconstruct
N_Img_Test = 100;
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

%% Plotting Normalized empirical entropy of hidden units, RMSE
h_fig = figure('Name', sprintf('MNIST: GBPRBM, H=%d: RMSE (test), Normalized Empirical Entropy of Hidden Units', opt.H), ...
               'NumberTitle', 'Off', 'Units', 'normalized', 'Position', [ 0.1432    0.4648    0.6792    0.4324]); 
subplot(131);
imagesc(H_Test);
xlabel('Test samples');
ylabel('Hidden unit index');
title(sprintf('Activations, RMSE (test) = %0.3f', RMSE_Test));

subplot(132);
[h_count, ~] = hist(H_Test',2);
p_h_m = h_count(2,:)/N_Samples_Test;
p_h_p = h_count(1,:)/N_Samples_Test;
bar(p_h_m);
set(gca, 'XLim', [0.5, opt.H+0.5]);
set(gca, 'YLim',[0 1]);
set(gca, 'View', [90 90]);
xlabel('Hidden unit index');
ylabel('p(h_j)');
title('Normalized histogram (Probability)');
Entropy = zeros(opt.H,1);
for j=1:opt.H
    if p_h_m(j)==1
       p_h_m(j)=1-eps; 
    elseif p_h_m(j)==0
       p_h_m(j)=1+eps; 
    end

    if p_h_p(j)==0
       p_h_p(j)=1+eps; 
    elseif p_h_p(j)==1
       p_h_p(j)=1-eps; 
    end
    Entropy(j) = -(p_h_p(j).*log2(p_h_p(j)) + p_h_m(j).*log2(p_h_m(j)));
end
% Normalizing empirical entropy of hidden units
NEEoHU = sum(Entropy)/opt.H;

subplot(133);
bar(Entropy);
set(gca, 'View', [90 90]);
set(gca, 'XLim', [0.5, opt.H+0.5]);
set(gca, 'YLim',[0 1]);
xlabel('Hidden unit index');
ylabel('Bits');
title(sprintf('Sum of Empirical Entropies of Hidden Units = %0.2f, Max = %d,\nNormalized Empirical Entropy of Hidden Units = %0.2f', ...
               sum(Entropy), opt.H, NEEoHU));    
drawnow;
FileName = fullfile(opt.DirSave, sprintf( 'NEEoHU,H=%02d', opt.H));
saveas(h_fig, [FileName '.png'], 'png');
hgsave(h_fig, [FileName '.fig']);

% Plot original and reconstructed images
MNIST_Plot_Pairs_of_Images(testData(1:N_Img_Test,:)', v_mean(:,1:N_Img_Test),...
      sprintf('MNIST: GBPRBM, H=%d, RMSE (Test)=%f, Orig. & Recon. Images', opt.H, RMSE_Test));
FileName = fullfile(opt.DirSave, sprintf( 'Orig_Recon_Images,H=%02d', opt.H));
saveas(gcf,[FileName '.png'], 'png');
hgsave(gcf,[FileName '.fig']);

N_Plots = min(opt.H,100);    
MNIST_Visualize_Weights(opt.W(:,1:N_Plots), sprintf('MNIST: GBPRBM (H=%d) Model Weights', opt.H));

FileName = fullfile(opt.DirSave, sprintf('Weights,H=%d', opt.H));
saveas(gcf,[FileName '.png'], 'png');
hgsave(gcf,[FileName '.fig']);