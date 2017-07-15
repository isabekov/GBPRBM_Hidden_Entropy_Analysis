function [RMSE_Test, NEEoIHU_Test, H_vec] = MNIST_GBPRBM_Test_Multiple_Models
% This function tests performance of multiple Gaussian-Bipolar RBM models
% by measuring RMSE and Normalized Empirical Entropy of Individual Hidden Units.
% Parameters of the models are loaded from "*.mat" files located at 
% "$Model_Dir" directory. File names are taken according to the number of
% hidden units. Test data is loaded from "MNIST" directory.

% Output: figure with the following plots:
%          number of hidden units vs. RMSE
%          number of hidden units vs. NEEoIHU

% Number of hidden units for multiple models
H_vec = [5 10 15 32 49 75 100 128 144 256 400 512 784 900 1024 1500 2048];
H_vec_len = length(H_vec);

% Root-Mean Square Error
RMSE_Test = zeros(1, H_vec_len);

% Normalized Empirical Entropy of Individual Hidden Units
NEEoIHU_Test   = zeros(1, H_vec_len);

% Load test data
load(fullfile('MNIST', 'MNIST_Test_Medal_Normalized.mat'));
Model_Dir = 'MNIST_GBPRBM_Model';
for n=1:H_vec_len
    FileName = fullfile(Model_Dir, sprintf('MNIST_GBPRBM_Model,H=%i.mat', H_vec(n)));
    opt = load(FileName);        
    % Testing
    [RMSE_Test(n), NEEoIHU_Test(n)] = MNIST_GBPRBM_Test(H_vec(n), opt, testData, testLabels);
    % Close popped up figures
    close(gcf);
    close(gcf);
    close(gcf);
end

h_fig = figure('Name','MNIST: GBPRBM Number of Hidden Units vs. RMSE and Normalized Empirical Entropy of Individual Hidden Units', 'NumberTitle', 'Off');
subplot(121);
plot(H_vec, RMSE_Test, 'sq-','LineWidth',  2);
axis tight;
title('RMSE (test)');
xlabel('# of Hidden Units');

subplot(122);
plot(H_vec, NEEoIHU_Test, 'sq-','LineWidth', 2);
axis tight;
title('NEEoIHU (test)');
xlabel('# of Hidden Units');

FileName = fullfile(Model_Dir, 'Test_Multiple_Models_RMSE_NEEoIHU');
save([FileName '.mat'], 'RMSE_Test', 'NEEoIHU_Test', 'H_vec');
saveas(h_fig, [FileName '.png'], 'png');
hgsave(h_fig, [FileName '.fig']);