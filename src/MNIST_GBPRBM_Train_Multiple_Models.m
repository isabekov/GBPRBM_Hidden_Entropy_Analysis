function MNIST_GBPRBM_Train_Multiple_Models
% This function trains multiple Gaussian-Bipolar RBM models
% Parameters of the models are stored in "*.mat" files in "MNIST_GBPRBM_Model"
% directory. File names are generated according to the number of
% hidden units. Training data is loaded from "MNIST" directory.

% Load MNIST training data
load(fullfile('MNIST', 'MNIST_Train_Medal_Normalized.mat'));
%trainData = trainData(1:5000,:);

% Number of hidden units for multiple models
H_vec = [5 10 15 32 49 75 100 128 144 256 400 512 784 900 1024 1500 2048];
H_vec_len = length(H_vec);

for n=1:H_vec_len
    %% Training
    MNIST_GBPRBM_Train(H_vec(n), trainData);    
end