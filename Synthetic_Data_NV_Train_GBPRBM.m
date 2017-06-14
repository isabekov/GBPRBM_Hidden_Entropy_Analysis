function Synthetic_Data_NV_Train_GBPRBM
% This function trains multiple Gaussian-Bipolar RBM models
% Parameters of the models are stored in "*.mat" files in
% "Synth_Data_GBPRBM_Model_V=%d" directory, where %d specifies the number
% of visible units of the data. File names are generated according to the
% number of hidden units. Training data is loaded from "Synth_Data" directory.

% Number of visible units for different data
V_vec = 2.^(3:7);
V_vec_len = length(V_vec);
% Number of hidden units for multiple models
H_vec = [4 8 16 32 48 64 72 80 100 120 136 144 148 169 192 212 230 248 256];
H_vec_len = length(H_vec);
for k=1:V_vec_len
    for n=1:H_vec_len
        % Load Synth_Data training data
        load(fullfile('Synth_Data', sprintf('Synth_Data_Train_V=%d.mat', V_vec(k))));
        %trainData = trainData(1:5000,:);
        % Training
        Synthetic_Data_GBPRBM_Train(H_vec(n), trainData);
    end
end