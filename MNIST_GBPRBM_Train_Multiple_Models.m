function opt = MNIST_GBPRBM_Train_Multiple_Models
% Load MNIST training data
load(fullfile('MNIST', 'MNIST_Train_Medal_Normalized.mat'));
%trainData = trainData(1:5000,:);
H_vec = [10 49 100 144 256 512 784 1024 2048];
H_vec_len = length(H_vec);
RMSE_Test_vec = zeros(1,H_vec_len);
NEEoHU_vec = zeros(1,H_vec_len);
% Load MNIST test data
load(fullfile('MNIST', 'MNIST_Test_Medal_Normalized.mat'));
for n=1:H_vec_len
    %% Training
    [opt ] = MNIST_GBPRBM_Train(H_vec(n), trainData);
    %% Testing
    [RMSE_Test_vec(n), NEEoHU_vec(n)] = MNIST_GBPRBM_Test(H_vec(n), opt, testData);
end
h_fig = figure('Name','Number of hidden units vs. RMSE and Normalized Empirical Entropy of Hidden Units', 'NumberTitle', 'Off');
subplot(211); plot(H_vec, RMSE_Test_vec);
title('RMSE (test)');

subplot(212); plot(H_vec, NEEoHU_vec);
xlabel('# of hidden units');
title('Normalized Empirical Entropy of Hidden Units');
FileName = fullfile(opt.DirSave, 'NEEoHU');
saveas(h_fig, [FileName '.png'], 'png');
hgsave(h_fig, [FileName '.fig']);
