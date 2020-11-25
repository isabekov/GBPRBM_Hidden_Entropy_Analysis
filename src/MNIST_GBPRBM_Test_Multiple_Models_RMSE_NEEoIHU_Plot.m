function MNIST_GBPRBM_Test_Multiple_Models_RMSE_NEEoIHU_Plot
Model_Dir = 'MNIST_GBPRBM_Model';
FileName = fullfile(Model_Dir, 'Test_Multiple_Models_RMSE_NEEoIHU.mat');
load(FileName);
% Number of visible units
V = 784;

h_fig = figure('Name','MNIST: GBPRBM H/V vs. RMSE and NEEoIHU', 'NumberTitle', 'Off');
subplot(121);
plot(H_vec/V, RMSE_Test, 'sq-','LineWidth',  2);
axis tight;
title('RMSE (test)');
xlabel('H/V');

subplot(122);
plot(H_vec/V, NEEoIHU_Test, 'sq-','LineWidth', 2);
axis tight;
title('NEEoIHU (test)');
xlabel('H/V');

FileName = fullfile(Model_Dir, 'Test_Multiple_Models_RMSE_NEEoIHU');
save([FileName '.mat'], 'RMSE_Test', 'NEEoIHU_Test', 'H_vec');
saveas(h_fig, [FileName '.png'], 'png');
hgsave(h_fig, [FileName '.fig']);