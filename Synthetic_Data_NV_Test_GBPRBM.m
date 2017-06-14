function [RMSE_Test, NEEoIHU_Test, V_vec, H_vec] = Synthetic_Data_NV_Test_GBPRBM
% This function tests performance of multiple Gaussian-Bipolar RBM models
% by measuring RMSE and Normalized Empirical Entropy of Individual Hidden
% Units. Parameters of the models are loaded from "*.mat" files located at
% "Synth_Data_GBPRBM_Model_V=%d" directory, where %d is the number of
% visible units which characterizes the data. File names of the models are
% taken according to the number of hidden units. Test data is loaded from
% "Synth_Data" directory.

% Output: figure with the following plots:
%          number of hidden units vs. RMSE
%          number of hidden units vs. NEEoIHU

% Number of visible units for different data
V_vec = 2.^(3:7)';
V_vec_len = length(V_vec);
% Number of hidden units for multiple models
H_vec = [4 8 16 32 48 64 72 80 100 120 136 144  148 169 192 212 230 248 256];
H_vec_len = length(H_vec);
RMSE_Test = zeros(V_vec_len, H_vec_len);
NEEoIHU_Test = zeros(V_vec_len, H_vec_len);
for k=1:V_vec_len
    for n=1:H_vec_len
        Model_Dir = sprintf('Synth_Data_GBPRBM_Model_V=%d', V_vec(k));
        FileName = fullfile(Model_Dir, sprintf('Synth_Data_GBPRBM_Model,H=%i.mat', H_vec(n)));
        opt = load(FileName);
        % Load Synth_Data test data
        load(fullfile('Synth_Data', sprintf('Synth_Data_Test_V=%d.mat', V_vec(k))));
        % testData = testData(1:1000,:);
        % Testing
        [RMSE_Test(k,n), NEEoIHU_Test(k,n)] = Synthetic_Data_GBPRBM_Test(H_vec(n), opt, testData);
        close(gcf);
    end
end

h_fig = figure('NumberTitle', 'Off', 'Name', 'GBPRBM: Number of Hidden Units vs. RMSE and Normalized Empirical Entropy of Individual Hidden Units');
subplot(221);
LineSpec={'k-.s', '-*', 'm--o', 'c-.d' , '-<', '--<'};
for i=1:V_vec_len
    plot(H_vec, RMSE_Test(i,:), LineSpec{i}, 'LineWidth', 2);
    hold on;
end
YLim = get(gca, 'YLim');
axis tight;
set(gca, 'YLim', YLim);
legend([repmat('V = ', V_vec_len, 1), num2str(V_vec)]);
title('RMSE (test)');
xlabel('# of Hidden Units');

subplot(222);
for i=1:V_vec_len
    plot(H_vec, NEEoIHU_Test(i,:), LineSpec{i}, 'LineWidth', 2);
    hold on;
end
YLim = get(gca, 'YLim');
axis tight;
set(gca, 'YLim', YLim);
legend([repmat('V = ', V_vec_len, 1), num2str(V_vec)]);
title('NEEoIHU (test)');
xlabel('# of Hidden Units');

subplot(223);
LineSpec={'k-.s', '-*', 'm--o', 'c-.d' , '-<', '--<'};
for i=1:V_vec_len
    plot(H_vec/V_vec(i), RMSE_Test(i,:), LineSpec{i}, 'LineWidth', 2);
    hold on;
end
YLim = get(gca, 'YLim');
axis tight;
set(gca, 'YLim', YLim);
legend([repmat('V = ', V_vec_len, 1), num2str(V_vec)]);
title('RMSE (test)');
xlabel('H/V');

subplot(224);
for i=1:V_vec_len
    plot(H_vec/V_vec(i), NEEoIHU_Test(i,:), LineSpec{i}, 'LineWidth', 2);
    hold on;
end
YLim = get(gca, 'YLim');
axis tight;
set(gca, 'YLim', YLim);
legend([repmat('V = ', V_vec_len, 1), num2str(V_vec)]);
title('NEEoIHU (test)');
xlabel('H/V');

FileName = fullfile('Synth_Data', 'Test_Multiple_Models_RMSE_NEEoIHU');
save([FileName '.mat'], 'RMSE_Test', 'NEEoIHU_Test', 'H_vec', 'V_vec');
saveas(h_fig, [FileName '.png'], 'png');
hgsave(h_fig, [FileName '.fig']);

% surface(repmat(V_vec, 1, H_vec_len), repmat(H_vec, V_vec_len, 1), RMSE_Test);
% set(gca, 'View', [129 22])
% title('RMSE (test)');
% xlabel('# of Visible Units');
% ylabel('# of Hidden Units');
%
% subplot(122);
% surface(repmat(V_vec, 1, H_vec_len), repmat(H_vec, V_vec_len, 1), NEEoIHU_Test);
% set(gca, 'View', [129 22])
% title('NEEoIHU (test)');
% xlabel('# of Visible Units');
% ylabel('# of Hidden Units');