function Synthetic_Data_NV_Test_GBPRBM_RMSE_NEEoIHU_Plot
FileName = fullfile('Synth_Data', 'Test_Multiple_Models_RMSE_NEEoIHU.mat');
load(FileName);
[V_vec_len, H_vec_len]  = size(RMSE_Test);
V_vec_len = V_vec_len - 1;
V_vec = V_vec(2:end);
RMSE_Test = RMSE_Test(2:end,:);
NEEoIHU_Test = NEEoIHU_Test(2:end,:);
h_fig = figure('NumberTitle', 'Off', 'Position',  [500   264   844   324], ...
        'Name', 'GBPRBM: Number of Hidden Units vs. RMSE and Normalized Empirical Entropy of Individual Hidden Units');
subplot(121);
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

subplot(122);
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
axes('Position', [0.6 0.2 0.2 0.2]);
for i=1:V_vec_len
    plot(H_vec/V_vec(i), NEEoIHU_Test(i,:), LineSpec{i}, 'LineWidth', 2);
    hold on;
end
set(gca, 'XLim', [0 1.5]);
FileName = fullfile('Synth_Data', 'Test_Multiple_Models_RMSE_NEEoIHU');
saveas(h_fig, [FileName '.png'], 'png');
hgsave(h_fig, [FileName '.fig']);
