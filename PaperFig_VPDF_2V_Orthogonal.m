function PaperFig_VPDF_2V_Orthogonal
opt.V = 2;
opt.H = 2;
% Orthogonal weights
opt.W = [9   2; ...
         2  -9];
% Visible bias
opt.b_v = [11 12]';
% Standard deviations for visible units
opt.sigma_v = [2 2]';
% Covariance matrix
opt.S = diag(opt.sigma_v.^2);
% Inverse of the "covariance matrix"
opt.Sinv = diag(1./(opt.sigma_v.^2));
% Hidden bias
opt.b_h =  (-opt.b_v'*opt.Sinv*opt.W)';
% Span ratio for all axes
opt.VRatio = [0.5 0.5];
% Number of samples used to quantize X-axis
opt.Nv1 = 150 ;
% Number of samples used to quantize Y-axis
opt.Nv2 = 151;
% Do not plot decision boundaries for p(h_j|v)
opt.plot_decision_regions = false;
opt.put_labels = true;
opt.FontSize = 10;
opt.Color = {[0, 0.6, 0.6],'red', [0, 0.7, 0], [0.8706, 0.4902, 0], 'green'};

figure('Name', 'Orthogonal Weights', ...
   'Units', 'normalized', 'Position', [ 0.1018    0.3919    0.6691    0.4753], ...
   'NumberTitle', 'Off');
opt.h_axes = subplot(121);
opt.h_axes_pdf = opt.h_axes;
opt = GBPRBM_Plot_PDF_of_Visible_Units_2V(opt);
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(a) ' str]);
%set(opt.h_cbar, 'Location', 'WestOutside');


opt.plot_type = 'imagesc';
opt.HRange = [10 10];
opt.N_Samples = [250 251];
opt.h_axes = subplot(122);
opt.h_axes_ent = opt.h_axes;
opt = GBPRBM_Plot_HEntropy_vs_HBias_2H_Experiment(opt);
str = get(opt.h_title, 'String');
str = strrep(str, ' (experiment)', '');
set(opt.h_title, 'String', ['(b) ' str]);
set(opt.h_cbar, 'YTick', [0.02 0.5 1 1.5 1.98]);
set(opt.h_cbar, 'YTickLabel', [0 0.5 1 1.5 2]);

opt.b_h_md =  (-opt.b_v'*opt.Sinv*opt.W)';
[x1 y1] = dsxy2figxy(opt.h_axes_ent, opt.b_h_md(1), opt.b_h_md(2));

Pos = get(opt.h_axes_pdf, 'Position');
x2 = (Pos(1) + 1.3*Pos(3));
y2 = (Pos(2) + 0.6*Pos(4));
annotation('arrow',[x1 x2], [y1 y2], 'LineWidth', 2);