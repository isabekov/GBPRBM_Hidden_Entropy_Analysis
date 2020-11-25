function PaperFig_HEntropy_1H
% Number of hidden units
opt.H = 1;

opt.V = 3;      
% Model geometry
opt.W = [10 -6 1]';
% Visible bias
opt.b_v = [8 5 3]';
% Standard deviations for visible units
opt.sigma_v = 1.5*[1 1 1]';
% Covariance matrix
opt.S = diag(opt.sigma_v.^2);
% Inverse of the "covariance matrix"
opt.Sinv = diag(1./(opt.sigma_v.^2));
opt.b_h_range = 5;
unpack_struct(opt);

figure('Name', 'Hidden entropy vs. Hidden bias for H=1', ...
       'Units', 'pixels', 'NumberTitle', 'Off', ...
       'Position', [ 189.6000  158.4000  500.0000  140.800]);

opt.plot_type = 'imagesc';
opt.HRatio = 0.6;
opt.N_Samples = 250;

% ================= Empirical Computation =====================
opt.h_axes = subplot(121);
opt.h_axes_emp = opt.h_axes;
opt = GBPRBM_Plot_HEntropy_vs_HBias_1H_Experiment(opt);
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(a) ' str]);

% ================= Analytical Computation ====================
opt.h_axes = subplot(122);
opt.h_axes_thr = opt.h_axes;
opt = GBPRBM_Plot_HEntropy_vs_HBias_1H_Theory(opt);
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(b) ' str]);

Pos = get(opt.h_axes_emp, 'Position');
set(opt.h_axes_emp, 'Position', [Pos(1) 1.05*Pos(2) Pos(3) 1*Pos(4)]);

Pos = get(opt.h_axes_thr, 'Position');
set(opt.h_axes_thr, 'Position', [Pos(1) 1.05*Pos(2) Pos(3) 1*Pos(4)]);
