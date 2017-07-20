function PaperFig_HEntropy_2H5_3H
% Number of hidden units
opt.H = 3;
% Number of visible units
opt.V = 2;
% Model geometry
opt.W = [10 6 2; ...
       -6  4 -2];
% Visible bias
opt.b_v = [8 5]';

% Standard deviations for visible units
%opt.sigma_v = [1.5 1.5]';
opt.sigma_v = 1.5*ones(1,opt.V);% Covariance matrix
% Covariance matrix
opt.S = diag(opt.sigma_v.^2);
% Inverse of the "covariance matrix"
opt.Sinv = diag(1./(opt.sigma_v.^2));
% Span ratio for all axes
opt.VRatio = [0.5 0.5];
opt.Nv1 = 150;
opt.Nv2 = 151;
opt.Nv2 = 152;
opt.Color = {[0, 0.6, 0.6],'red', [0, 0.7, 0], [0.8706, 0.4902, 0], 'green'};
opt.to_draw_lines = 0;
opt.plot_decision_regions = true;
opt.put_labels = true;
figure; opt.Color_Map = jet; close;
unpack_struct(opt);

opt.HRatio = [0.6 0.6 0.6];
opt.N_Samples = [150 151 152];
opt.Entropy_Level = 0.85;
opt.Alpha_Transparency = 0.8;
opt.View_Elevation = 38;
opt.View_Azimuth = 142;

%% ==================== Hidden entropy 2H5 =================
fxd = 2;
% ======== Mapping =========
% fxd = 1 => a = 2; b = 3;   
% fxd = 2 => a = 1; b = 3;   
% fxd = 3 => a = 1; b = 2;   
Dim = 1:H;
idx =ismember(Dim,fxd);
ab = num2cell(Dim(~idx));
[a b] = deal(ab{:});
h_const = zeros(H,1);
h_const(a) =+1;
h_const(b) =+1;   
opt.fxd = fxd;
opt.a = a;
opt.b = b;
opt.h_const = h_const;
opt.b_h_fxd = -W(:,fxd)'*Sinv*(b_v + W*h_const);
opt = GBPRBM_Plot_HEntropy_vs_HBias_2H5_Analysis(opt);
str = get(opt.h_title1, 'String');
set(opt.h_title1, 'String', ['(a) ' str]);
str = get(opt.h_title2, 'String');
set(opt.h_title2, 'String', ['(b) ' str]);
set(opt.h_axes_1, 'XDir', 'reverse');
set(opt.h_axes_2, 'XDir', 'reverse');

hXLabel = get(opt.h_axes_1, 'XLabel');
str = get(hXLabel, 'String');
set(hXLabel, 'String', ['$' str '$']);
set(hXLabel, 'Interpreter', 'LaTeX');

hXLabel = get(opt.h_axes_2, 'XLabel');
str = get(hXLabel, 'String');
set(hXLabel, 'String', ['$' str '$']);
set(hXLabel, 'Interpreter', 'LaTeX');

hYLabel = get(opt.h_axes_1, 'YLabel');
str = get(hYLabel, 'String');
set(hYLabel, 'String', ['$' str '$']);
set(hYLabel, 'Interpreter', 'LaTeX');

% Attention! Z-label is used instead of Y-label
hYLabel = get(opt.h_axes_2, 'ZLabel');
str = get(hYLabel, 'String');
set(hYLabel, 'String', ['$' str '$']);
set(hYLabel, 'Interpreter', 'LaTeX');

%% ==================== Hidden entropy 3H =================
opt = GBPRBM_Plot_HEntropy_vs_HBias_3H_Analysis(opt);
str = get(opt.h_title1, 'String');
set(opt.h_title1, 'String', ['(c) ' str]);
str = get(opt.h_title2, 'String');
set(opt.h_title2, 'String', ['(d) ' str]);
cameratoolbar;
set(gcf, 'Toolbar', 'figure');