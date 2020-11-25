function PaperFig_HEntropy_2H_VPDF_V2
% Number of visible units
opt.V = 2;
% Number of hidden units
opt.H = 2;
% Model geometry
%opt.W = [-9  2 ; ...
         %2  -4 ];          
opt.W = [10 6; ...
       -6 4];         
    
opt.b_v = [8 5]';
% Standard deviations for visible units
opt.sigma_v = [1.5 1.5]';
% Covariance matrix
opt.S = diag(opt.sigma_v.^2);
% Inverse of the "covariance matrix"
opt.Sinv = diag(1./(opt.sigma_v.^2));
    % Span ratio for all axes
opt.VRatio = [0.2 0.5];

opt.Nv1 = 350 ;
opt.Nv2 = 350;
figure; opt.Color_Map = jet; close;

opt.Color = {[0, 0.6, 0.6],'red', [0, 0.7, 0], [0.8706, 0.4902, 0], 'green'};
opt.to_draw_lines = 1;
opt.plot_decision_regions = true;
opt.put_labels = true;
unpack_struct(opt);  



opt.HRatio = [1 1];
% Number of samples used for quantization of b_h(1) and b_h(2)          
opt.N_Samples = [100 101];
opt.plot_type = 'imagesc';
figure('Name', 'Hidden Entropy vs. Hidden Bias (Empirical and Theorethical)', ...
   'Units', 'normalized', 'Position', [  0.1003    0.0417    0.6369    0.8268], ...
   'NumberTitle', 'Off');


% ================= Empirical Computation =====================
opt.h_axes = subplot(221);
opt.h_axes_ent_exp = opt.h_axes;
opt = GBPRBM_Plot_HEntropy_vs_HBias_2H_Experiment(opt);
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(a) ' str]);
hold on;
% =============== Theoretical Computation ====================

opt.h_axes = subplot(222);
opt = GBPRBM_Plot_HEntropy_vs_HBias_2H_Theory(opt);
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(b) ' str]);

opt.b_h = zeros(opt.H,1);
opt.b_h = (-opt.b_v'*opt.Sinv*opt.W)';
opt.b_h_md = opt.b_h;



% ================= Most Distant Components ===================
opt.LabelPosition = [0.82, 0.85];
opt.LineWidth = 2;
opt.h_axes = subplot(223);
opt.h_axes_md = opt.h_axes;
opt = GBPRBM_Plot_PDF_of_Visible_Units_2V(opt);
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(c) ' str]);
opt.h_title_md = opt.h_title;

% ================= 3 Active Components ===================
h = [0, -1]';
opt.b_h(1) = -W(:,1)'*Sinv*(b_v + W*h); %opt.b_h(1) = -17.250000;
h = [+1, 0]';
opt.b_h(2) = -W(:,2)'*Sinv*(b_v + W*h); %opt.b_h(2) = 15.500000;
opt.b_h_3conf = opt.b_h;
opt.h_axes = subplot(224);
opt.h_axes_3conf = opt.h_axes;

opt.LabelPosition = [0.82, 0.85];

opt = GBPRBM_Plot_PDF_of_Visible_Units_2V(opt);
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(d) ' str]);

[x1 y1] = dsxy2figxy(opt.h_axes_ent_exp, opt.b_h_md(1), opt.b_h_md(2));
%set(opt.h_title_md, 'Units', 'normalized');
Pos = get(opt.h_title_md, 'Position');
[x2 y2] = dsxy2figxy(opt.h_axes_md, Pos(1), Pos(2));
annotation('arrow',[x1 x2], [y1 y2], 'LineWidth', 2);


[x1 y1] = dsxy2figxy(opt.h_axes_ent_exp, opt.b_h_3conf(1), opt.b_h_3conf(2));
%set(opt.h_title_md, 'Units', 'normalized');
Pos = get(opt.h_title_md, 'Position');
[x2 y2] = dsxy2figxy(opt.h_axes_3conf, Pos(1), Pos(2));
annotation('arrow',[x1 x2], [y1 y2], 'LineWidth', 2);