function PaperFig_VPDF_1V_2V_3V
% Number of hidden units
opt.H = 3;
% Weight matrix
W_r = [10  6  2; ...
       -6  4 -2;
        1  3  5];     
% Visible bias
b_v_r = [8 5 3]';
% Standard deviations for visible units
sigma_v_r = 1.5*ones(opt.H,1);

% Span ratio for all axes
opt.VRatio = [0.1 0.1 0.1];

%  Number of samples used to quantize V1,V2 and V3 axes
opt.Nv1 = 150 ;
opt.Nv2 = 151;
opt.Nv3 = 152;

% Colormap for 2-dimensional plot of p(v1,v2)
opt.Color_Map = jet(64);

% Colors for drawing geometry of the model
opt.Color = {[0, 0.3, 0.6],'red', [0, 0.7, 0], [0.8706, 0.4902, 0], 'green'};

opt.plot_decision_regions = true;
opt.put_labels = true;
unpack_struct(opt);  

opt.plot_type = 'imagesc';
figure('Name', 'PDF of Visible Units for V=1,2', 'Units', 'pixels',...
      'Position', [ 373.6000,   23.2000,  401.6000,  506.4000], ...
      'NumberTitle', 'Off');

%% ================= V = 1 =====================
opt.h_axes = subplot(211);
set(opt.h_axes, 'Units', 'normalized');
opt.V = 1;
opt.W = W_r(1,:);
opt.b_v = b_v_r(1);
opt.sigma_v = sigma_v_r(1);
% Covariance matrix
opt.S = diag(opt.sigma_v.^2);
% Inverse of the "covariance matrix"
opt.Sinv = diag(1./(opt.sigma_v.^2));
% Extract variables from structure
unpack_struct(opt);
% Hidden bias
opt.b_h = zeros(opt.H,1);
 h = [0, +1, +1;
     -1,  0, +1;
     -1, -1, 0]';
opt.b_h(1) = -W(:,1)'*Sinv*(b_v + W*h(:,1)); %opt.b_h(1) =  -71.1111;
opt.b_h(2) = -W(:,2)'*Sinv*(b_v + W*h(:,2)); %opt.b_h(2) =         0;
opt.b_h(3) = -W(:,3)'*Sinv*(b_v + W*h(:,3)); %opt.b_h(3) =    7.1111;
 
% Positions of the labels of the decision regions for p(h_j|v) for j=1...H
opt.LabelPosition = [0.87, 0.42, 0.67];

% Plot probability of visible unit, p(v1)
opt = GBPRBM_Plot_PDF_of_Visible_Units_1V(opt);
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(a) ' str]);
hold on;


%% =============== V = 2 ====================
opt.h_axes = subplot(212);
set(opt.h_axes, 'Units', 'normalized');
opt.V = 2;
opt.W = W_r(1:2,:);
opt.b_v = b_v_r(1:2);
opt.sigma_v = sigma_v_r(1:2);
% Covariance matrix
opt.S = diag(opt.sigma_v.^2);
% Inverse of the "covariance matrix"
opt.Sinv = diag(1./(opt.sigma_v.^2));
% Extract variables from structure
unpack_struct(opt);
% Positions of the labels of the decision regions for p(h_j|v) for j=1...H
opt.LabelPosition = [0.77, 0.72, 0.67];
% Hidden bias
opt.b_h = zeros(opt.H,1);
h = [0, +1, +1;
     -1,  0, +1;
     -1, -1, 0]';
opt.b_h(1) = -W(:,1)'*Sinv*(b_v + W*h(:,1)); %opt.b_h(1) = -52.4444;
opt.b_h(2) = -W(:,2)'*Sinv*(b_v + W*h(:,2)); %opt.b_h(2) = -16.0000;
opt.b_h(3) = -W(:,3)'*Sinv*(b_v + W*h(:,3)); %opt.b_h(3) =  13.3333;

% Plot probability of visible units, p(v1,v2)
opt = GBPRBM_Plot_PDF_of_Visible_Units_2V(opt);
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(b) ' str]);
hold on;

%% =============== V = 3 ====================
figure('Name', 'PDF of Visible Units for V=3', ...
       'Units', 'normalized', ... %'Position', [ 1.1559   -0.1628    0.2884    0.8880], ...
       'NumberTitle', 'Off');
opt.h_axes = axes;
set(opt.h_axes, 'Units', 'normalized');
opt.V = 3;
opt.W = W_r(1:3,:);
opt.b_v = b_v_r(1:3);
opt.sigma_v = sigma_v_r(1:3);
% Covariance matrix
opt.S = diag(opt.sigma_v.^2);
% Inverse of the "covariance matrix"
opt.Sinv = diag(1./(opt.sigma_v.^2));
% Extract variables from structure
unpack_struct(opt);
% Hidden bias
opt.b_h = zeros(opt.H,1);
h = [0, +1, +1;
     -1,  0, +1;
     -1, -1, 0]';
opt.b_h(1) = -W(:,1)'*Sinv*(b_v + W*h(:,1)); %opt.b_h(1) = -57.3333;
opt.b_h(2) = -W(:,2)'*Sinv*(b_v + W*h(:,2)); %opt.b_h(2) = -25.3333;
opt.b_h(3) = -W(:,3)'*Sinv*(b_v + W*h(:,3)); %opt.b_h(3) =  15.5556;
  
% Number of sample to be drawn
opt.N_Samples = 1000;

% Scatter samples drawn from the probability of visible units, p(v1,v2,v3)
opt = GBPRBM_Scatter_PDF_of_Visible_Units_3V(opt);
cameratoolbar;
% View angle
view(23, 26);
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(c) ' str]);
hold on;
set(gcf, 'Renderer', 'OpenGL');
% Transparency
alpha(0.5);
set(gca, 'Projection', 'orthographic');

