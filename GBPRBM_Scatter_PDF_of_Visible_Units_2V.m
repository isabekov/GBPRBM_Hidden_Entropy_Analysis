function opt = GBPRBM_Scatter_PDF_of_Visible_Units_2V(varargin)
% This function generates samples from probability density function p(v1,v2)
% of the Gaussian-Bipolar Restricted Boltzmann Machine (GBPRBM) with 
% two visible units. It uses scatter plot to visualize the samples.

% Usage:
% Just run without arguments:
%             GBPRBM_Scatter_Probability_of_Visible_Units_2V
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Scatter_Probability_of_Visible_Units_2V(opt)


% Author: Altynbek Isabekov
% E-mail: aisabekov [at] ku.edu.tr

%% ================== Setting Parameters of the GBPRBM Model ==============
if ~isempty(varargin)
    if length(varargin) >=1
        if isstruct(varargin{1})
            opt = varargin{1};
            % Unpack structure
            unpack_struct(opt);
            if V~=2
                error('Number of visible units is not equal to 2. Quitting.')
            end
            if size(W,1)~=2
                error('Matrix W is not of size 2xH. Quitting.')
            end
            if length(sigma_v)~=2
                error('Vector sigma_v is not of size 2x1. Quitting.')
            end     
            if length(b_v)~=2
                error('Vector b_v is not of size 2x1. Quitting.')
            end
        end
    end
else
    %% Usage: pass structure "opt" as an argument to this function
    % Number of visible units
    opt.V = 2;
    % Number of hidden units
    opt.H = 3;
    % Model geometry
    opt.W = [9 2  3;
            2 -6 3];
    % Visible bias
    opt.b_v = [8 12]';
    % Standard deviations for visible units
    opt.sigma_v = [2 2]';
    % Covariance matrix
    opt.S = diag(opt.sigma_v.^2);
    % Inverse of the "covariance matrix"
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    % W = [ 10.0000   -1.3892    6.8200;
    %          0    7.8785   -7.3135  ];
    %     h_x = [0; x2; x3];
    %     b_h_act_x = -W(:,1)'*Sinv*(b_v + alpha*W*h_x);
    %     h_y = [y1; 0; y3];
    %     b_h_act_y = -W(:,2)'*Sinv*(b_v + alpha*W*h_y);
    %     h_z = [z1; z2; 0];
    %     b_h_act_z = -W(:,3)'*Sinv*(b_v + alpha*W*h_z);
    % Span ratio for all axes
    opt.VRatio = [0.4 0.4];
    % Number of dots to be scattered
    opt.N_Samples = 5000;
    % Number of samples used to quantize X-axis
    opt.Nv1 = 250;
    % Number of samples used to quantize Y-axis
    opt.Nv2 = 251;
    % Put labels
    opt.put_labels = true;    
    opt.Color = {'red', 'green', 'blue', 'magenta'};
    % Unpack structure
    unpack_struct(opt);  
    b_h =  (-b_v'*Sinv*W)';
    b_h(3) = -20.250000;
    b_h(1) = -17.250000; 
    b_h(2) = 15.500000;
    opt.b_h =  b_h;
end

if exist('h_axes','var')==1
   axes(h_axes); 
else
   figure('Name', 'GBPRBM Probability of visible units, p(v1,v2)', 'NumberTitle', 'Off');
end

if exist('plot_decision_regions','var')==0
    plot_decision_regions = true;
end

if exist('FontSize','var')==0
    FontSize = 10;
end

%% ======= Determination of the range of  X- and Y-axes ======
% Range is necessary to set proper XLim and YLim values.
% If all variables listed in "Var" exist, no need to determine boundaries -
% they are already given.
Var = {'XLimLow ', 'XLimHigh', 'YLimLow', 'YLimHigh'};
% If all variables listed in "Var" exist, no need to determine boundaries
% They are already given.
if  ~(sum(ismember(Var, who)) == length(Var))   
    % If variables listed in "Var" do not exist then calculated boundaries
    [opt] = GBPRBM_Visible_Units_Span(opt);
    % Cast calculated span for each of the axes
    XLimLow  = opt.v(1).AxesLimMin;
    opt.XLimLow = XLimLow;
    XLimHigh = opt.v(1).AxesLimMax;
    opt.XLimHigh = XLimHigh;
    YLimLow  = opt.v(2).AxesLimMin;
    opt.YLimLow = YLimLow;
    YLimHigh = opt.v(2).AxesLimMax;
    opt.YLimHigh = YLimHigh;
end

% Calculate probability of the hidden units
[opt.p_h opt.TT] = GBPRBM_Calculate_Probability_of_Hidden_Units(opt);

%% Inverse CDF Sampling
v1 = GBPRBM_Inverse_CDF_Sampling_X(opt);
opt.v1 = v1;
v2 = GBPRBM_Inverse_CDF_Sampling_Y_given_X(opt);

% Area of each marker
S = 1.5*ones(1,N_Samples);
% Renderer is set to Z-buffer
set(gcf,'Renderer','zbuffer');
% Scatter dots distributed according to p(v1,v2,v3)
scatter(v1,v2, S,'filled', 'Marker', '.', 'MarkerFaceColor', 'blue', ...
                     'MarkerEdgeColor', 'blue');
opt.v1 = v1;
opt.v2 = v2;
xlabel('v_1', 'FontSize', FontSize);
ylabel('v_2', 'FontSize', FontSize);
opt.h_title = title('p(v_1,v_2)', 'FontSize', FontSize);
set(gca, 'XLim', [XLimLow XLimHigh]);
set(gca, 'YLim', [YLimLow YLimHigh]);
set(gca, 'PlotBoxAspectRatioMode', 'manual');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
set(gca, 'Box', 'on');
daspect([1 1 1]);

%% Plot geometry of the model (location of the centroids)
GBPRBM_Visible_Units_Plot_Geometry(opt);

%% Plot decision regions for p(h_j|v) for j=1...H
if plot_decision_regions==true
    GBPRBM_Plot_Decision_Region_for_Visible_Units(opt);
end