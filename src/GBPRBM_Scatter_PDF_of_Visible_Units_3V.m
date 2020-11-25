function opt = GBPRBM_Scatter_PDF_of_Visible_Units_3V(varargin)
% This function generates samples from probability density function p(v1,v2,v3)
% of the Gaussian-Bipolar Restricted Boltzmann Machine (GBPRBM) with
% three visible units. It uses scatter plot to visualize the samples.

% Usage:
% Just run without arguments:
%             GBPRBM_Scatter_Probability_of_Visible_Units_3V
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Scatter_Probability_of_Visible_Units_3V(opt)


% Author: Altynbek Isabekov
% E-mail: aisabekov [at] ku.edu.tr

%% ================== Setting Parameters of the GBPRBM Model ==============
if ~isempty(varargin)
    if length(varargin) >=1
        if isstruct(varargin{1})
            opt = varargin{1};
            % Unpack structure
            unpack_struct(opt);
            if V~=3
                error('Number of visible units is not equal to 2. Quitting.')
            end
            if size(W,1)~=3
                error('Matrix W is not of size 3xH. Quitting.')
            end
            if length(sigma_v)~=3
                error('Vector sigma_v is not of size 3x1. Quitting.')
            end   
            if length(b_v)~=3
                error('Vector b_v is not of size 3x1. Quitting.')
            end            
        end
    end
else
    %% Usage: pass structure "opt" as an argument to this function
    % Number of visible units
    opt.V = 3;
    % Number of hidden units
    opt.H = 3;
    % Model geometry
    opt.W = [9 2  1;
            2 -6 1;
            -3 4 5];
 
    % Visible bias
    opt.b_v = [8 12 14]';
    % Standard deviations for visible units
    opt.sigma_v = [2 2 2]';
    % Covariance matrix
    opt.S = diag(opt.sigma_v.^2);
    % Inverse of the "covariance matrix"
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    % W = [ 10.0000   -1.3892    6.8200;
    %          0    7.8785   -7.3135  ];
    opt.b_h =  (-opt.b_v'*opt.Sinv*opt.W)';
    % Span ratio for all axes
    opt.VRatio = [0.4 0.4 0.4];
    % Number of dots to be scattered
    opt.N_Samples = 5000;
    % Number of samples used to quantize X-axis
    opt.Nv1 = 190;
    % Number of samples used to quantize Y-axis
    opt.Nv2 = 180;
    % Number of samples used to quantize Z-axis
    opt.Nv3 = 170;
    % Camera position
    opt.View_Elevation = 38;
    opt.View_Azimuth = 233;
    % Transparency value of the surfaces
    opt.Alpha_Transparency = 0.5;
    % Put labels
    opt.put_labels = true;
    opt.plot_decision_regions = true;
    opt.Color = {'red', 'green', 'blue', 'magenta'};
    % Hidden bias
    opt.b_h =  [-11.000000, 5.500000,-25.500000]';
    % Unpack structure
    unpack_struct(opt);        
end

if exist('h_axes','var')==1
   axes(h_axes); 
else
   figure('Name', 'GBPRBM Probability of visible units, p(v1,v2,v3)', ...
          'Units', 'normalized',  'NumberTitle', 'Off', ...
          'Position', [0.2438    0.0964    0.5220    0.7708]);
end

if exist('plot_decision_regions','var')==0
    plot_decision_regions = true;
end

if exist('FontSize','var')==0
    FontSize = 10;
end

if exist('Alpha_Transparency', 'var')==0
    Alpha_Transparency = 0.5;
end
Var = {'View_Azimuth', 'View_Elevation'};
% Reset azimuth and elevation values if they don't exist
if  (sum(ismember(Var, who)) ~= length(Var)) 
    View_Azimuth = 138;
    View_Elevation = 36;                
end

%% ======= Determination of the range of  X-, Y- and Z-axes ======
% Range is necessary to set proper XLim, YLim and ZLim values.
% If all variables listed in "Var" exist, no need to determine boundaries -
% they are already given.
Var = {'XLimLow ', 'XLimHigh', 'YLimLow', 'YLimHigh', 'ZLimLow', 'ZLimHigh'};
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
    ZLimLow  = opt.v(3).AxesLimMin;
    opt.ZLimLow = ZLimLow;
    ZLimHigh = opt.v(3).AxesLimMax;
    opt.ZLimHigh = ZLimHigh;
end

% Calculate probability of the hidden units
[opt.p_h opt.TT] = GBPRBM_Calculate_Probability_of_Hidden_Units(opt);

%% Inverse CDF Sampling
v1 = GBPRBM_Inverse_CDF_Sampling_X(opt);
opt.v1 = v1;
v2 = GBPRBM_Inverse_CDF_Sampling_Y_given_X(opt);
opt.v2 = v2;
v3 = GBPRBM_Inverse_CDF_Sampling_Z_given_XY(opt);
opt.v3 = v3;
% Area of each marker
S = 1.5*ones(1,N_Samples);

% Scatter dots distributed according to p(v1,v2,v3)
scatter3(v1,v2,v3, S,'filled', 'Marker', 'o', 'MarkerFaceColor', 'blue', ...
                     'MarkerEdgeColor', 'blue');
camproj('perspective');
% Change camera position
view(View_Azimuth, View_Elevation);
h = xlabel('v_1', 'FontSize', FontSize);
% Setting normalized units is necessary in OpenGL rendering mode
set(h, 'Units', 'normalized');
h = ylabel('v_2', 'FontSize', FontSize);
% Setting normalized units is necessary in OpenGL rendering mode
set(h, 'Units', 'normalized');
h = zlabel('v_3', 'FontSize', FontSize);
% Setting normalized units is necessary in OpenGL rendering mode
set(h, 'Units', 'normalized');
opt.h_title = title('p(v_1,v_2,v_3)', 'FontSize', FontSize);
% Setting normalized units is necessary in OpenGL rendering mode
set(opt.h_title, 'Units', 'normalized');

daspect([1 1 1]);
camlight('headlight');
set(gca, 'Box', 'on');
set(gca, 'XLim', [XLimLow XLimHigh]);
set(gca, 'YLim', [YLimLow YLimHigh]);
set(gca, 'ZLim', [ZLimLow ZLimHigh]);
set(gca, 'FontSize', FontSize);

%% Plot geometry of the model (location of the centroids)
GBPRBM_Visible_Units_Plot_Geometry(opt);

%% Plot decision regions for p(h_j|v) for j=1...H
if plot_decision_regions==true
    GBPRBM_Plot_Decision_Region_for_Visible_Units(opt);
end

% Renderer is set to OpenGL
set(gcf,'Renderer','OpenGL');
% Set transparency value for the surfaces (works only with OpenGL renderer)
alpha(Alpha_Transparency);