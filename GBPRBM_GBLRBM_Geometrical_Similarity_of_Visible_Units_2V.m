function GBPRBM_GBLRBM_Geometrical_Similarity_of_Visible_Units_2V(varargin)
% This function plots similar geometries of the Gaussian-Bipolar Restricted
% Boltzmann Machine (GBPRBM) and the Gaussian-Bernoulli Restricted Boltzmann 
% Machine (GBLRBM). Conversion of the model parameters for GBLRBM is provided.

% Usage:
% Just run without arguments:
%             GBLRBM_GBPRBM_Geometrical_Similarity_of_Visible_Units_2V
% Or pass a structure with parameters, look down below for the details
%             GBLRBM_GBPRBM_Geometrical_Similarity_of_Visible_Units_2V(opt)


% Author: Altynbek Isabekov
% E-mail: aisabekov [at] ku.edu.tr

%% ================== Setting Parameters of the GBPRBM Model ==============
if ~isempty(varargin)
    if length(varargin) >=1
        if isstruct(varargin{1})
            opt = varargin{1};
            V = varargin{1}.V;
            H = varargin{1}.H;
            W = varargin{1}.W;
            if V~=2
                error('Number of visible units is not equal to 2. Quitting.')
            end
            if size(W,1)~=2
                error('Matrix W is not of size 2xH. Quitting.')
            end
            b_v = varargin{1}.b_v;       
        end
    end
else
    % Number of visible units
    opt.V = 2;
    % Number of hidden units
    opt.H = 3;
    % Model geometry
    opt.W = [8 2  -3;
            4 -6  +2];
    % Visible bias
    opt.b_v = [8 12]';
    % Standard deviations for visible units
    opt.sigma_v = [2 3]';
    % Covariance matrix
    opt.S = diag(opt.sigma_v.^2);
    % Inverse of the "covariance matrix"
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    % Span ratio for all axes
    opt.VRatio = [0.2 0.2];
    opt.FontSize = 10;
    % Unpack structure
    unpack_struct(opt);   
end

%% Gaussian-Bipolar Restricted Boltzmann Machine 
% Calculating limits of the axes

opt = GBPRBM_Visible_Units_Span(opt);

% Add some more space beyond original XLim and YLim range
XLimLow  = opt.v(1).min - VRatio(1)*opt.v(1).range;
XLimHigh = opt.v(1).max + VRatio(1)*opt.v(1).range;

YLimLow  = opt.v(2).min - VRatio(2)*opt.v(2).range;
YLimHigh = opt.v(2).max + VRatio(2)*opt.v(2).range;


figure('Name', 'GBPRBM and GBLRBM', 'NumberTitle', 'Off', 'Units', 'normalized', ...
       'Position', [ 0.1347    0.2253    0.7430    0.5182]);
subplot(121);
set(gca, 'LooseInset', [0,0,0,0]);
set(gca, 'XLim', [XLimLow XLimHigh]);
set(gca, 'YLim', [YLimLow YLimHigh]);
set(gca,'DataAspectRatio',[1 1 1]);
set(gca, 'PlotBoxAspectRatio',[1 1 1]);
title('Gaussian-Bipolar RBM', 'FontSize', FontSize);
set(gca, 'YDir', 'normal');
set(gca, 'FontSize', FontSize);
xlabel('v_1', 'FontSize', FontSize);
ylabel('v_2', 'FontSize', FontSize);

%% Gaussian-Bipolar Restricted Boltzmann Machine
GBPRBM_Visible_Units_Plot_Geometry(opt);

%% Gaussian-Bernoulli Restricted Boltzmann Machine
%opt.b_v = b_v + (W*(-ones(H,1)) - W*ones(H,1))/2;
opt.b_v = b_v + W*(-ones(H,1));
opt.W = W*2;
subplot(122);
set(gca, 'LooseInset', [0,0,0,0]);
set(gca, 'XLim', [XLimLow XLimHigh]);
set(gca, 'YLim', [YLimLow YLimHigh]);
set(gca,'DataAspectRatio',[1 1 1]);
set(gca, 'PlotBoxAspectRatio',[1 1 1]);

title('Gaussian-Bernoulli RBM', 'FontSize', FontSize);
set(gca, 'YDir', 'normal');
set(gca, 'FontSize', FontSize);
xlabel('v_1', 'FontSize', FontSize);
ylabel('v_2', 'FontSize', FontSize);

GBLRBM_Visible_Units_Plot_Geometry(opt);

