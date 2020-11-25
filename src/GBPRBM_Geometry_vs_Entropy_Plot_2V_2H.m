function GBPRBM_Geometry_vs_Entropy_Plot_2V_2H(varargin)
% This function plots geometry of the centroid locations in p(v1,v2) of the
% Gaussian-Bipolar Restricted Boltzmann Machine (GBPRBM) with two visible 
% units. Hidden entropy as a function of hidden bias is also plotted for
% the given parameters.

% Usage:
% Just run without arguments:
%             GBPRBM_Geometry_vs_Entropy_Plot_2V_2H
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Geometry_vs_Entropy_Plot_2V_2H(opt)


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
                error('Vector b_v is not of of size 2x1. Quitting.')
            end            
        end
    end
else
    %% Usage: pass structure "opt" as an argument to this function
    % Number of visible units
    opt.V = 2;
    % Visible bias
    opt.b_v = [10 10]';
    % Number of hidden units
    opt.H = 2;
    % Standard deviations for visible units
    opt.sigma_v = [2 2]';
    % Covariance matrix
    opt.S = diag(opt.sigma_v.^2);
    % Inverse of the "covariance matrix"
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    theta = -20/360*2*pi;
    w1= [7 4]';
    w2 = [4 -7]'; 
    
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
    opt.W = [w1 R*w2];
    % Weight matrix of size (V x H)
    % Span ratio for all axes
    opt.VRatio = [0.4 0.1];
    % Span ratio for all axes
    opt.HRatio = [1.2 1.2];   
    % Number of samples used for quantization of b_h(1) and b_h(2)
    opt.N_Samples = [150 151];
    opt.FontSize = 10;
    opt.plot_type = 'imagesc';
    unpack_struct(opt); 
end

if exist('FontSize','var')==0
    FontSize = 10;
end

%% ======= Determination of the range of  X- and Y-axes ======
% Range is necessary to set proper XLim and YLim values.
% If all variables listed in "Var" exist, no need to determine boundaries -
% they are already given.
Var = {'XLimLow ', 'XLimHigh', 'YLimLow', 'YLimHigh'};
if  ~(sum(ismember(Var, who)) == length(Var))   
    % If variables listed in "Var" do not exist then calculated boundaries
    [opt] = GBPRBM_Visible_Units_Span(opt);
    % Cast calculated span for each of the axes
    XLimLow  = opt.v(1).AxesLimMin;
    XLimHigh = opt.v(1).AxesLimMax;
    YLimLow  = opt.v(2).AxesLimMin;
    YLimHigh = opt.v(2).AxesLimMax;
end

figure('Name', 'Model Geometry and Hidden Entropy', 'NumberTitle', 'Off');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [4 4 19 7]);

%% Plot geometry of the centroids' locations
opt.h_axes = subplot(121);
set(gca, 'XLim', [XLimLow XLimHigh]);
set(gca, 'YLim', [YLimLow YLimHigh]);
set(gca, 'PlotBoxAspectRatioMode', 'manual');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
set(gca, 'DataAspectRatioMode', 'manual');
set(gca, 'DataAspectRatio', [1 1 1]);
title('Geometry of p(v_1,v_2)', 'FontSize', FontSize);
set(gca, 'YDir', 'normal');
xlabel('v_1', 'FontSize', FontSize);
ylabel('v_2', 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);
GBPRBM_Visible_Units_Plot_Geometry(opt);
set(gca, 'box', 'on');


%% Plot hidden entropy as a function of hidden bias
opt.h_axes=subplot(122);
set(gca, 'PlotBoxAspectRatioMode', 'manual');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
set(gca,'DataAspectRatio',[1 1 1]);
set(gca, 'PlotBoxAspectRatio',[1 1 1]);
GBPRBM_Plot_HEntropy_vs_HBias_2H_Experiment(opt);