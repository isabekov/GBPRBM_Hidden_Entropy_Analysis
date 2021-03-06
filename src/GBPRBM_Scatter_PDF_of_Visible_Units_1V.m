function opt = GBPRBM_Scatter_PDF_of_Visible_Units_1V(varargin)
% This function generates samples from probability density function p(v1)
% of the Gaussian-Bipolar Restricted Boltzmann Machine (GBPRBM) with 
% a single visible unit. Histogram is used to visualize the samples.

% Usage:
% Just run without arguments:
%             GBPRBM_Scatter_Probability_of_Visible_Units_1V
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Scatter_Probability_of_Visible_Units_1V(opt)


% Author: Altynbek Isabekov
% E-mail: aisabekov [at] ku.edu.tr

%% ================== Setting Parameters of the GBPRBM Model ==============
if ~isempty(varargin)
    if length(varargin) >=1
        if isstruct(varargin{1})
            opt = varargin{1};
            % Unpack structure
            unpack_struct(opt);
            if V~=1
                error('Number of visible units is not equal to 1. Quitting.')
            end
            if size(W,1)~=1
                error('Matrix W is not of size 1xH. Quitting.')
            end
            if length(sigma_v)~=1
                error('sigma_v is not scalar. Quitting.')
            end   
            if length(b_v)~=1
                error('b_v is not scalar. Quitting.')
            end
        end
    end
else
    %% Usage: pass structure "opt" as an argument to this function
    %  Structure "opt" can have the following values
    % Number of visible units
    opt.V = 1;
    % Number of hidden units
    opt.H = 3;
    % Model geometry
    opt.W = [9 2 3];
    % Visible bias
    opt.b_v = 8;
    % Standard deviations for the visible units
    opt.sigma_v = 2;
    % Covariance matrix
    opt.S = diag(opt.sigma_v.^2);
    % Inverse of the "covariance matrix"
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    % Span ratio for all axes
    opt.VRatio = 0.4;
    % Number of samples used to quantize X-axis
    opt.Nv1 = 190;
    % Put labels (needed for GBPRBM_Visible_Units_Plot_Geometry)
    opt.put_labels = true;
    % Colors for drawing geometry of the model
    opt.Color = {[0, 0.3, 0.6],'red', [0, 0.7, 0], [0.8706, 0.4902, 0], 'green'};
    opt.FontSize = 10;
    % Hidden bias
    opt.b_h = zeros(opt.H,1); % preallocation
    % Entropy = 2 bits   
    % opt.b_h(1) = -6.750000;
    % opt.b_h(2) = -7.000000;
    % opt.b_h(3) = -14.250000; 
    
    % Entropy = 2 bits   
    % opt.b_h(1) = -15.750000; 
    % opt.b_h(2) = 2.000000;
    %opt.b_h(3) = -14.250000; 
    
    % Entropy = 2 bits   
    % opt.b_h(1) = -6.750000;
    % opt.b_h(2) = -10.000000;
    % opt.b_h(3) = -11.250000;    
    
    % Entropy = 1.585 bits     
    opt.b_h(1) = -6.750000;
    opt.b_h(2) = -10.000000;    
    opt.b_h(3) = -11.250000;
     
    % Number of dots to be scattered
    opt.N_Samples = 5000;
    % Number of samples used to quantize X-axis
    opt.Nv1 = 250;
    % Number of bins for histogram
    opt.Nbins = 30;
   
    % Unpack structure
    unpack_struct(opt);    
end

if exist('h_axes','var')==1
   axes(h_axes); 
else
   figure('Name', 'GBPRBM Probability of visible unit, p(v1)', 'NumberTitle', 'Off');
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
Var = {'XLimLow ', 'XLimHigh'};
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
end

% Calculate probability of the hidden units
[opt.p_h opt.TT] = GBPRBM_Calculate_Probability_of_Hidden_Units(opt);

%% Inverse CDF Sampling
v = GBPRBM_Inverse_CDF_Sampling_X(opt);

% Renderer is set to Z-buffer
%set(gcf,'Renderer','zbuffer');
% Scatter dots distributed according to p(v)
[N,X] = hist(v, Nbins);
Xd = diff(X);
binwidth = Xd(ceil(Nbins/2));
bar(X,N/(sum(N)*binwidth), 'EdgeColor', 'blue');
set(gca, 'XLim', [XLimLow XLimHigh]);
set(gca, 'FontSize', FontSize);
xlabel('v_1', 'FontSize', FontSize);
hold on;
opt.h_title = title('p(v_1)', 'FontSize', FontSize);
set(gca, 'XLim', [XLimLow XLimHigh]);

%% Plot geometry of the model (location of the centroids)
GBPRBM_Visible_Units_Plot_Geometry(opt);

%% Plot decision regions for p(h_j|v) for j=1...H
if plot_decision_regions==true
    GBPRBM_Plot_Decision_Region_for_Visible_Units(opt);
end