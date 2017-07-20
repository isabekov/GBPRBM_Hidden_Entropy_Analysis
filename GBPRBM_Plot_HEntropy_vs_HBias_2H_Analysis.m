function opt = GBPRBM_Plot_HEntropy_vs_HBias_2H_Analysis(varargin)
% This function plots hidden entropy of the Gaussian-Bipolar Restricted
% Boltzmann Machine (GBPRBM) with two hidden units as a function of
% hidden bias using both empirical evaluation and analytical solution
% of the high hidden entropy regions.

% Usage: 
% Just run without arguments:
%             GBPRBM_Plot_HEntropy_vs_HBias_2H_Analysis
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Plot_HEntropy_vs_HBias_2H_Analysis(opt)


% Author: Altynbek Isabekov
% E-mail: aisabekov [at] ku.edu.tr

%% ================== Setting Parameters of the GBPRBM Model ==============
if ~isempty(varargin)
    if length(varargin) >=1
        if isstruct(varargin{1})
            opt = varargin{1};
            % Unpack structure
            unpack_struct(opt); 
            if H~=2
                error('Number of hidden units is not equal to 2. Quitting.')
            end
            if sum(size(W)== [V H]) ~= 2
                error('Matrix W is not of size VxH. Quitting.')
            end
            if length(sigma_v)~= V
                error('Vector sigma_v is not of length V. Quitting.')
            end 
            if length(b_v)~= V
                error('Vector b_v is not of length V. Quitting.')
            end                  
        end
    end
 else
    %% Usage: pass structure "opt" as an argument to this function
    % Number of visible units
    opt.V = 2;
    % Visible bias
    opt.b_v = [8 5]';
    % Standard deviations for visible units
    opt.sigma_v = [2 3]';
    % Covariance-like matrice
    opt.S = diag(opt.sigma_v.^2);
    % Covariance-like matrice's inverse
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    % Number of hidden units
    opt.H = 2;
    % Weight matrix of size (V x H)
    opt.W = [-9  2 ; ...
              2  -4 ];   
    % Number of samples used for quantization of b_h(1) and b_h(2)          
    opt.N_Samples = [100 101];
    opt.Color_Map = jet(64);
    opt.FontSize = 10;    
    opt.plot_type = 'imagesc';
    % Span ratio for all axes
    opt.HRatio = [1.2 1.2];    
    unpack_struct(opt); 
end

figure('Name', 'Hidden Entropy vs. Hidden Bias for H=2 (Experiment and Theory)', ...
   'Units', 'normalized', 'Position', [0.0198, 0.3021, 0.5488, 0.4036], ...
   'NumberTitle', 'Off');

% ================= Empirical Computation =====================
opt.h_axes = subplot(121);
opt = GBPRBM_Plot_HEntropy_vs_HBias_2H_Experiment(opt);
hold on;
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(a) ' str]);

% ================= Theoretical Computation ===================
opt.h_axes = subplot(122);
opt = GBPRBM_Plot_HEntropy_vs_HBias_2H_Theory(opt);
str = get(opt.h_title, 'String');
set(opt.h_title, 'String', ['(b) ' str]);
