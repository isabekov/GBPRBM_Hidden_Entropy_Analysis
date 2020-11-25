function opt = GBPRBM_Plot_HEntropy_vs_HBias_1H_Analysis(varargin)
% This function plots hidden entropy of the Gaussian-Bipolar Restricted
% Boltzmann Machine (GBPRBM) with a single hidden unit as a function of
% hidden bias using both empirical evaluation and analytical solution
% of the high hidden entropy regions.

% Usage: 
% Just run without arguments:
%             GBPRBM_Plot_HEntropy_vs_HBias_1H_Analysis
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Plot_HEntropy_vs_HBias_1H_Analysis(opt)


% Author: Altynbek Isabekov
% E-mail: aisabekov [at] ku.edu.tr

%% ================== Setting Parameters of the GBPRBM Model ==============
if ~isempty(varargin)
    if length(varargin) >=1
        if isstruct(varargin{1})
            opt = varargin{1};
            % Unpack structure
            unpack_struct(opt); 
            if H~=1
                error('Number of hidden units is not equal to 1. Quitting.')
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
    % If the model parameters are not given, use the following example
    % Number of visible units
    opt.V = 2;
    % Visible bias
    opt.b_v = [10 12]';
    % Standard deviations for visible units
    opt.sigma_v = [2 3]';
    % Covariance-like matrice (although, it's not)
    opt.S = diag(opt.sigma_v.^2);
    % Covariance-like matrice's inverse
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    % Number of hidden units
    opt.H = 1;
    % Weight matrix of size (V x H)
    opt.W = [-9  ; ...
              2  ];          
    opt.N_Samples = 100;
    % Range for b_h is [-W'*Sinv*b_v - b_h_range, -W'*Sinv*b_v + b_h_range]
    opt.b_h_range = 5;
    unpack_struct(opt); 
end

figure('Name', 'Hidden Entropy vs. Hidden Bias (Experiment and Theory)', ...
   'Units', 'normalized', 'Position', [0.0198, 0.3021, 0.5488, 0.4036], ...
   'NumberTitle', 'Off');


% ================= Empirical Computation =====================
opt.h_axes = subplot(121);
opt = GBPRBM_Plot_HEntropy_vs_HBias_1H_Experiment(opt);
hold on;

% ================= Theoretical Computation ===================
opt.h_axes = subplot(122);
GBPRBM_Plot_HEntropy_vs_HBias_1H_Theory(opt);
