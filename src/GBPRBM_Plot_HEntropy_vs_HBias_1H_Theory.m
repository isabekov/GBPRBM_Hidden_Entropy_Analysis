function opt = GBPRBM_Plot_HEntropy_vs_HBias_1H_Theory(varargin)
% This function plots hidden entropy of the Gaussian-Bipolar Restricted
% Boltzmann Machine (GBPRBM) with a single hidden unit as a function of
% hidden bias using analytical solution of the high hidden entropy regions.

% Usage:
% Just run without arguments:
%             GBPRBM_Plot_HEntropy_vs_HBias_1H_Theory
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Plot_HEntropy_vs_HBias_1H_Theory(opt)


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
                error('Number of hidden units is not equal to 3. Quitting.')
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
    % Number of hidden units
    opt.H = 1;
    % Visible bias
    opt.b_v = [10 10]';
    % Standard deviations for visible units
    opt.sigma_v = [2 3]';
    % Covariance-like matrice (although, it's not)
    opt.S = diag(opt.sigma_v.^2);
    % Covariance-like matrice's inverse
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    % Weight matrix of size (V x H)
    opt.W = [9 2]';
    % Number of samples for b_h
    opt.N_Samples = 150;
    % Range for b_h is [-W'*Sinv*b_v - b_h_range, -W'*Sinv*b_v + b_h_range]
    opt.b_h_range = 5;
    % Unpack structure
    unpack_struct(opt);      
end

if exist('h_axes','var')==1
   axes(h_axes); 
else
   figure('Name', 'Hidden entropy vs. Hidden bias for H=1 (Theory)', 'NumberTitle', 'Off');
end

%% Determining boundaries
b_h_bound= -W'*Sinv*b_v;
b_h_lim = zeros(1,2);
b_h_lim(1) = b_h_bound - b_h_range;
% Upper limit for b_h{j}
b_h_lim(2) = b_h_bound + b_h_range;
% Sampling, b_h{j} is a row vector
b_h = linspace(b_h_lim(1), b_h_lim(2), N_Samples)';
c = b_h - b_h_bound;
Entropy = (c.*tanh(-c) + log(2*cosh(c)))/log(2);
plot(b_h, Entropy, 'LineWidth', 2);
axis tight;
set(gca,'YLim',[0 1.1]);
set(gca, 'YDir', 'normal');
xlabel('b_h(1)');
opt.h_title = title('Hidden entropy [bits] (theory)');
h = line(b_h_bound*[1 1], get(gca, 'YLim'));
set(h, 'LineWidth',2 , 'Color', 'red', 'LineStyle', '--');

text(b_h_bound,diff(get(gca, 'YLim')/2), ...
    '${\bf b}_h=-{\bf W}^T{\bf \Sigma}^{-1}{\bf b}_v$', ...
    'BackgroundColor', 'white','Interpreter', 'LaTeX');

c = b_h(1) - b_h_bound;
% Activated configuration of the hidden unit
h = 2*((c > 0) - 0.5);
text(b_h(1), 0, ['[' sprintf('%+i',h') ']'], ...
    'BackgroundColor', 'yellow', 'FontSize', 8, ...
    'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Left');

c = b_h(end) - b_h_bound;
% Activated configuration of the hidden unit
h = 2*((c > 0) - 0.5);
text(b_h(end),0, ['[' sprintf('%+i',h') ']'], ...
    'BackgroundColor', 'yellow', 'FontSize', 8, ...
    'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Right');
