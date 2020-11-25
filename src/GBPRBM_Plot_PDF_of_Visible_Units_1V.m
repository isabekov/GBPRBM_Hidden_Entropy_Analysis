function opt = GBPRBM_Plot_PDF_of_Visible_Units_1V(varargin)
% This function plots probability density function p(v1) of the
% Gaussian-Bipolar Restricted Boltzmann Machine (GBPRBM) with a single
% visible unit.

% Usage: 
% Just run without arguments:
%             GBPRBM_Plot_Probability_of_Visible_Units_1V
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Plot_Probability_of_Visible_Units_1V(opt)


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
if  ~(sum(ismember(Var, who)) == length(Var))   
    % If variables listed in "Var" do not exist then calculated boundaries
    [opt] = GBPRBM_Visible_Units_Span(opt);
    % Cast calculated span for each of the axes
    XLimLow  = opt.v(1).AxesLimMin;
    XLimHigh = opt.v(1).AxesLimMax;
end
% Linearly spaced coordinate vector v{1} should be of column type
v = linspace(XLimLow, XLimHigh,  Nv1)';

%% ================ Computation of partition function Z ===============
% Truth table for the hidden units
% Superimpose these two matrices
%         Values: "0" and "-1"            Values: "+1" and "0"
TT = double(dec2bin(0:(2^H)-1)-'1') +  double(dec2bin(0:(2^H)-1)-'0');
% Initial sum
Z=0;
% For every configuration of the hidden layer  (2^H different configurations)
for g=TT' % Hidden vector "g" of size: (H x 1)
    Z = Z + exp(1/2*(2*b_v + W*g)'*Sinv*W*g + b_h'*g);
end
Z = Z * sqrt((2*pi)^V*det(S));
% ============= End of computation of the partition function Z ============

%% ======================== Computation of p(v) ========================
% Three different representations of the p(v) equation are given below:
% Ugly version of the p(v) equation in ASCII
%
%         1        /  1          T   ___ (-1)         \
% p(v) =  —   exp | — — (v - b_v)    \        (v - b_v) | *
%         Z        \  2              /                / 
%                                    ———  
%      ___              ___  (-1)
%      \        //  T   \               T \    \
%    * /   exp ||  v    /        W + b_h   | h  |
%      ———      \\      ———               /    /
%       h

% LaTeX version of the same p(v) equation
% p(\mathbf{v}) = \frac{1}{Z} \exp \left(-\frac{1}{2}(\mathbf{v}-
%                  \mathbf{b}_v)^T \mathbf{\Sigma}^{-1}(\mathbf{v}-
%                  \mathbf{b}_v)\right)
%                  \sum_{\mathbf{h}}  
%                  \exp\left(\left( \mathbf{v}^T \mathbf{\Sigma}^{-1}\mathbf{W}
%                  +\mathbf{b}_h^T\right)\mathbf{h} \right)

% MATLAB pseudo-code for the same p(v) equation, given visible vector "v" 
% for h=TT' % for all hidden configuration of "h"
%  p_v = 1/Z*exp(-1/2*(v - b_v)'*Sinv*(v - b_v))*exp((v'*Sinv*W + b_h')*h)
% end

% Algorithm: 
% 1) For every hidden configuration "g" (there are 2^H of them) compute
%    the exponential energy:
%    ExpE = 1/Z*exp(-1/2(v - b_v)'*Sinv*(v - b_v) + v'*Sinv*(W*g) + (b_h'*g)).
% 2) Accumulate all ExpE.
% 3) Normalize "ExpE" with previously calculated partition function "Z" to obtain p(v).

ExpE = zeros(Nv1,1);
% For every configuration of the hidden layer  (2^H different configurations)
for g=TT'   % Hidden vector g of size: (H x 1)   
   ExpE = ExpE + exp(-(v - b_v).^2/(2*(sigma_v)^2) + v./(sigma_v)^2*(W*g) + b_h'*g);    
end

% Probability of observing visible unit
p_v = 1/Z.*ExpE;

clear('ExpE');
plot(v, p_v);
set(gca, 'YDir', 'normal');
xlabel('v_1', 'FontSize', FontSize);

opt.h_title = title('p(v_1)', 'FontSize', FontSize);
set(gca, 'XLim', [XLimLow XLimHigh]);


%% Plot geometry of the model (location of the centroids)
GBPRBM_Visible_Units_Plot_Geometry(opt);

%% Plot decision regions for p(h_j|v) for j=1...H
if plot_decision_regions==true
   GBPRBM_Plot_Decision_Region_for_Visible_Units(opt);
end