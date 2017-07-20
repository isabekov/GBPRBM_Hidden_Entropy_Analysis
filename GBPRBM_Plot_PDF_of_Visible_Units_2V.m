function opt = GBPRBM_Plot_PDF_of_Visible_Units_2V(varargin)
% This function plots probability density function p(v1,v2) of the
% Gaussian-Bipolar Restricted Boltzmann Machine (GBPRBM) with two visible 
% units.

% Usage:
% Just run without arguments:
%             GBPRBM_Plot_Probability_of_Visible_Units_2V
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Plot_Probability_of_Visible_Units_2V(opt)


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
            if ~(exist('plot_type', 'var'))
                plot_type = 'imagesc'; 
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
            2 -6  3];
    % Visible bias
    opt.b_v = [8 12]';
    % Standard deviations for visible units
    opt.sigma_v = [2 2]';
    % Covariance matrix
    opt.S = diag(opt.sigma_v.^2);
    % Inverse of the "covariance matrix"
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    % Span ratio for all axes
    opt.VRatio =  [0.4 0.4];
    % Number of samples used to quantize X-axis
    opt.Nv1 = 190;
    % Number of samples used to quantize Y-axis
    opt.Nv2 = 180;
    % Plot type
    opt.plot_type = 'imagesc';     
    % Put labels (needed for GBPRBM_Visible_Units_Plot_Geometry)
    opt.put_labels = true;
    opt.Color = {[0, 0.6, 0.6],'red', [0, 0.7, 0], [0.8706, 0.4902, 0], 'green'};    
    opt.plot_decision_regions = true;
    opt.FontSize = 10;
    % Unpack structure
    unpack_struct(opt);   
   
    b_h =  (-b_v'*Sinv*W)';
    b_h(3) = -20.250000;
    b_h(1) = -17.250000; 
    b_h(2) = 15.500000;
    opt.b_h =  b_h;
    opt.j = 2;
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
if  ~(sum(ismember(Var, who)) == length(Var))   
    % If variables listed in "Var" do not exist then calculated boundaries
    [opt] = GBPRBM_Visible_Units_Span(opt);
    % Cast calculated span for each of the axes
    XLimLow  = opt.v(1).AxesLimMin;
    XLimHigh = opt.v(1).AxesLimMax;
    YLimLow  = opt.v(2).AxesLimMin;
    YLimHigh = opt.v(2).AxesLimMax;
end

% PDF p(v) will be computed for each point on the coordinate grid "v". 
v = cell(1,2);

% Linearly spaced coordinate vectors v{1} and v{2} should be of column type
v{1} = linspace(XLimLow, XLimHigh,  Nv1)';
v{2} = linspace(YLimLow, YLimHigh,  Nv2)';
% Truth table for the hidden units
% Superimpose these two matrices
%         Values: "0" and "-1"            Values: "+1" and "0"
TT = double(dec2bin(0:(2^H)-1)-'1') +  double(dec2bin(0:(2^H)-1)-'0');

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
% 1) Compute the part of the p(v) equation which is only dependent on
%    visible units. It will be denoted as "A" and will be used multiple times
%    in subsequent calculations. 
%    Formula: A = 1/Z*exp(-1/2(v - b_v)'*Sinv*(v - b_v))
% 2) Compute the part of the p(v) equation which dependends both on
%    visible vectors and hidden parameters (W and b_h in this case).
%    It will be denoted as "B" and will be used multiple times in subsequent
%    calculations.
%    Formula: B = v'*Sinv*W + b_h' 
% 3) Calculate the sum of all terms exp(B*h), evaluated at all configurations
%    of the hidden vector "h". There are 2^H such configurations.
% 4) Multiply term "A" with previously calculated sum of exponentials to 
%    obtain p(v).
%
% The goal is to avoid redundant calculations, considering that arithmetical
% operationgs in double-precision floating-point format are costly.
% The approach is "to compute once, use many times".

% Number of points in the "v1" coordinate grid
Lv1 = length(v{1});
% Number of points in the "v2" coordinate grid
Lv2 = length(v{2});

% A template for the "A" part, independent of the hidden units
Vis_Only = cell(1,2);             
% A template for the "B" part, which is a function of both visible units "v"
% and hidden parameters ("W", "b_h")
Vis_Hid = cell(1,2);

for i=1:V
    % The lengths of v{1} and v{2} are different!
    Vis_Only{i} = -(v{i}-b_v(i)).^2/(2*(sigma_v(i))^2);   
    % Vis_Hid is a matrix of size Lv2 x H, obtained through an outer product
    Vis_Hid{i}  = v{i}./(sigma_v(i))^2*W(i,:); 
end
% Create a grid
% Replicate column vector Vis_Only{1} along dimension of v2
% Replicate row    vector Vis_Only{2} along dimension of v1
Asum = (  repmat(Vis_Only{1} ,  1 , Lv2) ...
        + repmat(Vis_Only{2}', Lv1,  1 ));

% The "B" part, which is a function of both visible units "v"
% and hidden parameters ("W", "b_h")
B = zeros(Lv1,Lv2,H);
for j=1:H
    B(:,:,j) = repmat(Vis_Hid{1}(:,j),  1 ,Lv2) + ...
             + repmat(Vis_Hid{2}(:,j)',Lv1, 1 ) + b_h(j);
end

a_max = -Inf;
% For every configuration of the hidden layer  (2^H different configurations)
for g=TT' % Hidden vector "g" of size: (H x 1)
    a = 1/2*(2*b_v + W*g)'*Sinv*W*g + b_h'*g;
    if a > a_max
       a_max = a; 
    end   
end

% Initial sum
ExpB=0;
% For every configuration of the hidden layer  (2^H different configurations)
for g=TT'   % Hidden vector g of size: (H x 1)
   SumTmp = zeros(Lv1,Lv2);
   for j=1:H
       SumTmp = SumTmp + B(:,:,j)*g(j);
   end
   ExpB = ExpB + exp(Asum + SumTmp - a_max);
end
%% ================ Computation of partition function Z ===============

% Initial sum
Z=0;
% For every configuration of the hidden layer  (2^H different configurations)
for g=TT' % Hidden vector "g" of size: (H x 1)
    Z = Z + exp(1/2*(2*b_v + W*g)'*Sinv*W*g + b_h'*g - a_max);
end
Z = Z * sqrt((2*pi)^V*det(S));
% ============= End of computation of the partition function Z ============
% Probability of observing visible units, p(v) = p(v1,v2) calculated
% for given values of v1 and v2
p_v = 1/Z*ExpB;

clear('A', 'B', 'ExpB', 'SumTmp');

if strcmp(plot_type, 'imagesc')
    imagesc(v{1},v{2},p_v');
    opt.h_cbar = colorbar;
    daspect([1 1 1]);
    set(gca, 'PlotBoxAspectRatio',[1 1 1]);
    colorbar;
elseif strcmp(plot_type, 'surface')
    surf(v{1},v{2},p_v');
    axis tight;
end

opt.h_title = title('p(v_1,v_2)', 'FontSize', FontSize);
xlabel('v_1', 'FontSize', FontSize);
ylabel('v_2', 'FontSize', FontSize);
set(gca, 'XLim', [XLimLow XLimHigh]);
set(gca, 'YLim', [YLimLow YLimHigh]);
set(gca, 'YDir', 'normal');
set(gca, 'Box', 'on');
set(gca, 'FontSize', FontSize);

%% Plot geometry of the model (location of the centroids)
GBPRBM_Visible_Units_Plot_Geometry(opt);

%% Plot decision regions for p(h_j|v) for j=1...H
if plot_decision_regions==true
    GBPRBM_Plot_Decision_Region_for_Visible_Units(opt);
end
