function opt = GBPRBM_Plot_HEntropy_vs_HBias_2H5_Analysis(varargin)
% This function plots hidden entropy of the Gaussian-Bipolar Restricted
% Boltzmann Machine (GBPRBM) with three hidden units as a function of
% hidden biases with one hidden bias set to a constant value. It uses both
% empirical evaluation and analytical solution of the high hidden entropy regions.

% Usage:
% Just run without arguments:
%             GBPRBM_Plot_HEntropy_vs_HBias_2H5_Analysis
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Plot_HEntropy_vs_HBias_2H5_Analysis(opt)


% Author: Altynbek Isabekov
% E-mail: aisabekov [at] ku.edu.tr

%% ================== Setting Parameters of the GBPRBM Model ==============
if ~isempty(varargin)
    if length(varargin) >=1
        if isstruct(varargin{1})
            opt = varargin{1};
            % Unpack structure
            unpack_struct(opt); 
            if H~=3
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
    % Visible bias
    opt.b_v = [10 12]';
    % Standard deviations for visible units
    opt.sigma_v = [2 3]';
    % Covariance-like matrice
    opt.S = diag(opt.sigma_v.^2);
    % Covariance-like matrice's inverse
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    % Number of hidden units
    opt.H = 3;
    % Camera position
    opt.View_Elevation = 38;
    opt.View_Azimuth = 33;
    % Weight matrix of size (V x H)
    opt.W = [-9   2  8; ...
              2  -4  0];   
    opt.Color_Map = jet(64);
    % Span ratio for all axes
    opt.HRatio = [0.8 0.8 0.8];
    opt.FontSize = 10;   
    % Index of the fixed hidden bias
    fxd = 2;
    % Number of samples used for quantization of b_h(x) and b_h(y)
    % where x and y are indices of the hidden bias vector corresponding
    % the hidden biases which are not fixed
    opt.N_Samples = [250 251];
    
    unpack_struct(opt); 

    % ======== Mapping =========
    % fxd = 1 => a = 2; b = 3;   
    % fxd = 2 => a = 1; b = 3;   
    % fxd = 3 => a = 1; b = 2;   
    Dim = 1:H;
    idx =ismember(Dim,fxd);
    ab = num2cell(Dim(~idx));
    [a b] = deal(ab{:});
    h_const = zeros(H,1);
    h_const(a) =+1;
    h_const(b) =-1;    
    b_h_fxd = -W(:,fxd)'*Sinv*(b_v + W*h_const);    
end

if exist('FontSize','var')==0
    FontSize = 10;
end

figure('Name', 'Hidden entropy vs. Hidden Bias (H=3, one bias is fixed)', ...
       'Units', 'normalized', 'Position', [ 0.2343    0.1589    0.4993    0.4688], ...
       'NumberTitle', 'Off');
   
%% Determining boundaries
% Index terms
r = (2.^fliplr([0:H-2]))/2; %#ok<*NBRAK>
s = sum(r) + 1;
b_h_bound = cell(1,H); 
TT = double(dec2bin(0:(2^(H-1))-1)-'1') +  double(dec2bin(0:(2^(H-1))-1)-'0');
for j=1:H
   b_h_bound{j} = zeros(1,2^(H-1));
   for h=TT' 
      idx = r*h + s;
        if j==1
            g = [0; h];
        elseif j==H
            g = [h; 0];
        else
            g = [h(1); 0; h(2)];
        end
      b_h_bound{j}(idx) = -W(:,j)'*Sinv*(b_v + W*g);
   end
end

b_h_lim = cell(1,H);

if isempty(varargin)||~isfield(varargin{1}, 'b_h_lim')
    if (W(:,1)'*W(:,2)) == 0
        range = 10;
        for j=1:H
            b_h_lim{j}(1) = -W(:,j)'*Sinv*b_v - range;
            b_h_lim{j}(2) = -W(:,j)'*Sinv*b_v + range;
        end
    else
        for j=1:H
            b_h_lim{j} = zeros(1,2);
            range = max(b_h_bound{j}) - min(b_h_bound{j}); 
            % Lower limit for b_h{j}
            b_h_lim{j}(1) = min(b_h_bound{j}) - HRatio(j)*range;
            % Upper limit for b_h{j}
            b_h_lim{j}(2) = max(b_h_bound{j}) + HRatio(j)*range;    
        end
    end
else
    b_h_lim = varargin{1}.b_h_lim;
end

%% Sampling, b_h{j} is a row vector
b_h_x = linspace(b_h_lim{a}(1), b_h_lim{a}(2), N_Samples(1))';
b_h_y = linspace(b_h_lim{b}(1), b_h_lim{b}(2), N_Samples(2))';

% Index terms
r = (2.^fliplr([0:H-1]))/2;
s = sum(r) + 1;

% Truth table for hidden layer, e.g. for H=3, TT is given as:
% MSB LSB
%  h1 h2 h3 -> Order is important!
% [-1 -1 -1;
%  -1 -1 +1;
%  -1 +1 -1;
%   .......
%  +1 +1 +1]
TT = double(dec2bin(0:(2^H)-1)-'1') +  double(dec2bin(0:(2^H)-1)-'0');
B = zeros(1,2^H);

exp_b_h_x = exp(b_h_x);
exp_b_h_y = exp(b_h_y);
Nm = zeros(length(b_h_x),length(b_h_y),2^H);
p_h = zeros(length(b_h_x),length(b_h_y),2^H);
Zn  = zeros(length(b_h_x),length(b_h_y));
% For every configuration of the hidden layer  (2^H different configurations)
for h=TT' 
   x = h(a);    
   y = h(b); 
   h_fxd = h(fxd);      
   % Decimal index, e.g. for h=[-1 -1]', idx = 1;  for h=[+1 +1]', idx = 4; 
   idx = r*h + s;
   B(idx) = exp(1/2*(2*b_v + W*h)'*Sinv*W*h + b_h_fxd*h_fxd);
   Nm(:,:,idx) = B(idx) * exp_b_h_x.^x * exp_b_h_y'.^y;
   Zn = Zn + Nm(:,:,idx);
end
% Normalization
for k = 1:2^H
    p_h(:,:,k) = Nm(:,:,k)./Zn;
end
clear Nm Zn exp_b_h_x exp_b_h_y;

Entropy = zeros(length(b_h_x),length(b_h_y));
for h=TT'
    idx = r*h + s;
    Entropy = Entropy - p_h(:,:,idx).*log2(p_h(:,:,idx));
end
clear p_h;

opt.h_axes_1 = subplot(121);
imagesc(b_h_x,b_h_y,Entropy');
h_cbar = colorbar;
set(h_cbar, 'Location', 'SouthOutside');
colormap(Color_Map);
set(h_cbar, 'XTick', [0, log2(2), log2(3), 0.99*log2(4)]);
set(h_cbar, 'XTickLabel', {'0', '1', 'log2(3)', '2'});

set(gca, 'PlotBoxAspectRatio', [1 1 1]);
set(gca, 'DataAspectRatioMode', 'manual');
set(gca, 'DataAspectRatio',[1 1 1]);
set(gca, 'FontSize', FontSize);

set(gca, 'YDir', 'normal');
xlabel(sprintf('b_h(%i)', a), 'FontSize', FontSize);
ylabel(sprintf('b_h(%i)', b), 'FontSize', FontSize);
opt.h_title1 = title('Hidden entropy [bits] (experiment)', 'FontSize', FontSize);
PosExp = get(gca, 'Position');


%% ================== Entropy = 1.585 bits (A is fixed) ===================
hold on; 
opt.h_axes_2 = subplot(122);
set(gca, 'XLim', b_h_lim{1});
set(gca, 'YLim', b_h_lim{2});
set(gca, 'ZLim', b_h_lim{3});
set(gca, 'YDir', 'normal');
set(gca, 'Box', 'on');
set(gcf, 'Renderer', 'painters');

view(View_Azimuth, View_Elevation);
xlabel('b_h(1)', 'FontSize', FontSize);
ylabel('b_h(2)', 'FontSize', FontSize);
zlabel('b_h(3)', 'FontSize', FontSize);
opt.h_title2 = title('Hidden entropy [bits] (theory)', 'FontSize', FontSize);

set(gca, 'PlotBoxAspectRatio', [1 1 1]);
set(gca, 'DataAspectRatioMode', 'manual');
set(gca, 'DataAspectRatio',[1 1 1]);
set(gca, 'FontSize', FontSize);
switch fxd
       case 1
           view(90,0);
       case 2
           view(0,0);
       case 3
           view(0,90);
end
opt.fxd = fxd;
opt.a = a;
opt.b = b;
opt.b_h_lim = b_h_lim;
opt.h_const = h_const;
opt.b_h_fxd = b_h_fxd;
opt = GBPRBM_Determine_High_HEntropy_Regions_3H(opt);        
GBPRBM_Plot_High_HEntropy_Regions_3H(opt);

PosOrig = get(gca, 'Position');
set(gca, 'Position', [PosOrig(1) PosExp(2:4)]);
h_cbar_cp = copyobj(h_cbar, gcf);
pause(1);
Pos = get(h_cbar, 'Position');
set(h_cbar_cp, 'Position', [PosOrig(1) Pos(2:4)]);
