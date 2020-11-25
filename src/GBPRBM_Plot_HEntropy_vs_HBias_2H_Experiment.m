function opt = GBPRBM_Plot_HEntropy_vs_HBias_2H_Experiment(varargin)
% This function plots hidden entropy of the Gaussian-Bipolar Restricted
% Boltzmann Machine (GBPRBM) with two hidden units as a function of
% hidden bias using empirical evaluation of the hidden entropy function.

% Usage:
% Just run without arguments:
%             GBPRBM_Plot_HEntropy_vs_HBias_2H_Experiment
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Plot_HEntropy_vs_HBias_2H_Experiment(opt)


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
    opt.b_v = [10 10]';
    % Standard deviations for visible units
    opt.sigma_v = [2 3]';
    % Covariance-like matrice
    opt.S = diag(opt.sigma_v.^2);
    % Covariance-like matrice's inverse
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    % Number of hidden units
    opt.H = 2;

    % Weight matrix of size (V x H)
    opt.W = [9  2 ; ...
             2 -6 ];
    % Number of samples used for quantization of b_h(1) and b_h(2)
    opt.N_Samples = [150 151];
    opt.plot_type = 'imagesc';
    opt.FontSize = 10;
    % Span ratio for all axes
    opt.HRatio = [1 1];    
    % Unpack structure
    unpack_struct(opt); 
end

if exist('h_axes','var')==1
   axes(h_axes); 
else
   figure('Name', 'Hidden entropy vs. Hidden bias for H=2 (Experiment)', 'NumberTitle', 'Off');
end

if exist('FontSize','var')==0
    FontSize = 10;
end

%% Determining boundaries
% Index terms
r = (2.^fliplr([0:H-2]))/2; %#ok<NBRAK>
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
        end
      b_h_bound{j}(idx) = -W(:,j)'*Sinv*(b_v + W*g);
   end
end
% Plot ratios for h1 and h2
b_h_lim = cell(1,H);
b_h     = cell(1,H);

if isempty(varargin)||~isfield(varargin{1}, 'b_h_lim')
    if (W(:,1)'*W(:,2)) == 0
        for j=1:H
            b_h_lim{j}(1) = -W(:,j)'*Sinv*b_v - HRange(j);
            b_h_lim{j}(2) = -W(:,j)'*Sinv*b_v + HRange(j);
        end
    else
        HRange = zeros(H,1);
        for j=1:H
            b_h_lim{j} = zeros(1,2);            
            HRange(j) = max(b_h_bound{j}) - min(b_h_bound{j}); 
            % Lower limit for b_h{j}
            b_h_lim{j}(1) = min(b_h_bound{j}) - HRatio(j)*HRange(j);
            % Upper limit for b_h{j}
            b_h_lim{j}(2) = max(b_h_bound{j}) + HRatio(j)*HRange(j);    
        end
    end
else
    b_h_lim = varargin{1}.b_h_lim;
end

for j=1:H
    % Sampling, b_h{j} is a row vector
    b_h{j} = linspace(b_h_lim{j}(1), b_h_lim{j}(2), N_Samples(j))';
end
% Index terms
r = (2.^fliplr([0:H-1]))/2; %#ok<NBRAK>
s = sum(r) + 1;

% Truth table for hidden layer, e.g. for H=2, TT is given as:
% MSB LSB
%  h1 h2 Order is important!
% [-1 -1;
%  -1 +1;
%  +1 -1;
%  +1 +1]
TT = double(dec2bin(0:(2^H)-1)-'1') +  double(dec2bin(0:(2^H)-1)-'0');
B = zeros(1,2^H);

exp_b_h{1} = exp(b_h{1});
exp_b_h{2} = exp(b_h{2});
Nm = zeros(length(b_h{1}),length(b_h{2}),2^H);
p_h = zeros(length(b_h{1}),length(b_h{2}),2^H);
Zn  = zeros(length(b_h{1}),length(b_h{2}));
% For every configuration of the hidden layer  (2^H different configurations)
for h=TT' 
   % Decimal index, e.g. for h=[-1 -1]', idx = 1;  for h=[+1 +1]', idx = 4; 
   idx = r*h + s;
   B(idx) = exp(1/2*(2*b_v + W*h)'*Sinv*W*h);
   Nm(:,:,idx) = B(idx)*exp_b_h{1}.^h(1)*exp_b_h{2}'.^h(2);
   Zn = Zn + Nm(:,:,idx);
end
% Normalization
for k = 1:2^H
    p_h(:,:,k) = Nm(:,:,k)./Zn;
end
clear Nm;

Entropy = zeros(length(b_h{1}),length(b_h{2}));
for h=TT'
    idx = r*h + s;
    Entropy = Entropy - p_h(:,:,idx).*log2(p_h(:,:,idx));
end

switch plot_type 
    case 'imagesc'
        imagesc(b_h{1},b_h{2},Entropy');
    case 'surface'
        surf(b_h{1},b_h{2},Entropy');
end
h_cbar = colorbar;
YTick = [0.001, 1, 0.99*log2(3)];
set(h_cbar, 'YTick', YTick);
YTickLabel = {'0','1', 'log2(3)'};
set(h_cbar, 'YTickMode', 'manual');
set(h_cbar, 'YTickLabel', YTickLabel);
set(h_cbar, 'FontSize', FontSize);
opt.h_cbar = h_cbar;

set(gca, 'PlotBoxAspectRatioMode', 'manual');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
set(gca, 'DataAspectRatioMode', 'manual');
set(gca, 'DataAspectRatio',[1 1 1]);
set(gca, 'FontSize', FontSize);

set(gca, 'YDir', 'normal');
xlabel('b_h(1)', 'FontSize', FontSize);
ylabel('b_h(2)', 'FontSize', FontSize);
opt.h_title = title('Hidden entropy [bits] (experiment)', 'FontSize', FontSize);

% Hidden bias which activates two most distant Gaussian components in p(v)
b_h = (-b_v'*Sinv*W)';
line(b_h(1)*[1 1], b_h(2)*[1 1], 'Color', 'magenta', ...
                  'LineWidth', 2, 'Marker', 'sq');
text(b_h(1),b_h(2), '${\bf b}_h=-{\bf W}^T{\bf \Sigma}^{-1}{\bf b}_v$', ...
   'BackgroundColor', 'white','Interpreter', 'LaTeX', 'FontSize', FontSize, ...
   'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'left');              
