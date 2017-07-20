function opt = GBPRBM_Plot_HEntropy_vs_HBias_2H_Theory(varargin)
% This function plots hidden entropy of the Gaussian-Bipolar Restricted
% Boltzmann Machine (GBPRBM) with two hidden units as a function of
% hidden bias using analytical solution of the high hidden entropy regions.

% Usage:
% Just run without arguments:
%             GBPRBM_Plot_HEntropy_vs_HBias_2H_Theory
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Plot_HEntropy_vs_HBias_2H_Theory(opt)


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
    opt.W = [9   2; ...
             2  -6];   
    opt.FontSize = 10;
    opt.Color_Map = jet(64);
    % Span ratio for all axes
    opt.HRatio = [0.8 0.8];
    unpack_struct(opt);    
end
if exist('h_axes','var')==1
   axes(h_axes); 
else
   figure('Name', 'Hidden entropy vs. Hidden bias for H=2 (Theory)', ...
          'Units', 'normalized', 'NumberTitle', 'Off');
end

if exist('FontSize','var')==0
    FontSize = 10;
end

%% Determining boundaries
if exist('b_h_lim','var')==1
    disp('Using provided axis limits.');
else      
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
            end
          b_h_bound{j}(idx) = -W(:,j)'*Sinv*(b_v + W*g);
       end
    end
    b_h_lim = cell(1,H);
    for j=1:H
        b_h_lim{j} = zeros(1,2);
        HRange = max(b_h_bound{j}) - min(b_h_bound{j}); 
        % Lower limit for b_h{j}
        b_h_lim{j}(1) = min(b_h_bound{j}) - HRatio(j)*HRange;
        % Upper limit for b_h{j}
        b_h_lim{j}(2) = max(b_h_bound{j}) + HRatio(j)*HRange;    
    end
end

set(gca, 'XLim', b_h_lim{1});
set(gca, 'YLim', b_h_lim{2});

set(gca, 'YDir', 'normal');
set(gcf, 'Renderer', 'painters');

xlabel('b_h(1)', 'FontSize', FontSize);
ylabel('b_h(2)', 'FontSize', FontSize);

opt.h_title = title('Hidden entropy [bits] (theory)', 'FontSize', FontSize);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
set(gca, 'DataAspectRatioMode', 'manual');
set(gca, 'DataAspectRatio',[1 1 1]);
set(gca, 'Box', 'on');
set(gca, 'FontSize', FontSize);
% Center point (most distant configurations separation point)
b_h =  (-b_v'*Sinv*W)';
line(b_h(1)*ones(1,2), b_h(2)*ones(1,2), ...
    'Color', 'blue', 'LineWidth', 2, 'Marker', 'sq');
hold on;

for fxd = 1:H
    % ======== Mapping =========
    % fxd = 1 => a = 2;   
    % fxd = 2 => a = 1;        
    switch fxd
        case 1
            a = 2;
        case 2
            a = 1;
    end
    % ==========================
    % For a = -1 and a = +1
    for h_param=[-1 +1]
        opt.fxd = fxd;
        opt.a = a;
        opt.b_h_lim = b_h_lim;
        h_const = zeros(H,1);
        h_const(a) = h_param(1);        
        opt.h_const = h_const;
        b_h_fxd = -W(:,fxd)'*Sinv*(b_v + W*h_const);
        opt.b_h_fxd = b_h_fxd;
        opt = GBPRBM_Determine_High_HEntropy_Regions_2H(opt);        
        GBPRBM_Plot_High_HEntropy_Regions_2H(opt);
    end
end
h_cbar = colorbar;
colormap(Color_Map);
set(h_cbar, 'YTick', [1, log2(2)/log2(3)*size(Color_Map,1), size(Color_Map,1)]);
set(h_cbar, 'YTickLabel', {0, log2(2), 'log2(3)'});
set(h_cbar, 'FontSize', FontSize);
% Hidden bias which activates two most distant Gaussian components in p(v)
b_h = (-b_v'*Sinv*W)';
C = Color_Map(round(log2(2)/log2(3)*size(Color_Map,1)),:);
line(b_h(1)*[1 1], b_h(2)*[1 1], 'Color', C, ...
                  'LineWidth', 2, 'Marker', 'sq');
text(b_h(1),b_h(2), '${\bf b}_h=-{\bf W}^T{\bf \Sigma}^{-1}{\bf b}_v$', ...
   'BackgroundColor', 'white','Interpreter', 'LaTeX', 'FontSize', FontSize, ...
   'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'left');


function opt = GBPRBM_Determine_High_HEntropy_Regions_2H(opt)
% Input variables:
% opt.W
% opt.b_v
% opt.Sinv
% opt.h_const
% opt.fxd
% opt.a
% opt.b
% opt.b_h_lim

if isstruct(opt)
    unpack_struct(opt);   
end
%% ================== Hidden entropy = log2(3) = 1.585 bits (A is fixed) ===================
fprintf('\n=========================================================\n');
fprintf('Hidden entropy vs. hidden bias (H=2), with one bias fixed:\n')
fprintf('b_h(%i) = %f\n', fxd, b_h_fxd);
fprintf('---------------------------------\n');
fprintf('Hidden entropy = log2(3) bits if\n');

h_fxd = -sign(h_const(a)*W(:,a)'*Sinv*W(:,fxd));
fprintf('    h_j=%+i corresponding to\n',h_fxd);
b_h_a = -W(:,a)'*Sinv*(b_v + h_fxd*W(:,fxd));
fprintf('    b_h(%i) = %f,\n', a, b_h_a);

%% =================== Hidden entropy = 1 bit ====================================
fprintf('---------------------------------\n');
fprintf('Hidden entropy = 1 bit if\n'); 

if h_const(a) == +1           
   fprintf('    b_h(%i) > %f.\n', a, b_h_a);
   A_range = [b_h_a, b_h_lim{a}(2)];   
else
   fprintf('    b_h(%i) < %f.\n', a, b_h_a);
   A_range = [b_h_lim{a}(1), b_h_a]; 
end    
    
% Output argument
opt.A_range = A_range;
opt.b_h_a = b_h_a;


function GBPRBM_Plot_High_HEntropy_Regions_2H(opt)
% Inputs:
% opt.A_range = A_range;
% opt.b_h_a = b_h_a;

if isstruct(opt)
    unpack_struct(opt);   
end
C = Color_Map(round(log2(2)/log2(3)*size(Color_Map,1)),:);
switch fxd
   case 1
       h = line(b_h_fxd*ones(1,2), A_range);                    
   case 2
       h = line(A_range, b_h_fxd*ones(1,2));     
end     
set(h, 'LineWidth', 2, 'Color', C);  
hold on;

%% ================== Hidden entropy = 1.585 bits (A is fixed) ===================
C = Color_Map(round(log2(3)/log2(3)*size(Color_Map,1)),:);
switch fxd
       case 1
           h = line(b_h_fxd*ones(1,2), b_h_a*ones(1,2));
       case 2
           h = line(b_h_a*ones(1,2), b_h_fxd*ones(1,2));
end
set(h, 'LineWidth', 2, 'Marker', 'sq', ...
       'MarkerEdgeColor',C, 'MarkerFaceColor',C);
