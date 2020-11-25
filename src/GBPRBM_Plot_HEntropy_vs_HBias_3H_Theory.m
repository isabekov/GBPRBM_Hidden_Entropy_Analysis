function opt = GBPRBM_Plot_HEntropy_vs_HBias_3H_Theory(varargin)
% This function plots hidden entropy of the Gaussian-Bipolar Restricted
% Boltzmann Machine (GBPRBM) with three hidden units as a function of
% hidden bias using analytical solution of the high hidden entropy regions.

% Usage:
% Just run without arguments:
%             GBPRBM_Plot_HEntropy_vs_HBias_3H_Theory
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Plot_HEntropy_vs_HBias_3H_Theory(opt)


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
    opt.W = [-9  2  8; ...
              2  -4   0];   
    opt.FontSize = 10;
    opt.Color_Map = jet(64);
    % Span ratio for all axes
    opt.HRatio = [0.8 0.8 0.8];
    unpack_struct(opt);    
end

if exist('h_axes','var')==1
   axes(h_axes); 
else
   figure('Name', 'Hidden entropy vs. Hidden bias for H=3 (Theory)', ...
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
            else
                g = [h(1); 0; h(2)];
            end
          b_h_bound{j}(idx) = -W(:,j)'*Sinv*(b_v + W*g);
       end
    end
    b_h_lim = cell(1,H);
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

set(gca, 'XLim', b_h_lim{1});
set(gca, 'YLim', b_h_lim{2});
set(gca, 'ZLim', b_h_lim{3});
set(gca, 'YDir', 'normal');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
set(gca, 'DataAspectRatioMode', 'manual');
set(gca, 'DataAspectRatio',[1 1 1]);
set(gca, 'FontSize', FontSize);
view(View_Azimuth, View_Elevation);
camproj('perspective');
set(gcf, 'Renderer', 'painters');

xlabel('b_h(1)', 'FontSize', FontSize);
ylabel('b_h(2)', 'FontSize', FontSize);
zlabel('b_h(3)', 'FontSize', FontSize);
opt.h_title = title('Hidden entropy [bits] (theory)');


% Hidden bias which activates two most distant Gaussian components in p(v)
b_h =  (-b_v'*Sinv*W)';
line(b_h(1)*ones(1,2), b_h(2)*ones(1,2), b_h(3)*ones(1,2), ...
    'Color', 'magenta', 'LineWidth', 2, 'Marker', 'o');
text(b_h(1),b_h(2), '${\bf b}_h=-{\bf W}^T{\bf \Sigma}^{-1}{\bf b}_v$', ...
        'BackgroundColor', 'white','Interpreter', 'LaTeX', 'FontSize', FontSize);
for fxd = 1:H
    % ======== Mapping =========
    % fxd = 1 => a = 2; b = 3;   
    % fxd = 2 => a = 1; b = 3;   
    % fxd = 3 => a = 1; b = 2;   
    Dim = 1:H;
    idx =ismember(Dim,fxd);
    ab = num2cell(Dim(~idx));
    [a b] = deal(ab{:});
    % ==========================
    % Truth table: TT = [-1 -1], [-1 +1], [+1 -1], [+1 +1]
    TT = double(dec2bin(0:(2^(H-1))-1)-'1') +  double(dec2bin(0:(2^(H-1))-1)-'0');
    for h_param=TT'
        opt.fxd = fxd;
        opt.a = a;
        opt.b = b;
        opt.b_h_lim = b_h_lim;
        h_const = zeros(H,1);
        h_const(a) = h_param(1);
        h_const(b) = h_param(2);
        opt.h_const = h_const;
        b_h_fxd = -W(:,fxd)'*Sinv*(b_v + W*h_const);
        opt.b_h_fxd = b_h_fxd;
        opt = GBPRBM_Determine_High_HEntropy_Regions_3H(opt);        
        GBPRBM_Plot_High_HEntropy_Regions_3H(opt);
    end   
end

for fxd = 1:H
    opt.fxd = fxd;
    % ======== Mapping =========
    % fxd = 1 => a = 2; b = 3;   
    % fxd = 2 => a = 1; b = 3;   
    % fxd = 3 => a = 1; b = 2;   
    Dim = 1:H;
    idx =ismember(Dim,fxd);
    ab = num2cell(Dim(~idx));
    [a b] = deal(ab{:});
    opt.a = a;
    opt.b = b;
    %% One-Bit Hidden Entropy Region, p(h_j,h_a,h_b) = p(h_j,-h_a,-h_b)
    for h_fxd = [-1 1]
        opt.h_fxd = h_fxd;        
        opt = GBPRBM_Determine_OBHER_Hm1_Antipodes_3H(opt);        
    end
end
Color_Map = opt.Color_Map;
h_lines = zeros(1,3);
C = Color_Map(round(log2(2)/log2(4)*size(Color_Map,1)),:);
h=findall(gca, 'Type','Patch', 'FaceColor', C);
h_lines(1) = h(1);
C = Color_Map(round(log2(3)/log2(4)*size(Color_Map,1)),:);
h=findall(gca, 'Type','Line', 'Color', C);
h_lines(2) = h(1);
C = Color_Map(round(log2(4)/log2(4)*size(Color_Map,1)),:);
h=findall(gca, 'Type','Line', 'MarkerFaceColor', C);
h_lines(3) = h(1);
legend(h_lines, {'1 bit', 'log_2(3) bits', '2 bits'});
