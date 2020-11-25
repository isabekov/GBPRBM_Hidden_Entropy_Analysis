function opt = GBPRBM_Plot_HEntropy_vs_HBias_3H_Experiment(varargin)
% This function plots hidden entropy of the Gaussian-Bipolar Restricted
% Boltzmann Machine (GBPRBM) with three hidden units as a function of
% hidden bias using empirical evaluation of the hidden entropy function.

% Usage:
% Just run without arguments:
%             GBPRBM_Plot_HEntropy_vs_HBias_3H_Experiment
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Plot_HEntropy_vs_HBias_3H_Experiment(opt)


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
    opt.V = 3;
    % Number of hidden units
    opt.H = 3;
    % Model geometry
    opt.W = [9 2  1;
            2 -6 1;
            -3 4 5];
    % Visible bias
    opt.b_v = [8 12 10]';
    % Standard deviations for visible units
    opt.sigma_v = [2 3 4]';
    % Covariance matrix
    opt.S = diag(opt.sigma_v.^2);
    % Inverse of the "covariance matrix"
    opt.Sinv = diag(1./(opt.sigma_v.^2));
    % Span ratio for all axes
    opt.HRatio = [0.8 0.8 0.8];
    % Number of samples used for quantization of b_h(1), b_h(2) and b_h(3)
    opt.N_Samples = [150 151 152];
    % Entropy level at which to plot an isosurface
    opt.Entropy_Level = 0.95;
    % Camera position
    opt.View_Elevation = 38;
    opt.View_Azimuth = -212;
    opt.FontSize = 10;
    % Transparency value of the surfaces
    opt.Alpha_Transparency = 0.5;
    % Unpack structure
    unpack_struct(opt);   
end

if exist('h_axes','var')==1
   axes(h_axes); 
else
   figure('Name', 'Hidden entropy vs. Hidden bias for H=3 (Experiment)', 'NumberTitle', 'Off');
end

if exist('FontSize','var')==0
    FontSize = 10;
end

if exist('Alpha_Transparency', 'var')==0
    Alpha_Transparency = 0.5;
end

Var = {'View_Azimuth', 'View_Elevation'};
% Reset azimuth and elevation values if they don't exist
if  (sum(ismember(Var, who)) ~= length(Var)) 
    View_Azimuth = 138;
    View_Elevation = 36;                
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
opt.b_h_lim = b_h_lim;

b_h = cell(1,H);
for j=1:H
% Sampling, b_h{j} is a row vector
    b_h{j} = linspace(b_h_lim{j}(1), b_h_lim{j}(2), N_Samples(j))';
end
% Index terms
r = (2.^fliplr([0:H-1]))/2; %#ok<NBRAK>
s = sum(r) + 1;
TT = double(dec2bin(0:(2^H)-1)-'1') +  double(dec2bin(0:(2^H)-1)-'0');
B = zeros(1,2^H);

L1 = length(b_h{1});
L2 = length(b_h{2});
L3 = length(b_h{3});

exp_b_h{1} = exp(b_h{1});
exp_b_h{2} = exp(b_h{2});
exp_b_h{3} = exp(b_h{3});

%Nm = zeros(L2,L1,L3,2^H);
p_h = cell(1,2^H);
for k=1:2^H
    % Memeory preallocation
    p_h{k} = zeros(L2,L1,L3);    
end
% Memeory preallocation
Zn  = zeros(L2,L1,L3);

% Truth table for hidden layer, e.g. for H=3, TT is given as:
% MSB LSB
%  h1 h2 h3 -> Order is important!
% [-1 -1 -1;
%  -1 -1 +1;
%  -1 +1 -1;
%   .......
%  +1 +1 +1]
% For every configuration of the hidden layer  (2^H different configurations)
for h=TT' 
   % Decimal index, e.g. for h=[-1 -1]', idx = 1;  for h=[+1 +1]', idx = 4; 
   idx = r*h + s;
   B(idx) = exp(1/2*(2*b_v + W*h)'*Sinv*W*h);
   % C: matrix of size L1 x L2
   C = B(idx)*exp_b_h{2}.^h(2)*exp_b_h{1}'.^h(1);
   % Replicate floor matrix C (dimensions 1,2, i.e. "x", "y")
   % along dimension 3 (i.e. "z")
   D = repmat(C,[1,1,L3]);
   % Replicate column vector exp_b_h{3} for every element on the floor
   % matrix
   E = repmat(exp_b_h{3}.^h(3), [1, L2, L1]);
   % Rotate tensor
   E = shiftdim(E,1);
   % Elementwise product of two tensors
   p_h{idx} = D.*E;
   Zn = Zn + p_h{idx};
end
clear('C', 'D', 'E');
%% Normalization
for k = 1:2^H
    p_h{k} = p_h{k}./Zn;
end

clear('Zn');
Entropy = zeros(L2,L1,L3);
for h=TT'
    idx = r*h + s;
    Entropy = Entropy - p_h{idx}.*log2(p_h{idx});
end
clear p_h;

% Plotting
[X,Y,Z] = meshgrid(b_h{1},b_h{2},b_h{3});

p = patch(isosurface(X,Y,Z,Entropy,Entropy_Level));
isonormals(X,Y,Z,Entropy,p);
set(p,'FaceColor','green','EdgeColor','none');
view(View_Azimuth, View_Elevation);
cameratoolbar;
lighting gouraud;
camlight('headlight');
axis tight;
%set(gca, 'YDir', 'normal');
set(gcf,'Renderer','OpenGL');
%set(gcf,'Renderer','zbuffer');
set(gca, 'FontSize', FontSize);
camproj('perspective');
daspect([1 1 1]);

% Add transparency
alpha(Alpha_Transparency);
h = xlabel('b_h(1)', 'FontSize', FontSize);
% Setting normalized units is necessary in OpenGL rendering mode
set(h, 'Units', 'normalized');

h = ylabel('b_h(2)', 'FontSize', FontSize);
% Setting normalized units is necessary in OpenGL rendering mode
set(h, 'Units', 'normalized');

h = zlabel('b_h(3)', 'FontSize', FontSize);
% Setting normalized units is necessary in OpenGL rendering mode
set(h, 'Units', 'normalized');

opt.h_title = title(['Hidden entropy at ' num2str(Entropy_Level) ' bits (experiment)'], 'FontSize', FontSize);
% Setting normalized units is necessary in OpenGL rendering mode
set(opt.h_title, 'Units', 'normalized');
% Hidden bias which activates two most distant Gaussian components in p(v)
b_h = (-b_v'*Sinv*W)';

%line(b_h(1)*[1 1], b_h(2)*[1 1], get(gca, 'ZLim'));
%line(b_h(1)*[1 1], get(gca, 'YLim'), b_h(3)*[1 1] );
%line(get(gca, 'XLim'), b_h(2)*[1 1], b_h(3)*[1 1] );

line(b_h(1)*[1 1], b_h(2)*[1 1], b_h(3)*[1 1], 'Color', 'magenta', ...
                 'LineWidth', 2, 'Marker', 'sq');
text(b_h(1),b_h(2), '${\bf b}_h=-{\bf W}^T{\bf \Sigma}^{-1}{\bf b}_v$', ...
        'BackgroundColor', 'white','Interpreter', 'LaTeX', 'FontSize', FontSize);