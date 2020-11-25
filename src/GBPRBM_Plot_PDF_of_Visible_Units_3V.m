function opt = GBPRBM_Plot_PDF_of_Visible_Units_3V(varargin)
% This function plots probability density function p(v1,v2,v3) of the
% Gaussian-Bipolar Restricted Boltzmann Machine (GBPRBM) with three visible
% units. The pdf is evaluated at certain points and plotted as an isosurface.

% Usage: 
% Just run without arguments:
%             GBPRBM_Plot_Probability_of_Visible_Units_3V
% Or pass a structure with parameters, look down below for the details
%             GBPRBM_Plot_Probability_of_Visible_Units_3V(opt)


% Author: Altynbek Isabekov
% E-mail: aisabekov [at] ku.edu.tr

%% ================== Setting Parameters of the GBPRBM Model ==============
if ~isempty(varargin)
    if length(varargin) >=1
        if isstruct(varargin{1})
            opt = varargin{1};
            % Unpack structure
            unpack_struct(opt); 
            if V~=3
                error('Number of visible units is not equal to 3. Quitting.')
            end
            if size(W,1)~=3
                error('Matrix W is not of size 3xH. Quitting.')
            end
            if length(sigma_v)~=3
                error('Vector sigma_v is not of size 3x1. Quitting.')
            end            
            if length(b_v)~=3
                error('Vector b_v is not of size 3x1. Quitting.')
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
    % W = [ 10.0000   -1.3892    6.8200;
    %          0    7.8785   -7.3135  ];
    %     h_x = [0; x2; x3];
    %     b_h_act_x = -W(:,1)'*Sinv*(b_v + alpha*W*h_x);
    %     h_y = [y1; 0; y3];
    %     b_h_act_y = -W(:,2)'*Sinv*(b_v + alpha*W*h_y);
    %     h_z = [z1; z2; 0];
    %     b_h_act_z = -W(:,3)'*Sinv*(b_v + alpha*W*h_z);
    % Span ratio for all axes
    opt.VRatio = [0.8 0.8 0.8];
    % Number of samples used to quantize X-axis
    opt.Nv1 = 190;
    % Number of samples used to quantize Y-axis
    opt.Nv2 = 180;
    % Number of samples used to quantize Z-axis
    opt.Nv3 = 170;
    % Camera position
    opt.View_Elevation = 38;
    opt.View_Azimuth = -212;
    % Put labels (needed for GBPRBM_Visible_Units_Plot_Geometry)
    opt.put_labels = true;    
    % Transparency value of the surfaces
    opt.Alpha_Transparency = 0.5;
    % Colors for drawing geometry of the model
    opt.Color = {[0, 0.3, 0.6],'red', [0, 0.7, 0], [0.8706, 0.4902, 0], 'green'};
    % Ratio of the pdf's maximum value at which isosurface will be plotted
    opt.Probability_Level = [0.0005 0.2  0.5 0.9];
    % Unpack structure
    unpack_struct(opt);    
    h = [+1 -1];    
    b_h =  (-b_v'*Sinv*W)';
    b_h(1) = -W(:,1)'*Sinv*(b_v + (h(2))*W(:,2)); 
    b_h(2) = -W(:,2)'*Sinv*(b_v + (h(1))*W(:,1));
    opt.b_h =  b_h;
end

if exist('h_axes','var')==1
   axes(h_axes); 
else
   figure('Name', 'GBPRBM Probability of visible units, p(v1,v2,v3)', 'NumberTitle', 'Off');
end

if exist('plot_decision_regions','var')==0
    plot_decision_regions = true;
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

if exist('Probability_Level', 'var')==0
    Probability_Level = [0.0005 0.2  0.5 0.9];
end
%% ======= Determination of the range of  X-, Y- and Z-axes ======
% Range is necessary to set proper XLim, YLim and ZLim values.
% If all variables listed in "Var" exist, no need to determine boundaries -
% they are already given.
Var = {'XLimLow ', 'XLimHigh', 'YLimLow', 'YLimHigh', 'ZLimLow', 'ZLimHigh'};
if  ~(sum(ismember(Var, who)) == length(Var)) 
    % If variables listed in "Var" do not exist then calculated boundaries
    [opt] = GBPRBM_Visible_Units_Span(opt);
    % Cast calculated span for each of the axes
    XLimLow  = opt.v(1).AxesLimMin;
    XLimHigh = opt.v(1).AxesLimMax;
    YLimLow  = opt.v(2).AxesLimMin;
    YLimHigh = opt.v(2).AxesLimMax;
    ZLimLow  = opt.v(3).AxesLimMin;
    ZLimHigh = opt.v(3).AxesLimMax;
end

% PDF p(v) will be computed for each point on the coordinate grid "v". 
v = cell(1,V);
% Linearly spaced coordinate vectors v{1} and v{2} should be of column type
v{1} = linspace(XLimLow, XLimHigh,  Nv1)';
v{2} = linspace(YLimLow, YLimHigh,  Nv2)';
v{3} = linspace(ZLimLow, ZLimHigh,  Nv3)';
% Truth table for the hidden neurons
% Superimpose these two matrices
%         Values: "0" and "-1"            Values: "+1" and "0"
TT = double(dec2bin(0:(2^H)-1)-'1') +  double(dec2bin(0:(2^H)-1)-'0');

%% ======================== Computation of p(v) ========================
a_max = -Inf;
% For every configuration of the hidden layer  (2^H different configurations)
for g=TT' % Hidden vector "g" of size: (H x 1)
    a = 1/2*(2*b_v + W*g)'*Sinv*W*g + b_h'*g;
    if a > a_max
       a_max = a; 
    end   
end
% Number of points in the "v1" coordinate grid
Lv1 = length(v{1});
% Number of points in the "v2" coordinate grid
Lv2 = length(v{2});
% Number of points in the "v2" coordinate grid
Lv3 = length(v{3});
% A template for the "A" part, independent of the hidden units
Vis_Only = cell(1,V);             
% A template for the "B" part, which is a function of both visible units "v"
% and hidden parameters ("W", "b_h")
Vis_Hid = cell(1,V);

for i=1:V
    % The lengths of v{1} and v{2} are different!
    Vis_Only{i} = -(v{i}-b_v(i)).^2/(2*(sigma_v(i))^2);   
    % Vis_Hid is a matrix of size Lv2 x H, obtained through an outer product
    Vis_Hid{i}  = v{i}./(sigma_v(i))^2; 
end
% Create a grid
% Replicate column vector Vis_Only{1} along dimension v2 and v3
% Replicate row    vector Vis_Only{2} along dimension v1 and v3
Asum = (  repmat(Vis_Only{1} , [ 1 , Lv2, Lv3]) + ...
        + repmat(Vis_Only{2}', [Lv1,  1 , Lv3]));

Temp = repmat(Vis_Only{3}, [1,Lv1,  Lv2]);
Asum = Asum + shiftdim(Temp, 1);

% The "B" part, which is a function of both visible units "v"
% and hidden parameters ("W", "b_h")

% Initial sum
ExpB=0;
% For every configuration of the hidden layer  (2^H different configurations)
for g=TT'   % Hidden vector g of size: (H x 1)
   Coef = W*g;
   Bsum = repmat(Vis_Hid{1}*Coef(1) , [ 1 , Lv2, Lv3]) + ...
        + repmat(Vis_Hid{2}'*Coef(2), [Lv1, 1 ,  Lv3 ]);
   Temp = repmat(Vis_Hid{3}*Coef(3), [1,   Lv1, Lv2]);
   ExpB = ExpB + exp(Asum + Bsum + shiftdim(Temp, 1) + b_h'*g - a_max);
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

% Plotting
[X,Y,Z] = meshgrid(v{1},v{2},v{3});
Level = max(max(max(p_v)))*Probability_Level;
L = length(Level);
Color = mat2cell(rand(L,3), ones(1,L), 3);
Color{1} = [ 0    0.6000    0.6000] ; %'black';
Color{2} = 'red';
Color{3} = 'green';
Color{4} = 'blue';

for k=1:L
    fv = isosurface(X,Y,Z,ipermute(p_v,[2 1 3]),Level(k));
    p = patch(fv);
    isonormals(X,Y,Z,ipermute(p_v,[2 1 3]),p);
    set(p,'FaceColor',Color{k},'EdgeColor','none')    
    hold on;
end

hold off;
str = sprintf(['p(v_1,v_2,v_3) sliced at ' repmat('%1.2e,',1,L-1) ' and %1.2e'], Level);
title(str, 'FontSize', FontSize);
xlabel('v_1', 'FontSize', FontSize);
ylabel('v_2', 'FontSize', FontSize);
zlabel('v_3', 'FontSize', FontSize);
set(gca, 'XLim', [XLimLow XLimHigh]);
set(gca, 'YLim', [YLimLow YLimHigh]);
set(gca, 'ZLim', [ZLimLow ZLimHigh]);

%% Plot geometry of the model (location of the centroids)
GBPRBM_Visible_Units_Plot_Geometry(opt);

%% Plot decision regions for p(h_j|v) for j=1...H
if plot_decision_regions==true
   GBPRBM_Plot_Decision_Region_for_Visible_Units(opt);
end

% Renderer
set(gcf,'Renderer','OpenGL');
%camproj('perspective');
% Add transparency to the plot (works only with OpenGL renderer)
alpha(Alpha_Transparency);
lighting gouraud;
camlight('headlight');
view(View_Azimuth, View_Elevation);
daspect([1 1 1]);
